import logging
import pathlib
from os import PathLike
from typing import Iterable

import hydra
import polars as pl
from omegaconf import DictConfig


def get_barcode_counts(lf: pl.LazyFrame, cfg: DictConfig) -> pl.LazyFrame:
    """Count cells per barcode and flag barcodes meeting the threshold.

    Groups by ``meta_barcode``, counts occurrences, and adds a
    ``barcode_ok`` column (non-null when count >= ``cfg.bc_threshold``).
    Retains the first ``meta_aa_changes`` and ``meta_variant_type``
    per barcode group.

    Args:
        lf: Cell-level lazy frame containing ``meta_barcode``,
            ``meta_aa_changes``, and ``meta_variant_type`` columns
            (as produced by ``filter_columns``).
        cfg: Hydra config supplying ``bc_threshold``.

    Returns:
        Lazy frame with one row per barcode, including ``count``,
        ``meta_aa_changes``, ``meta_variant_type``, and ``barcode_ok``.
    """
    barcode_lf = (
        lf.group_by("meta_barcode")
        .agg(
            [
                pl.len().alias("count"),
                pl.col("meta_aa_changes").first(),
                pl.col("meta_variant_type").first(),
            ]
        )
        .with_columns(
            pl.when(pl.col("count") >= cfg.bc_threshold)
            .then(pl.col("count"))
            .otherwise(None)
            .alias("barcode_ok")
        )
    )

    return barcode_lf


def get_barcodes_per_variant(cells_lf: pl.LazyFrame, cfg: DictConfig) -> pl.LazyFrame:
    """Count distinct barcodes per variant and flag variants meeting threshold.

    Groups by ``meta_aa_changes``, counts barcodes, and adds a
    ``variant_barcode_count_ok`` column (non-null when barcode count
    >= ``cfg.variant_bc_threshold``).

    Args:
        cells_lf: Cell-level lazy frame containing ``meta_barcode``
            and ``meta_aa_changes`` columns (as produced by
            ``filter_columns``).
        cfg: Hydra config supplying ``variant_bc_threshold``.

    Returns:
        Lazy frame with one row per variant, including
        ``barcode_count`` and ``variant_barcode_count_ok``.
    """
    barcodes_per_variant_lf = (
        cells_lf.group_by("meta_aa_changes")
        .agg(
            [
                pl.col("meta_barcode").n_unique().alias("barcode_count"),
            ]
        )
        .with_columns(
            pl.when(pl.col("barcode_count") >= cfg.variant_bc_threshold)
            .then(pl.col("barcode_count"))
            .otherwise(None)
            .alias("variant_barcode_count_ok")
        )
    )

    return barcodes_per_variant_lf


def add_qc_queries(
    lf: pl.LazyFrame, cfg: DictConfig
) -> tuple[pl.LazyFrame, pl.LazyFrame, pl.LazyFrame]:
    """Apply edit-distance, barcode-level, and variant-level QC filters.

    Filters are applied sequentially:
    1. Retain rows with ``editDistance`` <= ``cfg.edit_distance_threshold``.
    2. Remove barcodes below ``cfg.bc_threshold`` cell count.
    3. Remove variants below ``cfg.variant_bc_threshold``
       barcode count.

    Args:
        lf: Cell-level lazy frame to filter.
        cfg: Hydra config supplying QC thresholds and column names.

    Returns:
        A tuple of (filtered_lf, barcode_count_lf,
        variants_per_barcode_lf) where the latter two contain
        the intermediate QC summary frames.
    """
    logging.info("Adding edit distance QC query")
    lf = lf.filter(pl.col("meta_edit_distance") <= cfg.edit_distance_threshold)

    logging.info("Adding Barcode Level QC query")
    barcode_count_lf = get_barcode_counts(lf, cfg)
    lf = lf.join(
        barcode_count_lf.filter(pl.col("barcode_ok").is_not_null()).select(
            "meta_barcode"
        ),
        on="meta_barcode",
        how="inner",
    )

    logging.info("Adding Variant Level QC query")
    variants_per_barcode_lf = get_barcodes_per_variant(lf, cfg)
    lf = lf.join(
        variants_per_barcode_lf.filter(
            pl.col("variant_barcode_count_ok").is_not_null()
        ).select("meta_aa_changes"),
        on="meta_aa_changes",
        how="inner",
    )

    return lf, barcode_count_lf, variants_per_barcode_lf


def read_file(cell_file_path: pathlib.Path) -> pl.LazyFrame:
    """Read a single cell file into a lazy frame.

    Supports CSV and Parquet formats (extensions ``.csv``,
    ``.parquet``, ``.parq``, ``.pq``). Adds two metadata columns:
    ``meta_source_file`` (file path as a string) and
    ``meta_origin_file_idx`` (row index within the source file).

    Args:
        cell_file_path: Path to the cell data file.

    Returns:
        Lazy frame of the file contents with metadata columns.
    """
    logging.info("Scanning file %s", cell_file_path)

    if cell_file_path.suffix == ".csv":
        lf = pl.scan_csv(cell_file_path)
    elif cell_file_path.suffix in [".parquet", ".parq", ".pq"]:
        lf = pl.scan_parquet(cell_file_path)

    lf = lf.with_columns(
        pl.lit(str(cell_file_path)).alias("meta_source_file"),
        pl.row_index().alias("meta_origin_file_idx"),
    )

    return lf


def combine_cell_files(cell_files: Iterable[PathLike]) -> pl.LazyFrame:
    """Read and concatenate multiple cell files into one lazy frame.

    Args:
        cell_files: Iterable of paths to cell data files.

    Returns:
        Concatenated lazy frame of all input files.
    """
    return pl.concat([read_file(pathlib.Path(cell_file)) for cell_file in cell_files])


def filter_columns(lf: pl.LazyFrame, cfg: DictConfig) -> pl.LazyFrame:
    """Retain only the columns needed for QC and output.

    This function can be used to reduce memory usage by dropping
    unnecessary columns before writing output files.

    Args:
        lf: Lazy frame containing all columns.
        cfg: Hydra config supplying the barcode column name.

    Returns:
        Lazy frame with only the necessary columns retained.
    """
    lf = lf.with_columns(
        pl.col("aaChanges").alias("meta_aa_changes"),
        pl.col("variantType").alias("meta_variant_type"),
        pl.col("editDistance").alias("meta_edit_distance"),
        pl.col(cfg.barcode_col_name).alias("meta_barcode"),
    )

    schema_names = lf.collect_schema().names()
    cell_profiler_columns = [
        col for col in schema_names if len(col) > 0 and col[0].isupper() and "_" in col
    ]
    meta_columns = [col for col in schema_names if col.startswith("meta_")]

    return lf.select(pl.col(meta_columns + cell_profiler_columns))


def configure_logging(log_path: PathLike) -> None:
    """Configure root logger to write to a file and stdout.

    Sets log level to INFO and formats messages as
    ``timestamp - level - message``.

    Args:
        log_path: Path to the log file (created if absent).
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )


def log_config(cfg: DictConfig) -> None:
    """Log all config key-value pairs at INFO level.

    Args:
        cfg: Hydra config to log.
    """
    logging.info("Configuration:")
    for key, value in cfg.items():
        logging.info(f"  {key}: {value}")


@hydra.main(version_base=None, config_path="conf", config_name="config")
def main(cfg: DictConfig) -> None:
    """CLI entrypoint driven by Hydra configuration.

    Required overrides (no default in config.yaml):
        cell_files: path or list of paths to cell data files.
        output_dir: directory to write output files.

    Optional overrides (defaults set in config.yaml):
        output_root (default: null — no prefix on output filenames)
        bc_threshold (default: 10)
        variant_bc_threshold (default: 4)
        barcode_col_name (default: "upBarcode")

    Example usage::

        fisseq-qc-filter \\
            cell_files=[a.parquet,b.parquet] \\
            output_dir=/results
    """
    output_dir = pathlib.Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{cfg.output_root}." if cfg.output_root is not None else ""
    log_path = output_dir / f"{prefix}.fisseq_qc_filter.log"
    configure_logging(log_path)
    log_config(cfg)

    cell_files = (
        list(cfg.cell_files)
        if hasattr(cfg.cell_files, "__iter__") and not isinstance(cfg.cell_files, str)
        else [cfg.cell_files]
    )

    combined_lf = combine_cell_files(cell_files)
    combined_lf = filter_columns(combined_lf, cfg)
    combined_lf, barcode_count_lf, variants_per_barcode_lf = add_qc_queries(
        combined_lf, cfg
    )

    logging.info("Writing output files to %s", output_dir)
    for name, lf in [
        ("filtered_cells", combined_lf),
        ("barcode_counts", barcode_count_lf),
        ("variants_per_barcode", variants_per_barcode_lf),
    ]:
        logging.info("Writing %s", name)
        lf.sink_parquet(output_dir / f"{prefix}{name}.parquet")
