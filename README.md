# fisseq-qc-filter

A CLI tool for quality-filtering CellProfiler cell profiles from FISSEQ
experiments. It removes low-confidence barcode reads, rare barcodes, and
under-represented variants, then writes clean Parquet outputs for downstream
analysis.

## Installation

Requires Python >= 3.13. Install with [uv](https://docs.astral.sh/uv/):

```bash
uv sync
```

This installs the `fisseq-qc-filter` command into the project's virtual
environment.

## Usage

```bash
fisseq-qc-filter cell_files=<path> output_dir=<path> [overrides...]
```

### Required arguments

| Argument | Description |
|---|---|
| `cell_files` | Path to a single cell profile file, or a list of paths |
| `output_dir` | Directory to write output files (created if absent) |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `output_root` | `null` | Prefix for output filenames. If set, files are named `{output_root}.filtered_cells.parquet` etc. |
| `bc_threshold` | `10` | Minimum number of cells a barcode must appear in to be retained |
| `variant_bc_threshold` | `4` | Minimum number of passing barcodes a variant must have to be retained |
| `edit_distance_threshold` | `1` | Maximum edit distance from the reference barcode to accept |
| `barcode_col_name` | `upBarcode` | Name of the barcode column in the input files |

### Examples

Filter a single file:

```bash
fisseq-qc-filter \
    cell_files=/data/plate1_cells.parquet \
    output_dir=/results/plate1
```

Filter multiple files together:

```bash
fisseq-qc-filter \
    'cell_files=[/data/plate1.parquet,/data/plate2.parquet]' \
    output_dir=/results/combined
```

Override QC thresholds:

```bash
fisseq-qc-filter \
    cell_files=/data/plate1.parquet \
    output_dir=/results/strict \
    bc_threshold=20 \
    variant_bc_threshold=8 \
    edit_distance_threshold=0
```

Add a prefix to output filenames:

```bash
fisseq-qc-filter \
    cell_files=/data/plate1.parquet \
    output_dir=/results \
    output_root=plate1
# writes: plate1.filtered_cells.parquet, plate1.barcode_counts.parquet, ...
```

## Input format

Input files must be CSV or Parquet (`.csv`, `.parquet`, `.parq`, `.pq`).
Each row represents one cell. The following columns are expected:

| Column | Description |
|---|---|
| `editDistance` | Edit distance of the decoded barcode from the nearest reference |
| `aaChanges` | Amino acid change(s) for the variant carried by this cell |
| `variantType` | Classification of the variant (e.g. missense, synonymous) |
| `upBarcode` (configurable) | The upstream barcode sequence |

CellProfiler feature columns are identified by starting with an uppercase
letter and containing an underscore (e.g. `Cells_AreaShape_Area`). All other
columns are dropped before output.

### A word on file formats

In general, **using CSV input files will result in a 5x-20x increase in runtime** due to parsing overhead.
The QC filtering is run using lazy evaluation, meaning that the pipeline scans the input files as it goes rather than reading the whole file into memory at once.
This significantly decreases the memory footprint.
However, in the case of CSV files this means that input file has to be constantly parsed as the pipeline runs, leading to a significant runtime penalty.

## QC pipeline

Filters are applied in this order:

1. **Edit distance filter** — retains cells where `editDistance <=
   edit_distance_threshold`. Removes cells where the barcode could not be
   confidently decoded.

2. **Barcode-level filter** — counts how many cells carry each barcode.
   Barcodes with fewer than `bc_threshold` cells are dropped. This removes
   rare or spurious barcode integrations.

3. **Variant-level filter** — counts how many passing barcodes each variant
   has. Variants represented by fewer than `variant_bc_threshold` barcodes
   are dropped. This ensures each variant has enough independent measurements
   for reliable downstream statistics.

## Outputs

Three Parquet files are written to `output_dir`:

| File | Contents |
|---|---|
| `filtered_cells.parquet` | Cells passing all QC filters, with CellProfiler features and metadata columns |
| `barcode_counts.parquet` | Per-barcode cell counts and QC pass/fail status |
| `variants_per_barcode.parquet` | Per-variant barcode counts and QC pass/fail status |

If `output_root` is set, filenames are prefixed with `{output_root}.`.

Metadata columns added to `filtered_cells.parquet`:

| Column | Description |
|---|---|
| `meta_source_file` | Path to the input file the cell came from |
| `meta_source_file_idx` | Row index of the cell within its source file |
| `meta_barcode` | Barcode sequence |
| `meta_aa_changes` | Amino acid changes |
| `meta_edit_distance` | Edit distance at time of filtering |
