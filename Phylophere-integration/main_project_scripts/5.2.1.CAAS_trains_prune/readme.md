# CAAS Train Filtering and Parameter Sweep

This repository contains two scripts designed to filter CAAS (Clustering Analysis of Allelic Sites) data based on clustering density and perform a parameter sweep to evaluate different filtering thresholds:

1. **`filter_caas_clusters-param.py`**: A Python script that processes a CAAS input file, identifies high-density clusters of positions, and outputs a filtered file with positions marked as "Good" or "Discarded".
2. **`seq_test.sh`**: A Bash script that runs the Python script with a grid of parameters (`minlen` and `maxcaas`), tracks runtime, and generates a summary of discarded positions per output file.

## Purpose

The `filter_caas_clusters-param.py` script processes a tab-separated `.caas` file containing gene and position data, identifying clusters where the density of positions exceeds a specified threshold (`maxcaas`) over a minimum interval length (`minlen`). Positions in high-density clusters are marked as "Discarded", while others are marked as "Good". The script includes input validation, verbose logging, and efficient pruning logic.

The `seq_test.sh` script automates running the Python script with multiple combinations of `minlen` and `maxcaas`, organizes output files, and generates a summary file (`discarded_summary.tsv`) with the count of discarded positions and runtime for each parameter set.

## Requirements

- **Python 3.6+** (for `filter_caas_clusters-param.py`):
  - `pandas`: For data manipulation and CSV handling.
  - `numpy`: For efficient sorting and unique position extraction.
- **Bash** (for `seq_test.sh`):
  - Standard Unix utilities: `awk`, `find`, `bc`.
- **Input File**: A tab-separated `.caas` file with at least two columns: `Gene` (string) and `Position` (integer).

Install Python dependencies:
```bash
pip install pandas numpy
```

## Usage

### 1. `filter_caas_clusters-param.py`

**Purpose**: Filters positions in a `.caas` file based on clustering density.

**Command**:
```bash
./filter_caas_clusters-param.py -i input.caas [-c maxcaas] [-l minlen] [-v]
```

**Arguments**:
- `-i, --inputfile`: Path to the input `.caas` file (required).
- `-c, --maxcaas`: Maximum density threshold (0 to 1, default: 0.7).
- `-l, --minlen`: Minimum interval length (≥1, default: 3).
- `-v, --verbose`: Enable debug-level logging (optional).

**Input File Format**:
- Tab-separated file with at least `Gene` and `Position` columns.
- Example:
  ```
  Gene    Position
  GeneA   100
  GeneA   101
  GeneB   200
  ```

**Output**:
- A tab-separated file named `<input>.filtered.minlenX.maxcaasY.tsv` with columns: `Gene`, `Position`, `ClusteringFlag` ("Good" or "Discarded").
- A log file named `<input>.minlenX.maxcaasY.log` with processing details (DEBUG level if `--verbose` is used).

**Example**:
```bash
./filter_caas_clusters-param.py -i data.caas -c 0.5 -l 5 -v
```
- Output: `data.filtered.minlen5.maxcaas50.tsv`, `data.minlen5.maxcaas50.log`

### 2. `seq_test.bash`

**Purpose**: Runs `filter_caas_clusters-param.py` with a grid of `minlen` and `maxcaas` values, tracks runtime, and summarizes discarded positions.

**Command**:
```bash
./seq_test.bash input_file.caas [output_dir]
```

**Arguments**:
- `input_file.caas`: Path to the input `.caas` file (required).
- `output_dir`: Directory to store output files and logs (optional; defaults to current directory).

**Parameter Grid**:
- `minlen`: `[3, 10]`
- `maxcaas`: `[0.5, 0.6, 0.7, 0.8]`
- Edit `seq_test.bash` to modify these values.

**Output**:
- Filtered `.tsv` files for each parameter combination.
- Log files for each run (in `output_dir/logs` if specified).
- A summary file `discarded_summary.tsv` with columns: `File`, `DiscardedCount`, `Runtime`.

**Example**:
```bash
./seq_test.bash data.caas output
```
- Runs filtering for all combinations of `minlen` and `maxcaas`.
- Outputs `.tsv` and `.log` files to `output/` and `output/logs/`.
- Creates `output/discarded_summary.tsv` with counts and runtimes.

**Sample `discarded_summary.tsv`**:
```
File                            DiscardedCount  Runtime
data.filtered.minlen3.maxcaas50.tsv  10          00:00:05
data.filtered.minlen3.maxcaas60.tsv  8           00:00:04
```

## Key Features

- **Dynamic Pruning**: Identifies clusters with density ≥ `maxcaas` over spans ≥ `minlen`, merging overlapping intervals.
- **Input Validation**: Ensures the input file is tab-separated, has required columns, and contains valid integer positions.
- **Logging**: Outputs detailed logs (INFO or DEBUG) to both file and console.
- **Parameter Sweep**: Automates testing multiple parameter combinations and summarizes results.
- **Efficient Processing**: Uses `numpy` for sorting and `pandas` for data manipulation.

## Notes

- The Python script skips redundant interval checks to optimize performance without losing information.
- Debug logs include discarded positions per gene when `--verbose` is used.
- The Bash script tracks runtime for each parameter set and reports total runtime on exit.
- Modify the `MINLEN_VALUES` and `MAXCAAS_VALUES` arrays in `seq_test.bash` to adjust the parameter grid.

## Example Workflow

1. Prepare a `.caas` file (e.g., `data.caas`).
2. Run the parameter sweep:
   ```bash
   ./seq_test.bash data.caas results
   ```
3. Check outputs in `results/`:
   - Filtered `.tsv` files (e.g., `data.filtered.minlen3.maxcaas70.tsv`).
   - Log files in `results/logs/`.
   - Summary in `results/discarded_summary.tsv`.

## Troubleshooting

- **Input file errors**: Ensure the `.caas` file is tab-separated with `Gene` and `Position` columns.
- **Missing script**: Verify `filter_caas_clusters-param.py` is in the same directory as `seq_test.bash`.
- **Non-integer positions**: The Python script will exit with an error if the `Position` column contains non-integer values.
