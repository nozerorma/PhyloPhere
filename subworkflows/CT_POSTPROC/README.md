# CT Post-Processing Workflow

This module provides automated post-processing and characterization of CAAS (Context-Aware Amino Acid Substitutions) discovery results from the PhyloPhere CT workflow.

## Overview

The CT Post-Processing workflow consists of two main components:

1. **Cluster Filtering** (`ctpp_clustfilter.nf`): Identifies and filters high-density clusters of CAAS positions using parameter sweep or single-parameter modes
2. **Characterization** (`ctpp_characterization.nf`): Generates visual reports and statistical summaries of CAAS discovery results

## File Structure

```
subworkflows/CT_POSTPROC/
├── ctpp_clustfilter.nf          # Cluster filtering subworkflow
├── ctpp_characterization.nf     # Characterization & reporting subworkflow
└── local/
    ├── filter_caas_clusters-param.py    # Python filtering script
    └── 5.2.CT_postproc.Rmd             # R markdown reporting skeleton
```

## Cluster Filtering

### Modes

#### Exploratory Mode (default)
Performs a parameter sweep using all combinations of:
- **minlen**: {2, 3, 4, 10} (minimum cluster interval length)
- **maxcaas**: {0.6, 0.7, 0.8} (maximum CAAS density threshold)

Creates 12 parameter combinations with outputs organized in subdirectories:
```
{outdir}/ct_postproc/filter_exploratory/
├── minlen2_maxcaas60/
│   └── discovery.filtered.minlen2.maxcaas60.tsv
├── minlen2_maxcaas70/
│   └── discovery.filtered.minlen2.maxcaas70.tsv
├── minlen2_maxcaas80/
│   └── discovery.filtered.minlen2.maxcaas80.tsv
├── minlen3_maxcaas60/
...
└── discarded_summary.tsv        # Consolidated summary
```

#### Filter Mode
Runs a single filtering pass with user-specified parameters:
- `--filter_minlen`: Target minlen value
- `--filter_maxcaas`: Target maxcaas value

Output structure:
```
{outdir}/ct_postproc/filter_selected/
├── discovery.filtered.minlenX.maxcaasY.tsv
└── discarded_summary.tsv
```

### Output Format

Each filtered CAAS file is tab-separated with columns:
- **Gene**: Gene identifier (from discovery input)
- **Position**: CAAS position number
- **ClusteringFlag**: "Good" (retained) or "Discarded" (in high-density cluster)

### Consolidated Summary (`discarded_summary.tsv`)

Tab-separated file with columns:
- **File**: Input CAAS filename
- **DiscardedCount**: Number of positions marked as "Discarded"
- **Minlen**: Minimum cluster interval length (exploratory mode only)
- **Maxcaas**: Maximum density threshold value (exploratory mode only)

## Characterization & Reporting

Generates HTML reports and supplementary data files analyzing CAAS discovery results.

### Inputs

1. **Discovery file** (`.tab`): Original CAAS discovery output (Gene, Position, [other columns])
2. **Filter summary** (`discarded_summary.tsv`): Output from cluster filtering
3. **Gene length file** (tab-separated, no header):
   - Column 1: Gene name
   - Column 2: Gene length (integer)
4. **Filter directory**: Directory containing filtered results

### Output

Reports generated in `{outdir}/ct_postproc/reports/`:
- HTML report with visualizations
- TSV data tables with gene categories and statistics
- PNG plots for parameter analysis and distributions

## Usage Examples

### Exploratory Mode (Parameter Sweep)

```bash
nextflow run main.nf \
  --ct_postproc \
  --discovery_input "results/discovery.tab" \
  --gene_length_file "data/gene_lengths.tsv" \
  --generate_reports
```

Default parameters:
- `--caas_postproc_mode exploratory`
- `--minlen_values "2,3,4,10"`
- `--maxcaas_values "0.6,0.7,0.8"`
- `--outdir "Out"`

### Filter Mode (Single Parameter Set)

```bash
nextflow run main.nf \
  --ct_postproc \
  --caas_postproc_mode filter \
  --discovery_input "results/discovery.tab" \
  --filter_minlen 3 \
  --filter_maxcaas 0.7 \
  --gene_length_file "data/gene_lengths.tsv" \
  --generate_reports
```

### Cluster Filtering Only (No Reports)

```bash
nextflow run main.nf \
  --ct_postproc \
  --discovery_input "results/discovery.tab" \
  --generate_reports false
```

## Configuration

Default configuration in `conf/ct_postproc.config`:

```groovy
params {
    // Input/Output
    discovery_input = ""                    // Path to discovery.tab file (REQUIRED)
    gene_length_file = ""                   // Path to gene_length TSV file (REQUIRED if generate_reports=true)
    postproc_outdir = "${params.outdir}/ct_postproc"
    
    // Processing mode
    caas_postproc_mode = "exploratory"      // "exploratory" or "filter"
    
    // Exploratory mode parameters
    minlen_values = "2,3,4,10"              // Comma-separated minlen values
    maxcaas_values = "0.6,0.7,0.8"          // Comma-separated maxcaas values
    
    // Filter mode parameters
    filter_minlen = 3
    filter_maxcaas = 0.7
    
    // Processing options
    generate_reports = true
    verbose = false
}
```

## Algorithm Details

### Cluster Detection (filter_caas_clusters-param.py)

For each gene, the algorithm evaluates all possible position intervals `[start, end]` where:
- Interval span ≥ `minlen`
- Density = count_positions / span ≥ `maxcaas`

Positions within any interval exceeding the density threshold are marked as "Discarded".

**Example:**
```
Positions: [10, 11, 12, 50]
minlen=3, maxcaas=0.7

Interval [10, 12]: span=3, count=3, density=1.0 ≥ 0.7 → Discard [10, 11, 12]
Interval [10, 50]: span=41, count=4, density=0.098 < 0.7 → Keep
```

### Outlier Detection (Characterization)

Two complementary methods:

1. **Extreme Genes** (top 1% by CAAS density): Density = CAAS_count / Gene_length
2. **Dubious Genes** (IQR outliers): CAAS count > Q3 + 3×IQR

Genes are categorized as:
- **Both**: Extreme AND Dubious
- **Extreme**: Top 1% density only
- **Dubious**: IQR outlier only
- **Normal**: Neither

## Citation

Original CAAS methodology: [Insert publication reference]

Python filtering script based on: [seq_test.sh / filter_caas_clusters-param.py]

R markdown report template: 5.2.CT_postproc.Rmd (original manuscript analysis)

## Support

For issues or questions regarding the CAAS post-processing workflow, please refer to:
- Main pipeline documentation: [PhyloPhere README]
- Script documentation: inline comments in Python/R scripts
- Configuration examples: `conf/ct_postproc.config`
