# CAAP Mode in PhyloPhere

## Overview

CAAP (Convergent Amino Acid Properties) mode extends the classical CAAS detection by identifying convergence based on physicochemical properties rather than exact amino acid identity. This allows detection of functional convergence even when different amino acids with similar properties are substituted.

## What's New

### Grouping Schemes

CAAP implements 5 grouping schemes (GS0-GS4) that classify amino acids by different property sets:

- **GS0 (Identity)**: Each amino acid is its own group
  - Validates against classical CAAS (should produce identical results)
  - 20 groups: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y

- **GS1 (Physicochemical)**: 6 groups based on chemical properties
  - CV | AGPS | NDQE | RHK | ILMFWY | T

- **GS2 (Charge/Size)**: 7 groups based on charge and molecular size
  - C | AGV | DE | NQHW | RK | ILFP | YMTS

- **GS3 (Hydrophobicity)**: 6 groups based on hydrophobic properties
  - C | AGPST | NDQE | RHK | ILMV | FWY

- **GS4 (Fine-grained)**: 12 groups for detailed biochemical classification
  - C | AILV | ST | NQ | DE | RH | G | P | K | M | F | WY

### Output Format

#### Discovery Output
When `--caap_mode` is enabled, discovery output includes:
- **Mode column**: "CAAP" for property-based hits
- **CAAP_Group column**: The grouping scheme (GS0-GS4) that detected the convergence
- **One row per matching scheme**: A position may appear multiple times if detected by different schemes

Example:
```
Gene    Mode    CAAP_Group    Trait          Position    Substitution    ...
BRCA1   CAAP    GS0           config_np.tab  745         VEG/EEE        ...
BRCA1   CAAP    GS1           config_np.tab  745         132/333        ...
BRCA1   CAAP    GS2           config_np.tab  745         BBC/CCC        ...
```

#### Bootstrap Output
Bootstrap tests ALL grouping schemes for each discovered position:
```
Position        CAAP_Group    Count    Cycles    EmpiricalPval
BRCA1@745       GS0          1232     100000    0.01232
BRCA1@745       GS1          1232     100000    0.01232
BRCA1@745       GS2          1228     100000    0.01228
```

## Usage

### In Nextflow

Add the `--caap_mode` parameter to your Nextflow command:

```bash
nextflow run main.nf \
    --ct_tool "discovery,bootstrap" \
    --alignment "Data/protein_alignments/**/*" \
    --traitfile "Data/CAAS_traitfiles/traitfile.tab" \
    --paired_mode true \
    --max_conserved "1" \
    --caap_mode true \
    --outdir "Out/Phylophere_CAAP"
```

### In Standalone CAAStools

```bash
# Discovery
./ct discovery \
    -a alignment.phy \
    -t traitfile.tab \
    -o output.tsv \
    --paired_mode \
    --max_conserved 1 \
    --caap_mode

# Bootstrap
./ct bootstrap \
    -a alignment.phy \
    -t traitfile.tab \
    -s resample_dir/ \
    -o bootstrap.tsv \
    --discovery discovery_output.tsv \
    --paired_mode \
    --max_conserved 1 \
    --caap_mode
```

## Configuration

### In nextflow.config

```groovy
params {
    caap_mode = false  // Set to true to enable CAAP mode
}
```

### Override at runtime

```bash
nextflow run main.nf --caap_mode true
```

## Interpretation

### Scheme-Specific Results

Different grouping schemes may detect different positions:
- **GS0 results = Classical CAAS**: Validates the implementation
- **GS1-GS4 results**: May find additional convergence based on property classes
- **Multiple schemes per position**: Indicates robust convergence across different property classifications

### Example Interpretation

Position 745 detected by GS0, GS1, GS2:
- **GS0**: Exact amino acid convergence (VEG → EEE)
- **GS1**: Physicochemical group convergence (groups 1,3,2 → 3,3,3)
- **GS2**: Charge/size convergence (groups B,B,C → C,C,C)

This suggests the convergence is driven by both identity and broader physicochemical properties.

### Bootstrap Interpretation

Higher bootstrap support for certain schemes indicates:
- Property-level convergence is more consistent across trait permutations
- Functional constraints may be better captured by specific grouping schemes

## Technical Details

### Overlap Logic
- Uses Counter-based string overlap with character multiplicity
- "AAA/ASS" = 1 overlap (min of 3 A's and 1 A)
- "AAA/AAS" = 2 overlap (min of 3 A's and 2 A's)

### Pair-Aware Mode
- Conserved pairs tracked per grouping scheme
- Format: "overlap:pair_id1,pair_id2,..."
- Different schemes may show different conserved pair patterns

### P-value Calculation
- Hypergeometric test uses population-level group diversity
- All species in alignment recodified to group codes
- Frequency distribution calculated across entire population (~200 species)
- Different schemes produce different p-values based on group frequencies

## Files Modified

### PhyloPhere Nextflow Pipeline
- `nextflow.config`: Added `caap_mode` parameter
- `conf/modules.config`: Added `--caap_mode` flag to DISCOVERY and BOOTSTRAP processes
- `subworkflows/CT/ct_discovery.nf`: Passes args through (no change needed)
- `subworkflows/CT/ct_bootstrap.nf`: Passes args through (no change needed)

### CAAStools Core Modules
- `modules/caap_id.py`: New module implementing CAAP detection
- `modules/disco.py`: Modified to support CAAP mode
- `modules/boot.py`: Modified to support CAAP bootstrap
- `ct`: Added `--caap_mode` command-line flag

## Validation

GS0 results should exactly match classical CAAS results to validate correctness:
```bash
# Run classical CAAS
./ct discovery -a alignment.phy -t traits.tab -o classical.out

# Run CAAP mode
./ct discovery -a alignment.phy -t traits.tab -o caap.out --caap_mode

# Extract GS0 results and compare
grep "GS0" caap.out > caap_gs0.out
# These should have identical positions and p-values to classical.out
```

## References

- Original CAAS method: Muntané et al. (2018)
- Pair-aware implementation: Ramon et al. (in prep)
- CAAP extension: Ramon et al. (in prep)
