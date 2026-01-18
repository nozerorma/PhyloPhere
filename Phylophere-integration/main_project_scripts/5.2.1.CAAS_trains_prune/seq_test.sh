#!/bin/bash

# Usage: ./seq_test.sh input_file.caas [output_dir]

set -euo pipefail

INPUT="${1:-}"
OUTDIR="${2:-}"
SCRIPT="filter_caas_clusters-param.py"
SUMMARY_FILE="discarded_summary.tsv"

# ─── Start Total Runtime Timer ────────────────────────────────────────────────
SECONDS=0
trap 'echo -e "\nTotal runtime for all operations: $(date -u -d "@$SECONDS" +%H:%M:%S)"' EXIT

# ─── Validate Input ───────────────────────────────────────────────────────────
if [[ -z "$INPUT" ]]; then
  echo "Usage: $0 input_file.caas [output_dir]"
  exit 1
fi

if [[ ! -f "$INPUT" ]]; then
  echo "Error: Input file '$INPUT' not found."
  exit 1
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "Error: Python script '$SCRIPT' not found in current directory."
  exit 1
fi

# ─── Parameter Grid ───────────────────────────────────────────────────────────
MINLEN_VALUES=(3 4 10)
MAXCAAS_VALUES=(0.5 0.6 0.7 0.8)
# MINLEN_VALUES=(3)
# MAXCAAS_VALUES=(0.7)

# ─── Runtime Tracking ─────────────────────────────────────────────────────────
declare -A runtimes

# ─── Run Filtering ────────────────────────────────────────────────────────────
for minlen in "${MINLEN_VALUES[@]}"; do
  for maxcaas in "${MAXCAAS_VALUES[@]}"; do
    echo "▶ Running filter: minlen=$minlen, maxcaas=$maxcaas"
    
    t0=$(date +%s)
    python3 "$SCRIPT" -i "$INPUT" -l "$minlen" -c "$maxcaas" --verbose
    t1=$(date +%s)
    
    elapsed=$((t1 - t0))
      hrs=$((elapsed / 3600))
      mins=$(((elapsed % 3600) / 60))
      secs=$((elapsed % 60))
      printf "    ⏱ Runtime: %02d:%02d:%02d (hh:mm:ss)\n" "$hrs" "$mins" "$secs"

      # Capture output filename
      maxcaas_str=$(echo "$maxcaas * 100 / 1" | bc)
      outfile=$(find . -maxdepth 1 -type f -name "*filtered.minlen${minlen}.maxcaas${maxcaas_str}.tsv" -print -quit)
      if [[ -z "$outfile" ]]; then
        echo "⚠️  Warning: Could not locate output file for minlen=$minlen maxcaas=$maxcaas (expected maxcaas${maxcaas_str})"
        continue
      fi
      outfile=$(basename "$outfile")
      runtimes["$outfile"]=$(printf "%02d:%02d:%02d" "$hrs" "$mins" "$secs")
  done
done

# ─── Move Output (Optional) ───────────────────────────────────────────────────
if [[ -n "$OUTDIR" ]]; then
  mkdir -p "$OUTDIR/logs"
  mv -- *.tsv "$OUTDIR"
  echo "✔ Output files moved to: $OUTDIR"
  mv -- *.log "$OUTDIR/logs"
  echo "✔ Log files moved to: $OUTDIR"
else
  OUTDIR="."  # fallback to current directory
  echo "ℹ No output directory specified. Output files remain in: $(pwd)"
fi

# ─── Count Unique Discarded Positions and Save ────────────────────────────────
echo -e "\n🧮 Discarded Position Counts (unique Gene@Position per file):"
echo -e "File\tDiscardedCount\tRuntime" > "$OUTDIR/$SUMMARY_FILE"

for file in "$OUTDIR"/*.filtered.minlen*; do
  if [[ -f "$file" ]]; then
    echo "Processing file: $file"
    basefile=$(basename "$file")
    runtime="${runtimes[$basefile]:-00:00:00}"

    # Count unique discarded positions
    count=$(awk -F'\t' '
      function trim(s) { sub(/^[ \t\r\n]+/, "", s); sub(/[ \t\r\n]+$/, "", s); return s }

      BEGIN { gene_col=0; pos_col=0; flag_col=0 }

      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if (trim($i) == "Gene") gene_col = i;
          if (trim($i) == "Position") pos_col = i;
          if (trim($i) == "ClusteringFlag") flag_col = i;
        }
        if (flag_col == 0) {
          print "ERROR: Missing ClusteringFlag column" > "/dev/stderr";
          exit 1;
        }
        next
      }

      trim($flag_col) == "Discarded" {
        n++
      }

      END { print n+0 }
    ' "$file")

    echo -e "$basefile\t$count\t$runtime"
    echo -e "$basefile\t$count\t$runtime" >> "$OUTDIR/$SUMMARY_FILE"
  fi
done

echo -e "\n📄 Discarded counts & runtimes saved to: $OUTDIR/$SUMMARY_FILE"
echo "✅ All operations completed successfully."
echo "You can now check the output, logs, and runtime summary in: $OUTDIR"
echo "Thank you for using the CAAS train filtering script! 🚂✨"
