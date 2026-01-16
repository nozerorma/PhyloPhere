#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ---------------------------------------
# Part 1: Trait → Hypothesis → Analysis
# ---------------------------------------
for trait_dir in */; do
  trait=${trait_dir%/}
  for hyp_path in "$trait"/*/; do
    hyp=${hyp_path%/}
    hyp_name=$(basename "$hyp")

    for ana_path in "$hyp"/*/; do
      ana=${ana_path%/}
      ana_name=$(basename "$ana")

      if [[ "$ana_name" != "discovery" && "$ana_name" != "bootstrap" ]]; then
        continue
      fi
      
      # Output directory and file for this trait/analysis/hypothesis
      dir_out="$trait/$ana_name"
      mkdir -p "$dir_out"
      out_file="$dir_out/${hyp_name}.caas"
      > "$out_file"  # truncate before writing

      header_saved=false
      # Gather all .output files under this analysis
      find "$ana" -type f -name "*.output" | sort | while read -r file; do
        if ! $header_saved; then
          cat "$file" > "$out_file"
          header_saved=true
        else
          tail -n +2 "$file" >> "$out_file"
        fi
      done
    done
  done
done

# -------------------------------------------------
# Part 2: Merge all Hypothesis.caas into Analysis.caas
# -------------------------------------------------
for trait_dir in */; do
  trait=${trait_dir%/}

  # For each analysis folder under this trait
  for ana_dir in "$trait"/*/; do
    ana=${ana_dir%/}
    ana_name=$(basename "$ana")

    if [[ "$ana_name" != "discovery" && "$ana_name" != "bootstrap" ]]; then
      continue
    fi

    temp_file="$trait/${ana_name}.temp"
    final_file="$trait/${ana_name}.caas"

    > "$temp_file"
    header_saved=false

    # Concatenate each Hypothesis.caas in this analysis folder
    find "$ana" -maxdepth 1 -type f -name "*.caas" | sort | while read -r f; do
      if ! $header_saved; then
        cat "$f" > "$temp_file"
        header_saved=true
      else
        tail -n +2 "$f" >> "$temp_file"
      fi
    done

    # Remove duplicated header lines (assumes header starts with "Gene")
    awk '
      BEGIN { seen=0 }
      {
        if ($1 == "Gene") {
          if (seen == 0) {
            print; seen=1
          }
        } else {
          print
        }
      }
    ' "$temp_file" > "$final_file"

    rm "$temp_file"
  done
done

