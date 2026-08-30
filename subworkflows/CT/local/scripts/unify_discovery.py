#!/usr/bin/env python3
"""
Unify multi-hypothesis CAAStools discovery outputs.
Deduplicates candidate hits by (gene, position, caap_group), keeping 1 representative
row per unique candidate site while tracking n_hypotheses and supporting_hypotheses.
"""

import sys
import re
from collections import defaultdict

def extract_h_num(trait_str):
    m = re.search(r'traitfile_H(\d+)\.tab', str(trait_str))
    if m:
        return int(m.group(1))
    return float('inf')

def unify_discovery(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        with open(output_file, 'w') as out:
            out.write("gene\tmode\tcaap_group\ttrait\tposition\tcaas\tpattern\tffgn\tfbgn\tgfg\tgbg\tmfg\tmbg\tffg\tfbg\tms\tn_hypotheses\tsupporting_hypotheses\n")
        return

    header_cols = lines[0].split('\t')
    
    # Identify column indices
    col_idx = {col: i for i, col in enumerate(header_cols)}
    
    gene_i = col_idx.get('gene', 0)
    pos_i = col_idx.get('position', 4)
    caap_i = col_idx.get('caap_group', 2)
    trait_i = col_idx.get('trait', 3)
    
    # Group rows by (gene, position, caap_group)
    grouped_rows = defaultdict(list)
    for line in lines[1:]:
        cols = line.split('\t')
        if len(cols) <= max(gene_i, pos_i, caap_i):
            continue
        key = (cols[gene_i], cols[pos_i], cols[caap_i])
        grouped_rows[key].append(cols)

    # Build new header
    new_header = list(header_cols)
    if 'n_hypotheses' not in new_header:
        ms_pos = col_idx.get('ms', len(new_header) - 1) + 1 if 'ms' in col_idx else len(new_header)
        new_header.insert(ms_pos, 'n_hypotheses')
        new_header.insert(ms_pos + 1, 'supporting_hypotheses')
    
    out_lines = ['\t'.join(new_header)]

    for key, row_list in grouped_rows.items():
        # Sort rows by hypothesis number (prefer H1 if present)
        row_list.sort(key=lambda r: extract_h_num(r[trait_i]))
        
        # Collect all supporting hypotheses
        h_nums = []
        traits_seen = []
        for r in row_list:
            t = r[trait_i]
            if t not in traits_seen:
                traits_seen.append(t)
            h_num = extract_h_num(t)
            if h_num != float('inf'):
                h_str = f"H{h_num}"
                if h_str not in h_nums:
                    h_nums.append(h_str)
            else:
                if t not in h_nums:
                    h_nums.append(t)
        
        # Sort H numbers numerically
        h_nums_sorted = sorted(h_nums, key=lambda x: int(x[1:]) if x.startswith('H') and x[1:].isdigit() else 999999)
        
        n_hyp = str(len(traits_seen))
        supp_hyp = ",".join(h_nums_sorted)
        
        # Take representative row (first row after sorting by H number)
        rep_row = list(row_list[0])
        
        # Insert n_hypotheses and supporting_hypotheses
        ms_pos_row = col_idx.get('ms', len(rep_row) - 1) + 1 if 'ms' in col_idx else len(rep_row)
        rep_row.insert(ms_pos_row, n_hyp)
        rep_row.insert(ms_pos_row + 1, supp_hyp)
        
        out_lines.append('\t'.join(rep_row))

    with open(output_file, 'w') as out:
        out.write('\n'.join(out_lines) + '\n')
    print(f"Unified {len(lines)-1} raw rows into {len(grouped_rows)} unique candidate sites in {output_file}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: unify_discovery.py <input_discovery.tab> <output_unified_discovery.tab>")
        sys.exit(1)
    unify_discovery(sys.argv[1], sys.argv[2])
