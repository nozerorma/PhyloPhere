#!/usr/bin/env python3
"""
Aggregator — adapted for CT_ACCUMULATION in PhyloPhere.

Changes vs. original to_integrate version:
  - read_species_list(): reads 3-col headerless traitfile (species, trait, pair);
    assigns contrast = 1 for every entry (single group accumulation mode).
  - read_metadata_caas(): reads filtered_discovery.tsv produced by CT_POSTPROC.
    Columns are the unique source of truth:
      Gene, Position, tag, caas, is_significant, Pvalue, pvalue_boot, Pattern,
      convergence_description, convergence_mode, CAAP_Group, amino_encoded,
      is_conserved_meta, conserved_pair, sig_hyp, sig_perm, sig_both,
      top_change_type, bottom_change_type, change_side, low_confidence_nodes,
      asr_is_conserved, comments, ..., Trait
    No fallback to legacy formats.
"""

import argparse
import os
import glob
import gc
import csv

import numpy as np
import logging
from Bio import AlignIO
from collections import defaultdict


# --------------------------
# Helpers
# --------------------------

def _as_bool(x):
    return str(x).strip().lower() in {'true', 't', '1', 'yes', 'y'}


def natural_sort_key(chromosome):
    if not chromosome:
        return (9999, '')
    s = str(chromosome).strip()
    s = s.replace('CHR', 'chr').replace('Chr', 'chr')
    s2 = s.replace('chr', '')
    specials = {'x': 23, 'y': 24, 'm': 25, 'mt': 25}
    key = s2.lower()
    leading = ''
    for ch in key:
        if ch.isdigit():
            leading += ch
        else:
            break
    if leading:
        return (int(leading), key[len(leading):])
    if key in specials:
        return (specials[key], '')
    return (1000, key)


# --------------------------
# Core I/O and Aggregation
# --------------------------

def read_species_list(species_file):
    """Read the CT traitfile (3-col, no header): species, trait, pair.

    All species are assigned to contrast group 1 (single-group accumulation mode).
    """
    logging.info(f"Reading species list from {species_file} (3-col, no header expected)")
    def _default_species_entry():
        return {'contrast': set(), 'trait': set(), 'pair': set(),
                'trait_by_contrast': {}, 'trait_by_pair': {}}
    species_data = defaultdict(_default_species_entry)
    with open(species_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                # Try whitespace split as fallback
                parts = line.split()
            if len(parts) >= 3:
                species = parts[0]
                try:
                    trait = int(parts[1])
                    pair  = int(parts[2])
                except ValueError:
                    logging.warning(f"Skipping malformed traitfile line: {line[:100]}")
                    continue
                contrast = 1  # fixed — single accumulation group
                entry = species_data[species]
                entry['contrast'].add(contrast)
                entry['trait'].add(trait)
                entry['pair'].add(pair)
    logging.info(f"Loaded {len(species_data)} species (all assigned to contrast group 1)")
    return species_data


def read_genomic_info(genomic_file):
    logging.info(f"Reading genomic info from {genomic_file}")
    genes = []
    with open(genomic_file) as f:
        headers = f.readline().strip().split('\t')
        gene_idx = headers.index('gene')
        chr_idx  = headers.index('chr')
        start_idx = headers.index('start')
        end_idx   = headers.index('end')
        msa_length_idx = headers.index('length')
        for line in f:
            parts = line.strip().split('\t')
            try:
                genes.append({
                    'gene': parts[gene_idx],
                    'chr':  parts[chr_idx],
                    'start': int(float(parts[start_idx])),
                    'end':   int(float(parts[end_idx])),
                    'msa_length': int(float(parts[msa_length_idx]))
                })
            except (ValueError, IndexError) as e:
                logging.warning(f"Skipping malformed line in genomic info: {line.strip()[:100]} — {e}")
                continue
    sorted_genes = sorted(genes, key=lambda x: (natural_sort_key(x['chr']), x['start']))
    logging.info(f"Processed {len(sorted_genes)} genes from genomic info")
    return sorted_genes


def read_bg_info(bg_file):
    logging.info(f"Reading CAAS background from {bg_file}")
    genes = []
    with open(bg_file) as f:
        for line in f:
            gene_name = line.strip()
            if gene_name:
                genes.append(gene_name)
    logging.info(f"Loaded {len(genes)} background genes")
    return genes


<<<<<<< HEAD
def read_metadata_caas(metadata_file):
=======
def read_metadata_caas(metadata_file, use_all_mrca_filter=False):
>>>>>>> asr
    """Read CAAS metadata from a filtered_discovery.tsv file.

    Source-of-truth columns (tab-separated):
      Gene, Position, tag, caas, is_significant, Pvalue, pvalue_boot, Pattern,
      convergence_description, convergence_mode, CAAP_Group, amino_encoded,
      is_conserved_meta, conserved_pair, sig_hyp, sig_perm, sig_both,
      top_change_type, bottom_change_type, change_side, low_confidence_nodes,
      asr_is_conserved, comments, ..., Trait

    Returns: dict[group][gene][msa_pos] = {tag, Pattern, amino_encoded, Pvalue,
                                           pvalue_boot, isSignificant}
    """
    logging.info(f"Reading metadata CAAS from {metadata_file if metadata_file else 'None'}")
    metadata = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

    if not metadata_file:
        return metadata

    with open(metadata_file) as f:
        raw_header = f.readline().strip()
        if not raw_header:
            logging.warning(f"Metadata CAAS file is empty: {metadata_file} — skipping")
            return metadata
        sep = '\t' if '\t' in raw_header else ','
        h = raw_header.split(sep)

        if 'Gene' not in h:
            logging.warning(f"Metadata CAAS file has no 'Gene' column (header: {raw_header[:120]}) — skipping")
            return metadata

        gene_idx          = h.index('Gene')
        pos_idx           = h.index('Position')
        tag_idx           = h.index('tag')
        pattern_idx       = h.index('Pattern')
        amino_idx         = h.index('amino_encoded') if 'amino_encoded' in h else None
        pvalue_idx        = h.index('Pvalue')
        sig_idx           = h.index('sig_both')
        pboot_idx         = h.index('pvalue_boot') if 'pvalue_boot' in h else None
        group_idx         = h.index('CAAP_Group')
        conserved_idx     = h.index('is_conserved_meta') if 'is_conserved_meta' in h else None
        asr_conserved_idx = h.index('asr_is_conserved') if 'asr_is_conserved' in h else None
<<<<<<< HEAD
=======
        asr_root_conserved_idx = h.index('asr_root_conserved') if 'asr_root_conserved' in h else None
>>>>>>> asr

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(sep)
            try:
                gene      = parts[gene_idx].strip()
                msa_pos   = int(float(parts[pos_idx].strip()))
                tag       = parts[tag_idx].strip()
                pattern   = parts[pattern_idx].strip()
                amino_conv = parts[amino_idx].strip() if amino_idx is not None else ''
                try:
                    caas_pval = float(parts[pvalue_idx])
                except Exception:
                    caas_pval = float('inf')
                is_sig = _as_bool(parts[sig_idx])
                pboot  = None
                if pboot_idx is not None:
                    try:
                        pboot = float(parts[pboot_idx])
                    except Exception:
                        pass
                group = parts[group_idx].strip()
                if not group or group in ('NA', 'na', 'N/A'):
                    group = '1'
            except (IndexError, ValueError) as e:
                logging.warning(f"Skipping malformed meta_caas line: {line[:120]} — {e}")
                continue

            # Exclude GS0 entries.
            if group.upper() == 'GS0':
                continue

            # Exclude dubious-conserved positions: is_conserved_meta=TRUE but asr_is_conserved=FALSE.
            if conserved_idx is not None and asr_conserved_idx is not None:
                _is_cons = _as_bool(parts[conserved_idx]) if conserved_idx < len(parts) else False
                _is_asr  = _as_bool(parts[asr_conserved_idx]) if asr_conserved_idx < len(parts) else True
                if _is_cons and not _is_asr:
                    continue
<<<<<<< HEAD
=======
                if use_all_mrca_filter and asr_root_conserved_idx is not None:
                    _is_root = _as_bool(parts[asr_root_conserved_idx]) if asr_root_conserved_idx < len(parts) else True
                    if _is_cons and not _is_root:
                        continue
>>>>>>> asr

            metadata[group][gene][msa_pos] = {
                'tag': tag,
                'Pattern': pattern,
                'AminoConv': amino_conv,
                'Pvalue': caas_pval,
                'Pvalue.boot': pboot,
                'isSignificant': is_sig,
            }

    logging.debug(f"Meta-CAAS loaded. Groups: {list(metadata.keys())}")
    return metadata


def get_combined_species(species_data, group_combo):
    groups = [int(g) for g in str(group_combo).split(',')]
    return {sp for sp in species_data if any(g in species_data[sp]['contrast'] for g in groups)}


# Alignment ops

def load_alignment(filepath, fmt):
    try:
        return AlignIO.read(filepath, fmt)
    except Exception as e:
        raise ValueError(f"Invalid alignment file: {str(e)}")


def validate_species(alignment, required_species):
    present_species = {rec.id for rec in alignment}
    missing = required_species - present_species
    if missing:
        raise ValueError(f"Missing species in alignment: {', '.join(sorted(missing))}")


def calculate_conservation(alignment, group_species=None):
    seq_len = alignment.get_alignment_length()
    general_counts = [defaultdict(int) for _ in range(seq_len)]
    for record in alignment:
        for pos in range(seq_len):
            aa = record.seq[pos]
            if aa != '-':
                general_counts[pos][aa] += 1
    results = {}
    for pos in range(seq_len):
        g_counts = general_counts[pos]
        g_total  = sum(g_counts.values()) or 1
        max_val  = max(g_counts.values(), default=0)
        cons_idx = (max_val / g_total) * 100
        results[pos] = {'cons_idx': cons_idx}
    return results


def calculate_masked_positions(alignment, target_species):
    target_set = set(target_species)
    seq_len = alignment.get_alignment_length()
    masked = set()
    for pos in range(seq_len):
        for rec in alignment:
            if rec.id in target_set and rec.seq[pos] == '-':
                masked.add(pos)
                break
    return masked


def aggregate(args):
    logging.info("Starting background aggregation for CT accumulation randomizations...")

    species_data   = read_species_list(args.species_list)
    gene_info_list = read_genomic_info(args.genomic_info)
    bg_list        = read_bg_info(args.bg_caas)
    metadata_dict  = read_metadata_caas(
        args.metadata_caas,
        use_all_mrca_filter=getattr(args, 'use_all_mrca_filter', False),
    ) if args.metadata_caas else {}

    # Filter gene_info_list to background genes only
    if bg_list:
        bg_set = set(bg_list)
        gene_info_list = [g for g in gene_info_list if g['gene'] in bg_set]
        logging.info(f"Filtered to {len(gene_info_list)} genes present in background list")

    # Assign global integer positions
    gene_offsets = {}
    accumulated_position = 0
    for gene_info in gene_info_list:
        gene_offsets[gene_info['gene']] = accumulated_position
        accumulated_position += gene_info['msa_length']
    logging.info(f"Global positions assigned. Total positions: {accumulated_position}")

    # All species treated as one accumulation group (single-group mode)
    all_species = set(species_data.keys())
    logging.info(f"Total species for accumulation: {len(all_species)}")

    # Output: single enriched global CSV
    # Columns: gene, position, chr, start, end, msa_pos, cons_idx, masked, iscaas
    aggregated_filename = f"{args.output_prefix}_global.csv"
    aggregated_fieldnames = ['gene', 'position', 'chr', 'start', 'end', 'msa_pos', 'cons_idx', 'masked', 'iscaas']
    aggregated_file = open(aggregated_filename, 'w', newline='')
    aggregated_writer = csv.DictWriter(aggregated_file, fieldnames=aggregated_fieldnames)
    aggregated_writer.writeheader()

    overall_caas_cons = []

    # Use '*' (no extension filter) to match flat alignment directory — same as CT discovery
    # Split on the first '.' to extract the gene symbol, handling multi-dot names like
    # GENE.Species.filter2.phy produced by the CT discovery step.
    alignment_files = {
        os.path.basename(f).split('.')[0]: f
        for f in glob.glob(os.path.join(args.alignment_dir, '*'))
        if os.path.isfile(f)
    }
    logging.info(f"Found {len(alignment_files)} alignment files in {args.alignment_dir!r}.")
    if not alignment_files:
        logging.error(
            f"No alignment files found in {args.alignment_dir!r}. "
            "Check --alignment-dir path and --alignment-format."
        )

    genes_written = 0
    for gene_info in gene_info_list:
        gene_name = gene_info['gene']
        if gene_name not in alignment_files:
            logging.warning(f"No alignment file for gene {gene_name}, skipping")
            continue
        try:
            alignment = AlignIO.read(alignment_files[gene_name], args.alignment_format)
            seq_len   = alignment.get_alignment_length()
            if seq_len != gene_info['msa_length']:
                logging.warning(
                    f"MSA length mismatch for {gene_name}: observed {seq_len} "
                    f"vs metadata {gene_info['msa_length']}"
                )
            general_cons = calculate_conservation(alignment)
            masked_pos   = calculate_masked_positions(alignment, all_species)

            for msa_pos in range(seq_len):
                global_pos = gene_offsets[gene_name] + msa_pos
                is_caas    = any(
                    msa_pos in metadata_dict.get(g, {}).get(gene_name, {})
                    for g in metadata_dict
                ) if metadata_dict else False
                aggregated_writer.writerow({
                    'gene':     gene_name,
                    'msa_pos':  msa_pos,
                    'position': global_pos,
                    'chr':      gene_info['chr'],
                    'start':    gene_info['start'],
                    'end':      gene_info['end'],
                    'cons_idx': general_cons[msa_pos]['cons_idx'],
                    'masked':   (msa_pos in masked_pos),
                    'iscaas':   is_caas,
                })
                if is_caas:
                    overall_caas_cons.append(general_cons[msa_pos]['cons_idx'])

            genes_written += 1
            del alignment
            gc.collect()

        except Exception as e:
            logging.error(f"Error processing gene {gene_name}: {str(e)}")
            continue

    aggregated_file.close()

    if genes_written == 0:
        logging.error(
            f"No genes were written to {aggregated_filename}. "
            "Verify --alignment-dir contains files whose basenames before the first '.' "
            "match gene names from --genomic-info."
        )
    else:
        logging.info(f"Aggregation wrote {genes_written} genes to {aggregated_filename}.")

    # Deciles
    decile_file = f"{args.output_prefix}_deciles.csv"
    if overall_caas_cons:
        decile_thresholds = np.percentile(overall_caas_cons, list(range(10, 100, 10)))
        logging.info(f"Computed decile thresholds: {decile_thresholds}")
    else:
        decile_thresholds = []
        logging.warning("No CAAS cons_idx values recorded; deciles cannot be computed.")
    with open(decile_file, 'w', newline='') as df:
        writer = csv.writer(df)
        writer.writerow(['Decile', 'Threshold'])
        for i, threshold in enumerate(decile_thresholds, start=1):
            writer.writerow([i, threshold])
    logging.info(f"Decile thresholds written to: {decile_file}")
    logging.info("Aggregation complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CT_ACCUMULATION: background aggregator')
    parser.add_argument('-a', '--alignment-dir',   required=True, help='Directory containing alignment files')
    parser.add_argument('-f', '--alignment-format', default='phylip-relaxed')
    parser.add_argument('-i', '--genomic-info',    required=True, help='Gene genomic info TSV (gene, chr, start, end, length)')
    parser.add_argument('-s', '--species-list',    required=True, help='Species traitfile (3-col, no header: species trait pair)')
    parser.add_argument('-m', '--metadata-caas',   help='Meta-CAAS file (original or global_meta_caas.tsv format)')
    parser.add_argument('-b', '--bg-caas',         help='Cleaned background gene list (one gene per line, no header)')
    parser.add_argument('-o', '--output-prefix',   required=True, help='Prefix for output files')
    parser.add_argument('--use-all-mrca-filter', action='store_true',
                        help='Exclude conserved-meta rows failing asr_root_conserved')
    parser.add_argument('--log-level', default='INFO')

    args = parser.parse_args()
    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s - %(message)s')
    aggregate(args)
