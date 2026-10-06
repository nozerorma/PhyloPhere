#!/usr/bin/env python3
# concatenate.py — Build the global position table (one row per alignment column) for the accumulation null.
# PhyloPhere | subworkflows/CT_ACCUMULATION/local/src/aggregation/

"""
Aggregator: writes <prefix>_global.csv, the table the randomization phase samples from.

Every alignment column of every background gene gets one row with a global integer
position, a conservation value, a mask flag (column gapped in any species of the
traitfile) and an iscaas flag (the column is a CAAS of the filtered discovery table in
any group). Only the keys (group, gene, msa_pos) of the discovery table are used, never
its per-position values.

Called by:  CT_ACCUMULATION Nextflow process (ctacc_run.nf CT_ACCUMULATION_AGGREGATE → main.py --tool aggregate)
Inputs:     alignment directory, genomic-info TSV (gene, chr, start, end, length), traitfile
            (file or directory), filtered_discovery.tsv from CT_POSTPROC, cleaned background
            gene list, optional directory of Valdar <gene>.entropy.tsv files
Outputs:    <output-prefix>_global.csv with columns gene, position, chr, start, end, msa_pos,
            cons_idx, masked, iscaas
"""

# ── Standard library ──────────────────────────────────────────────────────────
import argparse
import os
import glob
import gc
import csv
from pathlib import Path

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import logging
from Bio import AlignIO
from collections import defaultdict


# ── Helpers ───────────────────────────────────────────────────────────────────


def natural_sort_key(chromosome):
    """Sort key for chromosome names: numeric first, then X, Y, M/MT, then anything else.

    Accepts names with or without a 'chr' prefix; empty names sort last.
    """
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


# ── Core I/O and aggregation ──────────────────────────────────────────────────


def read_species_list(species_file):
    """Read the CT traitfile or directory of traitfiles (3-col, no header): species, trait, pair.

    Returns a dict keyed by species. Only the keys are used downstream (the species whose
    gaps mask a column); every species is assigned to contrast group 1. A directory is read as
    its traitfile_H*.tab files, else its *.tab files except traitfile_fop.tab, else all of its files.
    """
    logging.info(f"Reading species list from {species_file} (3-col, no header expected)")
    def _default_species_entry():
        return {'contrast': set(), 'trait': set(), 'pair': set(),
                'trait_by_contrast': {}, 'trait_by_pair': {}}
    species_data = defaultdict(_default_species_entry)

    species_path = Path(species_file)
    files_to_read = []
    if species_path.is_dir():
        h_files = sorted(species_path.glob("traitfile_H*.tab"))
        if h_files:
            files_to_read = h_files
        else:
            files_to_read = [
                f for f in sorted(species_path.glob("*.tab"))
                if f.name != "traitfile_fop.tab"
            ]
        if not files_to_read:
            files_to_read = [f for f in sorted(species_path.glob("*")) if f.is_file()]
    elif species_path.is_file():
        files_to_read = [species_path]

    for fpath in files_to_read:
        with open(fpath, "r", encoding="utf-8-sig") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 3:
                    parts = line.split()
                if len(parts) >= 3:
                    species = parts[0]
                    try:
                        trait = int(parts[1])
                        pair  = int(parts[2])
                    except ValueError:
                        logging.warning(f"Skipping malformed traitfile line in {fpath}: {line[:100]}")
                        continue
                    contrast = 1  # one accumulation group: every species counts
                    entry = species_data[species]
                    entry['contrast'].add(contrast)
                    entry['trait'].add(trait)
                    entry['pair'].add(pair)
    logging.info(f"Loaded {len(species_data)} species across {len(files_to_read)} traitfile(s) (all assigned to contrast group 1)")
    return species_data


def read_genomic_info(genomic_file):
    """Read the genomic-info TSV (columns gene, chr, start, end, length; located by header name).

    Returns the genes with coordinates as dicts, sorted by chromosome (natural order) and start.
    Genes without start or end are dropped: they have no genomic position to order them by.
    """
    logging.info(f"Reading genomic info from {genomic_file}")
    genes = []
    with open(genomic_file) as f:
        headers = [c.strip().lower() for c in f.readline().strip().split('\t')]
        gene_idx = headers.index('gene')
        chr_idx  = headers.index('chr')
        start_idx = headers.index('start')
        end_idx   = headers.index('end')
        msa_length_idx = headers.index('length')
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            parts = line_str.split('\t')
            if len(parts) <= max(gene_idx, chr_idx, start_idx, end_idx, msa_length_idx):
                continue
            start_str = parts[start_idx].strip()
            end_str = parts[end_idx].strip()
            if not start_str or not end_str:
                logging.debug(f"Skipping unlocalized gene in genomic info (no coordinates): {parts[gene_idx]}")
                continue
            try:
                genes.append({
                    'gene': parts[gene_idx],
                    'chr':  parts[chr_idx],
                    'start': int(float(start_str)),
                    'end':   int(float(end_str)),
                    'msa_length': int(float(parts[msa_length_idx]))
                })
            except (ValueError, IndexError) as e:
                logging.warning(f"Skipping malformed line in genomic info: {line_str[:100]} — {e}")
                continue
    sorted_genes = sorted(genes, key=lambda x: (natural_sort_key(x['chr']), x['start']))
    logging.info(f"Processed {len(sorted_genes)} genes from genomic info")
    return sorted_genes


def read_bg_info(bg_file):
    """Read the background gene list (one gene per line, no header)."""
    logging.info(f"Reading CAAS background from {bg_file}")
    genes = []
    with open(bg_file) as f:
        for line in f:
            gene_name = line.strip()
            if gene_name:
                genes.append(gene_name)
    logging.info(f"Loaded {len(genes)} background genes")
    return genes


def read_metadata_caas(metadata_file):
    """Read the CAAS positions of a filtered_discovery.tsv (tab or comma separated).

    Columns are matched case-insensitively: Gene, Position (or msa_pos), caap_group (or caap,
    group), and optionally tag (or tag_support), convergence_type and amino_encoded. A row with
    no group value goes to group '1'. A file without a Gene or Position column is skipped with
    a warning and yields an empty table.

    Returns: dict[group][gene][msa_pos] = {tag, convergence_type, caas}, where caas holds the
    amino_encoded value. Only the (group, gene, msa_pos) keys are consumed downstream.
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
        h = [c.strip() for c in raw_header.split(sep)]

        col_map = {col.lower(): idx for idx, col in enumerate(h)}

        gene_idx = col_map.get('gene')
        if gene_idx is None:
            logging.warning(f"Metadata CAAS file has no 'Gene' column (header: {raw_header[:120]}) — skipping")
            return metadata

        pos_idx = col_map.get('position')
        if pos_idx is None:
            pos_idx = col_map.get('msa_pos')
        if pos_idx is None:
            logging.warning(f"Metadata CAAS file has no 'Position' column (header: {raw_header[:120]}) — skipping")
            return metadata

        tag_idx = col_map.get('tag')
        if tag_idx is None:
            tag_idx = col_map.get('tag_support')

        convergence_idx = col_map.get('convergence_type')
        amino_idx = col_map.get('amino_encoded')
        group_idx = col_map.get('caap_group')
        if group_idx is None:
            group_idx = col_map.get('caap')
        if group_idx is None:
            group_idx = col_map.get('group')

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(sep)
            try:
                gene = parts[gene_idx].strip()
                msa_pos = int(float(parts[pos_idx].strip()))
                tag = parts[tag_idx].strip() if tag_idx is not None and tag_idx < len(parts) else ''
                convergence = parts[convergence_idx].strip() if convergence_idx is not None and convergence_idx < len(parts) else ''
                amino_conv = parts[amino_idx].strip() if amino_idx is not None and amino_idx < len(parts) else ''
                group = parts[group_idx].strip() if group_idx is not None and group_idx < len(parts) else '1'
                if not group or group in ('NA', 'na', 'N/A'):
                    group = '1'
            except (IndexError, ValueError) as e:
                logging.warning(f"Skipping malformed meta_caas line: {line[:120]} — {e}")
                continue

            metadata[group][gene][msa_pos] = {
                'tag': tag,
                'convergence_type': convergence,
                'caas': amino_conv,
            }

    logging.debug(f"Meta-CAAS loaded. Groups: {list(metadata.keys())}")
    return metadata


# ── Alignment operations ──────────────────────────────────────────────────────


def calculate_conservation(alignment, group_species=None):
    """Per-column conservation: percentage of the non-gap residues that are the majority residue.

    Returns {column index (0-based): {'cons_idx': value}}; a column of gaps only, or an empty
    alignment, gets 0.0. group_species is accepted for call compatibility and is not used.
    """
    seq_len = alignment.get_alignment_length()
    if len(alignment) == 0:
        return {pos: {'cons_idx': 0.0} for pos in range(seq_len)}

    # One row per sequence, one uint8 ASCII code per column
    seqs = np.array([list(str(rec.seq).encode('ascii')) for rec in alignment], dtype=np.uint8)

    # Gap is '-' (ASCII 45)
    gap_char = 45
    results = {}
    for pos in range(seq_len):
        col = seqs[:, pos]
        non_gap = col[col != gap_char]
        g_total = len(non_gap)
        if g_total == 0:
            cons_idx = 0.0
        else:
            counts = np.bincount(non_gap)
            max_val = counts.max() if len(counts) > 0 else 0
            cons_idx = (max_val / g_total) * 100
        results[pos] = {'cons_idx': cons_idx}
    return results


def calculate_masked_positions(alignment, target_species):
    """Set of 0-based columns with a gap in at least one of the target species."""
    target_set = set(target_species)
    seq_len = alignment.get_alignment_length()
    target_recs = [rec for rec in alignment if rec.id in target_set]
    if not target_recs:
        return set()

    seqs = np.array([list(str(rec.seq).encode('ascii')) for rec in target_recs], dtype=np.uint8)
    masked_cols = np.any(seqs == 45, axis=0)
    return set(np.where(masked_cols)[0])


# ── Aggregation ───────────────────────────────────────────────────────────────


def aggregate(args):
    """Write <output_prefix>_global.csv for the genes of the genomic-info table.

    Global positions are assigned to every background gene in genomic order before the
    alignments are read, so a gene without an alignment file leaves a gap in the numbering.
    Reads from args: alignment_dir, alignment_format, genomic_info, species_list, bg_caas,
    metadata_caas, entropy_dir (optional), output_prefix.
    """
    logging.info("Starting background aggregation for CT accumulation randomizations...")

    species_data   = read_species_list(args.species_list)
    gene_info_list = read_genomic_info(args.genomic_info)
    bg_list        = read_bg_info(args.bg_caas)
    metadata_dict  = read_metadata_caas(args.metadata_caas) if args.metadata_caas else {}

    # Map gene names to their <gene>.entropy.tsv (the per-clade tables are not used)
    entropy_files = {}
    if getattr(args, 'entropy_dir', None) and os.path.isdir(args.entropy_dir):
        logging.info(f"Scanning entropy directory: {args.entropy_dir}")
        for f in os.listdir(args.entropy_dir):
            if f.endswith('.entropy.tsv') and not f.endswith('.clade_entropy.tsv'):
                gene = f.split('.')[0]
                entropy_files[gene] = os.path.join(args.entropy_dir, f)
        logging.info(f"Found {len(entropy_files)} matching entropy files")

    # Keep the background genes only (an empty list keeps all genes)
    if bg_list:
        bg_set = set(bg_list)
        gene_info_list = [g for g in gene_info_list if g['gene'] in bg_set]
        logging.info(f"Filtered to {len(gene_info_list)} genes present in background list")

    # Global positions: each gene starts where the previous one ends
    gene_offsets = {}
    accumulated_position = 0
    for gene_info in gene_info_list:
        gene_offsets[gene_info['gene']] = accumulated_position
        accumulated_position += gene_info['msa_length']
    logging.info(f"Global positions assigned. Total positions: {accumulated_position}")

    # Every species of the traitfile belongs to the one accumulation group
    all_species = set(species_data.keys())
    logging.info(f"Total species for accumulation: {len(all_species)}")

    # Output columns: gene, position (global), chr, start, end, msa_pos (0-based column in the gene), cons_idx, masked, iscaas
    aggregated_filename = f"{args.output_prefix}_global.csv"
    aggregated_fieldnames = ['gene', 'position', 'chr', 'start', 'end', 'msa_pos', 'cons_idx', 'masked', 'iscaas']
    aggregated_file = open(aggregated_filename, 'w', newline='')
    aggregated_writer = csv.DictWriter(aggregated_file, fieldnames=aggregated_fieldnames)
    aggregated_writer.writeheader()

    # No extension filter: the alignment directory is flat. The gene symbol is the text before
    # the first '.', so multi-dot names such as GENE.Species.filter2.phy map to GENE.
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

            # Valdar variability per column (the file's position is 1-based)
            entropy_vals = {}
            entropy_loaded = False
            if getattr(args, 'entropy_dir', None):
                if gene_name in entropy_files:
                    try:
                        with open(entropy_files[gene_name], 'r') as ef:
                            reader = csv.DictReader(ef, delimiter='\t')
                            for row in reader:
                                pos = int(row['position'])
                                var = float(row['variability'])
                                entropy_vals[pos - 1] = var
                        entropy_loaded = True
                    except Exception as e:
                        logging.warning(f"Error reading entropy file for {gene_name}: {e}")
                else:
                    logging.warning(f"No entropy file found for gene {gene_name}")

            general_cons = calculate_conservation(alignment)
            masked_pos   = calculate_masked_positions(alignment, all_species)

            for msa_pos in range(seq_len):
                global_pos = gene_offsets[gene_name] + msa_pos
                is_caas    = any(
                    msa_pos in metadata_dict.get(g, {}).get(gene_name, {})
                    for g in metadata_dict
                ) if metadata_dict else False

                # Conservation value and mask
                is_masked = (msa_pos in masked_pos)
                if getattr(args, 'entropy_dir', None):
                    # Variability mode: a column without a value gets an empty cons_idx
                    if entropy_loaded and msa_pos in entropy_vals:
                        cons_val = entropy_vals[msa_pos]
                    else:
                        cons_val = "" # empty value
                        is_masked = True # a column without a value is masked so the null excludes it
                else:
                    cons_val = general_cons[msa_pos]['cons_idx']

                aggregated_writer.writerow({
                    'gene':     gene_name,
                    'msa_pos':  msa_pos,
                    'position': global_pos,
                    'chr':      gene_info['chr'],
                    'start':    gene_info['start'],
                    'end':      gene_info['end'],
                    'cons_idx': cons_val,
                    'masked':   is_masked,
                    'iscaas':   is_caas,
                })

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
    parser.add_argument('--entropy-dir',           help='Directory containing Valdar entropy (.entropy.tsv) files')
    parser.add_argument('--log-level', default='INFO')

    args = parser.parse_args()
    numeric_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s - %(message)s')
    aggregate(args)
