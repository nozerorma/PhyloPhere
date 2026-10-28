#!/usr/bin/env python3
# position_sources.py: readers of position sources projected to CAAS coordinates (module, not a program).
# PhyloPhere | subworkflows/POSENRICH/local/src/

"""Read the position sources of POSENRICH and project them to CAAS columns.

Positions are zero-based columns of the CAAS trimmed alignment. MAP prot_ali_col and ori_codon_col are
one-based; hg38_aa_pos is the explicit reference protein residue. Pfam requires an explicit coordinate_space:
CAAS_0based with alignment_positions, or reference_protein_1based with verified endpoints. UCR requires
CAAS_0based, trimmed_alignment_1based or verified raw_alignment_1based coordinates. FUBAR hg38_aa_pos indexes
the reference protein; FADE and custom positions already use CAAS space. Plain text and gzip inputs are accepted.
Records that cannot be projected go to the `pending` list when one is given, never into a collection.
"""

import csv
import glob
import gzip
import os
import re
import sys

import pandas as pd


# ── Input handling ──────────────────────────────────────────────────────────────────────


def resolve_path(path):
    """Return an existing path for `path`, accepting a `.gz` sibling.

    Returns None when neither `path` nor `<path>.gz` exists. os.path.exists() follows
    symlinks, so a dangling symlink (a Nextflow stage-in whose target is not reachable)
    resolves to None instead of failing later on open.
    """
    if path is None:
        return None
    if os.path.exists(path):
        return path
    if not path.endswith(".gz") and os.path.exists(path + ".gz"):
        return path + ".gz"
    return None


def open_maybe_gz(path, mode="rt"):
    """Open a plain or gzip-compressed text file transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)



def load_universe_genes(cleaned_background_path):
    """Set of genes of the universe file (one per line), or None when absent."""
    if not cleaned_background_path or not os.path.exists(cleaned_background_path):
        return None
    genes = set()
    with open(cleaned_background_path, 'r') as f:
        for line in f:
            g = line.strip()
            if g:
                genes.add(g)
    print(f"Loaded {len(genes)} active genes from universe filter.")
    return genes


def parse_map_file(path):
    """Read explicit named coordinates; never infer residue numbers by counting rows."""
    selected_cols, col_to_genomic, residue_to_col, raw_to_col = [], {}, {}, {}
    coords, seen_cols = [], set()
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'status', 'prot_ali_col', 'ori_codon_col', 'hg38_aa_pos', 'hg38_nt_coord'}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"MAP missing columns: {sorted(required - set(reader.fieldnames or []))}")
        for row in reader:
            if row['status'] != 'selected':
                continue
            col = int(row['prot_ali_col']) - 1
            if col < 0 or col in seen_cols:
                raise ValueError('MAP has invalid or repeated trimmed columns')
            selected_cols.append(col)
            seen_cols.add(col)
            raw_to_col[int(row['ori_codon_col'])] = col
            residue = row['hg38_aa_pos']
            if residue not in ('NA', '', '-'):
                residue = int(residue)
                if residue in residue_to_col:
                    raise ValueError('MAP has repeated reference residues')
                residue_to_col[residue] = col
            genomic = row['hg38_nt_coord']
            if re.fullmatch(r'chr[0-9XYM]+:\d+', genomic):
                col_to_genomic[col] = genomic
                coords.append(int(genomic.split(':')[1]))
    strand = '-' if sum((b > a) - (b < a) for a, b in zip(coords, coords[1:])) < 0 else '+'
    return selected_cols, col_to_genomic, residue_to_col, strand, raw_to_col


def build_map_cache(map_dir, universe_genes):
    """Parse the MAP file of every gene of the universe into {gene: parsed fields}.

    The gene is the file name up to the first dot. All genes are kept when there is no
    universe. Files that fail to parse are skipped with a warning.
    """
    map_cache = {}
    if not map_dir or not os.path.exists(map_dir):
        print("No valid map_dir provided; skipping MAP coordinate caching.")
        return map_cache
    map_files = glob.glob(os.path.join(map_dir, "*.map.tsv"))
    print(f"Scanning {len(map_files)} MAP files...")
    seen = set()
    for path in map_files:
        gene = os.path.basename(path).split('.')[0]
        if gene in seen:
            raise ValueError(f"Several MAP files for gene {gene}; choose one reference")
        seen.add(gene)
        filename = os.path.basename(path)
        gene = filename.split('.')[0]
        if universe_genes is not None and gene not in universe_genes:
            continue
        try:
            selected_cols, col_to_genomic, residue_to_col, strand, raw_to_col = parse_map_file(path)
            map_cache[gene] = {
                'selected_cols': set(selected_cols),
                'col_to_genomic': col_to_genomic,
                'residue_to_col': residue_to_col,
                'strand': strand,
                'raw_to_col': raw_to_col,
                'map_file': filename
            }
        except Exception as e:
            print(f"Warning: failed to parse MAP file for {gene}: {e}", file=sys.stderr)
    print(f"Cached coordinates mapping for {len(map_cache)} genes.")
    return map_cache


# ── Writing ───────────────────────────────────────────────────────────────────



# ── Source readers ───────────────────────────────────────────────────────────────────


def column_runs(columns):
    """Runs of consecutive CAAS columns of a sorted list of ints, as [(start, end), ...] (inclusive)."""
    runs = []
    for c in columns:
        if runs and c == runs[-1][1] + 1:
            runs[-1][1] = c
        else:
            runs.append([c, c])
    return [(a, b) for a, b in runs]


def window_term(genomic, size=1_000_000):
    """(term name, description) of the genomic window of a 'chrN:coordinate' string, or None.

    Windows are chromosome-separated, `size` bp wide, with one-based inclusive ends. The name counts the window start in
    Mb, kb or bp, whichever divides `size` first, so 1 Mb windows are chr1_12M.
    """
    match = re.match(r'(chr[0-9XYM]+):(\d+)', genomic)
    if not match:
        return None
    chrom, pos_nt = match.group(1), int(match.group(2))
    bin_num = (pos_nt - 1) // size
    start, end = bin_num * size + 1, (bin_num + 1) * size
    unit, letter = (1_000_000, 'M') if size % 1_000_000 == 0 else (1_000, 'k') if size % 1_000 == 0 else (1, 'bp')
    return (f"{chrom}_{(bin_num * size) // unit}{letter}",
            f"Genomic window on {chrom} from {start} to {end} bp (1-based, inclusive)")


def collect_pfam_members(domain_file, map_cache, active_genes, pending=None, regions=None):
    """Pfam domain and clan terms: (pfam_terms, pfam_clan_terms).

    All overlapping domain and clan annotations are preserved. Each dict maps
    term -> (description, sorted "Gene:Position" member list). Without a domain file both are empty.
    When `pending` is a list, records that cannot be mapped to CAAS coordinates are appended to it
    instead of being dropped silently. When `regions` is a list, one record per run of consecutive CAAS
    columns of each domain and clan hit is appended to it (collection_id, term_id, term_label, gene,
    start, end, source_record_id, mapping_status, detail).
    """
    pfam_terms, pfam_clan_terms = {}, {}
    if not domain_file:
        return pfam_terms, pfam_clan_terms
    dom_cols = ['gene', 'pfam_id', 'target_name', 'description',
                'clan_acc', 'clan_name', 'ali_start', 'ali_end']
    df_dom = pd.read_csv(domain_file, sep='\t', keep_default_na=False)
    missing = set(dom_cols) - set(df_dom.columns)
    if missing:
        raise ValueError(f'Pfam table missing columns: {sorted(missing)}')
    if 'coordinate_space' not in df_dom:
        raise ValueError(f'Pfam table {domain_file} has no coordinate_space; '
                         'convert verified source coordinates with build_reference_tables.py '
                         'project-ppga or regenerate the annotation')
    for row in df_dom.itertuples(index=False):
        gene = str(row.gene).split('.')[0]
        if gene not in active_genes:
            continue

        pfam_id = str(row.pfam_id)
        clan_acc = str(row.clan_acc)
        clan_name = str(row.clan_name)
        target_name = str(row.target_name)
        desc = str(row.description) if row.description is not None else ''

        try:
            ali_start = int(row.ali_start)
            ali_end = int(row.ali_end)
        except (ValueError, TypeError):
            continue

        map_entry = map_cache.get(gene)
        residue_to_col = map_entry['residue_to_col'] if map_entry else {}
        if row.coordinate_space == 'CAAS_0based':
            columns = sorted({int(c) for c in str(row.alignment_positions).split(',') if c != ''})
            if any(c < 0 for c in columns):
                raise ValueError('Negative Pfam alignment position')
            if map_entry and not set(columns).issubset(map_entry['selected_cols']):
                raise ValueError(f'Pfam positions outside MAP for {gene}')
            mapping_status = 'direct_caas'
            if pending is not None and not map_entry:
                pending.append({'source': 'pfam', 'record_id': f"pfam:{gene}:{pfam_id}:{ali_start}-{ali_end}",
                                'gene': gene, 'raw_position': f"{ali_start}-{ali_end}",
                                'reason': 'unvalidated_no_map'})
        elif row.coordinate_space == 'reference_protein_1based':
            # A trimmed endpoint must not discard surviving interior residues.
            columns = sorted(c for residue, c in residue_to_col.items()
                             if ali_start <= residue <= ali_end)
            mapping_status = 'projected_via_map'
            if pending is not None and not map_entry:
                pending.append({'source': 'pfam', 'record_id': f"pfam:{gene}:{pfam_id}:{ali_start}-{ali_end}",
                                'gene': gene, 'raw_position': f"{ali_start}-{ali_end}",
                                'reason': 'no_map_reference_space'})
        else:
            raise ValueError(f'Unknown Pfam coordinate_space: {row.coordinate_space}')
        if not columns:
            continue
        hit_members = [f"{gene}:{c}" for c in columns]
        record_id = f"pfam:{gene}:{pfam_id}:{ali_start}-{ali_end}"

        if pfam_id not in pfam_terms:
            pfam_terms[pfam_id] = (f"{target_name} ({desc})", [])
        pfam_terms[pfam_id][1].extend(hit_members)

        detail = f"ali_start={ali_start};ali_end={ali_end};coordinate_space={row.coordinate_space}"
        has_clan = bool(clan_acc and clan_acc != 'NA')
        if regions is not None:
            for start, end in column_runs(columns):
                regions.append({'collection_id': 'pfam_domains', 'term_id': pfam_id, 'term_label': target_name,
                                'gene': gene, 'start': start, 'end': end, 'source_record_id': record_id,
                                'mapping_status': mapping_status, 'detail': detail})
                if has_clan:
                    regions.append({'collection_id': 'pfam_clans', 'term_id': clan_acc, 'term_label': clan_name,
                                    'gene': gene, 'start': start, 'end': end, 'source_record_id': record_id,
                                    'mapping_status': mapping_status, 'detail': detail})

        if has_clan:
            if clan_acc not in pfam_clan_terms:
                pfam_clan_terms[clan_acc] = (clan_name, [])
            pfam_clan_terms[clan_acc][1].extend(hit_members)

    for term in pfam_terms:
        pfam_terms[term] = (pfam_terms[term][0], sorted(set(pfam_terms[term][1])))
    for term in pfam_clan_terms:
        pfam_clan_terms[term] = (pfam_clan_terms[term][0], sorted(set(pfam_clan_terms[term][1])))
    return pfam_terms, pfam_clan_terms


def collect_ucr_members(ucr_file, map_cache, active_genes, pending=None):
    """UCR core/flank column sets and per-position region and variability.

    ucr_positions.tsv aggregates three detection methods (absolute, relative,
    sliding; bin/detect_ucr.py), each with its own window boundaries over the same
    conservation track, so a position can be core under one method and flank under
    another. Only the sliding method is read, which merges each conserved stretch
    into fewer, more contiguous blocks than the other two. Core and flank remain
    separate layers that may overlap. 'core' wins over flank in ucr_region_char
    whatever the chunk order. Returns (core_cols, flank_cols, region_char, variability_char),
    the first two as {gene: set of columns} and the others as {(gene, col): value}.
    When `pending` is a list, rows that cannot be projected to CAAS coordinates
    are appended to it.
    """
    gene_ucr_core_cols, gene_ucr_flank_cols = {}, {}
    ucr_region_char, ucr_variability_char = {}, {}
    if not ucr_file:
        return gene_ucr_core_cols, gene_ucr_flank_cols, ucr_region_char, ucr_variability_char
    # ucr_positions.tsv can be very large: only the needed columns are read, in
    # chunks, and rows are iterated with itertuples (iterrows is far slower).
    header = pd.read_csv(ucr_file, sep='\t', nrows=0).columns
    if 'coordinate_space' not in header:
        raise ValueError(f'UCR table {ucr_file} has no coordinate_space; '
                         'declare verified coordinates or convert the PPGA source with '
                         'build_reference_tables.py project-ppga')
    cols = ['gene', 'position', 'region_type', 'method', 'variability']
    if 'coordinate_space' in header:
        cols.append('coordinate_space')
    reader = pd.read_csv(ucr_file, sep='\t', usecols=cols,
                         dtype={'gene': str, 'region_type': str, 'method': str},
                         chunksize=500_000)
    for chunk in reader:
        chunk = chunk[chunk['method'] == 'sliding']
        for row in chunk.itertuples(index=False):
            gene = str(row.gene).split('.')[0]
            if gene not in active_genes:
                continue
            try:
                pos_residue = int(row.position)
            except (ValueError, TypeError):
                continue
            map_entry = map_cache.get(gene)
            record_id = f"ucr:{gene}:{pos_residue}:{row.region_type}:{row.method}"
            space = row.coordinate_space
            if space not in ('raw_alignment_1based', 'trimmed_alignment_1based', 'CAAS_0based'):
                raise ValueError(f'Unknown UCR coordinate_space: {space}')
            if space in ('CAAS_0based', 'trimmed_alignment_1based'):
                col = pos_residue - (space == 'trimmed_alignment_1based')
                if col < 0 or (map_entry and col not in map_entry['selected_cols']):
                    raise ValueError(f'UCR position outside MAP: {gene}:{col}')
            else:
                if not map_entry:
                    if pending is not None:
                        pending.append({'source': 'ucr', 'record_id': record_id, 'gene': gene,
                                        'raw_position': pos_residue, 'reason': 'no_map'})
                    continue
                col = map_entry['raw_to_col'].get(pos_residue)
            if col is None:
                if pending is not None:
                    pending.append({'source': 'ucr', 'record_id': record_id, 'gene': gene,
                                    'raw_position': pos_residue, 'reason': 'raw_position_not_in_map'})
                continue
            region = str(row.region_type)
            key = (gene, col)
            if region == 'core':
                gene_ucr_core_cols.setdefault(gene, set()).add(col)
                ucr_region_char[key] = 'core'
                ucr_variability_char[key] = row.variability
            elif region in ('flank_up', 'flank_down'):
                gene_ucr_flank_cols.setdefault(gene, set()).add(col)
                if ucr_region_char.get(key) != 'core':
                    ucr_region_char[key] = region
                    ucr_variability_char[key] = row.variability
    return gene_ucr_core_cols, gene_ucr_flank_cols, ucr_region_char, ucr_variability_char


def collect_fubar_members(fubar_file, map_cache, active_genes, pending=None):
    """FUBAR positive/purifying column sets and per-position calls.

    Unlike the two hit dicts, fubar_char records every position FUBAR tested, so
    'neutral' (tested, not significant) differs from a position absent from
    fubar_sites.tsv. Returns (positive_cols, purifying_cols, fubar_char). When `pending` is a
    list, rows that cannot be projected to CAAS coordinates are appended to it.
    """
    gene_pos_sel_cols, gene_neg_sel_cols, fubar_char = {}, {}, {}
    if not fubar_file:
        return gene_pos_sel_cols, gene_neg_sel_cols, fubar_char
    header = pd.read_csv(fubar_file, sep='\t', nrows=0).columns
    cols = ['gene', 'is_pos_hit', 'is_neg_hit', 'hg38_aa_pos']
    if 'coordinate_space' in header:
        cols.append('coordinate_space')
    if 'alignment_position' in header:
        cols.append('alignment_position')
    reader = pd.read_csv(fubar_file, sep='\t', usecols=cols,
                         dtype={'gene': str}, chunksize=500_000)
    for chunk in reader:
        for row in chunk.itertuples(index=False):
            gene = str(row.gene).split('.')[0]
            if gene not in active_genes:
                continue
            try:
                is_pos = int(row.is_pos_hit)
                is_neg = int(row.is_neg_hit)
            except (ValueError, TypeError):
                is_pos, is_neg = 0, 0

            space = getattr(row, 'coordinate_space', 'reference_protein_1based')
            try:
                pos_residue = int(row.hg38_aa_pos)
            except (ValueError, TypeError):
                if space != 'CAAS_0based':
                    continue
                pos_residue = 'NA'

            map_entry = map_cache.get(gene)
            record_id = (f"fubar:{gene}:CAAS:{row.alignment_position}" if space == 'CAAS_0based'
                         else f"fubar:{gene}:{pos_residue}")
            if space not in ('reference_protein_1based', 'CAAS_0based'):
                raise ValueError(f'Unknown FUBAR coordinate_space: {space}')
            if space == 'CAAS_0based':
                col = int(row.alignment_position)
                if col < 0 or (map_entry and col not in map_entry['selected_cols']):
                    raise ValueError(f'FUBAR position outside MAP: {gene}:{col}')
            else:
                if not map_entry:
                    if pending is not None:
                        pending.append({'source': 'fubar', 'record_id': record_id, 'gene': gene,
                                        'raw_position': pos_residue, 'reason': 'no_map'})
                    continue
                col = map_entry['residue_to_col'].get(pos_residue)
            if col is None:
                if pending is not None:
                    pending.append({'source': 'fubar', 'record_id': record_id, 'gene': gene,
                                    'raw_position': pos_residue, 'reason': 'reference_residue_not_in_map'})
                continue
            if is_pos == 1:
                gene_pos_sel_cols.setdefault(gene, set()).add(col)
                fubar_char[(gene, col)] = 'positive'
            elif is_neg == 1:
                gene_neg_sel_cols.setdefault(gene, set()).add(col)
                fubar_char[(gene, col)] = 'purifying'
            else:
                fubar_char[(gene, col)] = 'neutral'
    return gene_pos_sel_cols, gene_neg_sel_cols, fubar_char


def collect_fade_sites(fade_file, map_cache, active_genes, pending=None):
    """FADE selected sites (positions already in CAAS space) as {gene: set of columns}.

    FADE upstream already applies its configured BF threshold; no new cutoff here.
    When `pending` is a list, rows that cannot be validated against the MAP are
    appended to it.
    """
    result = {}
    if fade_file:
        frame = pd.read_csv(fade_file, usecols=['gene', 'position'])
        for row in frame.itertuples(index=False):
            gene, pos = str(row.gene).split('.')[0], int(row.position)
            if gene not in active_genes:
                continue
            if pos < 0:
                if pending is not None:
                    pending.append({'source': 'fade', 'record_id': f"fade:{gene}:{pos}", 'gene': gene,
                                    'raw_position': pos, 'reason': 'negative_position'})
                continue
            if gene in map_cache and pos not in map_cache[gene]['selected_cols']:
                if pending is not None:
                    pending.append({'source': 'fade', 'record_id': f"fade:{gene}:{pos}", 'gene': gene,
                                    'raw_position': pos, 'reason': 'position_outside_map'})
                continue
            result.setdefault(gene, set()).add(pos)
    return result


def collect_genomic_windows(map_cache, window_size=1_000_000):
    """Genomic windows of `window_size` bp over the MAP genomic coordinates as terms.

    Windows are chromosome-separated with one-based inclusive ends. This is genomic
    proximity, distinct from protein proximity.
    """
    gen_terms = {}
    for gene, entry in map_cache.items():
        col_to_genomic = entry['col_to_genomic']
        for col, genomic in col_to_genomic.items():
            term = window_term(genomic, window_size)
            if term is None:
                continue
            term_name, desc = term
            if term_name not in gen_terms:
                gen_terms[term_name] = (desc, [])
            gen_terms[term_name][1].append(f"{gene}:{col}")

    for term in gen_terms:
        gen_terms[term] = (gen_terms[term][0], sorted(set(gen_terms[term][1])))
    return gen_terms


def collect_custom_members(custom_file, map_cache, active_genes, strict=True, pending=None):
    """Custom marker terms: TSV without header, gene, position, term, optional
    description; '#' lines ignored. Positions must already be validated CAAS
    coordinates. With strict=True (default) an invalid row raises; otherwise it is
    appended to `pending` and skipped, so a registry build can report it instead of
    aborting.
    """
    custom_terms, members = {}, []
    if not (custom_file and os.path.exists(custom_file)):
        return custom_terms, members
    with open_maybe_gz(custom_file, 'rt') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            gene = fields[0]
            if gene not in active_genes:
                continue
            pos = int(fields[1])
            term_name = fields[2]
            record_id = f"custom:{gene}:{pos}:{term_name}"
            if pos < 0:
                if strict:
                    raise ValueError('Custom positions must be zero-based CAAS positions')
                if pending is not None:
                    pending.append({'source': 'custom', 'record_id': record_id, 'gene': gene,
                                    'raw_position': pos, 'reason': 'negative_position'})
                continue
            if gene in map_cache and pos not in map_cache[gene]['selected_cols']:
                if strict:
                    raise ValueError(f'Custom position outside MAP: {gene}:{pos}')
                if pending is not None:
                    pending.append({'source': 'custom', 'record_id': record_id, 'gene': gene,
                                    'raw_position': pos, 'reason': 'position_outside_map'})
                continue
            desc = fields[3] if len(fields) > 3 else "Custom annotation term"

            if term_name not in custom_terms:
                custom_terms[term_name] = (desc, [])
            custom_terms[term_name][1].append(f"{gene}:{pos}")
            members.append({'collection_id': 'custom', 'term_id': term_name, 'gene': gene,
                            'position': pos, 'source_record_id': f"custom:{gene}:{pos}:{term_name}",
                            'term_label': term_name, 'mapping_status': 'direct_caas',
                            'detail': desc})

    for term in custom_terms:
        custom_terms[term] = (custom_terms[term][0], sorted(set(custom_terms[term][1])))
    return custom_terms, members
