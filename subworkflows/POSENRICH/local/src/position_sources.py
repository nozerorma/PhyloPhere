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
import json
import glob
import gzip
import os
import re
import sys

import pandas as pd
from pathlib import Path
from posenrich_contract import read_table, identifier, position
sys.path.insert(0, str(Path(__file__).resolve().parents[3]/'VEP/local/src'))
import vep_common as vr


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
    from posenrich_universe import load_universe_genes as load_genes
    return set(load_genes(cleaned_background_path))


def parse_map_file(path, gene=None):
    """Validate alignment columns independently of verified biological coordinates."""
    path = Path(path)
    mapping = vr.load_map(path)
    gene = gene if gene is not None else path.name.removesuffix('.map.tsv')
    reference = vr.load_reference(path, gene, mapping)
    raw_to_col, raw_seen = {}, set()
    rows = read_table(path)
    if 'ori_codon_col' not in rows:
        raise ValueError('MAP missing ori_codon_col')
    for row in rows.to_dict('records'):
        raw = row['ori_codon_col']
        if not raw.isdigit() or int(raw) < 1:
            raise ValueError('Invalid raw alignment column')
        if int(raw) in raw_seen:
            raise ValueError('Duplicate original MAP column')
        raw_seen.add(int(raw))
        if row['status'] == 'selected':
            if int(raw) in raw_to_col:
                raise ValueError('Duplicate raw MAP column')
            raw_to_col[int(raw)] = int(row['prot_ali_col']) - 1
    residue_to_col, genomic = {}, {}
    if reference:
        for col, row in mapping.items():
            residue = row['aa_position']
            if residue is not None:
                if residue > len(reference['protein_sequence']):
                    raise ValueError('MAP residue outside verified protein')
                residue_to_col[residue] = col
                codon = reference['codons'].get(str(residue))
                if codon:
                    genomic[col] = codon['positions'][0]
    return list(mapping), genomic, residue_to_col, reference['strand'] if reference else 'NA', raw_to_col


def build_map_cache(map_dir, universe_genes):
    cache = {}
    if not map_dir:
        return cache
    if not Path(map_dir).is_dir():
        raise ValueError('MAP directory does not exist')
    genes = universe_genes if universe_genes is not None else {
        p.name.removesuffix('.map.tsv') for p in Path(map_dir).glob('*.map.tsv')}
    for gene in sorted(genes):
        path = vr.map_path(gene, map_dir)
        if path is None:
            continue
        selected, genomic, residue, strand, raw = parse_map_file(path, gene)
        reference = vr.load_reference(path, gene, vr.load_map(path))
        cache[gene] = dict(selected_cols=set(selected), col_to_genomic=genomic,
                           residue_to_col=residue, strand=strand, raw_to_col=raw,
                           reference_status='verified' if reference else 'unverified_reference',
                           reference_identity=({k: reference[k] for k in ('protein_id','transcript_id','gene_id',
                                               'annotation_release','assembly','species')} if reference else None),
                           map_file=path.name)
    return cache


# ── Writing ───────────────────────────────────────────────────────────────────



# ── Source readers ───────────────────────────────────────────────────────────────────


def verify_source_reference(path, gene, cache):
    manifest = Path(str(path)+'.manifest.json')
    if not cache or cache.get('reference_status') != 'verified':
        return False
    if not manifest.is_file():
        return False
    meta = json.loads(manifest.read_text())
    identity = meta.get('reference_identities', {}).get(gene)
    if identity is None:
        return False
    if identity != cache['reference_identity']:
        raise ValueError(f'Annotation and MAP reference identities differ for {gene}')
    return True


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
    df_dom = read_table(domain_file)
    missing = set(dom_cols) - set(df_dom.columns)
    if missing:
        raise ValueError(f'Pfam table missing columns: {sorted(missing)}')
    if 'coordinate_space' not in df_dom:
        raise ValueError(f'Pfam table {domain_file} has no coordinate_space; '
                         'convert verified source coordinates with build_reference_tables.py '
                         'project-ppga or regenerate the annotation')
    seen_pfam = {}
    for row in df_dom.itertuples(index=False):
        gene = identifier(row.gene)
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
        except (ValueError, TypeError) as exc:
            raise ValueError('Invalid Pfam coordinates') from exc
        if not 1 <= ali_start <= ali_end:
            raise ValueError('Invalid Pfam interval')

        pfam_terms.setdefault(pfam_id, (f'{target_name} ({desc})', []))
        if clan_acc and clan_acc != 'NA':
            pfam_clan_terms.setdefault(clan_acc, (clan_name, []))
        pk = (gene, pfam_id, ali_start, ali_end, row.coordinate_space)
        record = tuple(str(v) for v in row)
        if pk in seen_pfam:
            if seen_pfam[pk] != record:
                raise ValueError('Conflicting Pfam duplicate')
            continue
        seen_pfam[pk] = record
        map_entry = map_cache.get(gene)
        residue_to_col = map_entry['residue_to_col'] if map_entry else {}
        if row.coordinate_space == 'CAAS_0based':
            columns = sorted({position(c) for c in str(row.alignment_positions).split(',') if c != ''})
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
            if map_entry and not verify_source_reference(domain_file, gene, map_entry):
                if pending is not None:
                    pending.append(dict(source='pfam', record_id=f'{gene}:{pfam_id}', gene=gene,
                                        raw_position=f'{ali_start}-{ali_end}', reason='unverified_reference'))
                continue
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
                         dtype=str, keep_default_na=False,
                         chunksize=500_000)
    seen_ucr = {}
    for chunk in reader:
        chunk = chunk[chunk['method'] == 'sliding']
        for row in chunk.itertuples(index=False):
            gene = identifier(row.gene)
            if gene not in active_genes:
                continue
            try:
                pos_residue = int(row.position)
            except (ValueError, TypeError) as exc:
                raise ValueError('Invalid UCR coordinate') from exc
            if not str(row.position).isdigit():
                raise ValueError('UCR coordinate must be integer')
            uk = (gene, pos_residue, row.region_type, row.method)
            uv = (row.coordinate_space, row.variability)
            if uk in seen_ucr and seen_ucr[uk] != uv:
                raise ValueError('Conflicting UCR duplicate')
            seen_ucr[uk] = uv
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
            variability = float(row.variability)
            if not 0 <= variability <= 1:
                raise ValueError('Invalid UCR variability')
            key = (gene, col)
            if key in ucr_variability_char and float(ucr_variability_char[key]) != variability:
                raise ValueError('Conflicting UCR variability at projected position')
            if region not in ('core', 'flank_up', 'flank_down'):
                raise ValueError('Unknown UCR region type')
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
    if 'coordinate_space' not in header:
        raise ValueError('FUBAR requires explicit coordinate_space')
    cols = ['gene', 'is_pos_hit', 'is_neg_hit', 'hg38_aa_pos']
    if 'coordinate_space' in header:
        cols.append('coordinate_space')
    if 'alignment_position' in header:
        cols.append('alignment_position')
    reader = pd.read_csv(fubar_file, sep='\t', usecols=cols,
                         dtype=str, keep_default_na=False, chunksize=500_000)
    seen_fubar = {}
    for chunk in reader:
        for row in chunk.itertuples(index=False):
            gene = identifier(row.gene)
            if gene not in active_genes:
                continue
            try:
                is_pos = int(row.is_pos_hit)
                is_neg = int(row.is_neg_hit)
            except (ValueError, TypeError) as exc:
                raise ValueError('Invalid FUBAR call') from exc
            if str(row.is_pos_hit) not in ('0','1') or str(row.is_neg_hit) not in ('0','1') or is_pos + is_neg > 1:
                raise ValueError('Invalid/conflicting FUBAR call')

            space = getattr(row, 'coordinate_space', None)
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
                col = position(row.alignment_position)
                if col < 0 or (map_entry and col not in map_entry['selected_cols']):
                    raise ValueError(f'FUBAR position outside MAP: {gene}:{col}')
            else:
                if not map_entry:
                    if pending is not None:
                        pending.append({'source': 'fubar', 'record_id': record_id, 'gene': gene,
                                        'raw_position': pos_residue, 'reason': 'no_map'})
                    continue
                if not verify_source_reference(fubar_file, gene, map_entry):
                    if pending is not None:
                        pending.append(dict(source='fubar', record_id=record_id, gene=gene,
                                            raw_position=pos_residue, reason='unverified_reference'))
                    continue
                col = map_entry['residue_to_col'].get(pos_residue)
            if col is None:
                if pending is not None:
                    pending.append({'source': 'fubar', 'record_id': record_id, 'gene': gene,
                                    'raw_position': pos_residue, 'reason': 'reference_residue_not_in_map'})
                continue
            key = gene, col
            call = is_pos, is_neg
            if key in seen_fubar and seen_fubar[key] != call:
                raise ValueError('Conflicting FUBAR duplicate')
            seen_fubar[key] = call
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
        frame = read_table(fade_file, sep=',')
        if not {'gene', 'position'} <= set(frame):
            raise ValueError('FADE source missing gene/position')
        for row in frame.itertuples(index=False):
            gene, pos = identifier(row.gene), position(row.position)
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
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) not in (3,4):
                raise ValueError('Malformed custom source row')
            gene = identifier(fields[0])
            if gene not in active_genes:
                continue
            pos = int(fields[1])
            term_name = identifier(fields[2])
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
