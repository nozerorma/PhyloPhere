#!/usr/bin/env python3
# build_position_gmt.py — Position-level gene sets (GMT) and characterization tables for POSENRICH.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""Build positional collections and annotations, without running an enrichment test.

Output IDs are Gene:Position, with zero-based CAAS trimmed-alignment positions.
MAP prot_ali_col and ori_codon_col are one-based; hg38_aa_pos is the explicit
reference protein residue. Pfam requires explicit coordinate_space: CAAS_0based
with alignment_positions, or reference_protein_1based with verified endpoints. UCR position indexes the raw alignment (one-based). FUBAR hg38_aa_pos
indexes the reference protein; FADE and custom positions already use CAAS space.

Retained collections: Pfam domains/clans, genomic 1 Mb windows, custom features,
and global UCR/FUBAR/FADE characterization layers. Orthogroup-derived collections
are retired. Legacy CLI flags for their databases are ignored with a warning;
COSMIC, PrimateAI and Ensembl VEP remain report annotations through the VEP module.
All inputs are optional. Plain text and gzip inputs are accepted.
"""

# ── Dependencies ──────────────────────────────────────────────────────────────
import os
import sys
import argparse
import pandas as pd
import glob
import re
import gzip
import csv


# ── Input handling ────────────────────────────────────────────────────────────


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


def validate_required_inputs(args):
    """Normalize the optional inputs of args in place.

    A NO_FILE* sentinel or "None" becomes None; a path that does not exist becomes None
    with a warning; a `.gz` sibling replaces the path. map_dir is cleared with a warning
    when it is not a directory.
    """
    optional = {
        "--gene_ensembl_file": "gene_ensembl_file",
        "--domain_variability_file": "domain_variability_file",
        "--ucr_positions_file": "ucr_positions_file",
        "--fubar_sites_file": "fubar_sites_file",
        "--egg_members_file": "egg_members_file",
        "--egg_annotations_file": "egg_annotations_file",
        "--cosmic_db": "cosmic_db",
        "--pai3d_db": "pai3d_db",
        "--cleaned_background": "cleaned_background",
        "--custom_marker_file": "custom_marker_file",
        "--fade_sites_top_file": "fade_sites_top_file",
        "--fade_sites_bottom_file": "fade_sites_bottom_file",
    }

    # map_dir is optional and only needed for genomic coordinates
    if args.map_dir:
        if os.path.basename(args.map_dir).startswith("NO_") or not os.path.isdir(args.map_dir):
            print(f"WARNING: optional input --map_dir was not provided or not found (skipping genomic coordinate mapping): {args.map_dir}", file=sys.stderr)
            args.map_dir = None

    for flag, attr in optional.items():
        given = getattr(args, attr, None)
        if given and (os.path.basename(str(given)).startswith("NO_") or str(given) == "None"):
            setattr(args, attr, None)
            continue
        if given:
            resolved = resolve_path(given)
            if resolved is None:
                print(f"WARNING: optional input {flag} was provided but does not "
                      f"exist (skipping): {given}", file=sys.stderr)
                setattr(args, attr, None)
            else:
                setattr(args, attr, resolved)


# ── CLI ───────────────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(description="Build position-level GMT files for FCS enrichment analysis.")
    parser.add_argument("--gene_ensembl_file", required=False, default=None, help="Path to ensembl_genes.output containing human_protein_id column")
    parser.add_argument("--domain_variability_file", required=False, default=None, help="Path to domain_variability.tsv")
    parser.add_argument("--ucr_positions_file", required=False, default=None, help="Path to ucr_positions.tsv")
    parser.add_argument("--fubar_sites_file", required=False, default=None, help="Path to fubar_sites.tsv")
    parser.add_argument("--egg_members_file", required=False, default=None, help="Path to eggNOG members.tsv")
    parser.add_argument("--egg_annotations_file", required=False, default=None, help="Path to eggNOG annotations.tsv")
    parser.add_argument("--map_dir", required=False, default=None, help="Optional directory containing <GENE>*.map.tsv files")
    parser.add_argument("--cosmic_db", required=False, default=None, help="Path to Cosmic_MutantCensus_v104_GRCh38.tsv.gz")
    parser.add_argument("--pai3d_db", required=False, default=None, help="Path to PrimateAI-3D.hg38.txt.gz")
    parser.add_argument("--cleaned_background", required=False, default=None, help="Optional list of genes tested (universe filter)")
    parser.add_argument("--custom_marker_file", required=False, default=None, help="Optional custom marker file (gene, position, term, desc)")
    parser.add_argument("--fade_sites_top_file", required=False, default=None, help="Optional fade_sites_top.csv (gene,position,max_bf,target_aa) from FADE_TABLES")
    parser.add_argument("--fade_sites_bottom_file", required=False, default=None, help="Optional fade_sites_bottom.csv (gene,position,max_bf,target_aa) from FADE_TABLES")
    parser.add_argument("--output_dir", required=True, help="Output directory for generated GMT files")
    return parser.parse_args()


# ── Loading ───────────────────────────────────────────────────────────────────


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


def load_ensembl_mapping(gene_ensembl_file):
    """Return (ENSP to gene, gene to ENSP) from the gene/human_protein_id table.

    The ENSP version suffix is dropped; rows with an empty or NA protein ID are skipped.
    """
    if not gene_ensembl_file or not os.path.exists(gene_ensembl_file):
        print("No valid gene_ensembl_file provided; skipping Ensembl mapping.")
        return {}, {}
    # Only the two columns needed are read; the vectorized split is faster than a row loop.
    df = pd.read_csv(gene_ensembl_file, sep='\t',
                     usecols=['gene', 'human_protein_id'], dtype=str)
    df = df.dropna(subset=['human_protein_id'])
    df = df[df['human_protein_id'] != 'NA']
    df['ensp_clean'] = df['human_protein_id'].str.split('.').str[0]
    ensp_to_gene = dict(zip(df['ensp_clean'], df['gene']))
    gene_to_ensp = dict(zip(df['gene'], df['ensp_clean']))
    print(f"Loaded {len(ensp_to_gene)} Ensembl protein-to-gene mappings.")
    return ensp_to_gene, gene_to_ensp


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


def write_gmt(output_path, terms):
    """Write {term: (description, members)} as GMT, sorted by term; empty terms are omitted."""
    with open(output_path, 'w') as f:
        for term_name, (desc, members) in sorted(terms.items()):
            if members:
                member_str = "\t".join(members)
                f.write(f"{term_name}\t{desc}\t{member_str}\n")
    print(f"Wrote {len(terms)} terms to {output_path}")


# ── Collection builders ───────────────────────────────────────────────────────
# Each builder returns the GMT-ready terms plus the normalized member records
# (collection_id, term_id, gene, position, source_record_id, mapping_status) that
# posenrich_registry.py writes to the common registry. GMT stays a compatible
# export; the normalized tables are the canonical record.


class _DiscardMembers(list):
    """Annotation-only exports do not retain the normalized registry copy in RAM."""
    def append(self, item):
        pass


def collect_pfam_members(domain_file, map_cache, active_genes, pending=None, include_members=True):
    """Pfam domain/clan terms, per-position characterization and normalized members.

    All overlapping domain and clan annotations are preserved. Returns
    (pfam_terms, pfam_clan_terms, pfam_char, members): the first two map
    term -> (description, sorted "Gene:Position" member list), pfam_char maps
    (gene, col) -> (domain target names, clan names). Without a domain file every
    structure is empty. When `pending` is a list, records that cannot be mapped to
    CAAS coordinates are appended to it instead of being dropped silently.
    """
    pfam_terms, pfam_clan_terms, pfam_char, members = {}, {}, {}, []
    if not include_members: members = _DiscardMembers()
    if not domain_file:
        return pfam_terms, pfam_clan_terms, pfam_char, members
    dom_cols = ['gene', 'pfam_id', 'target_name', 'description',
                'clan_acc', 'clan_name', 'ali_start', 'ali_end']
    df_dom = pd.read_csv(domain_file, sep='\t', keep_default_na=False)
    missing = set(dom_cols) - set(df_dom.columns)
    if missing:
        raise ValueError(f'Pfam table missing columns: {sorted(missing)}')
    if 'coordinate_space' not in df_dom:
        raise ValueError(
            f'Pfam table {domain_file!r} has no coordinate_space. '
            'Regenerate with compute_domain_variability.py on the productive trimmed CAAS alignments, '
            'or use migrate_pfam_coordinates.py with the original hmmscan query and domtblout '
            'only when the query exactly matches that reference. '
            'Do not label ungapped query endpoints as reference_protein_1based without verification.')
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

        for c in columns:
            domains, clans = pfam_char.setdefault((gene, c), (set(), set()))
            domains.add(target_name)
            if clan_name and clan_name != 'NA':
                clans.add(clan_name)
            members.append({'collection_id': 'pfam_domains', 'term_id': pfam_id, 'gene': gene,
                            'position': c, 'source_record_id': record_id,
                            'term_label': target_name, 'mapping_status': mapping_status,
                            'detail': f"ali_start={ali_start};ali_end={ali_end};coordinate_space={row.coordinate_space}"})

        if clan_acc and clan_acc != 'NA':
            if clan_acc not in pfam_clan_terms:
                pfam_clan_terms[clan_acc] = (clan_name, [])
            pfam_clan_terms[clan_acc][1].extend(hit_members)
            for c in columns:
                members.append({'collection_id': 'pfam_clans', 'term_id': clan_acc, 'gene': gene,
                                'position': c, 'source_record_id': record_id,
                                'term_label': clan_name, 'mapping_status': mapping_status,
                                'detail': f"ali_start={ali_start};ali_end={ali_end};coordinate_space={row.coordinate_space}"})

    for term in pfam_terms:
        pfam_terms[term] = (pfam_terms[term][0], sorted(set(pfam_terms[term][1])))
    for term in pfam_clan_terms:
        pfam_clan_terms[term] = (pfam_clan_terms[term][0], sorted(set(pfam_clan_terms[term][1])))
    return pfam_terms, pfam_clan_terms, pfam_char, members


def collect_ucr_members(ucr_file, map_cache, active_genes, pending=None, include_members=True):
    """UCR core/flank column sets, per-position region/variability and members.

    ucr_positions.tsv aggregates three detection methods (absolute, relative,
    sliding; bin/detect_ucr.py), each with its own window boundaries over the same
    conservation track, so a position can be core under one method and flank under
    another. Only the sliding method is read, which merges each conserved stretch
    into fewer, more contiguous blocks than the other two. Core and flank remain
    separate annotation layers that may overlap: each is tested against the
    background on its own, so they need not be disjoint. 'core' wins over flank in
    ucr_region_char whatever the chunk order. Members carry the single term
    'global' of collection 'ucr_core'/'ucr_flank'. When `pending` is a list, rows
    that cannot be projected to CAAS coordinates are appended to it.
    """
    gene_ucr_core_cols, gene_ucr_flank_cols = {}, {}
    ucr_region_char, ucr_variability_char, members = {}, {}, []
    if not include_members: members = _DiscardMembers()
    if not ucr_file:
        return gene_ucr_core_cols, gene_ucr_flank_cols, ucr_region_char, ucr_variability_char, members
    # ucr_positions.tsv can be very large: only the needed columns are read, in
    # chunks, and rows are iterated with itertuples (iterrows is far slower).
    header = pd.read_csv(ucr_file, sep='\t', nrows=0).columns
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
            space = getattr(row, 'coordinate_space', 'raw_alignment_1based')
            if space not in ('raw_alignment_1based', 'CAAS_0based'):
                raise ValueError(f'Unknown UCR coordinate_space: {space}')
            if space == 'CAAS_0based':
                col = pos_residue
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
                members.append({'collection_id': 'ucr_core', 'term_id': 'global', 'gene': gene,
                                'position': col,
                                'source_record_id': record_id,
                                'term_label': 'UCR_core', 'mapping_status': 'direct_caas' if space == 'CAAS_0based' else 'projected_via_map',
                                'detail': f"region=core;raw_position={pos_residue};variability={row.variability}"})
            elif region in ('flank_up', 'flank_down'):
                gene_ucr_flank_cols.setdefault(gene, set()).add(col)
                if ucr_region_char.get(key) != 'core':
                    ucr_region_char[key] = region
                    ucr_variability_char[key] = row.variability
                members.append({'collection_id': 'ucr_flank', 'term_id': 'global', 'gene': gene,
                                'position': col,
                                'source_record_id': record_id,
                                'term_label': 'UCR_flank', 'mapping_status': 'direct_caas' if space == 'CAAS_0based' else 'projected_via_map',
                                'detail': f"region={region};raw_position={pos_residue};variability={row.variability}"})
    return gene_ucr_core_cols, gene_ucr_flank_cols, ucr_region_char, ucr_variability_char, members


def collect_fubar_members(fubar_file, map_cache, active_genes, pending=None, include_members=True):
    """FUBAR positive/purifying column sets, per-position calls and members.

    Unlike the two hit dicts, fubar_char records every position FUBAR tested, so
    'neutral' (tested, not significant) differs from a position absent from
    fubar_sites.tsv. members only carry the significant hits, the single term
    'global' of collection 'fubar_positive'/'fubar_purifying'. When `pending` is a
    list, rows that cannot be projected to CAAS coordinates are appended to it.
    """
    gene_pos_sel_cols, gene_neg_sel_cols, fubar_char, members = {}, {}, {}, []
    if not include_members: members = _DiscardMembers()
    if not fubar_file:
        return gene_pos_sel_cols, gene_neg_sel_cols, fubar_char, members
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
                members.append({'collection_id': 'fubar_positive', 'term_id': 'global', 'gene': gene,
                                'position': col, 'source_record_id': record_id,
                                'term_label': 'FUBAR_positive', 'mapping_status': 'direct_caas' if space == 'CAAS_0based' else 'projected_via_map',
                                'detail': f"hg38_aa_pos={pos_residue};call=positive"})
            elif is_neg == 1:
                gene_neg_sel_cols.setdefault(gene, set()).add(col)
                fubar_char[(gene, col)] = 'negative'
                members.append({'collection_id': 'fubar_purifying', 'term_id': 'global', 'gene': gene,
                                'position': col, 'source_record_id': record_id,
                                'term_label': 'FUBAR_purifying', 'mapping_status': 'direct_caas' if space == 'CAAS_0based' else 'projected_via_map',
                                'detail': f"hg38_aa_pos={pos_residue};call=negative"})
            else:
                fubar_char[(gene, col)] = 'neutral'
            # Coverage collection: every tested position, so "neutral" (tested, not
            # significant) stays distinguishable from absent from fubar_sites.tsv.
            members.append({'collection_id': 'fubar_tested', 'term_id': 'global', 'gene': gene,
                            'position': col, 'source_record_id': record_id,
                            'term_label': 'FUBAR_tested', 'mapping_status': 'direct_caas' if space == 'CAAS_0based' else 'projected_via_map',
                            'detail': f"hg38_aa_pos={pos_residue};call={fubar_char[(gene, col)]}"})
    return gene_pos_sel_cols, gene_neg_sel_cols, fubar_char, members


def collect_fade_sites(fade_file, map_cache, active_genes, collection_id, pending=None, include_members=True):
    """FADE selected sites (positions already in CAAS space) as columns and members.

    FADE upstream already applies its configured BF threshold; no new cutoff here.
    When `pending` is a list, rows that cannot be validated against the MAP are
    appended to it.
    """
    result, members = {}, []
    if not include_members: members = _DiscardMembers()
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
            members.append({'collection_id': collection_id, 'term_id': 'global', 'gene': gene,
                            'position': pos, 'source_record_id': f"fade:{gene}:{pos}",
                            'term_label': collection_id, 'mapping_status': 'direct_caas',
                            'detail': ''})
    return result, members


def collect_genomic_windows(map_cache, include_members=True):
    """1 Mb genomic bins over the MAP genomic coordinates, and members.

    Bins are chromosome-separated with one-based inclusive ends. This is genomic
    proximity, distinct from protein proximity.
    """
    gen_terms, members = {}, []
    if not include_members: members = _DiscardMembers()
    for gene, entry in map_cache.items():
        col_to_genomic = entry['col_to_genomic']
        for col, genomic in col_to_genomic.items():
            match = re.match(r'(chr[0-9XYM]+):(\d+)', genomic)
            if match:
                chrom = match.group(1)
                pos_nt = int(match.group(2))
                bin_num = (pos_nt - 1) // 1000000
                bin_start = bin_num * 1000000 + 1
                bin_end = (bin_num + 1) * 1000000

                term_name = f"{chrom}_{bin_num}M"
                desc = f"Genomic bin on {chrom} from {bin_start} to {bin_end} bp (1-based, inclusive)"

                if term_name not in gen_terms:
                    gen_terms[term_name] = (desc, [])
                gen_terms[term_name][1].append(f"{gene}:{col}")
                members.append({'collection_id': 'genomic_windows_1mb', 'term_id': term_name,
                                'gene': gene, 'position': col,
                                'source_record_id': f"map:{gene}:{col}",
                                'term_label': term_name, 'mapping_status': 'projected_via_map',
                                'detail': f"coord={genomic}"})

    for term in gen_terms:
        gen_terms[term] = (gen_terms[term][0], sorted(set(gen_terms[term][1])))
    return gen_terms, members


def build_char_layers(gene_ucr_core_cols, gene_ucr_flank_cols, gene_pos_sel_cols,
                      gene_neg_sel_cols, gene_fade_top_cols, gene_fade_bottom_cols):
    """Global characterization layers: one term per layer over all genes.

    Building them does not select or run an inferential test.
    """
    def _global_positions(gene_cols):
        members = []
        for g, cols in gene_cols.items():
            members.extend(f"{g}:{c}" for c in cols)
        return sorted(set(members))

    return {
        "UCR_core":        ("Ultra-conserved core positions",       _global_positions(gene_ucr_core_cols)),
        "UCR_flank":       ("UCR flanking positions",               _global_positions(gene_ucr_flank_cols)),
        "FUBAR_positive":  ("FUBAR positive-selection sites (FDR)",  _global_positions(gene_pos_sel_cols)),
        "FUBAR_purifying": ("FUBAR purifying-selection sites",       _global_positions(gene_neg_sel_cols)),
        "FADE_top_sig":    ("FADE directional selection sites, top (BF >= threshold)",    _global_positions(gene_fade_top_cols)),
        "FADE_bottom_sig": ("FADE directional selection sites, bottom (BF >= threshold)", _global_positions(gene_fade_bottom_cols)),
    }


def write_characterization_table(path, pfam_char, ucr_region_char, ucr_variability_char,
                                 fubar_char, fade_top_keys, fade_bottom_keys):
    """One row per annotated (Gene, Position): Pfam domain/clan, UCR region and its
    variability, and the FUBAR call, flattened for direct Gene/Position joins.
    The rows are the union of the keys of all sources; a position missing from one
    source has blanks in its columns. Returns the number of rows written.
    """
    char_keys = set(pfam_char) | set(ucr_region_char) | set(fubar_char) | fade_top_keys | fade_bottom_keys
    with open(path, 'w') as f:
        f.write("Gene\tPosition\tpfam_domain\tpfam_clan\tucr_region\tposition_variability\tfubar_selection\tfade_top_sig\tfade_bottom_sig\n")
        for gene, col in sorted(char_keys):
            pfam_domain, pfam_clan = (';'.join(sorted(values)) for values in pfam_char.get((gene, col), (set(), set())))
            ucr_region = ucr_region_char.get((gene, col), '')
            variability = ucr_variability_char.get((gene, col), '')
            fubar_selection = fubar_char.get((gene, col), '')
            f.write(f"{gene}\t{col}\t{pfam_domain}\t{pfam_clan}\t{ucr_region}\t{variability}\t{fubar_selection}\t{int((gene, col) in fade_top_keys)}\t{int((gene, col) in fade_bottom_keys)}\n")
    return len(char_keys)


def collect_custom_members(custom_file, map_cache, active_genes, strict=True, pending=None, include_members=True):
    """Custom marker terms: TSV without header, gene, position, term, optional
    description; '#' lines ignored. Positions must already be validated CAAS
    coordinates. With strict=True (default) an invalid row raises; otherwise it is
    appended to `pending` and skipped, so a registry build can report it instead of
    aborting.
    """
    custom_terms, members = {}, []
    if not include_members: members = _DiscardMembers()
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


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    """Build every GMT and table that the available inputs allow."""
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Resolve the inputs before any work: a dangling path is reported once, up front,
    # and a `.gz` sibling (eggNOG is stored compressed) replaces the plain name.
    validate_required_inputs(args)
    for retired in ('egg_members_file', 'egg_annotations_file', 'cosmic_db', 'pai3d_db'):
        if getattr(args, retired):
            print(f"WARNING: --{retired} ignored: orthogroup-derived collections retired", file=sys.stderr)

    universe_genes = load_universe_genes(args.cleaned_background)
    ensp_to_gene, gene_to_ensp = load_ensembl_mapping(args.gene_ensembl_file)
    map_cache = build_map_cache(args.map_dir, universe_genes)

    active_genes = set(universe_genes) if universe_genes is not None else set(map_cache.keys())
    if not active_genes and universe_genes is None:
        if universe_genes:
            active_genes = set(universe_genes)
        else:
            active_genes = set(ensp_to_gene.values())

    # 1. Pfam domains and clans (all overlapping domain and clan annotations kept)
    print("Compiling PFAM Domains & Clans...")
    pfam_terms, pfam_clan_terms, pfam_char, _ = collect_pfam_members(
        args.domain_variability_file, map_cache, active_genes, include_members=False)
    if args.domain_variability_file:
        write_gmt(os.path.join(args.output_dir, "pfam_domains.gmt"), pfam_terms)
        write_gmt(os.path.join(args.output_dir, "pfam_clans.gmt"), pfam_clan_terms)

    # 2. UCR positions per gene, split by region_type (core, flank_up, flank_down).
    print("Loading UCR positions...")
    gene_ucr_core_cols, gene_ucr_flank_cols, ucr_region_char, ucr_variability_char, _ = \
        collect_ucr_members(args.ucr_positions_file, map_cache, active_genes, include_members=False)

    # 3. FUBAR selection positions per gene, split by the sign of selection
    print("Loading FUBAR selection positions...")
    gene_pos_sel_cols, gene_neg_sel_cols, fubar_char, _ = collect_fubar_members(
        args.fubar_sites_file, map_cache, active_genes, include_members=False)

    # FADE upstream already applies its configured BF threshold. No new cutoff here.
    gene_fade_top_cols, _ = collect_fade_sites(
        args.fade_sites_top_file, map_cache, active_genes, 'fade_top_sig', include_members=False)
    gene_fade_bottom_cols, _ = collect_fade_sites(
        args.fade_sites_bottom_file, map_cache, active_genes, 'fade_bottom_sig', include_members=False)

    # 5. Genomic locations (1 Mbp chromosome bins). Needs genomic coordinates, so it
    # writes nothing without a map_dir. It iterates map_cache and not active_genes,
    # which can hold genes without a MAP entry.
    print("Compiling Genomic Locations (1 Mbp bins)...")
    gen_terms, _ = collect_genomic_windows(map_cache, include_members=False)
    write_gmt(os.path.join(args.output_dir, "genomic_locations.gmt"), gen_terms)

    # Global characterization sets are written separately from the local GMTs.
    # Building them does not select or run an inferential test.
    print("Compiling characterization layers (global overlap sets)...")
    char_layers = build_char_layers(gene_ucr_core_cols, gene_ucr_flank_cols, gene_pos_sel_cols,
                                    gene_neg_sel_cols, gene_fade_top_cols, gene_fade_bottom_cols)
    write_gmt(os.path.join(args.output_dir, "characterization_layers.tsv"), char_layers)

    # 6.6 Position characterization: one row per (Gene, Position) with the Pfam domain and
    # clan, UCR region and its variability, and the FUBAR call, flattened for direct
    # Gene/Position joins. The rows are the union of the keys of all sources; a position
    # missing from one source has blanks in its columns.
    print("Compiling position characterization table...")
    fade_top_keys = {(g, c) for g, cols in gene_fade_top_cols.items() for c in cols}
    fade_bottom_keys = {(g, c) for g, cols in gene_fade_bottom_cols.items() for c in cols}
    n_char = write_characterization_table(
        os.path.join(args.output_dir, "position_characterization.tsv"),
        pfam_char, ucr_region_char, ucr_variability_char, fubar_char, fade_top_keys, fade_bottom_keys)
    print(f"Wrote {n_char} rows to position_characterization.tsv")

    # 7. Custom features: tab-separated gene, position, term, optional description; "#" lines ignored
    if args.custom_marker_file and os.path.exists(args.custom_marker_file):
        print("Compiling Custom Features...")
        custom_terms, _ = collect_custom_members(args.custom_marker_file, map_cache, active_genes, include_members=False)
        write_gmt(os.path.join(args.output_dir, "custom_features.gmt"), custom_terms)


if __name__ == "__main__":
    main()
