#!/usr/bin/env python3
# build_position_gmt.py — Position-level gene sets (GMT) and characterization tables for POSENRICH.
# PhyloPhere | subworkflows/ENRICHMENT/local/src/

"""
BuildPositionGmt: builds the gene sets of alignment positions that posenrich_enrich.py
tests, plus the position annotation tables of the position report.

A position is identified as Gene:Column, where Column is the protein alignment column of
the MAP file (status "selected" rows of <GENE>*.map.tsv; the column of CAAS's Position).
A GMT line is term, description and the member position IDs, tab-separated. Sources:
Pfam domains and clans, 1 Mbp genomic bins, eggNOG orthogroups (all positions of the
member genes, and the subsets in UCR core or flank, under positive or purifying selection
by FUBAR, at COSMIC mutation sites and at PrimateAI-3D pathogenic sites), and custom
markers. Domains are mapped from ungapped reference residue numbers to columns through
the MAP files. Genes outside the universe (--cleaned_background) are dropped; when no MAP
file is available the active genes are the universe, or else the genes of the Ensembl map.

Every input is optional: a missing, dangling or NO_FILE* input is skipped with a warning
and the GMTs that need it are not written.

Called by:  POSENRICH_BUILD_GMT Nextflow process (posenrich.nf → build_position_gmt.py)
Inputs:     --gene_ensembl_file, --domain_variability_file, --ucr_positions_file, --fubar_sites_file,
            --egg_members_file, --egg_annotations_file, --map_dir, --cosmic_db, --pai3d_db,
            --cleaned_background, --custom_marker_file (a .gz sibling of a path is accepted)
            --fade_sites_top_file, --fade_sites_bottom_file are accepted but not read
Outputs:    <output_dir>/*.gmt  pfam_domains, pfam_clans, genomic_locations, orthogroups,
                ucr_core_orthogroups, ucr_flank_orthogroups, selection_pos_orthogroups,
                selection_neg_orthogroups, cosmic_orthogroups, pai3d_orthogroups, custom_features
            characterization_layers.tsv  global position sets in GMT layout (UCR core and flank,
                FUBAR positive and purifying, FADE top and bottom), tested by overlap and not ranked
            position_characterization.tsv  Gene, Position, pfam_domain, pfam_clan, ucr_region,
                position_variability, fubar_selection
            cosmic_coverage_genes.txt, pai3d_coverage_genes.txt  genes the database could annotate
"""

# ── Dependencies ──────────────────────────────────────────────────────────────
import os
import sys
import argparse
import pandas as pd
import glob
import re
import gzip


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
    }

    # map_dir is optional and only needed for genomic coordinates
    if args.map_dir:
        if args.map_dir.startswith("NO_FILE") or not os.path.isdir(args.map_dir):
            print(f"WARNING: optional input --map_dir was not provided or not found (skipping genomic coordinate mapping): {args.map_dir}", file=sys.stderr)
            args.map_dir = None

    for flag, attr in optional.items():
        given = getattr(args, attr, None)
        if given and (str(given).startswith("NO_FILE") or str(given) == "None"):
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
    parser.add_argument("--fade_sites_top_file", required=False, default=None, help="Optional fade_sites_top.csv (gene,position,max_bf,target_aa) from FADE_JSON_TO_CSV")
    parser.add_argument("--fade_sites_bottom_file", required=False, default=None, help="Optional fade_sites_bottom.csv (gene,position,max_bf,target_aa) from FADE_JSON_TO_CSV")
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
    """Parse one <GENE>.map.tsv; returns (selected_cols, col_to_genomic, residue_to_col, strand).

    Only rows with status "selected" count. Columns used: 2 status, 4 protein alignment
    column, 5 human amino acid (NA at a gap), 6 hg38 codon coordinate (chrN:pos).
    residue_to_col maps the 1-based index of the ungapped human residue to its column.
    The strand is inferred from the trend of the genomic coordinates along the columns.
    """
    selected_cols = []
    col_to_genomic = {}
    residue_to_col = {}
    non_gap_counter = 0
    coords = []

    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            status = fields[1]
            prot_col = fields[3]
            hg38_aa = fields[4]
            hg38_nt = fields[5] if len(fields) > 5 else 'NA'

            if status == 'selected':
                col = int(prot_col) # Map to prot_ali_col
                selected_cols.append(col)
                if hg38_nt != 'NA':
                    col_to_genomic[col] = hg38_nt
                    if ':' in hg38_nt:
                        try:
                            pos = int(hg38_nt.split(':')[1])
                            coords.append(pos)
                        except ValueError:
                            pass
                if hg38_aa != 'NA':
                    non_gap_counter += 1
                    residue_to_col[non_gap_counter] = col

    # Strand: the sign of the coordinate steps along increasing columns
    is_minus = False
    if len(coords) > 1:
        diffs = [coords[i+1] - coords[i] for i in range(len(coords)-1)]
        sign_sum = sum(1 if d > 0 else -1 for d in diffs if d != 0)
        is_minus = (sign_sum < 0)
    strand = '-' if is_minus else '+'

    return selected_cols, col_to_genomic, residue_to_col, strand


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
    for path in map_files:
        filename = os.path.basename(path)
        gene = filename.split('.')[0]
        if universe_genes is not None and gene not in universe_genes:
            continue
        try:
            selected_cols, col_to_genomic, residue_to_col, strand = parse_map_file(path)
            map_cache[gene] = {
                'selected_cols': selected_cols,
                'col_to_genomic': col_to_genomic,
                'residue_to_col': residue_to_col,
                'strand': strand
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


def write_gene_list(output_path, genes):
    """Write the genes, sorted, one per line."""
    with open(output_path, 'w') as f:
        for g in sorted(genes):
            f.write(f"{g}\n")
    print(f"Wrote {len(genes)} coverage genes to {output_path}")


# ── External databases ────────────────────────────────────────────────────────


def build_genomic_to_pos(map_cache):
    """Lookup (chrom, nucleotide position) -> {(gene, column)}, one entry per codon position.

    Shared by the coordinate-keyed databases (COSMIC, PAI3D). The three nucleotides of a
    codon follow the strand: forward from the stored coordinate on "+", backward on "-".
    """
    genomic_to_pos = {}
    for gene, cache in map_cache.items():
        strand = cache['strand']
        for col, genomic in cache['col_to_genomic'].items():
            match = re.match(r'(chr[0-9XYM]+):(\d+)', genomic)
            if match:
                chrom = match.group(1)
                coord = int(match.group(2))
                if strand == '+':
                    codon_positions = [coord, coord + 1, coord + 2]
                else:
                    codon_positions = [coord, coord - 1, coord - 2]

                for nt_pos in codon_positions:
                    key = (chrom, nt_pos)
                    if key not in genomic_to_pos:
                        genomic_to_pos[key] = set()
                    genomic_to_pos[key].add((gene, col))
    return genomic_to_pos


def scan_external_positions(db_path, chr_col_name, pos_col_name, genomic_to_pos,
                             filter_col_name=None):
    """Stream a gzip TSV keyed by genomic (chrom, pos) against `genomic_to_pos`.

    Returns (coverage_cols, filtered_cols):
      - coverage_cols: gene -> {column} for every row matching genomic_to_pos, whatever
        filter_col_name says. It answers which genes the database could have annotated at
        all, and defines the coverage-gene list used to restrict the enrichment background.
      - filtered_cols: gene -> {column} for the rows whose filter_col_name value contains
        "pathogenic" (case-insensitive). Without filter_col_name (COSMIC has no
        pathogenicity call; every somatic mutation counts) it equals coverage_cols.
    """
    coverage_cols = {}
    filtered_cols = {}
    matched = 0
    with gzip.open(db_path, 'rt') as f:
        header_line = f.readline().strip()
        header_cols = header_line.split('\t')
        try:
            chr_col = header_cols.index(chr_col_name)
            pos_col = header_cols.index(pos_col_name)
        except ValueError:
            print(f"Error: {db_path} is missing {chr_col_name}/{pos_col_name} columns.",
                  file=sys.stderr)
            return coverage_cols, filtered_cols
        filter_col = header_cols.index(filter_col_name) if filter_col_name else None

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) <= max(chr_col, pos_col):
                continue
            chrom = fields[chr_col]
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            try:
                pos_nt = int(fields[pos_col])
            except ValueError:
                continue

            key = (chrom, pos_nt)
            if key not in genomic_to_pos:
                continue

            is_filtered_match = True
            if filter_col is not None and filter_col < len(fields):
                is_filtered_match = 'pathogenic' in fields[filter_col].lower()

            for gene, col in genomic_to_pos[key]:
                coverage_cols.setdefault(gene, set()).add(col)
                matched += 1
                if is_filtered_match:
                    filtered_cols.setdefault(gene, set()).add(col)

    print(f"Mapped {matched} {os.path.basename(db_path)} rows to protein alignment columns.")
    return coverage_cols, filtered_cols


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    """Build every GMT and table that the available inputs allow."""
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Resolve the inputs before any work: a dangling path is reported once, up front,
    # and a `.gz` sibling (eggNOG is stored compressed) replaces the plain name.
    validate_required_inputs(args)

    universe_genes = load_universe_genes(args.cleaned_background)
    ensp_to_gene, gene_to_ensp = load_ensembl_mapping(args.gene_ensembl_file)
    map_cache = build_map_cache(args.map_dir, universe_genes)

    active_genes = set(map_cache.keys())
    if not active_genes:
        if universe_genes:
            active_genes = set(universe_genes)
        else:
            active_genes = set(ensp_to_gene.values())

    # 1. Pfam domains and clans
    print("Compiling PFAM Domains & Clans...")
    pfam_terms = {}
    pfam_clan_terms = {}
    # (gene, column) -> (pfam_domain, pfam_clan), for position_characterization.tsv. A
    # position covered by several domain hits keeps the last one read; the field is
    # descriptive only.
    pfam_char = {}

    if os.path.exists(args.domain_variability_file):
        dom_cols = ['gene', 'pfam_id', 'target_name', 'description',
                    'clan_acc', 'clan_name', 'ali_start', 'ali_end']
        df_dom = pd.read_csv(args.domain_variability_file, sep='\t',
                             usecols=dom_cols)
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
            start_col = residue_to_col.get(ali_start)
            end_col = residue_to_col.get(ali_end)
            selected_cols = map_entry['selected_cols'] if map_entry else []
            
            if start_col is not None and end_col is not None:
                columns = [c for c in selected_cols if start_col <= c <= end_col]
                members = [f"{gene}:{c}" for c in columns]
                
                if pfam_id not in pfam_terms:
                    pfam_terms[pfam_id] = (f"{target_name} ({desc})", [])
                pfam_terms[pfam_id][1].extend(members)

                for c in columns:
                    pfam_char[(gene, c)] = (target_name, clan_name if clan_name and clan_name != 'NA' else '')

                if clan_acc and clan_acc != 'NA':
                    if clan_acc not in pfam_clan_terms:
                        pfam_clan_terms[clan_acc] = (clan_name, [])
                    pfam_clan_terms[clan_acc][1].extend(members)

        for term in pfam_terms:
            pfam_terms[term] = (pfam_terms[term][0], sorted(list(set(pfam_terms[term][1]))))
        for term in pfam_clan_terms:
            pfam_clan_terms[term] = (pfam_clan_terms[term][0], sorted(list(set(pfam_clan_terms[term][1]))))

        write_gmt(os.path.join(args.output_dir, "pfam_domains.gmt"), pfam_terms)
        write_gmt(os.path.join(args.output_dir, "pfam_clans.gmt"), pfam_clan_terms)

    # 2. UCR positions per gene, split by region_type (core, flank_up, flank_down).
    #    ucr_positions.tsv aggregates three detection methods (absolute, relative,
    #    sliding; bin/detect_ucr.py), each with its own window boundaries over the same
    #    conservation track, so a position can be core under one method and flank under
    #    another. Only the sliding method is read, which merges each conserved stretch
    #    into fewer, more contiguous blocks than the other two. Core and flank remain
    #    separate annotation layers that may overlap: each is tested against the
    #    background on its own, so they need not be disjoint.
    print("Loading UCR positions...")
    gene_ucr_core_cols = {}    # region_type == core
    gene_ucr_flank_cols = {}   # region_type in {flank_up, flank_down}
    # (gene, column) -> region label ('core', 'flank_up' or 'flank_down') and the
    # variability of the position, for position_characterization.tsv. 'core' wins over
    # flank whenever both occur, whatever the chunk order.
    ucr_region_char = {}
    ucr_variability_char = {}
    if os.path.exists(args.ucr_positions_file):
        # ucr_positions.tsv can be very large: only the needed columns are read, in
        # chunks, and rows are iterated with itertuples (iterrows is far slower).
        reader = pd.read_csv(args.ucr_positions_file, sep='\t',
                             usecols=['gene', 'position', 'region_type', 'method', 'variability'],
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
                if not map_entry:
                    continue
                col = map_entry['residue_to_col'].get(pos_residue)
                if col is None:
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

    # 3. FUBAR selection positions per gene, split by the sign of selection
    print("Loading FUBAR selection positions...")
    gene_pos_sel_cols = {}   # positive selection (FDR)
    gene_neg_sel_cols = {}   # purifying selection
    # (gene, column) -> 'positive', 'negative' or 'neutral', for
    # position_characterization.tsv. Unlike the two hit dicts above, it records every
    # position FUBAR tested, so "neutral" (tested, not significant) differs from a
    # position absent from fubar_sites.tsv (absent from this dict).
    fubar_char = {}
    if os.path.exists(args.fubar_sites_file):
        reader = pd.read_csv(args.fubar_sites_file, sep='\t',
                             usecols=['gene', 'is_pos_hit',
                                      'is_neg_hit', 'hg38_aa_pos'],
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

                try:
                    pos_residue = int(row.hg38_aa_pos)
                except (ValueError, TypeError):
                    continue

                map_entry = map_cache.get(gene)
                if not map_entry:
                    continue
                col = map_entry['residue_to_col'].get(pos_residue)
                if col is None:
                    continue
                if is_pos == 1:
                    gene_pos_sel_cols.setdefault(gene, set()).add(col)
                    fubar_char[(gene, col)] = 'positive'
                elif is_neg == 1:
                    gene_neg_sel_cols.setdefault(gene, set()).add(col)
                    fubar_char[(gene, col)] = 'negative'
                else:
                    fubar_char[(gene, col)] = 'neutral'

    # 3.5 FADE directional selection positions per gene. The dicts are only initialized
    #     here: --fade_sites_top_file and --fade_sites_bottom_file are not read, so the
    #     FADE_top_sig and FADE_bottom_sig layers come out empty.
    print("Loading FADE directional selection positions...")
    gene_fade_top_cols = {}
    gene_fade_bottom_cols = {}

    # 4. COSMIC and PAI3D positions. Both are keyed by genomic coordinate, so they need
    #    the MAP coordinates (genomic_to_pos). Each also yields the list of genes it
    #    could have annotated at all (*_coverage_genes.txt), which posenrich_enrich.py
    #    uses to restrict the background of cosmic_orthogroups and pai3d_orthogroups
    #    instead of diluting the test with genes the database never observed.
    genomic_to_pos = None
    if (args.cosmic_db and os.path.exists(args.cosmic_db)) or \
       (args.pai3d_db and os.path.exists(args.pai3d_db)):
        genomic_to_pos = build_genomic_to_pos(map_cache)

    print("Loading COSMIC mutation positions...")
    gene_cosmic_cols = {}
    if args.cosmic_db and os.path.exists(args.cosmic_db):
        print(f"Streaming {args.cosmic_db}...")
        gene_cosmic_cols, _ = scan_external_positions(
            args.cosmic_db, 'CHROMOSOME', 'GENOME_START', genomic_to_pos)
        write_gene_list(os.path.join(args.output_dir, "cosmic_coverage_genes.txt"),
                         gene_cosmic_cols.keys())

    # 4b. PAI3D pathogenic positions per gene. Unlike COSMIC (every somatic mutation
    #     counts), PAI3D predicts the pathogenicity of each variant, so GMT membership is
    #     limited to variants whose `prediction` contains "pathogenic" (the substring
    #     rule of 14.Position_enrichment_report.Rmd's is_pathogenic; no score cutoff).
    #     Coverage includes every matched variant, since whether PAI3D could see a gene
    #     is a coverage question, not a pathogenicity one.
    print("Loading PAI3D pathogenicity positions...")
    gene_pai3d_cols = {}
    if args.pai3d_db and os.path.exists(args.pai3d_db):
        print(f"Streaming {args.pai3d_db}...")
        gene_pai3d_coverage, gene_pai3d_cols = scan_external_positions(
            args.pai3d_db, 'chr', 'pos', genomic_to_pos, filter_col_name='prediction')
        write_gene_list(os.path.join(args.output_dir, "pai3d_coverage_genes.txt"),
                         gene_pai3d_coverage.keys())

    # genomic_to_pos is the largest structure of the script (one entry per codon position
    # of the whole universe) and the two scans above are its only readers. It is released
    # here so it does not stay in memory through the steps below.
    genomic_to_pos = None

    # 5. Genomic locations (1 Mbp chromosome bins). Needs genomic coordinates, so it
    # writes nothing without a map_dir. It iterates map_cache and not active_genes,
    # which can hold genes without a MAP entry.
    print("Compiling Genomic Locations (1 Mbp bins)...")
    gen_terms = {}
    for gene, entry in map_cache.items():
        col_to_genomic = entry['col_to_genomic']
        for col, genomic in col_to_genomic.items():
            match = re.match(r'(chr[0-9XYM]+):(\d+)', genomic)
            if match:
                chrom = match.group(1)
                pos_nt = int(match.group(2))
                bin_num = pos_nt // 1000000
                bin_start = bin_num * 1000000
                bin_end = bin_start + 1000000
                
                term_name = f"{chrom}_{bin_num}M"
                desc = f"Genomic bin on {chrom} from {bin_start} to {bin_end} bp"
                
                if term_name not in gen_terms:
                    gen_terms[term_name] = (desc, [])
                gen_terms[term_name][1].append(f"{gene}:{col}")

    for term in gen_terms:
        gen_terms[term] = (gen_terms[term][0], sorted(list(set(gen_terms[term][1]))))
        
    write_gmt(os.path.join(args.output_dir, "genomic_locations.gmt"), gen_terms)

    # 6. eggNOG orthogroups: the baseline GMT (all positions of the member genes) and the
    #    restricted GMTs (UCR core and flank, positive and purifying selection, COSMIC, PAI3D).
    print("Compiling eggNOG Orthogroups...")
    ortho_terms = {}
    ortho_ucr_core_terms = {}
    ortho_ucr_flank_terms = {}
    ortho_pos_sel_terms = {}
    ortho_neg_sel_terms = {}
    ortho_cosmic_terms = {}
    ortho_pai3d_terms = {}
    
    # Orthogroup descriptions: annotations column 2 (id) and 4 (description)
    ortho_descs = {}
    if args.egg_annotations_file and os.path.exists(args.egg_annotations_file):
        with open_maybe_gz(args.egg_annotations_file, 'rt') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    ortho_descs[fields[1]] = fields[3]

    if args.egg_members_file and os.path.exists(args.egg_members_file):
        with open_maybe_gz(args.egg_members_file, 'rt') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                og_id = fields[1]
                members_list = fields[4].split(',')
                desc = ortho_descs.get(og_id, "No functional annotation")
                
                # Member genes of the orthogroup: members column 5, entries "<taxid>.<ENSP>"
                og_genes = []
                for m in members_list:
                    parts = m.split('.', 1)
                    if len(parts) >= 2:
                        ensp_id = parts[1]
                        gene = ensp_to_gene.get(ensp_id)
                        if gene and gene in active_genes:
                            og_genes.append(gene)
                
                if og_genes:
                    full_members = []
                    ucr_core_members = []
                    ucr_flank_members = []
                    pos_sel_members = []
                    neg_sel_members = []
                    cosmic_members = []
                    pai3d_members = []

                    for g in og_genes:
                        g_entry = map_cache.get(g)
                        if not g_entry:
                            continue
                        for col in g_entry['selected_cols']:
                            pos_id = f"{g}:{col}"
                            full_members.append(pos_id)
                            if col in gene_ucr_core_cols.get(g, set()):
                                ucr_core_members.append(pos_id)
                            if col in gene_ucr_flank_cols.get(g, set()):
                                ucr_flank_members.append(pos_id)
                            if col in gene_pos_sel_cols.get(g, set()):
                                pos_sel_members.append(pos_id)
                            if col in gene_neg_sel_cols.get(g, set()):
                                neg_sel_members.append(pos_id)
                            if col in gene_cosmic_cols.get(g, set()):
                                cosmic_members.append(pos_id)
                            if col in gene_pai3d_cols.get(g, set()):
                                pai3d_members.append(pos_id)

                    def _add(store, members, suffix):
                        if members:
                            if og_id not in store:
                                store[og_id] = (f"{desc}{suffix}", [])
                            store[og_id][1].extend(members)

                    _add(ortho_terms,           full_members,      "")
                    _add(ortho_ucr_core_terms,  ucr_core_members,  " [UCR core sites]")
                    _add(ortho_ucr_flank_terms, ucr_flank_members, " [UCR flank sites]")
                    _add(ortho_pos_sel_terms,   pos_sel_members,   " [positive selection sites]")
                    _add(ortho_neg_sel_terms,   neg_sel_members,   " [purifying selection sites]")
                    _add(ortho_cosmic_terms,    cosmic_members,    " [COSMIC somatic mutation sites]")
                    _add(ortho_pai3d_terms,     pai3d_members,     " [PAI3D pathogenic sites]")

        def _dedup(store):
            for t in store:
                store[t] = (store[t][0], sorted(set(store[t][1])))

        for _store in (ortho_terms, ortho_ucr_core_terms, ortho_ucr_flank_terms,
                       ortho_pos_sel_terms, ortho_neg_sel_terms, ortho_cosmic_terms,
                       ortho_pai3d_terms):
            _dedup(_store)

        write_gmt(os.path.join(args.output_dir, "orthogroups.gmt"), ortho_terms)
        write_gmt(os.path.join(args.output_dir, "ucr_core_orthogroups.gmt"), ortho_ucr_core_terms)
        write_gmt(os.path.join(args.output_dir, "ucr_flank_orthogroups.gmt"), ortho_ucr_flank_terms)
        write_gmt(os.path.join(args.output_dir, "selection_pos_orthogroups.gmt"), ortho_pos_sel_terms)
        write_gmt(os.path.join(args.output_dir, "selection_neg_orthogroups.gmt"), ortho_neg_sel_terms)

        if args.cosmic_db and os.path.exists(args.cosmic_db):
            write_gmt(os.path.join(args.output_dir, "cosmic_orthogroups.gmt"), ortho_cosmic_terms)

        if args.pai3d_db and os.path.exists(args.pai3d_db):
            write_gmt(os.path.join(args.output_dir, "pai3d_orthogroups.gmt"), ortho_pai3d_terms)

    # map_cache (comparable in size to genomic_to_pos), the orthogroup terms and the
    # per-gene COSMIC and PAI3D dicts are already written and not read below, which
    # reuses only gene_ucr_*_cols, gene_*_sel_cols and gene_fade_*_cols. They are released.
    map_cache = None
    ortho_terms = ortho_ucr_core_terms = ortho_ucr_flank_terms = None
    ortho_pos_sel_terms = ortho_neg_sel_terms = None
    ortho_cosmic_terms = ortho_pai3d_terms = ortho_descs = None
    gene_cosmic_cols = gene_pai3d_cols = gene_pai3d_coverage = None

    # 6.5 Characterization layers: global position sets for overlap tests, not for ranked
    #     FCS (they are too large and would dominate a Wilcoxon test). The `.tsv`
    #     extension keeps them out of the `*.gmt` glob; posenrich_enrich.py reads them
    #     through --characterization.
    print("Compiling characterization layers (global overlap sets)...")
    def _global_positions(gene_cols):
        members = []
        for g, cols in gene_cols.items():
            members.extend(f"{g}:{c}" for c in cols)
        return sorted(set(members))

    char_layers = {
        "UCR_core":        ("Ultra-conserved core positions",       _global_positions(gene_ucr_core_cols)),
        "UCR_flank":       ("UCR flanking positions",               _global_positions(gene_ucr_flank_cols)),
        "FUBAR_positive":  ("FUBAR positive-selection sites (FDR)",  _global_positions(gene_pos_sel_cols)),
        "FUBAR_purifying": ("FUBAR purifying-selection sites",       _global_positions(gene_neg_sel_cols)),
        "FADE_top_sig":    ("FADE directional selection sites, top (BF >= threshold)",    _global_positions(gene_fade_top_cols)),
        "FADE_bottom_sig": ("FADE directional selection sites, bottom (BF >= threshold)", _global_positions(gene_fade_bottom_cols)),
    }
    write_gmt(os.path.join(args.output_dir, "characterization_layers.tsv"), char_layers)

    # 6.6 Position characterization: one row per (Gene, Position) with the Pfam domain and
    # clan, UCR region and its variability, and the FUBAR call, flattened for direct
    # Gene/Position joins. The rows are the union of the keys of all sources; a position
    # missing from one source has blanks in its columns.
    print("Compiling position characterization table...")
    char_keys = set(pfam_char) | set(ucr_region_char) | set(fubar_char)
    with open(os.path.join(args.output_dir, "position_characterization.tsv"), 'w') as f:
        f.write("Gene\tPosition\tpfam_domain\tpfam_clan\tucr_region\tposition_variability\tfubar_selection\n")
        for gene, col in sorted(char_keys):
            pfam_domain, pfam_clan = pfam_char.get((gene, col), ('', ''))
            ucr_region = ucr_region_char.get((gene, col), '')
            variability = ucr_variability_char.get((gene, col), '')
            fubar_selection = fubar_char.get((gene, col), '')
            f.write(f"{gene}\t{col}\t{pfam_domain}\t{pfam_clan}\t{ucr_region}\t{variability}\t{fubar_selection}\n")
    print(f"Wrote {len(char_keys)} rows to position_characterization.tsv")

    # 7. Custom features: tab-separated gene, position, term, optional description; "#" lines ignored
    if args.custom_marker_file and os.path.exists(args.custom_marker_file):
        print("Compiling Custom Features...")
        custom_terms = {}
        with open(args.custom_marker_file, 'r') as f:
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
                desc = fields[3] if len(fields) > 3 else "Custom annotation term"
                
                if term_name not in custom_terms:
                    custom_terms[term_name] = (desc, [])
                custom_terms[term_name][1].append(f"{gene}:{pos}")

        for term in custom_terms:
            custom_terms[term] = (custom_terms[term][0], sorted(list(set(custom_terms[term][1]))))
            
        write_gmt(os.path.join(args.output_dir, "custom_features.gmt"), custom_terms)


if __name__ == "__main__":
    main()
