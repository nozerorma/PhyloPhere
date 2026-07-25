#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import glob
import re
import gzip


def resolve_path(path):
    """Return an existing path for `path`, transparently accepting a `.gz`
    sibling. Returns None if neither the given path nor `<path>.gz` exists.
    `os.path.exists()` follows symlinks, so a DANGLING symlink (e.g. a Nextflow
    stage-in whose target only lives on the HPC) resolves to None here — that is
    exactly the failure mode that used to silently drop whole databases."""
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
    """Fail loudly (non-zero exit) if any REQUIRED input is missing or a dangling
    symlink, listing every offender at once. Optional inputs (cosmic, background,
    custom markers) are only warned about. Resolves `.gz` siblings in place so
    downstream code can open the returned paths directly.

    Rationale: the builder previously guarded every database with a silent
    `if os.path.exists(...)`, so a broken stage-in produced a partial GMT set
    (e.g. only genomic_locations) with exit code 0 — indistinguishable from a
    deliberately omitted input. Required inputs must break the run visibly."""
    required = {
        "--gene_ensembl_file": "gene_ensembl_file",
        "--domain_variability_file": "domain_variability_file",
        "--ucr_positions_file": "ucr_positions_file",
        "--fubar_sites_file": "fubar_sites_file",
        "--egg_members_file": "egg_members_file",
        "--egg_annotations_file": "egg_annotations_file",
    }
    optional = {
        "--cosmic_db": "cosmic_db",
        "--pai3d_db": "pai3d_db",
        "--cleaned_background": "cleaned_background",
        "--custom_marker_file": "custom_marker_file",
    }

    missing = []
    for flag, attr in required.items():
        given = getattr(args, attr)
        resolved = resolve_path(given)
        if resolved is None:
            missing.append(f"  {flag} {given}")
        else:
            setattr(args, attr, resolved)

    # map_dir is a directory, not a file.
    if not (args.map_dir and os.path.isdir(args.map_dir)):
        missing.append(f"  --map_dir {args.map_dir} (directory not found)")

    for flag, attr in optional.items():
        given = getattr(args, attr)
        if given:
            resolved = resolve_path(given)
            if resolved is None:
                print(f"WARNING: optional input {flag} was provided but does not "
                      f"exist (skipping): {given}", file=sys.stderr)
                setattr(args, attr, None)
            else:
                setattr(args, attr, resolved)

    if missing:
        print("ERROR: the following REQUIRED position-GMT inputs are missing or "
              "are dangling symlinks:", file=sys.stderr)
        for m in missing:
            print(m, file=sys.stderr)
        print("Refusing to build a partial GMT set. Check that these paths exist "
              "on the machine running the job (a symlink to an HPC-only path is "
              "dangling here).", file=sys.stderr)
        sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(description="Build position-level GMT files for FCS enrichment analysis.")
    parser.add_argument("--gene_ensembl_file", required=True, help="Path to ensembl_genes.output containing human_protein_id column")
    parser.add_argument("--domain_variability_file", required=True, help="Path to domain_variability.tsv")
    parser.add_argument("--ucr_positions_file", required=True, help="Path to ucr_positions.tsv")
    parser.add_argument("--fubar_sites_file", required=True, help="Path to fubar_sites.tsv")
    parser.add_argument("--egg_members_file", required=True, help="Path to 9443_members.tsv")
    parser.add_argument("--egg_annotations_file", required=True, help="Path to 9443_annotations.tsv")
    parser.add_argument("--map_dir", required=True, help="Directory containing <GENE>*.map.tsv files")
    parser.add_argument("--cosmic_db", required=False, help="Path to Cosmic_MutantCensus_v104_GRCh38.tsv.gz")
    parser.add_argument("--pai3d_db", required=False, help="Path to PrimateAI-3D.hg38.txt.gz")
    parser.add_argument("--cleaned_background", required=False, help="Optional list of genes tested (universe filter)")
    parser.add_argument("--custom_marker_file", required=False, help="Optional custom marker file (gene, position, term, desc)")
    parser.add_argument("--fade_sites_top_file", required=False, help="Optional fade_sites_top.csv (gene,position,max_bf,target_aa) from FADE_JSON_TO_CSV")
    parser.add_argument("--fade_sites_bottom_file", required=False, help="Optional fade_sites_bottom.csv (gene,position,max_bf,target_aa) from FADE_JSON_TO_CSV")
    parser.add_argument("--output_dir", required=True, help="Output directory for generated GMT files")
    return parser.parse_args()

def load_universe_genes(cleaned_background_path):
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
    # Only the two columns we need; vectorized split beats row-wise iteration.
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

    # Determine strand trend
    is_minus = False
    if len(coords) > 1:
        diffs = [coords[i+1] - coords[i] for i in range(len(coords)-1)]
        sign_sum = sum(1 if d > 0 else -1 for d in diffs if d != 0)
        is_minus = (sign_sum < 0)
    strand = '-' if is_minus else '+'

    return selected_cols, col_to_genomic, residue_to_col, strand

def build_map_cache(map_dir, universe_genes):
    map_cache = {}
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

def write_gmt(output_path, terms):
    with open(output_path, 'w') as f:
        for term_name, (desc, members) in sorted(terms.items()):
            if members:
                member_str = "\t".join(members)
                f.write(f"{term_name}\t{desc}\t{member_str}\n")
    print(f"Wrote {len(terms)} terms to {output_path}")


def write_gene_list(output_path, genes):
    with open(output_path, 'w') as f:
        for g in sorted(genes):
            f.write(f"{g}\n")
    print(f"Wrote {len(genes)} coverage genes to {output_path}")


def build_genomic_to_pos(map_cache):
    """Codon-level genomic coordinate -> {(gene, col), ...} lookup, shared by
    every external coordinate-keyed database (COSMIC, PAI3D, ...)."""
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
      - coverage_cols: gene -> {col, ...} for EVERY row matching genomic_to_pos,
        regardless of filter_col_name. This answers "which genes could this
        external database structurally have annotated at all" — used to build
        the coverage-gene list for background restriction, independent of
        whatever functional filter defines GMT membership.
      - filtered_cols: gene -> {col, ...} restricted to rows where
        filter_col_name's value contains "pathogenic" (case-insensitive). If
        filter_col_name is None (COSMIC has no pathogenicity call — every
        somatic mutation counts), filtered_cols == coverage_cols.
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

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Fail loudly on missing/dangling required inputs BEFORE doing any work,
    # and resolve any `.gz` siblings (e.g. eggNOG ships compressed locally).
    validate_required_inputs(args)

    universe_genes = load_universe_genes(args.cleaned_background)
    ensp_to_gene, gene_to_ensp = load_ensembl_mapping(args.gene_ensembl_file)
    map_cache = build_map_cache(args.map_dir, universe_genes)

    active_genes = set(map_cache.keys())

    # 1. PFAM Domains & Clans
    print("Compiling PFAM Domains & Clans...")
    pfam_terms = {}
    pfam_clan_terms = {}
    
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

            residue_to_col = map_cache[gene]['residue_to_col']
            start_col = residue_to_col.get(ali_start)
            end_col = residue_to_col.get(ali_end)
            
            if start_col is not None and end_col is not None:
                columns = [c for c in map_cache[gene]['selected_cols'] if start_col <= c <= end_col]
                members = [f"{gene}:{c}" for c in columns]
                
                if pfam_id not in pfam_terms:
                    pfam_terms[pfam_id] = (f"{target_name} ({desc})", [])
                pfam_terms[pfam_id][1].extend(members)
                
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

    # 2. Load UCR Positions per gene, SPLIT by region_type.
    #    region_type ∈ {core, flank_up, flank_down}. ucr_positions.tsv aggregates
    #    THREE independent UCR-detection methods (absolute / relative / sliding —
    #    see detect_ucr.py upstream); reading all three indiscriminately made
    #    "core" and "flank" massively overlapping (a position can be core under
    #    one method's window and flank under another's), since each method
    #    computes its own window boundaries over the same underlying
    #    conservation track. Restricted to `sliding` only, which empirically has
    #    the least residual self-overlap of the three (window-merging
    #    consolidates each conserved stretch into fewer, more contiguous blocks
    #    per gene than the other two methods). Core and flank are kept as
    #    independent, potentially-overlapping annotation layers — not a forced
    #    partition: each is Fisher-tested against the background on its own, so
    #    there is no statistical requirement that they be disjoint.
    print("Loading UCR positions...")
    gene_ucr_core_cols = {}    # region_type == core
    gene_ucr_flank_cols = {}   # region_type in {flank_up, flank_down}
    if os.path.exists(args.ucr_positions_file):
        # ucr_positions.tsv can be very large (~2 GB): read only the needed
        # columns in chunks and iterate with itertuples (row-wise iterrows over
        # this file is pathologically slow).
        reader = pd.read_csv(args.ucr_positions_file, sep='\t',
                             usecols=['gene', 'position', 'region_type', 'method'],
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

                col = map_cache[gene]['residue_to_col'].get(pos_residue)
                if col is None:
                    continue
                region = str(row.region_type)
                if region == 'core':
                    gene_ucr_core_cols.setdefault(gene, set()).add(col)
                elif region in ('flank_up', 'flank_down'):
                    gene_ucr_flank_cols.setdefault(gene, set()).add(col)

    # 3. Load FUBAR Selection Positions per gene, SPLIT by selection sign.
    #    Positive-selection sites (is_pos_hit_fdr) are expected to be ENRICHED for
    #    trait-associated change; purifying sites (is_neg_hit) quantify constraint.
    #    Kept as separate layers rather than a single "selection" set.
    print("Loading FUBAR selection positions...")
    gene_pos_sel_cols = {}   # positive selection (FDR)
    gene_neg_sel_cols = {}   # purifying selection
    if os.path.exists(args.fubar_sites_file):
        # fubar_sites.tsv is per-codon genome-wide (~750 MB): read only the
        # needed columns in chunks and iterate with itertuples.
        reader = pd.read_csv(args.fubar_sites_file, sep='\t',
                             usecols=['gene', 'is_pos_hit_fdr',
                                      'is_neg_hit', 'hg38_aa_pos'],
                             dtype={'gene': str}, chunksize=500_000)
        for chunk in reader:
            for row in chunk.itertuples(index=False):
                gene = str(row.gene).split('.')[0]
                if gene not in active_genes:
                    continue
                try:
                    is_pos = int(row.is_pos_hit_fdr)
                    is_neg = int(row.is_neg_hit)
                except (ValueError, TypeError):
                    is_pos, is_neg = 0, 0

                if is_pos != 1 and is_neg != 1:
                    continue
                try:
                    pos_residue = int(row.hg38_aa_pos)
                except (ValueError, TypeError):
                    continue

                col = map_cache[gene]['residue_to_col'].get(pos_residue)
                if col is None:
                    continue
                if is_pos == 1:
                    gene_pos_sel_cols.setdefault(gene, set()).add(col)
                if is_neg == 1:
                    gene_neg_sel_cols.setdefault(gene, set()).add(col)

    # 3.5 Load FADE directional selection positions per gene (from
    #     FADE_JSON_TO_CSV's fade_sites_{top,bottom}.csv: gene,position,max_bf,
    #     target_aa, already restricted to sites clearing the classic BF>=100
    #     bar). Unlike FUBAR's hg38_aa_pos (a 1-indexed residue number needing
    #     residue_to_col translation), FADE's `position` is already 0-indexed
    #     in the same coordinate space CAAS's own Position column uses directly
    #     (confirmed empirically: joining FADE and CAAS positions with no
    #     translation gives the expected ~70x position-level enrichment,
    #     matching posenrich_enrich.py's own "obs-scores Position is prot_ali_col
    #     space" convention) -- so no map_cache lookup here, same as gene:position
    #     custom markers below.
    print("Loading FADE directional selection positions...")
    gene_fade_top_cols = {}
    gene_fade_bottom_cols = {}
    for fade_file, fade_cols in ((args.fade_sites_top_file, gene_fade_top_cols),
                                  (args.fade_sites_bottom_file, gene_fade_bottom_cols)):
        if not fade_file or not os.path.exists(fade_file):
            continue
        fade_df = pd.read_csv(fade_file, dtype={'gene': str})
        for row in fade_df.itertuples(index=False):
            gene = row.gene
            if gene not in active_genes:
                continue
            try:
                pos = int(row.position)
            except (ValueError, TypeError):
                continue
            fade_cols.setdefault(gene, set()).add(pos)

    # 4. Load COSMIC and PAI3D positions per gene. Both are external,
    #    coordinate-keyed databases with incomplete real-world coverage (unlike
    #    the internal deterministic computations above), so besides GMT
    #    membership we also emit which active genes each database could
    #    structurally have annotated at all (coverage_genes.txt) — consumed by
    #    posenrich_enrich.py to restrict the enrichment background for these
    #    two GMTs specifically, instead of diluting the test with genes/positions
    #    COSMIC/PAI3D never had a chance to observe.
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

    # 4b. Load PAI3D pathogenic positions per gene. Unlike COSMIC (every
    #     somatic mutation counts as evidence), PAI3D is a per-variant
    #     pathogenicity predictor, so GMT membership is restricted to variants
    #     the `prediction` column calls pathogenic (matching the same
    #     substring-match convention as 11.Scoring_report.Rmd's is_pathogenic
    #     — no score-cutoff fallback). Coverage (for background restriction)
    #     is tracked separately and includes every matched variant regardless
    #     of predicted pathogenicity, since "could PAI3D see this gene" is a
    #     coverage question, not a pathogenicity one.
    print("Loading PAI3D pathogenicity positions...")
    gene_pai3d_cols = {}
    if args.pai3d_db and os.path.exists(args.pai3d_db):
        print(f"Streaming {args.pai3d_db}...")
        gene_pai3d_coverage, gene_pai3d_cols = scan_external_positions(
            args.pai3d_db, 'chr', 'pos', genomic_to_pos, filter_col_name='prediction')
        write_gene_list(os.path.join(args.output_dir, "pai3d_coverage_genes.txt"),
                         gene_pai3d_coverage.keys())

    # 5. Genomic Locations (1 Mbp Chromosome Bins)
    print("Compiling Genomic Locations (1 Mbp bins)...")
    gen_terms = {}
    for gene in active_genes:
        col_to_genomic = map_cache[gene]['col_to_genomic']
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

    # 6. Orthogroups (eggNOG) - Baseline + restrictive GMTs (UCR core/flank,
    #    positive/purifying selection, COSMIC, PAI3D).
    print("Compiling eggNOG Orthogroups...")
    ortho_terms = {}
    ortho_ucr_core_terms = {}
    ortho_ucr_flank_terms = {}
    ortho_pos_sel_terms = {}
    ortho_neg_sel_terms = {}
    ortho_cosmic_terms = {}
    ortho_pai3d_terms = {}
    
    # Load descriptions
    ortho_descs = {}
    if os.path.exists(args.egg_annotations_file):
        with open_maybe_gz(args.egg_annotations_file, 'rt') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    ortho_descs[fields[1]] = fields[3]

    if os.path.exists(args.egg_members_file):
        with open_maybe_gz(args.egg_members_file, 'rt') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                og_id = fields[1]
                members_list = fields[4].split(',')
                desc = ortho_descs.get(og_id, "No functional annotation")
                
                # Identify member genes
                og_genes = []
                for m in members_list:
                    if m.startswith('9606.ENSP'):
                        parts = m.split('.')
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
                        for col in map_cache[g]['selected_cols']:
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

    # 6.5 Characterization layers — global functional layers for the report's
    #     hypergeometric/Fisher overlap test (NOT ranked by FCS: they are far too
    #     large and would dominate a Wilcoxon test). Written with a `.tsv`
    #     extension so the FCS `*.gmt` glob never picks them up; consumed
    #     separately by the enrich script's characterization step.
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

    # 7. Custom Features
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
