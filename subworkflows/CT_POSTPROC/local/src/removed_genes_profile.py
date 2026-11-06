#!/usr/bin/env python3
# removed_genes_profile.py: per-gene and per-position profile of the genes of each gene-filter category
# (Normal, Extreme, Dubious, Both), for report 8.1.
"""
Reduces the large reference tables (variability, Pfam, UCR, FUBAR) to small tables restricted to the candidate genes.
Every reference is optional; a missing one drops only its columns.

Inputs:  --discovery     pooled discovery TSV (Gene, Position, side), 0-based alignment columns (CAAS_0based)
         --gene-stats    gene_stats.tsv (Gene, category, length, n_CAAS)
         --clusters      selected_clusters.tsv (Gene, Position, clustering_flag); 'Discarded' = cluster train
         --gene-variability  gene_variability.tsv of the PPGA run (per gene, whole alignment)
         --entropy-dir   per-gene <gene>.prot.entropy.tsv (variability, g per column; position = Position + 1)
         --domains / --ucr / --fubar  position references with coordinate_space CAAS_0based
Outputs: gene_profile.tsv, position_profile.tsv, pfam_by_category.tsv (in --out-dir)

Coordinates: the reference tables must declare coordinate_space CAAS_0based; any other space stops the run, so that
a column is never compared across spaces.
"""
import argparse, os, re, sys
import numpy as np
import pandas as pd

CHUNK = 2_000_000


SPECIES = re.compile(r'\.[A-Z][a-z]+_[a-z]+(_[A-Za-z0-9]+)*$')


def make_resolver(genes):
    """Maps a reference-table gene to the gene name of the run, or None. The tables of one family name genes as
    'A1BG', 'A1BG.Homo_sapiens' or 'A1BG.Homo_sapiens.prot'; the first form present among the run genes wins."""
    cache = {}

    def resolve(g):
        if g not in cache:
            c = g[:-5] if g.endswith('.prot') else g
            cache[g] = next((x for x in (g, c, SPECIES.sub('', c)) if x in genes), None)
        return cache[g]
    return resolve


def check_space(df, path):
    if 'coordinate_space' not in df or set(df['coordinate_space'].dropna()) != {'CAAS_0based'}:
        sys.exit(f'{path}: coordinate_space must be CAAS_0based')


def read_chunks(path, usecols, resolve):
    """Stream a large table, keeping the rows of the candidate genes (renamed to the run gene names)."""
    for ch in pd.read_csv(path, sep='\t', usecols=usecols, chunksize=CHUNK, dtype={'gene': str}):
        check_space(ch, path) if 'coordinate_space' in ch else None
        ch['gene'] = ch['gene'].map(resolve)
        yield ch[ch['gene'].notna()]


def main():
    a = argparse.ArgumentParser()
    a.add_argument('--discovery', required=True); a.add_argument('--gene-stats', required=True)
    a.add_argument('--clusters', default=None); a.add_argument('--gene-variability', default=None)
    a.add_argument('--entropy-dir', default=None); a.add_argument('--domains', default=None)
    a.add_argument('--ucr', default=None); a.add_argument('--fubar', default=None)
    a.add_argument('--out-dir', required=True)
    args = a.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    stats = pd.read_csv(args.gene_stats, sep='\t', usecols=['Gene', 'category', 'length', 'n_CAAS', 'stratum'])
    disc = pd.read_csv(args.discovery, sep='\t', usecols=['Gene', 'Position', 'side'], dtype={'Gene': str})
    disc = disc.drop_duplicates(['Gene', 'Position', 'side'])
    if args.clusters:
        cl = pd.read_csv(args.clusters, sep='\t', usecols=['Gene', 'Position', 'clustering_flag'])
        disc = disc.merge(cl, on=['Gene', 'Position'], how='left')
        disc['train'] = disc['clustering_flag'].eq('Discarded')
    else:
        disc['train'] = False
    disc = disc.merge(stats[['Gene', 'category']], on='Gene', how='left')
    if disc['category'].isna().any():
        sys.exit('discovery genes missing from gene_stats')
    genes = set(stats['Gene'])
    resolve = make_resolver(genes)
    pos = disc[['Gene', 'Position', 'side', 'train', 'category']].copy()
    key = pos['Gene'] + '|' + pos['Position'].astype(str)

    # ── Variability: per column (entropy files) and per gene (gene_variability.tsv) ──
    gene = stats.set_index('Gene')
    if args.gene_variability:
        gv = pd.read_csv(args.gene_variability, sep='\t')
        gv['Gene'] = gv['gene'].map(resolve)
        gv = gv.dropna(subset=['Gene']).set_index('Gene')[['n_positions', 'pct_variable', 'mean_variability', 'mean_gap_fraction']]
        gene = gene.join(gv.rename(columns={'n_positions': 'n_columns', 'mean_variability': 'gene_mean_variability',
                                            'mean_gap_fraction': 'gene_mean_gap'}))
    if args.entropy_dir:
        var, gap, rest = {}, {}, {}
        cols_rest = {}
        by_gene = {g: set(d) for g, d in pos.groupby('Gene')['Position']}
        for g in genes:
            f = os.path.join(args.entropy_dir, f'{g}.prot.entropy.tsv')
            if not os.path.exists(f):
                continue
            e = pd.read_csv(f, sep='\t', usecols=['position', 'variability', 'g'])
            e['position'] -= 1                      # entropy is 1-based, discovery 0-based
            caas = by_gene.get(g, set())
            is_c = e['position'].isin(caas)
            cols_rest[g] = (e.loc[~is_c, 'variability'].mean(), e.loc[~is_c, 'g'].mean())
            for p, v, gg in zip(e.loc[is_c, 'position'], e.loc[is_c, 'variability'], e.loc[is_c, 'g']):
                var[f'{g}|{p}'] = v; gap[f'{g}|{p}'] = gg
        pos['variability'] = key.map(var); pos['gap_fraction'] = key.map(gap)
        r = pd.DataFrame(cols_rest, index=['mean_variability_noncaas', 'mean_gap_noncaas']).T
        gene = gene.join(r)

    # ── Pfam: per column membership ──
    if args.domains:
        dom = pd.read_csv(args.domains, sep='\t', dtype=str)
        check_space(dom, args.domains)
        dom['gene'] = dom['gene'].map(resolve)
        dom = dom[dom['gene'].notna()].copy()
        rows = []
        for r in dom.itertuples():
            for p in str(r.alignment_positions).split(','):
                if p != '':
                    rows.append((r.gene, int(p), r.pfam_id, r.target_name, r.clan_name))
        cols = pd.DataFrame(rows, columns=['Gene', 'Position', 'pfam_id', 'target_name', 'clan_name'])
        hit = cols.drop_duplicates(['Gene', 'Position'])
        pos['in_pfam'] = (pos['Gene'] + '|' + pos['Position'].astype(str)).isin(hit['Gene'] + '|' + hit['Position'].astype(str))
        gene['n_pfam_columns'] = hit.groupby('Gene').size()
        gene['n_pfam_domains'] = dom.groupby('gene').size()
        # Descriptive domain/clan table: genes carrying it and CAAS positions inside it, by gene category.
        cat = stats.set_index('Gene')['category']
        cols['category'] = cols['Gene'].map(cat)
        cp = cols.merge(pos[['Gene', 'Position', 'train']], on=['Gene', 'Position'], how='inner')
        out = []
        for label, kcols in (('pfam', ['pfam_id', 'target_name']), ('clan', ['clan_name'])):
            g1 = cols.groupby(kcols + ['category'])['Gene'].nunique().rename('genes')
            g2 = cp.groupby(kcols + ['category']).size().rename('caas_positions')
            g3 = cp[cp['train']].groupby(kcols + ['category']).size().rename('train_positions')
            t = pd.concat([g1, g2, g3], axis=1).fillna(0).reset_index()
            t.insert(0, 'level', label)
            t['term'] = t[kcols[-1]] if label == 'clan' else t['pfam_id'] + ' ' + t['target_name']
            out.append(t[['level', 'term', 'category', 'genes', 'caas_positions', 'train_positions']])
        pd.concat(out).to_csv(os.path.join(args.out_dir, 'pfam_by_category.tsv'), sep='\t', index=False)

    # ── UCR: per column region (core wins over flank) ──
    if args.ucr:
        reg, ncore = {}, {}
        want = set(key)
        for ch in read_chunks(args.ucr, ['gene', 'position', 'region_type', 'coordinate_space'], resolve):
            ch['k'] = ch['gene'] + '|' + ch['position'].astype(str)
            for k, r in zip(ch['k'], ch['region_type']):
                if k in want and reg.get(k) != 'core':
                    reg[k] = 'core' if r == 'core' else 'flank'
            c = ch[ch['region_type'] == 'core'].drop_duplicates('k')
            for g, n in c.groupby('gene').size().items():
                ncore[g] = ncore.get(g, 0) + n      # a gene is split across chunks only at a boundary: duplicates are rare
        pos['ucr_region'] = key.map(reg).fillna('none')
        gene['n_ucr_core_columns'] = pd.Series(ncore)

    # ── FUBAR: per site calls ──
    if args.fubar:
        acc, rows = {}, []
        want = set(key)
        for ch in read_chunks(args.fubar, ['gene', 'alignment_position', 'prob_pos', 'prob_neg', 'is_pos_hit', 'is_neg_hit',
                                           'coordinate_space'], resolve):
            for g, d in ch.groupby('gene'):
                s = acc.setdefault(g, [0, 0, 0])
                s[0] += len(d); s[1] += int(d['is_pos_hit'].sum()); s[2] += int(d['is_neg_hit'].sum())
            k = ch['gene'] + '|' + ch['alignment_position'].astype(str)
            rows.append(ch.loc[k.isin(want), ['gene', 'alignment_position', 'prob_pos', 'prob_neg', 'is_pos_hit', 'is_neg_hit']])
        fu = pd.concat(rows).rename(columns={'gene': 'Gene', 'alignment_position': 'Position'})
        pos = pos.merge(fu.drop_duplicates(['Gene', 'Position']), on=['Gene', 'Position'], how='left')
        f = pd.DataFrame(acc, index=['n_fubar_sites', 'n_fubar_pos', 'n_fubar_neg']).T
        gene = gene.join(f)

    # ── Gene-level aggregates of the CAAS positions ──
    pos['category'] = pos['category'].astype(str)
    g = pos.groupby('Gene')
    gene['n_positions'] = g.size()
    gene['n_train'] = g['train'].sum()
    for c in ('variability', 'gap_fraction'):
        if c in pos:
            gene[f'caas_mean_{c}'] = g[c].mean()
            gene[f'train_mean_{c}'] = pos[pos['train']].groupby('Gene')[c].mean()
    gene.rename_axis('Gene').reset_index().to_csv(os.path.join(args.out_dir, 'gene_profile.tsv'), sep='\t', index=False)
    pos.to_csv(os.path.join(args.out_dir, 'position_profile.tsv'), sep='\t', index=False)
    print(f'genes {len(gene)} positions {len(pos)}', file=sys.stderr)


if __name__ == '__main__':
    main()
