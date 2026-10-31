#!/usr/bin/env python3
"""Strict reduction of a complete single-partition FADE JSON to compact tables."""
import argparse
import json
import math
import sys

AA_LETTERS = 'ACDEFGHIKLMNPQRSTVWY'
GENES_COLUMNS = (['direction', 'gene', 'n_sites_tested', 'n_sig_site_aa', 'max_bf', 'top_target_aa', 'top_site']
                 + [f'sigcount_{aa}' for aa in AA_LETTERS])
SITES_COLUMNS = ['direction', 'gene', 'site', 'position', 'max_bf', 'target_aa']


def number(value):
    # HyPhy can serialize a saturated positive BF as an infinity string.
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise ValueError(f'Invalid Bayes factor: {value!r}')
    x = float(value)
    if math.isnan(x) or x < 0:
        raise ValueError(f'Invalid Bayes factor: {value!r}')
    return x


def reduce_json(js, gene, direction, threshold, expected_sites=None):
    if direction not in {'top', 'bottom'} or not math.isfinite(threshold) or threshold <= 0:
        raise ValueError('Invalid FADE direction or BF threshold')
    mle = js.get('MLE') or {}
    headers = mle.get('headers') or []
    bf_columns = [i for i,h in enumerate(headers) if isinstance(h, list) and h and h[0] in {'BayesFactor', 'BayesFactor[bias>0]'}]
    if len(bf_columns) != 1:
        raise ValueError('FADE JSON requires a unique BayesFactor header')
    col = bf_columns[0]
    content = mle.get('content')
    if not isinstance(content, dict) or set(content) != set(AA_LETTERS):
        raise ValueError('FADE JSON must contain all 20 target residues')
    lengths = set()
    for aa in AA_LETTERS:
        if not isinstance(content[aa], dict) or set(content[aa]) != {'0'}:
            raise ValueError('FADE requires exactly one partition, indexed 0')
        rows = content[aa]['0']
        if not isinstance(rows, list) or not rows or any(not isinstance(r,list) or len(r)!=len(headers) for r in rows):
            raise ValueError(f'Incomplete FADE rows for {aa}')
        lengths.add(len(rows))
    declared_sites = (js.get('input') or {}).get('number of sites')
    if len(lengths)!=1 or (expected_sites is not None and lengths!={expected_sites}):
        raise ValueError('FADE site coverage differs from the alignment')
    n = next(iter(lengths))
    if declared_sites is not None and declared_sites != n:
        raise ValueError('FADE input site count disagrees with MLE rows')
    best, best_site, best_aa = -1, None, None
    sig_counts, site_hits = {}, {}
    for aa in AA_LETTERS:
        sig_counts[aa] = 0
        for site, row in enumerate(content[aa]['0'], 1):
            bf = number(row[col])
            if bf > best:
                best, best_site, best_aa = bf, site, aa
            if bf >= threshold:
                sig_counts[aa] += 1
                if site not in site_hits or bf > site_hits[site][0]:
                    site_hits[site] = bf, aa
    hit = best >= threshold
    gene_row = [direction,gene,n,sum(sig_counts.values()),best,best_aa if hit else 'NA',best_site if hit else 'NA']
    gene_row += [sig_counts[aa] for aa in AA_LETTERS]
    sites = [[direction,gene,s,s-1,*site_hits[s]] for s in sorted(site_hits)]
    return gene_row, sites


def _line(values):
    def fmt(v):
        if isinstance(v,float):
            return 'Inf' if math.isinf(v) else format(v,'.17g')
        return str(v)
    return '\t'.join(fmt(v) for v in values)+'\n'


def write_tables(genes_path, sites_path, genes, sites):
    from pathlib import Path
    Path(genes_path).write_text(_line(GENES_COLUMNS)+''.join(_line(r) for r in genes))
    Path(sites_path).write_text(_line(SITES_COLUMNS)+''.join(_line(r) for r in sites))


def main():
    if len(sys.argv)==3 and sys.argv[1]=='--print-header':
        print(_line({'genes':GENES_COLUMNS,'sites':SITES_COLUMNS}[sys.argv[2]]),end='');return
    parser=argparse.ArgumentParser(description=__doc__)
    for flag in ['json','gene','direction','genes-out','sites-out']:
        parser.add_argument('--'+flag,required=True)
    parser.add_argument('--bf-threshold',type=float,default=100)
    parser.add_argument('--expected-sites',type=int)
    args=parser.parse_args()
    row,sites=reduce_json(json.loads(open(args.json).read()),args.gene,args.direction,args.bf_threshold,args.expected_sites)
    write_tables(args.genes_out,args.sites_out,[row],sites)


if __name__=='__main__':
    main()
