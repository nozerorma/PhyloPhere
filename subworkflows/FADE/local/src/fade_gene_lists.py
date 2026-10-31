#!/usr/bin/env python3
"""Literal gene identifiers and valid empty FADE backgrounds/hit lists."""
import argparse
import csv
from pathlib import Path
from fade_json_to_tables import number
from fade_merge_tables import SUMMARY, integer


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary',required=True);p.add_argument('--direction',required=True,choices=['top','bottom'])
    p.add_argument('--threshold',type=float);a=p.parse_args()
    rows=[];seen=set()
    with open(a.summary,newline='') as fh:
        reader=csv.DictReader(fh,delimiter='\t')
        if reader.fieldnames!=SUMMARY:raise ValueError('Invalid FADE summary header')
        for r in reader:
            gene=r['gene']
            if not gene or gene in seen:raise ValueError('Empty or duplicate gene identifier')
            seen.add(gene)
            if a.threshold is not None and (number(r['max_bf'])>=a.threshold)!=(integer(r['n_sig_site_aa'])>0):raise ValueError('FADE summary threshold mismatch')
            rows.append(r)
    Path('background.txt').write_text(''.join(r['gene']+'\n' for r in rows))
    Path(f'fade_{a.direction}_significant.txt').write_text(''.join(r['gene']+'\n' for r in rows if integer(r['n_sig_site_aa'])>0))


if __name__=='__main__':main()
