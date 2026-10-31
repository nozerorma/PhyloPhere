#!/usr/bin/env python3
"""Write genomic candidate SNVs; protein HGVS is deliberately not used offline."""
from vep_common import read_targets


def build(targets_file,output):
    data=read_targets(targets_file)
    if data.get('schema')!='phylophere.vep_targets.v1':raise ValueError('Invalid VEP targets')
    variants={v['variant_id']:v for t in data['targets'] if t['status']=='ready' for v in t['variants']}
    with open(output,'w') as f:
        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for vid,v in sorted(variants.items(),key=lambda kv:(kv[1]['chrom'],kv[1]['pos'],kv[0])):
            f.write(f"{v['chrom']}\t{v['pos']}\t{vid}\t{v['ref']}\t{v['alt']}\t.\t.\t.\n")
    return data,variants
