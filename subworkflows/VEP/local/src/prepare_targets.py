#!/usr/bin/env python3
"""Validate SCORING and verified MAP references before any external annotation."""
import argparse
from pathlib import Path
from vep_common import digest, pack_targets, prepare, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores',required=True);p.add_argument('--maps',required=True)
    p.add_argument('--proteins');p.add_argument('--output',default='vep_targets.json')
    p.add_argument('--fingerprint',required=True)
    a=p.parse_args()
    targets=prepare(a.scores,a.maps,a.proteins)
    inputs={'scores':digest(a.scores),'content_fingerprint':a.fingerprint,
            'maps':{x.name:digest(x) for x in sorted(Path(a.maps).iterdir()) if x.is_file() and (x.name.endswith('.map.tsv') or x.name.endswith('.map.tsv.vep.json'))}}
    if a.proteins:inputs['proteins']=digest(a.proteins)
    write_json(a.output,pack_targets(targets,inputs))
    write_json('validated_inputs.json',dict(schema='phylophere.vep_inputs.v1',inputs=inputs,
               target_count=len(targets),references={t['Gene']:t['reference'] for t in targets if t['reference']}))

if __name__=='__main__':main()
