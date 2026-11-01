#!/usr/bin/env python3
"""Validate SCORING and verified MAP references before any external annotation."""
import argparse
from pathlib import Path
from vep_common import pack_targets, prepare, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores',required=True);p.add_argument('--maps',required=True)
    p.add_argument('--proteins');p.add_argument('--alignments');p.add_argument('--output',default='vep_targets.json')
    p.add_argument('--fingerprint',required=True)
    a=p.parse_args()
    targets=prepare(a.scores,a.maps,a.proteins,a.alignments)
    inputs={'scores':str(a.scores),'maps':str(a.maps),'proteins':a.proteins,'alignments':a.alignments,'source_fingerprint':a.fingerprint}
    write_json(a.output,pack_targets(targets,inputs))
    write_json('validated_inputs.json',dict(schema='phylophere.vep_inputs.v1',inputs=inputs,
               target_count=len(targets),references={t['Gene']:t['reference'] for t in targets if t['reference']}))

if __name__=='__main__':main()
