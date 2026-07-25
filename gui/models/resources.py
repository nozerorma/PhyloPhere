#!/usr/bin/env python3
# resources.py — Resources tab config: local/slurm profile-level resource caps.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
nextflow.config's `local` and `slurm` profiles set different params.max_memory /
max_cpus / max_time ceilings (local: 12.GB/8/5.day; slurm: 128.GB/128/960.h) — these
are the actual per-profile knobs worth exposing, unlike the per-process-label values
in conf/resources.config and conf.slurm/resources.config, which were confirmed
byte-identical between the two profiles and are pipeline-internal, not user-tunable
per run. Exposed as --max_cpus/--max_memory/--max_time overrides on the generated
`nextflow run` invocation.
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass


@dataclass(kw_only=True)
class ResourcesConfig:
    # Defaults mirror nextflow.config's `local` profile.
    local_max_cpus: str = "8"
    local_max_memory: str = "12.GB"
    local_max_time: str = "5.day"

    # Defaults mirror nextflow.config's `slurm` profile.
    slurm_max_cpus: str = "128"
    slurm_max_memory: str = "128.GB"
    slurm_max_time: str = "960.h"
