#!/usr/bin/env python3
# resources.py — Resources tab config: local/slurm profile-level resource caps.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
nextflow.config's `local` and `slurm` profiles set different params.max_memory /
max_cpus / max_time ceilings (local: 64.GB/32/5.day; slurm: 128.GB/128/960.h) —
these are the profile-level knobs, exposed as --max_cpus/--max_memory/--max_time
overrides on the generated `nextflow run` invocation.

process_overrides is a different, finer-grained knob: per-process cpus/memory
overrides (was previously considered pipeline-internal and out of scope — see
conf/resources.config's ~50 withName/withLabel blocks — but per-run tuning turned
out to matter across very different machines: the shared Marvin SLURM cluster vs.
a 32cpu/64GB workstation vs. an 8cpu/16GB laptop). The three checked-in presets
(conf/resources.config.{slurm,local,local_lowspec}) are loaded wholesale via the
Resources tab's three preset buttons (see gui/resource_presets.py) into this list,
then freely hand-edited per row. Rendered into a `-c`-loaded override config
generated alongside the run scripts (see run_single.sh.j2) rather than as
--flags, since Nextflow only reads per-process resource directives from a config
file, never from the command line.
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass, field


@dataclass(kw_only=True)
class ProcessResourceOverride:
    """One `withName`/`withLabel` process-selector resource override row."""

    selector_type: str = "withName"  # "withName" | "withLabel"
    selector: str = ""
    cpus: str = ""  # blank = leave conf/resources.config's own value untouched
    memory: str = ""  # e.g. "8.GB"; blank = leave untouched


@dataclass(kw_only=True)
class ResourcesConfig:
    # Defaults mirror nextflow.config's `local` profile.
    local_max_cpus: str = "32"
    local_max_memory: str = "64.GB"
    local_max_time: str = "5.day"

    # Defaults mirror nextflow.config's `slurm` profile.
    slurm_max_cpus: str = "128"
    slurm_max_memory: str = "128.GB"
    slurm_max_time: str = "960.h"

    # Per-process cpus/memory overrides (see module docstring above).
    process_overrides: list[ProcessResourceOverride] = field(default_factory=list)
