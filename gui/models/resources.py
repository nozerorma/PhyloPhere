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
overrides, empty by default. conf/resources.config is the only source of per-process
values, so a row here is a deliberate deviation from it. The Resources tab's button
loads the current conf defaults into this list (see gui/resource_defaults.py) to edit.
Machine size is not expressed per process: local_max_* / slurm_max_* are ceilings that
Nextflow enforces on every request (process.resourceLimits), retries included.
Overrides are rendered into a `-c`-loaded config generated alongside the run scripts
(see run_single.sh.j2) rather than as --flags, since Nextflow only reads per-process
resource directives from a config file, never from the command line.
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

    # SLURM executor knobs (nextflow.config's executor.$slurm block) -- no local
    # equivalent, since the local executor has no submission queue to throttle.
    # Defaults mirror that block's own comment: the lab account's SLURM QOS caps
    # concurrent CPUs at 100 cluster-wide, shared across everything the lab runs,
    # so queueSize/submitRateLimit are deliberately conservative -- raising them
    # without also checking that ceiling risks the QOSMaxCpuPerUserLimit rejection
    # storm documented there.
    slurm_queue_size: str = "8"
    slurm_submit_rate_limit: str = "30/1min"
    slurm_exit_read_timeout: str = "4h"

    # Per-process cpus/memory overrides (see module docstring above).
    process_overrides: list[ProcessResourceOverride] = field(default_factory=list)
