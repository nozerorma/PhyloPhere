#!/usr/bin/env python3
# resources.py — Resources tab config: resource ceilings of the local and slurm profiles.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Resources: ceilings of the `local` and `slurm` profiles, SLURM executor settings and
per-process resource overrides.

The `local` and `slurm` profiles of nextflow.config set separate params.max_memory,
max_cpus and max_time ceilings (local: 64.GB/32/5.day; slurm: 128.GB/128/960.h). They
are passed as --max_cpus/--max_memory/--max_time on the generated `nextflow run`
command. Nextflow enforces them on every request (process.resourceLimits), retries
included, so machine size is not set per process.

process_overrides is a finer knob: per-process cpus and memory, empty by default.
conf/resources.config is the only source of per-process values, so a row here is a
deliberate deviation from it. The button of the Resources tab fills the list with the
current defaults (gui/resource_defaults.py) as a starting point. The rows are rendered
into a config loaded with `-c` and generated beside the run scripts (run_single.sh.j2),
because Nextflow reads per-process resource directives from config files only, not from
command-line flags.

Imported by: gui/models/project.py, gui/resource_defaults.py, gui/widgets/resource_table/,
gui/widgets/tabs/resources_tab.py
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

    # SLURM executor settings (executor.$slurm block of nextflow.config); the local
    # executor has no submission queue to throttle. The defaults follow the comment of
    # that block: the lab account's QOS caps concurrent CPUs at 100 cluster-wide, shared
    # by everything the lab runs, so queueSize and submitRateLimit are conservative.
    # Raising them without checking that cap risks QOSMaxCpuPerUserLimit rejections of
    # the submissions.
    slurm_queue_size: str = "8"
    slurm_submit_rate_limit: str = "30/1min"
    slurm_exit_read_timeout: str = "4h"

    # Per-process cpus and memory overrides (see the module docstring).
    process_overrides: list[ProcessResourceOverride] = field(default_factory=list)
