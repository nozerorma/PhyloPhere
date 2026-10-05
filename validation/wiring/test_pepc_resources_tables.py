"""pepc_resources_tables.py reads the report numbers from a run's trace, Nextflow log and RESAMPLE log."""
import importlib.util
import os
from pathlib import Path

import pytest

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
spec = importlib.util.spec_from_file_location("prt", ROOT / "validation/tier1/scripts/pepc_resources_tables.py")
prt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prt)

HEADER = "task_id\thash\tnative_id\tname\tstatus\texit\tsubmit\tduration\trealtime\t%cpu\tpeak_rss\tpeak_vmem\trchar\twchar\n"
STATS = ("Workflow completed > WorkflowStats[succeededCount=2; failedCount=0; ignoredCount=0; cachedCount=0; succeedDuration=1h 2m 3s; "
         "peakRunning=3; peakCpus=12; peakMemory=20 GB; ]")


def _run(tmp_path, tag, slow):
    results, work_run, work = tmp_path / tag / "results", tmp_path / tag / "run", tmp_path / tag / "work"
    (results / "pipeline_info").mkdir(parents=True)
    work_run.mkdir(parents=True)
    (work / "ab" / "cdef01").mkdir(parents=True)
    rows = [("1", "ab/cdef01", "CT:RESAMPLE (x)", "2026-10-06 00:00:01.000", "3m 50s", "609.5%", "2.5 GB"),
            ("2", "11/222222", "CAAS_CORE:CAAS_CORE_BATCHED (b)", "2026-10-06 00:04:00.000", slow, "31.5%", "3.2 GB"),
            ("3", "33/444444", "FADE:FADE_BATCHED (f)", "2026-10-06 00:05:00.000", "46.2s", "100.8%", "500 MB")]
    (results / "pipeline_info/execution_trace.txt").write_text(HEADER + "".join(
        f"{i}\t{h}\t1\t{n}\tCOMPLETED\t0\t{s}\t1s\t{rt}\t{cpu}\t{rss}\t1 GB\t1 MB\t1 MB\n" for i, h, n, s, rt, cpu, rss in rows))
    (work_run / ".nextflow.log").write_text("Oct-06 00:00:00.100 [main] DEBUG nextflow.cli.CmdRun - x\nOct-06 00:00:00.200 [main] DEBUG nextflow.Session - Run name: run_a1\n"
                                            f"Oct-06 00:10:30.900 [main] DEBUG n.trace.WorkflowStatsObserver - {STATS}\nOct-06 00:10:31.500 [main] INFO x - bye\n")
    (work / "ab/cdef01/.command.log").write_text(
        "[START] 2026-10-06 00:00:10 Harvesting pool (strategy: auto)\n"
        "[INFO] 2026-10-06 00:00:40 Pool filled entirely from Tier 1 (1000/1000) in 1012 draws\n"
        "[INFO] 2026-10-06 00:02:00 Design matching: 887/1000 candidates reach >= 100 FOP hypotheses (88.7%)\n"
        "[INFO] 2026-10-06 00:02:20 Design matching: 1021/1147 candidates reach >= 100 FOP hypotheses (89.0%)\n"
        "[START] 2026-10-06 00:02:20 FOP mirror harvest for 1000 accepted cycles (max_fop=100, 8 worker(s), batch=1000)\n"
        "[COMPLETE] 2026-10-06 00:03:50 FOP mirror: 1000 cycles\n")
    return {"results": str(results), "work_run": str(work_run), "work": str(work)}


def test_durations_and_sizes_are_read_in_seconds_and_bytes():
    assert prt.seconds("1h 2m 3s") == 3723 and prt.seconds("4m 28s") == 268 and prt.seconds("46.2s") == 46.2 and prt.seconds("245ms") == 0.245
    assert prt.nbytes("3.2 GB") == pytest.approx(3.2e9) and prt.nbytes("16.6 MB") == pytest.approx(16.6e6) and prt.nbytes("0") == 0.0


def test_the_tables_carry_the_numbers_of_the_files(tmp_path):
    g, p = _run(tmp_path, "g", "6m 33s"), _run(tmp_path, "p", "4m 26s")
    totals = prt.totals(g, p)
    assert "2026-10-06 00:00:00 → 00:10:31 (10m 31s)" in totals and "1h 2m 3s" in totals and "20 GB" in totals
    assert "3.2 GB (`CAAS_CORE:CAAS_CORE_BATCHED`)" in totals
    cost = prt.cost_centres(g, p, n=2).splitlines()
    assert cost[2].startswith("| `CAAS_CORE:CAAS_CORE_BATCHED` | 6m 33s | 31.5 % | 4m 26s |")      # the longest task first
    phases = prt.resample_phases(g, p)
    assert "| pool harvest | 30 s (1012 draws) |" in phases
    assert "1m 20s + 20 s top-up (1147 candidates, 89.0 % reach 100 hypotheses)" in phases
    assert "1m 30s (1000 cycles × 100 hypotheses, 8 workers)" in phases
    trace = prt.side_by_side(g, p).splitlines()
    assert trace[2].startswith("| `CT:RESAMPLE` | 3m 50s | 609.5% | 2.5 GB |")                    # ordered by submit time


def test_a_process_missing_from_one_run_is_marked_not_run(tmp_path):
    g, p = _run(tmp_path, "g", "6m 33s"), _run(tmp_path, "p", "4m 26s")
    trace = Path(p["results"]) / "pipeline_info/execution_trace.txt"
    trace.write_text("".join(l for l in trace.read_text().splitlines(True) if "FADE_BATCHED" not in l))
    assert "| `FADE:FADE_BATCHED` | 46.2s | 100.8% | 500 MB | not run |  |  |" in prt.side_by_side(g, p)
