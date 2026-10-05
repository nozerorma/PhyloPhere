#!/usr/bin/env python3
"""The tables of the Tier 1 PEPC resource report, from the files two runs leave behind.

For each run (a results directory holding `pipeline_info/execution_trace.txt`, and the run's work directory holding
`.nextflow.log` and the task directories) it prints the markdown of: the totals, the dominant cost centres, the RESAMPLE phase
timings and the per-process trace side by side. The numbers are read, never typed.

    python3 validation/tier1/scripts/pepc_resources_tables.py \\
        --geno validation/tier1/output/pepc/results/c4_complete validation/tier1/output/pepc/work/c4_complete_nxf_run validation/tier1/output/pepc/work/c4_complete \\
        --pheno validation/tier1/output/pepc/results/c4_phenotypic_complete ... 
"""
import argparse
import csv
import re
from datetime import datetime
from pathlib import Path

UNITS = {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}


def seconds(text):
    """Seconds of a Nextflow duration ('1h 2m 3s', '4m 28s', '46.2s', '245ms')."""
    total = 0.0
    for value, unit in re.findall(r"([\d.]+)\s*(ms|h|m|s)", text or ""):
        total += float(value) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600}[unit]
    return total


def nbytes(text):
    m = re.match(r"([\d.]+)\s*([KMGT]?B)", text or "")
    return float(m.group(1)) * UNITS[m.group(2)] if m else 0.0


def read_trace(results):
    path = Path(results) / "pipeline_info/execution_trace.txt"
    return list(csv.DictReader(open(path), delimiter="\t"))


def stats(work_run):
    """Launch and completion time, run name and the WorkflowStats of a run's .nextflow.log."""
    text = (Path(work_run) / ".nextflow.log").read_text(errors="ignore")
    ts = lambda s: datetime.strptime(s, "%b-%d %H:%M:%S.%f").replace(microsecond=0, year=datetime.now().year)
    stamps = re.findall(r"^(\w{3}-\d\d \d\d:\d\d:\d\d\.\d+) \[", text, re.M)       # first and last line of the log
    done = re.search(r"Workflow completed > WorkflowStats\[(.*?)\]", text)
    fields = dict(kv.strip().split("=", 1) for kv in done.group(1).split(";") if "=" in kv)
    return {"run": re.search(r"Run name: (\S+)", text).group(1), "launch": ts(stamps[0]), "done": ts(stamps[-1]), "stats": fields}


def hms(delta):
    s = int(round(delta.total_seconds()))
    return f"{s // 60}m {s % 60:02d}s"


def totals(geno, pheno):
    g, p = stats(geno["work_run"]), stats(pheno["work_run"])
    gt, pt = read_trace(geno["results"]), read_trace(pheno["results"])
    fmt = lambda s: f"{s['launch']:%Y-%m-%d %H:%M:%S} → {s['done']:%H:%M:%S} ({hms(s['done'] - s['launch'])})"
    big = lambda t: max(t, key=lambda r: nbytes(r["peak_rss"]))
    lines = ["| | genotypic | phenotypic |", "|---|---|---|",
             f"| Nextflow run name | `{g['run']}` | `{p['run']}` |",
             f"| Launch → completion (`.nextflow.log`) | {fmt(g)} | {fmt(p)} |",
             f"| Processes succeeded / failed | {g['stats']['succeededCount']} / {g['stats']['failedCount']} | {p['stats']['succeededCount']} / {p['stats']['failedCount']} |",
             f"| Cached tasks | {g['stats']['cachedCount']} | {p['stats']['cachedCount']} |",
             f"| Nextflow `succeedDuration` (cumulative task time) | {g['stats']['succeedDuration']} | {p['stats']['succeedDuration']} |",
             f"| `peakRunning` / `peakCpus` | {g['stats']['peakRunning']} / {g['stats']['peakCpus']} | {p['stats']['peakRunning']} / {p['stats']['peakCpus']} |",
             f"| `peakMemory` (requested, not resident) | {g['stats']['peakMemory']} | {p['stats']['peakMemory']} |",
             f"| Largest single-task `peak_rss` | {big(gt)['peak_rss']} (`{big(gt)['name'].split(' (')[0]}`) | {big(pt)['peak_rss']} (`{big(pt)['name'].split(' (')[0]}`) |"]
    return "\n".join(lines)


def short(name):
    return name.split(" (")[0]


def cost_centres(geno, pheno, n=6):
    gt, pt = read_trace(geno["results"]), read_trace(pheno["results"])
    gp = {short(r["name"]): r for r in gt}
    pp = {short(r["name"]): r for r in pt}
    top = sorted(gp, key=lambda k: -seconds(gp[k]["realtime"]))[:n]
    lines = ["| process | genotypic realtime | %cpu | phenotypic realtime | %cpu |", "|---|---|---|---|---|"]
    for k in top:
        q = pp.get(k)
        lines.append(f"| `{k}` | {gp[k]['realtime']} | {gp[k]['%cpu'].replace('%', ' %')} | {q['realtime'] if q else 'not run'} | {q['%cpu'].replace('%', ' %') if q else ''} |")
    return "\n".join(lines)


def resample_phases(geno, pheno):
    def one(run):
        row = next(r for r in read_trace(run["results"]) if short(r["name"]) == "CT:RESAMPLE")
        d = next(Path(run["work"]).glob(f"{row['hash'].split('/')[0]}/{row['hash'].split('/')[1]}*"))
        log = (d / ".command.log").read_text()
        stamp = lambda m: datetime.strptime(m, "%Y-%m-%d %H:%M:%S")
        t0 = stamp(re.search(r"\[START\] (\S+ \S+) Harvesting pool", log).group(1))
        pool = re.findall(r"\[INFO\] (\S+ \S+) Pool filled entirely from Tier 1 \((\d+)/\d+\) in (\d+) draws", log)
        match = re.findall(r"\[INFO\] (\S+ \S+) Design matching: (\d+)/(\d+) candidates reach >= \d+ FOP hypotheses \(([\d.]+)%\)", log)
        mirror = re.search(r"\[START\] (\S+ \S+) FOP mirror harvest for (\d+) accepted cycles \(max_fop=(\d+), (\d+) worker", log)
        done = re.search(r"\[COMPLETE\] (\S+ \S+) FOP mirror", log)
        secs = lambda a, b: int(round((stamp(b) - stamp(a)).total_seconds()))
        fmt = lambda s: f"{s // 60}m {s % 60:02d}s" if s >= 60 else f"{s} s"
        first = secs(pool[0][0], match[0][0])
        topup = secs(match[0][0], match[1][0]) if len(match) > 1 else 0
        return {"pool": f"{fmt(secs(row and f'{t0:%Y-%m-%d %H:%M:%S}', pool[0][0]))} ({pool[0][2]} draws)",
                "match": f"{fmt(first)}" + (f" + {fmt(topup)} top-up ({match[1][2]} candidates, {match[1][3]} % reach {mirror.group(3)} hypotheses)" if len(match) > 1 else f" ({match[0][2]} candidates, {match[0][3]} %)"),
                "mirror": f"{fmt(secs(mirror.group(1), done.group(1)))} ({mirror.group(2)} cycles × {mirror.group(3)} hypotheses, {mirror.group(4)} workers)"}
    g, p = one(geno), one(pheno)
    return "\n".join(["| phase | genotypic | phenotypic |", "|---|---|---|",
                      f"| pool harvest | {g['pool']} | {p['pool']} |",
                      f"| design matching: FOP harvest of every candidate | {g['match']} | {p['match']} |",
                      f"| FOP mirror | {g['mirror']} | {p['mirror']} |"])


def side_by_side(geno, pheno):
    gt, pt = read_trace(geno["results"]), read_trace(pheno["results"])
    gp = {short(r["name"]): r for r in sorted(gt, key=lambda r: r["submit"])}
    pp = {short(r["name"]): r for r in sorted(pt, key=lambda r: r["submit"])}
    names = list(gp) + [k for k in pp if k not in gp]
    lines = ["| process | geno realtime | geno %cpu | geno peak_rss | pheno realtime | pheno %cpu | pheno peak_rss |", "|---|---|---|---|---|---|---|"]
    cell = lambda r: (r["realtime"], r["%cpu"], r["peak_rss"]) if r else ("not run", "", "")
    for k in names:
        a, b = cell(gp.get(k)), cell(pp.get(k))
        lines.append(f"| `{k}` | {a[0]} | {a[1]} | {a[2]} | {b[0]} | {b[1]} | {b[2]} |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("geno", "pheno"):
        ap.add_argument(f"--{k}", nargs=3, required=True, metavar=("RESULTS", "WORK_RUN", "WORK"))
    a = ap.parse_args()
    geno = dict(zip(("results", "work_run", "work"), a.geno))
    pheno = dict(zip(("results", "work_run", "work"), a.pheno))
    for title, fn in (("Totals", totals), ("Dominant cost centres", cost_centres), ("RESAMPLE phases", resample_phases), ("Per-process trace", side_by_side)):
        print(f"### {title}\n\n{fn(geno, pheno)}\n")


if __name__ == "__main__":
    main()
