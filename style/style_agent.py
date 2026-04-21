#!/usr/bin/env python3
# style_agent.py — PhyloPhere pipeline style and coherence enforcement agent.
# PhyloPhere | style/

"""
StyleAgent: Audits Python, Nextflow, and R/Rmd modules in the PhyloPhere pipeline
for style consistency, correct headers, and code quality.

Runs in dry-run mode by default — prints issues without modifying files.
Use --apply to write auto-fixable changes (formatting + header injection).

Called by: developer / CI
Inputs:    PhyloPhere project root (inferred from script location)
Outputs:   Annotated report to stdout; modified files when --apply is used

Usage:
    python style/style_agent.py [OPTIONS]

Options:
    --apply           Apply all auto-fixable changes
    --lang LANG       Restrict to: python | nextflow | r
    --path PATH       Restrict to a specific file or directory
    --no-format       Skip formatter checks (black / isort / styler)
    --no-lint         Skip linter checks (flake8 / lintr / vulture)
    --no-headers      Skip header presence checks
    --no-structure    Skip structural checks (docstrings, process blocks, Rmd frontmatter)
    --verbose         Show formatter diffs inline
    --summary         Print only the final summary line
"""

import argparse
import ast
import re
import subprocess
import sys
from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import List, Optional, Tuple

# ── Constants ─────────────────────────────────────────────────────────────────

ROOT = Path(__file__).resolve().parent.parent

EXCLUDE_DIRS = {
    "work", ".git", ".history", "__pycache__", ".eggs",
    "build", "dist", "Docker(deprecated)", "apptainer", "pipeline_info",
}

# The PhyloPhere ASCII art signature — presence checked for "full" headers.
_PHYLO_ASCII_SIGNATURE = "██████╗"
_PHYLO_BANNER_SIGNATURE = "PHYLOPHERE:"

FULL_HEADER_TEMPLATE_NF = """\
#!/usr/bin/env nextflow

/*
#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: {filename}
#
*/"""

FULL_HEADER_TEMPLATE_PY = """\
#!/usr/bin/env python3
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: {filename}
#"""

FULL_HEADER_TEMPLATE_R = """\
#!/usr/bin/env Rscript
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: {filename}
#"""

# ── Enumerations & Data Classes ───────────────────────────────────────────────


class FileRole(Enum):
    NF_MAIN     = auto()  # main.nf — full PhyloPhere header
    NF_HELP     = auto()  # workflows/help.nf — full PhyloPhere header
    NF_WORKFLOW = auto()  # workflows/*.nf — slim header + DSL2 required
    NF_PROCESS  = auto()  # subworkflows/**/*.nf — slim header + process checks
    PY_MAIN     = auto()  # main.py / *_main.py / disambiguation_main.py — full header
    PY_MODULE   = auto()  # standalone *.py called directly from NF — slim header
    PY_PKG      = auto()  # src/**/*.py or __init__.py — slim header, no shebang
    R_SCRIPT    = auto()  # *.R called from NF via Rscript — slim header + shebang
    R_UTILITY   = auto()  # obj/*.R sourced by other R scripts — slim header, no shebang
    R_REPORT    = auto()  # *.Rmd — frontmatter field check


@dataclass
class Issue:
    file: Path
    category: str   # header | format | lint | dead-code | structure
    severity: str   # error | warning | info
    message: str
    line: Optional[int] = None
    diff: Optional[str] = None
    fixable: bool = False


@dataclass
class FileResult:
    path: Path
    role: FileRole
    issues: List[Issue] = field(default_factory=list)

    @property
    def has_issues(self) -> bool:
        return bool(self.issues)

    def add(self, issue: Optional[Issue]) -> None:
        if issue is not None:
            self.issues.append(issue)

    def extend(self, issues: List[Issue]) -> None:
        self.issues.extend(issues)


# ── File Classification ───────────────────────────────────────────────────────


def classify_file(path: Path) -> FileRole:
    """Determine a file's role in the pipeline, which drives header and structure rules."""
    rel = path.relative_to(ROOT)
    parts = set(rel.parts)
    name = path.name

    if path.suffix == ".nf":
        if name == "main.nf" and rel.parent == Path("."):
            return FileRole.NF_MAIN
        if name == "help.nf" and "workflows" in parts:
            return FileRole.NF_HELP
        if "workflows" in parts and "subworkflows" not in parts:
            return FileRole.NF_WORKFLOW
        return FileRole.NF_PROCESS

    if path.suffix == ".py":
        if name in ("main.py", "__main__.py") or name.endswith("_main.py"):
            return FileRole.PY_MAIN
        if "src" in parts or name == "__init__.py":
            return FileRole.PY_PKG
        return FileRole.PY_MODULE

    if path.suffix == ".Rmd":
        return FileRole.R_REPORT

    if path.suffix == ".R":
        if "obj" in parts:
            return FileRole.R_UTILITY
        return FileRole.R_SCRIPT

    raise ValueError(f"Cannot classify: {path}")


# ── Header Utilities ──────────────────────────────────────────────────────────


def _slim_expected(path: Path, comment: str) -> Tuple[str, str]:
    """Return the two expected slim header lines (filename line + PhyloPhere line)."""
    rel = str(path.parent.relative_to(ROOT)).replace("\\", "/")
    return (
        f"{comment} {path.name} — [TODO: describe purpose]",
        f"{comment} PhyloPhere | {rel}/",
    )


def _has_slim_header(lines: List[str], name: str, comment: str) -> bool:
    """True if the two-line slim header signature is present in the first 10 lines."""
    sig1_prefix = f"{comment} {name} —"
    sig2_prefix = f"{comment} PhyloPhere |"
    for i, line in enumerate(lines[:10]):
        if line.rstrip().startswith(sig1_prefix):
            if i + 1 < len(lines) and lines[i + 1].rstrip().startswith(sig2_prefix):
                return True
    return False


def _has_full_header(content: str) -> bool:
    return _PHYLO_ASCII_SIGNATURE in content and _PHYLO_BANNER_SIGNATURE in content


def check_slim_header(path: Path, comment: str) -> Optional[Issue]:
    lines = path.read_text(errors="replace").splitlines()
    if _has_slim_header(lines, path.name, comment):
        return None
    l1, l2 = _slim_expected(path, comment)
    return Issue(
        file=path,
        category="header",
        severity="error",
        message=f"Missing slim header. Expected:\n    {l1}\n    {l2}",
        fixable=True,
    )


def check_full_header(path: Path, lang: str) -> Optional[Issue]:
    content = path.read_text(errors="replace")
    if _has_full_header(content):
        return None
    template = {"nf": FULL_HEADER_TEMPLATE_NF, "py": FULL_HEADER_TEMPLATE_PY, "r": FULL_HEADER_TEMPLATE_R}[lang]
    filled = template.format(filename=path.name)
    return Issue(
        file=path,
        category="header",
        severity="error",
        message=f"Missing full PhyloPhere header. See archetype or add:\n{filled}",
        fixable=False,
    )


def inject_slim_header(path: Path, comment: str, has_shebang: bool) -> None:
    """Inject the slim header after the shebang line, replacing any existing slim header."""
    content = path.read_text(errors="replace")
    lines = content.splitlines(keepends=True)
    name = path.name
    l1, l2 = _slim_expected(path, comment)

    out: List[str] = []
    start = 0

    if has_shebang and lines and lines[0].startswith("#!"):
        out.append(lines[0])
        start = 1

    # Eat blank lines right after shebang
    while start < len(lines) and lines[start].strip() == "":
        start += 1

    # Remove any pre-existing slim header lines (our two-line format only)
    sig1 = f"{comment} {name} —"
    sig2 = f"{comment} PhyloPhere |"
    for _ in range(4):
        if start < len(lines):
            s = lines[start].rstrip()
            if s.startswith(sig1) or s.startswith(sig2):
                start += 1
            else:
                break

    # Eat blank separator
    if start < len(lines) and lines[start].strip() == "":
        start += 1

    out.append(f"{l1}\n")
    out.append(f"{l2}\n")
    out.append("\n")
    out.extend(lines[start:])
    path.write_text("".join(out))


# ── Shell Utility ─────────────────────────────────────────────────────────────


def _run(cmd: List[str]) -> Tuple[int, str]:
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT))
    return result.returncode, (result.stdout + result.stderr).strip()


def _tool_ok(name: str) -> bool:
    rc, _ = _run(["which", name])
    return rc == 0


# ── Python Checks ─────────────────────────────────────────────────────────────


def check_python_format(path: Path, apply: bool) -> List[Issue]:
    issues: List[Issue] = []

    if _tool_ok("black"):
        if apply:
            rc, out = _run(["black", "--quiet", str(path)])
            if rc != 0:
                issues.append(Issue(path, "format", "error", f"black failed:\n{out}"))
        else:
            rc, out = _run(["black", "--check", "--diff", str(path)])
            if rc != 0:
                issues.append(Issue(path, "format", "warning",
                                    "Needs black formatting", diff=out, fixable=True))
    else:
        issues.append(Issue(path, "format", "info", "black not found — install with: pip install black"))

    if _tool_ok("isort"):
        if apply:
            _run(["isort", "--quiet", str(path)])
        else:
            rc, out = _run(["isort", "--check", "--diff", str(path)])
            if rc != 0:
                issues.append(Issue(path, "format", "warning",
                                    "Needs isort import ordering", diff=out, fixable=True))
    else:
        issues.append(Issue(path, "format", "info", "isort not found — install with: pip install isort"))

    return issues


def check_python_lint(path: Path) -> List[Issue]:
    if not _tool_ok("flake8"):
        return [Issue(path, "lint", "info", "flake8 not found — install with: pip install flake8")]
    rc, out = _run(["flake8", str(path)])
    issues: List[Issue] = []
    for raw in out.splitlines():
        # flake8 output: path:line:col: CODE message
        m = re.match(r".+?:(\d+):\d+:\s+(.+)", raw)
        if m:
            issues.append(Issue(path, "lint", "warning", m.group(2), line=int(m.group(1))))
        elif raw.strip():
            issues.append(Issue(path, "lint", "warning", raw))
    return issues


def check_python_deadcode(path: Path) -> List[Issue]:
    if not _tool_ok("vulture"):
        return []
    whitelist = ROOT / "style" / "vulture_whitelist.py"
    cmd = ["vulture", str(path)]
    if whitelist.exists():
        cmd.append(str(whitelist))
    rc, out = _run(cmd)
    return [Issue(path, "dead-code", "info", line) for line in out.splitlines() if line.strip()]


def check_python_docstring(path: Path, role: FileRole) -> Optional[Issue]:
    """Check for module docstring and expected sections."""
    try:
        tree = ast.parse(path.read_text(errors="replace"))
    except SyntaxError as e:
        return Issue(path, "structure", "error", f"Syntax error — cannot parse: {e}")

    if not tree.body:
        return None

    first = tree.body[0]
    has_docstring = (
        isinstance(first, ast.Expr)
        and isinstance(first.value, ast.Constant)
        and isinstance(first.value.value, str)
    )
    if not has_docstring:
        return Issue(path, "structure", "warning",
                     "Missing module-level docstring. See archetype for expected sections.")

    doc = first.value.value
    if role == FileRole.PY_MODULE:
        missing = [s for s in ("Called by:", "Inputs:", "Outputs:") if s not in doc]
        if missing:
            return Issue(path, "structure", "info",
                        f"Module docstring missing sections: {', '.join(missing)}")
    elif role == FileRole.PY_MAIN:
        missing = [s for s in ("Called by:", "Usage:") if s not in doc]
        if missing:
            return Issue(path, "structure", "info",
                        f"Main docstring missing sections: {', '.join(missing)}")
    return None


# ── Nextflow Checks ───────────────────────────────────────────────────────────

_NF_PROCESS_RE = re.compile(
    r'\bprocess\s+(\w+)\s*\{((?:[^{}]|\{[^{}]*\})*)\}', re.DOTALL
)


def check_nextflow_structure(path: Path, role: FileRole) -> List[Issue]:
    content = path.read_text(errors="replace")
    issues: List[Issue] = []

    # DSL2 declaration only required in top-level workflow files
    if role in (FileRole.NF_MAIN, FileRole.NF_HELP, FileRole.NF_WORKFLOW):
        if "nextflow.enable.dsl" not in content:
            issues.append(Issue(path, "structure", "error",
                                "Missing 'nextflow.enable.dsl = 2'"))

    for m in _NF_PROCESS_RE.finditer(content):
        pname, body = m.group(1), m.group(2)
        if "input:" not in body:
            issues.append(Issue(path, "structure", "warning",
                                f"Process {pname}: missing 'input:' block"))
        if "output:" not in body:
            issues.append(Issue(path, "structure", "warning",
                                f"Process {pname}: missing 'output:' block"))
        if not any(kw in body for kw in ("script:", "shell:", "exec:")):
            issues.append(Issue(path, "structure", "warning",
                                f"Process {pname}: missing 'script:' / 'shell:' block"))
        if "label" not in body:
            issues.append(Issue(path, "structure", "warning",
                                f"Process {pname}: missing 'label' directive — needed for resource allocation"))
    return issues


# ── R / Rmd Checks ────────────────────────────────────────────────────────────

_REQUIRED_RMD_FIELDS = ("title", "author", "output", "params")
_REQUIRED_HTML_OPTS  = ("self_contained", "toc", "code_folding", "theme")


def check_rmd_frontmatter(path: Path) -> List[Issue]:
    try:
        import yaml
    except ImportError:
        return [Issue(path, "structure", "info",
                      "PyYAML not installed — skipping Rmd frontmatter check (pip install pyyaml)")]

    content = path.read_text(errors="replace")
    if not content.startswith("---"):
        return [Issue(path, "structure", "error",
                      "Rmd missing YAML frontmatter — file must start with '---'")]

    end = content.find("\n---", 3)
    if end == -1:
        return [Issue(path, "structure", "error", "Rmd YAML frontmatter block not closed")]

    try:
        fm = yaml.safe_load(content[3:end]) or {}
    except Exception as e:
        return [Issue(path, "structure", "error", f"Rmd YAML parse error: {e}")]

    issues: List[Issue] = []
    for f in _REQUIRED_RMD_FIELDS:
        if f not in fm:
            issues.append(Issue(path, "structure", "warning",
                                f"Rmd YAML missing required field: '{f}'"))

    html_opts = (fm.get("output") or {})
    if isinstance(html_opts, dict):
        html_opts = html_opts.get("html_document") or {}
        if isinstance(html_opts, dict):
            for opt in _REQUIRED_HTML_OPTS:
                if opt not in html_opts:
                    issues.append(Issue(path, "structure", "info",
                                        f"html_document missing recommended option: '{opt}'"))
    return issues


def check_r_format(path: Path, apply: bool) -> List[Issue]:
    dry_flag = "on" if not apply else "off"
    rc, out = _run([
        "Rscript", "--vanilla", "-e",
        f"if (!requireNamespace('styler', quietly=TRUE)) quit(status=2); "
        f"res <- styler::style_file('{path}', dry='{dry_flag}'); "
        f"if (!apply(res['changed'], 1, isTRUE)) quit(status=0) else quit(status=1)",
    ])
    if rc == 2:
        return [Issue(path, "format", "info",
                      "styler not installed — skipping R format check (install.packages('styler'))")]
    if rc == 1:
        return [Issue(path, "format", "warning",
                      "styler found style issues" + (" — applied" if apply else " — run with --apply to fix"),
                      fixable=True)]
    return []


def check_r_lint(path: Path) -> List[Issue]:
    rc, out = _run([
        "Rscript", "--vanilla", "-e",
        f"if (!requireNamespace('lintr', quietly=TRUE)) quit(status=2); "
        f"lints <- lintr::lint('{path}'); "
        f"if (length(lints) > 0) {{ print(lints); quit(status=1) }}",
    ])
    if rc == 2:
        return [Issue(path, "lint", "info",
                      "lintr not installed — skipping R lint check (install.packages('lintr'))")]
    issues: List[Issue] = []
    for line in out.splitlines():
        m = re.match(r".*?:(\d+):\d+: \w+: (.+)", line)
        if m:
            issues.append(Issue(path, "lint", "warning", m.group(2), line=int(m.group(1))))
        elif line.strip() and not line.startswith("<"):
            issues.append(Issue(path, "lint", "warning", line))
    return issues


# ── Per-File Processing ───────────────────────────────────────────────────────


def process_python(path: Path, opts: argparse.Namespace) -> FileResult:
    role = classify_file(path)
    result = FileResult(path=path, role=role)
    has_shebang = role != FileRole.PY_PKG

    if not opts.no_headers:
        if role == FileRole.PY_MAIN:
            result.add(check_full_header(path, "py"))
        else:
            issue = check_slim_header(path, "#")
            if issue and issue.fixable and opts.apply:
                inject_slim_header(path, "#", has_shebang)
            else:
                result.add(issue)

    if not opts.no_structure:
        result.add(check_python_docstring(path, role))

    if not opts.no_format:
        result.extend(check_python_format(path, opts.apply))

    if not opts.no_lint:
        result.extend(check_python_lint(path))
        result.extend(check_python_deadcode(path))

    return result


def process_nextflow(path: Path, opts: argparse.Namespace) -> FileResult:
    role = classify_file(path)
    result = FileResult(path=path, role=role)

    if not opts.no_headers:
        if role in (FileRole.NF_MAIN, FileRole.NF_HELP):
            result.add(check_full_header(path, "nf"))
        else:
            issue = check_slim_header(path, "//")
            if issue and issue.fixable and opts.apply:
                inject_slim_header(path, "//", has_shebang=True)
            else:
                result.add(issue)

    if not opts.no_structure:
        result.extend(check_nextflow_structure(path, role))

    return result


def process_r(path: Path, opts: argparse.Namespace) -> FileResult:
    role = classify_file(path)
    result = FileResult(path=path, role=role)

    if not opts.no_headers:
        if role == FileRole.R_REPORT:
            if not opts.no_structure:
                result.extend(check_rmd_frontmatter(path))
        else:
            issue = check_slim_header(path, "#")
            if issue and issue.fixable and opts.apply:
                inject_slim_header(path, "#", has_shebang=(role == FileRole.R_SCRIPT))
            else:
                result.add(issue)

    if role != FileRole.R_REPORT:
        if not opts.no_format:
            result.extend(check_r_format(path, opts.apply))
        if not opts.no_lint:
            result.extend(check_r_lint(path))

    return result


# ── File Collection ───────────────────────────────────────────────────────────


def collect_files(suffixes: set, restrict: Optional[Path] = None) -> List[Path]:
    search = restrict or ROOT
    found = []
    for p in sorted(search.rglob("*")):
        if any(ex in p.parts for ex in EXCLUDE_DIRS):
            continue
        if p.suffix in suffixes and p.is_file():
            found.append(p)
    return found


# ── Terminal Formatting ───────────────────────────────────────────────────────

_C = {
    "error":   "\033[91m",
    "warning": "\033[93m",
    "info":    "\033[96m",
    "bold":    "\033[1m",
    "dim":     "\033[2m",
    "reset":   "\033[0m",
}


def _c(key: str, text: str) -> str:
    return f"{_C[key]}{text}{_C['reset']}"


# ── Reporting ─────────────────────────────────────────────────────────────────


def print_report(results: List[FileResult], verbose: bool, apply: bool) -> int:
    """Print grouped report. Returns 0 if clean, 1 if issues found."""
    buckets = {
        "Python":    [r for r in results if r.role.name.startswith("PY_")],
        "Nextflow":  [r for r in results if r.role.name.startswith("NF_")],
        "R / Rmd":   [r for r in results if r.role.name.startswith("R_")],
    }

    total_issues = 0
    total_files_with_issues = 0

    for lang, lang_results in buckets.items():
        if not lang_results:
            continue
        dirty = [r for r in lang_results if r.has_issues]
        print(f"\n{_c('bold', '═' * 64)}")
        print(f"{_c('bold', lang)}  —  {len(lang_results)} files checked, {len(dirty)} with issues")
        print(_c("bold", "═" * 64))

        for r in dirty:
            rel = r.path.relative_to(ROOT)
            print(f"\n  {_c('bold', str(rel))}  {_c('dim', '[' + r.role.name + ']')}")
            by_cat: dict = {}
            for issue in r.issues:
                by_cat.setdefault(issue.category, []).append(issue)

            for cat, cat_issues in by_cat.items():
                for issue in cat_issues:
                    prefix = f"    [{cat.upper()}]"
                    sev_line = f"  {_c(issue.severity, prefix)}  {issue.message}"
                    if issue.line:
                        sev_line += f"  {_c('dim', f'(line {issue.line})')}"
                    if issue.fixable and not apply:
                        sev_line += _c("dim", "  ✦ auto-fixable")
                    print(sev_line)
                    if verbose and issue.diff:
                        for dl in issue.diff.splitlines()[:25]:
                            print(f"      {_c('dim', dl)}")

        total_issues += sum(len(r.issues) for r in dirty)
        total_files_with_issues += len(dirty)

    mode = "APPLIED" if apply else "DRY-RUN"
    print(f"\n{_c('bold', '═' * 64)}")
    print(f"{_c('bold', f'SUMMARY  [{mode}]')}")
    print(_c("bold", "═" * 64))
    print(f"  Files with issues : {total_files_with_issues}")
    print(f"  Total issues      : {total_issues}")
    if not apply and total_issues:
        fixable = sum(1 for r in results for i in r.issues if i.fixable)
        print(f"  Auto-fixable      : {fixable}  (re-run with --apply)")
    print()
    return 1 if total_issues else 0


# ── Entry Point ───────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(
        description="PhyloPhere style enforcement agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Archetypes in style/archetypes/ show the expected structure for each file type.",
    )
    parser.add_argument("--apply", action="store_true",
                        help="Apply auto-fixable changes (formatting + header injection)")
    parser.add_argument("--lang", choices=["python", "nextflow", "r"],
                        help="Restrict check to a single language")
    parser.add_argument("--path", type=Path, metavar="PATH",
                        help="Restrict check to a specific file or directory")
    parser.add_argument("--no-format",    action="store_true", dest="no_format")
    parser.add_argument("--no-lint",      action="store_true", dest="no_lint")
    parser.add_argument("--no-headers",   action="store_true", dest="no_headers")
    parser.add_argument("--no-structure", action="store_true", dest="no_structure")
    parser.add_argument("--verbose",      action="store_true",
                        help="Show formatter diffs inline")
    parser.add_argument("--summary",      action="store_true",
                        help="Print only the final summary")
    opts = parser.parse_args()

    restrict = opts.path.resolve() if opts.path else None
    mode = "APPLYING" if opts.apply else "DRY-RUN"
    print(f"{_c('bold', 'PhyloPhere Style Agent')}  [{mode}]")
    print(f"Root: {ROOT}\n")

    all_results: List[FileResult] = []

    if opts.lang in (None, "python"):
        for f in collect_files({".py"}, restrict):
            all_results.append(process_python(f, opts))

    if opts.lang in (None, "nextflow"):
        for f in collect_files({".nf"}, restrict):
            all_results.append(process_nextflow(f, opts))

    if opts.lang in (None, "r"):
        for f in collect_files({".R", ".Rmd"}, restrict):
            all_results.append(process_r(f, opts))

    if opts.summary:
        total = sum(len(r.issues) for r in all_results)
        dirty = sum(1 for r in all_results if r.has_issues)
        print(f"Total: {total} issues across {dirty} files")
        sys.exit(1 if total else 0)
    else:
        sys.exit(print_report(all_results, opts.verbose, opts.apply))


if __name__ == "__main__":
    main()
