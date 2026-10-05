#!/usr/bin/env python3
"""Cross-check of every surface a pipeline parameter lives on.

A parameter is declared in conf/*.config, read by the pipeline code, written by the GUI into params.json, held by a field
of the GUI model, shown by a tab field and rendered by a template. Each link can break on its own, and a broken link is
silent: a tab field whose value no template renders does nothing, a declared parameter nobody reads is dead weight, a
parameter the GUI never writes cannot be set from the GUI.

`audit(root)` returns the findings as (code, name, detail) tuples:

  DECLARED_UNREAD      declared in conf/, read by no pipeline code
  READ_UNDECLARED      read as params.X, declared nowhere (no default: the read is null unless the caller sets it)
  GUI_UNDECLARED       written by the GUI, declared nowhere
  GUI_UNREAD           written by the GUI, read by no pipeline code
  NOT_IN_GUI           declared and read, never written by the GUI
  FIELD_NOT_RENDERED   a model field no template or context code uses
  FIELD_NOT_IN_TAB     a field of a module model with no tab field
  TAB_FIELD_NOT_IN_MODEL  a tab field naming nothing in its module model
  DEFAULT_MISMATCH     the GUI model default differs from the conf default of the parameter it feeds
  BLANK_BURNS_DEFAULT  the GUI writes "" for a parameter whose conf default is not empty, so a blank field replaces the default
  DUPLICATE_FIELD      one name held by two module models
  TEMPLATE_STALE_KEY   a key of a GUI template (gui/templates/*.json) that the model no longer has; the loader ignores it silently

Run as a script it prints the findings and exits 1 when any is not in the allowlist (a TSV of code, name, reason).
"""
import argparse
import ast
import collections
import dataclasses
import os
import re
import subprocess
import sys
from pathlib import Path

PROFILES = (None, "slurm", "local")
# the module models and the tab file that shows each
MODULE_TABS = {
    "modules.caas": "caas_tab.py", "modules.disambiguation": "disambiguation_tab.py", "modules.accumulation": "accumulation_tab.py",
    "modules.rer": "rerconverge_tab.py", "modules.fade": "fade_tab.py", "modules.scoring": "scoring_tab.py",
    "modules.vep": "vep_tab.py", "modules.enrichment": "enrichment_tab.py",
}
# a module model's own bookkeeping fields, shown by the shared tab frame and not by a field spec
FRAME_FIELDS = {"enabled", "extra_flags"}


def _tracked(root):
    out = subprocess.run(["git", "-C", str(root), "ls-files"], capture_output=True, text=True)
    if out.returncode == 0 and out.stdout.strip():
        return [p for p in out.stdout.split("\n") if p]
    return [str(p.relative_to(root)) for p in Path(root).rglob("*") if p.is_file()]


def declared_params(root):
    """{name: {values}} of params.* resolved by Nextflow, over the default, slurm and local profiles."""
    declared = collections.defaultdict(set)
    for profile in PROFILES:
        cmd = ["nextflow", "config", "-flat", str(root)] + (["-profile", profile] if profile else [])
        out = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, NXF_OFFLINE="true"))
        if out.returncode != 0:
            raise RuntimeError(f"nextflow config failed: {out.stderr[-500:]}")
        for line in out.stdout.splitlines():
            m = re.match(r"params\.([A-Za-z_]\w*)\s*=\s*(.*)$", line)
            if m:
                declared[m.group(1)].add(m.group(2).strip())
    return declared


def _strip_comments(text):
    """The code without // line comments (a // inside a string, as in a URL, has no space before it) and /* */ blocks."""
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(^|\s)//.*$", r"\1", text, flags=re.M)


def read_params(root, files):
    """{name: {files}} of the parameters pipeline code reads: params.X, params['X'] and params.containsKey('X')."""
    reads = collections.defaultdict(set)
    code = [f for f in files if (f.endswith((".nf", ".groovy")) or f == "nextflow.config" or f.startswith(("conf/", "lib/")))
            and not f.startswith(("validation/", "gui/", "style/", "archive/"))]
    for f in code:
        text = _strip_comments((Path(root) / f).read_text(errors="ignore"))
        names = re.findall(r"\bparams\.([A-Za-z_]\w*)\b(?!\s*\()", text)
        names += re.findall(r"\bparams\[\s*['\"]([A-Za-z_]\w*)['\"]\s*\]", text)
        names += re.findall(r"\bparams\.containsKey\(\s*['\"]([A-Za-z_]\w*)['\"]", text)
        for n in names:
            reads[n].add(f)
    return reads


def gui_params(root):
    """{param: (variable, expression)} of the params.json the single-run template writes; variable is the shell variable
    the value is read from (None when it is a template expression), expression the value text."""
    text = (Path(root) / "gui/generation/templates/run_single.sh.j2").read_text()
    start = text.index("<<PARAMS_EOF")
    block = text[start:text.index("PARAMS_EOF", start + 20)]
    out = {}
    for m in re.finditer(r'^\s*(?:\{%.*?%\})*\s*"([A-Za-z_]\w*)"\s*:\s*(.*?),?\s*(?:\{%.*?%\})*\s*$', block, re.M):
        expr = m.group(2)
        v = re.search(r"\$\{?(\w+)", expr)
        out[m.group(1)] = (v.group(1) if v and "{{" not in expr else None, expr)
    return out


def exports(root):
    """{shell variable: Jinja expression} of the batch template's `export VAR="{{ expr }}"` lines."""
    text = (Path(root) / "gui/generation/templates/sbatch_array.sh.j2").read_text()
    return {m.group(1): m.group(2).strip() for m in re.finditer(r'^\s*export\s+(\w+)="\{\{\s*(.*?)\s*\}\}"', text, re.M)}


def _project_config(root):
    """The ProjectConfig class of the tree at `root`: the gui package is imported afresh from there and not left in sys.modules."""
    def purge():
        for name in [m for m in sys.modules if m == "gui" or m.startswith("gui.")]:
            del sys.modules[name]

    purge()
    sys.path.insert(0, str(root))
    try:
        from gui.models.project import ProjectConfig
        return ProjectConfig
    finally:
        sys.path.remove(str(root))
        purge()


def model_fields(root):
    """[(container, field, default)] of every leaf field of the GUI project model."""
    ProjectConfig = _project_config(root)

    def walk(obj, path):
        for f in dataclasses.fields(obj):
            v = getattr(obj, f.name)
            if dataclasses.is_dataclass(v):
                yield from walk(v, path + [f.name])
            else:
                yield ".".join(path), f.name, v

    return list(walk(ProjectConfig(), []))


def tab_fields(root):
    """{container: {field names}} from the FieldSpec(name=...) calls of each module tab."""
    out = {}
    for container, fname in MODULE_TABS.items():
        tree = ast.parse((Path(root) / "gui/widgets/tabs" / fname).read_text())
        names = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and getattr(node.func, "id", "") == "FieldSpec":
                for kw in node.keywords:
                    if kw.arg == "name" and isinstance(kw.value, ast.Constant):
                        names.add(kw.value.value)
        # tabs with custom widgets (check boxes grouped into one control) bind the model field through the config object
        names |= set(re.findall(r"\bconfig\.(\w+)\b", (Path(root) / "gui/widgets/tabs" / fname).read_text()))
        out[container] = names
    return out


_ALIASES = {
    "modules.caas": "caas", "modules.disambiguation": "disambig|disambiguation", "modules.scoring": "scoring",
    "modules.enrichment": "enrichment", "modules.rer": "rer", "modules.fade": "fade", "modules.vep": "vep",
    "modules.accumulation": "accum|accumulation", "precomputed": "pc|precomputed", "runtime": "runtime", "general": "general",
}


def rendered(root, container, field):
    """True when a template or the context code uses container.field (an alias of the container counts in the context code)."""
    gen = Path(root) / "gui/generation"
    templates = "\n".join(p.read_text() for p in (gen / "templates").glob("*.j2"))
    if re.search(r"\b%s\.%s\b" % (re.escape(container), re.escape(field)), templates):
        return True
    context = (gen / "context.py").read_text() + (gen / "render.py").read_text()
    alias = _ALIASES.get(container, re.escape(container))
    full = r"project\.%s" % re.escape(container)
    return re.search(r"\b(?:%s|%s)\.%s\b|getattr\(\s*(?:%s|%s)\s*,\s*['\"]%s['\"]" % (alias, full, re.escape(field), alias, full, re.escape(field)), context) is not None


def stale_template_keys(root):
    """['file: path.to.key'] for every key of a GUI template that is not a field of the project model."""
    import json
    ProjectConfig = _project_config(root)

    def check(obj, data, path, out):
        names = {f.name: getattr(obj, f.name) for f in dataclasses.fields(obj)}
        for k, v in data.items():
            if k not in names:
                out.append(".".join(path + [k]))
            elif dataclasses.is_dataclass(names[k]) and isinstance(v, dict):
                check(names[k], v, path + [k], out)

    found = []
    for f in sorted(Path(root).glob("gui/templates/*.json")) + sorted(Path(root).glob("validation/tier1/templates/*.json")):
        keys = []
        check(ProjectConfig(), json.loads(f.read_text()), [], keys)
        found += [f"{f.name}: {k}" for k in keys]
    return found


def _norm(value):
    s = str(value).strip().strip("'\"")
    if s.lower() in ("true", "false"):
        return s.lower()
    try:
        return repr(float(s))
    except ValueError:
        return s


def _mapped_field(expr, var, export_map):
    """(container, field, plain) the GUI value of a parameter comes from; plain is False when a filter or condition acts on it."""
    src = export_map.get(var, "") if var else expr
    m = re.search(r"((?:modules\.\w+)|runtime|general|precomputed|resources)\.(\w+)", src)
    if not m:
        return None
    plain = re.fullmatch(r"\s*%s\.%s(\s*\|\s*default\([^)]*\))?\s*" % (re.escape(m.group(1)), re.escape(m.group(2))), src) is not None
    return m.group(1), m.group(2), plain


def audit(root):
    root = Path(root)
    files = _tracked(root)
    declared, reads, gui = declared_params(root), read_params(root, files), gui_params(root)
    export_map, fields, tabs = exports(root), model_fields(root), tab_fields(root)
    findings = []
    add = lambda code, name, detail="": findings.append((code, name, detail))
    D, R, G = set(declared), set(reads), set(gui)
    for n in sorted(D - R):
        add("DECLARED_UNREAD", n, f"declared {sorted(declared[n])[0]}")
    for n in sorted(R - D):
        add("READ_UNDECLARED", n, ", ".join(sorted(reads[n])[:2]))
    for n in sorted(G - D):
        add("GUI_UNDECLARED", n)
    for n in sorted((G & D) - R):
        add("GUI_UNREAD", n)
    for n in sorted((D & R) - G):
        add("NOT_IN_GUI", n, ", ".join(sorted(reads[n])[:2]))
    by_name = collections.defaultdict(list)
    for container, field, _ in fields:
        if container.startswith("modules."):
            by_name[field].append(container)
    for name, cs in sorted(by_name.items()):
        if len(cs) > 1 and name not in FRAME_FIELDS:
            add("DUPLICATE_FIELD", name, " and ".join(cs))
    for container, field, default in fields:
        if container and not rendered(root, container, field) and field not in ("schema_version", "flavor"):
            add("FIELD_NOT_RENDERED", f"{container}.{field}")
        if container in tabs and field not in tabs[container] and field not in FRAME_FIELDS:
            add("FIELD_NOT_IN_TAB", f"{container}.{field}")
    model_names = {c: {f for cc, f, _ in fields if cc == c} for c in tabs}
    for container, names in tabs.items():
        for n in sorted(names - model_names[container]):
            add("TAB_FIELD_NOT_IN_MODEL", f"{container}.{n}")
    for k in stale_template_keys(root):
        add("TEMPLATE_STALE_KEY", k)
    defaults = {(c, f): d for c, f, d in fields}
    for param, (var, expr) in sorted(gui.items()):
        mapped = _mapped_field(expr, var, export_map)
        if not mapped or not mapped[2] or param not in declared:
            continue
        container, field, _ = mapped
        if (container, field) not in defaults or defaults[(container, field)] == "":
            continue   # a blank model field defers to the conf default; BLANK_BURNS_DEFAULT checks that deferral
        conf_values = {_norm(v) for v in declared[param]}
        if any("$" in v or "projectDir" in v or "baseDir" in v for v in conf_values):
            continue
        if _norm(defaults[(container, field)]) not in conf_values:
            add("DEFAULT_MISMATCH", param, f"GUI {container}.{field} = {defaults[(container, field)]!r}, conf {sorted(declared[param])}")
    for param, (var, expr) in sorted(gui.items()):
        if param in declared and var and " else " not in export_map.get(var, "") and re.fullmatch(r'"\$\{%s:-\}"' % var, expr.strip().rstrip(",")):
            conf_values = {_norm(v) for v in declared[param]} - {"", "None", "null"}
            if conf_values:
                add("BLANK_BURNS_DEFAULT", param, f"template fallback is empty, conf {sorted(declared[param])}")
    return findings


def load_allowlist(path):
    allowed = {}
    if Path(path).exists():
        for line in Path(path).read_text().splitlines():
            if line.strip() and not line.startswith("#"):
                code, name, reason = (line.split("\t") + ["", ""])[:3]
                allowed[(code, name)] = reason
    return allowed


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[2]))
    ap.add_argument("--allowlist", default=str(Path(__file__).with_name("param_audit_allow.tsv")))
    args = ap.parse_args()
    findings = audit(args.root)
    allowed = load_allowlist(args.allowlist)
    open_ = [f for f in findings if (f[0], f[1]) not in allowed]
    for code, name, detail in findings:
        print(f"{'ok ' if (code, name) in allowed else 'OPEN'}\t{code}\t{name}\t{detail}")
    print(f"\n{len(findings)} findings, {len(open_)} not in the allowlist", file=sys.stderr)
    sys.exit(1 if open_ else 0)


if __name__ == "__main__":
    main()
