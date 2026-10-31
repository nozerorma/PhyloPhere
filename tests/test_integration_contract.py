"""Static contracts between main.nf, the workflows, conf/*.config and the help.

No Nextflow run: everything is read from the sources, so a change that breaks one of these
contracts fails here before it fails in a run.
"""
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = [ROOT / 'nextflow.config', *sorted((ROOT / 'conf').glob('*.config'))]
NF_FILES = [ROOT / 'main.nf', *sorted((ROOT / 'workflows').glob('*.nf')), *sorted((ROOT / 'subworkflows').rglob('*.nf'))]

# Parameters of the executor and the container profiles: not pipeline options, so not in the help.
INFRASTRUCTURE = {'max_cpus', 'max_memory', 'max_time', 'slurm_exit_read_timeout', 'slurm_queue', 'slurm_queue_size',
                  'slurm_submit_rate_limit', 'use_apptainer', 'use_singularity'}
# Read through a bare `params.<name>` that is not a parameter: Groovy map methods and a property of the FADE json.
NOT_PARAMETERS = {'containsKey', 'json'}
# Read behind a containsKey guard and declared nowhere on purpose: the former name of caas_map_dir in older project files.
COMPATIBILITY_ALIASES = {'vep_map_dir'}
# Declared and read only outside Nextflow code: the generated run scripts and bin/resolve_core_inputs.py.
READ_OUTSIDE_NEXTFLOW = {'auto_fetch_gmt', 'auto_generate_ensembl', 'ensembl_dataset', 'remove_work'}


def strip_comments(text):
    return re.sub(r'//[^\n]*', '', text)


def declared_params():
    """Names assigned inside a params {} scope of any config file."""
    names = set()
    for path in CONFIGS:
        inside, depth = False, 0
        for line in strip_comments(path.read_text()).splitlines():
            if re.match(r'\s*params\s*\{', line):
                inside, depth = True, 0
            if inside:
                m = re.match(r'\s*([A-Za-z_]\w*)\s*=', line)
                if m and depth <= 1:
                    names.add(m.group(1))
                depth += line.count('{') - line.count('}')
                if depth <= 0 and '}' in line:
                    inside = False
    return names


def read_params():
    """Names read as params.<name> or params['<name>'] in the Nextflow code and the configs."""
    names = set()
    for path in NF_FILES + CONFIGS + sorted((ROOT / 'lib').glob('*.groovy')):
        text = strip_comments(path.read_text())
        names |= set(re.findall(r'\bparams\.([A-Za-z_]\w*)', text))
        names |= set(re.findall(r'params\[[\'"](\w+)[\'"]\]', text))
    return {n for n in names if not n.endswith('_')} - NOT_PARAMETERS - COMPATIBILITY_ALIASES


def fade_dynamic_params():
    """fadeSettings() reads these by name (params[name]), which no grep of params.<name> sees."""
    text = (ROOT / 'subworkflows/FADE/fade_contract.nf').read_text()
    return set(re.findall(r"get\('(fade_\w+)'", text))


def help_flags():
    return set(re.findall(r'(?<![\w-])--([A-Za-z][A-Za-z0-9_]*)', (ROOT / 'workflows/help.nf').read_text()))


def test_every_help_flag_is_declared():
    missing = sorted(help_flags() - declared_params())
    assert not missing, f'documented in the help and declared in no config: {missing}'


def test_every_parameter_read_by_the_code_is_declared():
    missing = sorted(read_params() - declared_params())
    assert not missing, f'read as params.<name> and declared in no config: {missing}'


def test_every_active_parameter_is_documented():
    active = (read_params() | fade_dynamic_params()) & declared_params()
    missing = sorted(active - help_flags() - INFRASTRUCTURE)
    assert not missing, f'read by the code and absent from the help: {missing}'


def test_declared_parameters_have_a_reader():
    readers = read_params() | fade_dynamic_params() | READ_OUTSIDE_NEXTFLOW | INFRASTRUCTURE
    unread = sorted(declared_params() - readers)
    assert not unread, f'declared in a config and read by nothing: {unread}'


def test_no_script_assigns_a_parameter():
    # The first assignment of a params key wins: an assignment in a script is silently ignored.
    offenders = []
    for path in NF_FILES:
        for number, line in enumerate(strip_comments(path.read_text()).splitlines(), 1):
            if re.match(r'\s*params\.[A-Za-z_]\w*\s*=[^=]', line):
                offenders.append(f'{path.relative_to(ROOT)}:{number}')
    assert not offenders, offenders


# ── Calls of main.nf against the take:/input: of what they call ──────────────────────────────────

def parse_includes(text):
    """alias -> (original name, file) for every include of main.nf."""
    includes = {}
    for body, source in re.findall(r'include\s*\{([^}]*)\}\s*from\s*[\'"]([^\'"]+)[\'"]', text):
        for item in body.split(';'):
            item = item.strip()
            if not item:
                continue
            m = re.match(r'(\w+)(?:\s+as\s+(\w+))?$', item)
            name, alias = m.group(1), m.group(2) or m.group(1)
            path = (ROOT / source).with_suffix('.nf') if not source.endswith('.nf') else ROOT / source
            includes[alias] = (name, path)
    return includes


def declared_arity(name, path):
    """Number of take: lines of a workflow or input: declarations of a process; None for a function."""
    text = path.read_text()
    m = re.search(r'^workflow\s+%s\s*\{(.*?)^\}' % name, text, re.S | re.M)
    if m:
        take = re.search(r'take:\s*\n(.*?)^\s*(?:main|emit|work):', m.group(1), re.S | re.M)
        lines = strip_comments(take.group(1)).splitlines() if take else []
        return len([l for l in lines if l.strip()])
    m = re.search(r'^process\s+%s\s*\{(.*?)^\}' % name, text, re.S | re.M)
    if m:
        inputs = re.search(r'input:\s*\n(.*?)^\s*output:', m.group(1), re.S | re.M)
        lines = strip_comments(inputs.group(1)).splitlines() if inputs else []
        return len([l for l in lines if re.match(r'\s*(path|val|tuple|each|env|stdin)\b', l)])
    return None


def call_arities(text, names):
    """name -> list of argument counts of every call of that name in the text."""
    calls = {}
    for m in re.finditer(r'(?<![\w.])(%s)\(' % '|'.join(sorted(names, key=len, reverse=True)), text):
        i, depth, args, current = m.end(), 1, 0, False
        while depth:
            c = text[i]
            if c in '([{':
                depth += 1
            elif c in ')]}':
                depth -= 1
            elif c == ',' and depth == 1:
                args += 1
            if depth and not c.isspace() and not (c == ',' and depth == 1):
                current = True
            i += 1
        calls.setdefault(m.group(1), []).append(args + 1 if current else 0)
    return calls


def test_main_calls_match_the_signatures():
    text = strip_comments((ROOT / 'main.nf').read_text())
    includes = parse_includes(text)
    arities = {alias: declared_arity(name, path) for alias, (name, path) in includes.items()}
    callable_names = {alias for alias, n in arities.items() if n is not None}
    assert {'SCORING', 'ENRICHMENT', 'POSENRICH', 'CAAS_CORE', 'RER_MAIN'} <= callable_names, callable_names
    calls = call_arities(text, callable_names)
    wrong = {alias: (arities[alias], seen) for alias, seen in calls.items() if any(n != arities[alias] for n in seen)}
    assert not wrong, f'(declared, called with): {wrong}'
    assert {'SCORING', 'ENRICHMENT', 'POSENRICH'} <= set(calls)


# ── Resources and retries ────────────────────────────────────────────────────────────────────────

def process_definitions():
    procs = {}
    for path in NF_FILES:
        for m in re.finditer(r'^process\s+(\w+)\s*\{(.*?)^\}', path.read_text(), re.S | re.M):
            body = m.group(2)
            procs[m.group(1)] = {'labels': re.findall(r"^\s*label\s+['\"](\w+)['\"]", body, re.M),
                                 'body': body, 'file': path}
    return procs


def selector_blocks():
    """(kind, pattern, attributes) of every withLabel/withName block of resources.config."""
    text = (ROOT / 'conf/resources.config').read_text()
    blocks = []
    for m in re.finditer(r"with(Label|Name)\s*:\s*['\"]?([^'\"{\s]+(?:\s*\|\s*[^'\"{\s]+)*)['\"]?\s*\{", text):
        depth, j = 1, m.end()
        while depth:
            depth += (text[j] == '{') - (text[j] == '}')
            j += 1
        attrs = {a.group(1): a.group(2).split('//')[0].strip()
                 for a in re.finditer(r'^\s*(cpus|memory|time|errorStrategy|maxRetries)\s*=\s*(.+)$', text[m.end():j - 1], re.M)}
        blocks.append((m.group(1), m.group(2), attrs))
    return blocks


def matches(pattern, name):
    return any(re.fullmatch(p.strip(), name) for p in pattern.split('|'))


def effective(name, procs, blocks):
    attrs = {}
    labels = {pat: a for kind, pat, a in blocks if kind == 'Label'}
    for label in procs[name]['labels']:
        attrs.update(labels.get(label, {}))
    for kind, pat, a in blocks:
        if kind == 'Name' and matches(pat, name):
            attrs.update(a)
    return attrs


def test_every_process_has_resources():
    procs, blocks = process_definitions(), selector_blocks()
    bare = [n for n, p in procs.items()
            if not effective(n, procs, blocks).keys() & {'cpus', 'memory', 'time'}]
    assert not bare, f'processes with no label or withName request: {sorted(bare)}'


def test_resource_selectors_match_something():
    procs, blocks = process_definitions(), selector_blocks()
    used_labels = {l for p in procs.values() for l in p['labels']}
    # include aliases of main.nf name processes that exist once in the code
    aliases = set(parse_includes(strip_comments((ROOT / 'main.nf').read_text())))
    dead = []
    for kind, pat, _ in blocks:
        if kind == 'Label' and pat not in used_labels:
            dead.append(f'label {pat}')
        if kind == 'Name' and not any(matches(pat, n) for n in set(procs) | aliases):
            dead.append(f'name {pat}')
    assert not dead, dead


def test_attempt_scaling_has_a_retry():
    procs, blocks = process_definitions(), selector_blocks()
    inert = []
    for name in sorted(procs):
        attrs = effective(name, procs, blocks)
        scales = any('task.attempt' in attrs.get(k, '') for k in ('cpus', 'memory', 'time'))
        strategy = attrs.get('errorStrategy')
        own = re.search(r'^\s*errorStrategy\s+([^\n]+)', procs[name]['body'], re.M)
        retries = ('retry' in strategy) if strategy else bool(own and 'retry' in own.group(1))
        if scales and not retries:
            inert.append(name)
    assert not inert, f'memory or time scales with task.attempt but nothing retries: {inert}'
