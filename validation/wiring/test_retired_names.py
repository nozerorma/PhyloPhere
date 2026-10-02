"""No tracked file outside a short allowlist names a process or parameter of the former permulation-null layout."""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
RETIRED = re.compile(r"PERM_REPLAY_BATCHED|\bPERM_REPLAY\b|CAAS_PERMS_(?:DISAMBIGUATE|AGGREGATE|REBUILD|MERGE_DETAIL)|process_perm_replay"
                     r"|ct_perm_replay_batch_size|ct_disambig_perms_batch_size|CT_PERM_REPLAY_BATCH_SIZE|CT_DISAMBIG_PERMS_BATCH_SIZE"
                     r"|caas_b0_diagnostic")
# history and frozen artifacts: the archive, measurement records, generated run scripts and templates under validation/,
# the style archetypes, and the migration that names the retired project fields
ALLOWED = ("archive/", "validation/", "style/", "docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md", "gui/models/serialization.py")


def test_no_tracked_file_outside_the_allowlist_names_the_retired_layout():
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split("\0")
    hits = []
    for name in filter(None, files):
        if name.startswith(ALLOWED):
            continue
        try:
            text = (ROOT / name).read_text()
        except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
            continue
        hits += [f"{name}:{n}: {line.strip()[:100]}" for n, line in enumerate(text.splitlines(), 1) if RETIRED.search(line)]
    assert not hits, "\n".join(hits)
