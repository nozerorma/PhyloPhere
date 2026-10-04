"""permulations.R with no cycle to harvest (N = 0): the FOP mirror loop must not run on an empty pool."""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "subworkflows/CT/local/scripts/permulations.R"

pytestmark = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


def _batch_starts(n, size):
    """The script's own batch_starts(), extracted from its source and evaluated in R."""
    m = re.search(r"^\s*(batch_starts <- function.*)$", SCRIPT.read_text(), re.M)
    assert m, "permulations.R has no batch_starts()"
    out = subprocess.run(["Rscript", "-e", f'{m.group(1)}\ncat(batch_starts({n}L, {size}L), sep=",")'], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return [int(x) for x in out.stdout.split(",") if x]


@pytest.mark.parametrize("n,size,expected", [(0, 1000, []), (1, 1000, [1]), (1000, 1000, [1]), (1001, 1000, [1, 1001]), (2500, 1000, [1, 1001, 2001])])
def test_the_fop_batches_start_every_batch_size_cycles_and_an_empty_pool_has_none(n, size, expected):
    assert _batch_starts(n, size) == expected


def test_the_fop_mirror_loop_iterates_over_batch_starts_not_over_a_bare_seq():
    text = SCRIPT.read_text()
    assert "seq(1L, length(pool), by = FOP_BATCH)" not in text and "in batch_starts(length(pool), FOP_BATCH)" in text
