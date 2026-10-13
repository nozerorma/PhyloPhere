#!/usr/bin/env python3
# test_chunk_count.py: the number of replay chunks of a gene follows the memory budget of a worker.
# PhyloPhere | tests/ct_replay_chunking/
"""replay_chunk_count: no budget or no rows is one chunk; the count is the smallest that fits the usable memory of a
worker, grows with the export size and falls when the budget grows (a retry with more memory)."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "subworkflows" / "CT_DISAMBIGUATION" / "local"))
from src.utils.gene_wrapper import replay_chunk_count, _REPLAY_BASE_MB, _REPLAY_MB_PER_MROW, _REPLAY_MEM_HEADROOM  # noqa: E402


def main():
    assert replay_chunk_count(5_000_000, 0) == 1            # no budget: chunking off
    assert replay_chunk_count(0, 5000) == 1                 # no rows
    # a gene that fits is not split; one that does not is split into the smallest number of chunks that fits
    for rows in (0.5e6, 1.5e6, 4.65e6, 6.4e6):
        for budget in (2500, 5000, 10000, 20000):
            k = replay_chunk_count(rows, budget)
            usable = max(_REPLAY_MEM_HEADROOM * budget - _REPLAY_BASE_MB, 0.2 * budget)
            assert rows / 1e6 * _REPLAY_MB_PER_MROW / k <= usable + 1e-6, (rows, budget, k)
            if k > 1:
                assert rows / 1e6 * _REPLAY_MB_PER_MROW / (k - 1) > usable, (rows, budget, k)
    # monotone: more rows, more chunks; more budget, fewer chunks
    ks = [replay_chunk_count(r, 5000) for r in (1e6, 2e6, 4e6, 8e6)]
    assert ks == sorted(ks) and ks[0] < ks[-1], ks
    kb = [replay_chunk_count(6.4e6, b) for b in (5000, 10000, 20000, 40000)]
    assert kb == sorted(kb, reverse=True) and kb[-1] == 1, kb
    print("test_chunk_count: OK", ks, kb)


if __name__ == "__main__":
    main()
