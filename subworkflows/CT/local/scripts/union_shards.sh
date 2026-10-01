#!/usr/bin/env bash
# union_shards.sh DEST SRC... : unite per-gene shard directories into DEST.
# Genes are disjoint across the sources. Files are hard-linked when the first shard can be linked from here
# (same filesystem) and copied otherwise; DEST ends up holding every source's files, b0/ subdirectories merged.
set -euo pipefail

dest="$1"
shift

link=l
# -quit stops find at the first match: a pipe into `head` would kill find with SIGPIPE (status 141) when
# there are more shards than a pipe holds.
probe="$(find -L "$@" -type f -name '*.tsv.gz' -print -quit)"
if [[ -n "$probe" ]]; then
    if ! ln -L "$probe" "${dest}.link_probe" 2>/dev/null; then link=""; fi
    rm -f "${dest}.link_probe"
fi
mkdir -p "$dest"
for d in "$@"; do
    cp -a${link}L "$d"/. "$dest"/
done
