#!/usr/bin/env bash
# union_shards.sh DEST SRC... : unite per-gene shard directories into DEST.
# Files are hard-linked when the first shard can be linked from here (same filesystem) and copied otherwise; DEST ends
# up holding every source's files, b0/ subdirectories merged. A gene that several sources hold (one shard per cycle
# block) is joined into one shard by merge_block_shards.py, run with UNION_PYTHON (default python3).
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

# Shard paths (relative to their source) found in more than one source.
shared="${dest}.shared"
find -L "$@" -type f -name '*.tsv.gz' -printf '%P\n' | sort | uniq -d > "$shared"

# --remove-destination: a shard held by several sources is replaced by the last one, then joined below.
for d in "$@"; do
    cp -a${link}L --remove-destination "$d"/. "$dest"/
done

if [[ -s "$shared" ]]; then
    # UNION_PYTHON may be a command prefix (container entrypoint plus interpreter): left unquoted to split.
    ${UNION_PYTHON:-python3} "$(dirname "${BASH_SOURCE[0]}")/merge_block_shards.py" \
        --dest "$dest" --shared "$shared" "$@"
fi
rm -f "$shared"
