#!/usr/bin/env python3
"""Exact foreground annotation of a rooted tree, retaining alignment coordinates."""
import argparse
from pathlib import Path
import dendropy
from Bio import SeqIO
from alignment_io import read_alignment


def parse_species_file(path):
    rows = [line.strip() for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    if not rows or len(rows) != len(set(rows)):
        raise ValueError(f'Empty or duplicate species list: {path}')
    if any(any(c.isspace() for c in sp) or ',' in sp for sp in rows):
        raise ValueError(f'Invalid species identifier: {path}')
    return set(rows)


def label_internal_nodes(tree, fg_species, internal_mode='all_descendants'):
    """Deterministic, collision-free identifiers, identical across directions."""
    used = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    fg_labels, index = set(), 0
    for node in tree.postorder_internal_node_iter():
        if node.parent_node is None:
            node.label = None
            continue
        index += 1
        label = f'IntNode{index}'
        while label in used:
            label = '_' + label
        used.add(label)
        node.label = label
        if internal_mode == 'all_descendants' and all(leaf.taxon.label in fg_species for leaf in node.leaf_iter()):
            fg_labels.add(label)
    return fg_labels


def serialize(tree, foreground, internal_labels, label='Foreground'):
    # Serialize nodes structurally: a label can never match another tip's suffix.
    def quote(name):
        return "'" + name.replace("'", "''") + "'" if any(c in name for c in " ()[],:;'{}\t\n") else name
    def visit(node):
        if node.is_leaf():
            name = node.taxon.label
            text = quote(name)
            tagged = name in foreground
        else:
            text = '(' + ','.join(visit(child) for child in node.child_node_iter()) + ')'
            text += quote(node.label) if node.label else ''
            tagged = node.label in internal_labels
        if tagged and node.parent_node is not None:
            text += '{' + label + '}'
        if node.edge_length is not None and node.parent_node is not None:
            text += ':' + format(node.edge_length, '.17g')
        return text
    return visit(tree.seed_node) + ';\n'


def annotate(tree_path, fasta_path, fg_path, output, *, bg_path=None, fasta_out=None,
             internal_mode='all_descendants', label='Foreground'):
    if internal_mode not in {'all_descendants', 'none'}:
        raise ValueError('Invalid internal-node mode')
    fg = parse_species_file(fg_path)
    bg = parse_species_file(bg_path) if bg_path else None
    if bg is not None and fg & bg:
        raise ValueError('Foreground and background pools overlap')
    tree = dendropy.Tree.get(path=str(tree_path), schema='newick', preserve_underscores=True)
    if tree.is_rooted is False or len(tree.seed_node.child_nodes()) != 2:
        raise ValueError('FADE requires a rooted, bifurcating reference root; no automatic rerooting')
    leaves = list(tree.leaf_node_iter())
    taxa = [node.taxon.label if node.taxon else None for node in leaves]
    if any(not t for t in taxa) or len(taxa) != len(set(taxa)):
        raise ValueError('Empty or duplicate tree tips')
    records = read_alignment(fasta_path)
    allowed = set(taxa) & {r.id for r in records}
    if bg is not None:
        allowed &= fg | bg
    result = dict(input_taxa=len(records), alignment_columns=len(records[0].seq),
                  retained_taxa=sorted(allowed), foreground_taxa=sorted(fg & allowed),
                  excluded_taxa=sorted({r.id for r in records} - allowed), internal_nodes=internal_mode)
    reason = 'fewer_than_two_taxa' if len(allowed) < 2 else 'no_foreground' if not fg & allowed else None
    if reason:
        return dict(result, status='excluded', reason=reason)
    tree.retain_taxa_with_labels(allowed)
    tree.is_rooted = True
    fg_internal = label_internal_nodes(tree, fg, internal_mode)
    Path(output).write_text(serialize(tree, fg & allowed, fg_internal, label))
    if fasta_out:
        SeqIO.write([r for r in records if r.id in allowed], fasta_out, 'fasta')
    return dict(result, status='eligible', foreground_internal_branches=len(fg_internal))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--species-file', required=True)
    parser.add_argument('--bg-species-file')
    parser.add_argument('--fasta_out')
    parser.add_argument('--internal-nodes', default='all_descendants')
    parser.add_argument('--label', default='Foreground')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = annotate(args.tree, args.fasta, args.species_file, args.output,
                      bg_path=args.bg_species_file, fasta_out=args.fasta_out,
                      internal_mode=args.internal_nodes, label=args.label)
    print(result)


if __name__ == '__main__':
    main()
