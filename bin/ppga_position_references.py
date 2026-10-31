"""Explicit projections for the verified PPGA productive-alignment source contract.

PPGA UCR positions index its entropy input (trimmed MSA, 1-based). Pfam
endpoints index the ungapped human sequence of that same MSA. Original human
residue numbers in MAP are a separate coordinate system. This converter does
not infer the source space from sequence lengths: project-ppga explicitly
selects that producer contract, and validates the supplied bundle against it.
"""
import csv
import gzip
import hashlib
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

from map_domain_variability import read_fasta, compute_domain_stats


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def index_files(directory, suffixes):
    result = {}
    for p in sorted(directory.iterdir()):
        if not p.is_file() or not any(p.name.endswith(s) for s in suffixes):
            continue
        gene = next(p.name[:-len(s)] for s in sorted(suffixes,key=len,reverse=True) if p.name.endswith(s))
        if gene in result:
            raise ValueError(f'{gene}: ambiguous files in {directory}')
        result[gene] = p
    return result


class References:
    def __init__(self, base):
        self.base = base
        self.alignments = index_files(base/'PROT/bmge', ('.fa', '.fasta', '.faa', '.fas'))
        self.maps = index_files(base/'MAP', ('.map.tsv',))
        self.entropy = index_files(base/'PROT_VAR_RAW', ('.prot.entropy.tsv', '.entropy.tsv'))
        self.validated = set()

    @lru_cache(maxsize=8)
    def gene(self, gene):
        if gene not in self.alignments or gene not in self.maps:
            raise ValueError(f'{gene}: missing productive alignment or MAP')
        seqs = read_fasta(self.alignments[gene])
        widths = {len(seq) for _, seq in seqs}
        if len(widths) != 1:
            raise ValueError(f'{gene}: unequal productive alignment lengths')
        width = widths.pop()
        refs = [s for h, s in seqs if 'homo_sapiens' in h.lower() or 'hg38' in h.lower()]
        if len(refs) > 1:
            raise ValueError(f'{gene}: ambiguous human reference')
        with self.maps[gene].open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        selected = {}
        for row in rows:
            if row['status'] != 'selected':
                continue
            col = int(row['prot_ali_col']) - 1
            if col in selected:
                raise ValueError(f'{gene}: duplicate MAP column')
            selected[col] = row
        if sorted(selected) != list(range(width)):
            raise ValueError(f'{gene}: MAP does not cover every productive alignment column')
        if gene not in self.entropy:
            raise ValueError(f'{gene}: missing producer entropy table')
        entropy = {}
        with self.entropy[gene].open() as f:
            for row in csv.DictReader(f, delimiter='\t'):
                pos = int(row['position'])
                if pos in entropy:
                    raise ValueError(f'{gene}: duplicate entropy position')
                entropy[pos] = {k: float(row[k]) for k in ('variability', 'g')}
        if sorted(entropy) != list(range(1, width + 1)):
            raise ValueError(f'{gene}: entropy does not index the productive alignment')
        self.validated.add(gene)
        return selected, refs[0] if refs else None, entropy


def project_ppga(base, output):
    base, output = Path(base), Path(output)
    if output.exists():
        raise ValueError('output directory already exists; choose a new reference version')
    refs = References(base)
    sources = [('ucr', base/'reports/variability/ucr_positions.tsv'),
               ('fubar', base/'reports/positive_selection/fubar_sites.tsv')]
    domains = base/'reports/variability/domain_variability.tsv'
    domtbl = base/'HMMSCAN/hmmscan.domtblout'
    if domains.exists():
        sources.insert(0, ('pfam', domains))
        if not domtbl.is_file():
            raise ValueError('Pfam projection requires the preserved hmmscan domtblout')
    hits = set()
    if domains.exists():
        with domtbl.open() as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                v = line.split()
                if len(v) < 22:
                    raise ValueError('Malformed hmmscan domtblout')
                hits.add((v[3], v[1].split('.')[0], int(v[17]), int(v[18]), int(v[5])))
    for _, source in sources:
        if not source.is_file():
            raise ValueError(f'missing producer table: {source}')
    output.mkdir(parents=True)
    completed = []
    for kind, source in sources:
        counts = Counter()
        dest = output/({'pfam': 'domain_variability.tsv.gz', 'ucr': 'ucr_positions.tsv.gz',
                        'fubar': 'fubar_sites.tsv.gz'}[kind])
        partial = Path(str(dest)+'.partial')
        with source.open() as f, gzip.open(partial, 'wt', newline='') as out:
            reader = csv.DictReader(f, delimiter='\t')
            if 'coordinate_space' in (reader.fieldnames or []):
                raise ValueError('project-ppga expects the original producer table, not an already projected table')
            extra = {'pfam': ['alignment_positions'], 'ucr': ['source_position'],
                     'fubar': ['alignment_position']}[kind]
            writer = csv.DictWriter(out, fieldnames=reader.fieldnames+extra+['coordinate_space'],
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for row in reader:
                counts['source_rows'] += 1
                if kind == 'ucr' and row['method'] != 'sliding':
                    counts['other_methods_not_used_by_collections'] += 1
                    continue
                gene = row['gene']
                sites, reference, entropy = refs.gene(gene)
                if kind == 'ucr':
                    col = int(row['position']) - 1
                    if col not in sites:
                        raise ValueError(f'{gene}: UCR position outside productive alignment: {row["position"]}')
                    row['source_position'], row['position'] = row['position'], str(col)
                elif kind == 'pfam':
                    if reference is None:
                        raise ValueError(f'{gene}: Pfam requires an unambiguous human query reference')
                    # Original PPGA hmmscan removes '-' only. X consumes a query position.
                    columns = [i for i, aa in enumerate(reference) if aa != '-']
                    lo, hi = int(row['ali_start']), int(row['ali_end'])
                    hit = (row['gene'], row['pfam_id'].split('.')[0], lo, hi, len(columns))
                    if hit not in hits:
                        raise ValueError(f'{gene}: Pfam hit/query length does not match the preserved scan: {hit}')
                    if not 1 <= lo <= hi <= len(columns):
                        raise ValueError(f'{gene}: Pfam endpoints outside productive query')
                    cols = columns[lo-1:hi]
                    row['alignment_positions'] = ','.join(map(str, cols))
                    for key, value in compute_domain_stats(cols, entropy).items():
                        if key in row:
                            row[key] = value
                else:
                    col = int(row['site']) - 1
                    if col not in sites:
                        raise ValueError(f'{gene}: FUBAR site outside MAP {col}')
                    for key in ('ori_codon_col', 'hg38_aa_pos', 'hg38_nt_coord'):
                        if row.get(key, 'NA') != sites[col][key]:
                            raise ValueError(f'{gene}: FUBAR {key} mismatch at {col}')
                    row['alignment_position'] = str(col)
                row['coordinate_space'] = 'CAAS_0based'
                writer.writerow(row)
                counts['written'] += 1
        partial.rename(dest)
        completed.append((kind, source, dest, dict(counts)))
        print(kind, dict(counts), flush=True)
    # Bind references to the full inventory, including files without annotations.
    alignments = {p.name: sha(p) for p in refs.alignments.values()}
    maps = {p.name: sha(p) for p in refs.maps.values()}
    for kind, source, dest, counts in completed:
        inputs = [source] + ([domtbl] if kind == 'pfam' else [])
        meta = dict(status='validated', coordinate_space='CAAS_0based',
                    projection_version='ppga_productive_alignment_v2', counts=counts,
                    producer_coordinate_contract=('productive_ungapped_human_query_1based' if kind == 'pfam'
                       else 'productive_alignment_1based' if kind == 'ucr' else 'FUBAR_site_1based'),
                    sources={str(p.resolve()): sha(p) for p in inputs}, table_sha256=sha(dest),
                    alignment_dir=str((base/'PROT/bmge').resolve()), alignments=alignments,
                    map_dir=str((base/'MAP').resolve()), maps=maps,
                    interpretation='Coordinate repair from original PPGA annotations; no new CAAS fit',
                    pfam_query_policy='original PPGA removes dash only; X consumes a position; reference gaps excluded from membership')
        Path(str(dest)+'.manifest.json').write_text(json.dumps(meta, indent=2)+'\n')
    manifest = dict(status='complete', projection_version='ppga_productive_alignment_v2',
                    base=str(base.resolve()), validated_genes=len(refs.validated),
                    outputs={dest.name: {'sha256': sha(dest), 'counts': counts}
                             for _, _, dest, counts in completed})
    (output/'reference_bundle_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest
