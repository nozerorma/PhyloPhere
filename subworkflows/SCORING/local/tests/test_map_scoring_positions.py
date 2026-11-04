import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SCRIPT = ROOT / 'subworkflows/SCORING/local/src/map_scoring_positions.py'
VEP = ROOT / 'subworkflows/VEP/local/src'


def run(tmp_path, genes, maps):
    mdir = tmp_path / 'maps'
    mdir.mkdir(parents=True, exist_ok=True)
    for name in maps:
        (mdir / f'{name}.map.tsv').write_text(
            'status\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n'
            'selected\t1\t1\tchr1:1000000\nselected\t2\t2\tchr1:1000001\n')
    scores = tmp_path / 'scores.tsv'
    scores.write_text('Gene\tPosition\tside\n' + ''.join(f'{g}\t{p}\ttop\n' for g, p in genes))
    out = tmp_path / 'mapped.tsv'
    result = subprocess.run([sys.executable, str(SCRIPT), '--scores', str(scores), '--map-dir', str(mdir),
                             '--vep-src', str(VEP), '--output', str(out)], capture_output=True, text=True)
    return result, out


def test_species_suffix_literal_ids_and_missing(tmp_path):
    result, out = run(tmp_path, [('001.Homo_sapiens', 0), ('001.Homo_sapiens', 1), ('missing', 0)], ['001'])
    assert result.returncode == 0, result.stderr
    with out.open() as fh:
        rows = {(r['Gene'], r['Position']): r for r in csv.DictReader(fh, delimiter='\t')}
    assert rows['001.Homo_sapiens', '0']['hg38_nt_coord'] == 'chr1:1000000'
    assert rows['001.Homo_sapiens', '1']['hg38_nt_coord'] == 'chr1:1000001'
    assert rows['001.Homo_sapiens', '0']['map_gene'] == '001'
    assert rows['missing', '0']['mapping_status'] == 'no_map'


def test_ambiguous_alias_and_version_mismatch(tmp_path):
    result, _ = run(tmp_path, [('X', 0)], ['X.Homo_sapiens', 'X.Papio_anubis'])
    assert result.returncode != 0 and 'Ambiguous producer gene alias' in result.stderr
    result, out = run(tmp_path / 'version', [('X.3.Homo_sapiens', 0)], ['X.2.Homo_sapiens'])
    assert result.returncode == 0, result.stderr
    assert '\tno_map\t' in out.read_text()


def test_two_scored_genes_cannot_share_one_map(tmp_path):
    result, _ = run(tmp_path, [('X', 0), ('X.Homo_sapiens', 1)], ['X'])
    assert result.returncode != 0 and 'Ambiguous producer gene alias' in result.stderr


def test_shared_vep_map_lookup_accepts_producer_alias(tmp_path):
    sys.path.insert(0, str(VEP))
    import vep_common as vr
    mdir = tmp_path / 'maps'; mdir.mkdir()
    (mdir / 'X.1.map.tsv').write_text('')
    assert vr.map_path('X.1.Homo_sapiens', mdir).name == 'X.1.map.tsv'
    assert vr.map_path('X.2.Homo_sapiens', mdir) is None
    assert vr.map_path('X.Homo_sapiens', mdir) is None
    assert vr.map_path('X', mdir) is None
    assert vr.resolve_gene('Y.Homo_sapiens', {'Y.Papio_anubis'}) == 'Y.Homo_sapiens'
