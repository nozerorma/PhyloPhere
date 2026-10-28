"""Legacy Pfam migration requires matching original scan and productive reference."""
import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'bin'))
sys.path.insert(0, str(ROOT / 'subworkflows/POSENRICH/local/src'))
from migrate_pfam_coordinates import migrate
from position_sources import collect_pfam_members


def fixture(tmp_path):
    domains = tmp_path / 'legacy.tsv'
    domains.write_text('gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\n'
                       'G1\tPF00001\tDomain\tdesc\tCL1\tClan\t1\t2\n')
    query = tmp_path / 'reference_seqs.fa'
    query.write_text('>G1\nAS\n')
    domtbl = tmp_path / 'hmmscan.domtblout'
    parts = ['x'] * 22
    for i, value in {1: 'PF00001.1', 3: 'G1', 17: '1', 18: '2'}.items():
        parts[i] = value
    domtbl.write_text(' '.join(parts) + '\n')
    align = tmp_path / 'alignments'
    align.mkdir()
    (align / 'G1.fa').write_text('>Homo_sapiens\nAX-S\n')
    output = tmp_path / 'converted.tsv'
    return domains, query, domtbl, align, output


def test_verified_migration_preserves_gaps_and_original(tmp_path):
    args = fixture(tmp_path)
    before = args[0].read_bytes()
    # Ambiguous input must fail even before a MAP is available.
    with pytest.raises(ValueError, match='has no coordinate_space'):
        collect_pfam_members(str(args[0]), {}, {'G1'})
    manifest = migrate(*args)
    row = next(csv.DictReader(args[-1].open(), delimiter='\t'))
    assert row['alignment_positions'] == '0,3'
    assert row['coordinate_space'] == 'CAAS_0based'
    assert args[0].read_bytes() == before
    assert manifest['n_hits'] == 1 and manifest['references'][0]['alignment_sha256']
    terms, clans = collect_pfam_members(str(args[-1]), {}, {'G1'})
    assert terms['PF00001'][1] == ['G1:0', 'G1:3']
    assert [v[1] for v in clans.values()] == [['G1:0', 'G1:3']]
    with pytest.raises(ValueError, match='already exists'):
        migrate(*args)


@pytest.mark.parametrize('bad', ['sequence', 'hit', 'bounds', 'ambiguous_reference', 'ambiguous_alignment'])
def test_unverified_migration_creates_no_output(tmp_path, bad):
    args = fixture(tmp_path)
    if bad == 'sequence':
        args[1].write_text('>G1\nAT\n')
    elif bad == 'hit':
        args[2].write_text(args[2].read_text().replace('PF00001.1', 'PF99999.1'))
    elif bad == 'bounds':
        args[0].write_text(args[0].read_text().replace('\t1\t2\n', '\t1\t3\n'))
    elif bad == 'ambiguous_reference':
        with (args[3] / 'G1.fa').open('a') as f:
            f.write('>Homo_sapiens_isoform\nAS\n')
    else:
        (args[3] / 'G1.fasta').write_text('>Homo_sapiens\nAS\n')
    with pytest.raises(ValueError):
        migrate(*args)
    assert not args[-1].exists()
    assert not Path(str(args[-1]) + '.manifest.json').exists()


def test_producer_keeps_domain_when_clan_metadata_is_missing(tmp_path, monkeypatch):
    import compute_domain_variability as producer
    align = tmp_path / 'alignments'
    align.mkdir()
    (align / 'G1.fa').write_text('>Homo_sapiens\nAX-S\n')
    clans = tmp_path / 'clans.tsv'
    clans.write_text('')
    hmm = tmp_path / 'Pfam-A.hmm'
    hmm.write_text('test HMM')
    monkeypatch.setattr(producer, 'ensure_pfam_cache', lambda _: (str(hmm), str(clans)))
    def scan(argv, **kwargs):
        fields = ['x'] * 22
        for index, value in {0: 'NewDomain', 1: 'PF99999.1', 3: 'G1', 6: '1e-8', 11: '1e-10', 17: '1', 18: '2'}.items():
            fields[index] = value
        Path(argv[argv.index('--domtblout')+1]).write_text(' '.join(fields)+' new domain description\n')
    monkeypatch.setattr(producer.subprocess, 'run', scan)
    output = tmp_path / 'generated'
    monkeypatch.setattr(sys, 'argv', ['compute_domain_variability.py', '--alignment-dir', str(align),
                                    '--output-dir', str(output)])
    producer.main()
    with (output / 'domain_variability.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 1
    assert rows[0]['target_name'] == 'NewDomain' and rows[0]['description'] == 'new domain description'
    assert rows[0]['clan_acc'] == rows[0]['clan_name'] == 'NA'
    assert rows[0]['alignment_positions'] == '0,3'
