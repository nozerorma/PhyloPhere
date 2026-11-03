"""Regression coverage for the CT discovery-row overlap contract."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path

LOCAL = Path(__file__).resolve().parent
sys.path.insert(0, str(LOCAL))
sys.path.insert(0, str(LOCAL.parents[1] / 'CT' / 'local'))

from modules.caas_id import residue_overlap
from src.core.contract import write_meta
from src.core.discovery import DiscoveryRows
from src.utils.amino import shared_residue_count


class DiscoveryOverlapTest(unittest.TestCase):
    def test_matches_ct_with_repeated_and_unresolved_residues(self):
        for fg, bg in [('AAAC', 'AADD'), ('ACD', 'DCA'), ('AAAC', 'DDDD'),
                       ('aa-XC', 'A?AD'), ('AAAAC', 'ADDD')]:
            with self.subTest(fg=fg, bg=bg):
                self.assertEqual(shared_residue_count(fg, bg), residue_overlap(fg, bg))

    def test_discovery_and_meta_preserve_occurrence_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / 'gene.b0.discovery.tsv'
            path.write_text('gene\tposition\tcaas\tpattern\tis_conserved_meta\tshared_residue_count\n'
                            'gene\t1\tAAAC/AADD\t4\tTRUE\t2\n')
            rows = DiscoveryRows(root, root)
            try:
                self.assertEqual(rows['gene'][0]['shared_residue_count'], '2')
                write_meta(path, root / 'meta')
                with (root / 'meta' / 'meta_caas.tsv').open() as fh:
                    self.assertEqual(next(csv.DictReader(fh, delimiter='\t'))['shared_residue_count'], '2')
            finally:
                rows.close()
            path.write_text(path.read_text().replace('TRUE\t2', 'TRUE\t1'))
            with self.assertRaisesRegex(ValueError, r'gene\.b0\.discovery\.tsv:2:.*expected 2'):
                DiscoveryRows(root, root)


if __name__ == '__main__':
    unittest.main()
