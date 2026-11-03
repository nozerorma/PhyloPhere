import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('network', Path(__file__).parent / 'src' / 'build_domino_network.py')
network = importlib.util.module_from_spec(spec)
spec.loader.exec_module(network)


class IdentifierMappingTest(unittest.TestCase):
    def test_network_keeps_original_ids_for_multiple_species(self):
        background = {'GENE.Homo_sapiens', 'AC010547.4.Papio_anubis', 'OTHER.Mus_musculus'}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'links.gz'
            with gzip.open(path, 'wt') as out:
                out.write('protein1 protein2 combined_score\np1 p2 800\np2 p3 850\n')
            edges = network.build_edges(path, {'p1': 'GENE', 'p2': 'AC010547.4', 'p3': 'OTHER'}, background, 700)
        self.assertEqual(edges, {('AC010547.4.Papio_anubis', 'GENE.Homo_sapiens'): 800,
                                 ('AC010547.4.Papio_anubis', 'OTHER.Mus_musculus'): 850})

    def test_ambiguous_alias_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Ambiguous STRING symbol'):
            network.resolve_background_symbols({'GENE.Homo_sapiens', 'GENE.Papio_anubis'}, {'GENE'})

    def test_zero_coverage_is_an_error(self):
        with self.assertRaisesRegex(ValueError, 'No background identifiers match'):
            network.build_edges('unused', {'p': 'UNRELATED'}, {'GENE.Homo_sapiens'}, 700)


if __name__ == '__main__':
    unittest.main()
