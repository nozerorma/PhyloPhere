"""Regression checks for legacy inputs and retained scientific safeguards.
Run: conda run --no-capture-output -n phylophere python -B -m unittest discover -s checks -v
"""
import csv
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SELECTION = ROOT / 'subworkflows/SELECTION/local/src'
sys.path.insert(0, str(SELECTION))
from extract_extreme_species import load_species_sets
from annotate_tree_fg import annotate


class CompatibilityTests(unittest.TestCase):
    def test_empty_pools_are_biological_exclusions(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / 'pool').write_text('')
            self.assertEqual(load_species_sets(d / 'pool'), (set(), set()))
            (d / 'pool').write_text('A\t1\nA\t0\n')
            with self.assertRaises(ValueError):
                load_species_sets(d / 'pool')
            (d / 'pool').write_text('')
            (d / 'tree').write_text('((A:.1,B:.1):.1,(C:.1,D:.1):.1);')
            (d / 'fasta').write_text('>A\nACDE\n>B\nACDE\n>C\nACDE\n>D\nACDE\n')
            result = annotate(d / 'tree', d / 'fasta', d / 'pool', d / 'out')
            self.assertEqual((result['status'], result['reason']), ('excluded', 'no_foreground'))
            self.assertFalse((d / 'out').exists())
            (d / 'pool').write_text('A\n')
            (d / 'tree').write_text('(A:.1,B:.1,C:.1,D:.1);')
            with self.assertRaises(ValueError):
                annotate(d / 'tree', d / 'fasta', d / 'pool', d / 'out')

    def test_postproc_cli_accepts_missing_completion_without_fabricating_science(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            with (d / 'master.csv').open('w') as f:
                w = csv.writer(f)
                w.writerow(['gene', 'msa_pos', 'side', 'asr_path_score'])
                w.writerow(['NA', 0, 'top', '.8'])
            run = subprocess.run([sys.executable, '-B', str(ROOT / 'subworkflows/CT_POSTPROC/local/src/prepare_postproc_input.py'),
                                  '--input', str(d / 'master.csv'), '--completion', ''], cwd=d, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            saved = json.loads((d / 'quality_universe.json').read_text())
            self.assertIsNone(saved['scientific'])
            self.assertIsNone(saved['completion_sha256'])
            self.assertIn('NA', saved['genes'])

    def test_gui_report_metadata_are_optional(self):
        sys.path.insert(0, str(ROOT))
        from gui.generation.report_registry import Listing, _postproc_find_slots, _signif_significance_find_slots
        with tempfile.TemporaryDirectory() as d:
            listing = Listing(d)
            postproc = {s.key: s for s in _postproc_find_slots(listing)}
            significance = {s.key: s for s in _signif_significance_find_slots(listing)}
            self.assertFalse(postproc['configuration_file'].required)
            self.assertFalse(postproc['quality_universe_file'].required)
            self.assertTrue(postproc['discovery_file'].required)
            self.assertFalse(significance['gene_metadata_input'].required)
            self.assertTrue(significance['position_scores_input'].required)

    @unittest.skipUnless(shutil.which('hyphy'), 'HyPhy is required for the native FADE check')
    def test_native_fade_reuses_baseline_without_bulk_hashes(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / 'input.fa').write_text('>A\nACDEFGHIKLMNPQRSTVWY\n>X.A\nACDEYGHIKLMNPQRSTVWY\n>B\nVCDEFGHIKLMNPQRSTAWY\n>C\nACDEFGHIKLMNPQKSTVWY\n')
            (d / 'tree.nwk').write_text('((A:.1,X.A:.1):.1,(B:.1,C:.1):.1);')
            (d / 'top').write_text('A\n')
            (d / 'bottom').write_text('B\n')
            (d / 'manifest.tsv').write_text('g\ttop\tinput.fa\ttree.nwk\ttop\tNO_FILE\n' +
                                           'g\tbottom\tinput.fa\ttree.nwk\tbottom\tNO_FILE\n')
            settings = dict(background_scope='all', internal_nodes='all_descendants', model='LG',
                            method='Variational-Bayes', grid=5, concentration=.5, bf_threshold=125)
            (d / 'settings.json').write_text(json.dumps(settings))
            fade = ROOT / 'subworkflows/FADE/local/src'
            env = os.environ.copy()
            env.pop('PHYLOPHERE_VERIFY_CHECKSUMS', None)
            run = subprocess.run([sys.executable, '-B', str(fade / 'run_fade.py'), '--manifest', 'manifest.tsv',
                                  '--settings', 'settings.json', '--batch-id', 'batch'], cwd=d, env=env,
                                 text=True, capture_output=True, timeout=120)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            audit = json.loads((d / 'g.bottom.fade.audit.json').read_text())
            self.assertTrue(audit['baseline_reused'])
            self.assertIn('alignment', audit['input_identity'])
            self.assertNotIn('alignment', audit['input_sha256'])
            (d / 'expected.json').write_text(json.dumps([['g', 'top'], ['g', 'bottom']]))
            run = subprocess.run([sys.executable, '-B', str(fade / 'fade_merge_tables.py'), '--direction', 'top',
                                  '--genes', 'batch.fade_genes.tsv', '--sites', 'batch.fade_sites.tsv',
                                  '--audits', 'g.top.fade.audit.json', 'g.bottom.fade.audit.json',
                                  '--expected', 'expected.json', '--settings', 'settings.json'], cwd=d, env=env,
                                 text=True, capture_output=True, timeout=30)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            metadata = json.loads((d / 'fade_run_metadata.json').read_text())
            self.assertEqual(metadata['output_sha256'], {})
            self.assertEqual(metadata['tested'], 1)

    def test_native_rer_contracts_and_complete_reports(self):
        with tempfile.TemporaryDirectory() as d:
            env = os.environ.copy()
            env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', VROOM_THREADS='1')
            run = subprocess.run(['Rscript', str(ROOT / 'checks/workflow_compatibility.R'), str(ROOT), d],
                                 cwd=d, env=env, text=True, capture_output=True, timeout=240)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            print(run.stdout)


if __name__ == '__main__':
    unittest.main()
