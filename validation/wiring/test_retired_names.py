"""No tracked file outside a short allowlist names a process or parameter of the former permulation-null layout."""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2]))
RETIRED = re.compile(r"PERM_REPLAY_BATCHED|\bPERM_REPLAY\b|CAAS_PERMS_(?:DISAMBIGUATE|AGGREGATE|REBUILD|MERGE_DETAIL)|process_perm_replay"
                     r"|ct_perm_replay_batch_size|ct_disambig_perms_batch_size|CT_PERM_REPLAY_BATCH_SIZE|CT_DISAMBIG_PERMS_BATCH_SIZE"
                     r"|caas_b0_diagnostic")
# history and frozen artifacts: the archive, measurement records, generated run scripts and templates under validation/,
# the style archetypes, and the migration that names the retired project fields
ALLOWED = ("archive/", "validation/", "style/", "docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md", "gui/models/serialization.py")


def test_no_tracked_file_outside_the_allowlist_names_the_retired_layout():
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split("\0")
    hits = []
    for name in filter(None, files):
        if name.startswith(ALLOWED):
            continue
        try:
            text = (ROOT / name).read_text()
        except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
            continue
        hits += [f"{name}:{n}: {line.strip()[:100]}" for n, line in enumerate(text.splitlines(), 1) if RETIRED.search(line)]
    assert not hits, "\n".join(hits)


RETIRED_OBSERVED = re.compile(r"DISCOVERY_BATCHED|CONCAT_DISCOVERY|CONCAT_BACKGROUND|CT_DISAMBIGUATION_(?:RUN|SPLIT_GENES|MERGE|PLOTS)"
                              r"|\bCAAS_PERMULATION\b|run_ct_discovery_batch|process_discovery_batched|ct_discovery\.nf|\basr_ready\b"
                              r"|merge_disambiguation_batches|regenerate_disambiguation_plots"
                              r"|disambiguation_main|disambiguation_db|disambiguation_writers|disambiguation_json|generate_bulk_plots"
                              r"|split_meta_caas_by_genes|\bprocess_all_genes\b|process_single_gene|export_from_db|list_gene_caas_positions"
                              r"|ct_discovery_batch_size|ct_disambig_batch_size|\bct_disambig_asr_mode\b|CT_DISCOVERY_BATCH_SIZE|DISAMBIG_BATCH_SIZE|\bASR_MODE\b")
# the frozen records, the archive and the style archetypes keep their history; the migration names the retired project fields
ALLOWED_OBSERVED = ("archive/", "validation/", "style/", "docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md", "gui/models/serialization.py")


def test_no_tracked_file_names_a_process_of_the_former_observed_chain():
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split("\0")
    hits = []
    for name in filter(None, files):
        if name.startswith(ALLOWED_OBSERVED):
            continue
        try:
            text = (ROOT / name).read_text()
        except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
            continue
        hits += [f"{name}:{n}: {line.strip()[:100]}" for n, line in enumerate(text.splitlines(), 1) if RETIRED_OBSERVED.search(line)]
    assert not hits, "\n".join(hits)


ALLOWED_KERNEL = ("archive/", "validation/", "style/", "docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md")
RETIRED_KERNEL = re.compile(r"modules\.disco\b|\bdisco\.py\b|\bct discovery\b|perm_replay\.output|collapse_fop_hits_by_base|parse_discovery_positions"
                            r"|recovery_boot|simtrait_revive_from_dir|--progress[-_]log")


def test_no_tracked_file_names_the_scalar_discovery_or_the_perm_replay_counts_path():
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split("\0")
    hits = []
    for name in filter(None, files):
        if name.startswith(ALLOWED_KERNEL):
            continue
        try:
            text = (ROOT / name).read_text()
        except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
            continue
        hits += [f"{name}:{n}: {line.strip()[:100]}" for n, line in enumerate(text.splitlines(), 1) if RETIRED_KERNEL.search(line)]
    assert not hits, "\n".join(hits)


# the scorer's legacy inputs (metadata read from a file, trait file, diagnostics dump) and the modules nothing reaches;
# the replay-performance record and the frozen artifacts keep their history
ALLOWED_SCORER = ("archive/", "validation/", "style/", "docs/CT_DISAMBIGUATION_REPLAY_PERFORMANCE.md")
RETIRED_SCORER = re.compile(r"\bcaas_metadata_path\b|\bdiagnostics_dir\b|\bcaas_hypothesis_domain_asr\b|\btip_dump_file\b"
                            r"|list_gene_caas_entries|read_caas_metadata_table|\bparse_trait_pairs\b"
                            r"|analyze_(?:gene|caas_position)_biochemistry|convert_biochem_result_to_dict"
                            r"|node_identification|\basr_only\b|identify_convergence_nodes|validate_node_mapping|\bget_descendants\b"
                            r"|export_posteriors_to_jsonl|\bparse_paml_rst\b|build_consolidated_multiset|validate_taxids_in_tree"
                            r"|\blabel_nodes\b|\bhas_polytomies\b|get_root_to_tip_paths|build_tip_node_lookup|build_node_path_from_mapping"
                            r"|ct_disambig_run_diagnostics|\brun_diagnostics\b|\bvalidate_asr_inputs\b")


def test_no_tracked_file_names_a_legacy_input_of_the_scorer_or_an_unreachable_module():
    files = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split("\0")
    hits = []
    for name in filter(None, files):
        if name.startswith(ALLOWED_SCORER):
            continue
        try:
            text = (ROOT / name).read_text()
        except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
            continue
        hits += [f"{name}:{n}: {line.strip()[:100]}" for n, line in enumerate(text.splitlines(), 1) if RETIRED_SCORER.search(line)]
    assert not hits, "\n".join(hits)
