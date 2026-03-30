"""
ASR Reconstruction
==================

Ancestral state reconstruction using PAML.
"""

import logging
import subprocess
import shutil
import os
from pathlib import Path
from typing import Optional
from dataclasses import dataclass

from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo._io import write as phylo_write

logger = logging.getLogger(__name__)


# Updated ASRConfig to include asr_input_dir
# Mapping for substitution models
MODEL_SPECS = {
    "poisson": {"paml_model": 0, "aa_rate_file": None},
    "poison": {"paml_model": 0, "aa_rate_file": None},  # alias
    "proportional": {"paml_model": 1, "aa_rate_file": None},
    "dayhoff": {"paml_model": 3, "aa_rate_file": "dayhoff.dat"},
    "jtt": {"paml_model": 3, "aa_rate_file": "jtt.dat"},
    "wag": {"paml_model": 3, "aa_rate_file": "wag.dat"},
    "lg": {"paml_model": 3, "aa_rate_file": "lg.dat"},
}
DEFAULT_MODEL = "lg"


@dataclass
class ASRConfig:
    model: str = "lg"  # Substitution model (parameterized)
    fix_blength: int = 1  # Fix initial branch lengths using M0 (default: 1)
    threads: int = 1  # Number of threads
    compute_asr: bool = True  # Whether to compute ASR
    paml_binary: Optional[str] = None  # Path to PAML binary
    output_dir: Path = Path("output")  # Default output directory
    asr_input_dir: Optional[Path] = None  # Added asr_input_dir


class ASRReconstructor:
    """
    Main class for ancestral state reconstruction.

    Supports:
    - PAML (internal computation)
    - External ASR results (read pre-computed files)
    """

    def __init__(self, config: ASRConfig):
        """
        Initialize ASR reconstructor.

        Args:
            config: ASR configuration
        """
        self.config = config
        self.paml_binary = self._find_paml_binary()
        self._project_root = Path(__file__).resolve().parents[2]
        logger.info(f"PAML binary: {self.paml_binary}")
        logger.info(f"ASR model requested: {self.config.model}")
        self._validate_config()

    def _validate_config(self) -> None:
        """Validate ASR configuration."""
        if self.config.compute_asr:
            # Check if codeml (PAML binary) is available
            tool_path = shutil.which("codeml")
            if not tool_path:
                raise FileNotFoundError(
                    "PAML codeml binary not found in PATH. "
                    "Please install or provide --paml-binary with full path to codeml."
                )
            self.paml_binary = tool_path
            logger.info(f"Found PAML codeml at {tool_path}")
        else:
            if not self.config.asr_input_dir or not self.config.asr_input_dir.exists():
                raise ValueError(
                    "--asr-input-dir must be provided when --compute-asr is False"
                )
            logger.info(f"Using pre-computed ASR from {self.config.asr_input_dir}")

    def reconstruct_gene(
        self,
        gene: str,
        alignment: MultipleSeqAlignment,
        tree: Tree,
        output_dir: Path,
    ) -> Path:
        """
        Reconstruct ancestral states for a single gene.

        Args:
            gene: Gene identifier
            alignment: Multiple sequence alignment
            tree: Phylogenetic tree (pruned to alignment species)
            output_dir: Directory for output files

        Returns:
            Path to ASR results file (.rst for PAML)

        Raises:
            RuntimeError: If ASR computation fails
        """
        if not self.config.compute_asr:
            # Return path to pre-computed results
            return self._find_precomputed_asr(gene)

        # Compute ASR using PAML
        return self._run_paml_asr(gene, alignment, tree, output_dir)

    def _find_precomputed_asr(self, gene: str) -> Path:
        """
        Find pre-computed ASR file for a gene.

        Args:
            gene: Gene identifier

        Returns:
            Path to ASR file

        Raises:
            FileNotFoundError: If ASR file not found
        """
        if self.config.asr_input_dir is None:
            raise ValueError(
                "asr_input_dir is None, cannot search for pre-computed ASR files"
            )

        # Try different file patterns
        candidates = [
            self.config.asr_input_dir / f"asr_{gene}/rst",  # PAML
        ]

        for candidate in candidates:
            if candidate.exists():
                logger.debug(f"Found pre-computed ASR for {gene}: {candidate}")
                return candidate

        raise FileNotFoundError(
            f"No pre-computed ASR found for {gene} in {self.config.asr_input_dir}. "
            f"Tried: {[c.name for c in candidates]}"
        )

    def _run_paml_asr(
        self,
        gene: str,
        alignment: MultipleSeqAlignment,
        tree: Tree,
        output_dir: Path,
    ) -> Path:
        """
        Run PAML ancestral state reconstruction.

        Args:
            gene: Gene identifier
            alignment: Multiple sequence alignment
            tree: Phylogenetic tree
            output_dir: Output directory

        Returns:
            Path to rst file

        Raises:
            RuntimeError: If PAML fails
        """
        # Create gene-specific output directory
        gene_dir = output_dir / f"asr_{gene}"
        gene_dir.mkdir(parents=True, exist_ok=True)

        # Write alignment and tree to temporary files
        aln_file = gene_dir / "alignment_paml.phy"
        tree_file = gene_dir / "tree_paml.nwk"

        # Write in sequential Phylip format (required by PAML)
        # BioPython's writers sometimes have issues, so we write it manually
        self._write_sequential_phylip(alignment, aln_file)
        # For tree, use the BioPython tree writing function
        phylo_write(tree, tree_file, "newick")

        # PAML requires a blank line after the tree
        with open(tree_file, "a") as f:
            f.write("\n")

        # Create control file
        ctl_file = self._create_control_file(gene_dir)

        logger.debug(
            f"PAML inputs for {gene}: alignment={aln_file}, tree={tree_file}, ctl={ctl_file}"
        )

        try:
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(max(1, self.config.threads))
            logger.debug(
                "Running codeml for %s with OMP_NUM_THREADS=%s in %s",
                gene,
                env["OMP_NUM_THREADS"],
                gene_dir,
            )
            # Run codeml in the gene directory (PAML expects control file in working dir)
            result = subprocess.run(
                [str(self.paml_binary), str(ctl_file.name)],
                cwd=str(gene_dir),
                check=False,  # Don't raise on non-zero exit (PAML sometimes exits with 1 even on success)
                capture_output=True,
                text=True,
                timeout=3600,  # 1 hour timeout
                env=env,
            )

            # Check for .rst file (PAML may return exit code 1 even if successful)
            rst_file = gene_dir / "rst"
            if rst_file.exists():
                logger.info(f"  ✓ PAML completed (exit code {result.returncode})")
                rst_size = rst_file.stat().st_size
                logger.info(f"  ✓ RST file generated: {rst_file} ({rst_size} bytes)")
                return rst_file

            # If rst doesn't exist, check for actual errors
            if result.returncode != 0:
                logger.error(
                    f"PAML failed for {gene} with exit code {result.returncode}"
                )
                logger.error(f"STDOUT: {result.stdout}")
                logger.error(f"STDERR: {result.stderr}")
                raise RuntimeError(
                    f"PAML failed for {gene}: exit code {result.returncode}"
                )

            # No rst file and no error = unexpected
            raise RuntimeError(f"PAML did not produce rst file for {gene}")

        except subprocess.CalledProcessError as e:
            logger.error(f"PAML failed for {gene} with exit code {e.returncode}")
            logger.error(f"STDOUT: {e.stdout}")
            logger.error(f"STDERR: {e.stderr}")
            # Also try to read any PAML output files
            try:
                paml_out = gene_dir / "paml.out"
                if paml_out.exists():
                    with open(paml_out) as f:
                        logger.error(f"PAML output file: {f.read()}")
            except Exception:
                pass
            raise RuntimeError(f"PAML failed for {gene}: exit code {e.returncode}")
        except subprocess.TimeoutExpired:
            logger.error(f"PAML timeout for {gene} (>1 hour)")
            raise RuntimeError(f"PAML timeout for {gene}")

    def _create_control_file(self, output_dir: Path) -> Path:
        """
        Create PAML control file.

        Args:
            alignment_path: Path to alignment (phylip format)
            tree_path: Path to tree (newick format)
            output_dir: Output directory

        Returns:
            Path to created control file

        Mostly based on PAML documentation for codeml using their example for mouse lemurs.
        https://github.com/abacus-gene/paml/blob/master/examples/MouseLemurs/codeml.ctl
        """
        ctl_path = output_dir / "codeml.ctl"

        model_key = self._normalize_model_key(self.config.model)
        model_cfg = MODEL_SPECS.get(model_key, MODEL_SPECS[DEFAULT_MODEL])
        rate_file_path = None
        if model_cfg.get("aa_rate_file"):
            rate_file_path = self._resolve_rate_file(model_cfg["aa_rate_file"])
            logger.info(f"Using amino acid rate file: {rate_file_path}")
        logger.info(
            f"Using substitution model '{model_key}' (PAML model={model_cfg['paml_model']})"
        )

        # Create control file content
        # Note: files will be written as alignment_paml.phy and tree_paml.nwk
        content = f"""seqfile = alignment_paml.phy
treefile = tree_paml.nwk
outfile = rst

noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
verbose = 1  * 0: concise; 1: detailed, 2: too much
runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

# Evolutionary model parameters
model = {model_cfg['paml_model']}                        * models for codons:
                    * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                    * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
aaDist = 0  *  0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a        
aaRatefile = {rate_file_path if rate_file_path else ""} * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
            * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
            * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
            * 13:3normal>0

icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
Mgene = 0
            * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
            * AA: 0:rates, 1:separate

fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
    kappa = 2  * initial or fixed kappa
fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
    omega = .4 * initial or fixed omega, for codons or codon-based AAs
fix_alpha = 0  * 0: estimate gamma shape parameter; 1: fix it at alpha
    alpha = 1.0 * initial or fixed alpha, 0:infinity (constant rate)
    Malpha = 0  * different alphas for genes
    ncatG = 8  * # of categories in dG of NSsites models # indicate 8 category

clock = 0 * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
getSE = 0
RateAncestor = 1

Small_Diff = .5e-6
cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?
fix_blength = {self.config.fix_blength} * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
method = 0 * Optimization method 0: simultaneous; 1: one branch a time
"""

        # Write control file
        with open(ctl_path, "w") as f:
            f.write(content)

        if self.config.fix_blength == 2:
            logger.info(
                "Preserving supplied tree topology and branch lengths during ASR (fix_blength=2)"
            )
        logger.debug(f"Created PAML control file: {ctl_path}")
        return ctl_path

    def _write_sequential_phylip(
        self, alignment: MultipleSeqAlignment, output_path: Path
    ) -> None:
        """
        Write alignment in sequential Phylip format (required by PAML).

        Args:
            alignment: Multiple sequence alignment
            output_path: Path to write Phylip file
        """
        nseq = len(alignment)
        seq_length = alignment.get_alignment_length()

        with open(output_path, "w") as f:
            # Header line
            f.write(f"{nseq:>6} {seq_length:>6}\n")

            # Write sequences one per line in sequential format
            for record in alignment:
                # Name: up to 10 characters, padded
                name = record.id[:10].ljust(10)
                # Sequence on same line - keep gaps as '-', PAML handles them
                # Replace non-standard amino acids with X
                seq_str = (
                    str(record.seq)
                    .replace("O", "X")
                    .replace("U", "X")
                    .replace("B", "X")
                    .replace("Z", "X")
                    .replace("J", "X")
                )
                f.write(f"{name} {seq_str}\n")

    def _normalize_model_key(self, model_name: Optional[str]) -> str:
        if not model_name:
            return DEFAULT_MODEL
        token = str(model_name).lower().split("+")[0].strip()
        if token not in MODEL_SPECS:
            logger.warning(
                f"Unrecognized substitution model '{model_name}', falling back to '{DEFAULT_MODEL}'"
            )
            return DEFAULT_MODEL
        return token

    def _resolve_rate_file(self, filename: str) -> Path:
        search_paths = [
            self._project_root / "dat" / filename,
        ]
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            conda_path = Path(conda_prefix)
            search_paths.extend(
                [
                    conda_path / "dat" / filename,
                    conda_path / "share" / "paml" / "dat" / filename,
                ]
            )
        for path in search_paths:
            if path.exists():
                return path
        raise FileNotFoundError(
            f"Rate file '{filename}' not found. Searched: {', '.join(str(p) for p in search_paths)}"
        )

    # Implement _find_paml_binary method
    def _find_paml_binary(self) -> Path:
        """Find the PAML codeml binary."""
        if self.config.paml_binary:
            return Path(self.config.paml_binary)
        codeml_path = shutil.which("codeml")
        if codeml_path:
            return Path(codeml_path)
        raise FileNotFoundError(
            "PAML codeml binary not found. Please specify the path."
        )
