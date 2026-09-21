seqfile = alignment_paml.phy
treefile = tree_paml.nwk
outfile = rst

noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
verbose = 1  * 0: concise; 1: detailed, 2: too much
runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

# Evolutionary model parameters
model = 3                        * models for codons:
                    * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                    * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
aaDist = 0  *  0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a        
aaRatefile = /home/miguel/IBE-UPF/PhD/PhyloPhere/validation/tier1/output/pepc/work/c4_complete/98/8b691df47e09cc76f0cf1d4e558f92/dat/lg.dat * only used for aa seqs with model=empirical(_F)
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
fix_blength = 2 * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional
method = 0 * Optimization method 0: simultaneous; 1: one branch a time
