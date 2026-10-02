// Runs the chain main.nf builds from the permulation core: CAAS_CORE (batches), CAAS_CORE_OBSERVED (the observed
// files from the b_0 slices) and CAAS_CORE_MERGE; writes the paths of the observed files and of the null's
// caas_perms.rds to <outdir>/chain_paths.txt. Used as the main.nf of a project that links the tree under test.
include { CAAS_CORE; CAAS_CORE_OBSERVED; CAAS_CORE_MERGE } from './subworkflows/CT/caas_permulation'

workflow {
    def alignments = Channel.fromList(params.mini_alignments.split(',').collect { f -> tuple(file(f).baseName, file(f)) })
    def cfg = Channel.value(file(params.mini_cfg))
    def core = CAAS_CORE(alignments, Channel.empty(), cfg, Channel.value(file(params.mini_resample)), Channel.value(file(params.mini_tree)),
                         Channel.value(file(params.mini_fop_pairs)), Channel.value(file(params.mini_lengths)))
    def obs = CAAS_CORE_OBSERVED(core.b0_observed.collect(), cfg.collect().map { items -> items[0] })
    def merged = CAAS_CORE_MERGE(core.pos_detail.collect(), Channel.value(file('NO_FILE')), Channel.value(file(params.mini_lengths)))
    obs.discovery.mix(obs.background, obs.background_genes, obs.master_csv, obs.global_meta_caas, merged.perms)
        .map { f -> f.toString() }.collectFile(name: 'chain_paths.txt', storeDir: params.outdir, newLine: true)
}
