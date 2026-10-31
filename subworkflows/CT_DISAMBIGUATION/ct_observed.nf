#!/usr/bin/env nextflow

// All scientific file inputs are staged. task.ext makes the copied implementation
// fingerprint part of the task cache key, including shared readers and rate tables.
def observedSourceHash(local_dir, reader, extra_dir=null) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    def sources = []
    file(local_dir).toFile().eachFileRecurse(groovy.io.FileType.FILES) { f ->
        if (f.name.endsWith('.py') || f.name.endsWith('.dat')) sources << f
    }
    if (extra_dir && file(extra_dir).isDirectory()) {
        file(extra_dir).toFile().eachFileRecurse(groovy.io.FileType.FILES) { f ->
            if (f.name.endsWith('.py') || f.name.endsWith('.sh') || f.name == 'ct') sources << f
        }
    }
    sources.sort { it.path }.each { f -> md.update(f.name.getBytes('UTF-8')); md.update(f.bytes) }
    md.update(file(reader).bytes)
    md.update(file("${local_dir}/../ct_observed.nf").bytes)
    return md.digest().encodeHex().toString()
}

def observed_source_hash = observedSourceHash("${baseDir}/subworkflows/CT_DISAMBIGUATION/local", "${baseDir}/subworkflows/CT/local/modules/design_io.py")
def observedQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

// Nextflow hashes directory metadata, even in deep mode. Fingerprint nested
// scientific inputs explicitly while retaining their original staged layout.
def observedInputHash(path) {
    def root = file(path).toFile()
    def md = java.security.MessageDigest.getInstance('SHA-256')
    def files = []
    if (root.isDirectory()) {
        root.eachFileRecurse(groovy.io.FileType.FILES) { files << it }
    } else if (root.isFile()) {
        files << root
    }
    files.sort { it.path }.each { f ->
        md.update((root.isDirectory() ? root.toPath().relativize(f.toPath()).toString() : f.name).getBytes('UTF-8'))
        md.update((byte)0)
        f.withInputStream { stream ->
            byte[] buffer = new byte[1048576]
            int n
            while ((n = stream.read(buffer)) != -1) md.update(buffer, 0, n)
        }
        md.update((byte)0)
    }
    return md.digest().encodeHex().toString()
}

process CAAS_OBSERVED {
    tag 'caas_observed'
    label 'process_resample'
    // Hash regular file inputs by content; directory fingerprints arrive as values.
    cache 'deep'
    ext observed_source_hash: observed_source_hash
    publishDir path: "${params.outdir}", mode: 'copy', overwrite: true, pattern: 'ct_disambiguation'
    publishDir path: "${params.outdir}/meta_caas", mode: 'copy', overwrite: true, pattern: 'meta_caas'

    input:
    path discovery, stageAs: 'inputs/discovery/*'
    tuple val(design_hash), path(trait_file, stageAs: 'inputs/design/*')
    path tree_file, stageAs: 'inputs/tree/*'
    path hyp_pairs, stageAs: 'inputs/weights/*'
    path taxid_map, stageAs: 'inputs/taxid/*'
    tuple val(alignment_hash), path(alignments, stageAs: 'inputs/alignments/*')
    path ensembl, stageAs: 'inputs/ensembl/*'

    output:
    path 'ct_disambiguation', emit: results_dir
    path 'ct_disambiguation/caas_convergence_master.csv', emit: master_csv
    path 'meta_caas', emit: meta_caas_dir
    path 'meta_caas/meta_caas.tsv', emit: meta_caas

    script:
    def local_dir = "${baseDir}/subworkflows/CT_DISAMBIGUATION/local"
    def reader = "${baseDir}/subworkflows/CT/local/modules/design_io.py"
    def cache_dir = params.ct_disambig_asr_cache_dir ?: 'asr_cache'
    def args = ['--alignment-dir', alignments.toString(), '--tree', tree_file.toString(),
                '--discovery', discovery.toString(), '--design', trait_file.toString(), '--output-dir', 'shards',
                '--asr-model', (params.ct_disambig_asr_model ?: 'lg').toString(),
                '--posterior-threshold', (params.ct_disambig_posterior_threshold != null ? params.ct_disambig_posterior_threshold : 0.1).toString(),
                '--workers', task.cpus.toString(), '--max-tasks-per-child', (params.ct_disambig_max_tasks_per_child ?: 50).toString(),
                '--asr-cache-dir', cache_dir.toString(), '--ali-format', (params.ali_format ?: 'auto').toString()]
    if (!new File(hyp_pairs.toString()).name.startsWith('NO_')) args += ['--fop-pairs', hyp_pairs.toString()]
    if (!new File(taxid_map.toString()).name.startsWith('NO_')) args += ['--taxid-mapping', taxid_map.toString()]
    if (!new File(ensembl.toString()).name.startsWith('NO_')) args += ['--ensembl-genes-file', ensembl.toString()]
    def payload = groovy.json.JsonOutput.toJson([
        scoring: args,
        contract: ['--b0-dirs', 'shards', '--design', trait_file.toString(), '--discovery-file', discovery.toString(), '--output-dir', '.']])
    def run = (params.use_singularity || params.use_apptainer) ? '/usr/local/bin/_entrypoint.sh python3' : 'python3'
    """
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    echo '${task.ext.observed_source_hash}' > observed_source.sha256
    cp -R ${observedQuote(local_dir)}/. .
    cp ${observedQuote(reader)} src/core/_ct_design_io.py
    cat > observed_args.json <<'PHYLOPHERE_OBSERVED_ARGS'
${payload}
PHYLOPHERE_OBSERVED_ARGS
    ${run} run_observed.py observed_args.json
    """
}
