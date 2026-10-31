#!/usr/bin/env nextflow

def rerQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

def rerHash(paths) {
    def md = java.security.MessageDigest.getInstance('SHA-256')
    paths.each { original ->
        def p = original
        if (original.name != 'NO_FILE') p = original.toRealPath()
        if (p.name == 'NO_FILE' && !p.exists()) { md.update('NO_FILE'.bytes); return }
        assert p.exists() : "RER input not found: ${p}"
        def files = [p]
        // Directory iteration below uses a stable relative-name/content sequence.
        if (p.isDirectory()) {
            files=[]
            p.toFile().eachFileRecurse(groovy.io.FileType.FILES) { files << it.toPath() }
        }
        files.sort { it.toString() }.each { f ->
            md.update((p.isDirectory() ? p.relativize(f).toString() : p.name).getBytes('UTF-8'))
            md.update((byte)0)
            f.toFile().withInputStream { stream ->
                byte[] b = new byte[1048576]; int n
                while ((n=stream.read(b)) != -1) md.update(b,0,n)
            }
            md.update((byte)0)
        }
    }
    md.digest().encodeHex().toString()
}

def rerSources(names) {
    rerHash([file("${baseDir}/subworkflows/RERCONVERGE/rer_contract.nf"),
             file("${baseDir}/subworkflows/RERCONVERGE/local/rer_contract.R")] +
             names.collect { file("${baseDir}/subworkflows/RERCONVERGE/${it}") })
}

def rerSingle(ch, label) {
    ch.toList().map { xs ->
        assert xs.size()==1 : "RER: expected one ${label}, got ${xs.size()}"
        assert xs[0].isFile() || xs[0].isDirectory() : "RER: ${label} not found: ${xs[0]}"
        xs[0]
    }
}

def rerRequired(value, label) {
    assert value : "RER requires ${label}"
    def p=file(value)
    assert p.exists() : "RER ${label} not found: ${p}"
    Channel.value(p)
}

// Compute source fingerprints before dispatch, so they are explicit cache inputs.
def rerStage(stage) {
    def stages = [trait:['local/build_rer_trait.R','rer_trait.nf'],
        trees:['local/rer_master_tree.R','rer_trees.nf'],
        matrix:['local/rer_matrix.R','rer_matrix.nf'],
        mode:['local/rer_trait_mode.R','rer_trait_mode.nf'],
        continuous:['local/continuous_rer.R','local/rer_permulations.R','local/rer_batches.R','rer_cont.nf'],
        binary:['local/binary_rer.R','local/rer_permulations.R','local/rer_batches.R','rer_bin.nf'],
        tables:['local/rer_tables.R','rer_tables.nf'],
        lists:['local/rer_gene_lists.R','rer_gene_lists.nf'],
        report:['local/5.RERconverge_report.Rmd','local/render_rer.R','rer_report.nf']]
    assert stages.containsKey(stage)
    def h = rerSources(stages[stage])
    if (stage=='report') h += rerHash([file("${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R")])
    h
}
