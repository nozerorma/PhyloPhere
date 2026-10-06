#!/usr/bin/env nextflow
// rer_trait.nf — Trait vector for RERconverge and detection of its type (binary or continuous).
// PhyloPhere | subworkflows/RERCONVERGE/

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  RER_TRAIT: runs build_rer_trait.R, which reads the trait file, builds the named
 *  trait vector (and the optional count vectors for the Haldane-Anscombe logit) and
 *  classifies the trait as binary (two values, recoded to 0/1 if needed) or
 *  continuous. The type file routes the analysis to RER_BIN or RER_CONT.
 *
 *  Consumes:  trait file (species column params.sp_colname, trait column params.traitname)
 *  Produces:  <trait>.polished.output (RData), <trait>.trait_type.output ('binary' or
 *             'continuous')
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ── Trait vector ───────────────────────────────────────────────────────────────

process RER_TRAIT {
    tag "$my_traitfile"

    label 'process_low'


    publishDir path: "${params.outdir}/rerconverge/rer_traits/", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    input:
    path my_traitfile

    output:
    path "${params.traitname}.polished.output", emit: polished
    path "${params.traitname}.trait_type.output", emit: trait_type

    script:
    def outputName   = "${params.traitname}.polished.output"
    def typeOutName  = "${params.traitname}.trait_type.output"

    if (params.use_singularity || params.use_apptainer) {
        """
        echo "Using Singularity/Apptainer"
        /usr/local/bin/_entrypoint.sh Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/build_rer_trait.R' \\
        ${ my_traitfile } \\
        ${ params.sp_colname } \\
        ${ params.traitname } \\
        ${ outputName } \\
        ${ typeOutName } \\
        "${ params.n_trait ?: '' }" \\
        "${ params.c_trait ?: '' }" \\
        "${ params.rer_transform ?: 'auto' }"
        """
    } else {
        """
        echo "Running locally"
        Rscript \\
        '$baseDir/subworkflows/RERCONVERGE/local/build_rer_trait.R' \\
        ${ my_traitfile } \\
        ${ params.sp_colname } \\
        ${ params.traitname } \\
        ${ outputName } \\
        ${ typeOutName } \\
        "${ params.n_trait ?: '' }" \\
        "${ params.c_trait ?: '' }" \\
        "${ params.rer_transform ?: 'auto' }"
        """
    }

}
