#!/usr/bin/env nextflow
// Human-reference characterization of CAAS candidate states. Public channels retained.
include { VEP_PREPARE } from '../subworkflows/VEP/vep_prepare'
include { PRIMATEAI_MAP } from '../subworkflows/VEP/primateai'
include { COSMIC_MAP } from '../subworkflows/VEP/cosmic'
include { ENSEMBL_VEP_ANNOTATE } from '../subworkflows/VEP/ensembl_vep'
include { vepHash; vepSources } from '../subworkflows/VEP/vep_contract'

workflow VEP {
    take:
    position_scores_input

    main:
    def primateai_out = Channel.empty()
    def cosmic_out = Channel.empty()
    def ensembl_out = Channel.empty()
    def pai = params.vep_primateai_db ?: ''
    def cosmic = params.cosmic_db ?: ''
    def ensembl = params.vep_ensembl ?: false
    if (pai || cosmic || ensembl) {
        assert params.caas_map_dir && file(params.caas_map_dir).isDirectory() : 'VEP requires an existing --caas_map_dir'
        assert (params.vep_species ?: 'homo_sapiens') == 'homo_sapiens' : 'VEP reference contract requires homo_sapiens'
        assert (params.vep_assembly ?: 'GRCh38') == 'GRCh38' : 'VEP reference contract requires GRCh38'
        def maps = file(params.caas_map_dir)
        def proteins = params.gene_ensembl_file ? file(params.gene_ensembl_file, checkIfExists: true) : []
        assert !ensembl || proteins : 'Ensembl VEP requires --gene_ensembl_file'
        def sources = vepSources(['vep_prepare.nf','local/src/vep_common.py','local/src/prepare_targets.py'])
        // Check all configured sources before dispatching tasks; no bulk hashes at startup.
        [pai,cosmic].findAll { it }.each { name ->
            def db=file(name,checkIfExists:true)
            assert db.isFile() && db.size()>0 : "VEP database must be a nonempty file: ${name}"
        }
        def alignments = params.alignment ? file(params.alignment,checkIfExists:true) : []
        def scores = (position_scores_input ?: Channel.empty()).collect().ifEmpty { [] }.map { files ->
            assert files.size() == 1 : "VEP requires exactly one position_scores.tsv, received ${files.size()}"
            files[0]
        }
        // Directory inputs need an inventory of names/size/mtime for normal resume.
        // vepHash uses metadata for data files; it never reads database/MAP bytes here.
        def refs_key=vepHash([maps]+(proteins ? [proteins] : [])+(alignments ? [alignments] : []))
        VEP_PREPARE(sources+refs_key, scores, Channel.value(maps), Channel.value(proteins), Channel.value(alignments))
        if (pai) {
            def db = file(pai, checkIfExists: true)
            assert db.isFile() && db.size()>0 : 'PrimateAI database must be a nonempty file'
            def hash = vepSources(['primateai.nf','local/src/vep_common.py','local/src/map_database.py','local/src/map_to_primateai.py'])
            PRIMATEAI_MAP(hash,VEP_PREPARE.out.targets,Channel.value(db))
            primateai_out = PRIMATEAI_MAP.out.primateai_tsv
        }
        if (cosmic) {
            def db = file(cosmic, checkIfExists: true)
            assert db.isFile() && db.size()>0 : 'COSMIC database must be a nonempty file'
            def hash = vepSources(['cosmic.nf','local/src/vep_common.py','local/src/map_database.py','local/src/map_to_cosmic.py'])
            COSMIC_MAP(hash,VEP_PREPARE.out.targets,Channel.value(db))
            cosmic_out = COSMIC_MAP.out.cosmic_tsv
        }
        if (ensembl) {
            def mode = params.vep_mode ?: 'auto'
            assert mode in ['auto','offline','database'] : 'VEP mode must be auto, offline or database'
            def cache_path = file(params.vep_cache_dir ?: "${System.properties['user.home']}/.cache/phylophere/vep/homo_sapiens_GRCh38")
            if (mode == 'auto') mode = cache_path.resolve('homo_sapiens').isDirectory() ? 'offline' : 'database'
            log.info "VEP: resolved source ${mode}; resources are never downloaded during annotation"
            def cache_input = mode == 'offline' ? cache_path : []
            def settings = [mode:mode,host:params.vep_db_host ?: 'ensembldb.ensembl.org',port:params.vep_db_port ?: 3306]
            def hash = vepSources(['ensembl_vep.nf','local/src/run_ensembl.py','local/src/build_vep_variants.py','local/src/join_vep_output.py','local/src/resolve_ensembl.py','local/src/vep_common.py'])
            // Database source is external and mutable; do not reuse its annotation task.
            if (mode == 'database') hash += java.util.UUID.randomUUID().toString()
            ENSEMBL_VEP_ANNOTATE(hash,VEP_PREPARE.out.targets,Channel.value(cache_input),Channel.value(settings))
            ensembl_out = ENSEMBL_VEP_ANNOTATE.out.ensembl_vep_tsv
        }
    } else {
        log.info 'VEP: no annotation source enabled'
    }

    emit:
    primateai_tsv = primateai_out
    cosmic_tsv = cosmic_out
    ensembl_vep_tsv = ensembl_out
}
