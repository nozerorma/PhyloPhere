// Copy to the root of a stable snapshot and rewrite ../ imports to ./ before running.
nextflow.enable.dsl=2
include { DOMINO_MODULES } from '../subworkflows/ENRICHMENT/domino.nf'
workflow {
 DOMINO_MODULES(Channel.value(file(params.wp8_background)),Channel.value(file(params.wp8_genes)),
                700,0.3,0.05,params.wp8_string_db,'caas')
 DOMINO_MODULES.out.modules_dir.subscribe { d -> assert d.listFiles().size()==0; println 'EMPTY_DOMINO_PASS' }
}
