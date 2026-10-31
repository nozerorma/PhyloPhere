#!/usr/bin/env nextflow
include { selectionHash; selectionQuote; selectionPositive } from '../SELECTION/selection_contract'

def fadeSources(stage) {
    def common = [file("${baseDir}/subworkflows/FADE/fade_contract.nf"),
                  file("${baseDir}/subworkflows/SELECTION/selection_contract.nf")]
    def names = [run:['local/src/run_fade.py','local/src/fade_json_to_tables.py','local/src/FADE_phylophere.bf','fade_run.nf'],
        merge:['local/src/fade_merge_tables.py','local/src/fade_json_to_tables.py','fade_tables.nf'],
        lists:['local/src/fade_gene_lists.py','local/src/fade_merge_tables.py','local/src/fade_json_to_tables.py','fade_gene_lists.nf'],
        report:['local/6.FADE_report.Rmd','local/src/render_fade.R','fade_report.nf']]
    assert names.containsKey(stage)
    def paths = common + names[stage].collect { file("${baseDir}/subworkflows/FADE/${it}") }
    if (stage=='run') paths += [file("${baseDir}/workflows/fade.nf"),
        file("${baseDir}/subworkflows/SELECTION/local/src/annotate_tree_fg.py"),
        file("${baseDir}/subworkflows/SELECTION/local/src/alignment_io.py"),
        file("${baseDir}/subworkflows/SELECTION/local/src/prepare_alignments.py")]
    if (stage=='report') paths += [file("${baseDir}/subworkflows/TRAIT_ANALYSIS/local/src/report_style.R")]
    selectionHash(paths.sort { it.toString() })
}

def fadeSpecies(path) {
    def rows=path.readLines().collect { it.trim() }.findAll { it && !it.startsWith('#') }
    assert rows && rows.unique(false).size()==rows.size() : "Empty or duplicate FADE species list: ${path}"
    assert rows.every { !(it =~ /[\s,]/).find() } : "Invalid FADE species identifier: ${path}"
    rows.toSet()
}

def fadeJson(value) { selectionQuote(groovy.json.JsonOutput.toJson(value)) }

def fadeSettings() {
    def get = { name, fallback -> params.containsKey(name) && params[name] != null ? params[name] : fallback }
    def s = [model:get('fade_model','LG'), method:get('fade_method','Variational-Bayes'),
        background_scope:get('fade_background_scope','all'),internal_nodes:get('fade_internal_nodes','all_descendants'),
        bf_threshold:get('fade_bf_threshold',100), grid:get('fade_grid',20), concentration:get('fade_concentration',0.5),
        chains:get('fade_chains',5),chain_length:get('fade_chain_length',2000000),burn_in:get('fade_burn_in',1000000),
        samples:get('fade_samples',1000),traitname:get('traitname','unknown_trait')]
    assert s.model in ['LG','WAG','JTT','JC69','mtMet','mtVer','gcpREV','HIVBm','HIVWm','GTR'] : "Unsupported FADE model: ${s.model}; Dayhoff is not provided by this HyPhy interface"
    assert s.method in ['Variational-Bayes','Collapsed-Gibbs','Metropolis-Hastings'] : 'Invalid FADE method'
    assert s.background_scope in ['all','opposite'] : 'Invalid FADE background scope'
    assert s.internal_nodes in ['all_descendants','none'] : 'Invalid FADE internal-node mode'
    ['bf_threshold','concentration'].each { key ->
        // params arrive as strings from params.json; accept plain decimal text only
        if (s[key] instanceof CharSequence && s[key].toString() ==~ /[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?/) s[key] = s[key].toString().toBigDecimal()
        assert s[key] instanceof Number && Double.isFinite(s[key] as double) && s[key]>0 : "Invalid FADE ${key}"
    }
    assert s.concentration >= 0.001 && s.concentration <= 1 : 'FADE concentration must be between 0.001 and 1'
    ['grid','chains','chain_length','burn_in','samples'].each { key -> s[key]=selectionPositive(s[key], "FADE ${key}") }
    assert s.grid>=5 && s.grid<=50 : 'FADE grid must be between 5 and 50'
    if (s.method!='Variational-Bayes') {
        assert s.chains>=2 && s.chains<=20 : 'FADE chains must be between 2 and 20'
        assert s.chain_length>=5000 && s.chain_length<=50000000 : 'Invalid FADE chain length'
        assert s.burn_in>=s.chain_length.intdiv(20) && s.burn_in<=s.chain_length*0.95 : 'Invalid FADE burn-in'
        assert s.samples>=50 && s.samples<=s.chain_length-s.burn_in : 'Invalid FADE sample count'
    }
    s
}
