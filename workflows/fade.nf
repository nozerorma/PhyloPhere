#!/usr/bin/env nextflow
// FADE: one audited contract for individual genes and batches. Inputs come from
// SELECTION_PREP; zero-based alignment coordinates are retained for downstream.
include { FADE_RUN; FADE_BATCHED } from '../subworkflows/FADE/fade_run'
include { FADE_TABLES as FADE_TABLES_TOP; FADE_TABLES as FADE_TABLES_BOTTOM } from '../subworkflows/FADE/fade_tables'
include { FADE_REPORT as FADE_REPORT_TOP; FADE_REPORT as FADE_REPORT_BOTTOM } from '../subworkflows/FADE/fade_report'
include { FADE_GENE_LISTS as FADE_GENE_LISTS_TOP; FADE_GENE_LISTS as FADE_GENE_LISTS_BOTTOM } from '../subworkflows/FADE/fade_gene_lists'
include { fadeSources; fadeSettings; fadeSpecies } from '../subworkflows/FADE/fade_contract'
include { selectionHash; selectionSingle; selectionPositive } from '../subworkflows/SELECTION/selection_contract'

def fadeBatch(rows, id) {
    rows.each { row ->
        assert row[0] && row[0] == row[0].trim() && !(row[0] =~ /[\t\r\n\/\\]/).find() : 'Invalid FADE gene identifier'
        assert [row[2],row[3],row[4],row[5]].findAll { it!=null }.every { !(it.name =~ /[\t\r\n]/).find() } : 'Invalid FADE input filename'
    }
    def fastas=rows.collect { it[2] }.unique()
    def trees=rows.collect { it[3] }.unique()
    def species=rows.collectMany { [it[4],it[5]] }.findAll { it!=null }.unique()
    [fastas,trees,species].each { paths ->
        assert paths.collect { it.name }.unique(false).size()==paths.size() : 'FADE staging filename collision'
    }
    def manifest=rows.collect { row ->
        [row[0],row[1],"fastas/${row[2].name}","trees/${row[3].name}","species/${row[4].name}",row[5]!=null ? "species/${row[5].name}" : 'NO_FILE'].join('\t')
    }.join('\n')+'\n'
    tuple(id,manifest,fastas,trees,species,selectionHash([*fastas,*trees,*species]))
}

def split_gene_lists(ch) {
    [bg:ch.flatten().filter { it.name=='background.txt' }.collect(),interest:ch.flatten().filter { it.name!='background.txt' }.collect()]
}

workflow FADE {
    take:
    fasta_top_ch
    fasta_bottom_ch
    tree_ch
    top_species_ch
    bottom_species_ch

    main:
    assert params.fade_direction in ['top','bottom','both'] : 'Invalid fade_direction'
    def wanted=params.fade_direction=='both' ? ['top','bottom'] : [params.fade_direction]
    def settings=fadeSettings()
    def source_hash=fadeSources('run')
    def merge_hash=fadeSources('merge')
    def batchSize=selectionPositive(params.fade_batch_size!=null ? params.fade_batch_size : 1,'fade_batch_size')
    def tree=tree_ch.collect().map { selectionSingle(it,'FADE tree') }
    def top=top_species_ch.collect().map { selectionSingle(it,'FADE top species') }
    def bottom=bottom_species_ch.collect().map { selectionSingle(it,'FADE bottom species') }
    def pools=top.combine(bottom).map { fg,bg ->
        assert fadeSpecies(fg).intersect(fadeSpecies(bg)).isEmpty() : 'FADE species pools overlap'
        [fg,bg]
    }
    def topRows=fasta_top_ch.map { gid,f -> tuple(gid,'top',f) }
    def bottomRows=fasta_bottom_ch.map { gid,f -> tuple(gid,'bottom',f) }
    def rows=topRows.mix(bottomRows).filter { wanted.contains(it[1]) }.combine(tree).combine(pools)
        .map { gid,d,f,t,fg,bg -> tuple(gid,d,f,t,d=='top' ? fg : bg,settings.background_scope=='opposite' ? (d=='top' ? bg : fg) : null) }
    def expected=rows.map { gid,d,f,t,fg,bg -> [gid,d] }.collect(flat:false).ifEmpty([]).map { keys ->
        assert keys.unique(false).size()==keys.size() : 'Duplicate FADE gene/direction inputs'
        keys.sort { a,b -> a[0]<=>b[0] ?: a[1]<=>b[1] }
    }
    def genes
    def sites
    def audits
    if (batchSize>1) {
        def batches=rows.groupTuple(by:0).toSortedList { a,b -> a[0]<=>b[0] }.flatMap().collate(batchSize).map { batch ->
            def flat=[]
            batch.each { g -> (0..<g[1].size()).each { i -> flat << [g[0],g[1][i],g[2][i],g[3][i],g[4][i],g[5][i]] } }
            flat.sort { a,b -> a[0]<=>b[0] ?: a[1]<=>b[1] }
            fadeBatch(flat,'batch_'+java.security.MessageDigest.getInstance('SHA-256').digest(flat.collect { it[0]+':'+it[1] }.join('\n').bytes).encodeHex().toString().take(16))
        }
        def result=FADE_BATCHED(batches,settings,source_hash)
        genes=result.genes;sites=result.sites;audits=result.audits.flatten()
    } else {
        def single=rows.map { row ->
            def b=fadeBatch([row],'single')
            tuple(row[0],row[1],b[1],b[2],b[3],b[4],b[5])
        }
        def result=FADE_RUN(single,settings,source_hash)
        genes=result.genes;sites=result.sites;audits=result.audits.flatten()
    }
    def allGenes=genes.collect().ifEmpty([]).map { paths -> paths.sort { it.toString() } }
    def allSites=sites.collect().ifEmpty([]).map { paths -> paths.sort { it.toString() } }
    def allAudits=audits.collect().ifEmpty([]).map { paths -> paths.sort { it.toString() } }
    def report_top_ch=Channel.empty(), report_bottom_ch=Channel.empty()
    def summary_top_ch=Channel.empty(), summary_bottom_ch=Channel.empty()
    def site_top_ch=Channel.empty(), site_bottom_ch=Channel.empty()
    def csv_top_ch=Channel.empty(), csv_bottom_ch=Channel.empty()
    def bg_top_ch=Channel.empty(), bg_bottom_ch=Channel.empty()
    def sig_top_ch=Channel.empty(), sig_bottom_ch=Channel.empty()
    if (wanted.contains('top')) {
        def tables=FADE_TABLES_TOP('top',allGenes,allSites,allAudits,expected,settings,merge_hash)
        summary_top_ch=tables.summary_tsv;site_top_ch=tables.site_tsv;csv_top_ch=tables.sites_csv
        report_top_ch=FADE_REPORT_TOP('top',tables.gene_counts_tsv,tables.sites_csv,top,tables.metadata).report
        def lists=split_gene_lists(FADE_GENE_LISTS_TOP('top',tables.summary_tsv).gene_lists)
        bg_top_ch=lists.bg;sig_top_ch=lists.interest
    }
    if (wanted.contains('bottom')) {
        def tables=FADE_TABLES_BOTTOM('bottom',allGenes,allSites,allAudits,expected,settings,merge_hash)
        summary_bottom_ch=tables.summary_tsv;site_bottom_ch=tables.site_tsv;csv_bottom_ch=tables.sites_csv
        report_bottom_ch=FADE_REPORT_BOTTOM('bottom',tables.gene_counts_tsv,tables.sites_csv,bottom,tables.metadata).report
        def lists=split_gene_lists(FADE_GENE_LISTS_BOTTOM('bottom',tables.summary_tsv).gene_lists)
        bg_bottom_ch=lists.bg;sig_bottom_ch=lists.interest
    }
    emit:
    report_top=report_top_ch
    report_bottom=report_bottom_ch
    summary_top=summary_top_ch
    summary_bottom=summary_bottom_ch
    site_tsv_top=site_top_ch
    site_tsv_bottom=site_bottom_ch
    sites_csv_top=csv_top_ch
    sites_csv_bottom=csv_bottom_ch
    gene_lists_bg_top=bg_top_ch
    gene_lists_bg_bottom=bg_bottom_ch
    gene_lists_sig_top=sig_top_ch
    gene_lists_sig_bottom=sig_bottom_ch
}
