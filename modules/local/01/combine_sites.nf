
process COMBINE_SITES {
    cpus    1
    memory "${task.attempt * 8} GB" 
    time   '2 h'
    label  'R'
    publishDir "${params.outdir}/rds", pattern: '*.rds', mode: 'copy'

    input:
    path(jacusa_tables)
    path(redi_counts)

    output:
    path "${params.name}.stranded_edSites.bed.gz", emit: edsites_bed
    path "${params.name}.res_other.rds",          emit: res_other
    path "${params.name}.nSamp_site_counts.rds",  emit: samp_site_counts
    path "${params.name}.siteStats_rediJOIN.rds", emit: sites_redi_join
    path "${params.name}.all_siteStats.rds",      emit: all_site_stats

    script:
    """
    01_combine_sites.R \\
        $params.name \\
        $redi_counts \\
        $params.sample_thresh \\
        $params.alt_count_thresh \\
        $params.depth_thresh \\
        $jacusa_tables
    """
}