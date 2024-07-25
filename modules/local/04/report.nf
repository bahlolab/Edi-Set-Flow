
process REPORT {
    cpus    1
    memory "${task.attempt * 8} GB"
    time   '2 h'
    label  'R'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(rmd)
    path(res_other)
    path(sites_tagged_context)
    path(sites_stats_filt)
    path(sample_site_depth)

    output:
    path("${params.name}.report.html")

    script:
    """
    if [[ -L "$rmd" ]]; then cp --remove-destination `readlink $rmd` $rmd; fi
    04_report.R \\
        $params.name \\
        $res_other \\
        $sites_tagged_context \\
        $sites_stats_filt \\
        $sample_site_depth
    """
}