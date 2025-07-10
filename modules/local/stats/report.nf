
process REPORT {
    cpus   { params.analysis_threads  }
    memory { (4 + params.analysis_threads) * task.attempt + ' GB' }
    time   { 4 * task.attempt + ' h'  }
    label 'edisetr'
    publishDir "${params.outdir}/report",  mode: 'copy'
    
    input:
    path(edisites_shards)
    path(sample_metadata)
    path(report_qmd)
    path(edisetr_lib)


    output:
    path report,             emit: report
    path "$prefix.*.csv.gz", emit: data

    script:
    report = "EdiSetFlow.report.html"
    prefix = 'EdiSetFlow'

    edisites_glob = edisites_shards[0].name.replaceFirst(/shard_[0-9]+/, 'shard_*')

    qmd_params = params
        .findAll { k, _v -> k.startsWith('report_') }
        .collect { k, v ->
            "--execute-param ${k.replaceFirst(/report_/, '')}:${v instanceof String ? "\"$v\"" : "$v"}"
        }.join(" ")

    """
    export R_LIBS_USER="\$(pwd -P)/$edisetr_lib"
    mkdir HOME && export HOME="HOME"
    export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 FLEXIBLAS_NUM_THREADS=2

    quarto render $report_qmd \\
        --output $report \\
        --execute-param edisites:"$edisites_glob" \\
        --execute-param metadata:"$sample_metadata" \\
        --execute-param output:"$prefix" \\
        --execute-param threads:$task.cpus \\
        $qmd_params
    """
}