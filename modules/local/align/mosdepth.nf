
process MOSDEPTH {
    cpus     2
    memory { 4 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label  'mosdepth'
    tag    "$sample"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path("${sample}.quantized.bed.gz"), path("${sample}.quantized.bed.gz.csi"), emit: coverage
    tuple path("${sample}.mosdepth.summary.txt"), path("${sample}.mosdepth.global.dist.txt"), emit: qc

    script:
    """
    export MOSDEPTH_Q0=COV
    export MOSDEPTH_Q1=CALL
    mosdepth \\
        --fast-mode \\
        --no-per-base \\
        --quantize 1:${params.disco_filter?.site_min_depth ?: 10}: \\
        --threads $task.cpus \\
        --mapq $params.min_mapq \\
        $sample \\
        $bam
    """
}