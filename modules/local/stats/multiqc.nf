
process MULTIQC {
    cpus 2
    memory '4 GB'
    time '1 h'
    label 'multiqc'
    publishDir "${params.outdir}/report", mode: 'copy', pattern: "*.html"
    errorStrategy 'ignore'

    input:
        path(files)

    output:
        tuple path("${output}_data/multiqc_general_stats.txt"),
              path("${output}.html")

    script:
    output = "multiqc_report"
    """
    multiqc . \\
        --title "Edi-Set-Flow MultiQC Report" \\
        --filename ${output}.html \\
        --module fastp \\
        --module star \\
        --module rseqc \\
        --module mosdepth \\
        --module vep
    """
    // --module samtools
}


