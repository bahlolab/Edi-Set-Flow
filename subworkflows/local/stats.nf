
include { MULTIQC } from '../../modules/local/stats/multiqc.nf'
include { SCATTER } from '../../modules/local/stats/scatter.nf'
include { REPORT } from '../../modules/local/stats/report.nf'

workflow STATS {
    take:
    edisetr_lib
    multiqc
    annotated_vcf

    main:

    MULTIQC(
        multiqc
    )

    SCATTER(
        annotated_vcf
    )

    REPORT(
        SCATTER.out,
        file(params.input, checkIfExists:true),
        file("$projectDir/bin/esf_report.qmd", checkIfExists:true),
        edisetr_lib
    )

}