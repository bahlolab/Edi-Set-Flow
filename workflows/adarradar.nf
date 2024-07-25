
// run checks etc here

// import subworkflows
include { M00_PREPROCESS } from '../subworkflows/local/m00_preprocess.nf'
include { M01_AGGREGATE  } from '../subworkflows/local/m01_aggregate.nf'
include { M02_INTERSECT  } from '../subworkflows/local/m02_intersect.nf'
include { M03_SAMDEPTH   } from '../subworkflows/local/m03_samdepth.nf'
include { M04_REPORT     } from '../subworkflows/local/m04_report.nf'

// TODO: check 1 row per sample
// TODO: check has bam if has jacusa
input = WfAdarRadar
    .read_csv(file(params.input), required: ['sample'])

// create channels from input
fastqs = Channel
    .fromList(input)
    .filter { it.bam == null & it.jacusa == null }
    .filter { it.fastq1 != null }
    .map { 
        [ it.sample,  [file(it.fastq1, checkIfExists: true)] + 
          (it.fastq2 == null ? [] : [file(it.fastq2, checkIfExists: true)]) ]
    }

bams_no_jacusa = Channel
    .fromList(input)
    .filter { it.bam != null & it.jacusa == null }
    .map { [it.sample, file(it.bam, checkIfExists: true), file("${it.bam}.bai", checkIfExists: true)] }

bams_with_jacusa = Channel
    .fromList(input)
    .filter { it.bam != null & it.jacusa != null }
    .map { [it.sample, file(it.bam, checkIfExists: true), file("${it.bam}.bai", checkIfExists: true)] }

jacusa_results = Channel
    .fromList(input)
    .filter { it.jacusa != null & it.bam != null }
    .map { [it.sample, file(it.jacusa, checkIfExists: true)] }

workflow ADARRADAR {

    M00_PREPROCESS(
        fastqs,
        bams_no_jacusa
    )

    M01_AGGREGATE(
        jacusa_results.mix(
            M00_PREPROCESS.out.jacusa_results
        )
    )

    M02_INTERSECT(
        M01_AGGREGATE.out.res_other,
        M01_AGGREGATE.out.samp_site_counts,
        M01_AGGREGATE.out.sites_redi_join,
        M01_AGGREGATE.out.all_site_stats,
        M01_AGGREGATE.out.bounding_ensg,
        M01_AGGREGATE.out.gen_features_intersect,
        M01_AGGREGATE.out.rm_repeats_intersect
    )

    M03_SAMDEPTH(
        bams_with_jacusa.mix(M00_PREPROCESS.out.bams),
        M01_AGGREGATE.out.res_other,
        M02_INTERSECT.out.sites_stats_filt,
        M02_INTERSECT.out.sites_filt_bed
    )

    M04_REPORT(
        // Module 1 ouputs
        M01_AGGREGATE.out.res_other,
        M01_AGGREGATE.out.samp_site_counts,
        M01_AGGREGATE.out.sites_redi_join,
        M01_AGGREGATE.out.all_site_stats,
        M01_AGGREGATE.out.bounding_ensg,
        M01_AGGREGATE.out.gen_features_intersect,
        M01_AGGREGATE.out.rm_repeats_intersect,
        // module 2 outputs
        M02_INTERSECT.out.sites_tagged_context,
        M02_INTERSECT.out.sites_stats_filt,
        M02_INTERSECT.out.sites_filt_bed,
        // module 3 outputs
        M03_SAMDEPTH.out.sample_site_depth
    )

}

workflow.onComplete { 
    log.info( 
        workflow.success ? 
            "\nOutputs written to: ${file(params.outdir).toAbsolutePath()}\n\nADAR-RADAR completed sucessfully." : 
            "\nADAR-RADAR failed."
    )
}
