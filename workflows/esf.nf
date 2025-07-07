
// run checks etc here

/* ---------------- IMPORT WORKFLOWS ---------------- */
include { SETUP          } from '../subworkflows/local/setup.nf'
include { ALIGN          } from '../subworkflows/local/align.nf'
include { COUNT as DISCO } from '../subworkflows/local/count.nf'
include { COUNT          } from '../subworkflows/local/count.nf'
include { ANNOT          } from '../subworkflows/local/annot.nf'
include { STATS          } from '../subworkflows/local/stats.nf'

workflow ESF {

    /* ----------------- READ INPUTS ------------------- */

    // TODO: check 1 row per sample
    input = WfEdiSetFlow.read_csv(file(params.input, checkIfExists:true), required: ['sample_id'])

    fastqs = Channel
        .fromList(input)
        .filter { it.bam == null } // exclude samples which already have bams
        .filter { it.fastq1 != null } // TODO - raise error if sample has no inputs
        .map { 
            [ it.sample_id,  [file(it.fastq1, checkIfExists: true)] + 
              (it.fastq2 == null ? [] : [file(it.fastq2, checkIfExists: true)]) ]
        }

    // bams = Channel
    //     .fromList(input)
    //     .filter { it.bam != null }
    //     .map { [it.sample_id, file(it.bam, checkIfExists: true), file("${it.bam}.bai", checkIfExists: true)] }

    /* ----------------- RUN SUBWORKFLOWS ------------------- */
    SETUP()

    ALIGN(
        fastqs,
        SETUP.out.ref_genome,
        SETUP.out.gtf_uncompressed,
        SETUP.out.bwamem2_index,
        SETUP.out.star_index,
        SETUP.out.prot_cod_bed
    )

    DISCO(
        SETUP.out.ref_genome,
        SETUP.out.intervals,
        SETUP.out.target_regions,
        ALIGN.out.sample_bams,
        [
            all      : false,
            by_id    : false,
            filter   : params.disco_filter,
            keep_geno: false,
            min_depth: params.disco_filter?.site_min_depth ?: 10
        ]
    )

    COUNT(
        SETUP.out.ref_genome,
        SETUP.out.intervals,
        DISCO.out.sites_bed,
        ALIGN.out.sample_bams,
        [
            all      : true,
            by_id    : true,
            filter   : params.count_filter,
            keep_geno: ['NALT', 'NREF'],
            min_depth: 1
        ]
    )

    ANNOT(
        SETUP.out.ref_genome,
        SETUP.out.gtf_tabixed,
        SETUP.out.rediportal_db,
        SETUP.out.repeats,
        SETUP.out.intervals,
        COUNT.out.sites_vcf
    )

    multiqc = ALIGN.out.multiqc.flatten().toSortedList()

    STATS(
        SETUP.out.edisetr_lib,
        multiqc,
        ANNOT.out.annotated_vcf
    )

    workflow.onComplete = { 
        log.info( 
            workflow.success ? 
                "\nOutputs written to: ${file(params.outdir).toAbsolutePath()}\nEdiSetFlow completed sucessfully" : 
                "\nEdiSetFlow failed"
        )
    }
}


