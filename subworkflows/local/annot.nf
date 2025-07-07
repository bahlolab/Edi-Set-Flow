
include { VEP       } from '../../modules/local/annot/vep.nf'
include { MERGE_VEP } from '../../modules/local/annot/merge_vep.nf'
include { VCFANNO   } from '../../modules/local/annot/vcfanno.nf'

workflow ANNOT {
    take:
    ref_genome
    gtf_tabixed
    rediportal_db
    repeats
    intervals
    sites_vcf

    main:
    
    VEP(
        sites_vcf,
        intervals,
        ref_genome,
        gtf_tabixed
    )

    MERGE_VEP(
        VEP.out.toSortedList{ a, b -> a.name <=> b.name }
    )

    VCFANNO(
        MERGE_VEP.out,
        rediportal_db,
        repeats
    )

    emit:
    annotated_vcf = VCFANNO.out
}