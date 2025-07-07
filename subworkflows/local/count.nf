

include { WHERE    } from '../../modules/local/count/where.nf'
include { JACUSA2  } from '../../modules/local/count/jacusa2.nf'
include { TO_VCF   } from '../../modules/local/count/to_vcf.nf'
include { MERGE    } from '../../modules/local/count/merge.nf'
include { GATHER   } from '../../modules/local/count/gather.nf'

workflow COUNT {
    /* 
        - Discover editing sites with jacusa2
        - Convert to VCF
        - Merge and filter results
        - Emit set of candidate sites
    */
    take: 
    ref_fasta
    intervals
    target_regions // bed, index, chrom_list
    sample_bams // sample, bam, bai, strand, cov, cov_idx
    opts

    main:

    WHERE(
        sample_bams.map { it[0, 4, 5] }, // sample, cov, cov_idx
        target_regions,
        opts.all ? true : false
    )

    jacusa_input = sample_bams
        .map { it.take(4) } // sample, bam, bai, strand,
        .join(WHERE.out, by:0) // sample, bam, bai, strand, where

    JACUSA2(
        jacusa_input,
        opts.subMap(['all', 'min_depth']),
    )

    TO_VCF(
        JACUSA2.out,
        ref_fasta,
        opts.subMap(['all'])
    )

    vcf_idx_channel = TO_VCF.out
        .map{ it.drop(1) }
        .toSortedList()
        .map{ it.transpose() }
    

    MERGE(
        intervals,
        vcf_idx_channel,
        target_regions,
        opts.subMap(['filter', 'all', 'keep_geno', 'by_id'])
    )

    GATHER(
        MERGE.out.toSortedList{ a, b -> a.name <=> b.name }
    )

    emit:
    sites_vcf = GATHER.out.vcf
    sites_bed = GATHER.out.bed
}