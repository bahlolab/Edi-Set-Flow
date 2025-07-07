include { FASTP    } from '../../modules/local/align/fastp.nf'
include { STAR     } from '../../modules/local/align/star.nf'
include { BWAMEM2  } from '../../modules/local/align/bwamem2.nf'
include { SAMTOOLS } from '../../modules/local/align/samtools.nf'
include { MOSDEPTH } from '../../modules/local/align/mosdepth.nf'
include { STRAND   } from '../../modules/local/align/strand.nf'

workflow ALIGN {
    take:
    fastqs
    ref_fasta
    ref_gtf
    bwamem2_index
    star_index
    infer_strand_bed
    
    main:
    // multiqc
    multiqc = Channel.empty()

    if (params.fastp_mode) {
        FASTP(fastqs)
        fastqs = FASTP.out.fastq
        multiqc = multiqc.concat(FASTP.out.qc)
    }
    
    if (params.aligner == 'STAR') {
        STAR(
            fastqs,
            star_index,
            ref_gtf
        )

        aligned_bams = STAR.out.bam
        multiqc = multiqc.concat(STAR.out.qc)

    } else {
        BWAMEM2(
            fastqs,
            ref_fasta.map {it[0]},
            bwamem2_index

        )
        aligned_bams = BWAMEM2.out
    }

    SAMTOOLS(
        aligned_bams,
        ref_fasta
    )

    bamdir = file("${params.outdir}/bam").with{ it.mkdirs(); it}.toRealPath()
    

    MOSDEPTH(
        SAMTOOLS.out
    )
    multiqc = multiqc.concat(MOSDEPTH.out.qc)

    STRAND(
        SAMTOOLS.out,
        infer_strand_bed
    )
    multiqc = multiqc.concat(STRAND.out.qc)
    sample_strand = STRAND.out.strand.map { [it[0], it[1].text.strip()] } //sample,strand

    sample_bams = 
        SAMTOOLS.out
        .join(sample_strand, by: 0) //sample, bam, bai, strand
        .join(MOSDEPTH.out.coverage, by:0) // sample, bam, bai, strand, cov, cov_idx

    // message about strand inference
    sample_bams
        .map { sam, _bam, _bai, strand, _cov, _idx  -> [strand, sam ] }
        .groupTuple()
        .map { strand, list -> 
            "INFO: Inferred ${list.size()} experiment${list.size() > 1 ? 's' : ''} as $strand" 
        }

        // create bam manifest
    sample_bams
        .map {sample, bam, _bai, strand, _cov, _idx -> "$sample,$strand,$bamdir/${bam.fileName}" }
        .collectFile(
            seed:     'sample,strand,bam',
            name:     "bams.csv",
            storeDir: params.outdir,
            newLine:  true,
            sort:     true, 
            cache:    false
        )

    // check not mixed strands
    sample_bams
        .map { _sam, _bam, _bai, strand, _cov, _idx  -> strand }
        .unique()
        .toList()
        .map { strands ->
            if ( strands.any { it =~ /UNSTRANDED/ } && strands.any { it =~ /(?:FIRST|SECOND)STRAND/ }) {
                sleep(1)
                error("ERROR: Mixture of stranded and unstranded reads detected. Set params.force_unstranded to proceed")
            }
        }

    sample_bams = sample_bams
        .map { sam, bam, bai, strand, cov, cov_idx  ->
            [sam, bam, bai, strand.replaceAll('^(PAIRED|SINGLE)-', ''), cov, cov_idx ] 
        }

    emit:
    sample_bams = sample_bams // sample, bam, bai, strand, cov, cov_idx
    multiqc     = multiqc
}