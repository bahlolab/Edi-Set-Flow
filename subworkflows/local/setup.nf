

include { INSTALL_EDISETR } from '../../modules/local/setup/install_edisetr.nf'
include { REF_GENOME      } from '../../modules/local/setup/ref_genome.nf'
include { REPEATS         } from '../../modules/local/setup/repeats.nf'
include { REDIPORTAL      } from '../../modules/local/setup/rediportal.nf'
include { GTF             } from '../../modules/local/setup/gtf.nf'
include { GTF_TO_BED      } from '../../modules/local/setup/gtf_to_bed.nf'
include { GNOMAD          } from '../../modules/local/setup/gnomad.nf'
include { BWAMEM2_IDX     } from '../../modules/local/setup/bwamem2_idx.nf'
include { STAR_IDX        } from '../../modules/local/setup/star_idx.nf'
include { DBSNP           } from '../../modules/local/setup/dbsnp.nf'
include { TARGET_REGIONS  } from '../../modules/local/setup/target_regions.nf'
include { SPLINTER        } from '../../modules/local/setup/splinter.nf'


workflow SETUP {
    /*
       Retreival and preprocessing of reference datasets
    */
    main:

    INSTALL_EDISETR(
        params.edisetr_repo,
        params.edisetr_ver
    )

    REF_GENOME(
        params.ref_fasta_url
    )

    REDIPORTAL(
        params.rediportal_db_url
    )
    rediportal_db = REDIPORTAL.out.db.map { db, idx, head -> [db, idx, head.text.strip().split('\t').toList() ] }

    exclude_bed = params.exclude_bed ? 
        Channel.fromPath(params.exclude_bed, checkIfExists:true).first() :
        Channel.empty()

    if (params.exclude_dbsnp) {
        DBSNP(
            params.dbsnp_url
        )
        exclude_bed = exclude_bed.concat(DBSNP.out)
    }

    if (params.exclude_gnomad) {
        GNOMAD(
            params.gnomad_url_pattern,
            params.gnomad_ver,
            Channel.fromList(params.gnomad_chr),
            params.gnomad_min_af 
        )
        exclude_bed = exclude_bed.concat(GNOMAD.out)
    }

    REPEATS(
        params.repeats_bed_url
    )
    
    GTF(
        params.gtf_url
    )

    GTF_TO_BED(
            GTF.out.uncompressed
    )

    bwamem2_index = null
    if (params.aligner == 'BWAMEM2') {
        bwamem2_index = BWAMEM2_IDX(
            REF_GENOME.out.map{ it[0] }
        )
    }
    star_index = null
    if (params.aligner == 'STAR') {
        // star needs uncompressed gtf
        STAR_IDX(
            REF_GENOME.out.map{ it[0] },
            GTF.out.uncompressed,
        )
        star_index = STAR_IDX.out.map { it[0].getParent() } // only need directory
    }

    include_bed = params.include_bed ? 
        Channel.fromPath(params.include_bed, checkIfExists:true).first() :
        Channel.empty()

    if (params.include_gtf) {
        include_bed = include_bed.concat(GTF_TO_BED.out.flat_bed)
    }

    if (params.include_rediportal) {
        include_bed = include_bed.concat(REDIPORTAL.out.bed.map{ it[0] })
    }

    TARGET_REGIONS(
        include_bed.toSortedList{ a, b -> a.name <=> b.name },
        exclude_bed.toSortedList{ a, b -> a.name <=> b.name }
    )

    SPLINTER(
        TARGET_REGIONS.out.map { it[0] }
    )

    intervals = SPLINTER.out
        .splitText()
        .map { it.strip() }
        .toList()
        .flatMap { it.withIndex() }
        .map { reg, i ->
            def width = params.n_intervals.toString().size()
            [ String.format("%0${width}d", i+1), reg.split(',') ] 
        }

    

    emit:
    ref_genome       = REF_GENOME.out
    repeats          = REPEATS.out
    rediportal_db    = rediportal_db
    rediportal_bed   = REDIPORTAL.out.bed
    gtf_tabixed      = GTF.out.tabixed
    gtf_uncompressed = GTF.out.uncompressed
    star_index       = star_index
    bwamem2_index    = bwamem2_index
    target_regions   = TARGET_REGIONS.out
    intervals        = intervals
    prot_cod_bed     = GTF_TO_BED.out.prot_cod_bed
    edisetr_lib      = INSTALL_EDISETR.out
}



