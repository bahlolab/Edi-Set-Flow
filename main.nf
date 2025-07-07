#!/usr/bin/env nextflow

// TODO - mode to disable discovery and just use REDIPORTAL catalgue
// 
// params.mode = 'DISCOVERY' // either DISCOVERY, CATALOGUE, or BOTH

/* ------------------ INPUT PARAMS ------------------*/
// input sample manifest in csv format: sample_id, fastq1, fastq2
params.input  = null
// directory to store pipeline outputs
params.outdir = 'output'
// what outputs to publish in outdir
params.publish = ['bam', 'sample_bcf', 'annot_vcf', 'report', 'diff_sites']
// run setup to download and generate required resources first
params.resource_dir      = 'esf_resources'
// reference resources
params.ref_fasta_url     = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz'
params.repeats_bed_url   = 'https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38/hg38_2020_rmsk.bed'
params.gtf_url           = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz'
params.rediportal_db_url = 'http://srv00.recas.ba.infn.it/webshare/ATLAS/download/TABLE1_hg38_v3.txt.gz'

/* ---- regions to include (union) ---*/
params.include_gtf        = true
params.include_rediportal = true
// path to bed file with additional regions - should be gzipped
params.include_bed        = null
// expland intervals in include by this much
params.include_pad        = 100

/* ---- regions to exclude (setdiff) ---*/
// exclude common variation from gnomaAD (v4.1)
params.exclude_gnomad      = true
params.gnomad_url_pattern  = "https://storage.googleapis.com/gcp-public-data--gnomad/release/<VER>/vcf/joint/gnomad.joint.v<VER>.sites.chr<CHR>.vcf.bgz"
params.gnomad_ver          = '4.1'
params.gnomad_chr          = (1..22) + ['X','Y']
params.gnomad_min_af       = 0.001
// exclude common variation from dbSNP(v153)
params.exclude_dbsnp       = true
params.dbsnp_url           = 'http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp153Common.bb'
// exclude custom regions - should be gzipped
params.exclude_bed         = null

/* ------------- Preproc/alignment -------------------*/
// whether or not to preprocess fastqs with fastp
// either null (disabled), 'MINIMAL' or 'STRICT'
params.fastp_mode = null
// Aligner of choice, either STAR or BWAMEM2 (if using BWAMEM2, index files must exist at refbase)
params.aligner = 'BWAMEM2'
// bam filters
params.drop_secondary     = true
params.drop_supplementary = false
params.drop_duplicates    = true


/* ----------- Edi-site detection ------------------*/
// min mapq used by JACUSA2
params.min_mapq = 20

// TODO: set true if analyse mix of stranded and unstranded data, else will error
// params.force_unstranded   = false

// infer enzyme and strand based on ADAR=A>G, APOBEC=C>T
params.infer_enzyme  = true
// remove sites with no matching enzyme pattern
params.filter_enzyme = false

// site discovery filters
params.disco_filter = [
    site_min_depth   : 10  ,  // minimum depth for site discovery with JACUSA2
    site_min_af      : 0.10,  // minimum frequency of site detection
    site_min_ac      : 3   ,  // minimum number of samples with site
    site_min_med_vaf : 0.01,  // min median variant allele-freq in samples with alt allele
    site_min_pass    : 0.66,  // min proportion of detected samples that pass JACUSA site filter
]
// site counting filters
params.count_filter = [
    site_min_depth   : 10   ,  // minimum depth for site_min_dp_frac filter
    site_min_dp_frac : 0.20 ,  // minimum number of sites with this fraction of coverage 
    site_min_af      : 0.20 ,  // minimum frequency of site detection
    site_min_ac      : 5    ,  // minimum number of samples with site
    site_min_med_vaf : 0.0  ,  // min median variant allele-freq in samples with alt allele
    site_min_pass    : 0.66 ,  // min proportion of detected samples that pass JACUSA site filter
]

// TODO  implement forced unstranded analysis
// params.ignore_strand = false

/* ----------- paralellism controls ---------- */
// number of genomic intervals to parallelise across for site merging and VEP
params.n_intervals       = 10
// number of threads for final report / GLM fitting
params.analysis_threads  = 32

// dev params
params.edisetr_repo = "bahlolab/Edi-Set-Flow/edisetr"
params.edisetr_ver  = "25.07-beta.1"

params.report_params = [
    model: 'quasibinomial', // GLM family to run, one of 'linear', 'arcsine', 'quasibinomial', 'binomial'
    fixed_effects: 'sex,age,death', // covariates to include, comma sep
    min_med_dp:  10,        // min median depth per site
    min_med_vaf: 0.001,     // min median VAF per site
    max_med_vaf: 1.0,       // max median VAF per site
    adar_only: true,        // only A>I sites
    grp_min_med_dp: 10,     // min median depth per group per site for GLM fit
    grp_min_med_vaf: 0.001, // min median VAF per group per site for GLM fit
    grp_max_med_vaf: 0.999 // max median VAF per group per site for GLM fit
]

include  { ESF } from './workflows/esf'

workflow { 
    log.info(
"""\
\u001B[1;34m ___ ___ ___  \u001B[1;32m ___ ___ _____  \u001B[1;33m ___ _    _____      __
\u001B[1;34m| __|   \\_ _|\u001B[1;32m / __| __|_   _|\u001B[1;33m | __| |  / _ \\ \\    / /
\u001B[1;34m| _|| |) | |  \u001B[1;32m\\__ \\ _|  | | \u001B[1;33m  | _|| |_| (_) \\ \\/\\/ / 
\u001B[1;34m|___|___/___| \u001B[1;32m|___/___| |_|   \u001B[1;33m|_| |____\\___/ \\_/\\_/\u001B[0m\
"""
    )

    ESF()
}

