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

// reference resources (set by -profile hg38/mm10/mm39)
params.genome          = 'generic'
params.ref_fasta_url   = null // required
params.ucsc_rmsk_url   = null // required
params.gtf_url         = null // required
params.rediportal_url  = null // optional

/* ---- regions to include (union) ---*/
params.include_gtf        = true
params.include_rediportal = true
// path to bed file with additional regions - should be gzipped
params.include_bed        = null
// expland intervals in include by this much
params.include_pad        = 100

/* ---- regions to exclude (setdiff) ---*/

// exclude common variation from gnomaAD (set by -profile hg38)
params.exclude_gnomad      = true
params.gnomad_url_pattern  = null
params.gnomad_ver          = null
params.gnomad_chr          = null
params.gnomad_min_af       = null

// exclude common variation from dbSNP (set by -profile hg38)
params.exclude_dbsnp       = true
params.dbsnp_url           = null

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
params.analysis_threads  = 16


params.report_model           = 'quasibinomial' // GLM family to run, one of 'linear', 'arcsine', 'quasibinomial', 'binomial'
params.report_fixed_effects   = null     // covariates to include, comma sep
params.report_min_med_dp      =  10      // min median depth per site
params.report_min_med_vaf     = 0.001    // min median VAF per site
params.report_max_med_vaf     = 1.0      // max median VAF per site
params.report_adar_only       = true     // only A>I sites
params.report_grp_min_med_dp  = 10       // min median depth per group per site for GLM fit
params.report_grp_min_med_vaf = 0.001    // min median VAF per group per site for GLM fit
params.report_grp_max_med_vaf = 0.999    // max median VAF per group per site for GLM fit

// force local R package install (useful for dev)
params.install_edisetr        = false         

include  { ESF } from './workflows/esf'

workflow { 

    WfEdiSetFlow.preflight_checks(params as Map)

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

