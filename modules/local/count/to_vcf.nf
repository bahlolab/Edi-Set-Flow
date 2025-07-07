
process TO_VCF {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    tag    "$sample"
    label 'bcftools'
    publishDir "${params.outdir}/vcf", 
        mode: 'copy', 
        enabled: { params.publish.contains('sample_bcf') && !opts.all }

    input:
    tuple val(sample), path(jacusa_bed)
    tuple path(ref_fasta), path(ref_index)
    val opts

    output:
    tuple val(sample), path(output), path("${output}.csi")

    script:

    output = "${sample}.jacusa_counts.bcf"
    """    
    gzip -dc $jacusa_bed \\
        | jacusa_to_vcf.awk ${opts.all ? '-v ALL=1' : ''} ${params.infer_enzyme ? '-v ENZ=1' : ''} \\
        | bcftools view --threads $task.cpus -Oz -o tmp.vcf.gz
    
    bcftools reheader tmp.vcf.gz \\
        --threads $task.cpus \\
        --fai ${ref_fasta}.fai \\
        --samples <(printf "SAMPLE\t$sample\n") \\
        | bcftools sort --write-index -Ob -o $output
    
    rm tmp.vcf.gz
    """
}
