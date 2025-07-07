
process BWAMEM2 {
    cpus   { 10                        }
    memory { 40 * task.attempt + ' GB' }
    time   { 4  * task.attempt + ' h'  }
    label 'bwamem2'
    tag "$sample"

    input:
    tuple val(sample), path(fastqs)
    path bwa_ref
    path bwa_ref_files

    output:
    tuple val(sample), path(bam)

    script:
    bam = "${sample}.bwamem2.bam"
    """
    bwa-mem2 mem -t ${task.cpus} -M $bwa_ref $fastqs \\
        | samtools view --threads 2 -bS -o $bam -
    """
}