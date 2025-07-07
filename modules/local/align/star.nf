
process STAR {
    cpus   { 16                        }
    memory { 40 * task.attempt + ' GB' }
    time   { 4  * task.attempt + ' h'  }
    label 'star'
    tag "$sample"

    input:
    tuple val(sample), path(fastqs)
    path star_genome_dir
    path star_gtf

    output:
    tuple val(sample), path(bam), emit: bam
    path log_out, emit: qc

    script:
    bam = "${sample}.star.bam"
    log_out = "${sample}.Log.final.out"
    """
    STAR-plain \\
        --runThreadN $task.cpus \\
        --genomeDir $star_genome_dir \\
        --sjdbGTFfile $star_gtf \\
        --readFilesIn $fastqs \\
        --readFilesCommand zcat \\
        --outFileNamePrefix ${sample}. \\
        --outSAMtype BAM Unsorted \\
        --twopassMode Basic \\
        --outFilterMismatchNoverLmax 0.10 \\
        --outSAMmapqUnique 60
    mv ${sample}.Aligned.out.bam $bam
    """
}