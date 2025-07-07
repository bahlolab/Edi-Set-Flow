
process SAMTOOLS {
    cpus   { 4                         }
    memory { 32 * task.attempt + ' GB' }
    time   { 6  * task.attempt + ' h'  }
    label 'samtools'
    publishDir "${params.outdir}/bam", mode: 'copy', enabled: params.publish.contains('bam')
    tag "$sample"
    /*
        - Runs fixmate, sort, markdup and calmd
        - Drops mapQ=0 reads since never used downstream
        - Optionally drops secondary or supplementary reads
    */

    input:
    tuple val(sample), path(in_bam)
    tuple path(ref_fasta), path(ref_fasta_files)

    output:
    tuple val(sample), path(out_bam), path("${out_bam}.bai")

    script:
    out_bam = "${in_bam.name}.sorted.md.bam"
    
    filter_opt = params.drop_secondary
        ? (params.drop_supplementary ? "-F 2304"  : "-F 256")
        : (params.drop_supplementary ? "-F 2048"  : "")

    dup_opt = params.drop_duplicates ? "-r" : ""
    """
    samtools fixmate -u -@ $task.cpus -m $in_bam - \\
        | samtools view -u -q 1 $filter_opt \\
        | samtools sort -m 4G -@ $task.cpus -O BAM -o _sorted.bam
    
    samtools markdup -u -@ $task.cpus $dup_opt _sorted.bam - \\
        | samtools calmd -u -@ $task.cpus - $ref_fasta \\
        | samtools view -@ $task.cpus -b --write-index -o $out_bam##idx##${out_bam}.bai

    rm _sorted.bam
    """
}