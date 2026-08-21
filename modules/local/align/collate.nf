
process COLLATE {
    cpus   { 4                         }
    memory { 8  * task.attempt + ' GB' }
    time   { 4  * task.attempt + ' h'  }
    label 'samtools'
    tag "$sample"
    /*
        - Name-collates a pre-aligned (typically coordinate-sorted) BAM so it
          matches the read-grouped ordering emitted by STAR/BWAMEM2, which the
          downstream SAMTOOLS `fixmate -m` step requires.
    */

    input:
    tuple val(sample), path(in_bam)

    output:
    tuple val(sample), path(out_bam)

    script:
    out_bam = "${sample}.collated.bam"
    """
    samtools collate -@ $task.cpus -o $out_bam $in_bam
    """
}
