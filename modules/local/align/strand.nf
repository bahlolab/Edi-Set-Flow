
process STRAND {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 1 * task.attempt + ' h'  }
    label 'rseqc'
    tag "$sample"

    input:
    tuple val(sample), path(bam), path(bai)
    path(refseq_bed)

    output:
    tuple val(sample), path("${sample}.infer_strand.txt"), emit: strand
    path "${sample}.infer_experiment.txt", emit: qc

    script:
    """
    infer_experiment.py -i ${bam} -r ${refseq_bed} > ${sample}.infer_experiment.txt
    
    infer_strand.py ${sample}.infer_experiment.txt > ${sample}.infer_strand.txt
    """
}