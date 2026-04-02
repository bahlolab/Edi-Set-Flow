
process FASTERQDUMP {
    cpus   { 2                        }
    memory { 4 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label  'sratools'
    tag    "$sample"

    input:
    tuple val(sample), val(accession)

    output:
    tuple val(sample), path("${sample}.R*.fq.gz")

    script:
    """
    # redirect HOME to avoid SRA config/cache being written to user home dir
    export HOME=\$(mktemp -d)

    fasterq-dump \\
        --outdir . \\
        --temp . \\
        --split-3 \\
        --min-read-len 1 \\
        --threads $task.cpus \\
        ${accession}

    if [ -f "${accession}_1.fastq" ]; then
        gzip ${accession}_1.fastq &
        gzip ${accession}_2.fastq
        wait
        mv ${accession}_1.fastq.gz ${sample}.R1.fq.gz
        mv ${accession}_2.fastq.gz ${sample}.R2.fq.gz
    else
        gzip ${accession}.fastq
        mv ${accession}.fastq.gz ${sample}.R1.fq.gz
    fi
    """
}
