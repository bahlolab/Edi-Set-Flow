
process GTF_TO_BED {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 1 * task.attempt + ' h'  }
    label 'bedops'

    input:
    path(gtf)

    output:
    path(flat_bed), emit: flat_bed
    path(prot_cod_bed), emit: prot_cod_bed

    script:
    flat_bed     = gtf.name.replaceFirst(/\.gtf(\.(?:gz|bgz))?$/, '') + ".flat.bed.gz"
    prot_cod_bed = gtf.name.replaceFirst(/\.gtf(\.(?:gz|bgz))?$/, '') + ".prot_cod.bed"
    """
    (zcat $gtf || cat $gtf) \\
        | gtf_to_bed.awk \\
        | bedops --merge - \\
        | bgzip --threads $task.cpus \\
        > $flat_bed

    (zcat $gtf || cat $gtf) \\
        | awk 'BEGIN {FS=OFS="\\t"} \$3 == "transcript" && \$9 ~ /transcript_type "protein_coding"/' \\
        | gtf_to_bed.awk \\
        > $prot_cod_bed
    """
}