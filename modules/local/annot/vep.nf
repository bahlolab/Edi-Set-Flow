def vep_cmd(gtf, fasta, cpus){
    """vep \\
            --fork $cpus \\
            --format vcf \\
            --vcf \\
            --dont_skip \\
            --numbers \\
            --symbol \\
            --no_stats \\
            --pick \\
            --gtf $gtf \\
            --fasta $fasta \\
            --output_file stdout"""
}

process VEP {
    cpus   { 4                        }
    memory { 8 * task.attempt + ' GB' }
    time   { 6 * task.attempt + ' h'  }
    label 'vep'
    tag   "$i"
    /*
        -annotate VCF using VEP, but match strand of annotations
        -split gtf into a forward strand and a reverse strand
    */
    input:
    tuple path(vcf), path(vcf_index)
    tuple val(i), val(region)
    tuple path(ref_fasta), path(ref_fasta_index)
    tuple path(both_gtf), path(both_index), path(plus_gtf), path(plus_index), path(minus_gtf), path(minus_index)

    output:
    path(vep)

    script:
    pref = vcf.name.replaceFirst(/\.vcf\.gz$/, ".part_${i}")
    vep  = pref + '.vep.vcf.gz'
    """
    tabix -h $vcf ${region.join(' ')} \\
        | awk -F'\t' '/^#/ {print; next} \$8 ~ /STRAND=\\./' \\
        | ${vep_cmd(both_gtf, ref_fasta, task.cpus)} \\
        | bgzip --threads $task.cpus \\
        > $vep

    tabix -h $vcf ${region.join(' ')} \\
        | awk -F'\t' '/^#/ {print; next} \$8 ~ /STRAND=\\+/' \\
        | ${vep_cmd(plus_gtf, ref_fasta, task.cpus)} \\
        | awk '!/^#/' \\
        | bgzip --threads $task.cpus \\
        >> $vep

    tabix -h $vcf ${region.join(' ')} \\
        | awk -F'\t' '/^#/ {print; next} \$8 ~ /STRAND=-/' \\
        | ${vep_cmd(minus_gtf, ref_fasta, task.cpus)} \\
        | awk '!/^#/' \\
        | bgzip --threads $task.cpus \\
        >> $vep
    """
}
