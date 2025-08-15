
process GNOMAD {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 12 * task.attempt + ' h' }
    tag "chr$chr"
    storeDir "${params.resource_dir}/${params.genome}/gnomAD"

    input:
    val url_pattern
    val ver
    val chr
    val min_af

    output:
    path output

    script:
    url = url_pattern.replaceAll('<VER>', ver).replaceAll('<CHR>', chr.toString())
    output = "gnomad_v${ver}_chr${chr}.common_${min_af}.bed.bgz"
    """
    bcftools view -Ou --threads $task.cpus -i "AF_joint>$min_af" -v snps $url \\
        | bcftools query -f '%CHROM\\t%POS\\n' \\
        | uniq \\
        | awk 'BEGIN{OFS="\\t"}{print \$1, \$2-1, \$2}' \\
        | bgzip --threads $task.cpus > $output
    """
}

