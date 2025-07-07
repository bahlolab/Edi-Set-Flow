
process GATHER {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 4 * task.attempt + ' h'  }
    label 'bcftools'
    /*
        - Gather VCF filters
        - We can use a naive merge (much faster) since records are already
    */
    input:
    path(vcfs)

    output:
    tuple path(vcf), path("${vcf}.csi"), emit: vcf
    tuple path(bed), path("${bed}.tbi"), path('chrom_list.txt'), emit: bed

    script:
    pref = vcfs[0].name.replaceFirst(/\.part_[0-9]+/, '.gathered').replaceFirst(/\.vcf\.gz$/, '')
    vcf  = "${pref}.vcf.gz"
    bed  = "${pref}.bed.bgz"
    """
    bcftools concat $vcfs --threads $task.cpus --naive -Oz -o $vcf
    bcftools index --threads $task.cpus $vcf

    bcftools query -f '%CHROM\\t%POS\\t%POS\\t%ID\\t.\\t%STRAND\\n' $vcf \\
        | awk 'BEGIN{OFS="\\t"} {\$2=\$2-1; print}' \\
        | uniq \\
        | bgzip --threads $task.cpus > $bed
        
    tabix --threads $task.cpus -0 -p bed $bed

    bgzip -@ $task.cpus -cd $bed | cut -f1 | uniq > chrom_list.txt
    """
}