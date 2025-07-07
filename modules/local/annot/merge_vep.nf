
process MERGE_VEP {
    cpus   { 4                        }
    memory { 4 * task.attempt + ' GB' }
    time   { 4 * task.attempt + ' h'  }
    label 'bcftools'
    /*
        - Sorts inputs
        - Concatentate in a single VCF
    */
    input:
    path(vcfs)

    output:
    tuple path(vcf), path("${vcf}.csi"), emit: vcf

    script:
    vcf = vcfs[0].name.replaceFirst(/\.part_[0-9]+/, '.gathered')
    """
    printf "%s\\n" $vcfs | xargs -P$task.cpus -I{} bcftools sort {} -W -Ob -o {}.sort.bcf

    bcftools concat *.sort.bcf -a -W --threads $task.cpus -Oz -o $vcf

    rm *.sort.bcf* || true
    """
}