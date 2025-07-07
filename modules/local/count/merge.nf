
process MERGE {
    cpus     2
    memory { 8 * task.attempt + ' GB' }
    time   { 6 * task.attempt + ' h'  }
    label  'bcftools'
    tag    "$i"
    /*
        - merge VCFs from each sample at specified region
        - add summary stats and filter
    */

    input:
    tuple val(i), val(region)
    tuple path(vcfs), path(indices)
    path targets
    val opts
    
    output:
    path(vcf)

    script:
    vcf = "jacusa.sites.part_${i}.${opts.all ? 'count': 'disco'}.vcf.gz"

    id_cmd = opts.by_id ? 
        "tabix $targets ${region.join(' ')} | cut -f4 > ids.txt" :
        ''

    filters = [] +
        (opts.filter?.site_min_af      ? ["AF>=$opts.filter.site_min_af"] : []) +
        (opts.filter?.site_min_ac      ? ["AC>=$opts.filter.site_min_ac"] : []) +
        (opts.filter?.site_min_med_vaf ? ["MEDIAN(VAF)>=$opts.filter.site_min_med_vaf"] : []) +
        (opts.filter?.site_min_pass    ? ["PASS_FREQ>=$opts.filter.site_min_pass"] : []) +
        (opts.filter?.site_min_dp_frac ? ["(COUNT((NREF+NALT)>=${opts.filter?.site_min_depth ?: 10})/N_SAMPLES)>=$opts.filter.site_min_dp_frac"] : []) +
        (params.filter_enzyme          ? ["(ADAR=1 || APOBEC=1)"] : []) +
        (opts.by_id                    ? ["ID=@ids.txt"] : [])

 
    filter_cmd = filters ?  "| bcftools view -Ou --include '${filters.join(' && ')}'" : ''

    write_cmd = "| bcftools " +
        (
            opts.keep_geno ?
            "annotate --remove ^${opts.keep_geno.collect{"FMT/$it"}.join(',')} " :
            "view --drop-genotypes "
        ) + 
        "--threads $task.cpus --write-index -Oz -o $vcf"
    
    """
    $id_cmd

    bcftools merge -Ou $vcfs \\
        --regions ${region.join(',')} \\
        --threads $task.cpus \\
        --missing-to-ref \\
        --merge id \\
        | bcftools +fill-tags -Ou -- \\
            -t 'AC,AF,PASS_FREQ:1=COUNT(PASS=1&&GT="1")/COUNT(PASS>=0&&GT="1")' \\
        $filter_cmd \\
        $write_cmd
    """
}