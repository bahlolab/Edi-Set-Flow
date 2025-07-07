
process WHERE {
    cpus   { 2                        }
    memory { 8 * task.attempt + ' GB' }
    time   { 1 * task.attempt + ' h'  }
    label 'bedops'
    tag "$sample"
    /*
        - intersect target regions with covered regions for each sample
    */

    input:
    tuple val(sample), path(sample_cov), path(sample_cov_idx)
    tuple path(targets), path(targets_idx), path(chrom_list)
    val all

    output:
    tuple val(sample), path(output)
    
    script:
    output = "${sample}.${all ? "all" : "callable"}.bed.bgz"
    
    flt_cmd = all ? "" : "| awk '/CALL\$/'"
    
    """
    while IFS= read -r CHROM; do
        echo "at \$CHROM"
        
        tabix $sample_cov \$CHROM \\
            $flt_cmd \\
            | bedops --intersect - <(tabix $targets \$CHROM ) \\
            | bgzip --threads $task.cpus \\
            >> $output \\
            || [[ \$? -eq 141 ]] # allow sigpipe from bedops closing early

    done < $chrom_list
    """
}