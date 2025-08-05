
process REPEATS {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label 'bcftools'
    storeDir "${params.resource_dir}/${params.genome}"
    
    input:
    val url

    output:
    tuple path(output), path("${output}.tbi")

    script:
    read_cmd = url.startsWith('http') ? 'wget -qO-' : 'cat'
    output = url.tokenize('/').last().replaceFirst(/\.gz$/, '') + '.bgz'
    
    """
    $read_cmd $url \\
        | gzip -dc \\
        | awk 'BEGIN{FS=OFS="\\t"} {print \$6, \$7, \$8, \$12, \$13, \$11}' \\
        | sort -k1,1 -k2,2n -k3,3n \\
        | bgzip --threads $task.cpus -c > $output

    tabix --threads $task.cpus -0 -p bed $output
    """
}