
process REF_GENOME {
    cpus   { 4                        }
    memory { 4 * task.attempt + ' GB' }
    time   { 4 * task.attempt + ' h'  }
    label 'samtools'
    storeDir "${params.resource_dir}"

    input:
    val url

    output:
    tuple path(output), path("${output}.*")

    script:
    read_cmd = url.startsWith('http') ? 'wget -qO-' : 'cat'
    output = url.tokenize('/').last().replace('.gz', '').replace('.bgz', '') + '.bgz'
    
    """
    $read_cmd $url \\
        | bgzip -@ $task.cpus -dc \\
        | bgzip -@ $task.cpus > $output
        
    samtools faidx -@ $task.cpus $output
    """
}