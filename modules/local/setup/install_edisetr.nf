
process INSTALL_EDISETR {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 1 * task.attempt + ' h'  }
    label 'edisetr'
    storeDir "${params.resource_dir}"
    
    input:
    val(repo)
    val(ver)

    output:
    path(output)

    script:
    output = "edisetr-$ver-lib"
    """
    mkdir $output
    install_edisetr.R $repo@$ver $output $task.cpus
    """
}