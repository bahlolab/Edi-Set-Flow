
process INSTALL_EDISETR {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 1 * task.attempt + ' h'  }
    label 'edisetr'
    tag "$hash"
    
    input:
    tuple path(path), val(hash)

    output:
    path(output)

    script:
    output = "edisetr-$hash-lib"
    """
    mkdir $output
    R -q -e \\
        "remotes::install_local(
            path = '$path',
            lib = '$output',
            force = TRUE,
            dependencies = TRUE,
            upgrade = 'never',
            build_vignettes = FALSE            
         )"
    """
}