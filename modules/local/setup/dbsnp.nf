
process DBSNP {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label 'bigbedtobed'
    storeDir "${params.resource_dir}"
    
    input:
    val url

    output:
    path(output)

    script:
    output = url.tokenize('/').last().replace('.bb', '') + '.bed.bgz'
    """
    bigBedToBed $url stdout \\
        | awk 'BEGIN { FS = OFS = "\\t" } \$14 == "snv" {print \$1, \$2, \$3 }' \\
        | gzip > $output
    """
}