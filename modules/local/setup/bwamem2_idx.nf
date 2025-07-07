
process BWAMEM2_IDX {
    cpus   { 2                         }
    memory { 96 * task.attempt + ' GB' }
    time   { 8  * task.attempt + ' h'  }
    label 'bwamem2'
    storeDir "${params.resource_dir}/$genome_dir"

    input:
    path(ref_fasta)

    output:
    path(outputs)

    script:
    base = ref_fasta.name.replaceFirst(/\.(?:fa|fasta)(\.(?:gz|bgz))?$/, '')
    genome_dir = 'bwamem2_' + base    
    outputs = 
        ['0123', 'amb', 'ann', 'bwt.2bit.64', 'pac']
        .collect { "${ref_fasta.name}.${it}" }
    """
    bwa-mem2 index $ref_fasta
    """
}