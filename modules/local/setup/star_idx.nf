
process STAR_IDX {
    cpus   { 16                        }
    memory { 48 * task.attempt + ' GB' }
    time   { 8  * task.attempt + ' h'  }
    label 'star'
    storeDir "${params.resource_dir}/${params.genome}"

    input:
    path ref_fasta
    path ref_gtf

    output:
    path outputs

    script:
    base = ref_fasta.name.replaceFirst(/\.(?:fa|fasta)(\.(?:gz|bgz))?$/, '')
    genome_dir = 'star_' + base
    
    outputs = [
        'chrLength.txt',
        'exonGeTrInfo.tab',
        'genomeParameters.txt',
        'sjdbInfo.txt',
        'chrNameLength.txt',
        'exonInfo.tab',
        'Log.out',
        'sjdbList.fromGTF.out.tab',
        'chrName.txt',
        'geneInfo.tab',
        'SA',
        'sjdbList.out.tab',
        'chrStart.txt',
        'Genome',
        'SAindex',
        'transcriptInfo.tab'
    ].collect { "${genome_dir}/$it" }

    cmd1 =  ref_fasta.name ==~ /.*(\.gz|\.bgz|\.gzip)$/ ?
        "gzip -dc $ref_fasta > ref_fasta" :
        "ln -s $ref_fasta ref_fasta"
    """
    $cmd1

    STAR-plain \\
        --runThreadN $task.cpus \\
        --runMode genomeGenerate \\
        --genomeDir $genome_dir \\
        --genomeFastaFiles ref_fasta \\
        --sjdbGTFfile $ref_gtf
    """
}