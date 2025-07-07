process JACUSA2 {
    cpus     8
    memory { 40 * task.attempt + ' GB' }
    time   { 12  * task.attempt + ' h' }
    maxRetries 2
    errorStrategy 'retry' 
    label 'jacusa2'
    tag "$sample"
    /*
        - runs JACUSA2 to detect editing sites and applies filters
        - VCF output not used as it drops strand information
        - set maxIterations=0 for speed since score is not used
        - there is a bug in multithread mode leading to stochastic dropping of the last few results
            - implement parallel with bash instead
    */
    input:
    tuple val(sample), path(bam), path(bai), val(strand), path(where)
    val opts

    output:
    tuple val(sample), path(output)
    

    script:
    output = "${sample}.jacusa.bed.gz"
    """
    mkdir -p tmp
    export JAVA_TOOL_OPTIONS="-Djava.io.tmpdir=\$PWD/tmp"
    
    zcat $where \\
        | awk -v N=$task.cpus '{ file = sprintf("bed_part_%02d", (NR-1)%N); print > file }'
    
    pids=()
    for BED in bed_part_*; do
        JACUSA2 call-1 $bam \\
            ${opts.all ? '-A' :'' } \\
            -c ${opts.min_depth ?: 5} \\
            -b \$BED \\
            -P $strand \\
            -m 20 \\
            -q $params.min_mapq \\
            -a B,I,M,S,Y \\
            -u DirMult:maxIterations=0 \\
            -p 1 \\
            -r "\$BED.jacusa.bed" &
        pids+=("\$!")
    done

    wait "\${pids[@]}"

    (   
        head -n2 bed_part_00.jacusa.bed
        cat bed_part_*.jacusa.bed \\
            | grep -v '^#'
    ) | gzip > $output

    rm bed_part_* || true
    """
}
