
process SPLINTER {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label 'bedops'
    
    input:
    path in_bed

    output:
    path 'intervals.txt'

    script:
    """
    bgzip --threads $task.cpus -cd $in_bed \\
        | bedops --range 1:1001  --merge - \\
        | bedops --range 0:-1000 --merge - \\
        | splinter.awk -v N=$params.n_intervals \\
        > intervals.txt
    """
}