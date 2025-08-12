
process GTF {
    cpus   { 2                        }
    memory { 4 * task.attempt + ' GB' }
    time   { 4 * task.attempt + ' h'  }
    label 'bcftools'
    storeDir "${params.resource_dir}/${params.genome}"
    /*
        - retrieve and indexe GTF
        - split into +/- strands for use with VEP in stranded modes
    */

    input:
    val url

    output:
    tuple path(both), path("${both}.tbi"), 
          path(plus), path("${plus}.tbi"), 
          path(minus), path("${minus}.tbi"),
          emit: tabixed
    path uncmp, emit: uncompressed

    script:
    read_cmd = url.startsWith('http') ? 'wget -qO-' : 'cat'
    pref  = url.tokenize('/').last().replaceFirst(/\.gtf(\.b?gz)/, '')
    both  = pref + '.gtf.bgz'
    plus  = pref + '.plus.gtf.bgz'
    minus = pref + '.minus.gtf.bgz'
    uncmp = pref + '.gtf'
    """
    set +o pipefail

    $read_cmd $url \\
        | gzip -dc \\
        | head -n10 \\
        | grep '^#' > header

    set -o pipefail

    (
        cat header
        $read_cmd $url \\
            | gzip -dc \\
            | grep -v '^#' \\
            | sort -k1,1 -k4,4n -k5,5n
    ) | bgzip --threads $task.cpus -c > $both \\
        && tabix --threads $task.cpus -p gff $both

    bgzip -cd --threads $task.cpus $both \\
        | awk -F'\\t' '/^#/ { print; next } \$7=="+"' \\
        | bgzip --threads $task.cpus \\
        > $plus \\
        && tabix --threads $task.cpus -p gff $plus

    bgzip -cd --threads $task.cpus $both \\
        | awk -F'\\t' '/^#/ { print; next } \$7=="-"' \\
        | bgzip --threads $task.cpus \\
        > $minus \\
        && tabix --threads $task.cpus -p gff $minus

    bgzip -cd --threads $task.cpus $both > $uncmp
    """
}