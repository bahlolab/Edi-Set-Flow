
process REDIPORTAL {
    cpus   { 2                        }
    memory { 4 * task.attempt + ' GB' }
    time   { 4 * task.attempt + ' h'  }
    label 'bcftools'
    storeDir "${params.resource_dir}"
    
    input:
    val url

    output:
    tuple path(output), path("${output}.tbi"), path(header)    , emit: db
    tuple path(bed), path("${bed}.tbi"), path(chroms), emit: bed

    script:
    read_cmd = url.startsWith('http') ? 'wget -qO-' : 'cat'
    pref     = "REDIPORTAL." + url.tokenize('/').last().replace('.gz', '').replace('.bgz', '') 
    output   = pref + '.bgz'
    header   = pref + '.header'
    bed      = pref.replaceFirst(/\.txt$/, '') + '.bed.bgz'
    chroms   = pref.replaceFirst(/\.txt$/, '') + '.chrom_list'

    """
    set +o pipefail

    $read_cmd $url \\
        | gzip -dc \\
        | head -n1 \\
        | awk -F'\\t' 'BEGIN{OFS=FS}{print \$2,\$3,\$4,\$5,\$1,\$6,\$7,\$8,\$9,\$10,\$11,\$12}' \\
        | sed \$'s:\\tEd\\t:\\tAlt\\t:' > $header

    (
        echo -n '#'
        cat $header
        $read_cmd $url \\
            | bgzip -@ $task.cpus -cd \\
            | tail -n+2 \\
            | awk -F'\\t' 'BEGIN{OFS=FS}{print \$2,\$3,\$4,\$5,\$1,\$6,\$7,\$8,\$9,\$10,\$11,\$12}'
    ) | bgzip -@ $task.cpus > $output

    set -o pipefail
    
    tabix -@ $task.cpus -s1 -b2 -e2 -c '#' $output

    bgzip -@ $task.cpus -cd $output \\
        | tail -n+2 \\
        | awk -F'\\t' 'BEGIN{OFS=FS}{print \$1, \$2-1, \$2, \$1 "_" \$2 "_" \$3 "_" \$4 ":" \$6, ".", \$6 }' \\
        | bgzip -@ $task.cpus > $bed

    tabix --threads $task.cpus -0 -p bed $bed

    bgzip -@ $task.cpus -cd $bed | cut -f1 | uniq > $chroms
    """
}