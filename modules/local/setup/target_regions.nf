
process TARGET_REGIONS {
    cpus   { 4                        }
    memory { 4 * task.attempt + ' GB' }
    time   { 1 * task.attempt + ' h'  }
    label 'bedops'

    input:
    path include_bed
    path exclude_bed

    output:
    tuple path(output), path("${output}.tbi"), path('chrom_list.txt')

    script:
    output = "targets.bed.gz"
    range = params.include_pad ? "--range -$params.include_pad:$params.include_pad" : ''
    
    exclude_bed.size() == 0 ?
        """
        bgzip -cd --threads $task.cpus $include_bed \\
            | sort-bed - \\
            | bedops $range --merge - \\
            | bgzip --threads $task.cpus > $output
        
        tabix -0 -p bed $output

        bgzip --threads $task.cpus -cd $output \\
            | cut -f1 \\
            | uniq > chrom_list.txt
        """ :
        """
        bgzip -cd --threads $task.cpus $include_bed \\
            | sort-bed - \\
            | bedops $range --merge - \\
            > include_bed &
        
        bgzip -cd --threads $task.cpus $exclude_bed \\
            | sort-bed - \\
            | bedops --merge - \\
            > exclude_bed
           
        wait

        bedops --difference include_bed exclude_bed \\
            | bgzip --threads $task.cpus > $output

        rm include_bed exclude_bed
        
        tabix -0 -p bed $output

        bgzip --threads $task.cpus -cd $output \\
            | cut -f1 \\
            | uniq > chrom_list.txt
        """
}