process SCATTER {
    cpus    2
    memory '4 GB'
    time   '6 h'
    label  'bcftools'
    publishDir "${params.outdir}/edisites", mode: 'copy'
    /*
        - Convert to TSV with bcftools +split-vep
        - Scatter for parallel processing
    */

    input:
    tuple path(vcf), path(idx)

    output:
    path("edisites_shard_*.tsv.gz")

    script:
    info = [
        'ID',
        'ADAR',
        'PASS_FREQ',
        'REDI_ACC',
        'REDI_REP_TYPE',
        'REDI_REP_ID',
        'REP_ID',
        // VEP fields
        'Consequence',
        'IMPACT',
        'SYMBOL',
        'Gene',
        'Feature_type',
        'Feature',
        'BIOTYPE',
        'EXON',
        'INTRON',
        'cDNA_position',
        'CDS_position',
        'Protein_position',
        'Amino_acids',
        'Codons',
        'STRAND'
    ]
    info_hdr = info.join('\\t')
    info_qry = info.collect {"%$it"}.join('\\t')

    """
    SAM_NALT=\$( bcftools query -l $vcf | sed 's:^:NALT_:' | paste -sd '\\t' - )
    SAM_NREF=\$( bcftools query -l $vcf | sed 's:^:NREF_:' | paste -sd '\\t' - )
    HEADER="${info_hdr}\\t\$SAM_NALT\\t\$SAM_NREF"

    bcftools +split-vep $vcf -d -A tab -f '${info_qry}[\\t%NALT][\\t%NREF]\\n' \\
        | awk -v header="\$HEADER" '{
            i = NR % $params.analysis_threads
            if (!(i in seen)) {
                filename = sprintf("edisites_shard_%03d.tsv", i)
                print header > filename
                seen[i] = filename
            }
            print >> seen[i]
        }'
    
    bgzip --threads $task.cpus edisites_shard_*.tsv
    """
}