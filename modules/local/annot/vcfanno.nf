
process VCFANNO {
    cpus   { 2                        }
    memory { 2 * task.attempt + ' GB' }
    time   { 2 * task.attempt + ' h'  }
    label 'vcfanno'
    /*
        - Annotate VCF with vcfanno
        - Adding columns "REDI_ACC", "REDI_REP_TYPE", "REDI_REP_ID" from REDIPORTAL
    */
    
    input:
    tuple path(vcf_in),   path(index_in)
    tuple path(redi_bgz), path(redi_idx)
    tuple path(rep_bed),  path(rep_idx)

    output:
    tuple path(output), path("${output}.csi")

    script:
    rep_block = 
        """\
        [[annotation]]
        file    = "$rep_bed"
        columns = [4, 5, 6]
        names   = ["REP_CLASS", "REP_FAM", "REP_ID"]
        ops     = ["self", "self", "self"]
        """.stripIndent().replaceAll('\n', '\n    ')
    
    redi_block = 
        """\
        [[annotation]]
        file    = "$redi_bgz"
        columns = [5, 8, 10]
        names   = ["REDI_ACC", "REDI_REP_TYPE", "REDI_REP_ID"]
        ops     = ["self", "self", "self"]
        """.stripIndent().replaceAll('\n', '\n    ') 

    output  = vcf_in.name.replaceFirst(/.vcf\.gz$/, ".vcfanno.vcf.gz")

    """
    cat << 'VCFANNO_CONF_EOF' > vcfanno.conf
    $rep_block
    $redi_block
    VCFANNO_CONF_EOF

    vcfanno -p $task.cpus vcfanno.conf $vcf_in \\
        | bgzip --threads $task.cpus > $output
    
    tabix -p vcf -C $output
    """
}
