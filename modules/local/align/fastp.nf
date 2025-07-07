
process FASTP {
    cpus     6
    memory { 4 * task.attempt + ' GB' }
    time   { 1 * task.attempt + ' h' }
    tag    "$sample"
    label  'fastp'

    input:
    tuple val(sample), path(fastqs)

    output:
    tuple val(sample), path("${sample}.fastp.R{1,2}.fq.gz"), emit: fastq
    path "${sample}.fastp.json",                             emit: qc

    script:
    // Build the --in/--out flags depending on SE vs PE
    def fastp_io = fastqs.size() == 2 ?
        "--in1 ${fastqs[0]} --in2 ${fastqs[1]} " +
        "--out1 ${sample}.fastp.R1.fq.gz --out2 ${sample}.fastp.R2.fq.gz" :
        "--in1 ${fastqs[0]} --out1 ${sample}.fastp.R1.fq.gz"

    // Toggle minimal vs strict filtering
    def filter_opts = params.fastp_mode == 'STRICT' ?
        // Strict: quality & poly‐tail trimming + adapter removal
        """\
        --detect_adapter_for_pe \\
        --cut_front \\
        --cut_tail \\
        --qualified_quality_phred 15 \\
        --unqualified_percent_limit 40 \\
        --cut_mean_quality 20 \\
        --n_base_limit 5 \\
        --trim_poly_x \\
        --poly_x_min_len 7 \\
        --n_base_limit 5"""
        .stripIndent().replaceAll('\n', '\n        ') :
        // Minimal: adapter removal pairing only
        "        --disable_quality_filtering"

    """\
    fastp $fastp_io \\
        --thread $task.cpus \\
        --compression 6 \\
        --length_required 30 \\
        $filter_opts \\
        --dont_eval_duplication \\
        --json ${sample}.fastp.json

    sed -i 's;${fastqs.size()==2 ? fastqs[0] : fastqs};${sample};g' ${sample}.fastp.json
    """
}