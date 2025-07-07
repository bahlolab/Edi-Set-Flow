#!/usr/bin/awk -f
#
# to_vcf.awk — JACUSA2 -> VCF, splitting multiallelics but preserving depth
#               even when no ALT is observed

# infer_enzyme(ref, alt, orig_strand, name, enz_ref, enz_alt)
#   returns name if the REF/ALT on the given strand matches the enzyme,
#   and (if orig_strand==".") will update the global `strand` to "+" or "-".
function infer_enzyme(r, a, p_str, name, eref, ealt, res) {
    res = ""
    # plus-strand or ambiguous
    if (p_str=="+" || p_str==".") {
        if (r==eref && a==ealt) {
            if (p_str==".") strand = "+"
            return name
        }
    }
    # minus-strand or ambiguous
    if (p_str=="-" || p_str==".") {
        if (r==comp[eref] && a==comp[ealt]) {
            if (p_str==".") strand = "-"
            return name
        }
    }
    return res
}

BEGIN {
    FS = OFS = "\t"
    # complement map
    comp["A"]="T"; comp["T"]="A"
    comp["C"]="G"; comp["G"]="C"
    # base arrays for picking alt order
    split("ACGT", alt_plus, "")
    split("TGCA", alt_minus, "")
    # map strand+base → index
    for (i=1; i<=4; i++) {
        ref_index["+ " alt_plus[i]]  = i
        ref_index["- " alt_minus[i]] = i
    }

    # VCF headers
    print "##fileformat=VCFv4.2"
    print "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Strand relative to the reference genome\">"
    if (ENZ) print "##INFO=<ID=ADAR,Number=0,Type=Flag,Description=\"Inferred ADAR site\">"
    if (ENZ) print "##INFO=<ID=APOBEC,Number=0,Type=Flag,Description=\"Inferred APOBEC site\">"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype call (0/1)\">"
    print "##FORMAT=<ID=PASS,Number=1,Type=Integer,Description=\"Site passes JACUSA filters\">"
    print "##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency\">"
    print "##FORMAT=<ID=NALT,Number=1,Type=Integer,Description=\"Alternate base count\">"
    print "##FORMAT=<ID=NREF,Number=1,Type=Integer,Description=\"Reference base count\">"
    print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"

    fmt = "GT:PASS:VAF:NALT:NREF"
}

# skip any existing header lines
/^#/ { next }

{
    chrom       = $1
    pos         = $3
    orig_strand = $6
    bases       = $7
    site_flt    = $9
    orig_ref    = $10

    pass         = (site_flt=="*" ? "1" : "0")
    logic_strand = (orig_strand=="." ? "+" : orig_strand)

    # parse allele counts
    split(bases, ac, ",")

    # choose ref & alt arrays based on logic_strand
    if (logic_strand=="-") {
        ref = comp[orig_ref]
        for (i=1; i<=4; i++) alt[i] = alt_minus[i]
    } else {
        ref = orig_ref
        for (i=1; i<=4; i++) alt[i] = alt_plus[i]
    }

    # find index of the REF in the alt array
    ri = ref_index[logic_strand " " ref]

    # emit one VCF line per alt allele
    for (i=1; i<=4; i++) {
        if (alt[i]==ref || (!ALL && ac[i]==0)) continue
        strand      = orig_strand
        alt_base  = alt[i]
        gt        = (ac[i]>0 ? 1 : 0)
        adr       = ac[ri]
        ada       = ac[i]
        dp        = adr + ada
        vaf       = sprintf("%.5f", ada/(dp?dp:1))

        info = "STRAND=" strand
        if (ENZ) {
            # try ADAR, then APOBEC
            enz           = infer_enzyme(ref, alt_base, orig_strand, "ADAR"  ,"A", "G")
            if (!enz) enz = infer_enzyme(ref, alt_base, orig_strand, "APOBEC","C", "T")
            if (enz) info = enz ";STRAND=" strand
        }

        id        = chrom "_" pos "_" ref "_" alt_base ":" strand
        genotype  = gt ":" pass ":" vaf ":" ada ":" adr

        printf("%s\t%s\t%s\t%s\t%s\t.\t.\t%s\t%s\t%s\n",
               chrom, pos, id, ref, alt_base, info, fmt, genotype)
    }
}
