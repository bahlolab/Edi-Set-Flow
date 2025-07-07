#!/usr/bin/awk -f
#
# gtf2bed.awk — Convert a GTF into a BED file,
#
BEGIN {
    FS = OFS = "\t"
}
# skip header/comment lines
/^#/ { next }

{
    # name the relevant GTF columns
    chrom     = $1
    gtf_start = $4
    gtf_end   = $5
    strand    = $7
    attrs     = $9

    # pull out gene_id
    name = "."
    n = split(attrs, a, ";")
    for (i = 1; i <= n; i++) {
        sub(/^[ \t]+/, "", a[i])
        if (a[i] ~ /^gene_id/) {
            gsub(/gene_id[ \t]*"|"/, "", a[i])
            name = a[i]
            break
        }
    }

    # convert to BED5
    bed_start = gtf_start - 1    # zero-based
    bed_end   = gtf_end          # end stays 1-based

    print chrom, bed_start, bed_end, name, ".", strand
}
