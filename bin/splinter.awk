#!/usr/bin/awk -f
# splinter.awk: split a BED stream into N segments,
# collapsing each chromosome to one [min…max] per segment,
# and outputting each segment as a comma-separated list of 1-based intervals to stdout.

# function must be defined at top level in awk
function append_curr() {
    region = curr_chr ":" (curr_lo+1) "-" curr_hi
    if (line) {
        line = line "," region
    } else {
        line = region
    }
}

BEGIN {
    if (N == "" || N < 1) {
        print "Error: must call with -v N=<number_of_segments> (N>=1)" > "/dev/stderr"
        exit 1
    }
}

/^#/ { next }

# 1st pass: read & store
{
    n++
    chr[n]   = $1
    st[n]    = $2 + 0
    en[n]    = $3 + 0
    len[n]   = en[n] - st[n]
    total   += len[n]
}

END {
    if (n < 1) exit

    width = total / N
    f     = 1
    sum   = 0
    line  = ""

    for (i = 1; i <= n; i++) {
        d = len[i]
        if (i == 1) {
            curr_chr = chr[i]
            curr_lo  = st[i]
            curr_hi  = en[i]
            sum      = d
            continue
        }

        # chromosome change
        if (chr[i] != curr_chr) {
            append_curr()
            curr_chr = chr[i]
            curr_lo  = st[i]
            curr_hi  = en[i]
            sum      += d
            continue
        }

        # width exceed: finish segment
        if (sum + d > width && f < N) {
            append_curr()
            print line
            f = f + 1
            line = ""
            curr_chr = chr[i]
            curr_lo  = st[i]
            curr_hi  = en[i]
            sum      = d
            continue
        }

        # extend current collapse
        if (st[i] < curr_lo) curr_lo = st[i]
        if (en[i] > curr_hi) curr_hi = en[i]
        sum += d
    }

    # final flush
    append_curr()
    print line
}
