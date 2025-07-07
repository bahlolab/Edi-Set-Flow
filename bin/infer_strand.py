#!/usr/bin/env python3
import sys
import re

result = open(sys.argv[1]).read()
thresh = 0.75
pattern = re.compile(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?$')
frac = [float(match.group(0)) for line in result.splitlines() if (match := pattern.search(line))]
frac_fr = frac[1] / (frac[1] + frac[2])
frac_rf = frac[2] / (frac[1] + frac[2])

pair = "PAIRED" if "PairEnd" in result else "SINGLE"

strand = 'UNSTRANDED'
strand = 'FR-SECONDSTRAND' if frac_fr > thresh else strand
strand = 'RF-FIRSTSTRAND'  if frac_rf > thresh else strand

print(pair + '-' + strand)
