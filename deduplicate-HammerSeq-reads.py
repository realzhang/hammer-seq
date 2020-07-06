#!/usr/bin/env python3
# input is STDIN, that can be piped from the output of samtools view of mapped SAM file.
# zhangzhuqiang@ibp.ac.cn 
# last update: 2020.07.05

import sys
import re

total = 0
dup = 0
set_dup = set()
resultFH = open("dedup.report.txt", "w")
for line in sys.stdin:
    r1 = line
    r2 = sys.stdin.readline() # read in one line, including the endline;
    total += 1
    flag, seq1 = [ re.split(r"\s+", r1)[i] for i in [1,9] ]
    # print("Flag is %s"%flag)
    # print("Seq1 is %s"%seq1)
    seq2 = re.split(r"\s+", r2)[9]
    comb_seq = ( seq1+seq2 if flag == '0' else seq2+seq1 ) # swap the pairs if in the opposite direction.

    if( comb_seq in set_dup ):
        dup += 1
        # print("Found dup!")
        continue
    else:
        set_dup.add( comb_seq )
        # print("New!")
        print(r1, end='')
        print(r2, end='')

resultFH.write("Total reads:\t%d\nDuplicates:\t%d (%f%%)\n" % (total, dup, dup/total*100)  )
resultFH.close()
