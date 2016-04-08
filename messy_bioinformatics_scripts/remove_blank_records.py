#!/usr/common/usg/languages/python/2.7.4/bin/python

import sys, gzip

f1 = gzip.open(sys.argv[1], 'rb')
f2 = gzip.open(sys.argv[2], 'rb')

w1 = gzip.open(sys.argv[1].strip('.fq.gz').strip('.fastq.gz')+'.cln.fastq.gz', 'wb')
w2 = gzip.open(sys.argv[2].strip('.fq.gz').strip('.fastq.gz')+'.cln.fastq.gz', 'wb')

lines1 = [f1.readline() for x in range(4)]
lines2 = [f2.readline() for x in range(4)]

while lines1[0] != "" and lines2[0] != "":
    if lines1[1] != "\n" and lines2[1] != "\n":
        for i in range(4):
            w1.write(lines1[i])
            w2.write(lines2[i])

    lines1 = [f1.readline() for x in range(4)]
    lines2 = [f2.readline() for x in range(4)]

f1.close(); f2.close();
w1.close(); w2.close();

if lines1[0] != lines2[0]:
    print "The fastq files were not of equal length.  Please check files and try again."
    print lines1
    print lines2
    exit(1)

exit(0)
