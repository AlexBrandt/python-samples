import pysam
import sys
filename = sys.argv[1]
print reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(filename) ])
