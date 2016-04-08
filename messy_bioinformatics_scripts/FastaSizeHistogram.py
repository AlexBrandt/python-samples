# Plots fasta file contig size

from Bio import SeqIO
import argparse
import numpy as np
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
parser = argparse.ArgumentParser()
parser.add_argument("f", type=str,
        help="The fasta file to be analyzed.")

parser.add_argument("o", type=str,
        help="The output pdf.")

parser.add_argument("c", nargs="?", default=100, type=int,
        help="Bin counts (optional)")


parser.add_argument("lower", nargs="?", default=0, type=int,
        help="x lower bound")

parser.add_argument("upper", nargs="?", default=0, type=int,
        help="x upper bound")

args = parser.parse_args()
contig_length = []

fasta_handle = open(args.f, 'rU')
for record in SeqIO.parse(fasta_handle, "fasta") :
        contig_length.append(len(record.seq))
plt.hist(contig_length, args.c, linewidth=0)
# plt.xscale('log')
plt.xlabel('Contig length')
plt.ylabel('Counts')
plt.title('Contig Histogram (bins = %i)' % (args.c))
if args.lower or args.upper:
    plt.xlim([args.lower,args.upper])
plt.savefig(args.o, bbox_inches='tight',alpha=.8)

print "Min: ",
print np.amin(contig_length)
print "Max: ",
print np.amax(contig_length)
print "Mean: ",
print np.mean(contig_length)
print "Std. dev.: ",
print np.std(contig_length)
