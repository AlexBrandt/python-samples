from Bio import SeqIO
import sys

input_handle = open(sys.argv[1], "r")
output_handle = open(sys.argv[1].strip('.fa') + '.selected.fa', "w")
input_seq_iterator = SeqIO.parse(input_handle, "fasta")

short_sequences = [record for record in input_seq_iterator \
    if len(record.seq) >= 1000]

SeqIO.write(short_sequences, output_handle, "fasta")
output_handle.close()
