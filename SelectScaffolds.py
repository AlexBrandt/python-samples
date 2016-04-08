#!/usr/bin/python

import gzip
import argparse


def main():
    parser = argparse.ArgumentParser(
            description='Select scaffolds from .fa or .fa.gz.')
    parser.add_argument("query", help="comma separated scaffold names",
                        type=str)
    parser.add_argument("scaffolds_fasta", help="the reference .fa or .fa.gz",
                        type=str)
    parser.add_argument("-g", "--gzip", help="the file is gzipped",
                        action='store_true')
    args = parser.parse_args()
    targets = args.query.split(',')

    if args.gzip:
        file = gzip.open(args.scaffolds_fasta)
    else:
        file = open(args.scaffolds_fasta)

    active = False
    for line in file:
        if ">" in line:
            if line.strip()[1:] in targets:
                print line,
                active = True
            else:
                active = False
        else:
            if active:
                print line,

if __name__ == "__main__":
    main()
