#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import collections
from Bio import SeqIO


__author__ = 'Colin Anthony'


def main(infile, outpath, name):
    print(infile)
    # Todo: add code to de-multiplex a raw fastq file based on gene regions
    # todo: look at existing options seqmagic/fastqc/vsearch?


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input fasta file, with all the time points in one file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the prefix name of your outfile', required=True)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name

    main(infile, outpath, name)
