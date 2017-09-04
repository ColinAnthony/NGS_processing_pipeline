#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import collections
import argparse
from Bio import SeqIO
#

__author__ = 'colin'


def hapl_collapse(fastafile, seq_name, outfile):
    '''
    Haplotype a fasta file by collapsing identical sequences into a single sequence with a frequency count.
    :param fastafile: input fasta file
    :param seq_name: (str) prefix for sequence name
    :param outfile: the path to where the outfile is created
    :return: dictionary of key = sequence, value  = count of that sequence
    '''

    # import seqs to dict and get count of each unique sequence
    total = 0
    dct = collections.defaultdict(int)
    for seq_record in SeqIO.parse(open(fastafile), "fasta"):
        dct[str(seq_record.seq).replace("~", "-").upper()] += 1
        total += 1

    # list of sequences, sorted by count, max to min
    sorted_sequences_by_count = sorted(dct, key=dct.get, reverse=True)

    # clear any existing file previously created
    with open(outfile, "w") as handle:
        handle.write("")

    # initialize a counter to add an incremental number to each sequence
    for counter, sequence in enumerate(sorted_sequences_by_count):
        frq = round(dct[sequence] / float(total), 3)
        num = str(counter).zfill(3)
        n = seq_name + "_" + str(num) + "_" + str(frq)

        # write the haplotyped sequence to file
        with open(outfile, "a") as handle:
            handle.write('>{0}\n{1}\n'.format(n, sequence))


def main(infile, outpath):

    print("Collapsing to haplotypes on file: \n\t {}".format(infile))

    # set the outfile name
    name = os.path.split(infile)[-1]
    name = name.replace("_sep.fasta", "_hap.fasta")
    out = os.path.join(outpath, name)

    # get prefix for sequence names from infile name
    seq_name = "_".join(infile.split("_")[:4])

    # run the haplotyping function
    hapl_collapse(infile, seq_name, out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='collapses sequences by unique, with freq added to the name',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', type=str, required=True,
                        help='The input fata file')
    parser.add_argument('-o', '--outpath', required=True, type=str,
                        help='path for output.')

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath

    main(infile, outpath)
