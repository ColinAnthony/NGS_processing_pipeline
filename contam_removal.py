#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import collections
from Bio import SeqIO

__author__ = 'Colin Anthony'


def find_contam(inseq):
    '''
    :param inseq: (str) the sequence to blast
    :return: (bool) True or False depending on whether sequence is hiv or not
    '''

    # todo blast the seq, if not hiv
    blast_result = ''

    if blast_result is 'HIV':
        contam_result = False
    else:
        contam_result = True
    return contam_result


def main(read1, read2, outpath):

    cln_read1 = read1.replace("R1.fastq", "cln_R1.fastq")
    cln_read2 = read1.replace("R2.fastq", "cln_R2.fastq")
    read1_out = os.path.join(outpath, cln_read1)
    read2_out = os.path.join(outpath, cln_read2)

    # initialize dict with key = seq name, value = True
    contam_names_d = collections.defaultdict()
    read1_kept = 0
    read2_kept = 0
    read1_removed = 0
    read2_removed = 0

    # clear out any existing outfiles for read1
    with open(read1_out, 'w') as handle:
        handle.write("")

    # check each sequence to see if it is a contaminating sequence
    for seq_record_R1 in SeqIO.parse(open(read1), "fastq"):
        name = seq_record_R1.name
        seq = seq_record_R1.seq
        is_contam = find_contam(seq)

        # if is contam save name to remove from read2
        if is_contam:
            contam_names_d[name] = is_contam
            read1_removed += 1

        # if is not contam, append to list to write to file
        else:
            read1_kept += 1
            # todo use biowriter?
            with open(read1_out, 'a') as handle:
                handle.write(seq_record_R1)

    # clear out any existing outfiles for read2
    with open(read2_out, 'w') as handle2:
        handle2.write("")

    # process read2 sequences, removing the paired reads that were removed in read1
    for seq_record_R2 in SeqIO.parse(open(read2), "fastq"):
        name = seq_record_R2.name
        name_parts = name.split(' ')
        name_parts[-1] = name_parts[-1].replace("2:", "1:")
        lookup_name = " ".join(name_parts)
        if lookup_name in contam_names_d.keys():
            read2_removed += 1
        else:
            read2_kept += 1
            # todo use biowriter?
            with open(read2_out, 'a') as handle2:
                handle2.write(seq_record_R2)

    # check for errors. If read1 and read2 are not of equal length, downstream processes will fail
    if len(read1_kept) != len(read2_kept):
        print("{0} sequences kept from read1 file\n{1} sequences kept from read2\n".format(read1_kept, read2_kept))
        print("Read1 and Read2 must be the same length, check the raw files (/0raw/) "
              "some sequences may have been removed before running the pipeline")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-r1', '--read1', default=argparse.SUPPRESS, type=str,
                        help='The read1 (R1) fastq file', required=True)
    parser.add_argument('-r2', '--read2', default=argparse.SUPPRESS, type=str,
                        help='The read2 (R2) fastq file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be copied', required=True)

    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2
    outpath = args.outpath

    main(read1, read2, outpath)
