#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import argparse
import collections
from Bio import SeqIO
#

__author__ = 'Colin Anthony'


def gethxb2(dict):
    '''
    :param dict: a dictionary of your aligned input sequences. Must contain HXB2, with HXB2 in the header
    :return: the HXB2 sequence as a string
    '''
    found = False
    hxb2_seq = None
    hxb2_key = None
    for k in dict.keys():
        if "HXB2" in k.upper():
            found = True
            hxb2_key = k
            hxb2_seq = dict[k]
            print("Found hxb2 ref. seq. Its full name is: ", hxb2_key)
            break
    if not found:
        print("We could not find a sequence with 'HXB2' in its name. "
              "Please make sure you have an HXB2 ref. seq. included")

    return str(hxb2_key), str(hxb2_seq)


def main(infile, outpath, field):

    alignment = SeqIO.parse(infile, "fasta")
    d = collections.defaultdict(list)
    hxb2_key, hxb2_seq = gethxb2(d)
    del d[hxb2_key]

    # account for python zero indexing
    field -= 1

    # get unique field(s) for splitting by time point
    for sequence_obj in alignment:
        unique_field = sequence_obj.description.split("_")[0:field]
        new_name = "_".join(unique_field) + "_sep.fasta"
        out_file_name = os.path.join(outpath, new_name)
        d[out_file_name].append(sequence_obj)

    # write the grouped sequences to their outfiles
    for out_file, seq_objs in d.items():
        for seq_obj in seq_objs:
            with open(out_file, 'a') as handle:
                handle.write(">{0}\n{1}\n".format(seq_obj.name, str(seq_obj.seq)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots loop stats from csv file (produced by loop_stats.py)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', type=str, required=True,
                        help='The input fasta file')
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-f', '--field', type=int, default=4, required=False,
                        help="The field that differentiates your samples/time points (use the last field if multiple."
                             "(ie: 4 for 'CAP177_2000_004wpi_V3C4_GGGACTCTAGTG_28, or 2 for SVB008_SP_GGTAGTCTAGTG_231")

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    field = args.field

    main(infile, outpath, field)
