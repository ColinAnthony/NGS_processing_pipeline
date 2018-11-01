#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import argparse
import collections
from itertools import groupby


__author__ = 'Colin Anthony'


def py3_fasta_iter(fasta_name):
    """
    modified from Brent Pedersen: https://www.biostars.org/p/710/#1412
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(str(fasta_name), 'r')
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header_str = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header_str, seq)


def fasta_to_dct(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = py3_fasta_iter(file_name)
    for k, v in my_gen:
        new_key = k.replace(" ", "_")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[new_key] = str(v).replace("~", "_")

    return dct


def main(infile, outpath, field):

    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)

    alignment = fasta_to_dct(infile)
    d = collections.defaultdict(list)

    # get unique field(s) for splitting by time point
    for name, seq in alignment.items():
        if "HXB2" in name:
            print("found HXB2, excluding from output")
        else:
            unique_field = name.split("_")[0:field]
            new_name = "_".join(unique_field) + "_sep.fasta"
            out_file_name = os.path.join(outpath, new_name)
            d[out_file_name].append((name, str(seq)))

    # write the grouped sequences to their outfiles
    for out_file, seq_objs in d.items():
        for seq_obj in seq_objs:
            with open(out_file, 'a') as handle:
                handle.write(">{0}\n{1}\n".format(seq_obj[0], seq_obj[1]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots loop stats from csv file (produced by loop_stats.py)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--infile', type=str, required=True,
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
