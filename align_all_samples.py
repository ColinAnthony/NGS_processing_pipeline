#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import subprocess
from subprocess import DEVNULL
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


def fasta_to_dct_rev(file_name):
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
        dct[str(v).replace("~", "_")] = new_key

    return dct


def align_dna(in_file, out_file):
    """
    does a mafft multiple sequence alignment on a fasta file
    :param in_file: (str) path and filename of the file to align
    :param out_file: (str) path and nae of the aligned outfile
    :return: (dict) dictionary of aligned sequences for all conserved regions
    """
    cmd = 'mafft {0} > {1}'.format(in_file, out_file)
    print("Aligning the sequences, please wait")
    subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    # todo error message add mafft to path


def main(infile, outpath, name):
    outpath = os.path.abspath(outpath)
    infile = os.path.abspath(infile)
    outfile = os.path.join(outpath, name)
    # run mafft alignment on sequences
    align_dna(infile, outfile)

    print("Aligning completed")


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
