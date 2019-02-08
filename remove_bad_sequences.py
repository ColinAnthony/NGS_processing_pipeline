#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import argparse
from itertools import groupby
import os
import collections


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
        dct[new_key] = str(v).replace("~", "_").upper()

    return dct


def translate_dna(sequence):
    """
    :param sequence: a DNA string
    :return: a protein string from the forward reading frame with the least number of stop codons, if the lowest number is equal between two reading frames "*" is returned
    """
    codontable = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        '---': '-',
    }

    seq = sequence.upper()
    protein_sequence = []

    for n in range(0, len(seq), 3):
        if seq[n:n + 3] in codontable:
            protein_sequence.append(codontable[seq[n:n + 3]])

        else:
            protein_sequence.append("X")

    p1_stop = protein_sequence[:-1].count("*")
    # p1_err = protein_sequence.count("X")

    if p1_stop == 0:  # and p1_err == 0:
        return "".join(protein_sequence)
    else:
        return None


def degen_remove(d):
    """
    :param d: (dict) dictionary of sequence names and DNA sequences)
    :return: (dict) (dict) dictionary with sequences with degenerate bases removed
    """
    badseq = ['M', 'R', 'W', 'S', 'Y', 'K', 'B', 'D', 'H', 'V', 'N']
    degen = 0
    n_d = collections.defaultdict(str)
    bad_d = collections.defaultdict(str)
    for name, seq in d.items():
        found = False

        # only replace consecutive "N's" with "", over the first and last 4 bp of the read
        if "N" in seq[-4:]:
            replace_n_end = seq[-4:].replace("N", "")
            seq = seq[:-4] + replace_n_end

        if "N" in seq[:4]:
            replace_n_start = seq[:4].replace("N", "")
            seq = replace_n_start + seq[4:]

        for base in badseq:
            if base in seq and not found:
                bad_d[name] = seq
                degen += 1
                found = True
        if not found:
            n_d[name] = seq

    return n_d, bad_d, degen


def stops_remove(d, frame):
    """
    :param d: (dict) dictionary of sequence names and DNA sequences)
    :return: (dict) dictionary with sequences with stop codons removed
    """
    stops = 0
    good_d = collections.defaultdict(str)
    bad_d = collections.defaultdict(str)
    for name, seq in d.items():
        pseq = translate_dna(seq[frame:])
        if pseq is None:
            bad_d[name] = seq
            stops += 1
        else:
            good_d[name] = seq

    return good_d, bad_d, stops


def length_check(d):
    """
    :param d: (dict) dictionary of sequence names and DNA sequences)
    :return: (dict) dictionary with short sequences removed
    """
    short = 0
    n_d = collections.defaultdict(str)
    bad_d = collections.defaultdict(str)

    for name, seq in d.items():
        seq_check = seq.replace("-", "")
        if len(seq_check) < length:
            bad_d[name] = seq
            short += 1
        else:
            n_d[name] = seq

    return n_d, bad_d, short


def main(infile, outp, frame, stops, length, logfile):

    # set outfile names
    n = os.path.split(infile)[1]
    out_gd = n.replace(".fasta", '_clean.fasta')
    # out_bd = n.replace(".fasta", '_dirty.fa')
    outfile_good = os.path.join(outp, out_gd)
    # outfile_bad = os.path.join(outp, out_bd)
    d = fasta_to_dct(infile)
    inseq_no = len(d)

    # remove sequences with degenerate bases
    print('removing degenerates')
    cln1_d, bad_d1, degen_no = degen_remove(d)

    # remove sequences with stop codons
    stops_no = 0
    if stops:
        print('removing stops')
        # account for python zero indexing
        frame = int(frame) - 1
        cln2_d, bad_d2, stops_no = stops_remove(cln1_d, frame)
    else:
        cln2_d = cln1_d
        bad_d2 = {}

    # remove short sequences
    if length is None:
        print('removing short sequences')
        cln3_d = cln2_d
        bad_d3 = {}
        short_no = 0
    else:
        cln3_d, bad_d3, short_no = length_check(cln2_d)

    # get totals for kept and removed seqs
    kept = len(cln3_d)
    removed = inseq_no - kept

    # write the cleaned sequences to file
    with open(outfile_good, "w") as handle:
        for seq_name, seq in cln3_d.items():
            handle.write(">{0}\n{1}\n".format(seq_name, seq))

    ## merge the 'bad sequence' dictionaries
    # bad_part = dict(bad_d1, **bad_d2)
    # all_bad = dict(bad_part, **bad_d3)

    ## write the bad sequences to file
    # with open(outfile_bad, "w") as handle:
    #     for seq_name1, seq in all_bad.items():
    #         handle.write(">{0}\n{1}\n".format(seq_name1, seq))

    # calculate stats on removed sequences for each operation
    percent_degen = round((degen_no / inseq_no) * 100, 2)
    if stops:
        percent_stop = round((stops_no / inseq_no) * 100, 2)
    else:
        percent_stop = None
    if length is not None:
            percent_short = round((short_no / inseq_no) * 100, 2)
    else:
        percent_short = None
    percent_kept = round((kept / inseq_no) * 100, 2)
    percent_tot_rem = round((removed / inseq_no) * 100, 2)

    # write stats to log file
    with open(logfile, 'a') as handle:
        handle.write("\n{0}\nFile:                                    = {1}          \n".format(("-"*40), infile))
        handle.write("Number of input sequences                = {0}          \n".format(inseq_no))
        handle.write("Number of sequences with non ACGT bases  = {0} - ({1} %)\n".format(degen_no, percent_degen))
        if stops:
            handle.write(" of sequences with stop codons           = {0} - ({1} %)\n".format(stops_no, percent_stop))
        if length is not None:
            handle.write("Number of short sequences                = {0} - ({1} %)\n".format(short_no, percent_short))
        handle.write("Total number of sequences removed        = {0} - ({1} %)\n".format(removed, percent_tot_rem))
        handle.write("Total number of sequences kept           = {0} - ({1} %)\n".format(kept, percent_kept))
        handle.write("\n")
        if stops:
            handle.write("Note, any frame shift errors not introducing stop codons will remain\n")
            handle.write("\n")

    print("for file:", n)
    print("Sequences kept = {0} - ({1} %)".format(kept, percent_kept))
    print("-"*40 + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='removes sequences with stop codons and '
                                                 'degenerate bases (will remove most indel errors',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--infile', type=str,
                        help='The input .fastq file', required=True)
    parser.add_argument('-o', '--outpath', type=str,
                        help='The path to where the output files will be created', required=True)
    parser.add_argument('-f', '--frame', default=1,
                        help='The reading frame (1, 2 or 3)', required=False)
    parser.add_argument('-s', '--stops', default=False, action='store_true',
                        help='Remove sequences with stop codons?)', required=False)
    parser.add_argument('-l', '--length', type=int, default=250,
                        help='The minimum read length)', required=False)
    parser.add_argument('-lf', '--logfile', default=argparse.SUPPRESS, type=str,
                        help='The path and name of the log file', required=True)


    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    frame = args.frame
    stops = args.stops
    length = args.length
    logfile = args.logfile


    main(infile, outpath, frame, stops, length, logfile)
