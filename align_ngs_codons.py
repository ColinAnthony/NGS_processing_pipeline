#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import collections
from itertools import groupby
import regex
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA
# from Bio import pairwise2


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
        dct[new_key] = v.upper()

    return dct


def fasta_to_dct_rev(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(list)
    my_gen = py3_fasta_iter(file_name)
    for k, v in my_gen:
        new_key = k.replace(" ", "_")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[str(v).replace("~", "_").upper()].append(new_key)

    return dct


def prelim_align(sequence, reference, name):
    """
    Pairwise align sequence to reference, to find reading frame and frame-shift in/dels
    :param sequence: (str) a query DNA sequence
    :param reference: (str) a reference DNA sequence (must start in reading frame 1)
    :param name: the name of the query sequence
    :return: (str) aligned query sequence, (str) aligned ref sequence, (int) reading frame for query sequence
    """
    # # biopython pairwise (scoring: match,  miss-match, gap open, gap extend)
    # alignment = pairwise2.align.localms(sequence, reference, 3, -1, -8, -2)
    # start_end_positions = [(0, alignment[4]), (alignment[3], alignment[4])]

    # reference = consensus of subtype C
    alignment, score, start_end_positions = local_pairwise_align_ssw(DNA(sequence), DNA(reference), gap_open_penalty=8,
                                                                     gap_extend_penalty=2, match_score=4,
                                                                     mismatch_score=-1.5, substitution_matrix=None)

    # get the query and ref aligned seqs
    seq_align = str(alignment[0])
    ref_align = str(alignment[1])
    # print(">{0}\n{1}".format(name, seq_align))
    # print(">ref\n{}".format(ref_align))

    # get start position for reference and query
    ref_start = start_end_positions[1][0]
    seq_start = start_end_positions[0][0]

    # calculate reading frame (reference must start in frame 1)
    frame = (ref_start - seq_start) % 3

    # report align start pos if alignment starts after query seq start
    if start_end_positions[0][0] != 0:
        print("  seq    ,    ref")
        print("start/end, start/end")
        print(start_end_positions)
        sys.exit("Preliminary alignment issue\nQuery sequence was truncated during alignment\nExiting")

    return seq_align, ref_align, frame


def gap_padding(seq_align, ref_align, frame, sequence, name):
    """
    pads sequence with gaps for indels and to set reading frame to frame 1
    :param seq_align: (str) an aligned query sequence
    :param ref_align: (str) the corresponding aligned reference sequence
    :param frame: (int) the reading frame for the query sequence
    :param sequence: (str) the un-aligned query sequence
    :param name: (str) the query sequence name
    :return: (str) a gap-padded query sequence to correct for indels and set reading frame to frame 1
    """
    # convert query seq to list to allow mutability
    new_seq = list(sequence)

    # find frame shirt deletions (gap in query)
    for i in range(frame, len(seq_align)):
        if seq_align[i] == "-" and i > 2:
            if seq_align[i + 1] != "-" and seq_align[i - 1] != "-":
                new_seq.insert(i, "-")

    # find frame shirt insertions (gap in reference)
    for i in range(frame, len(ref_align)):
        if ref_align[i] == "-":
            if ref_align[i + 1] != "-" and ref_align[i - 1] != "-":
                gap_frame = i % 3
                if gap_frame == 0:
                    new_seq.insert(i + 1, "--")
                if gap_frame == 1:
                    new_seq.insert(i + 3, "--")
                if gap_frame == 2:
                    new_seq.insert(i + 2, "--")

    # pad query to be in reading frame 1
    if frame != 3:
        lead_gap = "-" * frame
        new_seq.insert(0, lead_gap)

    # convert to string and return
    new_seq = "".join(new_seq)
    # print(">{0}\n{1}".format(name, new_seq))
    return "".join(new_seq)


def translate_dna(sequence):
    """
    :param sequence: (str) a DNA sequence string
    :return: (str) a protein string from the forward reading frame 1
    """
    codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
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
    prot = []

    for n in range(0, len(seq), 3):
        if seq[n:n + 3] in codontable:
            residue = codontable[seq[n:n + 3]]
        else:
            residue = "X"

        prot.append(residue)

    return "".join(prot)


def split_regions(sequence):
    print("")


def main(infile, ref, outpath, name):

    # get absolute paths
    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)
    name = name + "_aligned.fasta"
    outfile = os.path.join(outpath, name)
    # print(outfile)

    # read in fasta file
    in_seqs_d = fasta_to_dct_rev(infile)
    reference = list(fasta_to_dct(ref).values())[0]

    # generate seq_code to seq name list lookup dictionary
    first_look_up_d = collections.defaultdict(list)
    first_seq_code_d = collections.defaultdict(str)
    for i, (seq, names_list) in enumerate(in_seqs_d.items()):
        unique_id = str(i).zfill(4)
        first_look_up_d[unique_id] = names_list
        first_seq_code_d[seq] = unique_id

    # translate sequences
    for seq, code in first_seq_code_d.items():
        s_name = first_look_up_d[code][0]
        # print(s_name)
        seq_align, ref_align, frame = prelim_align(seq, reference, s_name)
        padded_sequence = gap_padding(seq_align, ref_align, frame, seq, s_name)
        prot_seq = translate_dna(padded_sequence)

    print("Aligning completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input file', required=True)
    parser.add_argument('-r', '--ref', default=argparse.SUPPRESS, type=str,
                        help='The reference sequence file. Must be in reading frame 1', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path for the output file', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='The name for the output file', required=True)

    args = parser.parse_args()
    infile = args.infile
    ref = args.ref
    outpath = args.outpath
    name = args.name

    main(infile, ref, outpath, name)
