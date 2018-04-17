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
from skbio.alignment import StripedSmithWaterman
from skbio import DNA


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


def prelim_align(sequence, reference):
    """
    Pairwise align sequence to reference, to find reading frame and frame-shift in/dels
    :param sequence: (str) a DNA sequence
    :return: a sequence with gaps to account for frame shift indels and the reading frame
    """
    # reference = consensus of subtype C
    alignment, score, start_end_positions = local_pairwise_align_ssw(DNA(sequence), DNA(reference), gap_open_penalty=8,
                                                                     gap_extend_penalty=2, match_score=2,
                                                                     mismatch_score=-1, substitution_matrix=None)

    # get the query and ref aligned seqs
    seq_align = str(alignment[0])
    ref_align = str(alignment[1])

    # convert query seq to list to allow mutability
    new_seq = list(sequence)

    # get start position for reference and query
    ref_start = start_end_positions[1][0]
    seq_start = start_end_positions[0][0]

    # calculate reading frame (reference must start in frame 1)
    frame = (ref_start - seq_start) % 3

    # report align start pos if alignment starts after query seq start
    if start_end_positions[0][0] != 0:
        print("seq    ,    ref")
        print(start_end_positions)

    # find frame shirt deletions (gap in query)
    for i in range(frame, len(seq_align)):
        if seq_align[i] == "-":
            if seq_align[i+1] != "-" and seq_align[i - 1] != "-":
                new_seq[i] = "-"

    # find frame shirt insertions (gap in reference)
    for i in range(frame, len(ref_align)):
        if ref_align[i] == "-":
            if ref_align[i + 1] != "-" and ref_align[i - 1] != "-":
                new_seq.insert(i+2, "--")

    # pad query to be in reading frame 1
    if frame != 3:
        lead_gap = "-" * frame
        new_seq.insert(0, lead_gap)

    # convert to string and return
    new_seq = "".join(new_seq)

    return "".join(new_seq), frame


def translate_dna(sequence):
    """
    :param sequence: (str) a DNA sequence string
    :return: (str) a protein string from the forward reading frame 1
    """

    # protein_from_dna = sequence.translate()


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

    seq = sequence.upper().replace("-", "")

    prot = []
    prot_fr_1 = []
    prot_fr_2 = []
    prot_fr_3 = []

    for n in range(0, len(seq), 3):
        if seq[n:n + 3] in codontable:
            residue = codontable[seq[n:n + 3]]
        else:
            residue = "X"
        prot_fr_1.append(residue)

    for n in range(1, len(seq), 3):
        if seq[n:n + 3] in codontable:
            residue = codontable[seq[n:n + 3]]
        else:
            residue = "X"
        prot_fr_2.append(residue)

    for n in range(2, len(seq), 3):
        if seq[n:n + 3] in codontable:
            residue = codontable[seq[n:n + 3]]
        else:
            residue = "X"
        prot_fr_3.append(residue)

    # find stop codons
    if "*" not in prot_fr_1:
        prot = prot_fr_1
        return "".join(prot)
    elif "*" not in prot_fr_2:
        prot = prot_fr_2
        return "".join(prot)
    elif "*" not in prot_fr_3:
        prot = prot_fr_3
        return "".join(prot)


def split_regions(sequence):
    print("")


def main(infile, outpath, name):

    # get absolute paths
    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)
    name = name + "_aligned.fasta"
    outfile = os.path.join(outpath, name)
    # print(outfile)

    # read in fasta file
    in_seqs_d = fasta_to_dct_rev(infile)
    reference = "ATGAGAGTGAGGGGGATACTGAGGAATTGGCCACAATGGTGGATATGGGGCATCTTAGGCTTTTGGATGTTAATGATTTGTAGTGTGGTGGGAAACTTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAAGAAGCAAAAACTACTCTATTCTGTGCATCAGATGCTAAAGCATATGAGAAAGAAGTGCATAATGTCTGGGCTACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAATAGTTTTGGAAAATGTAACAGAAAATTTTAACATGTGGAAAAATGACATGGTGGATCAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAGTTGACCCCACTCTGTGTCACTTTAAATTGTACAAATGCTACTATTAATAATACTTACTACAATAATAGCATGAATGAAGAAATAAAAAATTGCTCTTTCAATACAACCACAGAAATAAGAGATAAGAAACAGAAAGCGTATGCACTTTTTTATAGACCTGATATAGTACCACTTAATGAGAATAATAGTGAGTATATATTAATAAATTGTAATACCTCAACCATAACACAAGCCTGTCCAAAGGTCACTTTTGACCCAATTCCTATACATTATTGTGCTCCAGCTGGTTATGCGATTCTAAAGTGTAATAATAAGACATTCAATGGGACAGGACCATGCAATAATGTCAGCACAGTACAATGTACACATGGAATTAAGCCAGTGGTATCAACTCAACTACTGTTAAATGGTAGCCTAGCAGAAGAAGAGATAATAATTAGATCTGAAAATCTGACAGACAATGCCAAAACAATAATAGTACATCTTAATGAATCTGTAGAAATTGTGTGTACAAGACCCAACAATAATACAAGAAAAAGTATAAGGATAGGACCAGGACAAACATTCTATGCAACAGGTGACATAATAGGAGACATAAGACAAGCACATTGTAACATTAGTAAAAAAAAATGGAATAAAACTTTAGAAAAGGTAAAGGAAAAATTAAAAGAACACTTCCCTAATAAAACAATAAAATTTGAACCATCCTCAGGAGGGGACCTAGAAATTACAACACATAGCTTTAATTGTAGAGGAGAATTTTTCTATTGCAATACATCAAAACTGTTTAATAATACATACAATAGTACAACAAATACAAATGCAACCATCACACTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGGAGGTAGGACGAGCAATGTATGCCCCTCCCATTGCAGGAAACATAACATGTAACTCAAATATCACAGGACTACTATTGACACGTGATGGAGGAAAAAATAACACAAATAACACAGAGACATTCAGACCTGGAGGAGGAAATATGAAGGACAATTGGAGAAGTGAATTATATAAATATAAAGTGGTAGAAATTAAGCCATTGGGAATAGCACCCACTAAGGCAAAAAGGAGAGTGGTGGAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTGTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCGGCGTCAATAACGCTGACGGTACAGGCCAGACAATTGTTGTCTGGTATAGTGCAACAGCAAAGCAATTTGCTGAGAGCTATAGAGGCGCAACAGCATATGTTGCAACTCACAGTCTGGGGCATTAAGCAGCTCCAGGCAAGAGTCCTGGCTATAGAAAGATACCTAAAGGATCAACAGCTCCTAGGAATTTGGGGCTGCTCTGGAAAACTCATCTGCACCACTGCTGTGCCTTGGAACTCCAGTTGGAGTAATAAATCTCAAGAAGATATTTGGGATAACATGACCTGGATGCAGTGGGATAGAGAAATTAGTAATTACACAAACACAATATACAGGTTGCTTGAAGAATCGCAAAACCAGCAGGAGAAAAATGAAAAAGATTTATTAGCATTGGACAGTTGGAAAAATCTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAAATATTCATAATGATAGTAGGAGGCTTGATAGGTTTAAGAATAATTTTTGCTGTGCTTTCTATAGTGAATAGAGTTAGGCAGGGATACTCACCTTTGTCGTTTCAGACCCTTACCCCAAACCCGAGGGGACCCGACAGGCTCGGAAGAATCGAAGAAGAAGGTGGAGAGCAAGACAGAGACAGATCCATTCGATTAGTGAGCGGATTCTTAGCACTTGCCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCAATTGAGAGACTTCATATTGATTGCAGCGAGAGCAGTGGAACTTCTGGGACGCAGCAGTCTCAGGGGACTACAGAGGGGGTGGGAAGCCCTTAAGTATCTGGGAAGTCTTGTGCAGTATTGGGGTCTGGAACTAAAAAAGAGTGCTATTAGTCTGCTTGATACCATAGCAATAGCAGTAGCTGAAGGAACAGATAGGATTATAGAATTAATACAAAGAATTTGTAGAGCTATCCGCAACATACCTAGAAGAATAAGACAGGGCTTTGAAGCAGCTTTGCTATAA"

    # generate seq_code to seq name list lookup dictionary
    first_look_up_d = collections.defaultdict(list)
    first_seq_code_d = collections.defaultdict(str)
    for i, (seq, names_list) in enumerate(in_seqs_d.items()):
        unique_id = str(i).zfill(4)
        first_look_up_d[unique_id] = names_list
        first_seq_code_d[seq] = unique_id

    unit_test_d = collections.defaultdict
    # translate sequences
    for seq, code in first_seq_code_d.items():
        s_name = first_look_up_d[code][0]
        padded_sequence, frame = prelim_align(seq, reference)

        prot_seq = translate_dna(padded_sequence, frame)


    print("Aligning completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path for the output file', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='The name for the output file', required=True)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name

    main(infile, outpath, name)
