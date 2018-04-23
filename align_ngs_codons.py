#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import collections
from itertools import groupby
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA
from skbio import Protein
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


def pairwise_align_dna(sequence, reference, name):
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


def posnumcalc(hxb2seq, start):
    """
    Calculates the positional numbering relative to hxb2
    :param hxb2seq: (str) hxb2 protein sequence
    :param start:  (int) start amino acid position for hxb2
    :return: (list) list of position numbers [1, 2, 3, 4, 4.01, 4.02, 5]
    """
    pos_num = []
    n = start
    s = 0.01
    m = len(hxb2seq) - 1
    for i, resi in enumerate(hxb2seq):
        if i == 0 and resi == '-':
            print("Can't start numbering. HXB2 starts with a gap. Use a longer HXB2 sequence for the numbering")
        if i != m:
            if resi != '-' and hxb2seq[i+1] != '-':
                pos_num.append(n)
                n += 1
            elif resi != '-' and hxb2seq[i+1] == '-':
                g = n
                pos_num.append(g)
            elif resi == '-' and hxb2seq[i+1] == '-':
                g = n + s
                pos_num.append(g)
                s += 0.01
            elif resi == '-' and hxb2seq[i+1] != '-':
                g = n + s
                pos_num.append(g)
                s = 0.01
                n += 1
        else:
            if resi != '-':
                pos_num.append(n)
            elif resi == '-':
                g = n + s
                pos_num.append(g)

    return pos_num


def find_cons_var_regions(prot_sequence):
    """
    function to identify variable and conserved region of a sequence (envelope)
    :param prot_sequence: (str) a protein query sequence
    :return: (dict) dict of cons regions, (dict) dict of var regions
    """
    hxb2_prot = "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVE" \
                "QMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQA" \
                "CPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKR" \
                "IRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNT" \
                "EGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREK" \
                "RAVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNA" \
                "SWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSP" \
                "LSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLN" \
                "ATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILLX"

    blosum62 = {
        'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1,
              'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': 0, '*': -4},
        'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2,
              'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': -1, 'Z': 0, 'X': -1,
              '*': -4},
        'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0,
              'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3, 'B': 3, 'Z': 0, 'X': -1, '*': -4},
        'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1,
              'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4},
        'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3,
              'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1, 'B': -3, 'Z': -3, 'X': -2,
              '*': -4},
        'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1,
              'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2, 'B': 0, 'Z': 3, 'X': -1, '*': -4},
        'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1,
              'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4},
        'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2,
              'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3, 'B': -1, 'Z': -2, 'X': -1,
              '*': -4},
        'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1,
              'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3, 'B': 0, 'Z': 0, 'X': -1, '*': -4},
        'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3,
              'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3, 'B': -3, 'Z': -3, 'X': -1, '*': -4},
        'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2,
              'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1, 'B': -4, 'Z': -3, 'X': -1, '*': -4},
        'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5,
              'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 1, 'X': -1, '*': -4},
        'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1,
              'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1, 'B': -3, 'Z': -1, 'X': -1, '*': -4},
        'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3,
              'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1, 'B': -3, 'Z': -3, 'X': -1, '*': -4},
        'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3,
              'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2, 'B': -2, 'Z': -1, 'X': -2,
              '*': -4},
        'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0,
              'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2, 'B': 0, 'Z': 0, 'X': 0, '*': -4},
        'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1,
              'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0, 'B': -1, 'Z': -1, 'X': 0, '*': -4},
        'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2,
              'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3, 'B': -4, 'Z': -3, 'X': -2,
              '*': -4},
        'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2,
              'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1, 'B': -3, 'Z': -2, 'X': -1, '*': -4},
        'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2,
              'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4, 'B': -3, 'Z': -2, 'X': -1, '*': -4},
        'B': {'A': -2, 'R': -1, 'N': 3, 'D': 4, 'C': -3, 'Q': 0, 'E': 1, 'G': -1, 'H': 0, 'I': -3, 'L': -4, 'K': 0,
              'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 4, 'Z': 1, 'X': -1, '*': -4},
        'Z': {'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -3, 'Q': 3, 'E': 4, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1,
              'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2, 'B': 1, 'Z': 4, 'X': -1, '*': -4},
        'X': {'A': 0, 'R': -1, 'N': -1, 'D': -1, 'C': -2, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1, 'K': -1,
              'M': -1, 'F': -1, 'P': -2, 'S': 0, 'T': 0, 'W': -2, 'Y': -1, 'V': -1, 'B': -1, 'Z': -1, 'X': -1, '*': -4},
        '*': {'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4,
              'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'Z': -4,
              'X': -4, '*': 1},
        '-': {'A': -4, 'R': -4, 'N': -4, 'D': -4, 'C': -4, 'Q': -4, 'E': -4, 'G': -4, 'H': -4, 'I': -4, 'L': -4,
              'K': -4, 'M': -4, 'F': -4, 'P': -4, 'S': -4, 'T': -4, 'W': -4, 'Y': -4, 'V': -4, 'B': -4, 'Z': -4,
              'X': -4, '*': 1}
    }

    blosum80 = {
        'A': {'A': 5, 'R': -2, 'N': -2, 'D': -2, 'C': -1, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -2, 'L': -2, 'K': -1,
              'M': -1, 'F': -3, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0, 'B': -2, 'Z': -1, 'X': -1, '*': -6},
        'R': {'A': -2, 'R': 6, 'N': -1, 'D': -2, 'C': -4, 'Q': 1, 'E': -1, 'G': -3, 'H': 0, 'I': -3, 'L': -3, 'K': 2,
              'M': -2, 'F': -4, 'P': -2, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': -2, 'Z': 0, 'X': -1,
              '*': -6},
        'N': {'A': -2, 'R': -1, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': -1, 'G': -1, 'H': 0, 'I': -4, 'L': -4, 'K': 0,
              'M': -3, 'F': -4, 'P': -3, 'S': 0, 'T': 0, 'W': -4, 'Y': -3, 'V': -4, 'B': 4, 'Z': 0, 'X': -1, '*': -6},
        'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -4, 'Q': -1, 'E': 1, 'G': -2, 'H': -2, 'I': -4, 'L': -5, 'K': -1,
              'M': -4, 'F': -4, 'P': -2, 'S': -1, 'T': -1, 'W': -6, 'Y': -4, 'V': -4, 'B': 4, 'Z': 1, 'X': -1, '*': -6},
        'C': {'A': -1, 'R': -4, 'N': -3, 'D': -4, 'C': 9, 'Q': -4, 'E': -5, 'G': -4, 'H': -4, 'I': -2, 'L': -2, 'K': -4,
              'M': -2, 'F': -3, 'P': -4, 'S': -2, 'T': -1, 'W': -3, 'Y': -3, 'V': -1, 'B': -4, 'Z': -4, 'X': -1,
              '*': -6},
        'Q': {'A': -1, 'R': 1, 'N': 0, 'D': -1, 'C': -4, 'Q': 6, 'E': 2, 'G': -2, 'H': 1, 'I': -3, 'L': -3, 'K': 1,
              'M': 0, 'F': -4, 'P': -2, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -3, 'B': 0, 'Z': 3, 'X': -1, '*': -6},
        'E': {'A': -1, 'R': -1, 'N': -1, 'D': 1, 'C': -5, 'Q': 2, 'E': 6, 'G': -3, 'H': 0, 'I': -4, 'L': -4, 'K': 1,
              'M': -2, 'F': -4, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 1, 'Z': 4, 'X': -1, '*': -6},
        'G': {'A': 0, 'R': -3, 'N': -1, 'D': -2, 'C': -4, 'Q': -2, 'E': -3, 'G': 6, 'H': -3, 'I': -5, 'L': -4, 'K': -2,
              'M': -4, 'F': -4, 'P': -3, 'S': -1, 'T': -2, 'W': -4, 'Y': -4, 'V': -4, 'B': -1, 'Z': -3, 'X': -1,
              '*': -6},
        'H': {'A': -2, 'R': 0, 'N': 0, 'D': -2, 'C': -4, 'Q': 1, 'E': 0, 'G': -3, 'H': 8, 'I': -4, 'L': -3, 'K': -1,
              'M': -2, 'F': -2, 'P': -3, 'S': -1, 'T': -2, 'W': -3, 'Y': 2, 'V': -4, 'B': -1, 'Z': 0, 'X': -1, '*': -6},
        'I': {'A': -2, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -3, 'E': -4, 'G': -5, 'H': -4, 'I': 5, 'L': 1, 'K': -3,
              'M': 1, 'F': -1, 'P': -4, 'S': -3, 'T': -1, 'W': -3, 'Y': -2, 'V': 3, 'B': -4, 'Z': -4, 'X': -1, '*': -6},
        'L': {'A': -2, 'R': -3, 'N': -4, 'D': -5, 'C': -2, 'Q': -3, 'E': -4, 'G': -4, 'H': -3, 'I': 1, 'L': 4, 'K': -3,
              'M': 2, 'F': 0, 'P': -3, 'S': -3, 'T': -2, 'W': -2, 'Y': -2, 'V': 1, 'B': -4, 'Z': -3, 'X': -1, '*': -6},
        'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -4, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -3, 'K': 5,
              'M': -2, 'F': -4, 'P': -1, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': -1, 'Z': 1, 'X': -1,
              '*': -6},
        'M': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -2, 'Q': 0, 'E': -2, 'G': -4, 'H': -2, 'I': 1, 'L': 2, 'K': -2,
              'M': 6, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -2, 'V': 1, 'B': -3, 'Z': -2, 'X': -1, '*': -6},
        'F': {'A': -3, 'R': -4, 'N': -4, 'D': -4, 'C': -3, 'Q': -4, 'E': -4, 'G': -4, 'H': -2, 'I': -1, 'L': 0, 'K': -4,
              'M': 0, 'F': 6, 'P': -4, 'S': -3, 'T': -2, 'W': 0, 'Y': 3, 'V': -1, 'B': -4, 'Z': -4, 'X': -1, '*': -6},
        'P': {'A': -1, 'R': -2, 'N': -3, 'D': -2, 'C': -4, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': -4, 'L': -3,
              'K': -1, 'M': -3, 'F': -4, 'P': 8, 'S': -1, 'T': -2, 'W': -5, 'Y': -4, 'V': -3, 'B': -2, 'Z': -2, 'X': -1,
              '*': -6},
        'S': {'A': 1, 'R': -1, 'N': 0, 'D': -1, 'C': -2, 'Q': 0, 'E': 0, 'G': -1, 'H': -1, 'I': -3, 'L': -3, 'K': -1,
              'M': -2, 'F': -3, 'P': -1, 'S': 5, 'T': 1, 'W': -4, 'Y': -2, 'V': -2, 'B': 0, 'Z': 0, 'X': -1, '*': -6},
        'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -2, 'K': -1,
              'M': -1, 'F': -2, 'P': -2, 'S': 1, 'T': 5, 'W': -4, 'Y': -2, 'V': 0, 'B': -1, 'Z': -1, 'X': -1, '*': -6},
        'W': {'A': -3, 'R': -4, 'N': -4, 'D': -6, 'C': -3, 'Q': -3, 'E': -4, 'G': -4, 'H': -3, 'I': -3, 'L': -2,
              'K': -4, 'M': -2, 'F': 0, 'P': -5, 'S': -4, 'T': -4, 'W': 11, 'Y': 2, 'V': -3, 'B': -5, 'Z': -4, 'X': -1,
              '*': -6},
        'Y': {'A': -2, 'R': -3, 'N': -3, 'D': -4, 'C': -3, 'Q': -2, 'E': -3, 'G': -4, 'H': 2, 'I': -2, 'L': -2, 'K': -3,
              'M': -2, 'F': 3, 'P': -4, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -2, 'B': -3, 'Z': -3, 'X': -1, '*': -6},
        'V': {'A': 0, 'R': -3, 'N': -4, 'D': -4, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -4, 'I': 3, 'L': 1, 'K': -3,
              'M': 1, 'F': -1, 'P': -3, 'S': -2, 'T': 0, 'W': -3, 'Y': -2, 'V': 4, 'B': -4, 'Z': -3, 'X': -1, '*': -6},
        'B': {'A': -2, 'R': -2, 'N': 4, 'D': 4, 'C': -4, 'Q': 0, 'E': 1, 'G': -1, 'H': -1, 'I': -4, 'L': -4, 'K': -1,
              'M': -3, 'F': -4, 'P': -2, 'S': 0, 'T': -1, 'W': -5, 'Y': -3, 'V': -4, 'B': 4, 'Z': 0, 'X': -1, '*': -6},
        'Z': {'A': -1, 'R': 0, 'N': 0, 'D': 1, 'C': -4, 'Q': 3, 'E': 4, 'G': -3, 'H': 0, 'I': -4, 'L': -3, 'K': 1,
              'M': -2, 'F': -4, 'P': -2, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3, 'B': 0, 'Z': 4, 'X': -1, '*': -6},
        'X': {'A': -1, 'R': -1, 'N': -1, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': -1, 'L': -1,
              'K': -1, 'M': -1, 'F': -1, 'P': -1, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': -1, 'B': -1, 'Z': -1,
              'X': -1, '*': -6},
        '*': {'A': -6, 'R': -6, 'N': -6, 'D': -6, 'C': -6, 'Q': -6, 'E': -6, 'G': -6, 'H': -6, 'I': -6, 'L': -6,
              'K': -6, 'M': -6, 'F': -6, 'P': -6, 'S': -6, 'T': -6, 'W': -6, 'Y': -6, 'V': -6, 'B': -6, 'Z': -6,
              'X': -6, '*': 1}
    }

    # index = +1 for 0 based indexing of positions
    hxb2_cons_regions_index_d = {"C1": (0, 135),
                                 "C2": (151, 394),
                                 "C3": (408, 463),
                                 "C4": (466, 856),}

    hxb2_var_regions_index_d = {"V1": (135, 151),
                                 "V2": (394, 408),
                                 "V3": (464, 466),}

    alignment, score, start_end_positions = local_pairwise_align_ssw(Protein(prot_sequence), Protein(hxb2_prot),
                                                                     gap_open_penalty=8, gap_extend_penalty=2,
                                                                     match_score=4, mismatch_score=-1.5,
                                                                     substitution_matrix=blosum80)

    # get the query and ref aligned seqs
    seq_align = str(alignment[0])
    ref_align = str(alignment[1])
    print(">{0}\n{1}".format(name, seq_align))
    print(">ref\n{}".format(ref_align))

    # get start position for reference and query
    ref_start = start_end_positions[1][0]
    seq_start = start_end_positions[0][0]
    if seq_start > 0:
        sys.exit("Query sequence was truncated during alignment\nExiting")

    # get hxb2 numbering
    hxb2_numbering = posnumcalc(ref_align, ref_start)
    print(hxb2_numbering)

    # extract the cons and var regions to dicts
    cons_d = collections.defaultdict(str)
    var_d = collections.defaultdict(str)

    for cons_region, index_tup in hxb2_cons_regions_index_d.items():
        start = index_tup[0]
        end = index_tup[1]
        if start in hxb2_numbering and end in hxb2_numbering:
            # check that indexing is correct! not trimming off a resi at the end...
            cons_region_slice = seq_align[start:end]
            cons_d[cons_region] = cons_region_slice
        elif start not in hxb2_numbering and end in hxb2_numbering:
            cons_region_slice = seq_align[:end]
            cons_d[cons_region] = cons_region_slice
        elif start in hxb2_numbering and end not in hxb2_numbering:
            cons_region_slice = seq_align[start:]
            cons_d[cons_region] = cons_region_slice
        else:
            sys.exit("HXB2 coordinates for {0} not found\nHXb2 numbering was {}".format(cons_region, hxb2_numbering))

    for var_region, index_tup in hxb2_var_regions_index_d.items():
        start = index_tup[0]
        end = index_tup[1]
        if start in hxb2_numbering and end in hxb2_numbering:
            # check that indexing is correct! not trimming off a resi at the end...
            var_region_slice = seq_align[start:end]
            var_d[var_region] = var_region_slice
        elif start not in hxb2_numbering and end in hxb2_numbering:
            var_region_slice = seq_align[:end]
            var_d[var_region] = var_region_slice
        elif start in hxb2_numbering and end not in hxb2_numbering:
            var_region_slice = seq_align[start:]
            var_d[var_region] = var_region_slice
        else:
            sys.exit("HXB2 coordinates for {0} not found\nHXb2 numbering was {}".format(var_region, hxb2_numbering))

    return cons_d, var_d


def write_cons_var_regions_to_file(cons_d, var_d):
    print("")
    cons_fn = ''
    with open(cons_fn, 'w') as handle:
        handle.write('')


    return cons_fn


def call_aligner(prot_seq_d):
    """
    Takes a dict of protein sequences, writes them to a temp file and aligns them with mafft.
    Aligned file is read back in and returned as a dictionary
    :param prot_seq_d: (dict) dict of protein sequences: key = sequence, value = ID code
    :return: (dict) dictionary of aligned protein sequences: key = sequence, value = ID code
    """
    region_aligned_d = collections.defaultdict()

    # write dict to temp file

    # align temp file

    # read aligned file back in as dict

    # remove temp file

    return region_aligned_d


def join_regions(cons_regions, var_regions):
    print("")
    joined_d = collections.defaultdict(int)

    return joined_d


def backtranslate(padded_dna_input, prot_align):
    print("")
    dna_align = ''
    resi_count = 0
    for resi in prot_align:
        if resi == '-':
            dna_align += '---'
        else:
            dna_align += padded_dna_input[(resi_count * 3):((resi_count * 3) + 3)]
            resi_count += 1

    return dna_align


def main(infile, ref, outpath, name):

    # get absolute paths
    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)
    name = name + "_aligned.fasta"
    outfile = os.path.join(outpath, name)

    get_script_path = os.path.realpath(__file__)
    script_folder = os.path.split(get_script_path)[0]
    script_folder = os.path.abspath(script_folder)
    hxb2_file = os.path.join(script_folder, "HXB2_seqs.fasta")
    hxb2 = list(fasta_to_dct(hxb2_file).values())[0]
    hxb2_trans = translate_dna(hxb2.replace("-", ""))

    # read in fasta file and reference
    in_seqs_d = fasta_to_dct_rev(infile)
    reference = list(fasta_to_dct(ref).values())[0]
    reference_trans = translate_dna(reference.replace("-", ""))

    # generate seq_code to seq name list lookup dictionary
    first_look_up_d = collections.defaultdict(list)
    first_seq_code_d = collections.defaultdict(str)
    for i, (seq, names_list) in enumerate(in_seqs_d.items()):
        unique_id = str(i).zfill(4)
        first_look_up_d[unique_id] = names_list
        first_seq_code_d[seq] = unique_id

    # translate sequences

    cons_regions_dct = collections.defaultdict()
    var_regions_dct = collections.defaultdict()
    prot_seq_d = collections.defaultdict(dict)
    for seq, code in first_seq_code_d.items():
        s_name = first_look_up_d[code][0]
        # print(s_name)
        seq_align, ref_align, frame = pairwise_align_dna(seq, reference, s_name)
        padded_sequence = gap_padding(seq_align, ref_align, frame, seq, s_name)
        prot_seq = translate_dna(padded_sequence)
        prot_ref = translate_dna(reference)
        cons_regions_dct, var_regions_dct = find_cons_var_regions(prot_seq) # split apart? use .update()

    write_cons_var_regions_to_file()

    align_prot_d = call_aligner(prot_seq_d)

    # for align_prot, padd_code_d in align_prot_d.items():
    #     dna_aligned = backtranslate(padded_sequence, align_prot)

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
