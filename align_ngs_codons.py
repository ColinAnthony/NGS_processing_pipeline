#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import collections
from itertools import groupby
import random
import string
import subprocess
from subprocess import DEVNULL
import regex
import seqanpy
import pandas as pd


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
        v = v.replace("-", "")
        new_key = k.replace(" ", "_")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[new_key] = v.upper()

    return dct


def fasta_to_dct_keep_gap(file_name):
    """
    :param file_name: The fasta formatted file to read from.
    :return: a dictionary of the contents of the file name given. Dictionary in the format:
    {sequence_id: sequence_string, id_2: sequence_2, etc.}
    """
    dct = collections.defaultdict(str)
    my_gen = py3_fasta_iter(file_name)
    for k, v in my_gen:
        v = v
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
        new_key = k.replace("-", "")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[str(v).replace(" ", "_").upper()].append(new_key)

    return dct


def pairwise_align_dna(sequence, reference, regex_complied):
    """
    Pairwise align sequence to reference, to find reading frame and frame-shift in/dels
    :param sequence: (str) a query DNA sequence
    :param reference: (str) a reference DNA sequence (must start in reading frame 1)
    :param regex_complied: (regex_obj) a compiled regex pattern
    :return: (str) aligned query sequence, (str) aligned ref sequence, (int) reading frame for query sequence
    """

    # do overlap pairwise alignment to not get truncated query sequence
    overlap = seqanpy.align_overlap(sequence, reference, band=-1, score_match=4, score_mismatch=-1, score_gapext=-2,
                                    score_gapopen=-15)
    overlap = list(overlap)

    seq_align = overlap[1]
    ref_align = overlap[2]
    # print(">sqseq1\n{}\n".format(seq_align))
    # print(">sqref1\n{}\n".format(ref_align))

    # get start position in the seq, if not starting at index 0
    if seq_align[0] == '-':
        seq_start = regex_complied.search(seq_align).end()
    else:
        seq_start = 0

    # get end position in the seq, if not starting at index 0
    if seq_align[-1] == '-':
        # reverse the string and find first non-gap character, mult match.end by -1 to get non-reversed index
        seq_end = (regex_complied.search(seq_align[::-1]).end()) * -1
    else:
        seq_end = None

    # ref start will be 0 for align_overlap
    ref_start = 0

    # calculate reading frame (reference must start in frame 0)
    frame = (seq_start - ref_start) % 3

    # truncate the overlap alignment to the region of interest
    seq_align = seq_align[seq_start:seq_end]
    ref_align = ref_align[seq_start:seq_end]
    # print(">sqseq2\n{}\n".format(seq_align))
    # print(">sqref2\n{}\n".format(ref_align))
    # input("enter")
    return seq_align, ref_align, frame


def gap_padding(seq_align, ref_align, frame, regex_complied):
    """
    pads sequence with gaps for indels and to set reading frame to frame 1
    :param seq_align: (str) an aligned query sequence
    :param ref_align: (str) the corresponding aligned reference sequence
    :param frame: (int) the reading frame for the query sequence
    :param regex_complied: (regex_obj) a compiled regex pattern
    :return: (str) a gap-padded query sequence to correct for indels and set reading frame to frame 1
    """
    # convert query seq to list to allow mutability
    new_seq = list(seq_align)
    indel_gap_fix_master = []

    # get index and len of all del gaps to fix (gaps in query)
    all_gap_positions_seq = regex_complied.finditer(seq_align)

    for gap_obj in all_gap_positions_seq:
        gap_start = gap_obj.start()
        gap = gap_obj.captures()[0]
        gap_len = len(gap)
        gap_shift = gap_len % 3
        if gap_shift == 1:
            # changed this from 2 to 1 testing
            new_gap = "-" * 1
        elif gap_shift == 2:
            # changed this from 1 to 2 testing
            new_gap = "-" * 2
        else:
            new_gap = ""

        gap_len_in_seq = gap_len
        indel_gap_fix_master.append((gap_start, gap_len_in_seq, new_gap))

    # get index and len of all ins gaps to fix (gaps in ref)
    all_gap_positions_ref = regex_complied.finditer(ref_align)

    for gap_obj in all_gap_positions_ref:
        gap_start = gap_obj.start()
        gap = gap_obj.captures()[0]
        gap_len = len(gap)
        gap_shift = gap_len % 3
        if gap_shift == 1:
            new_gap = "-" * 2
        elif gap_shift == 2:
            new_gap = "-" * 1
        else:
            new_gap = ""
        if new_gap != "":
            gap_len_in_seq = 0
            indel_gap_fix_master.append((gap_start, gap_len_in_seq, new_gap))

    # sort the list of all gaps to insert by gap start pos
    indel_gap_fix_master = sorted(indel_gap_fix_master)
    # print(indel_gap_fix_master)

    # fix all gaps in order
    index_adjust = 0
    for gap_tuple in indel_gap_fix_master:
        # print("index", gap_tuple[0])
        idx = gap_tuple[0] + index_adjust
        # print("new index", idx)
        new_gap = gap_tuple[2]
        new_gap_len = len(new_gap)
        old_gap_len = gap_tuple[1]
        if old_gap_len == new_gap_len:
            continue
        if new_gap == "":
            # remove unnecessary gap, gap len multiple of 3
            # print("old gap :", old_gap_len)
            # print("removing:", to_remove)
            del new_seq[idx:(idx + old_gap_len)]
            index_adjust -= old_gap_len

        elif old_gap_len == 0:
            # add gap for frame-shift insertion
            # print("inerting:", new_gap)

            if new_gap_len == 1:
                # todo: should this be "idx+2" as well??
                new_seq.insert(idx, new_gap)
                index_adjust += 1
            elif new_gap_len == 2:
                new_seq.insert(idx+2, new_gap)
                index_adjust += 1
                # index_adjust += len(new_gap) -1

        else:
            # add gap for frame-shift deletion
            if old_gap_len == 1:
                # print("inerting:", new_gap)
                new_seq.insert(idx, new_gap)
                index_adjust += 1
            else:
                gap_change = old_gap_len - new_gap_len
                # print("old gap", old_gap_len)
                # print("truncating:", new_seq[idx:(idx + gap_change)])
                del new_seq[idx:(idx + gap_change)]
                if old_gap_len > new_gap_len:
                    index_adjust -= gap_change
                else:
                    index_adjust += gap_change

        # print("inx_adj :", index_adjust)

    # pad query to be in reading frame 1
    if frame != 0:
        lead_gap = "-" * frame
        new_seq.insert(0, lead_gap)

    # convert to string and return
    padded_seq = "".join(new_seq)
    # print(">padded\n{}".format(padded_seq))
    if seq_align.replace("-", "") != padded_seq.replace("-", ""):
        print("somethig went wrong, input and padded sequences are different")
        print(">input\n{1}".format(name, seq_align))
        print(">padded\n{1}".format(name, padded_seq))

    return padded_seq


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
    # convert terminal stop codon to X so that it doesn't get truncated by mafft
    if prot[-1] == "*":
        prot[-1] = "X"

    return "".join(prot)


def prot_pairwise_align(prot_sequence, ref_prot, regex_complied):
    """
    function to pairwise align protein sequence to ref
    :param prot_sequence: (str) a protein query sequence
    :param ref_prot: (str) a protein translation of the reference
    :param regex_complied: (regex_obj) a compiled regex pattern
    :return: (dict) dict of cons regions, (dict) dict of var regions
    """

    overlap = seqanpy.align_overlap(prot_sequence, ref_prot, band=-1, score_match=4, score_mismatch=-1, score_gapext=-2,
                                    score_gapopen=-15)
    overlap = list(overlap)

    seq_align = overlap[1]
    ref_align = overlap[2]
    # print(">seq_align\n{}".format(seq_align))
    # print(">ref_align\n{}".format(ref_align))

    if seq_align[0] == '-':
        seq_start = regex_complied.search(seq_align).end()
    else:
        seq_start = 0

    # get end position in the seq, if not starting at index 0
    if seq_align[-1] == "X":
        seq_end = None
    elif seq_align[-1] == '-':
        seq_end = (regex_complied.search(seq_align[::-1]).end()) * - 1
    else:
        seq_end = - 1

    seq_align_new = seq_align[seq_start:seq_end]
    ref_align_new = ref_align[seq_start:seq_end]

    start_end_positions = [(seq_start, seq_end), (seq_start, seq_end)]

    alignment_d = {"inseq": prot_sequence,
                   "query": seq_align_new,
                   "reference": ref_align_new,
                   "start_end_positions": start_end_positions,
                   }
    # print(">seq_align\n{}".format(seq_align_new[-5:]))
    # print(">ref_align\n{}".format(ref_align_new))

    return alignment_d
def set_cons_var_regions(ref_type):
    """
    contains the coordinates for the start and end of each conserved and variable region
    :ref_type: (str) the reference coordinates to use
    :return:
    """
    # todo: change this to regex to find region
    if ref_type != "HXB2":
        ref_cons_regions_index_d = {"C1": (0, 130),
                                    "C2": (156, 177),
                                    "C3": (185, 383),
                                    "C4": (395, 440),
                                    "C5": (447, 844),
                                    }

        ref_var_regions_index_d = {"V1": (130, 156),
                                   "V2": (177, 185),
                                   "V3": (383, 395),
                                   "V4": (440, 447),
                                   }

        full_order = ["C1", "V1", "C2", "V2", "C3", "V3", "C4", "V4", "C5"]

    else:
        ref_cons_regions_index_d = {"C1": (0, 131),
                                    "C2": (152, 184),
                                    "C3": (190, 394),
                                    "C4": (410, 459),
                                    "C5": (465, 856),
                                    }

        ref_var_regions_index_d = {"V1": (131, 152),
                                   "V2": (184, 190),
                                   "V3": (394, 410),
                                   "V4": (459, 465),
                                   }

        full_order = ["C1", "V1", "C2", "V2", "C3", "V3", "C4", "V4", "C5"]

    return ref_cons_regions_index_d, ref_var_regions_index_d, full_order
def posnumcalc(ref_seq, start):
    """
    Calculates the positional numbering relative to hxb2
    :param ref_seq: (str) hxb2 protein sequence
    :param start:  (int) start amino acid position for hxb2
    :return: (list) list of position numbers [1, 2, 3, 4, 4.01, 4.02, 5]
    """
    ref_resi_len = len(ref_seq.replace("-", ""))
    pos_nums = list(range(start, (start + ref_resi_len) + 1))

    return pos_nums
def get_cons_regions(prot_align_d, end_reg, ref_cons_regions_index_d, ref_numbering):
    """
    function to extract conserved regions of a sequence (envelope)
    :param prot_align_d: (dict) dict of pairwise alignment object
    :param end_reg: (str) end region key
    :param ref_cons_regions_index_d: (dict) of cons regions and their indexes in the reference prot seq
    :param ref_numbering: (list) list of reference sequence position numbers relative to start of coding region
    :return: (dict) dict of cons regions
    """
    ref_start = prot_align_d["start_end_positions"][1][0]
    seq_align = prot_align_d["query"]
    ref_align = prot_align_d["reference"]

    # get ref numbering
    # ref_numbering = posnumcalc(ref_align, ref_start)

    # extract the cons regions to dicts
    cons_d = collections.defaultdict(str)

    for cons_region, index_tup in ref_cons_regions_index_d.items():
        start = index_tup[0]
        end = index_tup[1]
        if start in ref_numbering and end in ref_numbering:
            if cons_region == end_reg:
                cons_region_slice = seq_align[ref_numbering.index(start):]
                cons_d[cons_region] = cons_region_slice
            else:
                cons_region_slice = seq_align[ref_numbering.index(start):ref_numbering.index(end)]
                cons_d[cons_region] = cons_region_slice.replace("-", "")

        elif start not in ref_numbering and end in ref_numbering:
            cons_region_slice = seq_align[:ref_numbering.index(end)]
            cons_d[cons_region] = cons_region_slice.replace("-", "")

        elif start in ref_numbering and end not in ref_numbering:
            cons_region_slice = seq_align[ref_numbering.index(start):]
            cons_d[cons_region] = cons_region_slice.replace("-", "")

        else:
            print("{} region not present".format(cons_region))

    return cons_d
def get_var_regions(prot_align_d, end_reg, ref_var_regions_index_d, ref_numbering):
    """
    function to extract conserved regions of a sequence (envelope)
    :param prot_align_d: (dict) dict of pairwise alignment object
    :param end_reg: (str) end region key
    :param ref_var_regions_index_d: (dict) of variable regions and their indexes in the reference prot seq
    :param ref_numbering: (list) list of reference sequence position numbers relative to start of coding region
    :return: (dict) dict of var regions
    """
    ref_start = prot_align_d["start_end_positions"][1][0]
    seq_align = prot_align_d["query"]
    ref_align = prot_align_d["reference"]

    # get hxb2 numbering
    # ref_numbering = posnumcalc(ref_align, ref_start)

    # extract the var regions to dicts
    var_d = collections.defaultdict(str)

    for var_region, index_tup in ref_var_regions_index_d.items():
        start = index_tup[0]
        end = index_tup[1]
        if start in ref_numbering and end in ref_numbering:
            if var_region == end_reg:
                var_region_slice = seq_align[ref_numbering.index(start):]
                var_d[var_region] = var_region_slice
            else:
                var_region_slice = seq_align[ref_numbering.index(start):ref_numbering.index(end)]
                var_d[var_region] = var_region_slice.replace("-", "")

        elif start not in ref_numbering and end in ref_numbering:
            var_region_slice = seq_align[:ref_numbering.index(end)]
            var_d[var_region] = var_region_slice.replace("-", "")

        elif start in ref_numbering and end not in ref_numbering:
            var_region_slice = seq_align[ref_numbering.index(start):]
            var_d[var_region] = var_region_slice.replace("-", "")
        else:
            print("{} region not present".format(var_region))

    return var_d
def find_start_end_regions(alignment_d, ref_cons_regions_index_d, ref_var_regions_index_d, full_order, bounds_lower):
    """
    function to take a pairwise align obj in dict format and find the start and end regions for cons and var boundaries
    :param alignment_d: (dict) dict of pairwise alignment object
    :param ref_cons_regions_index_d: (dict) of conserved regions
    :param ref_var_regions_index_d: (dict) of variable regions
    :param full_order: (list) a list of all the conserved and variable regions in sequential order
    :param bounds_lower: (int) start of reference sequence in amino acid numbering
    :return: (str) start and end region keys
    """
    # seq_align = alignment_d["query"]
    ref_align = alignment_d["reference"]
    ref_start = alignment_d["start_end_positions"][1][0]

    # get hxb2 numbering
    if bounds_lower is None:
        bounds_lower = 0
    ref_start_idx = ref_start + bounds_lower
    ref_numbering = posnumcalc(ref_align, ref_start_idx)

    cons_region_list = []
    var_region_list = []

    for cons_region, index_tup in ref_cons_regions_index_d.items():
        start = index_tup[0]
        end = index_tup[1]
        if start in ref_numbering or end in ref_numbering:
            cons_region_list.append(cons_region)
        else:
            print("{} region not present".format(cons_region))

    for var_region, index_tup in ref_var_regions_index_d.items():
        start = index_tup[0]
        end = index_tup[1]
        if start in ref_numbering or end in ref_numbering:
            var_region_list.append(var_region)
        else:
            print("{} region not present".format(var_region))

    # find start/end region
    start_cons = sorted(cons_region_list)[0]
    end_cons = sorted(cons_region_list)[-1]
    start_var = sorted(var_region_list)[0]
    end_var = sorted(var_region_list)[-1]

    # get start region
    if full_order.index(start_cons) > full_order.index(start_var):
        start = start_var
    else:
        start = start_cons

    # get end region
    if full_order.index(end_cons) > full_order.index(end_var):
        end = end_cons
    else:
        end = end_var

    return start, end, ref_numbering


def get_var_regions_dict(ref_type, gene_region, regions_path):
    """
    imports regex search string from reference csv file
    :param ref_type: (str) name of the reference
    :param gene_region: (str) name of the gene
    :param regions_path: (str) path to the reference csv file
    :return: (dict) key = variable gene region name , value = sequence string
    """
    regions_file = os.path.join(regions_path, "HIV_var_cons_regions.csv")
    data = pd.read_csv(regions_file, sep=',', header=0, parse_dates=True, na_values=[' '])
    df = pd.DataFrame(data)
    print("reference is: ", ref_type)
    print("gene region is: ", gene_region)
    ref_df = df.loc[df["reference_type"] == ref_type]
    ref_gene_df = ref_df.loc[df["gene"] == gene_region]
    var_regions_dict = dict(zip(ref_gene_df["gene_region"], ref_gene_df["sequence"]))

    return var_regions_dict


def find_cons_var_regions(prot_sequence, regions_dict):

    regions_index_d = collections.defaultdict(int)
    regions_sequences_d = collections.defaultdict(str)
    for var_reg_name, var_seq in regions_dict.items():
        region_key = var_reg_name.split("_")[-1]
        var_reg = var_reg_name.split("_")[0]
        error = 3
        pattern = "({0}){{e<{1}}}".format(var_seq, error)
        match = regex.search(pattern, prot_sequence, regex.BESTMATCH)
        if match is not None:
            if region_key == "start":
                slice_index = match.start()
            elif region_key == "end":
                slice_index = match.end()
            else:
                sys.exit("error in region name: {}\nshould end in 'start' or 'end'.".format(var_reg_name))
            regions_index_d[var_reg] = slice_index
        else:
            print(var_reg_name, "not found")
            slice_index = None
            regions_index_d[var_reg] = slice_index

    for region, reg_idx in sorted(regions_index_d.items(), key=lambda x: x[0].split("_")[0]):
        print("do something")

    return regions_index_d


def write_regions_to_file(region_dict, path_for_tmp_file):
    """
    function to write dictionary to fasta file
    :param region_dict: (2D dict) dictionary = {seq_code: {region: sequence}}
    :param path_for_tmp_file: (str) path to where temp files will be created
    :return: (list) list of file names that have been written to file.
    """
    # create unique temp file name prefix
    temp_file_prefix = "tmp_"
    temp_file_prefix += ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))
    temp_file = os.path.join(path_for_tmp_file, temp_file_prefix)

    # write dicts to file
    file_names = []

    # catch error when only 1 sequence is present
    if len(region_dict.keys()) < 2:
        sys.exit("must have more than 1 sequence to align\nexiting")

    for seq_code, region_d in region_dict.items():
        for region, seq in region_d.items():
            # complete file name
            region_fn = "{0}_{1}.fasta".format(temp_file, region)

            if region_fn not in file_names:
                file_names.append(region_fn)

            with open(region_fn, 'a') as handle:
                handle.write(">{0}\n{1}\n".format(str(seq_code), seq))

    return file_names


def call_aligner(file_names):
    """
    Takes a dict of protein sequences, writes them to a temp file and aligns them with mafft.
    Aligned file is read back in and returned as a dictionary
    :param file_names: (list) list of the files to align
    :return: (dict) dictionary of aligned protein sequences: key = sequence, value = ID code
    """
    region_aligned_d = collections.defaultdict(dict)

    for file in file_names:
        region = os.path.split(file)[-1].split("_")[2].replace(".fasta", "")
        outfile = file.replace(".fasta", "_aligned.fasta")
        cmd = "mafft {0} > {1}".format(file, outfile)
        subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        aligned_region_d = fasta_to_dct_keep_gap(outfile)
        os.unlink(file)
        os.unlink(outfile)
        region_aligned_d[region].update(aligned_region_d)

    return region_aligned_d


def pad_var_region_to_longest(var_regions_dct):
    """
    pads length variable regions to all have the same length
    :param var_regions_dct: (dict) dict of the different variable regions for each sequence
    :return: (dict) gap padded dict of the different variable regions for each sequence
    """
    # initialise dicts
    new_var_regions_dct = collections.defaultdict(dict)
    var_region_lens = collections.defaultdict(list)
    max_lens_d = collections.defaultdict(int)

    # get len of each var region
    for seq_code, var_regions_d in var_regions_dct.items():
        for var_region, var_seq in var_regions_d.items():
            var_region_lens[var_region].append(len(var_seq))

    # get max len for each var region
    for var_reg, var_len_list in var_region_lens.items():
        max_lens_d[var_reg] = max(var_len_list)

    for seq_code, var_regions_d in var_regions_dct.items():
        for var_region, var_seq in var_regions_d.items():
            max_len = max_lens_d[var_region]
            this_seq_len = len(var_seq)
            gaps_to_pad = max_len - this_seq_len
            half_point = this_seq_len // 2
            new_var_seq = var_seq[:half_point] + "-" * gaps_to_pad + var_seq[half_point:]
            new_var_regions_dct[seq_code][var_region] = new_var_seq

    return new_var_regions_dct


def join_regions(cons_regions, padded_var_regions, full_order):
    """
    function to join conserved and variable regions to re-create the full sequence
    :param cons_regions: (dict) dictionary of the conserved regions {"C1": {"code": "seq"}}
    :param padded_var_regions: dictionary of the variable  regions {"V1": {"code": "seq"}}
    :param full_order: (list) a list of all the conserved and variable regions in sequential order
    :return: (dict) dictionary of re-created sequences {"code": "seq"}
    """
    joined_d = collections.defaultdict(str)

    # find start/end region
    cons = []
    var = []

    for i in full_order:
        if i in cons_regions:
            cons.append(i)

    for i in full_order:
        if i in padded_var_regions:
            var.append(i)

    start_cons = cons[0]
    end_cons = cons[-1]
    start_var = var[0]
    end_var = var[-1]

    # get start region
    if full_order.index(start_cons) > full_order.index(start_var):
        start = full_order.index(start_var)
    else:
        start = full_order.index(start_cons)

    # get end region
    if full_order.index(end_cons) > full_order.index(end_var):
        end = full_order.index(end_cons)
    else:
        end = full_order.index(end_var)

    final_order = full_order[start:end + 1]

    # stitch the regions together in the right order
    for seq_region in final_order:
        if seq_region in cons_regions:
            region_d = cons_regions[seq_region]
            for seq_code, seq in region_d.items():
                joined_d[seq_code] += seq
        elif seq_region in padded_var_regions:
            region_d = padded_var_regions[seq_region]
            for seq_code, seq in region_d.items():
                joined_d[seq_code] += seq
        else:
            sys.exit("region not in conserved or var dictionary\ncan't re-create full-length sequence")

    return joined_d


def backtranslate(padded_dna_d, prot_align_d):
    """
    function to backtranslate aligned protein sequence to aligned DNA sequence
    :param padded_dna_d: (dict) of gap padded dna sequences in frame 1, with indels padded with gaps {code: padded_seq}
    :param prot_align_d: (dict) of protein sequences {code: prot_seq}
    :return: (dict) of aligned dna sequences
    """
    dna_align_d = collections.defaultdict(str)
    for code, prot_seq in prot_align_d.items():
        dna_seq = padded_dna_d[code]
        dna_align = ''
        resi_count = 0
        for resi in prot_seq:
            if resi == '-':
                dna_align += '---'
            else:
                dna_align += dna_seq[(resi_count * 3):((resi_count * 3) + 3)]
                resi_count += 1
        dna_align_d[code] = dna_align

        if dna_align.replace("-", "") != dna_seq.replace("-", ""):
            print("input and output sequences are not identical")
            print(">{0}_in\n{1}\n".format(code, dna_seq))
            print(">{0}_out\n{1}\n".format(code, dna_align))

    return dna_align_d


def main(infile, outpath, name, ref, gene, bounds_lower, bounds_upper, var_align):

    # get absolute paths
    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)
    out_name = name + "_aligned.fasta"
    bad_name = name + "_NOT_aligned.fasta"
    outfile = os.path.join(outpath, out_name)
    badfile = os.path.join(outpath, bad_name)

    # get the reference DNA and prot sequence sequence
    gene_region = ref + "_" + gene
    get_script_path = os.path.realpath(__file__)
    script_folder = os.path.split(get_script_path)[0]
    script_folder = os.path.abspath(script_folder)
    ref_file = os.path.join(script_folder, "{}_seqs.fasta".format(ref))

    # get gene sequence
    if bounds_lower:
        bl = ((bounds_lower * 3) - 3)
    else:
        bl = None
    if bounds_upper:
        bu = ((bounds_upper * 3) + 3)
    else:
        bu = None

    ref_seqs = fasta_to_dct(ref_file)
    reference = ref_seqs[gene_region][bl:bu]

    ref_prot_seq = translate_dna(reference.replace("-", ""))

    # read in fasta file and reference
    in_seqs_d = fasta_to_dct_rev(infile)

    # generate seq_code to seq name list lookup dictionary
    first_look_up_d = collections.defaultdict(list)
    first_seq_code_d = collections.defaultdict(str)
    for i, (seq, names_list) in enumerate(in_seqs_d.items()):
        unique_id = str(i).zfill(4)
        first_look_up_d[unique_id] = names_list
        first_seq_code_d[seq] = unique_id

    # initialize dictionaries to collect cons and var regions and gap padded sequences
    cons_regions_dct = collections.defaultdict(dict)
    var_regions_dct = collections.defaultdict(dict)
    padded_seq_dict = collections.defaultdict(str)

    # set the cons and var regions
    ref_cons_regions_index_d, ref_var_regions_index_d, full_order = set_cons_var_regions(ref)
    regex_complied_1 = regex.compile(r"(^[-]*)", regex.V1)
    regex_complied_2 = regex.compile(r"([-]+)", regex.V1)
    # get the sequences for variable region boundaries for the ref-gene_region
    var_regions_dict = get_var_regions_dict(ref, gene, script_folder)
    var_reg_regex_compiled_d = collections.OrderedDict()
    # todo compile not working :(
    for var_region_name, var_seq in var_regions_dict.items():
        regex_complied = regex.compile("({0}){{e<3}}".format(var_seq))
        # print(regex_complied)
        var_reg_regex_compiled_d[var_region_name] = regex_complied

    # open file for sequences that could not be properly translated (hence codon aligned)
    with open(badfile, 'w') as handle:
        for seq, code in first_seq_code_d.items():
            # get pairwise alignment for query to reference
            seq_align, ref_align, frame = pairwise_align_dna(seq, reference, regex_complied_1)
            # correct for reading frame and indels
            padded_sequence = gap_padding(seq_align, ref_align, frame, regex_complied_2)
            padded_seq_dict[code] = padded_sequence

            # translate query
            prot_seq = translate_dna(padded_sequence)
            # if the seq could not be translated, write to file and skip
            if prot_seq.count("*") > 2:
                print("error in getting seq into frame", prot_seq)
                names_list = first_look_up_d[code]
                # del first_look_up_d[code]
                del padded_seq_dict[code]
                for name_bad in names_list:
                    handle.write(">{0}\n{1}\n".format(name_bad, seq))
                continue

            # get prot pairwise alignment
            # prot_align_d = prot_pairwise_align(prot_seq, ref_prot_seq, regex_complied_1)

            # new regex calls
            cons_regions_dct[code] = find_cons_var_regions(prot_seq, var_regions_dict)

            # get start and end regions
            # start_region, end_region, ref_numbering = find_start_end_regions(prot_align_d, ref_cons_regions_index_d,
            #                                                                  ref_var_regions_index_d, full_order,
            #                                                                  bounds_lower)

            # extract conserved regions
            # cons_regions_dct[code] = get_cons_regions(prot_align_d, end_region, ref_cons_regions_index_d, ref_numbering)

            # extract variable regions
            # var_regions_dct[code] = get_var_regions(prot_align_d, end_region, ref_var_regions_index_d, ref_numbering)

    # pad the variable regions with '-', to the longest sequence
    new_var_regions_dct = pad_var_region_to_longest(var_regions_dct)

    # write the collected conserved regions to file and align
    tmp_cons_file_to_align = write_regions_to_file(cons_regions_dct, outpath)
    align_cons_prot_d = call_aligner(tmp_cons_file_to_align)

    # write the collected variable regions to file and align (optional)
    if var_align:
        tmp_var_file_to_align = write_regions_to_file(new_var_regions_dct, outpath)
        var_prot_d = call_aligner(tmp_var_file_to_align)
    else:
        # reformat dict for joining of regions
        var_prot_d = collections.defaultdict(dict)
        for seq_code, region_d in new_var_regions_dct.items():
            for region, seq in region_d.items():
                var_prot_d[region][seq_code] = seq

    # join the different variable and conserved regions together in the right order
    joined_regions_d = join_regions(align_cons_prot_d, var_prot_d, full_order)

    # back-translate the protein alignment to a dna alignment
    dna_aligned = backtranslate(padded_seq_dict, joined_regions_d)

    # write back-translated DNA alignment to file
    with open(outfile, 'w') as handle:
        for code, align_seq in dna_aligned.items():
            names_list = first_look_up_d[code]
            for seq_name in names_list:
                handle.write(">{0}\n{1}\n".format(seq_name, align_seq))

    print("Codon aligning completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Codon aligns NGS HIV-1 ENV sequences using mafft.'
                                                 'Conserved regions are aligned separately and variable length/'
                                                 'difficult to align regions can be aligned or padded with gaps '
                                                 'to the longest subsequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path for the aligned output file', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='The name for the output file, ie: sample/participant or study name', required=True)
    parser.add_argument('-r', '--ref', default="CONSENSUS_C", choices=["CONSENSUS_C", "HXB2", "CON_OF_CONS"], type=str,
                        action='store', help='The choice of reference sequence. Either consensus of subtype C,'
                                             'HXB2 or consensus of consensus', required=False)
    parser.add_argument('-g', '--gene', default="ENV", choices=["ENV", "GAG", "POL", "PRO", "NEF", "VIF", "VPR",
                                                                "REV", "VPU"], type=str,
                        help='The name for the gene region (ENV, GAG, POL, PRO, NEF, VIF, VPR, REV, VPU)',
                        required=False)
    parser.add_argument('-bl', '--bounds_lower', default=None, type=int,
                        help='The start amino acid positions of your data, relative to the reference, '
                             'if this is not a full sequence',
                        required=False)
    parser.add_argument('-bu', '--bounds_upper', default=None, type=int,
                        help='The end amino acid positions of your data, relative to the reference, '
                             'if this is not a full sequence',
                        required=False)
    parser.add_argument('-v', '--var_align', default=False, action="store_true",
                        help='Align the variable regions as well. May produce messy alignment', required=False)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name
    ref = args.ref
    gene = args.gene
    bounds_upper = args.bounds_upper
    bounds_lower = args.bounds_lower
    var_align = args.var_align

    main(infile, outpath, name, ref, gene, bounds_lower, bounds_upper, var_align)
