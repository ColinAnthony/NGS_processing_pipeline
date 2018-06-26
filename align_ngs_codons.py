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
import pandas as pd
import regex
import seqanpy
from pprint import pprint


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


def get_order(sub_regions):
    if sub_regions:
        if sub_regions == "GP120" or sub_regions == "GP160":
            full_order = ["C1", "V1", "C2", "V2", "C3", "V3", "C4", "V4", "C5"]
        elif sub_regions == "C0C1" or sub_regions == "C2C3" or sub_regions == "GP41":
            full_order = ["C1"]
        elif sub_regions == "C1C2" or sub_regions == "C3C5":
            full_order = ["C1", "V1", "C2", "V2", "C3"]
        elif sub_regions == "P17" or sub_regions == "P24":
            full_order = ["C1"]
        else:
            sys.exit("Incorrect Env region")
    else:
        full_order = None

    return full_order


def pairwise_align_dna(sequence, reference, regex_complied, gene):
    """
    Pairwise align sequence to reference, to find reading frame and frame-shift in/dels
    :param sequence: (str) a query DNA sequence
    :param reference: (str) a reference DNA sequence (must start in reading frame 1)
    :param regex_complied: (regex_obj) a compiled regex pattern
    :param gene: (str) the target gene (ENV, GAG, POL, etc...
    :return: (str) aligned query sequence, (str) aligned ref sequence, (int) reading frame for query sequence
    """
    # do overlap pairwise alignment to not get truncated query sequence
    if gene == "ENV":
        overlap = seqanpy.align_overlap(sequence, reference, band=-1, score_match=4, score_mismatch=-1, score_gapext=-3,
                                        score_gapopen=-14)
    else:
        # for other regions
        overlap = seqanpy.align_overlap(sequence, reference, band=-1, score_match=4, score_mismatch=-2, score_gapext=-3,
                                        score_gapopen=-14)
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
    if ref_align[0] == '-':
        ref_start = regex_complied.search(ref_align).end()
    else:
        ref_start = 0

    # calculate reading frame (reference must start in frame 0)
    frame = (seq_start - ref_start) % 3

    # truncate the overlap alignment to the region of interest
    seq_align = seq_align[seq_start:seq_end]
    ref_align = ref_align[seq_start:seq_end]
    # print(">sqseq2\n{}\n".format(seq_align))
    # print(">sqref2\n{}\n".format(ref_align))

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
        if gap_start != 0:
            ins_pos = gap_start % 3

            if ins_pos == 1:
                ins_pos = gap_start + 3
            elif ins_pos == 2:
                ins_pos = gap_start + 2
            else:
                # ins_pos must == 0
                ins_pos = gap_start + 1
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
                indel_gap_fix_master.append((ins_pos, gap_len_in_seq, new_gap))

    # sort the list of all gaps to insert by gap start pos
    indel_gap_fix_master = sorted(indel_gap_fix_master)

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
    # stop codons coded as "Z" as mafft removed "*" characters
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
                  'TAC': 'Y', 'TAT': 'Y', 'TAA': 'Z', 'TAG': 'Z',
                  'TGC': 'C', 'TGT': 'C', 'TGA': 'Z', 'TGG': 'W',
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
    ref_df = df.loc[df["reference_type"] == ref_type]
    ref_gene_df = ref_df.loc[df["gene"] == gene_region]
    var_regions_dict = dict(zip(ref_gene_df["gene_region"], ref_gene_df["sequence"]))
    error = dict(zip(ref_gene_df["gene_region"], ref_gene_df["error"]))
    return var_regions_dict, error


def get_ref_start_end(ref_type, sub_regions, regions_path):
    """
    imports reference sequence start and end positions
    :param ref_type: (str) name of the reference
    :param sub_regions: (str) name of the env subregion
    :param regions_path: (str) path to the reference csv file
    :return: (dict) key = variable gene region name , value = sequence string
    """
    regions_file = os.path.join(regions_path, "gene_sub_regions_start_end.csv")
    data = pd.read_csv(regions_file, sep=',', header=0, parse_dates=True, na_values=[' '])
    df = pd.DataFrame(data)
    ref_df = df.loc[df["reference_type"] == ref_type]
    ref_region = ref_df.loc[df["sub_region"] == sub_regions]
    ref_start = int(ref_region["ref_start"]) * 3
    ref_end = int(ref_region["ref_end"]) * 3

    return ref_start, ref_end


def find_var_region_boundaries(prot_sequence, regions_dict, sub_regions, errors_allowed):
    """
    uses regex to find the start and end coordinates of HIV1 variable regions
    :param prot_sequence: (str) a protein sequence
    :param regions_dict: (dict) key = variable region name, value = regex string
    :param sub_regions: (str/None) name of env sub-region if present, else None
    :param errors_allowed: (int) number of miss-matches allowed in regex
    :return:
    """
    regions_index_d = collections.defaultdict(int)
    if sub_regions == "C0C1" or sub_regions == "C2C3" or sub_regions == "GP41" or sub_regions == "P17" \
            or sub_regions == "P24" or not sub_regions:
        regions_index_d["None"] = None
        return regions_index_d
    else:
        for var_reg_name, var_seq in regions_dict.items():
            if sub_regions == "C1C2" and var_reg_name.split("_")[0].upper() in ["V3", "V4"]:
                continue
            elif sub_regions == "C3C5" and var_reg_name.split("_")[0].upper() in ["V1", "V2"]:
                continue

            else:
                region_key = var_reg_name.split("_")[-1].lower()
                # set error in regex
                error = errors_allowed[var_reg_name]
                pattern = "{0}{{e<{1}}}".format(var_seq, error)

                match = regex.search(pattern, prot_sequence, regex.BESTMATCH)

                # if failed to get regex match for start of V1 for C1C2 amplicon data, try shorter regex search pattern
                if match is None and var_reg_name == "V1_start" and sub_regions == "C1C2":
                    alt_pattern = r'(XL[NKIE]C[NRTSI]){e<1}'
                    match = regex.search(alt_pattern, prot_sequence, regex.BESTMATCH)
                # if failed to get regex match for end of V2 for C1C2 amplicon data, try shorter regex search pattern
                if match is None and var_reg_name == "V2_end" and sub_regions == "C1C2":
                    alt_pattern = r'(Y[RKIV]L[IT][NRS]CN){e<2}'
                    match = regex.search(alt_pattern, prot_sequence, regex.BESTMATCH)
                    if match is None and var_reg_name == "V2_end" and sub_regions == "C1C2":
                        alt_pattern = r'(Y[RKIV]L[IT]X){e<1}'
                        match = regex.search(alt_pattern, prot_sequence, regex.BESTMATCH)
                # if failed to get regex match for start of V3 for C3C5 amplicon data, try shorter regex search pattern
                if match is None and var_reg_name == "V3_start" and sub_regions == "C3C5":
                    alt_pattern = r'(FYC[ND]T[ST].LF[NTSKD]){e<3}'
                    match = regex.search(alt_pattern, prot_sequence, regex.BESTMATCH)
                # if failed to get regex match for start of V4 for C3C5 amplicon data, try shorter regex search pattern
                if match is None and var_reg_name == "V4_start" and sub_regions == "C3C5":
                    alt_pattern = r'(LI[LV][TVL]RDGG.){e<3}'
                    match = regex.search(alt_pattern, prot_sequence, regex.BESTMATCH)
                    if match is None and var_reg_name == "V4_start" and sub_regions == "C3C5":
                        alt_pattern = r'(L[IL][LV][TVL]RD.){e<3}'
                        match = regex.search(alt_pattern, prot_sequence, regex.BESTMATCH)
                # if failed to get regex match for end of V4 for C3C5 amplicon data, try shorter regex search pattern
                if match is None and var_reg_name == "V4_end" and sub_regions == "C3C5":
                    alt_pattern = r'(E[TIV]FR){e<1}'
                    match = regex.search(alt_pattern, prot_sequence, regex.BESTMATCH)

                if match is not None:
                    if region_key == "start":
                        slice_index = match.end()
                    elif region_key == "end":
                        slice_index = match.start()
                    else:
                        sys.exit("error in region name: {}\nshould end in 'start' or 'end'.".format(var_reg_name))
                    regions_index_d[var_reg_name] = slice_index
                else:
                    slice_index = "missing"
                    regions_index_d[var_reg_name] = slice_index
                    print(var_reg_name, "not found")

        return regions_index_d


def check_for_missing_regex(var_region_index_dct):
    """
    return true if an enrty in the dict has value == False
    :param var_region_index_dct: (dict) key = sequence code, value = {key=region, value= index or bool}
    :return: bool
    """
    value = False
    for code, var_idx_d in var_region_index_dct.items():
        for region, idx in var_idx_d.items():
            if idx == "missing":
                value = True
                return value

    return value


def get_cons_regions(prot_sequence, regions_indx_dict, sub_regions):
    """
    extract conserved regions from protein sequence
    :param prot_sequence: (str) protein sequence
    :param regions_indx_dict: (dict) key = var region name, value = index in prot_sequence
    :param sub_regions: (str) the Env regions present in the sequence, if any
    :return: (dict) key = conserved region name, value = (str) slice of prot_sequence
    """

    cons_reg_sequences_d = collections.defaultdict(str)

    if sub_regions == "C0C1" or sub_regions == "C2C3" or sub_regions == "GP41" or sub_regions == "P17" \
            or sub_regions == "P24" or not sub_regions:
        cons_reg_sequences_d[sub_regions] = prot_sequence
        return cons_reg_sequences_d

    else:
        if sub_regions == "GP120" or sub_regions == "GP160":
            for code, reg_idx_d in sorted(regions_indx_dict.items(), key=lambda x: x[0].split("_")[0]):
                cons_reg_sequences_d["C1"] = prot_sequence[:reg_idx_d["V1_start"]]
                cons_reg_sequences_d["C2"] = prot_sequence[reg_idx_d["V1_end"]:reg_idx_d["V2_start"]]
                cons_reg_sequences_d["C3"] = prot_sequence[reg_idx_d["V2_end"]:reg_idx_d["V3_start"]]
                cons_reg_sequences_d["C4"] = prot_sequence[reg_idx_d["V3_end"]:reg_idx_d["V4_start"]]
                cons_reg_sequences_d["C5"] = prot_sequence[reg_idx_d["V4_end"]:]

        elif sub_regions == "C1C2":
            for code, reg_idx_d in sorted(regions_indx_dict.items(), key=lambda x: x[0].split("_")[0]):
                cons_reg_sequences_d["C1"] = prot_sequence[:reg_idx_d["V1_start"]]
                cons_reg_sequences_d["C2"] = prot_sequence[reg_idx_d["V1_end"]:reg_idx_d["V2_start"]]
                cons_reg_sequences_d["C3"] = prot_sequence[reg_idx_d["V2_end"]:]

        elif sub_regions == "C3C5":

            for code, reg_idx_d in sorted(regions_indx_dict.items(), key=lambda x: x[0].split("_")[0]):
                cons_reg_sequences_d["C1"] = prot_sequence[:reg_idx_d["V3_start"]]
                cons_reg_sequences_d["C2"] = prot_sequence[reg_idx_d["V3_end"]:reg_idx_d["V4_start"]]
                if reg_idx_d["V4_end"] == len(prot_sequence):
                    cons_reg_sequences_d["C3"] = prot_sequence[reg_idx_d["V4_end"] - 2:]
                else:
                    cons_reg_sequences_d["C3"] = prot_sequence[reg_idx_d["V4_end"]:]

        else:
            sys.exit("something went wrong")

        return cons_reg_sequences_d


def get_var_regions(prot_sequence, regions_indx_dict, sub_regions):
    """
    extract variable regions from protein sequence
    :param prot_sequence: (str) protein sequence
    :param regions_indx_dict: (dict) key = var region name, value = index in prot_sequence
    :param sub_regions: (str) the Env regions present in the sequence, if any
    :return: (dict) key = conserved region name, value = (str) slice of prot_sequence
    """
    var_reg_sequences_d = collections.defaultdict(str)

    if sub_regions == "C0C1" or sub_regions == "C2C3" or sub_regions == "GP41" or sub_regions == "P17" \
            or sub_regions == "P24" or not sub_regions:
        var_reg_sequences_d[sub_regions] = ""
        return var_reg_sequences_d

    else:
        if sub_regions == "GP120" or sub_regions == "GP160":
            for code, reg_idx_d in sorted(regions_indx_dict.items(), key=lambda x: x[0].split("_")[0]):
                for region, reg_idx in reg_idx_d.items():
                    # if missing the regex, send to bad output
                    if reg_idx is None:
                        return None

                var_reg_sequences_d["V1"] = prot_sequence[reg_idx_d["V1_start"]:reg_idx_d["V1_end"]]
                var_reg_sequences_d["V2"] = prot_sequence[reg_idx_d["V2_start"]:reg_idx_d["V2_end"]]
                var_reg_sequences_d["V3"] = prot_sequence[reg_idx_d["V3_start"]:reg_idx_d["V3_end"]]
                var_reg_sequences_d["V4"] = prot_sequence[reg_idx_d["V4_start"]:reg_idx_d["V4_end"]]

        elif sub_regions == "C1C2":
            for code, reg_idx_d in sorted(regions_indx_dict.items(), key=lambda x: x[0].split("_")[0]):
                for region, reg_idx in reg_idx_d.items():
                    # if missing the regex, send to bad output
                    if reg_idx is None:
                        return None

                var_reg_sequences_d["V1"] = prot_sequence[reg_idx_d["V1_start"]:reg_idx_d["V1_end"]]
                var_reg_sequences_d["V2"] = prot_sequence[reg_idx_d["V2_start"]:reg_idx_d["V2_end"]]

        elif sub_regions == "C3C5":
            for code, reg_idx_d in sorted(regions_indx_dict.items(), key=lambda x: x[0].split("_")[0]):
                for region, reg_idx in reg_idx_d.items():
                    # if missing the regex, send to bad output
                    if reg_idx is None:
                        return None

                var_reg_sequences_d["V1"] = prot_sequence[reg_idx_d["V3_start"]:reg_idx_d["V3_end"]]
                var_reg_sequences_d["V2"] = prot_sequence[reg_idx_d["V4_start"]:reg_idx_d["V4_end"]]
        else:
            sys.exit("something went wrong")

        return var_reg_sequences_d


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


def call_aligner(file_names, var):
    """
    Takes a dict of protein sequences, writes them to a temp file and aligns them with mafft.
    Aligned file is read back in and returned as a dictionary
    :param file_names: (list) list of the files to align
    :param var: (bool) True if the dict is for the variable regions, False if not
    :return: (dict) dictionary of aligned protein sequences: key = sequence, value = ID code
    """
    region_aligned_d = collections.defaultdict(dict)

    for file in file_names:
        region = os.path.split(file)[-1].split("_")[2].replace(".fasta", "")
        outfile = file.replace(".fasta", "_aligned.fasta")
        if var:
            cmd = "mafft --amino --op 1 --ep 0.1 {0} > {1}".format(file, outfile)
        else:
            cmd = "mafft --amino {0} > {1}".format(file, outfile)
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


def join_regions(cons_regions, padded_var_regions, full_order, sub_regions):
    """
    function to join conserved and variable regions to re-create the full sequence
    :param cons_regions: (dict) dictionary of the conserved regions {"C1": {"code": "seq"}}
    :param padded_var_regions: dictionary of the variable  regions {"V1": {"code": "seq"}}
    :param full_order: (list) a list of all the conserved and variable regions in sequential order
    :param sub_regions: (str) the Env regions present in the sequence, if any
    :return: (dict) dictionary of re-created sequences {"code": "seq"}
    """
    joined_d = collections.defaultdict(str)

    if sub_regions == "C0C1" or sub_regions == "C2C3" or sub_regions == "GP41" or sub_regions == "P17" \
            or sub_regions == "P24" or not sub_regions:
        for seq_region, region_d in cons_regions.items():
            for seq_code, seq in region_d.items():
                joined_d[seq_code] = seq

        return joined_d

    else:
        for seq_region in full_order:
            if seq_region in cons_regions:
                region_d = cons_regions[seq_region]
                for seq_code, seq in region_d.items():
                    joined_d[seq_code] += seq
            elif seq_region in padded_var_regions:
                region_d = padded_var_regions[seq_region]
                for seq_code, seq in region_d.items():
                    joined_d[seq_code] += seq

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
            print("Input and output sequences are not identical\nSequence removed from output")
            print(">{}\n{}\n".format(code + "_align_prot", prot_seq))
            print(">{}\n{}\n".format(code + "_DNA_input", dna_seq))
            print(">{0}_out\n{1}\n".format(code + "_DNA_back-translated", dna_align))
            print(">prot_align_{}".format(code))
            print("-------------")
            del dna_align_d[code]

    return dna_align_d


def main(infile, outpath, name, ref, gene, var_align, sub_region, user_ref):

    # get absolute paths
    infile = os.path.abspath(infile)
    outpath = os.path.abspath(outpath)
    out_name = name + "_aligned.fasta"
    bad_name = name + "_NOT_aligned.fasta"
    outfile = os.path.join(outpath, out_name)
    badfile = os.path.join(outpath, bad_name)

    print("Starting codon alignment for file {}".format(infile))
    print("-------------------------------------")
    print("Output aligned file will be {}".format(outfile))
    print("Reference sequence is {}".format(ref))
    print("-------------------------------------")
    print("Gene region is {}".format(gene))
    print("Gene sub-region (if any) is {}".format(sub_region))
    print("-------------------------------------")
    # get the reference DNA sequence
    get_script_path = os.path.realpath(__file__)
    script_folder = os.path.split(get_script_path)[0]
    script_folder = os.path.abspath(script_folder)
    if not user_ref:
        gene_region = ref + "_" + gene
        ref_file = os.path.join(script_folder, "reference_sequences.fasta")
        ref_seqs = fasta_to_dct(ref_file)
        reference = ref_seqs[gene_region]
        hxb2_name = "HXB2_{}".format(gene)
        ref_hxb = ref_seqs[hxb2_name]
    else:
        ref_seqs = fasta_to_dct(user_ref)
        reference = ref_seqs[list(ref_seqs)[0]]

        ref_file = os.path.join(script_folder, "reference_sequences.fasta")
        ref_seqs = fasta_to_dct(ref_file)
        hxb2_name = "HXB2_{}".format(gene)
        ref_hxb = ref_seqs[hxb2_name]
        print("custom ref", reference)

    # read in fasta file and reference
    in_seqs_d = fasta_to_dct_rev(infile)

    # generate seq_code to seq name list lookup dictionary
    first_look_up_d = collections.defaultdict(list)
    first_seq_code_d = collections.defaultdict(str)

    # get internal reference
    longest_seq = ''
    seq_length = 200
    for i, (seq, names_list) in enumerate(in_seqs_d.items()):
        unique_id = str(i).zfill(4)
        first_look_up_d[unique_id] = names_list
        first_seq_code_d[seq] = unique_id
        seq_len = len(seq)
        seq_abundance = len(names_list)
        if seq_len > seq_length and seq_abundance > 3:
            seq_length = seq_len
            longest_seq = seq

    # add hxb2 to the alignment
    first_look_up_d[hxb2_name] = [hxb2_name]
    if sub_region:
        hxb2_start, hxb2_end = get_ref_start_end("HXB2", sub_region, script_folder)
        ref_hxb = ref_hxb[hxb2_start:hxb2_end]
        first_seq_code_d[ref_hxb] = hxb2_name
    else:
        first_seq_code_d[ref_hxb] = hxb2_name

    if not user_ref:
        # get reading frame of most abundant internal reference
        internal_ref_frame_1 = longest_seq
        internal_ref_frame_2 = "N" + longest_seq
        internal_ref_frame_3 = "NN" + longest_seq
        frame_1_tr = translate_dna(internal_ref_frame_1)
        frame_2_tr = translate_dna(internal_ref_frame_2)
        frame_3_tr = translate_dna(internal_ref_frame_3)
        frame_1_stops = frame_1_tr[:-1].count("Z")
        frame_2_stops = frame_2_tr[:-1].count("Z")
        frame_3_stops = frame_3_tr[:-1].count("Z")
        if frame_1_stops < 1:
            internal_reference = internal_ref_frame_1
        elif frame_2_stops < 1:
            internal_reference = internal_ref_frame_2
        elif frame_3_stops < 1:
            internal_reference = internal_ref_frame_3
        else:
            internal_reference = None

        if internal_reference is not None:
            reference = internal_reference
            user_ref = True

    # initialize dictionaries to collect cons and var regions and gap padded sequences
    cons_regions_dct = collections.defaultdict(dict)
    var_regions_dct = collections.defaultdict(dict)
    padded_seq_dict = collections.defaultdict(str)

    # set the cons and var regions
    full_order = get_order(sub_region)

    regex_complied_1 = regex.compile(r"(^[-]*)", regex.V1)
    regex_complied_2 = regex.compile(r"([-]+)", regex.V1)

    # get the sequences for variable region boundaries for the ref-gene_region
    if gene == "ENV":
        var_region_regex_dct, errors_allow = get_var_regions_dict(ref, gene, script_folder)
    else:
        var_region_regex_dct = {}
        errors_allow = None

    if not user_ref and sub_region:
        ref_start, ref_end = get_ref_start_end(ref, sub_region, script_folder)
        reference = reference[ref_start:ref_end]

    # compile the regex strings
    # var_reg_regex_compiled_d = collections.OrderedDict()
    # todo compile not working :(
    # for var_region_name, var_seq in var_region_regex_dct.items():
    #     regex_complied = regex.compile("({0}){{e<3}}".format(var_seq))
    #     # print(regex_complied)
    #     var_reg_regex_compiled_d[var_region_name] = regex_complied

    # open file for sequences that could not be properly translated (hence codon aligned)
    print("Processing sequences\n")
    bad_seq_counter = 0
    with open(badfile, 'w') as handle:
        for seq, code in first_seq_code_d.items():
            seq = seq.replace("-", "")
            var_region_index_dct = collections.defaultdict(dict)
            # get pairwise alignment for query to reference
            seq_align, ref_align, frame = pairwise_align_dna(seq, reference, regex_complied_1, gene)

            # pad indels with gaps
            padded_sequence = gap_padding(seq_align, ref_align, frame, regex_complied_2)

            # translate query
            prot_seq = translate_dna(padded_sequence)

            # if the seq could not be translated, write to file and skip
            if prot_seq[:-1].count("Z") > 1:
                print("error getting seq into frame", prot_seq)
                names_list = first_look_up_d[code]
                for name_bad in names_list:
                    bad_seq_counter += 1
                    handle.write(">{0}\n{1}\n".format(name_bad, seq))

                continue
            else:
                padded_seq_dict[code] = padded_sequence

            # get the var region boundaries, if any
            var_region_index_dct[code] = find_var_region_boundaries(prot_seq, var_region_regex_dct, sub_region,
                                                                    errors_allow)
            # if one or more of the var region boundaries was not found, write to file and skip
            missing_regex = check_for_missing_regex(var_region_index_dct)
            if missing_regex:
                pprint(var_region_index_dct[code])
                print("error finding one or more variable region boundary", prot_seq)
                names_list = first_look_up_d[code]
                del padded_seq_dict[code]
                for name_bad in names_list:
                    bad_seq_counter += 1
                    handle.write(">{0}\n{1}\n".format(name_bad, seq))
                continue

            # extract conserved regions
            cons_regions_dct[code] = get_cons_regions(prot_seq, var_region_index_dct, sub_region)

            # extract variable regions
            var_regions_dct[code] = get_var_regions(prot_seq, var_region_index_dct, sub_region)

    # write the collected conserved regions to file and align
    print("Aligning conserved regions sequences\n")
    tmp_cons_file_to_align = write_regions_to_file(cons_regions_dct, outpath)
    var = False
    align_cons_prot_d = call_aligner(tmp_cons_file_to_align, var)

    # write the collected variable regions to file and align (optional)
    if var_align:
        print("Aligning variable region sequences\n")
        tmp_var_file_to_align = write_regions_to_file(var_regions_dct, outpath)
        var = True
        var_prot_d = call_aligner(tmp_var_file_to_align, var)

    # pad the variable regions with '-', to the longest sequence
    else:
        print("Padding variable region sequences with gaps\n")
        new_var_regions_dct = pad_var_region_to_longest(var_regions_dct)

        # reformat dict for joining of regions
        var_prot_d = collections.defaultdict(dict)
        for seq_code, region_d in new_var_regions_dct.items():
            for region, seq in region_d.items():
                var_prot_d[region][seq_code] = seq

    # join the different variable and conserved regions together in the right order
    print("Joining conserved and variable regions\n")
    joined_regions_d = join_regions(align_cons_prot_d, var_prot_d, full_order, sub_region)

    # back-translate the protein alignment to a dna alignment
    print("Back-translating from protein to DNA alignment\n")
    dna_aligned = backtranslate(padded_seq_dict, joined_regions_d)

    # write back-translated DNA alignment to file
    print("Writing DNA alignment to outfile\n")
    with open(outfile, 'w') as handle:
        for code, align_seq in dna_aligned.items():
            names_list = first_look_up_d[code]
            for seq_name in names_list:
                handle.write(">{0}\n{1}\n".format(seq_name, align_seq))

    print("Total sequences not aligned: ", bad_seq_counter)
    print("Codon aligning completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Codon aligns NGS HIV-1 ENV sequences using mafft.'
                                                 'Conserved regions are aligned separately and variable length/'
                                                 'difficult to align sub-regions can be aligned or padded with gaps '
                                                 'to the longest sub-sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', default=argparse.SUPPRESS, type=str,
                        help='The input file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path for the aligned output file', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='The name for the output file, ie: sample/participant or study name', required=True)
    parser.add_argument('-r', '--ref', default="CONSENSUS_C", action='store',
                        choices=["CONSENSUS_C", "CONSENSUS_A1", "CONSENSUS_B", "HXB2", "CON_OF_CONS"], type=str,
                        help='The choice of reference sequence. Either consensus of subtype, '
                             'HXB2 or consensus of consensus', required=False)
    parser.add_argument('-g', '--gene', default="ENV", action="store",
                        choices=["ENV", "GAG", "POL", "PRO", "NEF", "VIF", "VPR", "REV", "VPU"], type=str,
                        help='The name for the gene region (ENV, GAG, POL, PRO, NEF, VIF, VPR, REV, VPU)',
                        required=False)
    parser.add_argument('-reg', '--regions', default=False, action="store",
                        choices=["C0C1", "C1C2", "C2C3", "C3C5", "GP41", "GP120", "GP160", "P17", "P24"], type=str,
                        help='the variable regions in your data', required=False)
    parser.add_argument('-v', '--var_align', default=False, action="store_true",
                        help='Align the variable regions as well. May produce messy alignment', required=False)
    parser.add_argument('-u', '--user_ref', default=False, type=str,
                        help='the path and file name for the custom DNA reference sequence, '
                             'must start in reading frame 1', required=False)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name
    ref = args.ref
    gene = args.gene
    var_align = args.var_align
    regions = args.regions
    user_ref = args.user_ref

    if gene == "ENV":
        if not regions:
            sys.exit("must use the -reg flag for ENV")

    main(infile, outpath, name, ref, gene, var_align, regions, user_ref)
