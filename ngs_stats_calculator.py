#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import collections
from glob import glob
from Bio import SeqIO

__author__ = 'Colin Anthony'


def fasta_to_dct(fn):
    '''
    converts a fasta file to a dictionary where key = seq name, value = sequence
    :param fn: a fasta file
    :return: a dictionary
    '''
    dct = collections.OrderedDict()
    for seq_record in SeqIO.parse(open(fn), "fasta"):
        dct[seq_record.description.replace(" ", "_")] = str(seq_record.seq).replace("~", "-").upper()
    return dct


def fastq_to_dct(fn):
    '''
    converts a fastq file to a dictionary where key = seq name, value = sequence
    :param fn: a fastq file
    :return: a dictionary
    '''
    dct = collections.OrderedDict()
    for seq_record in SeqIO.parse(open(fn), "fastq"):
        dct[seq_record.description.replace(" ", "_")] = str(seq_record.seq).replace("~", "-").upper()
    return dct


def main(inpath, outfile):

    print("Calculating sequencing depth and yield statistics")
    # initialize master dict to return
    stats_d = collections.defaultdict(list)

    all_names = collections.defaultdict(str)

    # calculate number of raw sequences
    stats_d["headers"].append("name")
    stats_d["headers"].append("raw_sequences")
    raw_files = os.path.join(inpath, "0raw", "*_R1.fastq")
    for raw_file in glob(raw_files):
        name = os.path.split(raw_file)[-1].replace("_R1.fastq", "")
        all_names[name] = name
        raw_d = fastq_to_dct(raw_file)
        total_raw = str(len(raw_d.keys()))
        stats_d[name].append(name)
        stats_d[name].append(total_raw)

    # calculate number of sequences after contam removal
    stats_d["headers"].append("contam_removed")
    contam_files = os.path.join(inpath, "1contam_removal", "*_cln_R1.fastq")
    for contam_file in glob(contam_files):
        name = os.path.split(contam_file)[-1].replace("_cln_R1.fastq", "")
        print("contram\n", contam_file)
        print(name)
        if name not in all_names.keys():
            print("Can't match name for merged file with parent file name")
            print(name)
            sys.exit()

        contam_rem_d = fastq_to_dct(contam_file)
        total_contam_rem = str(len(contam_rem_d.keys()))
        stats_d[name].append(total_contam_rem)

    # calculate number of merged sequences
    merged_files = os.path.join(inpath, "2consensus", "binned", "*", "*_mergePEAR", "*assembled.fastq")
    stats_d["headers"].append("merged_sequences")
    for merged_file in glob(merged_files):
        # get rid of the generic file name "merged.fastq.assembled.fastq"
        path_split1 = os.path.split(merged_file)[0]
        # get rid of the generic folder name "*_mergePEAR"
        path_split2 = os.path.split(path_split1)[0]
    # get the sample name
        name = os.path.split(path_split2)[-1]
        # Check you have the correct sample name
        if name not in all_names.keys():
            print("Can't match name for merged file with parent file name")
            print(name)
            sys.exit()

        merged_d = fastq_to_dct(merged_file)
        total_merged = str(len(merged_d.keys()))
        stats_d[name].append(total_merged)

    # calculate number of consensus sequences
    stats_d["headers"].append("consensus_sequences")
    consensus_files = os.path.join(inpath, "2consensus", "binned", "*", "*_buildConsensus", "*kept_buildConsensus.fastq")
    for consensus_file in glob(consensus_files):
        name = os.path.split(consensus_file)[-1].replace("_kept_buildConsensus.fastq", "")
        if name not in all_names.keys():
            print("Can't match name for consensus file with parent file name")
            print(name)
            sys.exit()
        consensus_d = fastq_to_dct(consensus_file)
        total_consensus = str(len(consensus_d.keys()))
        stats_d[name].append(total_consensus)

    # calculate number of cleaned sequences
    stats_d["headers"].append("cleaned_sequences")
    cleaned_files = os.path.join(inpath, "3cleaned", "*_clean.fasta")
    for cleaned_file in glob(cleaned_files):
        name = os.path.split(cleaned_file)[-1].replace("_clean.fasta", "")
        if name not in all_names.keys():
            print("Can't match name for cleaned file with parent file name")
            print(name)
            sys.exit()
        clean_d = fasta_to_dct(cleaned_file)
        total_clean = str(len(clean_d.keys()))
        stats_d[name].append(total_clean)
    # write the stats to the log file
    with open(outfile, 'w') as handle:
        # write the headers
        headers_to_write = ",".join(stats_d["headers"])
        handle.write(headers_to_write + "\n")
        del stats_d["headers"]

        # write the stats
        for sequence_file in stats_d.keys():
            lines_to_write = ",".join(stats_d[sequence_file])
            handle.write(lines_to_write + "\n")

    print("Stats calculations on your NGS samples are complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='calculate stats on sample processing',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--inpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-o', '--outfile', default=argparse.SUPPRESS, type=str,
                        help='The path and name for the stats output file', required=True)

    args = parser.parse_args()
    inpath = args.inpath
    outfile = args.outfile

    main(inpath, outfile)
