from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import collections
from glob import glob
from itertools import groupby
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
        new_key = k.replace(" ", "_")
        if new_key in dct.keys():
            print("Duplicate sequence ids found. Exiting")
            raise KeyError("Duplicate sequence ids found")
        dct[new_key] = str(v).replace("~", "_")

    return dct


def main(inpath, outfile):

    print("Calculating sequencing depth and yield statistics")
    # initialize master dict to return
    stats_d = collections.defaultdict(list)

    all_names = collections.defaultdict(str)

    # clear the stats_out file
    with open(outfile, 'w') as handle:
        handle.write("")

    # initialise the dictionary with headers for the outfile
    stats_d["headers"].append("name")
    stats_d["headers"].append("raw_sequences")
    stats_d["headers"].append("merged_sequences")
    stats_d["headers"].append("consensus_sequences")
    stats_d["headers"].append("cleaned_sequences")
    stats_d["headers"].append("no_contaminants")
    stats_d["headers"].append("prop_detection_1_percent_variant")
    stats_d["headers"].append("variant_freq_with_95_percent_detection")

    binned_folders = os.path.join(inpath, "1consensus", "binned", "*")
    cleaned_files = os.path.join(inpath, "2cleaned", "*_clean.fasta")
    contam_files = os.path.join(inpath, "3contam_removal", "*_good.fasta")

    bin_list = ['raw', 'merged', 'consensus']
    for binned_folder in glob(binned_folders):
        name = os.path.split(binned_folder)[-1]
        all_names[name] = "True"
        raw = os.path.join(binned_folder, "n001_fwd_loadData", "n001_fwd_loadData.csv")
        merged = os.path.join(binned_folder, "n019_mergePEAR", "n019_mergePEAR.csv")
        consensus = os.path.join(binned_folder, "n023_buildConsensus", "n023_buildConsensus.csv")
        stats_d[name].append(name)
        for index, item in enumerate([raw, merged, consensus]):
            data = pd.read_csv(item, sep=',', header=0, parse_dates=True, na_values=[' '])
            df = pd.DataFrame(data)
            headers = list(df)
            count = list(df[headers[2]])[0]

            stats_d[name].append(count)

    # calculate number of cleaned sequences
    for cleaned_file in glob(cleaned_files):
        name = os.path.split(cleaned_file)[-1].replace("_clean.fasta", "")
        if name not in all_names.keys():
            print("Can't match name for cleaned file with parent file name")
            print("name", name)
            print("not in", all_names.keys())
            sys.exit()

        clean_d = fasta_to_dct(cleaned_file)
        total_clean = str(len(clean_d.keys()))
        stats_d[name].append(total_clean)

    # calculate number of sequences after contam removal
    for contam_file in glob(contam_files):
        name = os.path.split(contam_file)[-1].replace("_good.fasta", "")
        if name not in all_names.keys():
            print("Can't match name for no_contam file with parent file name")
            print("name", name)
            print("not in", all_names.keys())
            sys.exit()

        contam_rem_d = fasta_to_dct(contam_file)
        total_contam_rem = str(len(contam_rem_d.keys()))
        stats_d[name].append(total_contam_rem)

        # specify target frequency to detect
        freq = 1

        # calculate probability of detection
        num_consensus_seqs = int(total_contam_rem)
        p = freq/100
        if num_consensus_seqs != 0:
            var_detection_limit = round((1 - ((1 - p) ** num_consensus_seqs))*100, 2)
            stats_d[name].append(str(var_detection_limit))
            var_freq_with_95perc_prob = round((1 - (0.05 ** (1 / num_consensus_seqs))) * 100, 3)
            stats_d[name].append(str(var_freq_with_95perc_prob))
        else:
            stats_d[name].append('NaN')
            stats_d[name].append("NaN")

    # write the stats to the log file
    with open(outfile, 'w') as handle:
        # write the headers
        headers_to_write = ",".join(stats_d["headers"])
        handle.write(headers_to_write + "\n")
        # write the stats
        for sequence_file in stats_d.keys():
            if sequence_file != "headers":
                stats_list = stats_d[sequence_file]
                # if sample failed at some point, add NaN's
                if len(stats_list) < len(stats_d["headers"]):
                    nan_to_add = 6 - len(stats_list)
                    for i in range(nan_to_add):
                        stats_list.append("NaN")
                lines_to_write = ",".join(str(x) for x in stats_list)
                handle.write(lines_to_write + "\n")

    print("Stats calculations on your NGS samples are complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='calculate stats on sample processing',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in', '--inpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-o', '--outfile', default=argparse.SUPPRESS, type=str,
                        help='The path and name for the stats output file', required=True)

    args = parser.parse_args()
    inpath = args.inpath
    outfile = args.outfile

    main(inpath, outfile)
