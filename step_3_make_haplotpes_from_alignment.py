#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import subprocess
from glob import glob


__author__ = 'Colin Anthony'


def main(infile, script_folder):

    # split into sample files
    alignment_path = os.path.split(infile)[0]
    parent_path = os.path.split(alignment_path)[0]
    haplo_outpath = os.path.join(parent_path, '5haplotype')
    split_by_unique = os.path.join(script_folder, 'split_fasta_into_subfiles.py')
    cmd6 = 'python3 {0} -i {1} -o {2}'.format(split_by_unique, infile, haplo_outpath)

    subprocess.call(cmd6, shell=True)

    # haplotype
    split_fasta_files = os.path.join(haplo_outpath, "*_sep.fasta")
    haplotyper = os.path.join(script_folder, 'haplotyper_freq.py')
    for split_fasta in glob(split_fasta_files):
        cmd7 = 'python3 {0} -i {1} -o {2}'.format(haplotyper, split_fasta, haplo_outpath)
        subprocess.call(cmd7, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='takes an alignment file and splits it into time points '
                                                 '(or participants/ samples etc...) and collapses identical sequences '
                                                 'in each file, adding a frequency to the name',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', type=str,
                        help='The path and name of the aligned fasta file', required=True)
    parser.add_argument('-sf', '--script_folder', default=argparse.SUPPRESS, type=str,
                        help='the path to the folder containing the pipeline scripts', required=True)

    args = parser.parse_args()
    infile = args.infile
    script_folder = args.script_folder

    main(infile, script_folder)
