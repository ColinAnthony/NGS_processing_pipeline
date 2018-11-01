#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import subprocess
from glob import glob


__author__ = 'Colin Anthony'


def main(infile, field):
    get_script_path = os.path.realpath(__file__)
    script_folder = os.path.split(get_script_path)[0]
    script_folder = os.path.abspath(script_folder)

    # split into sample files
    infile = os.path.abspath(infile)
    alignment_path = os.path.split(infile)[0]
    parent_path = os.path.split(alignment_path)[0]
    haplo_outpath = os.path.join(parent_path, '5haplotype')
    if not os.path.isdir(haplo_outpath):
        print("could not find the 5haplotype folder\n check your file naming structure conforms to the required format")
        sys.exit()
    split_by_unique = os.path.join(script_folder, 'split_fasta_into_subfiles.py')
    cmd6 = 'python3 {0} -in {1} -o {2} -f {3}'.format(split_by_unique, infile, haplo_outpath, field)

    subprocess.call(cmd6, shell=True)

    # haplotype
    split_fasta_files = os.path.join(haplo_outpath, "*_sep.fasta")
    haplotyper = os.path.join(script_folder, 'haplotyper_freq.py')
    for split_fasta in glob(split_fasta_files):
        cmd7 = 'python3 {0} -in {1} -o {2}'.format(haplotyper, split_fasta, haplo_outpath)
        subprocess.call(cmd7, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='takes an alignment file and splits it into time points '
                                                 '(or participants/ samples etc...) and collapses identical sequences '
                                                 'in each file, adding a frequency to the name',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-in', '--infile', type=str,
                        help='The path and name of the aligned fasta file', required=True)
    parser.add_argument('-f', '--field', type=int, default=4, required=False,
                        help="The field that differentiates your samples/time points (use the last field if multiple."
                             "(ie: 4 for 'CAP177_2000_004wpi_V3C4_GGGACTCTAGTG_28, or 2 for SVB008_SP_GGTAGTCTAGTG_231")

    args = parser.parse_args()
    infile = args.infile
    field = args.field

    main(infile, field)
