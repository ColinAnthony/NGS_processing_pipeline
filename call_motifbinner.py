#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import argparse
import sys
from glob import glob
import subprocess
from subprocess import DEVNULL


__author__ = 'Colin Anthony'


def get_primer_lens_score(primer):
    count_Ns = primer.count("N")
    if count_Ns == 0:
        print("could not find primer ID in primer")
        sys.exit()
    total = len(primer)
    non_Ns = total - count_Ns - 1
    if non_Ns < 0:
        print("Error in primer {}".format(primer))
    primer_lens = '1' + str(count_Ns) + "," + str(non_Ns)
    primer_score = str(count_Ns + int((non_Ns * 0.8)))
    return primer_lens, primer_score


def run_motifbinner(logfile, inpath, read1, read2, outpath, fwd_primer, cDNA_primer, name_prefix):
    fwd_read = os.path.join(inpath, read1)
    fwd_primer_lens, fwd_primer_score = get_primer_lens_score(fwd_primer)
    rev_read = os.path.join(inpath, read2)
    cDNA_primer_lens, cDNA_primer_score  = get_primer_lens_score(cDNA_primer)
    fwd_pid = 'NULL'
    rev_pid_fragment = 2
    cmd = 'MotifBinner2.R --fwd_file={0} ' \
          '--fwd_primer_seq ={1} ' \
          '--fwd_primer_lens={2} ' \
          '--fwd_primer_min_score={3} ' \
          '--rev_file={4} ' \
          '--rev_primer_seq={5} ' \
          '--rev_primer_lens={6} ' \
          '--rev_primer_min_score={7} ' \
          '--fwd_pid_in_which_fragment={8} ' \
          '--rev_pid_in_which_fragment={} ' \
          '--output_dir={} ' \
          '--base_for_names={} ' \
          '--ncpu=3 ' \
          '--min_read_length=290 ' \
          '--overlapping'.format(fwd_read, fwd_primer, fwd_primer_lens, fwd_primer_score, rev_read, cDNA_primer,
                                 cDNA_primer_lens, cDNA_primer_score, fwd_pid, rev_pid_fragment, outpath, name_prefix)

    if os.path.exists(logfile):
        with open(logfile, 'a') as handle:
            handle.write("MotifBinner2 commands:\n{}".format(cmd))

    subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)


def main(logfile, inpath, outpath, fwd_primer, cDNA_primer):

    for file in glob(inpath + '*R1.fastq'):
        read1 = file
        read2 = file.replace("R1.fastq", "R2.fastq")
        name_prefix = file.replace("R1.fastq", "")

        run_motifbinner(logfile, inpath, read1, read2, outpath, fwd_primer, cDNA_primer, name_prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--inpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be created', required=True)
    parser.add_argument('-f', '--fwd_primer', default=argparse.SUPPRESS, type=str,
                        help='The fwd primer for these samples', required=True)
    parser.add_argument('-r', '--cDNA_primer', default=argparse.SUPPRESS, type=str,
                        help='The cDNA primer for these samples', required=True)
    parser.add_argument('-l', '--logfile', default=argparse.SUPPRESS, type=str,
                        help='The path and name of the log file', required=True)

    args = parser.parse_args()
    inpath = args.inpath
    outpath = args.outpath
    fwd_primer = args.fwd_primer
    cDNA_primer = args.cDNA_primer
    logfile = args.logfile

    main(inpath, outpath, fwd_primer, cDNA_primer, logfile)
