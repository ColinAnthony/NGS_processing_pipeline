#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import argparse
import subprocess
from subprocess import DEVNULL
from glob import glob

__author__ = 'Colin Anthony'


def main(path, name, script_folder, gene_region, fwd_primer, cDNA_primer, frame, stops, length, envelope):

    # initialize the log file
    main_path = os.path.join(path, name)
    logfile = os.path.join(main_path, "logfile.txt")
    with open(logfile, 'w') as handle:
        handle.write("Initializing log file \n")

    # call MotifBinner
    gene_region_folder = os.path.join(main_path, gene_region)
    inpath = os.path.join(gene_region_folder, '1raw')
    cons_outpath = os.path.join(gene_region_folder, '2consensus')
    motifbinner = os.path.join(script_folder, 'call_motifbinner.py')
    cmd2 = 'python3 {0} -p {1} -n {2} -g {3}'.format(motifbinner, logfile, inpath, cons_outpath, fwd_primer, cDNA_primer)
    subprocess.call(cmd2, shell=True)

    # copy data from nested binned folders into 2consensus folder
    path_to_consensus = '2consensus/*/n028_cons_seqLength/*kept_cons_seqLength.fastq/'
    desired_path = cons_outpath
    for cons_file in glob(path_to_consensus):
        oldpath, oldname = os.path.split(cons_file)
        new_name = os.path.join(desired_path, oldname)
        os.rename(cons_file, new_name)

    # convert to fasta
    fastq_path = os.path.join(gene_region_folder, '2cons')
    for fastq in glob(fastq_path + '*.fastq'):
        fasta = fastq.replace("fastq", "fasta")
        cmd3 = 'seqmagick convert {0} {1}'.format(fastq, fasta)
        subprocess.call(cmd3, shell=True)

    # remove the fastq files
    remove_fastq = os.path.join(fastq_path, '*.fastq')
    os.remove(remove_fastq)

    # call remove bad sequences
    remove_bad_seqs = os.path.join(script_folder, 'remove_bad_sequences.py')
    clean_path = os.path.join(gene_region_folder, '3cleaned')
    for fasta_file in glob(fastq_path + '*.fasta'):
        cmd4 = 'python3 {0} -i {1} -o {2} -f {3} -s {4} -l {5] -lf {6}'.format(remove_bad_seqs, fasta_file, clean_path,
                                                                              frame, stops, length, logfile)
        subprocess.call(cmd4, shell=True)

    # call cat all cleaned files into one file
    all_clean_path = os.path.join(gene_region_folder, '3cleaned')
    clean_name = name + "_" + gene_region + ".fasta"
    all_cleaned_outname = os.path.join(all_clean_path, clean_name)
    with open(all_cleaned_outname, 'w') as outfile:
        for fasta_file in glob(clean_path + '*clean.fa'):
            with open(fasta_file) as infile:
                for line in infile:
                    outfile.write(line)

    # move concatenated file to 4aligned
    aln_path = os.path.join(gene_region_folder, '4aligned')
    move_file = os.path.join(aln_path, clean_name)
    os.rename(all_cleaned_outname, move_file)

    # call align all samples
    if envelope:
        align_all = os.path.join(script_folder, 'align_all_env_samples.py')
        infile = move_file
        # infile, outpath, name, loop, v_loop_align, dna, aligner
        inpath, fname = os.path.split(infile)
        fname = name.replace(".fasta", "_aligned.fasta")
        cmd = 'python3 {0}  -i {1} -o {2} -n {3}'.format(align_all, infile, aln_path, fname)
        subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    else:
        align_all = os.path.join(script_folder, 'align_all_samples.py')
        infile = move_file
        inpath, fname = os.path.split(infile)
        fname = name.replace(".fasta", "_aligned.fasta")
        cmd = 'python3 {0}  -i {1} -o {2} -n {3}'.format(align_all, infile, aln_path, fname)
        subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Call make_folders script BEFORE running this script '
                                                 'Then copy your data into the /1raw/ folder'
                                                 'This script runs the NGS data processing pipeline. '
                                                 'Run this script using screen or nohup',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--path', default=argparse.SUPPRESS, type=str,
                        help='The path where you created the folders', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the name of the participant', required=True)
    parser.add_argument('-g', '--gene_region', default=argparse.SUPPRESS, type=str,
                        help='the genomic region being sequenced', required=True)
    parser.add_argument('-sf', '--script_folder', default=argparse.SUPPRESS, type=str,
                        help='the path to the folder containing the pipeline scripts', required=True)
    parser.add_argument('-f', '--fwd_primer', default=argparse.SUPPRESS, type=str,
                        help='The fwd primer for these samples', required=True)
    parser.add_argument('-r', '--cDNA_primer', default=argparse.SUPPRESS, type=str,
                        help='The cDNA primer for these samples', required=True)
    parser.add_argument('-fr', '--frame', type=int,
                        help='The reading frame (1, 2 or 3)', default=1, required=False)
    parser.add_argument('-s', '--stops', default=False, action='store_true',
                        help='Remove sequences with stop codons?)', required=False)
    parser.add_argument('-l', '--length', type=int,
                        help='The minimum read length)', required=False)
    parser.add_argument('-e', '--envelope', default=False, action='store_true',
                        help='are these sequences from HIV envelope?)', required=False)

    args = parser.parse_args()
    path = args.path
    name = args.name
    script_folder = args.script_folder
    gene_region = args.gene_region
    fwd_primer = args.fwd_primer
    cDNA_primer = args.cDNA_primer
    frame = args.frame
    stops = args.stops
    length = args.length
    envelope = args.envelope

    main(path, name, script_folder, gene_region, fwd_primer, cDNA_primer, frame, stops, length, envelope)
