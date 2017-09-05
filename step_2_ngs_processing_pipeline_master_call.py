#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import sys
from shutil import copyfile
import argparse
import subprocess
from glob import glob
import regex
from Bio import SeqIO


__author__ = 'Colin Anthony'
#

def fasta_to_dct(fn):
    '''
    :param fn: (str)infile name
    :return: (dict) dictionary of names and sequences
    '''
    dct = {}
    for seq_record in SeqIO.parse(open(fn), "fasta"):
        dct[seq_record.description.replace(" ", "_").upper()] = str(seq_record.seq).replace("-", "").upper()
    return dct


def rename_sequences(raw_files_search):

    print("renaming raw files")
    for inf_R1 in raw_files_search:
        inf_R2 = inf_R1.replace("R1_001.fastq", "R2_001.fastq")
        outf_R1 = inf_R1.replace("-", "_")
        outf_R2 = inf_R2.replace("-", "_")

        outf_R1_rename = regex.sub("S[0-9][0-9]_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq", "R1.fastq", outf_R1)
        print(outf_R1)
        if outf_R1_rename == outf_R1:
            print("Unable to rename R1 file {0}\ncheck the file renaming regex".format(inf_R1))
            # check if they have already been renamed
            if outf_R1.split("_")[-1] == "R1.fastq":
                return
                # raise Exception("Already renamed?")
            else:
                sys.exit()
        else:
            os.rename(inf_R1, outf_R1_rename)

        outf_R2_rename = regex.sub("S[0-9][0-9]_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9].fastq", "R2.fastq", outf_R2)
        if outf_R2_rename == outf_R2:
            print("Unable to rename R2 file {0}\ncheck the file renaming regex".format(inf_R2))
            if outf_R2.split("_")[-1] == "R2.fastq":
                return
                # raise Exception("Already renamed?")
            else:
                sys.exit()
        else:
            os.rename(inf_R2, outf_R2_rename)


def call_motifbinner(raw_files, motifbinner, cons_outpath, counter, logfile):
    for file in raw_files:
        read1 = file
        read2 = file.replace("R1.fastq", "R2.fastq")
        name_prefix = os.path.split(file)[-1].replace("_R1.fastq", "")

        cmd1 = 'python3 {0} -r1 {1} -r2 {2} -o {3} -f {4} -r {5} -n {6} -c {7} -l {8}'.format(motifbinner,
                                                                                              read1,
                                                                                              read2,
                                                                                              cons_outpath,
                                                                                              fwd_primer,
                                                                                              cDNA_primer,
                                                                                              name_prefix,
                                                                                              counter,
                                                                                              logfile)

        subprocess.call(cmd1, shell=True)
        counter += 1


def delete_gaps(fasta_infiles):

    for fasta_file in fasta_infiles:
        temp_out = fasta_file.replace(".fasta", ".fasta.bak")
        new_fasta = temp_out.replace(".bak", "")
        cons_d = fasta_to_dct(fasta_file)
        with open(temp_out, 'w') as handle:
            for seq_name, seq in cons_d.items():

                handle.write('>{0}\n{1}\n'.format(seq_name, seq))

        os.remove(fasta_file)
        copyfile(temp_out, new_fasta)
        os.remove(temp_out)


def call_fasta_cleanup(contam_removed_fasta, remove_bad_seqs, clean_path, logfile, stops):

    for fasta_file in contam_removed_fasta:
        if stops:
            cmd4 = 'python3 {0} -i {1} -o {2} -f {3} -s {4} -l {5} -lf {6}'.format(remove_bad_seqs,
                                                                                   fasta_file,
                                                                                   clean_path,
                                                                                   frame,
                                                                                   stops,
                                                                                   length,
                                                                                   logfile)
        else:
            cmd4 = 'python3 {0} -i {1} -o {2} -f {3} -l {4} -lf {5}'.format(remove_bad_seqs,
                                                                            fasta_file,
                                                                            clean_path,
                                                                            frame,
                                                                            length,
                                                                            logfile)
        if os.path.exists(logfile):
            with open(logfile, 'a') as handle:
                handle.write("\nremove_bad_sequences commands:\n{}\n".format(cmd4))

        subprocess.call(cmd4, shell=True)


def call_contam_check(consensuses, contam_removal_script, contam_removed_path, gene_region, logfile):

    for consensus_file in consensuses:
        cmd3 = 'python3 {0} -i {1} -o {2} -g {3} -l {4}'.format(contam_removal_script,
                                                         consensus_file,
                                                         contam_removed_path,
                                                         gene_region,
                                                         logfile)

        subprocess.call(cmd3, shell=True)


def call_align(envelope, script_folder, to_align, aln_path, fname):

    if envelope is not None:
        align_all = os.path.join(script_folder, 'align_all_env_samples.py')

        cmd5 = 'python3 {0}  -i {1} -o {2} -l {3}'.format(align_all, to_align, aln_path, fname, envelope)
    else:
        align_all = os.path.join(script_folder, 'align_all_samples.py')
        cmd5 = 'python3 {0}  -i {1} -o {2} -n {3}'.format(align_all, to_align, aln_path, fname)

    subprocess.call(cmd5, shell=True)


def main(path, name, script_folder, gene_region, fwd_primer, cDNA_primer, frame, stops, length, envelope):

    # initialize the log file
    logfile = os.path.join(path, (gene_region + "_logfile.txt"))
    with open(logfile, 'w') as handle:
        handle.write("Log File,{0}_{1}\n".format(name, gene_region))

    # rename the raw sequences
    raw_fastq_inpath = os.path.join(path, '0raw')
    raw_files_search = os.path.join(raw_fastq_inpath, "*R1*.fastq")
    raw_files = glob(raw_files_search)
    if not raw_files:
        print("No raw files were found\n"
              "Check that files end with R1.fastq and R2.fastq")
        sys.exit()

    rename_sequences(raw_files)

    # run the call_MotifBinner script which will loop over fastq files in the target folder
    motifbinner = os.path.join(script_folder, 'call_motifbinner.py')
    rename_in = raw_files_search = os.path.join(raw_fastq_inpath, "*R1*.fastq")
    cons_outpath = os.path.join(path, '1consensus', 'binned')
    counter = 0
    call_motifbinner(rename_in, motifbinner, cons_outpath, counter, logfile)
    
    # check if the consensus files exist
    path_to_nested_consensuses = os.path.join(path, '1consensus/binned/*/*_buildConsensus/*_kept_buildConsensus.fastq')
    nested_consesnsuses = glob(path_to_nested_consensuses)
    if not nested_consesnsuses:
        print("No consensus sequences were found\n"
              "This is likely if MotifBinner was not able to complete\n"
              "Do your fastq sequences end in R1.fastq/R2.fastq?\n"
              "Check the primer sequences are correct\nand check the binning report in the 1consensus/binned/ folder")
        sys.exit()

    # copy data from nested binned folders into 1consensus folder
    print("Copy fastq files from nested folders to '1consensus' folder")
    consensus_path = os.path.join(path, '1consensus')
    for cons_file in nested_consesnsuses:
        old_path, old_name = os.path.split(cons_file)
        new_name = old_name.replace("_kept_buildConsensus", "")
        new_file = os.path.join(consensus_path, new_name)
        copyfile(cons_file, new_file)

    # convert copied fastq to fasta
    print("Converting fastq to fasta")
    cons_search_path = os.path.join(consensus_path, '*.fastq')
    consensuses = glob(cons_search_path)
    for fastq in consensuses:
        fasta = fastq.replace("fastq", "fasta")
        cmd2 = 'seqmagick convert {0} {1}'.format(fastq,
                                                  fasta)

        subprocess.call(cmd2, shell=True)

    # remove the copied fastq files
    print("Removing the copied fastq files")
    remove_fastq = cons_search_path
    for old_fastq_copy in glob(remove_fastq):
        os.remove(old_fastq_copy)

    # delete any gaps characters in the fasta sequences
    print("deleting gaps in consensus sequences")
    consensus_infiles = os.path.join(consensus_path, '*.fasta')
    consensus_infiles = glob(consensus_infiles)
    delete_gaps(consensus_infiles)

    # call remove bad sequences
    print("Removing 'bad' sequences")
    remove_bad_seqs = os.path.join(script_folder, 'remove_bad_sequences.py')
    clean_path = os.path.join(path, '2cleaned')
    if not consensus_infiles:
        print("Could not find consensus fasta files\n"
              "It is possible that something went wrong when copying the consensus sequences from the nested folders, "
              "to the 1consensus folder")
        sys.exit()

    call_fasta_cleanup(consensus_infiles, remove_bad_seqs, clean_path, logfile, stops)

    # remove contaminating sequences
    print("removing contaminating non-HIV sequences")
    contam_removal_script = os.path.join(script_folder, "contam_removal.py")
    clean_search = os.path.join(clean_path, "*clean.fasta")
    contam_removed_path = os.path.join(path, '3contam_removal')
    clean_files = glob(clean_search)

    if not clean_files:
        print("Could not find cleaned fasta files\n"
              "It is possible that there were no sequences remaining after removal of sequences with degenerate bases\n"
              "If you specified -s (remove sequences with stop codons, did you specify the correct reading frame?\n"
              "Do your sequences extend past the end of the reading frame")
    call_contam_check(clean_files, contam_removal_script, contam_removed_path, gene_region, logfile)

    # get the HXB2 sequence for the gene region
    hxb2_file = os.path.join(script_folder, "HXB2_seqs.fasta")
    hxb2 = fasta_to_dct(hxb2_file)
    hxb2_gene = "HXB2_" + gene_region.split("_")[0]
    hxb2_seq = hxb2[hxb2_gene]

    # cat all cleaned files into one file + the relevant HXB2 sequence
    print("merging all cleaned and contam removed fasta files into one file")
    all_clean_path = os.path.join(path, '3contam_removal')
    clean_name = name + "_" + gene_region + "_all.fasta"
    all_cleaned_outname = os.path.join(all_clean_path, clean_name)
    cleaned_files_search = os.path.join(contam_removed_path, '*_good.fasta')
    cleaned_files = glob(cleaned_files_search)
    if not cleaned_files:
        print("No cleaned fasta files were found\n"
              "Check that the fasta files still have sequences in them after the removal of bad sequences")
        sys.exit()

    with open(all_cleaned_outname, 'w') as outfile:
        outfile.write(">{0}\n{1}\n".format(hxb2_gene, hxb2_seq))
        for fasta_file in cleaned_files:
            with open(fasta_file) as infile:
                for line in infile:
                    outfile.write(line + "\n")

    # move concatenated file to 4aligned
    print("moving concatenated file to 4aligned folder")
    aln_path = os.path.join(path, '4aligned')
    move_file = os.path.join(aln_path, clean_name)
    copyfile(all_cleaned_outname, move_file)
    os.remove(all_cleaned_outname)

    # call alignment script
    print("Aligning the sequences")
    to_align = move_file
    inpath, fname = os.path.split(to_align)
    fname = fname.replace(".fasta", "_aligned.fasta")
    call_align(envelope, script_folder, to_align, aln_path, fname)

    # call funcion to calculate sequencing stats
    print("Calculating alignment stats")
    call_stats_calc = os.path.join(script_folder, 'ngs_stats_calculator.py')
    stats_outfname = (name + "_" + gene_region + '_sequencing_stats.csv')
    stats_outpath = os.path.join(path, stats_outfname)
    cmd6 = 'python3 {0} -i {1} -o {2}'.format(call_stats_calc, path, stats_outpath)
    subprocess.call(cmd6, shell=True)

    print("The sample processing has been completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Call make_folders script BEFORE running this script '
                                                 'Then copy your data into the /1raw/ folder'
                                                 'This script runs the NGS data processing pipeline. '
                                                 'It is a good idea to run this script using screen or nohup',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--path', default=argparse.SUPPRESS, type=str,
                        help='The path to the gene region subfolder (GAG_1/ or env_C1C2/ or POL_1/...)', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the prefix name of your outfile', required=True)
    parser.add_argument('-g', '--gene_region', default=argparse.SUPPRESS, type=str,
                        help='the genomic region being sequenced, '
                             'valid options: GAG_1/GAG_2/ENV_C1C2/POL_1/NEF_1 etc..', required=True)
    parser.add_argument('-sf', '--script_folder', default=argparse.SUPPRESS, type=str,
                        help='the path to the folder containing the pipeline scripts', required=True)
    parser.add_argument('-f', '--fwd_primer', default=argparse.SUPPRESS, type=str,
                        help='The fwd primer for these samples (eg: NNNNGGAAATATGGAAAGGAAGGAC)', required=True)
    parser.add_argument('-r', '--cDNA_primer', default=argparse.SUPPRESS, type=str,
                        help='The cDNA primer for these samples (eg: NNNNNNNNNNNTCTTCTAATACTGTATCATCTG)', required=True)
    parser.add_argument('-fr', '--frame', type=int,
                        help='The reading frame (1, 2 or 3)', required=False)
    parser.add_argument('-s', '--stops', default=False, action='store_true',
                        help='Remove sequences with stop codons?)', required=False)
    parser.add_argument('-l', '--length', type=int,
                        help='The minimum read length)', required=False)
    parser.add_argument('-e', '--envelope', default=None, nargs="+",
                        help='If your sequences are of HIV envelope, which V-loops are in the sequence?'
                             '(eg: V1 V2) (options include: V1, V2 , V3, V4, V5)', required=False)

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
