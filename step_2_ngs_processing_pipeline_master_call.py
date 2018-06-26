#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import os
import sys
from shutil import copyfile
from shutil import rmtree
from distutils.dir_util import copy_tree
import argparse
import subprocess
from glob import glob
import re
from itertools import groupby
import collections


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


def rename_sequences(raw_files_search):
    """
    rename raw files to cleaner name ending in R1.fastq or R2.fastq
    :param raw_files_search: list of files to rename
    :return:
    """
    print("renaming raw files")
    for inf_R1 in raw_files_search:
        inf_R2 = inf_R1.replace("R1_001.fastq", "R2_001.fastq")
        path = os.path.split(inf_R1)[0]
        inf_R1_name = os.path.split(inf_R1)[-1]
        inf_R2_name = os.path.split(inf_R2)[-1]

        outf_R1 = inf_R1_name.replace("-", "_")
        outf_R2 = inf_R2_name.replace("-", "_")
        outf_R1_rename = re.sub("S[0-9]+_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq", "R1.fastq", outf_R1)
        outf_R1_rename_with_path = os.path.join(path, outf_R1_rename)

        if outf_R1_rename == outf_R1:
            # check if they have already been renamed
            if outf_R1.split("_")[-1] == "R1.fastq":
                print("file was already in correct format?")
                return
            else:
                print("Unable to rename R1 file {0}\ncheck the file renaming regex".format(inf_R1))
                sys.exit()
        else:
            os.rename(inf_R1, outf_R1_rename_with_path)
        print(outf_R1_rename)

        outf_R2_rename = re.sub("S[0-9]+_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9].fastq", "R2.fastq", outf_R2)
        outf_R2_rename_with_path = os.path.join(path, outf_R2_rename)
        if outf_R2_rename == outf_R2:
            if outf_R2.split("_")[-1] == "R2.fastq":
                print("file was already in correct format?")
                return
            else:
                print("Unable to rename R2 file {0}\ncheck the file renaming regex".format(inf_R2))
                sys.exit()
        else:
            os.rename(inf_R2, outf_R2_rename_with_path)


def call_motifbinner(raw_files, motifbinner, cons_outpath, fwd_primer, cDNA_primer, nonoverlap, counter, logfile):
    """
    function to pass args to the script that calls the motifbinner2
    :param raw_files: (list) of all the read 1 files
    :param motifbinner: (str) call motifbinner script name
    :param cons_outpath: (str) desired outpath
    :param fwd_primer: (str) of fwd primer
    :param cDNA_primer: (str) of rev/cDNA primer
    :param nonoverlap: (bool) False for overlapping read 1 and 2, True of read 1 and 2 don't overlap
    :param counter: (int) count of number of times the script has been called (so that we only write to log once)
    :param logfile: (str) path and name of the log file
    :return:
    """

    if type(raw_files) is not list:
        raise TypeError('Expected list of raw files, got: ', raw_files)

    if nonoverlap:
        overlap_flag = "-v"
    else:
        overlap_flag = ""

    for file in raw_files:
        read1 = file
        read2 = file.replace("R1.fastq", "R2.fastq")
        name_prefix = os.path.split(file)[-1].replace("_R1.fastq", "")

        cmd1 = 'python3 {0} -r1 {1} -r2 {2} -o {3} -f {4} -r {5} -n {6} -c {7} -l {8} {9} '.format(motifbinner,
                                                                                              read1,
                                                                                              read2,
                                                                                              cons_outpath,
                                                                                              fwd_primer,
                                                                                              cDNA_primer,
                                                                                              name_prefix,
                                                                                              counter,
                                                                                              logfile,
                                                                                              overlap_flag)

        try:
            subprocess.call(cmd1, shell=True)
            counter += 1
        except Exception as e:
            raise e


def delete_gaps(fasta_infiles):
    """
    removes gap characters from binned consensus sequence fasta file
    :param fasta_infiles: list of consensus sequence fasta file
    :return: None
    """

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


def call_fasta_cleanup(consensus_fasta, remove_bad_seqs, clean_path, length, logfile):
    """
    function to pass args to remove bad sequences script
    :param consensus_fasta: (str) list of binned consensus sequence fasta files
    :param remove_bad_seqs: (str) name of script to run
    :param clean_path: (str) desired outpath
    :param length: (int) min length of sequence allowed
    :param logfile: the path and name of the log file
    :return:
    """

    for fasta_file in consensus_fasta:
        cmd4 = 'python3 {0} -i {1} -o {2} -l {3} -lf {4}'.format(remove_bad_seqs,
                                                                        fasta_file,
                                                                        clean_path,
                                                                        length,
                                                                        logfile)
        if os.path.exists(logfile):
            with open(logfile, 'a') as handle:
                handle.write("\nremove_bad_sequences commands:\n{}\n".format(cmd4))

        subprocess.call(cmd4, shell=True)


def call_contam_check(consensuses, contam_removal_script, contam_removed_path, gene_region, logfile):
    """
    function to pass args to contam check script
    :param consensuses: list of cleaned fasta files
    :param contam_removal_script: path to script
    :param contam_removed_path: output path location
    :param gene_region: the gene region
    :param logfile: the path and name of the log file
    :return: None
    """

    for consensus_file in consensuses:
        cmd3 = 'python3 {0} -i {1} -o {2} -g {3} -l {4}'.format(contam_removal_script,
                                                                consensus_file,
                                                                contam_removed_path,
                                                                gene_region,
                                                                logfile)

        subprocess.call(cmd3, shell=True)


def call_align(script_folder, to_align, aln_path, fname, ref, gene, sub_region, user_ref):
    """
    function to call the alignment script
    :param script_folder: (str) path to the folder containing the repo scripts
    :param to_align: (str) file to align
    :param aln_path: (str) the output path
    :param fname: (str) the prefix for the output file name
    :param ref: (str) the reference type
    :param gene_region: (str) the HIV gene region
    :param sub_region: (str) the HIV gene sub-region
    :param user_ref: (str) path to the user reference fasta file
    :return: None
    """
    align_function = os.path.join(script_folder, 'align_ngs_codons.py')
    if sub_region:
        reg = "-reg {}".format(sub_region)
    else:
        reg = ""
    if not user_ref:
        usr_ref = ""
    else:
        usr_ref = "-u {}".format(user_ref)

    cmd5 = 'python3 {0}  -i {1} -o {2} -n {3} -r {4} -g {5} -v {6} {7}'.format(align_function, to_align, aln_path, fname,
                                                                            ref, gene, reg, usr_ref)
    subprocess.call(cmd5, shell=True)


def main(path, name, gene_region, sub_region, fwd_primer, cDNA_primer, nonoverlap, length, run_step,
         run_only, user_ref):

    get_script_path = os.path.realpath(__file__)
    script_folder = os.path.split(get_script_path)[0]
    script_folder = os.path.abspath(script_folder)

    path = os.path.abspath(path)
    gene = gene_region.split("_")[0]
    # define logfile filename
    logfile = os.path.join(path, (gene_region + "_logfile.txt"))
    if not os.path.isfile(logfile):
        # initialize the log file
        with open(logfile, 'w') as handle:
            handle.write("Log File,{0}_{1}\n".format(name, gene_region))

    folders_to_make = ['0raw_temp', '1consensus_temp', '2cleaned_temp', '3contam_removal_temp']
    print("running pipeline from step:", run_step)
    initual_run_step = run_step
    if run_step <= 4:
        # make temp dirs
        print("making temp folders")

        for folder in folders_to_make:
            flder = os.path.join(path, folder)
            if os.path.isdir(flder):
                print("Deleting exisiting folders")
                rmtree(flder)
            os.makedirs(flder, exist_ok=True)

    new_data = os.path.join(path, "0new_data")

    # Step 1: rename the raw sequences
    if run_step == 1:
        # move files from new_data to 0raw_temp
        files_to_move = os.path.join(new_data, "*.fastq")
        move_folder = os.path.join(path, '0raw_temp')
        for file in glob(files_to_move):
            file_name = os.path.split(file)[-1]
            move_location = os.path.join(move_folder, file_name)
            copyfile(file, move_location)

        # do the renaming
        raw_fastq_inpath = os.path.join(path, '0raw_temp')
        raw_files_search = os.path.join(raw_fastq_inpath, "*R1*.fastq")
        raw_files = glob(raw_files_search)
        if not raw_files:
            print("No raw files were found\n"
                  "Check that files end with R1.fastq and R2.fastq")
            sys.exit()

        rename_sequences(raw_files)
        run_step += 1

        if run_only:
            # copy back to permanent folder, remove temp folder
            run_step = 10

    # Step 2: run the call_MotifBinner script which will loop over fastq files in the target folder
    if run_step == 2:
        move_folder = os.path.join(path, '0raw_temp')
        if initual_run_step == 2:
            files_to_move = os.path.join(new_data, "*.fastq")
            for file in glob(files_to_move):
                file_name = os.path.split(file)[-1]
                move_location = os.path.join(move_folder, file_name)
                copyfile(file, move_location)
        
        motifbinner = os.path.join(script_folder, 'call_motifbinner.py')
        rename_in_search = os.path.join(move_folder, "*_R1.fastq")
        rename_in = glob(rename_in_search)
        cons_outpath = os.path.join(path, '1consensus_temp', 'binned')
        counter = 0
        try:
            call_motifbinner(rename_in, motifbinner, cons_outpath, fwd_primer, cDNA_primer, nonoverlap, counter,
                             logfile)
            run_step += 1
        except Exception as e:
            print(e)
            print("now quitting..")
            sys.exit()

        # check if the consensus files exist
        nested_consensuses_path = os.path.join(path,
                                               '1consensus_temp/binned/*/*_buildConsensus/*_buildConsensus.fastq')
        nested_consesnsuses = glob(nested_consensuses_path)
        if not nested_consesnsuses:
            print("No consensus sequences were found\n"
                  "This is likely if MotifBinner was not able to complete\n"
                  "Do your fastq sequences end in R1.fastq/R2.fastq?\n"
                  "Check that the primer sequences are correct\n"
                  "Check the binning report in the 1consensus/binned/ folder")
            sys.exit()

        # copy data from nested binned folders into 1consensus folder
        print("Coping fastq files from nested folders to '1consensus_temp' folder")
        consensus_path = os.path.join(path, '1consensus_temp')
        for cons_file in nested_consesnsuses:
            old_path, old_name = os.path.split(cons_file)
            new_name1 = old_name.replace("_buildConsensus", "")
            new_name = new_name1.replace("_kept", "")
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
        consensus_path = os.path.join(path, '1consensus_temp')
        consensus_search = os.path.join(consensus_path, '*.fasta')
        consensus_infiles = glob(consensus_search)
        delete_gaps(consensus_infiles)

        if run_only:
            # copy back to permanent folder, remove temp folder
            run_step = 10

    # Step 3: call remove bad sequences
    if run_step == 3:
        move_folder = os.path.join(path, '1consensus_temp')
        if initual_run_step == 3:
            files_to_move = os.path.join(new_data, "*.fasta")
            for file in glob(files_to_move):
                file_name = os.path.split(file)[-1]
                move_location = os.path.join(move_folder, file_name)
                copyfile(file, move_location)

        if nonoverlap:
            print("move folder", move_folder)
            search_fwd_rev = os.path.join(move_folder, "*rev.fasta")
            for file in glob(search_fwd_rev):
                out = file + "_temp.fasta"
                print(out)
                print(file)
                cmd_rev_comp = 'seqmagick convert --reverse-complement {0} {1}'.format(file, out)
                subprocess.call(cmd_rev_comp, shell=True)
                os.unlink(file)
                os.rename(out, file)
            input("enter")
        print("Removing 'bad' sequences")
        remove_bad_seqs = os.path.join(script_folder, 'remove_bad_sequences.py')
        consensus_search = os.path.join(move_folder, '*.fasta')
        consensus_infiles = glob(consensus_search)
        clean_path = os.path.join(path, '2cleaned_temp')
        print(consensus_search)
        if not consensus_infiles:
            print("Could not find consensus fasta files\n"
                  "It is possible something went wrong when copying the consensus sequences from the nested folders "
                  "to the 1consensus folder")
            sys.exit()

        call_fasta_cleanup(consensus_infiles, remove_bad_seqs, clean_path, length, logfile)
        run_step += 1

        if run_only:
            # copy back to permanent folder, remove temp folder
            run_step = 10

    # Step 4: remove contaminating sequences
    if run_step == 4:
        move_folder = os.path.join(path, '2cleaned_temp')
        if initual_run_step == 4:
            files_to_move = os.path.join(new_data, "*.fasta")
            for file in glob(files_to_move):
                file_name = os.path.split(file)[-1]
                move_location = os.path.join(move_folder, file_name)
                copyfile(file, move_location)

        print("removing contaminating non-HIV sequences")
        contam_removal_script = os.path.join(script_folder, "contam_removal.py")
        clean_search = os.path.join(move_folder, "*clean.fasta")
        contam_removed_path = os.path.join(path, '3contam_removal_temp')
        clean_files = glob(clean_search)
        hxb2_region = {"GAG": "GAG", "POL": "POL", "PRO": "POL", "RT": "POL", "RNASE": "POL", "INT": "POL",
                       "ENV": "ENV", "GP120": "ENV", "GP41": "ENV", "NEF": "NEF",
                       "VIF": "VIF", "VPR": "VPR", "REV": "REV", "VPU": "VPU"}

        region_to_check = hxb2_region[gene]
        if not clean_files:
            print("Could not find cleaned fasta files\n"
                  "It is possible there were no sequences remaining after removal of sequences with degenerate bases\n"
                  )

        call_contam_check(clean_files, contam_removal_script, contam_removed_path, region_to_check, logfile)

        # copy back to permanent folder, remove temp folder
        run_step = 10

    # copy data to permanent folders, delete temp folders/files
    if run_step == 10:
        print('copying data from temp folders')
        for folder in folders_to_make:
            temp_folder = os.path.join(path, folder)
            perm_folder = temp_folder.replace("_temp", "")
            if not os.path.isdir(temp_folder):
                print("temp folder does not exist")
                sys.exit()
            if not os.path.isdir(perm_folder):
                print("folder does not exist")
                sys.exit()

            copy_tree(temp_folder, perm_folder)
            rmtree(temp_folder)

        # clear the new_data folder
        condition = True
        if run_only:
            print('removing temp files')
            while condition:
                response = input("remove files from 0new_data folder? (yes or no)")
                if response.lower() == "yes" or response.lower() == "y":
                    condition = False
                    new_files_to_remove = os.path.join(new_data, "*")
                    for file in glob(new_files_to_remove):
                        os.unlink(file)
                elif response.lower() == "no" or response.lower() == "n":
                    condition = False
                    print("not deleting files")
                else:
                    print("response not valid, please enter either: 'yes', or 'y' to delete or 'no', or 'n' to keep")
            sys.exit()

        else:
            new_files_to_remove = os.path.join(new_data, "*")
            for file in glob(new_files_to_remove):
                os.unlink(file)
            run_step = 5

    # Step 5: set things up to align sequences
    if run_step == 5:

        # cat all cleaned files into one file + the relevant HXB2 sequence
        print("merging all cleaned and contam removed fasta files into one file")
        all_clean_path = os.path.join(path, '3contam_removal')
        contam_removed_path = os.path.join(path, '3contam_removal')
        if nonoverlap:
            clean_name_fwd = name + "_" + gene_region + "_fwd_all.fasta"
            clean_name_rev = name + "_" + gene_region + "_rev_all.fasta"
            all_cleaned_outname_fwd = os.path.join(all_clean_path, clean_name_fwd)
            all_cleaned_outname_rev = os.path.join(all_clean_path, clean_name_rev)

            cleaned_files_search_fwd = os.path.join(contam_removed_path, '*fwd_good.fasta')
            cleaned_files_search_rev = os.path.join(contam_removed_path, '*rev_good.fasta')

            cleaned_files_fwd = glob(cleaned_files_search_fwd)
            cleaned_files_ref = glob(cleaned_files_search_rev)

            if not cleaned_files_fwd:
                sys.exit("No cleaned fwd-fasta files were found\n"
                      "Check that the fasta files still have sequences in them after the removal of bad sequences")
            if not cleaned_files_ref:
                sys.exit("No cleaned rev-fasta files were found\n"
                      "Check that the fasta files still have sequences in them after the removal of bad sequences")
        else:
            clean_name = name + "_" + gene_region + "_all.fasta"
            all_cleaned_outname = os.path.join(all_clean_path, clean_name)

            cleaned_files_search = os.path.join(contam_removed_path, '*_good.fasta')
            cleaned_files = glob(cleaned_files_search)

            if not cleaned_files:
                print("No cleaned fasta files were found\n"
                      "Check that the fasta files still have sequences in them after the removal of bad sequences")
                sys.exit()

        if nonoverlap:
            with open(all_cleaned_outname_fwd, 'w') as outfile:
                for fasta_file in cleaned_files_fwd:
                    with open(fasta_file) as infile:
                        for line in infile:
                            outfile.write(line + "\n")

            with open(all_cleaned_outname_rev, 'w') as outfile:
                for fasta_file in cleaned_files_ref:
                    with open(fasta_file) as infile:
                        for line in infile:
                            outfile.write(line + "\n")
        else:
            with open(all_cleaned_outname, 'w') as outfile:
                for fasta_file in cleaned_files:
                    with open(fasta_file) as infile:
                        for line in infile:
                            outfile.write(line + "\n")

        # move concatenated file to 4aligned
        print("moving concatenated file to 4aligned folder")
        aln_path = os.path.join(path, '4aligned')
        if nonoverlap:
            move_file_fwd = os.path.join(aln_path, clean_name_fwd)
            copyfile(all_cleaned_outname_fwd, move_file_fwd)
            os.unlink(all_cleaned_outname_fwd)
            move_file_rev = os.path.join(aln_path, clean_name_rev)
            copyfile(all_cleaned_outname_rev, move_file_rev)
            os.unlink(all_cleaned_outname_rev)
        else:
            move_file = os.path.join(aln_path, clean_name)
            copyfile(all_cleaned_outname, move_file)
            os.unlink(all_cleaned_outname)

        # call alignment script
        print("Aligning the sequences")
        if nonoverlap:
            for move_file in [move_file_fwd, move_file_rev]:
                to_align = move_file
                inpath, fname = os.path.split(to_align)
                fname = fname.replace(".fasta", "")
                ref = "CONSENSUS_C"

                # infile, outpath, name, ref, gene, var_align, sub_region, user_ref)
                call_align(script_folder, to_align, aln_path, fname, ref, gene_region, sub_region, user_ref)

                # translate alignment
                transl_name = fname.replace("_aligned.fasta", "_aligned_translated.fasta")
                cmd = "seqmagick convert --sort length-asc --upper --translate dna2protein --line-wrap 0 {0} {1}".format(fname, transl_name)
                subprocess.call(cmd, shell=True)

        else:
            to_align = move_file
            inpath, fname = os.path.split(to_align)
            fname = fname.replace(".fasta", "")
            ref = "CONSENSUS_C"

            # infile, outpath, name, ref, gene, var_align, sub_region, user_ref)
            call_align(script_folder, to_align, aln_path, fname, ref, gene, sub_region, user_ref)

            # translate alignment
            fname = to_align.replace(".fasta", "_aligned.fasta")
            transl_name = fname.replace("_aligned.fasta", "_aligned_translated.fasta")
            cmd = "seqmagick convert --sort length-asc --upper --translate dna2protein --line-wrap 0 {0} {1}".format(
                fname, transl_name)
            subprocess.call(cmd, shell=True)

        run_step += 1

        if run_only:
            sys.exit()

    # call funcion to calculate sequencing stats
    if run_step == 6:
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
                        help='The path to the gene region subfolder (GAG_1/ or ENV_C1C2/ or POL_1/...)', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the prefix name of your outfile. Usually the participant name', required=True)
    parser.add_argument('-g', '--gene_region', default="ENV", type=str,
                        choices=["ENV", "ENV_1", "ENV_2", "ENV_3", "ENV_4", "ENV_5", "GAG", "GAG_1", "GAG_2", "POL",
                                 "POL_1", "POL_2", "POL_3", "POL_4", "POL_5", "PRO", "NEF_1", "VIF_1", "VPR", "REV", "VPU"],
                        help='the genomic region being sequenced, '
                             'valid options: etc..', required=True)
    parser.add_argument('-reg', '--regions', default=False, action="store",
                        choices=["C0C1", "C1C2", "C2C3", "C3C5", "GP41", "GP120", "GP160", "P17", "P24"], type=str,
                        help='the variable regions in your data', required=False)
    parser.add_argument('-f', '--fwd_primer', default=argparse.SUPPRESS, type=str,
                        help='The fwd primer for these samples (eg: NNNNGGAAATATGGAAAGGAAGGAC)', required=False)
    parser.add_argument('-r', '--cDNA_primer', default=argparse.SUPPRESS, type=str,
                        help='The cDNA primer for these samples (eg: NNNNNNNNNNNTCTTCTAATACTGTATCATCTG)', required=False)
    parser.add_argument('-l', '--length', type=int,
                        help='The minimum read length)', required=False)
    parser.add_argument('-v', '--nonoverlap', default=False, action='store_true',
                        help="Use if reads don't overlap)", required=False)
    parser.add_argument('-u', '--user_ref', default=False, type=str,
                        help='the path and file name for the custom DNA reference sequence for codon aligning\n'
                             'must start in reading frame 1', required=False)
    parser.add_argument('-rs', '--run_step', type=int, default=1,
                        help='rerun the pipeline from a given step:\n'
                             '1 = step 1: rename raw files;\n'
                             '2 = step 2: run MotifBinner2;\n'
                             '3 = step 3: clean consensus sequences;\n'
                             '4 = step 4: remove contam sequences;\n'
                             '5 = step 5: align the sequences;\n'
                             '6 = step 6: calculate sequencing depth stats for each step of pipeline ', required=False)
    parser.add_argument('-ro', '--run_only', default=False, action='store_true',
                        help='run only the specified run_step)', required="--run_step" in sys.argv)

    args = parser.parse_args()
    path = args.path
    name = args.name
    gene_region = args.gene_region
    regions = args.regions
    fwd_primer = args.fwd_primer
    cDNA_primer = args.cDNA_primer
    nonoverlap = args.nonoverlap
    length = args.length
    run_step = args.run_step
    run_only = args.run_only
    user_ref = args.user_ref
    if gene_region == "ENV":
        if not regions:
            sys.exit("must use the -reg flag for ENV")

    main(path, name, gene_region, regions, fwd_primer, cDNA_primer, nonoverlap, length, run_step,
         run_only, user_ref)
