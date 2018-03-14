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

    print("renaming raw files")
    for inf_R1 in raw_files_search:
        print(inf_R1) # /data/hvtn503/breakthrough/GAG_2/0raw_temp/159400956_1110_013wpi_GAG_2_NN_S104_L001_R1_001.fastq
        inf_R2 = inf_R1.replace("R1_001.fastq", "R2_001.fastq")
        path = os.path.split(inf_R1)[0]
        inf_R1_name = os.path.split(inf_R1)[-1]  # 159400956_1110_013wpi_GAG_2_NN_S104_L001_R1_001.fastq
        inf_R2_name = os.path.split(inf_R2)[-1]

        outf_R1 = inf_R1_name.replace("-", "_")  # 159400956_1110_013wpi_GAG_2_NN_S104_L001_R1_001.fastq
        outf_R2 = inf_R2_name.replace("-", "_")
        outf_R1_rename = re.sub("S[0-9]+_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq", "R1.fastq", outf_R1)
        outf_R1_rename_with_path = os.path.join(path, outf_R1_rename)
        # /data/hvtn503/breakthrough/GAG_2/0raw_temp/159400956_1110_013wpi_GAG_2_NN_S104_L001_R1_001.fastq


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


def call_fasta_cleanup(contam_removed_fasta, remove_bad_seqs, clean_path, frame, length, logfile, stops):

    for fasta_file in contam_removed_fasta:
        if stops:
            cmd4 = 'python3 {0} -i {1} -o {2} -f {3} -s -l {4} -lf {5}'.format(remove_bad_seqs,
                                                                                   fasta_file,
                                                                                   clean_path,
                                                                                   frame,
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
        loops = " ".join(envelope)
        align_all = os.path.join(script_folder, 'align_all_env_samples.py')

        cmd5 = 'python3 {0}  -i {1} -o {2} -n {3} -l {4}'.format(align_all, to_align, aln_path, fname, loops)
    else:
        align_all = os.path.join(script_folder, 'align_all_samples.py')
        cmd5 = 'python3 {0}  -i {1} -o {2} -n {3}'.format(align_all, to_align, aln_path, fname)

    subprocess.call(cmd5, shell=True)


def main(path, name, gene_region, fwd_primer, cDNA_primer, nonoverlap, frame, stops, length, envelope, run_step,
         run_only):
    get_script_path = os.path.realpath(__file__)
    script_folder = os.path.split(get_script_path)[0]
    path = os.path.abspath(path)
    script_folder = os.path.abspath(script_folder)

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
            print(move_folder)

            for file in glob(move_folder):
                print(file)
                if "rev.fasta" in file.split("_"):
                    out = file + ".temp"
                    print(out)
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

        call_fasta_cleanup(consensus_infiles, remove_bad_seqs, clean_path, frame, length, logfile, stops)
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

        gene = gene_region.split("_")[0]
        region_to_check = hxb2_region[gene]
        if not clean_files:
            print("Could not find cleaned fasta files\n"
                  "It is possible there were no sequences remaining after removal of sequences with degenerate bases\n"
                  "If you specified -s (remove sequences with stop codons, did you specify the correct reading frame?\n"
                  "Do your sequences extend past the end of the reading frame")

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
        contam_removed_path = os.path.join(path, '3contam_removal')
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
        os.unlink(all_cleaned_outname)

        # call alignment script
        print("Aligning the sequences")
        to_align = move_file
        inpath, fname = os.path.split(to_align)
        fname = fname.replace(".fasta", "_aligned.fasta")
        call_align(envelope, script_folder, to_align, aln_path, fname)
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
                        help='the prefix name of your outfile', required=True)
    parser.add_argument('-g', '--gene_region', default=argparse.SUPPRESS, type=str,
                        help='the genomic region being sequenced, '
                             'valid options: GAG_1/GAG_2/ENV_C1C2/POL_1/NEF_1 etc..', required=True)
    parser.add_argument('-f', '--fwd_primer', default=argparse.SUPPRESS, type=str,
                        help='The fwd primer for these samples (eg: NNNNGGAAATATGGAAAGGAAGGAC)', required=False)
    parser.add_argument('-r', '--cDNA_primer', default=argparse.SUPPRESS, type=str,
                        help='The cDNA primer for these samples (eg: NNNNNNNNNNNTCTTCTAATACTGTATCATCTG)', required=False)
    parser.add_argument('-fr', '--frame', type=int,
                        help='The reading frame (1, 2 or 3)', required=False)
    parser.add_argument('-s', '--stops', default=False, action='store_true',
                        help='Remove sequences with stop codons?)', required=False)
    parser.add_argument('-l', '--length', type=int,
                        help='The minimum read length)', required=False)
    parser.add_argument('-e', '--envelope', type=str, default=None, nargs="+",
                        help='If your sequences are of HIV envelope, which V-loops are in the sequence?'
                             '(eg: V1 V2) (options include: V1, V2 , V3, V4, V5)', required=False)
    parser.add_argument('-v', '--nonoverlap', default=False, action='store_true',
                        help="Use if reads don't overlap)", required=False)
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
    fwd_primer = args.fwd_primer
    cDNA_primer = args.cDNA_primer
    nonoverlap = args.nonoverlap
    frame = args.frame
    stops = args.stops
    length = args.length
    envelope = args.envelope
    run_step = args.run_step
    run_only = args.run_only

    main(path, name, gene_region, fwd_primer, cDNA_primer, nonoverlap, frame, stops, length, envelope, run_step,
         run_only)
