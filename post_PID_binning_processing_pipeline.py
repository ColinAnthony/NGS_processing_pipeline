#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import subprocess
from subprocess import DEVNULL


__author__ = 'Colin Anthony'


def main(path, name, script_folder, gene_region):

    # call make folders script
    create_folders = os.path.join(script_folder, 'create_folders.py')
    cmd1 = 'python {0} -p {1} -n {2} -g {3}'.format(path, name, create_folders, gene_region)
    subprocess.call(cmd1, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    main_path = os.path.join(path, name)
    logfile = os.path.join(main_path, "logfile.txt")
    with open(logfile, 'w') as handle:
        handle.write("Initializing log file")
    # cd into 1raw
    cmd1 = 'cd {}'.format()
    subprocess.call()
    cmd2 =
    # call binner
    # pull data into folders
    # call remove bad sequences
    # call cat all files into one
    # call align all samples

    cmd = 'mafft {0} > {1}'.format(tmp_file_in, tmp_file_out)
    subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--path', default=argparse.SUPPRESS, type=str,
                        help='The path where the folders will be created', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the name of the participant', required=True)
    parser.add_argument('-g', '--gene_region', default=argparse.SUPPRESS, type=str,
                        help='the genomic region being sequenced', required=True)
    parser.add_argument('-f', '--script_folder', default=argparse.SUPPRESS, type=str,
                        help='the path to the folder containing the pipeline scripts', required=True)

    args = parser.parse_args()
    path = args.path
    name = args.name
    script_folder = args.script_folder
    gene_region = args.gene_region

    main(path, name, script_folder, gene_region)
