#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import subprocess
from subprocess import DEVNULL


__author__ = 'Colin Anthony'


def main(path, name, script_folder):

    # call make folders script
    # cd into top level
    logfile = os.path.join(path, "logfile.txt")
    with open(logfile, 'w') as handle:
        handle.write("Initializing log file")
    # cd into 1raw
    cmd1 =
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
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=list,
                        help='the name of the participant', required=True)
    parser.add_argument('-f', '--script_folder', default=argparse.SUPPRESS, type=list,
                        help='the path to the folder containing the pipeline scripts', required=True)

    args = parser.parse_args()
    path = args.path
    name = args.name
    script_folder = args.script_folder

    main(path, name, script_folder)
