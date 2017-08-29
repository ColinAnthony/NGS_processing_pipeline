#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse


__author__ = 'colin'


def main(your_path, fname):
    print(your_path)

    main_dir = os.path.join(your_path, fname)
    second_level_dirs = ['1raw',  '2consensus',  '3cleaned', '4aligned',  '5haplotype',  '6analysis']
    third_level_dirs = ['aa_frq', 'divergence', 'entropy', 'glycans', 'loops', 'tree']

    if not os.path.exists(main_dir):
        os.makedirs(main_dir)

    for folder in second_level_dirs:
        if folder == '6analysis':
            for nested_folder in third_level_dirs:
                set_path =  os.path.join(main_dir, folder)
                make_folder = os.path.join(set_path, nested_folder)
                os.makedirs(make_folder)
        else:
            make_folder = os.path.join(main_dir, folder)
            os.makedirs(make_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='create folder structure for NGS data cleanup/analysis',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-p', '--path', default=argparse.SUPPRESS, type=str,
                        help='The path where the folders will be created', required=True)
    parser.add_argument('-n', '--name', default=argparse.SUPPRESS, type=str,
                        help='the name of the participant', required=True)

    args = parser.parse_args()
    path = args.path
    name = args.name

    main(path, name)
