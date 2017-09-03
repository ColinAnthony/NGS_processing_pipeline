#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import sys
import argparse
import collections
import tempfile
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import regex



__author__ = 'Colin Anthony'


def fasta_to_dct(fn):
    '''
    converts a fasta file to a dictionary where key = seq name, value = sequence
    :param fn: a fasta file
    :return: a dictionary
    '''
    dct = collections.OrderedDict()
    for seq_record in SeqIO.parse(open(fn), "fasta"):
        dct[seq_record.description.replace(" ", "_").upper()] = str(seq_record.seq).replace("~", "-").upper()
    return dct


def blastn_seqs(query_sequence):
    '''
    :param inseq: (str) the sequence to blast
    :return: (bool) True or False depending on whether sequence is hiv or not
    '''
    tmp_dir = tempfile.gettempdir()
    tmpfile = os.path.join(tmp_dir, "tmp_blast_results.xml")

    blast_results = NCBIWWW.qblast('blastn', 'nr', query_sequence)
    # blastn_cline = NcbiblastnCommandline(query=query_sequence, db="nr", evalue=0.001, outfmt=5. perc_identity=50,
    #                                      out=tmpfile, max_target_seqs=3)
    # stdout, stderr = blastn_cline()

    with open(tmpfile, 'w') as handle:
        handle.write(blast_results.read())

    e_value_threshold = 0.04

    for blast_record in NCBIXML.parse(open(tmpfile)):
        # skip empty records
        if blast_record.alignments:
            for hit in blast_record.alignments:
                # todo is this really taking the top hit? is the xml sorted by e value?
                top = hit.hsps[0]
                if top.expect < e_value_threshold:
                    if "HIV-1" in hit.title.split(" "):
                        return True
                    else:
                        return False


def regex_contam(query_sequence, hxb2_region):
    '''
    :param inseq: (str) the sequence to blast
    :return: (bool) True or False depending on whether sequence is hiv or not
    '''
    # todo hardcode hxb2 regex string for GAG_2? for now
    # hxb2_region = ''
    error = int(len(query_sequence) * 0.6)
    hxb2_regex = "r'({0}){{{1}}}'".format(hxb2_region, error)
    query = str(query_sequence)
    match = regex.search(hxb2_regex, query, regex.BESTMATCH)
    if match is not None:
        print("no match")
        return True
    else:
        return False


def main(consensus, outpath):

    # initialize file names
    cln_cons_name = os.path.split(consensus)[-1]
    cln_cons = cln_cons_name.replace(".fasta", "contam_rem.fasta")
    consensus_out = os.path.join(outpath, cln_cons)
    contam_seqs = cln_cons_name.replace(".fasta", "contam_seqs.fasta")
    contam_out = os.path.join(outpath, contam_seqs)

    # clear out any existing outfiles for consensus
    with open(consensus_out, 'w') as handle:
        handle.write("")
    with open(contam_out, 'w') as handle:
        handle.write("")

    # get hxb2 seq for regex??
    hxb2_seq = "gag1"
    hxb2_seq = "gag2"

    # store all consensus seqs in a dict
    all_sequences_d = fasta_to_dct(consensus)

    # check each sequence to see if it is a contaminating sequence
    for name, seq in all_sequences_d.items():
        # is_contam = regex_contam(seq, hxb2_seq)
        is_contam = blastn_seqs(seq)

        # if is contam save to contam file
        if is_contam:
            with open(contam_out, 'a') as handle1:
                outstr = "{0}{1}\n{2}\n".format('>', name, seq)
                handle1.write(outstr)

        # if is not contam, to write to file
        else:
            with open(consensus_out, 'a') as handle2:
                outstr = "{0}{1}\n{2}\n".format('>', name, seq)
                handle2.write(outstr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--consensus', default=argparse.SUPPRESS, type=str,
                        help='The read1 (R1) fastq file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be copied', required=True)

    args = parser.parse_args()
    consensus = args.consensus
    read2 = args.read2
    outpath = args.outpath

    main(consensus, outpath)
