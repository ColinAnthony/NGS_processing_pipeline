#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
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
        dct[seq_record.description.replace(" ", "_").upper()] = str(seq_record.seq).replace("-", "").upper()
    return dct


def blastn_seqs(query_file, logfile):
    '''
    :param inseq: (str) the sequence to blast
    :return: (bool) True or False depending on whether sequence is hiv or not
    '''
    tmp_dir = tempfile.gettempdir()
    tmpfile = os.path.join(tmp_dir, "tmp_blast_results.xml")

    # blast_results = NCBIWWW.qblast('blastn', 'nr', query_file)
    blastn_cline = NcbiblastnCommandline(query=query_file, db="nt", evalue=0.001, outfmt=5, perc_identity=50,
                                         out=tmpfile, max_target_seqs=3)
    stdout, stderr = blastn_cline()
    print(stderr, "\n", stdout)
    # with open(tmpfile, 'w') as handle:
    #     handle.write(blast_results.read())

    e_value_threshold = 0.01

    for blast_record in NCBIXML.parse(open(tmpfile)):

        # skip empty records
        if blast_record.alignments:

            # take top hit
            top_hit = blast_record.alignments[0]

            # get hit title
            title = top_hit.title.split(" ")
            e_val = top_hit.hsps[0].expect
            if e_val < e_value_threshold:

                # is the hit to HIV-1?
                if "HIV-1" in title:
                    print("HIV")
        #            return False
                else:
                    print("not hiv", title)
        #         else:
        #             with open(logfile, 'a') as handle:
        #                 handle.write(seqname + "\ttop blastn hit:\t" + "_".join(title) + "\n")
        #             return True
        #
        # else:
        #     print("No blastn hits")
        #     with open(logfile, 'a') as handle:
        #         handle.write(seqname + "\ttop blastn hit:\t" + "None" + "\n")
        #     return True


def regex_contam(query_sequence, hxb2_region):
    '''
    :param query_sequence: (str) the sequence to blast
    :param hxb2_region: (str_ hxb2 sequence for the relevant gene region
    :return: (bool) True or False depending on whether sequence is hiv or not
    '''
    # todo hardcode hxb2 regex string for GAG_2? for now
    # hxb2_region = ''
    error = int(len(query_sequence) * 0.6)
    hxb2_regex = "r'({0}){{{1}}}'".format(hxb2_region, error)
    query = str(query_sequence)
    match = regex.search(hxb2_regex, query, regex.BESTMATCH)
    if match is not None:
        return False
    else:
        print("no match")
        return True


def main(consensus, outpath, logfile):

    # initialize file names
    cln_cons_name = os.path.split(consensus)[-1]
    cln_cons = cln_cons_name.replace(".fasta", "_good.fasta")
    consensus_out = os.path.join(outpath, cln_cons)
    contam_seqs = cln_cons_name.replace(".fasta", "_contam_seqs.fasta")
    contam_out = os.path.join(outpath, contam_seqs)

    with open(logfile, 'a') as handle:
        handle.write("Contam removal step:\nSequences identified as contaminants (if any):")

    # clear out any existing outfiles for consensus
    with open(consensus_out, 'w') as handle1:
        handle1.write("")
    with open(contam_out, 'w') as handle2:
        handle2.write("")

    # get hxb2 seq for regex??
    # gag1
    hxb2_seq = "GCGAAAAATTAGATAATTGGGAAAGAATTAAGTTAAGGCCAGGAGGAAAGAAACACTATATGCTAAAAC"
    # gag2
    hxb2_seq = "ACCAAATGAAAGACTGTACTGAGAGGCAGGCTAATTTTTTAGGGAAAATTTGGCCTTCCTACAAGGGGA"

    # store all consensus seqs in a dict
    all_sequences_d = fasta_to_dct(consensus)
    is_contam = blastn_seqs(consensus, logfile)

    # check each sequence to see if it is a contaminating sequence
    # for name, seq in all_sequences_d.items():
    #     # is_contam = regex_contam(seq, hxb2_seq)
    #     is_contam = blastn_seqs(seq, name, logfile)
    #     # if is contam, save to contam file
    #     if is_contam:
    #         print("Non HIV sequence found:\n", name, "\n", seq)
    #         with open(contam_out, 'a') as handle1:
    #             outstr = "{0}{1}\n{2}\n".format('>', name, seq)
    #             handle1.write(outstr)
    #
    #     # if is not contam, to write to good outfile
    #     else:
    #         print("HIV")
    #         with open(consensus_out, 'a') as handle2:
    #             outstr = "{0}{1}\n{2}\n".format('>', name, seq)
    #             handle2.write(outstr)
    #     input("enter")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--consensus', default=argparse.SUPPRESS, type=str,
                        help='The read1 (R1) fastq file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be copied', required=True)
    parser.add_argument('-l', '--logfile', default=argparse.SUPPRESS, type=str,
                        help='The path and name of the log file', required=True)

    args = parser.parse_args()
    consensus = args.consensus
    outpath = args.outpath
    logfile = args.logfile

    main(consensus, outpath, logfile)
