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
import sys

__author__ = 'Colin Anthony'


def fasta_to_dct(fn):
    '''
    converts a fasta file to a dictionary where key = seq name, value = sequence
    :param fn: a fasta file
    :return: a dictionary
    '''
    dct = collections.defaultdict(list)
    for seq_record in SeqIO.parse(open(fn), "fasta"):
        dct[seq_record.description.replace(" ", "_").upper()] = str(seq_record.seq).replace("-", "").upper()
        # dct[str(seq_record.seq).replace("-", "").upper()].append(seq_record.description.replace(" ", "_").upper())
    return dct


def blastn_seqs(infile, gene_region):
    '''
    :param name: sequence name
    :param sequence: (str) the sequence to blast
    :param logfile: the logfile path and name
    :return: (bool) True or False depending on whether sequence is not hiv
    '''

    # assign temp file
    tmp_dir = tempfile.gettempdir()
    tmp_file = os.path.join(tmp_dir, "tmp_blast_results.xml")

    # clear the tmp_file
    with open(tmp_file, 'w') as handle:
        handle.write("")

    target_gene = gene_region.upper().split("_")[0]
    # blast settings
    # format_fasta = ">{0}\n{1}".format(s_name, q_sequence)
    blastdb = "lanl_hiv_db"
    outformat = 5
    e_value = 0.0001
    threads = 4
    max_hits = 1

    print("running blast")
    ## run online blast
    # blast_results = NCBIWWW.qblast('blastn', 'nt', query=infile, entrez_query='"HIV-1"[organism]')
    ## write online blast results to file
    # with open(tmpfile, 'w') as handle:
    #     handle.write(blast_results.read())

    # run local blast
    blastn_cline = NcbiblastnCommandline(query=infile, db=blastdb, evalue=e_value, outfmt=outformat, perc_identity=50,
                                         out=tmp_file, num_threads=threads,  max_target_seqs=max_hits)

    stdout, stderr = blastn_cline() # stdin=format_fasta
    print("stderr", stderr)
    print("stout", stdout)

    good_records = collections.defaultdict(list)
    bad_records = collections.defaultdict(list)

    if os.path.isfile(tmp_file):
        all_blast_results = NCBIXML.parse(open(tmp_file))
    else:
        print("the blast failed to write an outfile")
        sys.exit()

    for blast_record in all_blast_results:
        # get query name
        query_seq_name = blast_record.query

        # was there a hit to something in the db?
        if blast_record.alignments:
            for alignment in blast_record.alignments:
                    title_name = alignment.title.split(" ")[0]
                    region = title_name.upper().split("_")[-1]
                    for hsp in alignment.hsps:
                        # get the e_value in case you want to store it
                        exp_value = hsp.expect

                    if region == target_gene:
                        good_records[query_seq_name] = "_hiv_" + region
                    else:
                        bad_records[query_seq_name] = "_hiv_" + region

        else:
            # no hit in db
            bad_records[query_seq_name] = "_not_hiv_" + "no_hit"

    return bad_records, good_records


def regex_contam(name, query_sequence):
    '''
    :param query_sequence: (str) the sequence to blast
    :param hxb2_region: (str_ hxb2 sequence for the relevant gene region
    :return: (bool) True or False depending on whether sequence is hiv or not
    '''

    # gag1
    # hxb2_seq = "GCGAAAAATTAGATAATTGGGAAAGAATTAAGTTAAGGCCAGGAGGAAAGAAACACTATATGCTAAAAC"
    # gag2
    # hxb2_seq = "ACCAAATGAAAGACTGTACTGAGAGACAGGCTAATTTTTTAGGGAAAATTTGGCCTTCCCACAAGGGGAGGCCAGGGAATTTCC"

    # error = int(1 - (len(query_sequence) * 0.7))
    # hxb2_regex = "r'({0}){{e<{1}}}'".format(hxb2_region, error)
    hxb2_regex = r'(ACCAAATGAAAGACTGTACTGAGAGACAGGCTAATTTTTTAGGGAAAATTTGGCCTTCCCACAAGGGGAGGCCAGGGAATTTCC){e<8}'
    match = regex.search(hxb2_regex, query_sequence, regex.BESTMATCH)

    if match is not None:
        return False
    else:
        return True


def main(consensus, outpath, gene_region, logfile):
    print(consensus)
    # initialize file names
    cln_cons_name = os.path.split(consensus)[-1]
    cln_cons = cln_cons_name.replace("_clean.fasta", "_good.fasta")
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

    # store all consensus seqs in a dict
    all_sequences_d = fasta_to_dct(consensus)

    # checck for contam
    contam, not_contam = blastn_seqs(consensus, gene_region)

    # check each sequence to see if it is a contaminant
    for name, seq in all_sequences_d.items():

        # if the sequence is not hiv, save to contam file
        if name in contam.keys():
            print("Non HIV sequence found:\n\t", name)
            new_name = name + contam[name]
            with open(contam_out, 'a') as handle1:
                outstr = ">{0}\n{1}\n".format(new_name, seq)
                handle1.write(outstr)

        # if is not contam, to write to good outfile
        elif name in not_contam.keys():
            with open(consensus_out, 'a') as handle2:
                outstr = ">{0}\n{1}\n".format(name, seq)
                handle2.write(outstr)

    else:
        print("not contam, not not contam. Something strange happened")
        # todo add logging

    print("contam check complete")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--consensus', default=argparse.SUPPRESS, type=str,
                        help='The read1 (R1) fastq file', required=True)
    parser.add_argument('-o', '--outpath', default=argparse.SUPPRESS, type=str,
                        help='The path to where the output file will be copied', required=True)
    parser.add_argument('-l', '--logfile', default=argparse.SUPPRESS, type=str,
                        help='The path and name of the log file', required=True)
    parser.add_argument('-g', '--gene_region', default=argparse.SUPPRESS, type=str,
                        help='the genomic region being sequenced, '
                             'valid options: GAG_1/GAG_2/ENV_C1C2/POL_1/NEF_1 etc..', required=True)

    args = parser.parse_args()
    consensus = args.consensus
    outpath = args.outpath
    gene_region = args.gene_region
    logfile = args.logfile

    main(consensus, outpath, gene_region, logfile)
