#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import os
import argparse
import collections
from itertools import groupby
# from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
import sys


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


def blastn_seqs(infile, gene_region, outpath):
    """
    :param infile: fasta file to blast
    :param gene_region: (str) the target gene region
    :param outpath: path to the outfile
    :return: (bool) True or False depending on whether sequence is not hiv
    """

    # assign temp file
    tmp_out_file = os.path.join(outpath, "tmp_blast.xml")
    target_gene = [gene_region.upper().split("_")[0]]

    # allow blast hit to pass if blasts to overlapping gene
    if "VIF" in target_gene:
        target_gene = ["VPR", "VIF"]
    if "POL" in target_gene:
        try:
            check_for_overlap = gene_region.upper().split("_")[1]
            if check_for_overlap.upper().split("_")[1] == "5":
                target_gene = ["INT", "VIF"]
        except IndexError as e:
            print(e, "\nno overlapping gene regions detected")
            pass

    # blast settings
    # format_fasta = ">{0}\n{1}".format(s_name, q_sequence)
    blastdb = "lanl_hiv_db"
    outformat = 5
    e_value = 0.000001
    threads = 4

    print("running blast")
    # # run online blast
    # blast_results = NCBIWWW.qblast('blastn', 'nt', query=infile, entrez_query='"HIV-1"[organism]')
    # # write online blast results to file
    # with open(tmp_out_file, 'w') as handle:
    #     handle.write(blast_results.read())

    # run local blast
    blastn_cline = NcbiblastnCommandline(query=infile, db=blastdb, evalue=e_value, outfmt=outformat, perc_identity=80,
                                         out=tmp_out_file, num_threads=threads)

    stdout, stderr = blastn_cline() # stdin=format_fasta
    print("stderr = ", stderr)
    print("stout = ", stdout)

    good_records = collections.defaultdict(list)
    bad_records = collections.defaultdict(list)

    if os.path.isfile(tmp_out_file):
        all_blast_results = NCBIXML.parse(open(tmp_out_file))
    else:
        print("the blast failed to write an outfile")
        sys.exit()

    for blast_record in all_blast_results:
        # get query name
        query_seq_name = blast_record.query.upper()

        # was there a hit to something in the db?
        if blast_record.alignments:
            found = False
            first_region = ''
            for i, alignment in enumerate(blast_record.alignments):
                title_name = alignment.title.split(" ")[0]
                region = title_name.upper().split("_")[-1]
                if i == 0:
                    first_region = title_name.upper().split("_")[-1]

                if region in target_gene:
                    found = True
                    good_records[query_seq_name] = "_hiv_" + region.lower()
                    break
            # todo: should this be indented?
            if not found:
                bad_records[query_seq_name] = "_hiv_" + first_region.lower()

        else:
            # no hit in db
            bad_records[query_seq_name] = "_not_hiv_" + "no_hit"

    os.unlink(tmp_out_file)

    return bad_records, good_records


def main(consensus, outpath, gene_region, logfile):
    print(consensus)
    # initialize file names
    cln_cons_name = os.path.split(consensus)[-1]
    cln_cons = cln_cons_name.replace("_clean.fasta", "_good.fasta")
    consensus_out = os.path.join(outpath, cln_cons)
    contam_seqs = cln_cons_name.replace("_clean.fasta", "_contam_seqs.fasta")
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
    contam, not_contam = blastn_seqs(consensus, gene_region, outpath)
    # set all output names to uppercase to ensure input > output names match
    contam_names = [x.upper() for x in contam.keys()]
    not_contam_names = [x.upper() for x in not_contam.keys()]

    # check each sequence to see if it is a contaminant
    for name, seq in all_sequences_d.items():
        # set input name to uppercase to ensure match to output name
        name = name.upper()
        # if the sequence is not hiv, save to contam file
        if name in contam_names:
            print("Non HIV sequence found:\n\t", name)
            new_name = name + contam[name]
            with open(contam_out, 'a') as handle1:
                outstr = ">{0}\n{1}\n".format(new_name, seq)
                handle1.write(outstr)

        # if is not contam, to write to good outfile
        elif name in not_contam_names:
            with open(consensus_out, 'a') as handle2:
                outstr = ">{0}\n{1}\n".format(name, seq)
                handle2.write(outstr)

        else:
            print("Input names did not match with output names. Something strange happened for: ", name)
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
