# !/usr/bin/python3
import csv
from subprocess import call
from subprocess import Popen
from subprocess import PIPE
import argparse
import os
import ntpath
import json
import logging
import time
from glob import glob

# Non-standard
import regex


__author__ = "Colin Anthony, Jon Ambler, David Matten"
__copyright__ = "Something"
__credits__ = ["Colin Anthony", "Jon Ambler", "David Matten"]
__license__ = "TBD"
__version__ = "0.2"
__maintainer__ = "Colin Anthony"
__email__ = ""
__status__ = "Testing"


def is_divisable(dividend, divisor):
    """
    Check if a number is divisible by another without remainders.
    :param dividend: https://en.wikipedia.org/wiki/Division_(mathematics)
    :param divisor: https://en.wikipedia.org/wiki/Division_(mathematics)
    :return: True or false.
    """
    if dividend % divisor == 0:
        return True
    else:
        return False


def make_seq_wild(input_sequence):
    """
    Replaces non-standard nucleotides with a '?' character for use in the regex function.
    :param input_sequence: A nucleotide sequence as a string.
    :return: The nucleotide sequence with replacements as a string.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    wild_seq = ''

    if len(input_sequence) <= 0:
        logging.warning("0 length sequence submitted to the make_seq_wild function.")

    for nuc in input_sequence:
        if nuc not in nucleotides:
            wild_seq += '?'
        else:
            wild_seq += nuc
    return wild_seq


def break_into_all_kmers(aString):
    """
    returns a list containing all possible substrings of all possible lengths for a given string.
    :param aString: The input string
    :return: A list of strings
    """
    length = len(aString)
    return [aString[i:j + 1] for i in range(length) for j in range(i, length)]


def make_kmer_key_list(list_of_strings, min_len=0, max_len=10000000, max_kmers=-1):
    """
    A kmer key list is a list of substrings of a given string that may be used to uniquely identify the string.
    Each substring is represented only once in the list of keys. Here, each string represents a primer, and the primer
    is mapped back to the primer dictionary based on its sequence.
    This function is called by add_kmer_keys.
    Why pass a list of all the primers? This is so the k-mers are unique to each primer,
    :param list_of_strings:
    :param min_len: The shortest length of sequence allowed for the k-mers
    :param max_len: The longest length of sequence allowed for the k-mers
    :param max_kmers: Limit the number of k-mers to be included. Longer k-mers are kept and shortest are discarded.
    :return: a dictionary where the keys are the string / primer and the values are a list of unique k-mers for that
    primer.
    """
    string_dict = {}
    unique_dict = {}
    for a_string in list_of_strings:
        sub_string_list = break_into_all_kmers(a_string)
        string_dict[a_string] = sub_string_list
        # Initialise a dictionary with an empty list
        unique_dict[a_string] = []

    for a_string in list_of_strings:
        for a_kmer in string_dict[a_string]:
            unique = True
            for another_string in list(string_dict.keys()):
                if another_string != a_string:
                    if a_kmer in string_dict[another_string]:
                        unique = False
            if unique and min_len <= len(a_kmer) <= max_len:
                unique_dict[a_string].append(a_kmer)

    if max_kmers > 0:
        shorter_unique_dict = {}
        for a_string in list(unique_dict.keys()):
            shorter_unique_dict[a_string] = []
            count = 0
            while count < max_kmers and len(unique_dict[a_string]) > 0:
                largest_kmer = max(unique_dict[a_string], key=len)
                shorter_unique_dict[a_string].append(largest_kmer)
                unique_dict[a_string].remove(largest_kmer)
                count += 1

            if count == 0:
                # Send error if no k-mers are found for the primer
                logging.error("No unique k-mers found for primer " + a_string)

        unique_dict = shorter_unique_dict

    return unique_dict


def make_primer_dict(primer_file):
    """
    This function creates the primer info dict containing information including the primer name, sequence,
    and a version of the primer with the non-standard nucleotide characters replaced to allow the use of
    the python regex package.
    :param primer_file:
    :param forward_index_len: This is the length of sequence that comes before the primer sequence in the read
    (in the case of the forward sequence)
    :param reverse_index_len: This is the length of sequence that comes before the primer sequence in the read
    (in the case of the reverse sequence)
    :param forward_init_seq:
    :param reverse_init_seq:
    :return: A dictionary containing the relevant information for the primers.
    """
    res_dict = {}
    with open(primer_file, 'r') as csvfile:
        primer_file_obj = csv.reader(csvfile, delimiter=',')
        # Skip header
        next(primer_file_obj, None)
        for row in primer_file_obj:
            sub_region = row[1].replace(' ', '')
            if sub_region == "None" or sub_region == None:
                sub_region = False

            res_dict[row[0]] = {
                'sub_region': sub_region,
                'fwd': row[2][int(row[3]):].replace(' ', ''),
                'rev': row[4][int(row[5]):].replace(' ', ''),
                'fwd_full': row[2].replace(' ', ''),
                'rev_full': row[4].replace(' ', ''),
                'overlap': row[6],
                'fwd_preseq': int(row[3]),
                'rev_preseq': int(row[5]),
                'fwd_wild': make_seq_wild(row[2][int(row[3]):].replace(' ', '')),
                'rev_wild': make_seq_wild(row[4][int(row[5]):].replace(' ', '')),
                'exact_matches_found': 0,
                'regex_matches_found': 0,
                'kmer_matches_found': 0,
                'blast_matches_found': 0,
            }

    number_of_regions = len(res_dict.keys())
    logging.info("Primer file parsed: " + str(number_of_regions) + " records found")

    return res_dict


def export_line_to_fastq(fileTag, headerLine, seqLine, plusLine, scoreLine, infast_name, out_dir, patient_list):
    """
    For the export of a fastq sequence, line by line.
    :param fileTag: the gene region???
    :param headerLine:
    :param seqLine:
    :param plusLine:
    :param scoreLine:
    :return:
    """
    out_file_name = infast_name.replace('multiplex', fileTag)
    out_file_path = out_dir + patient_list + '/' + fileTag + '/0new_data/'
    outfile = os.path.join(out_file_path, out_file_name)
    with open(outfile, "a+") as handle:
        handle.write(headerLine)
        handle.write(seqLine)
        handle.write(plusLine)
        handle.write(scoreLine)

    return True


def wildcard_seq_match(primer, sequence, error_rate):
    """
    Check if the primer sequence is found in the sequence, allowing for wildcards.
    :param primer: The primer to test for a match
    :param sequence: The sequence that will be compared to the primer.
    :param error_rate: The number of changes allowed between the sequences
    :return: True or false depending on if there is a match.
    """

    primer_pattern = r'(' + primer + '){e<' + str(error_rate) + '}'

    match = regex.search(primer_pattern, sequence, regex.BESTMATCH)

    if match is not None:
        return True
    else:
        return False


def add_kmer_keys(primerDict):
    """
    Add the primer keys to an existing primer dictionary.
    :param primerDict: A primer dict created by the make_primer_dict function.
    :return: a primer dictionary object with the unique k-mer keys added.
    """

    fwd_primer_list = []

    for a_forward_primer in list(primerDict.keys()):
        fwd_primer_list.append(primerDict[a_forward_primer]['fwd'])

    fwd_kmer_key_list = make_kmer_key_list(fwd_primer_list, min_len=10, max_len=20, max_kmers=30)

    for primer in list(primerDict.keys()):
        for kmer_seq in list(fwd_kmer_key_list.keys()):
            if primerDict[primer]['fwd'] == kmer_seq:
                primerDict[primer]['fwd_keys'] = fwd_kmer_key_list[kmer_seq]

    rev_primer_list = []

    for a_reverse_primer in list(primerDict.keys()):
        rev_primer_list.append(primerDict[a_reverse_primer]['rev'])

    rev_kmer_key_list = make_kmer_key_list(rev_primer_list, min_len=10, max_len=20, max_kmers=30)

    for primer in list(primerDict.keys()):
        for kmer_seq in list(rev_kmer_key_list.keys()):
            if primerDict[primer]['rev'] == kmer_seq:
                primerDict[primer]['rev_keys'] = rev_kmer_key_list[kmer_seq]

    return primerDict


def split_by_primers(fastq_file, primer_dict, orientation, infast_name, out_dir, patient_list, regex_error_rate,
                     make_sure=False):
    """
    Take an input fastq file and split it into individual fastq files, with the split based on the presence of
    a primer sequence specified in a dictionary.
    :param fastq_file:
    :param primer_dict:
    :param orientation: Whether these are forward (R1) or reverse (R2) reads. Options are 'fwd' or 'rev'.
    :param regex_error_rate: error rate for regex matching.
    :param make_sure: When set to True, a blast search is included in the process.
    :return:
    """
    line_number_quality = 1
    line_number_plus = 2
    line_number_sequence = 3
    line_number_header = 4

    # Get the shortest primer length
    shortest_primer_length = 1000000
    for aPrimer in list(primer_dict.keys()):
        if len(primer_dict[aPrimer][orientation]) < shortest_primer_length:
            shortest_primer_length = len(primer_dict[aPrimer][orientation])

    if make_sure:
        # todo this should be moved to a function of its own
        # Make a fasta file with all primers
        temp_fasta = open(out_dir + 'primerList.fasta', 'w')

        for geneRegion in primer_dict:
            header_line = '>' + geneRegion + '\n'
            seq_line = primer_dict[geneRegion][orientation] + '\n'
            temp_fasta.write(header_line)
            temp_fasta.write(seq_line)
        temp_fasta.close()

        # Make a blast database
        create_temp_blast_db(out_dir + 'primerList.fasta', out_dir + 'primers')

    for seq_line in open(fastq_file, 'r'):
        # Get the header
        if is_divisable(line_number_header, 4):
            current_header = seq_line

        # Get the sequence and match to a primer
        if is_divisable(line_number_sequence, 4):
            detected_primer = 'None'
            current_sequence = seq_line
            for geneRegion in list(primer_dict.keys()):
                # The geneRegion is the gene target of the primer
                seq_primer_region = seq_line[primer_dict[geneRegion][orientation + '_preseq']:
                                             len(primer_dict[geneRegion][orientation]) +
                                             primer_dict[geneRegion][orientation + '_preseq']]

                # Level 1: Check exact matches
                if primer_dict[geneRegion][orientation] == seq_primer_region:
                    detected_primer = geneRegion
                    primer_dict[geneRegion]['exact_matches_found'] += 1

                # Level 2: If no exact match, then look for regex matches
                elif wildcard_seq_match(primer_dict[geneRegion][orientation + '_wild'], seq_primer_region,
                                        regex_error_rate):
                    detected_primer = geneRegion
                    primer_dict[geneRegion]['regex_matches_found'] += 1

                # Level 3: If no wildmatch, look for kmer matches
                else:
                    for kmer in primer_dict[geneRegion][orientation + '_keys']:
                        if kmer in seq_primer_region:
                            detected_primer = geneRegion
                            primer_dict[geneRegion]['kmer_matches_found'] += 1

            # Finally, if there is still no match, use the blastDB
            if detected_primer == 'None' and make_sure:

                detected_primer = primer_blast_search('primers', current_sequence[:shortest_primer_length], out_dir)
                if detected_primer != 'None':
                    primer_dict[geneRegion]['blast_matches_found'] += 1

        # Get the plus sign
        if is_divisable(line_number_plus, 4):
            current_plus = seq_line

        # Get the quality scores
        if is_divisable(line_number_quality, 4):
            current_quality = seq_line

        # At the end of the record, write to the appropriate output file
        if is_divisable(line_number_quality, 4):
            export_line_to_fastq(detected_primer, current_header, current_sequence, current_plus, current_quality,
                                 infast_name, out_dir, patient_list)

        line_number_quality += 1
        line_number_plus += 1
        line_number_sequence += 1
        line_number_header += 1

    splitReport = open(orientation + '_splitReport.csv', 'w')
    splitReport.write('Region,Exact matches,Regex matches,Kmer matches,Blast matches\n')
    for a_gene_region in primer_dict.keys():
        out_line = a_gene_region + ',' + \
                   str(primer_dict[a_gene_region]['exact_matches_found']) + ',' + \
                   str(primer_dict[a_gene_region]['regex_matches_found']) + ',' + \
                   str(primer_dict[a_gene_region]['kmer_matches_found']) + ',' + \
                   str(primer_dict[a_gene_region]['blast_matches_found']) + '\n'
        splitReport.write(out_line)

        logging.info('Exact matches found: ' + str(primer_dict[a_gene_region]['exact_matches_found']))
        logging.info('Regex matches found: ' + str(primer_dict[a_gene_region]['regex_matches_found']))
        logging.info('K-mer matches found: ' + str(primer_dict[a_gene_region]['kmer_matches_found']))
        logging.info('Blast matches found: ' + str(primer_dict[a_gene_region]['blast_matches_found']))

        # Check to see if no matches were found
        if primer_dict[a_gene_region]['exact_matches_found'] == 0:
            logging.warning('No exact matches found')
        if primer_dict[a_gene_region]['regex_matches_found'] == 0:
            logging.warning('No regex matches found')
        if primer_dict[a_gene_region]['kmer_matches_found'] == 0:
            logging.warning('No k-mer matches found')

    splitReport.close()


def split_by_primers_matchpair(fastq_R1_file, fastq_R2_file, primer_dict, orientation, infast_R1_name, infast_R2_name,
                               out_dir, patient_list, regex_error_rate, make_sure=False):
    """
    Take an input fastq file and split it into individual fastq files, with the split based on the presence of
    a primer sequence specified in a dictionary. This version does matching on the R1 reads and finds the matching
    R2 read based on position in the fastq file.
    :param fastq_R1_file:
    :param fastq_R2_file:
    :param primer_dict:
    :param orientation: Whether these are forward (R1) or reverse (R2) reads. Options are 'fwd' or 'rev'.
    :param regex_error_rate: error rate for regex matching.
    :param make_sure: When set to True, a blast search is included in the process.
    :return:
    """
    line_number_quality = 1
    line_number_plus = 2
    line_number_sequence = 3
    line_number_header = 4

    # Get the shortest primer length
    shortest_primer_length = 1000000
    for aPrimer in list(primer_dict.keys()):
        if len(primer_dict[aPrimer][orientation]) < shortest_primer_length:
            shortest_primer_length = len(primer_dict[aPrimer][orientation])

    if make_sure:
        # todo this should be moved to a function of its own
        # Make a fasta file with all primers
        temp_fasta = open(out_dir + 'primerList.fasta', 'w')

        for geneRegion in primer_dict:
            header_line = '>' + geneRegion + '\n'
            seq_line = primer_dict[geneRegion][orientation] + '\n'
            temp_fasta.write(header_line)
            temp_fasta.write(seq_line)
        temp_fasta.close()

        # Make a blast database
        create_temp_blast_db(out_dir + 'primerList.fasta', out_dir + 'primers')

    for seq_line, R2_seq_line in zip(open(fastq_R1_file, 'r'), open(fastq_R2_file, 'r')):
        # Get the header
        if is_divisable(line_number_header, 4):
            current_header = seq_line
            current_header_R2 = R2_seq_line

        # Get the sequence and match to a primer
        if is_divisable(line_number_sequence, 4):
            detected_primer = 'None'
            current_sequence = seq_line
            current_sequence_R2 = R2_seq_line
            for geneRegion in list(primer_dict.keys()):
                # The geneRegion is the gene target of the primer
                seq_primer_region = seq_line[primer_dict[geneRegion][orientation + '_preseq']:
                                             len(primer_dict[geneRegion][orientation]) +
                                             primer_dict[geneRegion][orientation + '_preseq']]

                # Level 1: Check exact matches
                if primer_dict[geneRegion][orientation] == seq_primer_region:
                    detected_primer = geneRegion
                    primer_dict[geneRegion]['exact_matches_found'] += 1

                # Level 2: If no exact match, then look for regex matches
                elif wildcard_seq_match(primer_dict[geneRegion][orientation + '_wild'], seq_primer_region,
                                        regex_error_rate):
                    detected_primer = geneRegion
                    primer_dict[geneRegion]['regex_matches_found'] += 1

                # Level 3: If no wildmatch, look for kmer matches
                else:
                    for kmer in primer_dict[geneRegion][orientation + '_keys']:
                        if kmer in seq_primer_region:
                            detected_primer = geneRegion
                            primer_dict[geneRegion]['kmer_matches_found'] += 1

            # Finally, if there is still no match, use the blastDB
            if detected_primer == 'None' and make_sure:

                detected_primer = primer_blast_search('primers', current_sequence[:shortest_primer_length], out_dir)
                if detected_primer != 'None':
                    primer_dict[geneRegion]['blast_matches_found'] += 1

        # Get the plus sign
        if is_divisable(line_number_plus, 4):
            current_plus = seq_line
            current_plus_R2 = R2_seq_line

        # Get the quality scores
        if is_divisable(line_number_quality, 4):
            current_quality = seq_line
            current_quality_R2 = R2_seq_line

        # At the end of the record, write to the appropriate output file for the R1 file
        if is_divisable(line_number_quality, 4):
            export_line_to_fastq(detected_primer, current_header, current_sequence, current_plus, current_quality,
                                 infast_R1_name, out_dir, patient_list)

            export_line_to_fastq(detected_primer, current_header_R2, current_sequence_R2, current_plus_R2,
                                 current_quality_R2,
                                 infast_R2_name, out_dir, patient_list)

        line_number_quality += 1
        line_number_plus += 1
        line_number_sequence += 1
        line_number_header += 1

    splitReport = open(orientation + '_splitReport.csv', 'w')
    splitReport.write('Region,Exact matches,Regex matches,Kmer matches,Blast matches\n')
    for a_gene_region in primer_dict.keys():
        out_line = a_gene_region + ',' + \
                   str(primer_dict[a_gene_region]['exact_matches_found']) + ',' + \
                   str(primer_dict[a_gene_region]['regex_matches_found']) + ',' + \
                   str(primer_dict[a_gene_region]['kmer_matches_found']) + ',' + \
                   str(primer_dict[a_gene_region]['blast_matches_found']) + '\n'
        splitReport.write(out_line)

        logging.info('Exact matches found: ' + str(primer_dict[a_gene_region]['exact_matches_found']))
        logging.info('Regex matches found: ' + str(primer_dict[a_gene_region]['regex_matches_found']))
        logging.info('K-mer matches found: ' + str(primer_dict[a_gene_region]['kmer_matches_found']))
        logging.info('Blast matches found: ' + str(primer_dict[a_gene_region]['blast_matches_found']))

        # Check to see if no matches were found
        if primer_dict[a_gene_region]['exact_matches_found'] == 0:
            logging.warning('No exact matches found')
        if primer_dict[a_gene_region]['regex_matches_found'] == 0:
            logging.warning('No regex matches found')
        if primer_dict[a_gene_region]['kmer_matches_found'] == 0:
            logging.warning('No k-mer matches found')

    splitReport.close()


def create_temp_blast_db(fasta_filepath, db_identifier):
    """
    A command line call to create a local blast database from the given set of sequences.
    :param fasta_filepath:
    :param db_identifier:
    :return:
    """
    # Create a blast database for a genome
    database_name = db_identifier + "_BlastDB"
    return call(["makeblastdb", "-in", fasta_filepath, "-out", database_name, "-dbtype", "nucl"])


def primer_blast_search(db_identifier, sequence, out_dir):
    """
    This function is used to search an existing blast database for a sequence and is a wrapper for blastn.
    :param db_identifier: The identifier for the blast database.
    :param sequence: The query sequence as a string.
    :return: Returns the name of the sequence that is the best hit. In this case, the name of the primer.
    """
    # Search the selected database for the sequence
    # The input sequence must be in fasta format... which may prove a problem
    # Return a ordered result as a list of dictionaries... to be decided
    blast_db_directory = out_dir
    blast_database = blast_db_directory + db_identifier + "_BlastDB"

    # Make a fasta file with the sequence in it
    temp_dir = out_dir
    temp_fasta = open(temp_dir + 'temp.fa', 'w')
    temp_fasta.write(">temp_fasta_file\n")
    temp_fasta.write(sequence + "\n")
    temp_fasta.close()

    # Doing the blast
    #  blastn -db F11_BlastDB -query /Users/Admin/Dropbox/Programs/tools/Cell/Cell_core/blasta.fa
    # print "blasting"
    full_directory = temp_dir + "temp.fa"
    blast_result_raw = Popen(["blastn", "-db", blast_database, "-query", full_directory, "-outfmt", "7", "-word_size",
                              "15", "-evalue", "1000"], stdout=PIPE)
    blast_result = blast_result_raw.stdout.read()
    blast_result = blast_result.decode('ascii')

    # print len(blast_result)

    # Formatting into a list
    new_list = blast_result.split('\n')

    # Remove empty list items
    new_list = [_f for _f in new_list if _f]

    # make into dict & add to list
    List_of_result_lists = []
    for entry in new_list:
        if entry[0] != "#":
            split_entry = entry.split('\t')
            List_of_result_lists.append(split_entry)

    # Lists are in the form:
    # query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end,
    # evalue, bit score

    if len(List_of_result_lists) > 0:
        return List_of_result_lists[0][1]
    else:
        return "None"


def process_input_dir(fastq_directory_path):
    """
    This function finds all the fastq files in a given directory, and groups the based on the patient, sample,
    and the read pair. So R1 and R2 are grouped together, and samples belonging to the same patient are grouped.
    :param fastq_directory_path: A string that is the path to the directory containing the fastq files.
    :return: a dict
    """
    search_path = os.path.join(fastq_directory_path, "*.fastq")
    file_list = glob(search_path)

    infile_dict = {}

    for a_file in file_list:
        a_file = os.path.abspath(a_file)
        file_parts = a_file.split("_")
        if "None" in file_parts:
            continue
        if "R1" in a_file:
            path = os.path.split(a_file)[0]
            inf_R1_name = os.path.split(a_file)[-1]
            outf_R1 = inf_R1_name.replace("-", "_")
            outf_R1_rename = regex.sub("S[0-9]+_L[0-9][0-9][0-9]_R1_[0-9][0-9][0-9].fastq", "R1.fastq", outf_R1)
            outf_R1_rename_with_path = os.path.join(path, outf_R1_rename)

            if outf_R1_rename == outf_R1:
                # check if they have already been renamed
                if outf_R1.split("_")[-1] == "R1.fastq":
                    print("file was already in correct format?")
                else:
                    print()
                    raise ValueError("Unable to rename R1 file {0}\ncheck the file renaming regex".format(a_file))
            else:
                os.rename(a_file, outf_R1_rename_with_path)

            print(outf_R1_rename)

            orientation = "R1"
            patient = outf_R1_rename.split('_')[0]
            identifier = outf_R1_rename.replace("_R1.fastq", "")

        elif "R2" in a_file:
            path = os.path.split(a_file)[0]
            inf_R2_name = os.path.split(a_file)[-1]
            outf_R2 = inf_R2_name.replace("-", "_")
            outf_R2_rename = regex.sub("S[0-9]+_L[0-9][0-9][0-9]_R2_[0-9][0-9][0-9].fastq", "R2.fastq", outf_R2)
            outf_R2_rename_with_path = os.path.join(path, outf_R2_rename)
            if outf_R2_rename == outf_R2:
                if outf_R2.split("_")[-1] == "R2.fastq":
                    print("file was already in correct format?")
                else:
                    raise ValueError("Unable to rename R1 file {0}\ncheck the file renaming regex".format(a_file))
            else:
                os.rename(a_file, outf_R2_rename_with_path)

            orientation = "R2"
            patient = outf_R2_rename.split('_')[0]
            identifier = outf_R2_rename.replace("_R2.fastq", "")

        else:
            continue
        # if a_file[-5:] == 'fastq' or a_file[-2:] == 'fq':
        #
        #     # Identify patient tag
        #     patient = a_file.split('_')[0]
        #
        #     # Identify if R1 or R2
        #     orientation = a_file.split('_')[-3]
        #
        #     # Get filename without R1/2
        #     identifier = a_file.split('_')[:-4] + a_file.split('_')[-2:]
        #     identifier = ''.join(identifier)

        if patient not in infile_dict.keys():
            infile_dict[patient] = {identifier: {orientation: a_file}}

        elif identifier not in infile_dict[patient].keys():
            infile_dict[patient][identifier] = {orientation: a_file}

        elif orientation not in infile_dict[patient][identifier].keys():
            infile_dict[patient][identifier][orientation] = a_file

    logging.info(infile_dict)

    return infile_dict


def check_for_previous_runs(output_directory_to_check, parsed_input_file_dict):

    file_list = os.listdir(output_directory_to_check)

    # if output_directory_to_check[-1] != '/':
    #     output_directory_to_check = output_directory_to_check + '/'
    print(output_directory_to_check)
    output_directory_to_check = os.path.abspath(output_directory_to_check)
    print(output_directory_to_check)
    setup_log = os.path.join(output_directory_to_check, 'setup_log.csv')
    previous_run_dict = {}

    if not os.path.isfile(setup_log):

        logging.info('No previous runs detected.')

        print('Running first time setup')

        with open(setup_log, 'w') as handle:
            handle.write('Patient,Sample,R1,R2,status\n')

            for patient in parsed_input_file_dict.keys():

                for sample in parsed_input_file_dict[patient].keys():

                    line_string = patient + ',' + sample + ',' + parsed_input_file_dict[patient][sample]['R1'] + ',' + \
                    parsed_input_file_dict[patient][sample]['R2'] + ',Unprocessed' + '\n'

                    handle.write(line_string)

    else:

        print('Reading existing setup file')

        with open(setup_log, 'r') as handle:
            # setup_log = open(setup_log, 'r')

            completed_sample_run_list = []
            unprocessed_sample_run_list = []

            # for line in setup_log:
            for line in handle:
                line = line.strip('\n')
                line = line.split(',')

                patient = line[0]
                sample_identifier = line[1]
                R1_read = line[2]
                R2_read = line[3]
                status = line[4]

                if status == 'Complete':
                    completed_sample_run_list.append(sample_identifier)
                if status == 'Unprocessed':
                    unprocessed_sample_run_list.append(sample_identifier)

                if patient != 'Patient':

                    if patient not in previous_run_dict.keys():
                        previous_run_dict[patient] = {sample_identifier: {'R1': R1_read, 'R2': R2_read, 'status': status}}

                    elif line[1] not in previous_run_dict[line[0]].keys():
                        previous_run_dict[line[0]][line[1]] = {'R1': line[2], 'R2': line[3], 'status': status}

            logging.info('Unprocessed samples: ' + str(len(unprocessed_sample_run_list)))
            logging.info('Processed samples: ' + str(len(completed_sample_run_list)))
            # setup_log.close()

    return previous_run_dict, setup_log


def update_setuplog(status_to_update, log_file):
    # add path
    print(log_file)
    old_log = log_file.replace('setup_log.csv', 'setup_log_old.csv')
    os.rename(log_file, old_log)

    log_file_new = open(log_file, 'wt')
    log_file_old = open(old_log, 'rt')

    with log_file_old as lineIn:
        with log_file_new as lineOut:
            for line in lineIn:
                lineOut.write(line.replace('Unprocessed', status_to_update))

    os.unlink(old_log)
    log_file_new.close()

    return True


def main(config_file, main_pipeline, haplotype):
    """
    The main function for the pipeline
    :param config_file:
    :param main_pipeline:
    :param haplotype:
    :return:
    """

    print("Parsing the config file:")

    with open(config_file) as json_data_file:
        data = json.load(json_data_file)

    regex_error_rate = data["demiltiplexSettings"]["error_rate"]
    fwd_match_only = data["demiltiplexSettings"]["fwd_only"]
    should_do_blast = False
    if data["demiltiplexSettings"]["do_blast"] == "yes":
        should_do_blast = True

    # ----------------------- Parsing input files -----------------------
    infile_directory = data['input_data']['fastq_dir']
    if infile_directory[-1] != '/':
        infile_directory = infile_directory + '/'

    input_file_dict = process_input_dir(data['input_data']['fastq_dir'])

    prev_run_dict, log_file = check_for_previous_runs(data['input_data']['out_folder'], input_file_dict)

    patient_list = list(input_file_dict.keys())

    prev_run_patient_list = list(prev_run_dict.keys())

    # Catch to see if samples are added and being re-run
    if len(patient_list) != len(prev_run_patient_list) and len(prev_run_patient_list) > 0:
        print('The patient list from the last run and this one do not match.')
        quit()

    # ----------------------- Logging the run parameters -----------------------
    logging.debug(data)
    logging.info(input_file_dict)
    logging.info(data['input_data']['primer_csv'])
    logging.info(data['input_data']['out_folder'])
    logging.info(patient_list)

    # Make sure we are working in the right dir
    out_dir = data['input_data']['out_folder']

    if out_dir[-1] != '/':
        out_dir = out_dir + '/'

    # ----------------------- Logging info -----------------------
    timestr = time.strftime("%Y%m%d-%H%M%S")
    pipeline_logging_name = "Pipeline_{}.log".format(timestr)
    pipeline_logging_file = os.path.join(out_dir, pipeline_logging_name)
    logging.basicConfig(filename=pipeline_logging_file, level=logging.DEBUG)

    '''
    if len(patient_list) > 1:
        logging.warning('The patient list provided is either too long or not a list.')
        quit()
    '''

    # The actual running part.
    test_primer_dict = make_primer_dict(data['input_data']['primer_csv'])
    logging.debug(test_primer_dict)

    test_primer_dict = add_kmer_keys(test_primer_dict)

    logging.info("Primer data parsed. " + str(len(test_primer_dict.keys())) + " primer entries found.")

    # This is where we choose what to do if a run was already done.

    incomplete_run_detected = False
    demultiplex_fastq = False

    if len(prev_run_patient_list) == 0:
        # Then no previous run was detected, and everything must be done.
        create_file_structure = True
        demultiplex_fastq = True

    else:
        # The pipeline has been run, and so the file structure is in place.
        create_file_structure = False

        for a_patient in prev_run_dict.keys():
            for sample in prev_run_dict[a_patient].keys():
                print(prev_run_dict[a_patient][sample])
                if prev_run_dict[a_patient][sample]['status'] == 'Unprocessed':
                    incomplete_run_detected = True

    # Deciding to run demultiplexing again if incomplete run detected.
    #demultiplex_fastq = demultiplex

    if incomplete_run_detected:
        print('An incomplete run has been detected, please remove all files in the output directory and re-run the pipeline.')
        quit()

    run_main_pipe = main_pipeline
    make_haplotypes = haplotype

    print('Create file structure needed: ' + str(create_file_structure))
    print('Incomplete run detected: ' + str(incomplete_run_detected))

    if create_file_structure:
        print("Creating file structure")

        import step_1_create_folders

        # Create folders for each of the gene regions
        for gene_region in test_primer_dict.keys():
            step_1_create_folders.main('./', gene_region, patient_list)

        # Add extra one for reads where no genes matched
        step_1_create_folders.main('./', 'None', patient_list)

    if demultiplex_fastq:
        print("Demultiplexing fastq")

        for a_patient in input_file_dict.keys():

            for a_sample in input_file_dict[a_patient].keys():
                r1_file_path = input_file_dict[a_patient][a_sample]['R1']
                r2_file_path = input_file_dict[a_patient][a_sample]['R2']

                if fwd_match_only == "yes":
                    infast_R1_name = ntpath.basename(r1_file_path)
                    infast_R2_name = ntpath.basename(r2_file_path)

                    split_by_primers_matchpair(r1_file_path, r2_file_path,
                                                test_primer_dict, 'fwd', infast_R1_name, infast_R2_name, out_dir,
                                                a_patient,
                                                regex_error_rate, make_sure=should_do_blast)
                else:
                    infast_name = ntpath.basename(r1_file_path)
                    split_by_primers(r1_file_path, test_primer_dict, 'fwd',
                                    infast_name, out_dir, a_patient, regex_error_rate, make_sure=should_do_blast)

                    infast_name = ntpath.basename(r2_file_path)
                    split_by_primers(r2_file_path, test_primer_dict, 'rev',
                                    infast_name, out_dir, a_patient, regex_error_rate, make_sure=should_do_blast)

        update_complete = update_setuplog('Complete', log_file)
        print('De-multiplex complete: ' + str(update_complete))

    if run_main_pipe:
        print("Running main pipeline")

        import step_2_ngs_processing_pipeline_master_call

        for gene_region, gene_dict in test_primer_dict.items():
            print("Running pipeline for:", gene_region)

            # don't run if the gene_region is None: sequences that couldn't be assigned to a gene region
            if gene_region is not None or gene_region != "None":
                overlap = gene_dict['overlap']
                if overlap == "no":
                    nonoverlap = True
                else:
                    nonoverlap = False

                for a_patient_entry in patient_list:

                    # Adding the required parameters
                    path = out_dir + a_patient_entry + '/' + gene_region
                    sub_region = gene_dict['sub_region']
                    if not sub_region:
                        sub_region = False
                    else:
                        sub_region = sub_region
                    user_ref = False
                    # main(path, name, gene_region, sub_region, fwd_primer, cDNA_primer, nonoverlap, length, run_step,
                    #      run_only, cores)
                    # Calling step 2
                    try:
                        step_2_ngs_processing_pipeline_master_call.main(path,
                                                                        data['pipelineSettings']['out_prefix'],
                                                                        gene_region,
                                                                        sub_region,
                                                                        test_primer_dict[gene_region]['fwd_full'],
                                                                        test_primer_dict[gene_region]['rev_full'],
                                                                        nonoverlap,
                                                                        data['pipelineSettings']['min_read_length'],
                                                                        data['pipelineSettings']['run_step'],
                                                                        False,
                                                                        user_ref,
                                                                        data['pipelineSettings']['cores'])
                    except Exception as e:
                        # Todo is this try except actually necessary
                        print(e)
                        pass

        if make_haplotypes:
            print("Making haplotypes from alignment")

            import step_3_make_haplotpes_from_alignment

            step_3_make_haplotpes_from_alignment.main(
                data['haplotype_settings']['infile'],
                data['haplotype_settings']['field'],
            )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''A de-multiplexer tool Input for primers is csv separated by a comma. 
        Headings include: gene_name,fwd_sequence,fwd_PID_len,rev_sequence,rev_pid,overlapping
        overlapping = True if the the fwd and rev reads overlap, else overlapping = False
        ''', epilog="""Version 0.1""")

    parser.add_argument('-c', '--config_file', type=str, help='Configuration file for the run in JSON format',
                        required=True)
    parser.add_argument('-m', '--no_main_pipeline', default=True, action='store_false',
                        help='Do not run the main pipeline', required=False)
    parser.add_argument('-hap', '--no_haplotype', default=True, action='store_false',
                        help='Do not run the haplotyping pipeline', required=False)
    args = parser.parse_args()

    config_file = args.config_file
    main_pipeline = args.no_main_pipeline
    haplotype = args.no_haplotype

    main(config_file, main_pipeline, haplotype)
