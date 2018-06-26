#!/usr/bin/python3
from __future__ import print_function
from __future__ import division
import csv
import fnmatch
from subprocess import call
from subprocess import Popen
from subprocess import PIPE
import argparse
import os
import ntpath
import json


__author__ = 'Jon Ambler'


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


def make_seq_wild(inputSequence):
    """
    Replaces non-standard nucleotides with a '?' character for use in the regex function.
    :param inputSequence: A nucleotide sequence as a string.
    :return: The nucleotide sequence with replacements as a string.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    wild_seq = ''

    for nuc in inputSequence:
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
    Each substring is represented only once in the list of keys.
    :param list_of_strings:
    :param min_len:
    :param max_len:
    :param max_kmers:
    :return:
    """
    string_dict = {}
    unique_dict = {}
    for a_string in list_of_strings:
        sub_string_list = break_into_all_kmers(a_string)
        string_dict[a_string] = sub_string_list
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
                # print len(unique_dict[a_string])

                count += 1
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

    return res_dict


def export_line_to_fastq(fileTag, headerLine, seqLine, plusLine, scoreLine, infast_name, out_dir, patient_list):
    """
    For the export of a fastq sequence, line by line.
    :param fileTag:
    :param headerLine:
    :param seqLine:
    :param plusLine:
    :param scoreLine:
    :return:
    """
    out_file_name = infast_name.replace('multiplex', fileTag)
    out_file_path = out_dir + patient_list[0] + '/' + fileTag + '/0new_data/'

    # file_obj = open(out_file_path + out_file_name, "a+")
    # file_obj.write(headerLine)
    # file_obj.write(seqLine)
    # file_obj.write(plusLine)
    # file_obj.write(scoreLine)
    # file_obj.close()
    with open(out_file_path + out_file_name, "a+") as handle:
        handle.write(headerLine)
        handle.write(seqLine)
        handle.write(plusLine)
        handle.write(scoreLine)

    return True


def wildcard_seq_match(primer, sequence):
    """
    Check if the primer sequence is found in the sequence, allowing for wildcards.
    :param primer:
    :param sequence:
    :return:
    """

    # iupac_list = ['R', 'Y', 'S', 'W', 'K']
    # for char in iupac_list:
    #     primer = primer.replace(char, '?')

    # Make sure variable is empty
    lst = [sequence]
    # filtered = []

    filtered = fnmatch.filter(lst, primer)
    if len(filtered) == 1:
        return True
    else:
        return False


def add_kmer_keys(primerDict):
    """
    Add the primer keys to an existing primer dictionary.
    :param primerDict:
    :return:
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


def split_by_primers(fastq_file, primer_dict, orientation, infast_name, out_dir, patient_list, make_sure=True):
    """
    Take an input fastq file and split it into individual fastq files, with the split based on the presence of
    a primer sequence specified in a dictionary.
    :param fastq_file:
    :param primer_dict:
    :param orientation: Whether these are forward (R1) or reverse (R2) reads. Options are 'fwd' or 'rev'.
    :param make_sure: When set to True, a blast search is included in the process.
    :return:
    """
    line_number_quality = 1
    line_number_plus = 2
    line_number_sequence = 3
    line_number_header = 4

    if make_sure:
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

        # Get the shortest primer length
        shortest_primer_length = 1000000
        for aPrimer in list(primer_dict.keys()):
            if len(primer_dict[aPrimer][orientation]) < shortest_primer_length:
                shortest_primer_length = len(primer_dict[aPrimer][orientation])

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
                elif wildcard_seq_match(primer_dict[geneRegion][orientation + '_wild'], seq_primer_region):
                    detected_primer = geneRegion
                    primer_dict[geneRegion]['regex_matches_found'] += 1

                # Level 3: If no wildmatch, look for kmer matches
                else:
                    for kmer in primer_dict[geneRegion][orientation + '_keys']:
                        if kmer in seq_primer_region:
                            detected_primer = geneRegion
                            primer_dict[geneRegion]['kmer_matches_found'] += 1

            # Finally, if there is still no match, use the blastDB
            if detected_primer == 'None':

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


def main(config_file, output_dir, demultiplex, main_pipeline, haplotype):
    print("Parsing the config file:")

    with open(config_file) as json_data_file:
        data = json.load(json_data_file)
    print(data)

    print(data['input_data']['fwd_fastq_file'])
    print(data['input_data']['rev_fastq_file'])
    print(data['input_data']['primer_csv'])
    print(data['input_data']['out_folder'])
    print(data['input_data']['patient_list'])

    # Make sure we are working in the right dir
    out_dir = data['input_data']['out_folder']
    if out_dir[-1] != '/':
        out_dir = out_dir + '/'

    patient_list = data['input_data']['patient_list']
    if len(patient_list) > 1:
        print("The patient list provided is either too long or not a list.")
        quit()

    # The actual running part.
    test_primer_dict = make_primer_dict(data['input_data']['primer_csv'])
    print(test_primer_dict)

    test_primer_dict = add_kmer_keys(test_primer_dict)

    print("Primer data parsed. " + str(len(test_primer_dict.keys())) + " primer entries found.")

    os.chdir(output_dir)

    create_file_structure = True
    demultiplex_fastq = demultiplex
    run_main_pipe = main_pipeline
    make_haplotypes = haplotype

    if create_file_structure:
        print("Creating file structure")

        import step_1_create_folders

        for gene_region in test_primer_dict.keys():
            step_1_create_folders.main('./', gene_region, patient_list)

        # Add extra one for reads where no genes matched
        step_1_create_folders.main('./', 'None', patient_list)

    if demultiplex_fastq:
        print("Demultiplexing fastq")
        infast_name = ntpath.basename(data['input_data']['fwd_fastq_file'])
        split_by_primers(data['input_data']['fwd_fastq_file'], test_primer_dict, 'fwd',
                         infast_name, output_dir, patient_list)

        infast_name = ntpath.basename(data['input_data']['rev_fastq_file'])
        split_by_primers(data['input_data']['rev_fastq_file'], test_primer_dict, 'rev',
                         infast_name, output_dir, patient_list)

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

                # Adding the required parameters
                path = output_dir + patient_list[0] + '/' + gene_region
                sub_region = gene_dict['sub_region']
                if not sub_region:
                    sub_region = False
                else:
                    sub_region = sub_region
                user_ref = False

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
                                                                    )
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
                data['haplotype_settings']['script_folder']
            )

    # Paths for use in testing
    '''
    primer_file = '/Volumes/External/CIDRI_Data/HIV_group_Carolyn/Multiplex_regions_primers.csv'
    fastq_file = '/Volumes/External/CIDRI_Data/HIV_group_Carolyn/test_dataset/CAP206-2000-008wpi-multiplex_S17_L001_R1_001.fastq'
    out_dir = '/Volumes/External/CIDRI_Data/HIV_group_Carolyn/test_results/'	
    # For testing
    list_of_seqs = ['GGCTGTGGTATATAAAAATATTCATMATGA', 'GCCATAAGAAAAGCCATATTAGGAC', 'ATAAGACAGGGCTTTGAAGCAGC',
                    'AGCAGAGAGCTTCAGGTTCG', 'ACACATAGCTTTAATTGTRGAGGAGAATTT']

    '''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''A de-multiplexer tool Input for primers is csv separated by a comma. 
        Headings include: gene_name,fwd_sequence,fwd_PID_len,rev_sequence,rev_pid,overlapping
        overlapping = True if the the fwd and rev reads overlap, else overlapping = False
        ''',  epilog="""Version 0.1""")

    parser.add_argument('-c', '--config_file', type=str, help='Configuration file for the run in JSON format')
    parser.add_argument('-o', '--output_dir', type=str, help='Location to write the output fastq files')
    parser.add_argument('-d', '--demultiplex', default=True, action='store_false',
                        help='Do not run the demultiplexing step')
    parser.add_argument('-m', '--main_pipeline', default=True, action='store_false',
                        help='Do not run the main pipeline')
    parser.add_argument('-hap', '--haplotype', default=True, action='store_false',
                        help='Do not run the haplotyping pipeline')
    args = parser.parse_args()

    config_file = args.config_file
    output_dir = args.output_dir
    main_pipeline = args.main_pipeline
    haplotype = args.haplotype
    demultiplex = args.demultiplex

    main(config_file, output_dir, demultiplex, main_pipeline, haplotype)
