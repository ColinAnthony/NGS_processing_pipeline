#!/usr/bin/python
from __future__ import print_function
from __future__ import division
import tempfile
import subprocess
from subprocess import DEVNULL
import argparse
import os
import sys
import collections
import regex
from Bio import SeqIO
#

__author__ = 'Colin Anthony'


def translate_dna(sequence):
    '''
    :param sequence: (str) a DNA sequence string
    :return: (str) a protein string from the forward reading frame 1
    '''
    codontable = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
                  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
                  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
                  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
                  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
                  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
                  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
                  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                  'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
                  'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
                  '---': '-',
                  }
    seq = sequence.upper()
    prot = []
    for n in range(0,len(seq),3):
        if seq[n:n+3] in codontable:
            prot.append(codontable[seq[n:n+3]])
        else:
            prot.append("X")
    return "".join(prot)


def fasta_to_dct_rev(fn):
    '''Converts a fasta file to a dictionary
    :param fn: a fasta file
    :return: a dictionary key = sequence, value = list of sequence names with that sequence
    '''

    dct = collections.defaultdict(list)
    for seq_record in SeqIO.parse(open(fn), "fasta"):
        # TODO replace the list of lists, with a list of strings.
        dct[str(seq_record.seq).replace("-", "").upper()].append([seq_record.description.replace(" ", "_")])
    print("")
    return dct


def fasta_to_dct(fn):
    '''Converts a fasta file to a dictionary
    :param fn: a fasta file
    :return: a dictionary
    '''

    dct = collections.defaultdict(list)
    for seq_record in SeqIO.parse(open(fn), "fasta"):
        dct[seq_record.description.replace(" ", "_")] = str(seq_record.seq).upper()
    return dct


def slice_n_dice(seq, loop_list):
    '''
    :param seq: (str) DNA sequence
    :param loop_list (list) list of V-loops to strip out of sequence
    :return: (str) sliced conserved regions of sequence, sliced v-loop region
    '''
    start_v1 = r'(TTAAA[CT]TG[CT]ACCAAT){s<6}'
    end_v1 = r'(AATTGC[TA]CTTTCAAT[AG][TC]AAC[TC]ACA){s<6}'
    start_v2 = r'(TATA[AG]A[CT][CT]TGAT[AG]TA){s<6}'
    end_v2 = r'(TATA[GT]ATTAATAAATTG[CT]AA[CT]){s<6}'
    start_v3 = r'(AGAAAAAGTATAAGGATA){s<6}' # AAAACAATAATAGTACATCTC #r'AC[AG]GACAATGCTAAAAC[AC]){s<6}'
    end_v3 = r'(GGACCAGGACAAACATTC){s<6}' # GAATCTGTAGAGATTAATTGT #ATAGTACA[CG]CT[TG]AA[TC][AG][AC][AG]TCT){s<6}'
    start_v4 = r'(TGT[AG]GAGGAGAATTTTTCTATTG[CT]){s<6}'
    end_v4 = r'(TGCA[AG]AATAAAACAAATTATAAAC){s<6}'
    start_V5 = r'(GG[AG]CTA[AC]TATTG[AG][CT]ACGTGATGG[AT]GG[AC]){s<6}'
    end_V5 = r'(A[AC][ATC]G[AG]GAAATTC){s<6}'
    loop_regions = {'V1': [start_v1,
                           end_v1],
                    'V2': [start_v2,
                           end_v2],
                    'V3': [start_v3,
                           end_v3],
                    'V4': [start_v4,
                           end_v4],
                    'V5': [start_V5,
                           end_V5],
                    }
    # define loop regions start and end strings
    for region in loop_list:
        if region.upper() not in loop_regions.keys():
            print("Incorrect loop. \nValid options are 'V1, V2, V3, V4 and V5'")
            sys.exit()

    number_loops = len(loop_list)
    last_loop = number_loops - 1

    outseq = []
    outloop = []
    # if only one region, extract regions either side of V-loop
    if number_loops == 1:
        start = loop_regions[loop_list[0]][0]
        end = loop_regions[loop_list[0]][1]
        l_srt = regex.search(start, seq, regex.BESTMATCH)
        # if match found, slice match region out of sequence to prevent mismatches in later regex
        if l_srt is not None:
            c1 = seq[:l_srt.end()]
            seq = seq[l_srt.end():]
        else:
            tag = 'start ' + loop_list[0]
            print(tag, 'Not found')
            return None, None
        l_end = regex.search(end, seq, regex.BESTMATCH)
        # if match found, slice match region out of sequence to prevent mismatches in later regex
        if l_end is not None:
            v_loop = seq[:l_end.start()]
            c2 = seq[l_end.start():]
        else:
            tag = 'end ' + loop_list[0]
            print(tag, 'Not found')
            return None, None
        outseq = [c1, c2]
        outloop.append(v_loop)

    # if more than one loop specified, apply more complex extraction of conserved regions
    elif number_loops > 1:
        for i, lp in enumerate(loop_list):
            # for first loop
            if i == 0:
                start = loop_regions[lp][0]
                end = loop_regions[lp][1]
                l_srt = regex.search(start, seq, regex.BESTMATCH)
                if l_srt is not None:
                    c1 = seq[:l_srt.end()]
                    seq = seq[l_srt.end():]
                else:
                    tag = 'start ' + loop_list[0]
                    print(tag, 'Not found')
                    return None, None
                l_end = regex.search(end, seq, regex.BESTMATCH)
                if l_end is not None:
                    v_loop = seq[:l_end.start()]
                    seq = seq[l_end.start():]
                else:
                    tag = 'end ' + loop_list[0]
                    print(tag, 'Not found')
                    return None, None
                outseq.append(c1)
                outloop.append(v_loop)

            # for last loop
            elif i == last_loop:
                start = loop_regions[lp][0]
                end = loop_regions[lp][1]
                l_srt = regex.search(start, seq, regex.BESTMATCH)
                if l_srt is not None:
                    c3 = seq[:l_srt.end()]
                    seq = seq[l_srt.end():]
                else:
                    tag = 'start ' + loop_list[0]
                    print(tag, 'Not found')
                    return None, None
                l_end = regex.search(end, seq, regex.BESTMATCH)
                if l_end is not None:
                    v_loop = seq[:l_end.start()]
                    cn = seq[l_end.start():]
                else:
                    tag = 'end ' + loop_list[0]
                    print(tag, 'Not found')
                    return None, None
                outseq.append(c3)
                outseq.append(cn)
                outloop.append(v_loop)

            # for all loops in between
            else:
                start = loop_regions[lp][0]
                end = loop_regions[lp][1]
                l_srt = regex.search(start, seq, regex.BESTMATCH)
                if l_srt is not None:
                    c2 = seq[:l_srt.end()]
                    seq = seq[l_srt.end():]
                else:
                    tag = 'start ' + loop_list[0]
                    print(tag, 'Not found')
                    return None, None
                l_end = regex.search(end, seq, regex.BESTMATCH)
                if l_end is not None:
                    v_loop = seq[:l_end.start()]
                    seq = seq[l_end.start():]
                else:
                    tag = 'end ' + loop_list[0]
                    print(tag, 'Not found')
                    return None, None
                outseq.append(c2)
                outloop.append(v_loop)
    return outseq, outloop


def align_cons_dna(tmp_file, loop_s, nam, slave_dict):
    '''
    :param tmp_file: (srt) path to temp folder
    :param loop_s: (list) list of the variable loops to extract
    :param nam: (str) suffix for temp outfile name
    :param slave_dict (dict) dictionary containing the conserved regions f the sequence
    :return: (dict) dictionary of aligned sequences for all conserved regions
    '''
    # write each conserved region to a separate file
    for j in range(len(loop_s) + 1):
        for index, seq_list in slave_dict.items():
            with open(tmp_file + '_' + str(j) + '.fas', 'a') as handle1:
                handle1.write('>' + nam + '_' + str(index) + '\n' + str(seq_list[j]) + '\n')

    aln_seqs = collections.defaultdict(dict)
    print("Aligning the sequences, please wait")
    for j in range(len(loop_s)+1):
        a = tmp_file + '_' + str(j) + '.fas'
        b = a.replace('.fas', '.aln')
        cmd = 'mafft {0} > {1}'.format(a, b)
        subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        aln_seqs[str(j)] = fasta_to_dct(b)
        os.remove(a)
        os.remove(b)
    print("Aligning completed")
    new_aln_seqs = collections.defaultdict(dict)
    s_al = collections.OrderedDict(sorted(aln_seqs.items()))
    for region, dict_al in s_al.items():
        for seq_name, seqs_list in dict_al.items():
            indx = seq_name.split('_')[-1]
            new_aln_seqs[region][indx] = seqs_list
    return new_aln_seqs


def align_cons_prt(tmp_file, loop_s, nam, slave_dict, aligner):
    '''
    :param tmp_file: (srt) path to temp folder
    :param loop_s: (list) list of the variable loops to extract
    :param nam: (str) suffix for temp outfile name
    :param slave_dict (dict) dictionary containing the conserved regions f the sequence
    :return: (dict) dictionary of aligned sequences for all conserved regions
    '''
    for j in range(len(loop_s) + 1):
        with open(tmp_file + '_' + str(j) + '.fas', 'w') as handle1:
            handle1.write('')

    # write each conserved region to a separate file
    for j in range(len(loop_s) + 1):
        for index, seq_list in slave_dict.items():
            with open(tmp_file + '_' + str(j) + '.fas', 'a') as handle1:
                handle1.write('>' + nam + '_' + str(j) + '_' + str(index) + '\n'
                              + str(translate_dna(seq_list[j])) + '\n')

    aln_seqs = collections.defaultdict(dict)
    print("Aligning the sequences, please wait")
    for j in range(len(loop_s)+1):
        a = tmp_file + '_' + str(j) + '.fas'
        b = a.replace('.fas', '.aln')
        if aligner:
            cmd = 'muscle -in {0} -out {1} -seqtype protein'.format(a, b)
        else:
            cmd = 'mafft --amino {0} > {1}'.format(a, b)
        subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        aln_seqs[str(j)] = fasta_to_dct(b)
        os.remove(a)
        os.remove(b)
    print("Aligning completed")

    new_cons_d = collections.defaultdict(dict)
    s_aln_seqs = collections.OrderedDict(sorted(aln_seqs.items()))

    for cons_region, aln_d in s_aln_seqs.items():
        for index, prot in aln_d.items():
            backtrans = ''
            resi_count = 0
            index = index.split('_')[-1]
            cons_region = int(cons_region)
            dnaseq = slave_dict[index][cons_region]
            for resi in prot:
                if resi != '-':
                    backtrans += dnaseq[(resi_count * 3):((resi_count * 3) + 3)]
                    resi_count += 1
                else:
                    backtrans += '---'
            new_cons_d[str(cons_region)][index] = backtrans

    new_cons_d = collections.OrderedDict(sorted(new_cons_d.items()))
    return new_cons_d


def v_loop_handler(loop_dict):
    '''
    :param loop_dict: (dict) dictionary of extracted loop sequences
    :return: (dict) dictionary of extracted loop sequences padded with gaps in the middle
    '''
    length_d = collections.defaultdict(dict)
    max_l = collections.defaultdict(int)
    s = collections.OrderedDict(sorted(loop_dict.items()))
    for ix, loopseql in s.items():
        for cnter, lseq in enumerate(loopseql):
            length_d[ix][cnter] = len(lseq)
            if len(lseq) > max_l[cnter]:
                max_l[cnter] = len(lseq)
    # split loop in half and pad with gaps to longest seq
    new_loop_dict = collections.defaultdict(list)
    for ixd, loopseql in loop_dict.items():
        for c, lseq in enumerate(loopseql):
            l = len(lseq)
            gaps = '-' * (max_l[c] - l)
            left = int((l/3)//2)*3
            new_lseq = lseq[:left] + gaps + lseq[left:]
            new_loop_dict[ixd].append(new_lseq)
    return new_loop_dict


def v_loop_aligner(loop_dict, tmppath, aligner):
    '''
        :param loop_dict: (dict) dictionary of extracted loop sequences
        :param tmppath: (str) path to temp folder including temp outfile prefix
        :return: (dict) dictionary of aligned loop sequences
    '''
    ref_dict = collections.defaultdict(list)
    trans_dict = collections.defaultdict(list)
    num_loops = 0
    tmp_file = os.path.join(tmppath, 'vloop')

    # remove pre-existing files
    for ix, loopseq_l in loop_dict.items():
        for cnter, loops in enumerate(loopseq_l):
            with open(tmp_file + '_' + str(cnter) + '.fas', 'w') as handle1:
                handle1.write('')
    # translate loops from dna to prot
    for ix, loopseq_l in loop_dict.items():
        num_loops = len(loopseq_l)
        for cnter, loops in enumerate(loopseq_l):
            transl = translate_dna(loops)
            trans_dict[ix].append([cnter, transl])
            tmpfas = '>v_loopseq_{0}_{1}\n{2}\n'.format(cnter, ix, transl)
            with open(tmp_file + '_' + str(cnter) + '.fas', 'a') as handle1:
                handle1.write(tmpfas)
                ref_dict[ix].append([loops])
    # align prot sequence for each loop
    aln_loop = collections.defaultdict(dict)
    print("Aligning the sequences, please wait")
    for j in range(num_loops):
        a = tmp_file + '_' + str(j) + '.fas'
        b = a.replace('.fas', '.aln')
        if aligner:
            cmd = 'muscle -in {0} -out {1} -seqtype protein'.format(a, b)
        else:
            cmd = 'mafft --amino {0} > {1} '.format(a, b)
        subprocess.call(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
        aln_loop[str(j)] = fasta_to_dct(b)
        os.remove(a)
        os.remove(b)
    print("Aligning completed")
    # backtranslate protein align to dna
    new_loop_d = collections.defaultdict(list)
    s = collections.OrderedDict(sorted(aln_loop.items()))
    for cnt, dictn in s.items():
        for seq_n, prot in dictn.items():
            idx = seq_n.split('_')[-1]
            cnt = int(cnt)
            dnaseq = loop_dict[idx][cnt]
            backtrans = ''
            resi_count = 0
            for resi in prot:
                if resi == '-':
                    backtrans += '---'
                else:
                    backtrans += dnaseq[(resi_count*3):((resi_count*3) + 3)]
                    resi_count += 1
            new_loop_d[idx].append(backtrans)
    return new_loop_d


def patch_and_repair(aln_seqs, loop_dict):
    '''
    :param aln_seqs: (dict) dictionary of aligned conserverved regions
    :param loop_dict: (dict) dictionary of extracted loop sequences
    :return: (dict) dictionary of lookup indexes and merged conserved and variable regions
    '''
    master_universe = collections.defaultdict(list)
    for cons_region, al_dict in aln_seqs.items():
        for s_idx, sequence in al_dict.items():
            master_universe[s_idx].append([cons_region, sequence])
    for look_up_val, lst in master_universe.items():
        full_seq = []
        item = sorted(lst)
        loop_list = loop_dict[look_up_val]
        for i, ls in enumerate(item):
            full_seq.append(ls[1])
            if i < len(loop_list):
                full_seq.append(loop_list[i])
            else:
                continue
        full_seq = "".join(full_seq)
        master_universe[look_up_val] = full_seq
    return master_universe


def writeout_and_sanity_check(master_universe, finalout, master_ref_dict, l_bad, inf):
    '''
    :param master_universe: (dict) dictionary containing the index and final merged sequence
    :param finalout: (dict) dictionary containing the index and final merged sequence
    :param master_ref_dict: (dict) dictionary containing the index and list of sequence names
    :param l_bad: (int) number of sequences where v-loop not found
    :param inf: (str) path and name of original infule
    :return: (None) writes outfile and prints input/output sequence check
    '''

    final_out_fn = finalout + "_align.fasta"
    with open(final_out_fn, 'w') as handle:
        for inx, seq in master_universe.items():
            names_list = master_ref_dict[inx]
            for n in names_list:
                handle.write('>{0}\n{1}\n'.format(n[0], seq))

    # orig infile
    orig_dct = fasta_to_dct(inf)
    # new output dct name = final_out_fn
    output_dct = fasta_to_dct(final_out_fn)
    # check the sequences are identical to input
    identical_seq = 0
    diffs = 0
    total_inseqs = len(orig_dct)

    # compare the elements of orig_dct and output_dct.
    for seq_id, seq in output_dct.items():
        orig_seq = orig_dct[seq_id].replace("-", "")
        seq = seq.replace("-", "")
        if seq != orig_seq:
            print('>' + str(seq_id) + '_orig')
            print(orig_seq)
            print('>' + str(seq_id) + '_new')
            print(seq)
            print("\n")
            diffs += 1
        else:
            identical_seq += 1

    if identical_seq + l_bad + diffs == total_inseqs:
        sanity = True
    else:
        sanity = False
    if diffs == 0:
        happy = True
    else:
        happy = False
    print('Total Number of actual seqs=                                         {}'.format(total_inseqs))
    print('Number of seqs with loop not found=                                  {}'.format(l_bad))
    print('Number of sequences which are identical to input=                    {}'.format(identical_seq))
    print('Number of sequences differing between infile and outfile=            {}'.format(diffs))
    print('Sanity check (Numbers all add up):                            Sane = {}'.format(sanity))
    print('Happiness check (all out sequences are same as in sequences: Happy = {}'.format(happy))


def main(inf, outp, nam, loop_s, v_loop_align, dna, aligner):
    print(inf)
    # check that correct loop ID has been given
    loop_s = [x.upper() for x in loop_s]
    loop_s = sorted(loop_s)
    print(loop_s)
    # set outfile paths and names
    tmppath = tempfile.gettempdir()
    tmp_file = os.path.join(tmppath, nam)
    not_found = os.path.join(outp, nam + '_not_found')
    finalout = os.path.join(outp, nam)
    # read infile to dictionary, sequence as key, value as list of seq_names
    d = fasta_to_dct_rev(inf)

    # remove old output files from a previous run
    for j, region in enumerate(loop_s):
        with open(tmp_file + '_' + str(j) + '.fas', 'w') as handle1:
            handle1.write('')

    # remove old files from a previous run
    with open(not_found + '.fas', 'w') as handle2:
        handle2.write('')

    # dict to keep track of sequence and name list
    master_ref_dict = collections.OrderedDict()

    # dict of extracted loops and their index
    loop_dict = collections.defaultdict(list)

    # dict of extracted conserved regions and their index
    slave_dict = collections.defaultdict(dict)

    # counter for sequences with loops not found
    l_bad = 0

    for indx, (seq, namelist) in enumerate(d.items()):
        master_ref_dict[str(indx)] = namelist

        # slice conserved and v-loop regions from each sequence
        cons, var = slice_n_dice(seq, loop_s)
        if cons is None:
            # write out sequences with no match to loop region regex
            with open(not_found + '.fas', 'a') as handle2:
                for n in namelist:
                    l_bad += 1
                    handle2.write('>' + str(n[0]) + '_' + str(indx) + '\n' + str(seq) + '\n')
        else:
            # populated a dict with conserved regions and a dict with loop regions
            slave_dict[str(indx)] = cons
            loop_dict[str(indx)] = var

    # if can't find loops in any sequences, exit
    if len(slave_dict) < 2:
        print("loops not found in sequences. Can't align < 2 sequences\n(exiting)")

    # align the different conserved regions
    if dna is False:
        aln_seqs = align_cons_prt(tmp_file, loop_s, nam, slave_dict, aligner)
    else:
        # translate conserved regions to prot and align, then back translate to DNA
        aln_seqs = align_cons_dna(tmp_file, loop_s, nam, slave_dict)

    # handle the v-loops
    if v_loop_align is False:
        loop_dict = v_loop_handler(loop_dict)

    else:
        # align translated loops, then back translate to DNA
        loop_dict = v_loop_aligner(loop_dict, tmppath, aligner)

    # put conserved regions and loops back together
    master_universe = patch_and_repair(aln_seqs, loop_dict)

    # write final output, check that input sequences and names match the aligned output
    writeout_and_sanity_check(master_universe, finalout, master_ref_dict, l_bad, inf)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Aligns sequences with variable loops (HIV env sequences)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', type=str, default=argparse.SUPPRESS,
                        help='The input fasta file with all sequences from all time points', required=True)
    parser.add_argument('-o', '--outpath', type=str, default=argparse.SUPPRESS,
                        help='The path to where the output files will be created', required=True)
    parser.add_argument('-n', '--name', type=str, default=argparse.SUPPRESS,
                        help="The name of your outfile. Don't add a suffix (.fasta)", required=True)
    parser.add_argument('-l', '--loop', type=str, nargs='+', default=argparse.SUPPRESS,
                        help='The V-loops to handle. Valid options are v1, (v2), v3, v4 and v5', required=True)
    parser.add_argument('-v', '--v_loop_align', default=False, action='store_true',
                        help='Align the v-loops (rather than split them in half)', required=False)
    parser.add_argument('-d', '--dna', default=False, action='store_true',
                        help='Use this option if your sequences are not in coding space (reading frame = 1)'
                             '\nConserved regions will be aligned as amino acid translated '
                             'sequences by default.', required=False)
    parser.add_argument('-a', '--aligner', default=False, action='store_true',
                        help='Align with muscle instead of mafft)', required=False)

    args = parser.parse_args()
    infile = args.infile
    outpath = args.outpath
    name = args.name
    loop = args.loop
    v_loop_align = args.v_loop_align
    dna = args.dna
    aligner = args.aligner
    main(infile, outpath, name, loop, v_loop_align, dna, aligner)
