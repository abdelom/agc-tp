<<<<<<< HEAD
#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""
import time
import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
from tqdm import tqdm
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 10,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as  monfich:
        seq = ""
        for line in monfich:
            if line.startswith(">"):
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line[:-1]
        yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    generator_seq = read_fasta(amplicon_file, minseqlen)
    dict_seq = {}
    for seq in generator_seq:
        if seq not in dict_seq.keys():
            dict_seq[seq] = 1
        else:
            dict_seq[seq] += 1
    dict_seq = sorted(dict_seq.items(), key=lambda x: x[1], reverse=True)
    for elem in dict_seq:
        if elem[1] >= mincount:
            yield list(elem)


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size]
              for i in range(0, len_seq, chunk_size)
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    generator = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    list_otu = []
    nb_otu = 1
    boule = True
    for sequence1 in generator:
        if boule:
            list_otu.append(list(sequence1))
            boule = False
        else:
            index = 0
            flag = True
            while index < nb_otu :
                align = nw.global_align(sequence1[0], list_otu[index][0], gap_open=-1, gap_extend=-1,\
                matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
                if get_identity(align) > 97.0:
                    flag = False
                    break
                index += 1
            if flag:
                list_otu.append(list(sequence1))
                nb_otu +=1
    return list_otu

def fill(text, width=80):
    """Sp, lit text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as filout:
        for index, otu in enumerate(OTU_list):
            filout.write(f">OTU_{index + 1} occurrence:{otu[1]}\n")
            filout.write(fill(otu[0]) + '\n')

"""
def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    generator_kmer = cut_kmer(sequence, kmer_size)
    for kmer in generator_kmer:
        if kmer not in kmer_dict.keys():
            kmer_dict[kmer] = [id_seq]
        elif id_seq not in kmer_dict[kmer]:
            kmer_dict[kmer].append(id_seq)
    return kmer_dict
"""
def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    l_kmer = get_unique(cut_kmer(sequence, kmer_size))
    for kmer in l_kmer:
        try : 
            kmer_dict[kmer].append(id_seq)
        except: 
            kmer_dict[kmer] = [id_seq]
        #if kmer not in kmer_dict.keys():
         #   kmer_dict[kmer] = [id_seq]
        #else:
         #   kmer_dict[kmer].append(id_seq)
    return kmer_dict

"""
def search_mates(kmer_dict, sequence, kmer_size):
    generator_kmer = common(cut_kmer(sequence, kmer_size), kmer_dict.keys())
    list_tmp = []
    for k_mer in generator_kmer:
        list_tmp = list_tmp + kmer_dict[k_mer]
    best_mates = Counter(list_tmp).most_common(2)
    return [best_mates[0][0], best_mates[1][0]]
"""

def search_mates(kmer_dict, sequence, kmer_size):
    allfound = []
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer not in kmer_dict:
            continue
        allfound += kmer_dict[kmer]
    best_mates = Counter(allfound).most_common(2)
    return [seq_id for seq_id, count in best_mates]


def detect_chimera(perc_identity_matix):
    mean_std = statistics.mean([statistics.stdev(seg) for seg in perc_identity_matix])
    if mean_std < 5.0:
        return False
    flag = 0
    if perc_identity_matix[0][0] > perc_identity_matix[0][1]:
        flag = 1
    for seg in perc_identity_matix:
        if flag == 0 and seg[0] > seg[1]:
            return True
        if flag == 1 and seg[0] < seg[1]:
            return True
    return False

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    generator = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    i = 0
    not_chimeral = []
    k_mer_dict = {}
    counter = 0
    for seq in tqdm(generator, total = len(generator)):
        if i < 2:
            k_mer_dict = get_unique_kmer(k_mer_dict, seq[0], i, kmer_size)
            i += 1
            not_chimeral.append(seq[0])
            yield seq
        else:
            best = search_mates(k_mer_dict, seq[0], kmer_size)
            chunk_list = get_chunks(seq[0], chunk_size)
            chunk_seq_list = [get_chunks(not_chimeral[best[0]], chunk_size)]
            chunk_seq_list += [get_chunks(not_chimeral[best[1]], chunk_size)]
            perc_identity_matrix = [[] for c in range(len(chunk_list))]
            for j in range(len(chunk_seq_list)):
                for l,chunk in enumerate(chunk_list):
                    perc_identity_matrix[l].append(get_identity(
                                nw.global_align(chunk, chunk_seq_list[j][l],
                                    gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                        '../agc')) + "/MATCH")))
            if not detect_chimera(perc_identity_matrix):
                k_mer_dict = get_unique_kmer(k_mer_dict, seq[0], i, kmer_size)
                i += 1
                not_chimeral.append(seq[0])
                yield seq





#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    if args.output_file:
        OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
        write_OTU(OTU_list, args.output_file)

if __name__ == '__main__':
    main()
=======
#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as  monfich:
        seq = ""
        for line in monfich:
            if line.startswith(">"):
                if len(seq) > minseqlen:
                    yield seq
                seq = ""
            else:
                seq += line[:-1]
        yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    generator_seq = read_fasta(amplicon_file, minseqlen)
    dict_seq = {}
    for seq in generator_seq:
        if seq not in dict_seq.keys():
            dict_seq[seq] = 1
        else:
            dict_seq[seq] += 1
    dict_seq = sorted(dict_seq.items(), key=lambda x: x[1], reverse=True)
    for elem in dict_seq:
        if elem[1]> mincount:
            yield list(elem)


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
<<<<<<< HEAD
    generator = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    list_otu = [generator[0]]
    nb_otu = 1
    for sequence1 in generator[1:]:
        index = 0
        flag = True
        print("a")
        while index < nb_otu :
            align = nw.global_align(sequence1[0], list_otu[index][0], gap_open=-1, gap_extend=-1,\
            matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            print(get_identity(align))
            if get_identity(align) > 97.0:
                flag = False
                break
            index += 1
        if flag:
            list_otu.append(list(sequence1))
            nb_otu +=1
=======
    generator = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    list_otu = []
    flag = True
    generator = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    for sequence1 in generator:
        if flag:
            list_otu.append(sequence1)
            flag = False
        else:
            for index, sequence2 in enumerate(list_otu):
                align = nw.global_align(sequence1[0], sequence2[0], gap_open=-1, gap_extend=-1,\
                matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
                if get_identity(align) < 97.0:
                    list_otu.append(list(sequence1))
                elif get_identity(align) > 97.0 and sequence2[1] > sequence1[1]:
                    list_otu[index] = list(sequence1)
                break
>>>>>>> 14b12d7e047ec25e185a1d09e1b3d30a784cdf6a
    return list_otu

def fill(text, width=80):
    """Sp, lit text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as filout:
        for index, otu in enumerate(OTU_list):
<<<<<<< HEAD
            filout.write(f">OTU_{index + 1} occurrence:{otu[1]}\n")
=======
            filout.write(">OTU_{} occurrence:{}\n".format(index + 1, otu[1]))
>>>>>>> 14b12d7e047ec25e185a1d09e1b3d30a784cdf6a
            filout.write(fill(otu[0]) + '\n')


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    generator_kmer = cut_kmer(sequence, kmer_size)
    for kmer in generator_kmer:
        if kmer not in kmer_dict.keys():
            kmer_dict[kmer] = [id_seq]
        elif id_seq not in kmer_dict[kmer]:
            kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    generator_kmer = cut_kmer(sequence, kmer_size)
    list_tmp = []
    for k_mer in generator_kmer:
        if k_mer in kmer_dict.keys():
            list_tmp = list_tmp + kmer_dict[k_mer]
    best_mates = Counter(list_tmp).most_common(2)
    return [best_mates[0][0], best_mates[1][0]]


def detect_chimera(perc_identity_matix):
    mean_std = statistics.mean([statistics.stdev(seg) for seg in perc_identity_matix])
    if mean_std < 5.0:
        return False
    flag = 0
    if perc_identity_matix[0][0] > perc_identity_matix[0][1]:
        flag = 1
    for seg in perc_identity_matix:
        if flag == 0 and seg[0] > seg[1]:
            return True
        if flag == 1 and seg[0] < seg[1]:
            return True
    return False

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    generator = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    i = 0
    not_chimeral = []
    k_mer_dict = {}
    for seq in generator:
        if i < 2:
            k_mer_dict = get_unique_kmer(k_mer_dict, seq[0], i, kmer_size)
            i += 1
            print(i)
            not_chimeral.append(seq[0])
            yield seq
        else:
            best = search_mates(k_mer_dict, seq[0], kmer_size)
            chunk_list = get_chunks(seq[0], chunk_size)
            chunk_seq_list = [get_chunks(not_chimeral[best[0]], chunk_size)]
            chunk_seq_list += [get_chunks(not_chimeral[best[1]], chunk_size)]
            perc_identity_matrix = [[] for c in range(len(chunk_list))]
            for j in range(len(chunk_seq_list)):
                for l,chunk in enumerate(chunk_list):
                    perc_identity_matrix[l].append(get_identity(
                                nw.global_align(chunk, chunk_seq_list[j][l], 
                                    gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                        '../agc')) + "/MATCH")))
            if detect_chimera(perc_identity_matrix):
                k_mer_dict = get_unique_kmer(k_mer_dict, seq[0], i, kmer_size)
                i += 1
                not_chimeral.append(seq[0])
                yield seq





#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    if args.output_file:
        OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
        write_OTU(OTU_list, args.output_file)

if __name__ == '__main__':
    main()
>>>>>>> ea253b91ed4ae8e6c9208cfb266d40096f75d88f
