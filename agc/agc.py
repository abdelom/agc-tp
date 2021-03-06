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
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=10,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    """
    Keyword arguments:
    amplicon file, str: entry file's name .fasta.gz
    minseqlen, int: minimum lengh neccessary to conserve
    a sequence
    return value:
    sequence generator
    """
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
    """
    Keyword arguments:
    amplicon file, str: entry file's name .fasta.gz
    minseqlen, int: minimum lengh neccessary to conserve
    a sequence
    mincount, int: minimum count neccessary to conserve
    a sequence

    return value:
    sequence generator
    """
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
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size]
            for i in range(0, len_seq, chunk_size)
            if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers
    Keyword arguments:
    sequence, str: sequence to cur into kmer
    kmer_size, int: kmer length

    return value:
    kmer generator
    """
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]


def get_identity(alignment_list):
    """Prend en une liste de s??quences align??es au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux.
    Keyword arguments:
    alignement_list: liste de deux str
    un alignment
    return value:
    float, identity percentage between the two sequence
    """
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def fill(text, width=80):
    """Sp, lit text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    """sauvegarde des otus dans un fichier au format fasta
    Keyword arguments:
    otu_LIST, list: liste d'otu
    output_file, str: nom du fichier de sauvegarde
    """
    with open(output_file, "w") as filout:
        for index, otu in enumerate(OTU_list):
            filout.write(f">OTU_{index + 1} occurrence:{otu[1]}\n")
            filout.write(fill(otu[0]) + '\n')


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
    arguments:
    kmer_dict, dictionnaire: cl??: kmer, valeur: liste de s??quence contenat le kmer
    sequence, str: s??quence
    id_seq: id de la s??quence
    kmer_size, int: taille de kmer
    return value:
    kmer_dict, dictionnaire: cl??: kmer, valeur: liste de s??quence contenat le kmer
    """
    l_kmer = get_unique(cut_kmer(sequence, kmer_size))
    for kmer in l_kmer:
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """
    arguments:
    kmer_dict, dictionnaire: cl??: kmer, valeur: liste de s??quence contenat le kmer
    sequence, str: chunk
    id_seq: id de la s??quence
    kmer_size, int: taille de kmer
    return value:
    list, de  taille 2 avec l'id des deux s??quence partageant le plus de kmer avec
    le chunk
    """
    liste = []
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer not in kmer_dict:
            continue
        liste += kmer_dict[kmer]
    best = Counter(liste).most_common(2)
    return [best[0][0], best[1][0]]


def detect_chimera(perc_identity_matix):
    """
    d??termine si une s??quence est chm??rique ou non
    argument:
    perc_identity_matix, matrice: contient le taux d'itentit?? entre
    les chunk et les s??quences parentes
    return value:
    Bool??en, false, si la s??quence n'est pas chim??rique, true sinon
    """
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
    """retire les s??quences chim??rique
    Keyword arguments:
    amplicon file, str: entry file's name .fasta.gz
    minseqlen, int: minimum lengh neccessary to conserve
    a sequence
    mincount, int: minimum count neccessary to conserve
    a sequence
    chunk_size, int, la taille des chunk
    kmer_size, int: la taille des kmer
    return value:
    generateur de s??quences
    """
    generator = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    i = 0
    not_chimeral = []
    k_mer_dict = {}
    for seq in generator:
        if i < 2:
            k_mer_dict = get_unique_kmer(k_mer_dict, seq[0], i, kmer_size)
            i += 1
            not_chimeral.append(seq[0])
            yield seq
        else:
            best = []
            chunk_list = get_chunks(seq[0], chunk_size)
            for chunk in chunk_list:
                best += search_mates(k_mer_dict, chunk, kmer_size)
            best = Counter(best).most_common(2)
            chunk_seq_list = [get_chunks(not_chimeral[best[0][0]], chunk_size)]
            chunk_seq_list += [get_chunks(not_chimeral[best[1][0]], chunk_size)]
            perc_identity_matrix = [[] for c in range(len(chunk_list))]
            for j in range(len(chunk_seq_list)):
                for l, chunk in enumerate(chunk_list):
                    perc_identity_matrix[l].append(get_identity(
                        nw.global_align(chunk, chunk_seq_list[j][l],
                            gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),
                                '../agc')) + "/MATCH")))
            if not detect_chimera(perc_identity_matrix):
                k_mer_dict = get_unique_kmer(k_mer_dict, seq[0], i, kmer_size)
                i += 1
                not_chimeral.append(seq[0])
                yield seq


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """effectue les diff??rentes ??tapes de filtrage.
    une fois les chim??res enlev??s, la fonction r??alise une approche gloutone pour constitu??
    la lisste d'otu
    Keyword arguments:
    amplicon file, str: entry file's name .fasta.gz
    minseqlen, int: minimum lengh neccessary to conserve
    a sequence
    mincount, int: minimum count neccessary to conserve
    a sequence
    return value:
    list, la liste d'otu
    """
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
            while index < nb_otu:
                align = nw.global_align(sequence1[0], list_otu[index][0], gap_open=-1,\
                gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))
                if get_identity(align) > 97.0:
                    flag = False
                    break
                index += 1
            if flag:
                list_otu.append(list(sequence1))
                nb_otu += 1
    return list_otu


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
        start = time.time()
        OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen,\
        args.mincount, args.chunk_size, args.kmer_size)
        write_OTU(OTU_list, args.output_file)
        print(time.time() - start)

if __name__ == '__main__':
    main()
