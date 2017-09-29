# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
from nltk import ngrams, corpus
import string
import json
from pprint import pprint
#import enchant


# python Measure_ProteinCounts.py <class> <psg_dict_loc>
# python Measure_ProteinCounts.py Neg ./2of4brif.txt ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0


def main_body():
    global my_eng_dict
    parser = argparse.ArgumentParser(prog='Measure_ProteinCounts', usage='Measure_ProteinCounts.py <class> <path_engdict> <psg_dict_loc>', description='Script to count protein mentions')
    parser.add_argument('myclass', help='The proteins associated with the provided class are counted')
    parser.add_argument('path_engdict', help='Location of the english dictionary which is a list of english words')
    parser.add_argument('psg_dict_loc', help='The location where the posprot or negprot dict files are present')

    args = parser.parse_args()
    #eng_dict = enchant.Dict("en_US")

    my_eng_dict = {}
    create_english_dict(args.path_engdict)

    allowed_classes = ['Pos', 'Neg']
    if args.myclass not in allowed_classes:
        print "The provided class is not allowed: ", args.myclass, ". Class may only be either \'Pos\' or \'Neg\'"
        print "Please try again!"
        exit()

    protein_counts = {}
    for passagedictfile in ["openaccess_passage_dict.json", "train_passage_dict.json"]:

        fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
        if os.path.isfile(fullpath_passagedictfile):
            file_handle = open(fullpath_passagedictfile, "rb")
            passage_dict = json.load(file_handle)
            file_handle.close()
        else:
            print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
            exit()

        for pmid in passage_dict:
            for label in passage_dict[pmid]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                for prot in passage_dict[pmid][label][args.myclass]:
                    if prot not in protein_counts:
                        protein_counts[prot] = 1
                    else:
                        protein_counts[prot] += 1

    prot_counts_sorted = sorted(protein_counts, key=lambda k: protein_counts[k], reverse=True)
    for k in prot_counts_sorted:
        w_list = k.split()
        if len(w_list) > 1:
            all_eng = True
            for w in w_list:
                if not w.isdigit() and not is_englishword(w):   # Assumption: If a word is not a number AND not a common english word, then that word is part of a protein
                    all_eng = False; break
            if all_eng:
                print k, "\t:-:\t", protein_counts[k]
        else:
            if k.isdigit() or is_englishword(k):
                print k, "\t:-:\t", protein_counts[k]


def is_englishword(word):
    start_char = word[0]
    if start_char not in my_eng_dict:
        return False
    if word in my_eng_dict[start_char]:
        return True
    else:
        return False


def create_english_dict(path_engdict):
    global my_eng_dict
    alphabets = "abcdefghijklmnopqrstuvwxyz"
    if not os.path.isfile(path_engdict):
        print "Couldn't locate the english dictionary text file - " + path_engdict
        exit()

    for alphabet in list(alphabets):
        my_eng_dict[alphabet] = {}
        my_eng_dict[alphabet][alphabet] = 1

    # Assumption: There are no repetitions of words in text file and all words are in lowercase
    with open(path_engdict, 'rt') as f:
        for word in f:
            word = word.strip()
            start_char = word[0]
            if word in my_eng_dict[start_char]:
                print "Duplicate words present in english words text file"
                exit()
            if not word.islower():
                print "English word in text file is not in lowercase"
                exit()
            my_eng_dict[start_char][word] = 1


if __name__ == "__main__":
    main_body()