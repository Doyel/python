# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
import string
from pprint import pprint
import cPickle as pickle;


# python GetUniProtIdsforProteins.py <file_loc> <uniprotdict>
# python GetUniProtIdsforProteins.py /ua/ml-group/big-mechanism-project/OverRepresentationScores/top_common_english_words_proteins.txt ./DatumKBProteins_UniProtId.p

def main_body():
    global DatumKBProtein_UniprotID_dict
    parser = argparse.ArgumentParser(prog='GetUniProtIdsforProteins', usage='GetUniProtIdsforProteins.py <file_loc> <uniprotdict>', description='Script to get UniProtIds from common english words that are proteins')
    parser.add_argument('file_loc', help='The location where the list of common eng words selected by me are present')
    parser.add_argument('uniprotdict', help='The location of the DatumKB proteins to UniprotIDs dict')

    args = parser.parse_args()

    if os.path.isfile(args.uniprotdict):
        print "Found - DatumKBProteins_UniProtId.p!"
        DatumKBProtein_UniprotID_dict = pickle.load(open("DatumKBProteins_UniProtId.p", "rb"))
    else:
        print "The dictionary - DatumKBProteins_UniProtId.p was not found in provided directory. Please give the right location!"
        exit()

    engwords_uniprotID_dict = {}
    create_english_dict(args.file_loc, engwords_uniprotID_dict)

    fn = open(args.file_loc, 'w')
    for prot in engwords_uniprotID_dict:
        fn.write(prot + ",  " + engwords_uniprotID_dict[prot] + '\n')
    fn.close()


def create_english_dict(file_loc, engwords_uniprotID_dict):
    if not os.path.isfile(file_loc):
        print "Couldn't locate the file containing a list of common english words that are also protein names - " + file_loc
        exit()

    # Assumption: There are no repetitions of words in text file and all words are in lowercase
    with open(file_loc, 'rt') as f:
        for word in f:
            word = word.strip()
            engwords_uniprotID_dict[word] = DatumKBProtein_UniprotID_dict[word][1]


if __name__ == "__main__":
    main_body()
