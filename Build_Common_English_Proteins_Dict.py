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


# python Build_Common_English_Proteins_Dict.py <file_loc>
# python Build_Common_English_Proteins_Dict.py /ua/ml-group/big-mechanism-project/OverRepresentationScores/top_common_english_words_proteins.txt


def main_body():
    global my_eng_dict
    parser = argparse.ArgumentParser(prog='Build_Common_English_Dict', usage='Build_Common_English_Proteins_Dict.py <file_loc>', description='Script to build a dict for common english words that are also proteins')
    parser.add_argument('file_loc', help='The location where the list of common eng words selected by me are present')

    args = parser.parse_args()

    my_eng_dict = {}
    create_english_dict(args.file_loc)
    pickle.dump(my_eng_dict, open("common_english_words_proteins_dict.p", "wb"))


def create_english_dict(path_engdict):
    global my_eng_dict
    if not os.path.isfile(path_engdict):
        print "Couldn't locate the file containing a list of common english words that are also protein names - " + path_engdict
        exit()

    # Assumption: There are no repetitions of words in text file and all words are in lowercase
    with open(path_engdict, 'rt') as f:
        for word in f:
            word = word.strip()
            start_char = word[0]
            if start_char not in my_eng_dict:
                my_eng_dict[start_char] = {}
            if word in my_eng_dict[start_char]:
                print "Duplicate words present in english words text file"
                exit()
            if not word.islower():
                print "English word in text file is not in lowercase"
                exit()
            my_eng_dict[start_char][word] = 1


if __name__ == "__main__":
    main_body()