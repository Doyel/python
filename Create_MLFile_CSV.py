# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import collections;
import string;
import json
import cPickle as pickle;
from random import shuffle
from nltk.tokenize import StanfordTokenizer
from nltk.tokenize import sent_tokenize, word_tokenize


# python Create_MLFile_CSV.py <label> <datafiletype> <outputdir>
# python Create_MLFile_CSV.py RNAi Training ./datafiles/RNAi/


def main_body():
    parser = argparse.ArgumentParser(prog='Create_MLFile_CSV', usage='Create_MLFile_CSV.py <label> <datafiletype> <csv_loc>', description='Script to create labelled ML file')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('datafiletype', help='The type of the ML file that needs to be created - training, test or openaccess')
    parser.add_argument('csvLoc', help='The directory in which to store the csv files')
    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    allowed_filetype = ["Training", "Test", "OpenAccess"]
    if args.datafiletype not in allowed_filetype:
        print "The given file type is not allowed! File type can only be any one of \'Training\', \'Test\', \'OpenAccess\'"
        print "Please try again!"
        exit()

    label = args.label
    if args.datafiletype == "Training":
        outfilename = "vw_" + label + "_" + args.datafiletype + "_File_NoOpenAccess.txt"
        csv_file_name = label + "_train.csv"
    elif args.datafiletype == "Test":
        outfilename = "vw_" + label + "_" + args.datafiletype + "_File.txt"
        csv_file_name = label + "_test.csv"
    elif args.datafiletype == "OpenAccess":
        outfilename = "vw_" + label + "_Training_File_OpenAccess.txt"
        csv_file_name = label + "_openaccess.csv"

    csv_full_path = os.path.join(args.csvLoc, csv_file_name)
    if os.path.isfile(csv_full_path):
        print "Found - " + csv_file_name + "! Everything is perfect in the world!"
    else:
        print "Couldn't locate the csv file - " + csv_file_name + "! Please place it in the specified directory."
        exit()

    myrecords = []
    with open(csv_full_path, 'rt') as f1:
        for myline in f1:
            pieces = myline.strip().split("|_|")    # pmid |_| prot |_| class |_| tag |_| passage_text
            vw_label = pieces[2]
            vw_tag = pieces[3]
            lex_feat_str = create_lexical_features(pieces[4].strip())
            rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str + " |ProteinFeature " + pieces[1]
            myrecords.append(rec)

    shuffle(myrecords)
    outfile_full_path = os.path.join(args.csvLoc, outfilename)
    data_fn1 = open(outfile_full_path, 'w')
    for rc in myrecords:
        data_fn1.write(rc + '\n')
    data_fn1.close()


def create_lexical_features(paragraph):
    #tokens = StanfordTokenizer().tokenize(paragraph.strip())
    tokens = word_tokenize(paragraph.strip())
    #tokens = prune_lessfrequent_words(tokens)
    countdict = collections.Counter(tokens)
    cleanse_punct(countdict)
    mystr = ' '.join(['%s:%s' % (k, v) for k, v in countdict.iteritems()])
    return mystr


def cleanse_punct(countdict):
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)
    for word in countdict.keys():
        if (":" in word) or (len(word.strip()) == 0):
            del countdict[word]
            continue
        out = word.translate(remove_punctuation_map)
        if len(out.strip()) == 0:
            del countdict[word]


if __name__ == "__main__":
    main_body()
