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


# python Create_MLFile_OpenAccess.py <label> <datafiletype>
# python Create_MLFile_OpenAccess.py RNAi OpenAccess


def main_body():
    parser = argparse.ArgumentParser(prog='Create_MLFile_OpenAccess.py', usage='Create_MLFile_OpenAccess.py <label> <datafiletype>', description='Script to create labelled ML file')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('datafiletype', help='The type of the ML file that needs to be created - training, test or openaccess')    
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
        passagedictfile = "train_passage_dict.json"
    elif args.datafiletype == "Test":
        outfilename = "vw_" + label + "_" + args.datafiletype + "_File.txt"
        passagedictfile = "test_passage_dict.json"
    elif args.datafiletype == "OpenAccess":
        outfilename = "vw_" + label + "_Training_File_OpenAccess.txt"
        passagedictfile = "openaccess_passage_dict.json"
        
    if os.path.isfile(passagedictfile):
        print "Found - " + passagedictfile + "! Everything is perfect in the world!"
        file_handle = open(passagedictfile, "rb")
        passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in the current directory."
        exit()

    #print sorted(passage_dict.keys()); exit()    
    
    myrecords = []
    for pmid in sorted(passage_dict.keys()):
        print "PMID: ", pmid
        for myclass in passage_dict[pmid][label]:               # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            for prot in passage_dict[pmid][label][myclass]:     # passage_dict[pmid][label][myclass] is a dict
                for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:     # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                    vw_label = contg_psg_dict["class"]                              # contg_psg_dict represents one training instance
                    vw_tag = contg_psg_dict["tag"]
                    lex_feat_str = create_lexical_features(contg_psg_dict["textOfInterest"])
                    rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str
                    myrecords.append(rec)

    shuffle(myrecords)
    data_fn1 = open(outfilename, 'w')
    for rc in myrecords:
        data_fn1.write(rc + '\n')
    data_fn1.close()


def create_lexical_features(paragraph):
    tokens = StanfordTokenizer().tokenize(paragraph.strip())
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
