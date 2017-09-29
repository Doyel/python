# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
from collections import Counter;
import string, ntpath
import json
import cPickle as pickle;
from random import shuffle
from nltk.tokenize import StanfordTokenizer
from nltk.tokenize import sent_tokenize, word_tokenize


# python Create_MLFile.py <label> <datafiletype> <outputdir> <path_passage_dict>                       # [--feedback <Path to feedback passage dict>]    Optional
# python Create_MLFile.py RNAi OpenAccess ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/RNAi/ ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/openaccess_passage_dict.json
# python Create_MLFile.py RNAi Training ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/RNAi/ ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/train_passage_dict.json
# python Create_MLFile.py RNAi Test ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/RNAi/ ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/test_passage_dict.json
# python Create_MLFile.py RNAi FeedBack ./user_feedback_JSON/Apr4_2017/datafiles/RNAi/ ./user_feedback_JSON/Apr4_2017/openaccess_passage_dict_WithFeedback.json


def main_body():
    #global allowed_words
    #parser = argparse.ArgumentParser(prog='Create_MLFile.py', usage='Create_MLFile.py <label> <datafiletype> <outputdir> [--feedback <Path to feedback passage dict>]', description='Script to create labelled ML file')

    parser = argparse.ArgumentParser(prog='Create_MLFile_LexicalOnly.py',
                                     usage='Create_MLFile.py <label> <datafiletype> <outputdir> <passagedictPath>',
                                     description='Script to create labelled ML file')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('datafiletype', help='The type of the ML file that needs to be created - training, test or openaccess')
    parser.add_argument('outputdir', help='The directory in which to store the ML data files')
    parser.add_argument('passagedictPath', help='The directory in which the passage dict JSON file is located')
    #parser.add_argument('-f', '--feedback', help='Path to the feedback passage dict json file')  # Optional arg

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    allowed_filetype = ["Training", "Test", "OpenAccess", "FeedBack"]
    if args.datafiletype not in allowed_filetype:
        print "The given file type is not allowed! File type can only be any one of \'Training\', \'Test\', \'OpenAccess\'"
        print "Please try again!"
        exit()

    #allowed_words = pickle.load(open("allowed_words_dict_3.p", "rb"))
    #allowed_words = pickle.load(open("allowed_words_dict_1.p", "rb"))

    label = args.label
    passagedictfile = args.passagedictPath
    fname = ntpath.basename(args.passagedictPath).strip()

    if args.datafiletype == "Training":
        outfilename = "vw_" + label + "_" + args.datafiletype + "_File_NoOpenAccess.txt"
        exp_fname = "train_passage_dict.json"
    elif args.datafiletype == "Test":
        outfilename = "vw_" + label + "_" + args.datafiletype + "_File.txt"
        exp_fname = "test_passage_dict.json"
    elif args.datafiletype == "OpenAccess":
        outfilename = "vw_" + label + "_Training_File_OpenAccess.txt"
        exp_fname = "openaccess_passage_dict.json"
    elif args.datafiletype == "FeedBack":
        outfilename = "vw_" + label + "_Training_File_OpenAccess_WithFeedback.txt"
        exp_fname = "openaccess_passage_dict_WithFeedback.json"

    if fname != exp_fname:
        print "Incorrect passage dict file was provided! Expected file is: ", exp_fname
        exit()

    print "Datafile to be created: ", outfilename

    if os.path.isfile(passagedictfile):
        print "Found - " + fname + "! Everything is perfect in the world!"
        sys.stdout.flush()
        file_handle = open(passagedictfile, "rb")
        passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + fname + " in the given directory! Please give the correct path."
        exit()

    myrecords = []; num_rec_skipped = 0
    for pmid in passage_dict:
        print "PMID: ", pmid
        sys.stdout.flush()
        for myclass in passage_dict[pmid][label]:               # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            for prot in passage_dict[pmid][label][myclass]:     # passage_dict[pmid][label][myclass] is a dict
                for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:     # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                    vw_label = contg_psg_dict["class"]                              # contg_psg_dict represents one training instance
                    vw_tag = contg_psg_dict["tag"]
                    lex_feat_str = create_lexical_features(contg_psg_dict["textOfInterest"])
                    if len(lex_feat_str) == 0:
                        num_rec_skipped += 1
                        continue
                    if "weight" in contg_psg_dict:
                        rec = vw_label + " " + contg_psg_dict["weight"] + " " + vw_tag + "|LexicalFeatures " + lex_feat_str  # + " |ProteinFeature " + prot
                    else:
                        rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str #+ " |ProteinFeature " + prot
                    myrecords.append(rec)

    print "No. of records that were skipped due to frequency pruning: ", num_rec_skipped
    sys.stdout.flush()
    shuffle(myrecords)
    outfile_full_path = os.path.join(args.outputdir, outfilename)
    silentremove(outfile_full_path)
    data_fn1 = open(outfile_full_path, 'w')
    for rc in myrecords:
        data_fn1.write(rc + '\n')
    data_fn1.close()


def create_lexical_features(sentList):
    #tokens = StanfordTokenizer().tokenize(paragraph.strip())
    tokens = []
    for sent in sentList:
        tokens.extend(word_tokenize(sent.strip()))
    countdict = Counter(tokens)
    cleanse_punct(countdict)
    #keep_allowed_words_count(countdict)
    if not countdict:   # If all tokens are deleted due to frequency pruning
        return ""
    mystr = ' '.join(['%s:%s' % (k, v) for k, v in countdict.iteritems()])
    return mystr


def cleanse_punct(countdict):
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)
    for word in countdict.keys():
        if (":" in word) or ("|" in word) or (len(word.strip()) == 0):
            del countdict[word]
            continue
        out = word.translate(remove_punctuation_map)
        if len(out.strip()) == 0:
            del countdict[word]
            continue
        #if len(word) <= 3 and containsPunct(word):
            #print word
            #del countdict[word]     # To prune out tokens like n't, 's, 'll


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured


def keep_allowed_words_alphabet(countdict):
    for word in countdict.keys():
        start_char = word[0]
        try:
            if word not in allowed_words[start_char]:   # allowed_words[start_char] is a list
                del countdict[word]
        except KeyError:
            del countdict[word]


def keep_allowed_words_count(countdict):
    for word in countdict.keys():
        if word not in allowed_words:   # allowed_words is a dict
            del countdict[word]


def containsPunct(word):
    for c in set(string.punctuation):
        if c in word: return 1
    return 0


if __name__ == "__main__":
    main_body()
