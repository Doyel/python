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
from random import shuffle
import cPickle as pickle
from nltk.tokenize import sent_tokenize, word_tokenize
from Utilities import read_config_file


# python Create_MLFile.py <label> <datafiletype> <outputdir> <path_passage_dict> <selectedFeatJSON> <allowed_words>                      # [--feedback <Path to feedback passage dict>]    Optional
# python Create_MLFile.py RNAi OpenAccess ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/RNAi/ ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/openaccess_passage_dict.json Selected_DepParseFeatures.json allowed_words_dict_20.p
# python Create_MLFile.py RNAi Training ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/RNAi/ ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/train_passage_dict.json Selected_DepParseFeatures.json allowed_words_dict_20.p
# python Create_MLFile.py RNAi Test ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/RNAi/ ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/test_passage_dict.json Selected_DepParseFeatures.json allowed_words_dict_20.p
# python Create_MLFile.py RNAi FeedBack ./user_feedback_JSON/Apr4_2017/datafiles/RNAi/ ./user_feedback_JSON/Apr4_2017/openaccess_passage_dict_WithFeedback.json Selected_DepParseFeatures.json allowed_words_dict_20.p


def main_body():
    #parser = argparse.ArgumentParser(prog='Create_MLFile.py', usage='Create_MLFile.py <label> <datafiletype> <outputdir> [--feedback <Path to feedback passage dict>]', description='Script to create labelled ML file')

    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Create_MLFile.py',
                                     usage='Create_MLFile.py <label> <datafiletype> <outputdir> <passagedictPath> <selectedFeatJSON> <allowed_words>',
                                     description='Script to create labelled ML file')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('datafiletype', help='The type of the ML file that needs to be created - training, test or openaccess')
    parser.add_argument('outputdir', help='The directory in which to store the ML data files')
    parser.add_argument('passagedictPath', help='The directory in which the passage dict JSON file is located')
    parser.add_argument('selectedFeatJSON', help='The JSON file that contains the selected dep parse features')
    parser.add_argument('allowedWords', help='Pickle file that lists all words allowed to be lexical features')
    #parser.add_argument('-f', '--feedback', help='Path to the feedback passage dict json file')  # Optional arg

    args = parser.parse_args()
    pickle_init(args.allowedWords)

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

    if os.path.isfile(os.path.join(json_location, args.selectedFeatJSON)):
        print "Found - " + args.selectedFeatJSON + "! The world's order continues to be maintained!"
        file_handle = open(os.path.join(json_location, args.selectedFeatJSON), "rb")
        selecteddict = json.load(file_handle)
        file_handle.close()
    else:
        print "The dependency parse file containing the selected features - " + args.selectedFeatJSON + " was not found in the JSON directory. Exiting!"
        exit()

    myrecords = []; num_rec_skipped = 0
    for pmid in passage_dict:
        print "PMID: ", pmid
        sys.stdout.flush()
        for myclass in passage_dict[pmid][label]:               # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            for prot in passage_dict[pmid][label][myclass]:     # passage_dict[pmid][label][myclass] is a dict
                matched_prots = passage_dict[pmid][label][myclass][prot]["MatchedProtein"]
                for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:     # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                    sent_ids = getSentenceIds(contg_psg_dict, fname)
                    vw_label = contg_psg_dict["class"]                              # contg_psg_dict represents one training instance
                    vw_tag = contg_psg_dict["tag"]
                    lex_feat_str = create_lexical_features(contg_psg_dict["textOfInterest"])
                    depparse_feat_str = create_depParse_featuresStr(pmid, sent_ids, matched_prots, selecteddict[label])
                    if len(lex_feat_str) == 0:  # If there are no lexical features (i.e. no words) in a passage, then there could not be any viable depparse features extracted from that psg as well
                        num_rec_skipped += 1
                        continue
                    if "weight" in contg_psg_dict:
                        if len(depparse_feat_str) > 0:
                            rec = vw_label + " " + contg_psg_dict["weight"] + " " + vw_tag + "|LexicalFeatures " + lex_feat_str + " |DepParseFeatures " + depparse_feat_str
                        else:
                            rec = vw_label + " " + contg_psg_dict["weight"] + " " + vw_tag + "|LexicalFeatures " + lex_feat_str
                    else:
                        if len(depparse_feat_str) > 0:
                            rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str + " |DepParseFeatures " + depparse_feat_str
                        else:
                            rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str
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


def create_depParse_featuresStr(pmid, sent_ids, matched_prots, mydict):
    feat_str = []
    if pmid in mydict:
        for sentid in sent_ids: # Collect all high frequency dep parse features that may have been extracted from the
            sentid = str(sentid)    # sentences constituting the passage representing the current train/test instance
            if sentid in mydict[pmid]:  # Also we need to collect only those dep parse features in the sentences that
                for mtchprot in matched_prots:  # are sourced out of the proteins corresponding to the current protein
                    if mtchprot in mydict[pmid][sentid]:    # 'prot' in passage_dict
                        feat_str.extend(mydict[pmid][sentid][mtchprot])
            #else:
            #    print "The current sentence has no high frequency dep parse features extracted out of it!", sentid
    #else:
    #    print "None of the sentences in the article have any high frequency dep parse features extracted out of it", pmid
    feat_str = list(set(feat_str))
    if len(feat_str) > 0:
        return " ".join(feat_str)
    else:
        return ""


def create_lexical_features(sentList):
    #tokens = StanfordTokenizer().tokenize(paragraph.strip())
    tokens = []
    for sent in sentList:
        tokens.extend(word_tokenize(sent.strip()))
    countdict = Counter(tokens)
    cleanse_punct(countdict)
    keep_allowed_words_count(countdict)
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


def getSentenceIds(train_instance_dict, fname):
    sent_ids = []
    if fname.startswith("openaccess"):
        for para_dict in train_instance_dict["details"]:    # train_instance_dict["details"] is a list of paragraphs
            for sent_dict in para_dict["sentenceDetails"]:  # para_dict["sentenceDetails"] is a list of positive sentences within a paragraph
                sent_ids.append(sent_dict["absoluteId"])
    else:
        for sent_dict in train_instance_dict["sentenceDetails"]:
            sent_ids.append(sent_dict["index"])
    return sent_ids


def pickle_init(allowedWordsLoc):
    global allowed_words
    allowed_words = pickle.load(open(allowedWordsLoc, "rb"))
    # allowed_words = pickle.load(open("allowed_words_dict_1.p", "rb"))

if __name__ == "__main__":
    main_body()
