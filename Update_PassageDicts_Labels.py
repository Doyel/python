# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
from nltk import ngrams, corpus
import string
import cPickle as pickle;
import json
from pprint import pprint


# python Update_PassageDicts_Labels.py ./trigger_words_dict.p ./datafiles_BaseLine1 ./datafiles_BaseLine2
# python Update_PassageDicts_Labels.py ./trigger_words_dict.p ./datafiles_BaseLine1 ./datafiles_BaseLine2 -w


def main_body():
    parser = argparse.ArgumentParser(prog='Update_PassageDicts_Labels', usage='Update_PassageDicts_Labels.py <trigg_dict_loc> <psg_dict_loc> <output_dir> [-w]', description='Script to update both passage dicts wrt trigger words')
    parser.add_argument('trigg_dict_loc', help='Location of the dictionary that contains Mark\'s trigger words')
    parser.add_argument('psg_dict_loc', help='The location where both passage_dict files - training and OA, are present')
    parser.add_argument('output_dir', help='The location where the passage_dict files are saved')
    parser.add_argument('-w', action='store_true', help='Switch to indicate whether to update the weighted passagedicts')

    args = parser.parse_args()
    trigger_dict = pickle.load(open(args.trigg_dict_loc, "rb"))

    if args.w:
        listpassagedicts = ["openaccess_passage_dict_Weights.json", "train_passage_dict_Weights.json"]
    else:
        listpassagedicts = ["openaccess_passage_dict.json", "train_passage_dict.json"]


    for passagedictfile in listpassagedicts:

        fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
        if os.path.isfile(fullpath_passagedictfile):
            file_handle = open(fullpath_passagedictfile, "rb")
            passage_dict = json.load(file_handle)
            file_handle.close()
        else:
            print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
            exit()

        print "Processing passage file: ", passagedictfile
        for pmid in passage_dict:
            for label in trigger_dict:
                if label not in passage_dict[pmid]:
                    continue
                for prot in passage_dict[pmid][label]["Pos"]:   # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                    for contg_psg_dict in passage_dict[pmid][label]["Pos"][prot]["passageDetails"]:  
                        found = check_trigger_word(contg_psg_dict["textOfInterest"], trigger_dict[label])
                        if not found:
                            contg_psg_dict["class"] = "-1"  # Hence there may be negative instances inside passage_dict[pmid][label]["Pos"]

        output_filename = os.path.join(args.output_dir, passagedictfile)
        # The below scenario will happen if I give the same source and output dir for the passage dicts
        if (args.output_dir == args.psg_dict_loc) and os.path.isfile(output_filename):
            renamed_passagedictfile = passagedictfile.split(".")[0] + "_OriginalLabels." + passagedictfile.split(".")[1]
            print "Found an earlier / original version of the passage dict file: ", passagedictfile, " in the output dir: ", args.output_dir
            print "Renaming the earlier version to: ", renamed_passagedictfile
            os.rename(output_filename, os.path.join(args.output_dir, renamed_passagedictfile))
        json_fn = open(output_filename, 'w')
        json.dump(passage_dict, json_fn, indent=4, ensure_ascii=False)
        json_fn.close()


def check_trigger_word(sentList, trigger_wordlist):
    for word in trigger_wordlist:
        for sent in sentList:
            idx = sent.find(word)
            if idx != -1:
                #print "Found is true"
                return True
    return False


if __name__ == "__main__":
    main_body()
