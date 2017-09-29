# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
from nltk import ngrams, corpus
import string, ntpath
import cPickle as pickle;
import json
from pprint import pprint


# python Print_PositiveInstances.py RNAi ./user_feedback_JSON/Apr4_2017/openaccess_passage_dict_WithFeedback.json
# python Print_PositiveInstances.py RNAi ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/train_passage_dict.json


def main_body():
    parser = argparse.ArgumentParser(prog='Print_PositiveInstances', usage='Print_PositiveInstances.py <label> <psg_dict_loc>', description='Script to print positive instances for a given label from a passage dict')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('psg_dict_loc', help='The location where both passage_dict files - training and OA, are present')

    args = parser.parse_args()
    label = args.label

    fname = ntpath.basename(args.psg_dict_loc).strip()
    if os.path.isfile(args.psg_dict_loc):
        file_handle = open(args.psg_dict_loc, "rb")
        passage_dict = json.load(file_handle)
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + fname + "! Please place it in - ", args.psg_dict_loc
        exit()

    positive_instances = {}
    print "Processing passage file: ", fname
    for pmid in passage_dict:
        positive_instances[pmid] = {}
        positive_instances[pmid][label] = {}
        for prot in passage_dict[pmid][label]["Pos"]:   # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            positive_instances[pmid][label][prot] = []
            for contg_psg_dict in passage_dict[pmid][label]["Pos"][prot]["passageDetails"]:
                positive_instances[pmid][label][prot].extend(contg_psg_dict["textOfInterest"])  # contg_psg_dict["textOfInterest"] is a LIST of phrases

    drname = ntpath.dirname(args.psg_dict_loc)
    output_filename = os.path.join(drname, "Positive_Instances_Feedback_"+label+".json")
    json_fn = open(output_filename, 'w')
    json.dump(positive_instances, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


if __name__ == "__main__":
    main_body()
