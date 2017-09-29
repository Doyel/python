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


# python Update_PassageDicts_Weights.py RNAi 1 ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/datafiles_WithWeights_Baseline1
# python Update_PassageDicts_Weights.py RNAi 5 ./user_feedback_JSON/Apr4_2017 -f


def main_body():
    parser = argparse.ArgumentParser(prog='Update_PassageDicts_Weights', usage='Update_PassageDicts_Weights.py <psg_dict_loc>', description='Script to update both passage dicts with importance weights')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) to be processed')
    parser.add_argument('tot_weight', help='The total weight assigned to a prot-article pair')
    parser.add_argument('psg_dict_loc', help='The location where both passage_dict files - training and OA, are present')
    parser.add_argument('-f', action='store_true')

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    label = args.label
    total_wt = float(args.tot_weight)

    if args.f:
        listpassagedicts = ["openaccess_passage_dict_WithFeedback.json"]
    else:
        listpassagedicts = ["openaccess_passage_dict.json", "train_passage_dict.json"]

    temp_dict = {}
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
            if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
                continue
            temp_dict[pmid] = {}
            for myclass in passage_dict[pmid][label]:
                for prot in passage_dict[pmid][label][myclass]:   # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                    if prot not in temp_dict[pmid]:     # The same protein can have pos as well as neg instances after user gives feedback
                        temp_dict[pmid][prot] = {}
                        temp_dict[pmid][prot]["numInstances"] = 0

                    temp_dict[pmid][prot]["numInstances"] += len(passage_dict[pmid][label][myclass][prot]["passageDetails"])
                    if myclass == "Pos":
                        temp_dict[pmid][prot]["Pos"] = len(passage_dict[pmid][label][myclass][prot]["passageDetails"])
                    else:
                        temp_dict[pmid][prot]["Neg"] = len(passage_dict[pmid][label][myclass][prot]["passageDetails"])

        for pmid in temp_dict:
            for prot in temp_dict[pmid]:
                temp_dict[pmid][prot]["weight"] = format(total_wt / float(temp_dict[pmid][prot]["numInstances"]), '.2f')

        #pprint(temp_dict); exit()

        for pmid in passage_dict:
            for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label][myclass] is a dict
                    for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                        contg_psg_dict["weight"] = temp_dict[pmid][prot]["weight"]

        print "Making a copy of the original passage dict file: ", passagedictfile, "\thaving the phrase \"_NoWeights\""
        fullpath_renamed_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile.split(".")[0] + "_NoWeights." + passagedictfile.split(".")[1])
        os.rename(fullpath_passagedictfile, fullpath_renamed_passagedictfile)
        #output_filename = os.path.join(args.psg_dict_loc, passagedictfile.split(".")[0] + "_Weights." + passagedictfile.split(".")[1])
        print "Re-creating original passage dict file that contains weights: ", passagedictfile
        json_fn = open(fullpath_passagedictfile, 'w')
        json.dump(passage_dict, json_fn, indent=4, ensure_ascii=False)
        json_fn.close()

    if args.f:
        protpairs_fname = "PMID-Prot_Pairs_Feedback.json"
    else:
        protpairs_fname = "PMID-Prot_Pairs.json"
    output_filename = os.path.join(args.psg_dict_loc, protpairs_fname)
    json_fn = open(output_filename, 'w')
    json.dump(temp_dict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


if __name__ == "__main__":
    main_body()
