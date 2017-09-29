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


# python Update_PassageDicts_Weights.py 1 ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/datafiles_WithWeights_Baseline1 "RNAi|KO|omission"
# python Update_PassageDicts_Weights.py 5 ./user_feedback_JSON/Apr4_2017 RNAi|KO|omission -f


def main_body():
    parser = argparse.ArgumentParser(prog='Update_PassageDicts_Weights', usage='Update_PassageDicts_Weights.py <weight> <psg_dict_loc> <label_string>', description='Script to update both passage dicts with importance weights')
    parser.add_argument('tot_weight', help='The total weight assigned to a prot-article pair')
    parser.add_argument('psg_dict_loc', help='The location where both passage_dict files - training and OA, are present')
    parser.add_argument('label_string', help='The labels that need to be processed, provided as a | delimeted string')
    parser.add_argument('-f', action='store_true')

    args = parser.parse_args()
    total_wt = float(args.tot_weight)
    user_labels = str(args.label_string).split("|")

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    for mylabel in user_labels:
        if mylabel not in allowed_labels:
            print "The given label is not allowed: ", mylabel, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
            print "Please try again!"
            exit()

    if args.f:
        passagedictfile = "openaccess_passage_dict_WithFeedback.json"
        fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
        if os.path.isfile(fullpath_passagedictfile):
            file_handle = open(fullpath_passagedictfile, "rb")
            FB_passage_dict = json.load(file_handle)
            file_handle.close()
            print "Found the dictionary - ", passagedictfile
        else:
            print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
            exit()
        listpassagedicts = [FB_passage_dict]
        listpassagedictsfiles = ["openaccess_passage_dict_WithFeedback.json"]
    else:
        passagedictfile = "openaccess_passage_dict.json"
        fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
        if os.path.isfile(fullpath_passagedictfile):
            file_handle = open(fullpath_passagedictfile, "rb")
            OA_passage_dict = json.load(file_handle)
            file_handle.close()
            print "Found the dictionary - ", passagedictfile
        else:
            print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
            exit()

        passagedictfile = "train_passage_dict.json"
        fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
        if os.path.isfile(fullpath_passagedictfile):
            file_handle = open(fullpath_passagedictfile, "rb")
            NonOA_passage_dict = json.load(file_handle)
            file_handle.close()
            print "Found the dictionary - ", passagedictfile
        else:
            print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
            exit()
        listpassagedicts = [OA_passage_dict, NonOA_passage_dict]
        listpassagedictsfiles = ["openaccess_passage_dict.json", "train_passage_dict.json"]

    temp_dict = {}
    for label in user_labels:
        print "Label being processed: ", label
        temp_dict[label] = {}
        for passage_dict in listpassagedicts:
            for pmid in passage_dict:
                if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
                    continue
                temp_dict[label][pmid] = {}     # Both passage dicts have a different set of PMIDs
                for myclass in passage_dict[pmid][label]:   # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                    for prot in passage_dict[pmid][label][myclass]:
                        if prot not in temp_dict[label][pmid]:     # The same protein can have pos as well as neg instances after user gives feedback
                            temp_dict[label][pmid][prot] = {}
                            temp_dict[label][pmid][prot]["numInstances"] = 0
                        else:
                            if not args.f:
                                print "For Baseline models, control should not come here!"
                                print "In Baseline models, a protein can either be a pos protein or a neg protein but NOT both!"
                                exit()

                        temp_dict[label][pmid][prot]["numInstances"] += len(passage_dict[pmid][label][myclass][prot]["passageDetails"])
                        if myclass == "Pos":
                            temp_dict[label][pmid][prot]["Pos"] = len(passage_dict[pmid][label][myclass][prot]["passageDetails"])
                        else:
                            temp_dict[label][pmid][prot]["Neg"] = len(passage_dict[pmid][label][myclass][prot]["passageDetails"])

            for pmid in passage_dict:
                if pmid not in temp_dict[label]:
                    continue
                for prot in temp_dict[label][pmid]:
                    temp_dict[label][pmid][prot]["weight"] = format(total_wt / float(temp_dict[label][pmid][prot]["numInstances"]), '.2f')

            #pprint(temp_dict); exit()

            for pmid in passage_dict:
                if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
                    continue
                for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                    for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label][myclass] is a dict
                        for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                            if "weight" in contg_psg_dict:
                                print "The key \"weight\" is already present for this instance of the protein-article pair for predicate: ", label
                                print "Please check!"; exit()
                            contg_psg_dict["weight"] = temp_dict[label][pmid][prot]["weight"]

    for i, passagedictfile in enumerate(listpassagedictsfiles):
        fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
        print "Making a copy of the original passage dict file: ", passagedictfile, "\thaving the phrase \"_NoWeights\""
        fullpath_renamed_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile.split(".")[0] + "_NoWeights." + passagedictfile.split(".")[1])
        os.rename(fullpath_passagedictfile, fullpath_renamed_passagedictfile)
        #output_filename = os.path.join(args.psg_dict_loc, passagedictfile.split(".")[0] + "_Weights." + passagedictfile.split(".")[1])
        print "Re-creating original passage dict file containing weights: ", passagedictfile
        json_fn = open(fullpath_passagedictfile, 'w')
        json.dump(listpassagedicts[i], json_fn, indent=4, ensure_ascii=False)
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
