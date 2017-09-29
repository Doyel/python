# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import os
import ntpath
import json
import sys
reload(sys)
sys.setdefaultencoding("utf-8")
from random import shuffle
from Create_MLFile import create_lexical_features
from ParamTuning_InternalCV import writefile

# First tune weights, then tune parameters!!!
# python Create_Weighted_Training_File.py RNAi ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/datafiles_Baseline2 1 1
# python Create_Weighted_Training_File.py RNAi ./datafiles_TriggerHighlights_BaseLine2 1 1
# python Create_Weighted_Training_File.py KO ./datafiles_TriggerHighlights_BaseLine2 1 1
# python Create_Weighted_Training_File.py RNAi ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/datafiles_Baseline2 1 5 --feedback ./user_feedback_JSON/Apr4_2017/


def main_body():
    parser = argparse.ArgumentParser(prog='Create_Weighted_Training_File', usage='Create_Weighted_Training_File.py <label> <passageDict_dir> <num_folds> <optim_method> [--feedback <Path to feedback passage dict>]', description='Script to tune instance weights')
    parser.add_argument('label', help='Predicate being analyzed')
    parser.add_argument('passageDict_dir', help='Dir where all the original passage dicts are located')
    parser.add_argument('NonOA_Wgt', help='Weight to be distributed over a prot-PMID pair for a Non OpenAccess article')
    parser.add_argument('OA_Wgt', help='Weight to be distributed over a prot-PMID pair for an OpenAccess article')
    parser.add_argument('-f', '--feedback', help='Path to the feedback passage dict json file')  # Optional arg

    args = parser.parse_args()
    label = args.label
    listpassagedicts = [os.path.join(args.passageDict_dir, "train_passage_dict.json")]

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    if args.feedback is None:
        listpassagedicts.append(os.path.join(args.passageDict_dir, "openaccess_passage_dict.json"))
    else:
        listpassagedicts.append(os.path.join(args.feedback, "openaccess_passage_dict_WithFeedback.json"))

    for i, elem in enumerate(listpassagedicts):
        if os.path.isfile(elem):
            print "Found - " + ntpath.basename(elem)
            file_handle = open(elem, "rb")
            if i == 0:
                train_passage_dict = json.load(file_handle)
            else:
                OA_passage_dict = json.load(file_handle)
            file_handle.close()
        else:
            print "The passage dict file - " + ntpath.basename(elem) + " was not found in " + ntpath.dirname(elem)
            exit()

    print "Creating optimally weighted Training file"
    optimal_NoOA_wt = int(args.NonOA_Wgt)
    optimal_OA_wt = int(args.OA_Wgt)
    update_passagedict(label, train_passage_dict, optimal_NoOA_wt)
    train_instances = create_instances(label, train_passage_dict)
    update_passagedict(label, OA_passage_dict, optimal_OA_wt)
    OAtrain_instances = create_instances(label, OA_passage_dict)
    weight_str = getcomb_wgtstr(int(optimal_NoOA_wt), int(optimal_OA_wt))

    train_records = train_instances + OAtrain_instances
    shuffle(train_records)
    if args.feedback is None:
        trainfilename = "vw_" + label + "_Weights_" + weight_str + "_Training_File.txt"
        writefile(os.path.join(args.passageDict_dir, label, trainfilename), train_records)
        print "Training File created at: ", os.path.join(args.passageDict_dir, label)
    else:
        trainfilename = "vw_" + label + "_Weights_" + weight_str + "_WithFeedback_Training_File.txt"
        writefile(os.path.join(args.feedback, label, trainfilename), train_records)
        print "Training File created at: ", os.path.join(args.feedback, label)


def getcomb_wgtstr(weight1, weight2):
    mystr1 = str(weight1)
    mystr2 = str(weight2)
    return mystr1 + "_" + mystr2


def update_passagedict(label, passage_dict, total_wt):
    temp_dict = {}
    for pmid in passage_dict:
        if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
            continue
        temp_dict[pmid] = {}
        for myclass in passage_dict[pmid][label]:
            for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                if prot not in temp_dict[pmid]:  # The same protein can have pos as well as neg instances after user gives feedback
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

    for pmid in passage_dict:
        if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
            continue
        for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label][myclass] is a dict
                for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                    contg_psg_dict["weight"] = temp_dict[pmid][prot]["weight"]


def create_instances(label, passage_dict):
    myrecords = []
    for pmid in passage_dict:
        if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
            continue
        for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label][myclass] is a dict
                for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                    vw_label = contg_psg_dict["class"]  # contg_psg_dict represents one training instance
                    vw_tag = contg_psg_dict["tag"]
                    lex_feat_str = create_lexical_features(contg_psg_dict["textOfInterest"])
                    if "weight" in contg_psg_dict:
                        rec = vw_label + " " + contg_psg_dict["weight"] + " " + vw_tag + "|LexicalFeatures " + lex_feat_str  # + " |ProteinFeature " + prot
                    else:
                        rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str  # + " |ProteinFeature " + prot
                    myrecords.append(rec)

    return myrecords


if __name__ == "__main__":
    main_body()