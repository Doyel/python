# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import json
import math
from pprint import pprint;
import cPickle as pickle;

# python PredictedConfidences_JSON_Output.py RNAi <Location of vw_Predictions_file> ./../datafiles/SRI_Eval_09152017/datafiles_WithWeights_BaseLine2 ./../datafiles/SRI_Eval_09152017/datafiles_WithWeights_BaseLine2/RNAi/models/vw_RNAi_Individual_Confidences_Test.json -t

def main_body():
    global mydict
    parser = argparse.ArgumentParser(prog='PredictedConfidences_JSON_Output.py', usage='PredictedConfidences_JSON_Output.py <label> <vw_prediction_file> <passage_dict_loc> <confidences_filename> [-t]', description='Script to create confidences file')
    parser.add_argument('label', help='The class (in words) we are trying to predict')
    parser.add_argument('vw_predict_file', help='Prediction file created by vw from the test datafile')
    parser.add_argument('passage_dict_loc', help='Location of the passage dict')
    parser.add_argument('Out_filename', help='Path along with name of the Individual confidences json file')
    parser.add_argument('-t', action='store_true')

    args = parser.parse_args()
    label = args.label

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if label not in allowed_labels:
        print "The given label is not allowed: ", label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    if args.t:
        passagedictfile = "test_passage_dict.json"
    else:
        passagedictfile = "train_passage_dict.json"

    if os.path.isfile(os.path.join(args.passage_dict_loc, passagedictfile)):
        print "Found - " + passagedictfile + "! Everything is perfect in the world!"
        file_handle = open(os.path.join(args.passage_dict_loc, passagedictfile), "rb")
        passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
        file_handle.close()
    else:
        print "The dictionary - " + passagedictfile + " was not found in current directory."
        exit()

    prediction_dict = {}
    tags_passage_dict = {}
    my_json_dict = {}

    create_pred_dict(args.vw_predict_file, prediction_dict)
    create_tag_dict(passage_dict, tags_passage_dict, label)

    for pmid in prediction_dict:
        my_json_dict[pmid] = {}
        for prot in prediction_dict[pmid]:
            if check_for_unicode(prot):
                prot = prot.decode("utf-8")
            my_json_dict[pmid][prot] = []
            for tag in prediction_dict[pmid][prot]:
                temp_dict = {}
                temp_dict["Tag"] = tag
                temp_dict["Distance Label"] = tags_passage_dict[pmid][prot][tag]["Distance Label"]
                temp_dict["Passage"] = tags_passage_dict[pmid][prot][tag]["Passage"]
                temp_dict["Predicted Confidence"] = prediction_dict[pmid][prot][tag]
                temp_dict["Index List"] = tags_passage_dict[pmid][prot][tag]["Index List"]
                my_json_dict[pmid][prot].append(temp_dict)

    json_fn = open(args.Out_filename, 'w')
    json.dump(my_json_dict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


def create_pred_dict(pfile, tdict):
    with open(pfile, 'rt') as f1:
        for line in f1:
            if check_for_unicode(line):
                line = line.decode("utf-8")
            line_list = line.split()
            tag = line_list[1].strip()
            pmid = line_list[1].strip().split("_", 1)[0]
            prots = line_list[1].strip().split("--")[1:]
            prots = list(set(prots))    # prots is a distinct list of proteins extracted from the tag of an instance
            if pmid.strip() not in tdict:
                tdict[pmid.strip()] = {}

            for prt in prots:   # A tag will have multiple proteins only if they are positive
                prt = prt.replace("__", " ")
                if prt.strip() not in tdict[pmid.strip()]:
                    tdict[pmid.strip()][prt.strip()] = {}

                if tag in tdict[pmid.strip()][prt.strip()]:
                    print "Tag is present twice in the test data file!"
                    exit()
                tdict[pmid.strip()][prt.strip()][tag] = float(line_list[0].strip())


def create_tag_dict(passage_dict, tags_passage_dict, label):
    for pmid in passage_dict:
        tags_passage_dict[pmid] = {}
        for myclass in passage_dict[pmid][label]:               # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            for prot in passage_dict[pmid][label][myclass]:     # passage_dict[pmid][label][myclass] is a dict
                if prot not in tags_passage_dict[pmid]:
                    tags_passage_dict[pmid][prot] = {}
                for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:     # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                    if contg_psg_dict["tag"] in tags_passage_dict[pmid][prot]:                       # contg_psg_dict represents one training instance
                        print "Tag already present in dict: ", contg_psg_dict["tag"]
                        exit()
                    tags_passage_dict[pmid][prot][contg_psg_dict["tag"]] = {}
                    tags_passage_dict[pmid][prot][contg_psg_dict["tag"]]["Passage"] = contg_psg_dict["textOfInterest"]
                    tags_passage_dict[pmid][prot][contg_psg_dict["tag"]]["Distance Label"] = contg_psg_dict["class"]
                    tags_passage_dict[pmid][prot][contg_psg_dict["tag"]]["Index List"] = getIndexList(contg_psg_dict)


def check_for_unicode(mystring):
    try:
        mystring.decode('ascii')
    except UnicodeDecodeError:
        #print mystring
        return True
    else:
        return False


def getIndexList(my_psg_dict):
    mylist = []
    for sentence in my_psg_dict["sentenceDetails"]:
        mylist.append(sentence["index"])
    return mylist


if __name__ == "__main__":
    main_body()
