# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import math
from pprint import pprint;
from itertools import izip
import cPickle as pickle;


# python Generate_Confidences_BaseLine2_TestPassage.py RNAi ./datafiles_BaseLine2/RNAi/models/vw_RNAi_predictions_1_1.txt ./datafiles_BaseLine2/RNAi/vw_RNAi_Test_File.txt ./datafiles_BaseLine2/RNAi/models/vw_RNAi_Confidences_1_1.txt


def main_body():
    global confidence_dict
    parser = argparse.ArgumentParser(prog='Generate_Confidences_BaseLine2_TestPassage.py', usage='Generate_Confidences_BaseLine2_TestPassage.py <label> <vw_prediction_file> <loc_test_datafile> <confidences_filename>', description='Script to create confidences file')
    parser.add_argument('label', help='The class (in words) we are trying to predict')
    parser.add_argument('vw_predict_file', help='Prediction file created by vw from the test datafile')
    parser.add_argument('loc_test_datafile', help='Test datafile that contains the actual classes for each test instance')
    parser.add_argument('Out_filename', help='Name of the confidences file')

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    confidence_dict = {}
    confidence_dict["1"] = {}
    confidence_dict["-1"] = {}

    with open(args.vw_predict_file) as predfile, open(args.loc_test_datafile) as testfile:
        for pred, test in izip(predfile, testfile):
            predline = pred.strip(); testline = test.strip()
            test_label = testline.split()[0].strip()
            line_list = predline.split()
            pmid = line_list[1].strip().split("_", 1)[0].strip()
            prots = line_list[1].strip().split("--")[1:]
            populate_confidence_dict(test_label, pmid, list(set(prots)), float(line_list[0].strip()))

    test_pmids = list(set(confidence_dict["1"].keys() + confidence_dict["-1"].keys()))
    DatumKBdict = pickle.load(open("pubmedid_extras.p", "rb"))
    test_DatumKBdict = {}  # temp_dict is a flattened version of mydict.
    create_act_dict(test_pmids, test_DatumKBdict, DatumKBdict)   # It only contains those PMIDs that have extras datums of interest and are included in the Test file list.

    ctr = 0
    for pmid in test_DatumKBdict.keys():
        for prot in test_DatumKBdict[pmid]:
            if args.label.strip() in test_DatumKBdict[pmid][prot]:
                ctr += 1
    # Use this value as the denominator for Recall, not the positives present in the test file
    print "Actual positives in the Datum KB for label", args.label.strip(), "with respect to the test file: ", ctr

    records = []
    populate_confidence_file(records, "1")
    populate_confidence_file(records, "-1")

    conf_fn1 = open(args.Out_filename, 'w')
    for rec in records:
        conf_fn1.write(rec + '\n')
    conf_fn1.close()


def populate_confidence_dict(label, pmid, prots, conf_value):
    global confidence_dict
    if pmid.strip() in confidence_dict[label]:
        for prt in prots:
            if prt.strip() in confidence_dict[label][pmid.strip()]:
                confidence_dict[label][pmid.strip()][prt.strip()].append(conf_value)
            else:
                confidence_dict[label][pmid.strip()][prt.strip()] = [conf_value]
    else:
        confidence_dict[label][pmid.strip()] = {}
        for prt in prots:
            confidence_dict[label][pmid.strip()][prt.strip()] = [conf_value]


def populate_confidence_file(records, label):
    for pmid in confidence_dict[label]:
        for prot in confidence_dict[label][pmid]:
            noisyOR_Conf = compute_noisyOR(confidence_dict[label][pmid][prot])
            actual_label = label
            records.append(str(noisyOR_Conf) + "," + actual_label)


def create_act_dict(mypmids, tdict, mydict):
    print "Total no of pmids in the test file list: ", len(mypmids)
    for pmid in mypmids:
        if pmid in mydict:
            tdict[pmid] = {}          # tdict contains only those PMIDs that exist both in DatumKB and the test file
            for ent in mydict[pmid]:
                mylist = [];
                if len(mydict[pmid][ent]['reqstest']) > 0:
                    mylist.append("reqs"); mylist.extend(mydict[pmid][ent]["reqstest"])
                if len(mydict[pmid][ent]['dnreqstest']) > 0:
                    mylist.append("dnreqs"); mylist.extend(mydict[pmid][ent]["dnreqstest"])
                tdict[pmid][ent] = list(set(mylist))


def compute_noisyOR(mylist):
    prod = 1.0
    for num in mylist:
        prod = prod * float(1.0 - num)
    return float(1.0 - prod)


if __name__ == "__main__":
    main_body()