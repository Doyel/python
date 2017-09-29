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

# python Generate_Confidences.py <label> ./<vw_Predictions_file> <filename> ./pubmedid_extras.p
# python Generate_Confidences.py RNAi ./../datafiles_BaseLine1/RNAi/models/vw_RNAi_predictions_1_1_L1-01.txt ./../datafiles_BaseLine1/RNAi/models/vw_RNAi_Confidences_1_1_L1-01.txt ./pubmedid_extras.p ./datafiles_BaseLine1/RNAi/vw_RNAi_Test_Aggregated_Confidences.json

def main_body():
    global mydict
    parser = argparse.ArgumentParser(prog='Generate_Confidences.py',
                                     usage='Generate_Confidences.py <label> <vw_prediction_file> <confidences_filename> <DatumKB_dict> [-json <Name of the json file>]',
                                     description='Script to create confidences file')
    parser.add_argument('label', help='The class (in words) we are trying to predict')      
    parser.add_argument('vw_predict_file', help='Prediction file created by vw from the test datafile')
    parser.add_argument('Out_filename', help='Name of the confidences file')   
    parser.add_argument('DatumKB_dict', help='Dictionary created from the datum KB given by Mark')
    parser.add_argument('-json', nargs='?', help='Name of the JSON output file that stores aggregate confidences')
    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()
    
    mydict = pickle.load( open( args.DatumKB_dict, "rb" ) )
    confidence_dict = {}    
    
    generate_prediction_confs(args.vw_predict_file, confidence_dict)
    generate_actuals(args.label, confidence_dict)

    my_json = {}
    conf_fn1 = open(args.Out_filename, 'w'); ctr = 0
    for pmid in confidence_dict.keys():
        my_json[pmid] = {}
        for prot in confidence_dict[pmid]:
            my_json[pmid][prot] = {}
            my_json[pmid][prot]["Distance Label"] = confidence_dict[pmid][prot][1]
            my_json[pmid][prot]["Aggregated Confidence"] = confidence_dict[pmid][prot][0]
            conf_fn1.write(str(confidence_dict[pmid][prot][0]) + "," + confidence_dict[pmid][prot][1] + '\n')
            if confidence_dict[pmid][prot][1] == "1":   ctr += 1
    conf_fn1.close()

    print "Positive instances in the confidences file for label", args.label, "is: ", ctr

    if args.json is not None:
        json_fn = open(args.json, 'w')
        json.dump(my_json, json_fn, indent=4, ensure_ascii=False)
        json_fn.close()

    
def generate_actuals(label, actual_dict):
    temp_dict = {}  # temp_dict is a flattened version of mydict. It only contains those PMIDs that have extras datums of interest and are included in the Test file list.
    create_act_dict(actual_dict.keys(), temp_dict)  # actual_dict.keys() will contain all PMIDs present in the test file list
    for pmid in actual_dict.keys():
        if pmid not in temp_dict: # Article does not have any extras datum that are of interest (i.e. it does not have any 'reqs' and 'do not reqs' extras datums)
            for prot in actual_dict[pmid]:
               actual_dict[pmid][prot].append("-1")
            continue 
        for prot in actual_dict[pmid]: # Article does has extras datum that are of interest to us (i.e. it has 'reqs' and 'do not reqs' extras datums)
            if prot in temp_dict[pmid]: 
                if label.strip() in temp_dict[pmid][prot]:
                    actual_dict[pmid][prot].append("1")
                else:
                    actual_dict[pmid][prot].append("-1")
            else:   # Not all proteins detected in the article are entities of 'reqs' or 'dnreqs' extras datums
                actual_dict[pmid][prot].append("-1")
    
    ctr = 0
    print "Total no of pmids that have extras and are included in the test file list: ", len(temp_dict.keys())
    for pmid in temp_dict.keys():
        for prot in temp_dict[pmid]:
            if label.strip() in temp_dict[pmid][prot]:
                ctr += 1
    # Use this value as the denominator for Recall, not the positives present in the test file
    print "Actual positives in the Datum KB for label", label.strip(), "with respect to the test file: ", ctr
            

def create_act_dict(mypmids, tdict):
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


def generate_prediction_confs(predict_file, pred_dict):
    temp_dict = {}
    create_pred_dict(predict_file, temp_dict)
    for pmid in temp_dict.keys():
        pred_dict[pmid] = {}
        for prot in temp_dict[pmid]:
            if "__" in prot:
                myprot = prot.replace("__", " ")
            else:
                myprot = prot
            pred_dict[pmid][myprot] = [compute_noisyOR(temp_dict[pmid][prot])] 


def create_pred_dict(pfile, tdict):
    with open(pfile, 'rt') as f1:        
        for line in f1:
            line_list = line.split()
            pmid = line_list[1].strip().split("_", 1)[0]
            prots = line_list[1].strip().split("--")[1:]
            prots = list(set(prots))    # prots is a distinct list of proteins extracted from the tag of an instance
            if pmid.strip() in tdict:
                for prt in prots:   # A tag will have multiple proteins only if they are positive
                    if prt.strip() in tdict[pmid.strip()]:
                        tdict[pmid.strip()][prt.strip()].append(float(line_list[0].strip()))
                    else:
                        tdict[pmid.strip()][prt.strip()] = [float(line_list[0].strip())]
            else:
                tdict[pmid.strip()] = {}
                for prt in prots:
                    tdict[pmid.strip()][prt.strip()] = [float(line_list[0].strip())]


def compute_noisyOR(mylist):
    prod = 1.0
    for num in mylist:        
        prod = prod * float(1.0 - num)
    return float(1.0 - prod)


def sigmoid(num):
    return float(1.0 / (1.0 + math.exp(-num)))


if __name__ == "__main__":
    main_body()
