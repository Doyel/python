# -*- coding: utf-8 -*-
__author__ = 'Doyel'


import json;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import re;
import argparse;
import os;
import ntpath
import cPickle as pickle;
from pprint import pprint;
from Utilities import read_config_file

# python Create_MongoDBFile.py ./user_feedback_JSON/Jan18_2017/openaccess_passage_dict_WithFeedback.json Datums_MongoDB_OpenAccess_WithFeedback.json
# python Create_MongoDBFile.py ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/openaccess_passage_dict.json Datums_MongoDB_OpenAccess.json
# python Create_MongoDBFile.py ./datafiles_TriggerHighlights_BaseLine1/openaccess_passage_dict.json Datums_MongoDB_OpenAccess_WithTriggers.json

def __main__():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='JSON Converter', usage='Create_MongoDBFile.py <psg_dict_loc> <output_fname>', description='Script that reads passage_dict.json and creates a MongoDB JSON file')
    parser.add_argument('psg_dict_loc', help='The location where the OpenAccess passage_dict json file is present')
    parser.add_argument('output_fname', help='The name of the MongoDB database file that will be generated')
    args = parser.parse_args()

    # allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    # labels_toProcess = ["RNAi", "reqs"]
    labels_toProcess = ["RNAi", "KO", "omission", "reqs"]
    print "The MongoDB database file will be generated for the below labels!"
    print ", ".join(labels_toProcess)

    passagedictfile = ntpath.basename(args.psg_dict_loc)
    drname = ntpath.dirname(args.psg_dict_loc)
    fullpath_passagedictfile = os.path.join(drname, passagedictfile)
    if os.path.isfile(fullpath_passagedictfile):
        print "Found the dictionary - " + passagedictfile + " Everything is perfect in the world!"
        file_handle = open(fullpath_passagedictfile, "rb")
        passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
        exit()

    summary_dict = pickle.load( open( "summary_dict_Extras_OpenAccess.p", "rb" ) )
    dict_nameToid = pickle.load(open('dictionary_nameToUniProtId.p', 'rb'))

    test_labels = ['RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    type_labels = ['reqs', 'dnreqs']

    grand_list = []
    for pmid in passage_dict:
        curr_pmid_biblio = summary_dict[pmid]
        mydict = {}; mydict["Datums"] = []; datum_ctr = 0
        for label in labels_toProcess:
            if label not in passage_dict[pmid]: continue
            for prot in passage_dict[pmid][label]["Pos"]:
                mydatum = {}
                mydatum["confidence"] = 0.0
                mydatum["datum_id"] = pmid + "_" + str(datum_ctr)
                mydatum["timestamp"] = passage_dict[pmid][label]["Pos"][prot]["Timestamp"]
                mydatum["map"] = {}
                mydatum["map"]["TreatmentEntity"] = {}
                mydatum["map"]["TreatmentEntity"]["Entity"] = {}
                mydatum["map"]["TreatmentEntity"]["Entity"]["strings"] = prot
                mydatum["map"]["TreatmentEntity"]["Entity"]["UniprotId"] = dict_nameToid[prot]
                if label in test_labels:
                    mydatum["map"]["TreatmentTest"] = {}
                    mydatum["map"]["TreatmentTest"]["Highlight"] = passage_dict[pmid][label]["Pos"][prot]["HighlightObject"]
                    mydatum["map"]["TreatmentTest"]["Text"] = label
                elif label in type_labels:
                    mydatum["map"]["TreatmentType"] = {}
                    mydatum["map"]["TreatmentType"]["Highlight"] = passage_dict[pmid][label]["Pos"][prot]["HighlightObject"]
                    mydatum["map"]["TreatmentType"]["Text"] = label

                datum_ctr += 1
                mydict["Datums"].append(mydatum)    # There will be a datum for each unique protein-label pair in an article

        if len(mydict["Datums"]) == 0: continue
        mydict.update(curr_pmid_biblio)
        grand_list.append(mydict)

    mongodb_datafile = os.path.join(json_location, args.output_fname)
    json_fn = open(mongodb_datafile, 'w')
    json.dump(grand_list, json_fn, indent=4, ensure_ascii=False)


__main__();
