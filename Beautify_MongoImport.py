# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
import json
import copy
import ntpath
import datetime
from pprint import pprint;


# python Beautify_MongoImport.py ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0 ./user_feedback_JSON/Apr4_2017/user_edits.json
# python Beautify_MongoImport.py ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0 ./user_feedback_JSON/Apr4_2017/articles.json
# python Beautify_MongoImport.py ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0 ./user_feedback_JSON/Jan18_2017/user_edits_incremental.json

# openaccess_passage_dict.json contains the initial highlights.
# user_edits.json contains the most recent user changes to the highlights for any unique PMID-label-protein triple that were ever given feedback upon, since the beginning of time


def main_body():
    global pruned_words_list
    parser = argparse.ArgumentParser(prog='Beautify_MongoImport', usage='Beautify_MongoImport.py <psg_dict_loc> <mongoimportedjson>', description='Script to beautify the specified mongoimported JSON file')
    parser.add_argument('psg_dict_loc', help='The location where the OpenAccess passage_dict json file is present')
    parser.add_argument('mongoimportedjson', help='The name of the user feedback JSON file')
    args = parser.parse_args()

    # ["24106086", "20190815", "21629263",      "22833096", "16009723", "18583988",     "14517278", "19783983", "18411307",     "20890305", "18504304", "19234442",     "22245064", "20141835", "21573184"]
    ignore_PMIDS = []       #["20547488", "22433566", "20026654"]

    passagedictfile = "openaccess_passage_dict.json"
    fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
    if os.path.isfile(fullpath_passagedictfile):
        print "Found the dictionary - " + passagedictfile + " Everything is perfect in the world!"
        file_handle = open(fullpath_passagedictfile, "rb")
        passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
        exit()

    file_handle = open(args.mongoimportedjson, "rb")
    myfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
    file_handle.close()
    for mydict in myfile_data:
        mydict.pop("_id", None)

    # for line in file_handle:
    #     mydict = dict(json.loads(line.strip()))
    #     mydict.pop("_id", None)
    #     myfile_data.append(mydict)
    # file_handle.close()

    fname = ntpath.basename(args.mongoimportedjson).strip()
    drname = ntpath.dirname(args.mongoimportedjson)
    write_fname = os.path.join(drname,fname.split(".")[0]+"_Beautified."+fname.split(".")[1])
    silentremove(write_fname)
    json_fn = open(write_fname, 'w')
    json.dump(myfile_data, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()

    if fname != "user_edits.json" and fname != "user_edits_BIG_MECHANISM.json" and fname != "user_edits_TMP.json" and (not fname.startswith("user_edits_DUPL")):
        print "Filename is not one of user_edits.json, user_edits_BIG_MECHANISM.json, user_edits_TMP.json or user_edits_DUPL.json. Hence exiting at this point!"
        exit()

    mongo_passage_dict = {}     # Something that Matt was supposed to give me
    for mongo_doc in myfile_data:      # mongo_doc is a single document that corresponds to an unique PMID-label-protein triple
        if len(mongo_doc.keys()) > 1:
            print "Current document has multiple PMIDs. Please Check!"
            pprint(mongo_doc)
            exit()
        for pmid in mongo_doc:
            if pmid in ignore_PMIDS:
                print "Ignoring PMID: ", pmid, ". Please remove from ignore list if you want to collect userfeedback for this article!"
                continue
            if len(mongo_doc[pmid].keys()) > 1:
                print "Current document has multiple labels. Please Check!"
                pprint(mongo_doc[pmid])
                exit()
            for label in mongo_doc[pmid]:
                if len(mongo_doc[pmid][label].keys()) > 1:
                    print "Current document has multiple proteins. Please Check!"
                    pprint(mongo_doc[pmid][label])
                    exit()

                #if pmid == "18411307" and label == "RNAi" and mongo_doc[pmid][label].keys()[0] == "mdc1":
                #    print "#------------------------------------------------------------------#"
                #    print "Ignoring the triple - PMID: ", pmid, " Label: ", label, " Protein: ", mongo_doc[pmid][label].keys()[0]
                #    print "#------------------------------------------------------------------#"
                #    continue    # mongo_doc corresponds to an UNIQUE PMID-label-protein triple or else each triple having the same values for PMID-LABEL-PROTEIN gets ignored

                if pmid not in mongo_passage_dict:
                    mongo_passage_dict[pmid] = {}
                if label not in mongo_passage_dict[pmid]:
                    mongo_passage_dict[pmid][label] = {}

                for protein in mongo_doc[pmid][label]:
                    document_ts = getTimestamp(mongo_doc[pmid][label][protein])
                    if protein not in mongo_passage_dict[pmid][label]:
                        mongo_passage_dict[pmid][label][protein] = {}
                        mongo_passage_dict[pmid][label][protein]["Timestamp"] = document_ts
                    else:
                        print "There are multiple documents having the same PMID-label-protein triple"
                        print "PMID: ", pmid, " Label: ", label, " Protein: ", protein
                        print "My assumption that each mongo_doc corresponds to an unique PMID-label-protein triple is VIOLATED!"
                    if mongo_passage_dict[pmid][label][protein]["Timestamp"] <= document_ts:
                        # Assumption: All list members of the "HighlightObject" list object have yellow highlight color
                        mongo_passage_dict[pmid][label][protein]["HighlightObject"] = copy.deepcopy(mongo_doc[pmid][label][protein]["HighlightObject"])
                        mongo_passage_dict[pmid][label][protein]["ReadableTime"] = datetime.datetime.utcfromtimestamp(document_ts / 1000).strftime('%m-%d-%Y %H:%M:%S')

    # Update openaccess_passage_dict using mongo_passage_dict. openaccess_passage_dict is the original passage dict that contains the initial highlights.
    # mongo_passage_dict contains the most recent user feedback for any PMID-label-protein triple, that I do not want to ignore!
    for pmid in mongo_passage_dict:
        for label in mongo_passage_dict[pmid]:
            for protein in mongo_passage_dict[pmid][label]:
                if passage_dict[pmid][label]["Pos"][protein]["Timestamp"] < mongo_passage_dict[pmid][label][protein]["Timestamp"]:
                    passage_dict[pmid][label]["Pos"][protein]["HighlightObject"] = copy.deepcopy(mongo_passage_dict[pmid][label][protein]["HighlightObject"])
                    passage_dict[pmid][label]["Pos"][protein]["Timestamp"] = mongo_passage_dict[pmid][label][protein]["Timestamp"]

    write_fname = os.path.join(drname, "openaccess_passage_dict_WithFeedback.json")
    json_fn = open(write_fname, 'w')
    json.dump(passage_dict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()

    write_fname = os.path.join(drname, "SavedHighlights_WithFeedback.json")
    json_fn = open(write_fname, 'w')
    json.dump(mongo_passage_dict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured


def getTimestamp(document_protdict):
    return int(document_protdict["Timestamp"]["$numberLong"])


if __name__ == "__main__":
    main_body()
