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


# python Collect_Mark_YellowHighlights.py ./user_feedback_JSON/Apr4_2017/user_edits_BIG_MECHANISM.json  Mark_Yellow_Highlights_Extra_Apr11.json


def main_body():
    global pruned_words_list
    parser = argparse.ArgumentParser(prog='Collect_Mark_YellowHighlights', usage='Collect_Mark_YellowHighlights.py <mongoimportedjson> <outfile>', description='Script to collect all yellow highlights given by Mark')
    parser.add_argument('mongoimportedjson', help='The name of the user feedback JSON file')
    parser.add_argument('outfile', help='Name of the file in which the yellow highlights will be stored')
    args = parser.parse_args()

    processlabel = "RNAi"
    #PMIDS_toProcess = ["21629263", "22833096", "16009723", "18583988", "18411307", "20890305", "18504304", "19234442", "22245064", "20141835", "21573184", "20547488", "19783983", "22433566", "14517278", "20026654"]
    PMIDS_toProcess = ["24106086", "20190815"]
    drname = ntpath.dirname(args.mongoimportedjson)

    file_handle = open(args.mongoimportedjson, "rb")
    myfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
    file_handle.close()
    for mydict in myfile_data:
        mydict.pop("_id", None)

    mongo_passage_dict = {}     # Something that Matt was supposed to give me
    for mongo_doc in myfile_data:      # mongo_doc is a single document that corresponds to an unique PMID-label-protein triple
        if len(mongo_doc.keys()) > 1:
            print "Current document has multiple PMIDs. Please Check!"
            pprint(mongo_doc)
            exit()
        for pmid in mongo_doc:
            if pmid in PMIDS_toProcess:
                print "Processing: ", pmid

                if len(mongo_doc[pmid].keys()) > 1:
                    print "Current document has multiple labels. Please Check!"
                    pprint(mongo_doc[pmid])
                    exit()

                if mongo_doc[pmid].keys()[0] == processlabel:   # Only consider dicts with labels of my interest!
                    for label in mongo_doc[pmid]:
                        if len(mongo_doc[pmid][label].keys()) > 1:
                            print "Current document has multiple proteins. Please Check!"
                            pprint(mongo_doc[pmid][label])
                            exit()

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
                                exit()

                            mongo_passage_dict[pmid][label][protein]["YellowHighlights"] = []
                            for hilite in mongo_doc[pmid][label][protein]["HighlightObject"]:
                                if "rgb(255, 255, 123)" in hilite[0]:
                                    mongo_passage_dict[pmid][label][protein]["YellowHighlights"].append(hilite[1])

                            mongo_passage_dict[pmid][label][protein]["ReadableTime"] = datetime.datetime.utcfromtimestamp(document_ts / 1000).strftime('%m-%d-%Y %H:%M:%S')

    write_fname = os.path.join(drname, args.outfile)
    json_fn = open(write_fname, 'w')
    json.dump(mongo_passage_dict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


def getTimestamp(document_protdict):
    return int(document_protdict["Timestamp"]["$numberLong"])


if __name__ == "__main__":
    main_body()

