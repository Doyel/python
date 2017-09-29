# -*- coding: utf-8 -*-
__author__ = 'Doyel'


import json;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import re;
import argparse;
import os;
import cPickle as pickle;
from pprint import pprint;

def __main__():
    parser = argparse.ArgumentParser(prog='JSON Converter', usage='Convert_SRIfiles_JSON.py <dirpath>', description='Script that reads pre-processed JSON files sent by SRI and creates a single JSON file')
    parser.add_argument('dir_path', help='Path where the pre-processed JSON files are located')
    args = parser.parse_args()

    summary_dict = pickle.load( open( "summary_dict.p", "rb" ) )

    yes_subj_uniprot = 0; no_subj_uniprot = 0
    yes_trtmnt_uniprot = 0; no_trtmnt_uniprot = 0
    grand_list = []
    for myfile in os.listdir(args.dir_path):
        if myfile.startswith("."): continue     # Skip over hidden files        
        current_file = os.path.join(args.dir_path, myfile)
        file_handle = open(current_file, "rb")
        myfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
        file_handle.close();

        if len(myfile_data) == 0:   # If the file has no datums in it then continue to the next json file
            print "Skipping file: ", myfile; continue

        pmcid = myfile.split(".")[0]
        #print "Processing PMCID: ", pmcid
        curr_pmcid_biblio = summary_dict[pmcid]
        
        evid_dict = {}; ctr = 1;
        for elem_dict in myfile_data:
            
            #################--------------- The input to this script will be pre-processed JSON files ----------------###############

            #if len(elem_dict["evidence"]) > 1:  print "evidence has more than 1 element in its list"
            #if elem_dict["evidence"][0] not in evid_dict:
            #    evid_dict[elem_dict["evidence"][0]] = pmcid + "_" + str(ctr); ctr+=1;

            if "subject" in elem_dict["map"]:   # elem_dict["map"] is a dictionary. 'subject' is a key of elem_dict["map"]
                for subj in elem_dict["map"]["subject"]: # elem_dict["map"]["subject"] is a list of dictionaries. subj is a dictionary
                    if len(subj["Entity"]["strings"]) > 1:  print "Subject: ", subj["Entity"]["strings"]    # This print statement does not print anything
                    subj["Entity"]["strings"] = subj["Entity"]["strings"][0]    # Converting subj["Entity"]["strings"] from a list that always contains a single element to a string.
                    if "uniprotSym" in subj["Entity"]:     # subj["Entity"] is a dictionary
                        #print "Subject: ", subj["Entity"]["strings"], "\t UniprotId: ", subj["Entity"]["uniprotSym"]
                        yes_subj_uniprot += 1
                    else:
                        no_subj_uniprot += 1

            if "treatment" in elem_dict["map"]: # elem_dict["map"] is a dictionary. 'treatment' is a key of elem_dict["map"]
                for trtmnt in elem_dict["map"]["treatment"]: # elem_dict["map"]["treatment"] is a list of dictionaries. trmnt is a dictionary
                    if len(trtmnt["Entity"]["strings"]) > 1:  print "Treatment: ", trtmnt["Entity"]["strings"]  # This print statement does not print anything
                    trtmnt["Entity"]["strings"] = trtmnt["Entity"]["strings"][0]    # Converting trtmnt["Entity"]["strings"] from a list that always contains a single element to a string.
                    if "uniprotSym" in trtmnt["Entity"]:   # trtmnt["Entity"] is a dictionary
                        #print "Treatment: ", trtmnt["Entity"]["strings"], "\t UniprotId: ", trtmnt["Entity"]["uniprotSym"]
                        yes_trtmnt_uniprot += 1
                    else:
                        no_trtmnt_uniprot += 1

        #for elem_dict in myfile_data:
        #    elem_dict["evidence_id"] = evid_dict[elem_dict["evidence"][0]]

        mydict = {}
        mydict["PMCID"] = pmcid
        mydict["Datums"] = myfile_data
        mydict.update(curr_pmcid_biblio)
        grand_list.append(mydict)

        #exit()
    print "No. of datums that have uniprotids for subject strings: ", yes_subj_uniprot, "\t Out of: ", (yes_subj_uniprot + no_subj_uniprot)
    print "No. of datums that have uniprotids for treatment strings: ", yes_trtmnt_uniprot, "\t Out of: ", (yes_trtmnt_uniprot + no_trtmnt_uniprot)
    json_fn = open("Datums_SRI_PreProcessed.json", 'w')
    json.dump(grand_list, json_fn, indent=4, ensure_ascii=False)


__main__();
