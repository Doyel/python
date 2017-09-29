# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
from pprint import pprint;
import os;
import json;
import cPickle as pickle;

#python Extras_OpenAccess.py ./datums_03_2016.json ./PMCIDS_with_Extras_OpenAccess.txt

def main_body():
    parser = argparse.ArgumentParser(prog='Extras OpenAccess', usage='Extras_OpenAccess.py <datums_JSON_file> <pmcidfilelist>', description='Script to get datums having Extras fields of type reqs/do not req which are in OpenAccess')
    parser.add_argument('datums_file', help='newest version of the datums file given by Mark')    
    parser.add_argument('PMCIDfilelist', help='File listing the PMCIDs that need to be downloaded')
    args = parser.parse_args()

    summary_dict = pickle.load( open( "summary_dict_Extras_PMCID.p", "rb" ) )

    pmidlist = []
    with open(args.PMCIDfilelist, 'rt') as f1:
        for pmcid in f1:
            pmidlist.append(summary_dict[pmcid.strip()]["PMID"])

    file_handle = open(args.datums_file, "rb")
    myfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
    file_handle.close();

    myjson = []; tot_datums = 0;
    for elem in myfile_data:    # elem represents a datum and is a dictionary  
        if elem["source"]["pmid"].strip() not in pmidlist:
            continue      
        myjson.append(elem); tot_datums += 1;      

    json_fn = open("datums_OpenAccess.json", 'w')
    json.dump(myjson, json_fn, indent=4, ensure_ascii=False)
    print "The no. of total datums in JSON file: ", tot_datums

    pmid_fn1 = open("PubMedIDS_with_Extras_OpenAccess.txt", 'w');
    for pm in pmidlist:    
        pmid_fn1.write(pm+'\n')
    pmid_fn1.close()

    #---------------------------------------------------------------------------------------------------#

    my5json = []; my5pmidlist = pmidlist[:10]
    my5pmidlist.sort(); print my5pmidlist, "\n";
    for elem in myfile_data:    # elem represents a datum and is a dictionary  
        if elem["source"]["pmid"].strip() not in my5pmidlist:
            continue    
        if 'extras' in elem:
            if len(elem["extras"]) == 0:
                continue
            thisdatum_reqs = False; thisdatum_dnreqs = False
            for ext in elem["extras"]:  # ext is a list element which is actually a dictionary
                if ext["type"] == "reqs":                    
                    thisdatum_reqs = True;
                if ext["type"] == "does not req":                   
                    thisdatum_dnreqs = True;

            if thisdatum_reqs or thisdatum_dnreqs:
                my5json.append(elem);   
    
    newlist = sorted(my5json, key=lambda k: k['source']['pmid'])          
    json_fn5 = open("datums_my5PMCIDs.json", 'w')
    json.dump(newlist, json_fn5, indent=4, ensure_ascii=False)

    pmid_pmcid = pickle.load( open( "pmid_pmcid_Extras.p", "rb" ) )
    for pm in my5pmidlist:
        print pm, "  -->  ", pmid_pmcid[pm]
    
    print "The no. of total reqs/not reqs datums in the 5 articles: ", len(my5json), "   ", len(newlist)


if __name__ == "__main__":
    main_body()
