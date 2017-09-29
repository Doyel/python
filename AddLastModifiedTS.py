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


# python AddLastModifiedTS.py ./user_feedback_JSON/Apr4_2017/Datums_MongoDB_OpenAccess_WithFeedback_Apr12.json ./user_feedback_JSON/Apr4_2017/articles_BIG_MECHANISM_Beautified.json


def main_body():
    global pruned_words_list
    parser = argparse.ArgumentParser(prog='AddLastModifiedTS', usage='AddLastModifiedTS.py <MongoDBfile> <articlesjson>', description='Script to beautify the specified mongoimported JSON file')
    parser.add_argument('MongoDBfilepath', help='Filepath to the mongodb database file that contains the user feedback')
    parser.add_argument('articlesjson', help='Filepath to the earlier articles.json file that contains the last modified TS')
    args = parser.parse_args()

    # ["14517278", "19783983", "18411307", "20890305", "18504304", "19234442", "22245064", "20141835", "21573184", "20547488", "22433566", "20026654"]
    ignore_PMIDS = ["20026654"]
    
    MongoDBfile = ntpath.basename(args.MongoDBfilepath)
    drname = ntpath.dirname(args.MongoDBfilepath)
    if os.path.isfile(args.MongoDBfilepath):
        print "Found the dictionary - " + MongoDBfile + " Everything is perfect in the world!"
        file_handle = open(args.MongoDBfilepath, "rb")
        mongodb_dict = json.load(file_handle)  
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + MongoDBfile + "! Please place it in - ", args.MongoDBfilepath
        exit()
        
    articlesfile = ntpath.basename(args.articlesjson)
    if os.path.isfile(args.articlesjson):
        print "Found the dictionary - " + articlesfile + " Everything is perfect in the world!"
        file_handle = open(args.articlesjson, "rb")
        articles_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + articlesfile + "! Please place it in - ", args.articlesjson
        exit()
        
    lastupdated_dict = {}
    for art in articles_dict:
        if art["PMID"] in ignore_PMIDS: 
            print "Ignoring the PMID: ", art["PMID"]
            continue
        if "lastUpdated" in art:
            #print art["lastUpdated"]
            lastupdated_dict[art["PMID"]] = dict(art["lastUpdated"])
    
    for art in mongodb_dict:
        if art["PMID"] in lastupdated_dict:
            art["lastUpdated"] = lastupdated_dict[art["PMID"]]  

    mongodb_datafile = os.path.join(drname, MongoDBfile.split(".")[0]+"_TS."+MongoDBfile.split(".")[1])
    json_fn = open(mongodb_datafile, 'w')
    json.dump(mongodb_dict, json_fn, indent=4, ensure_ascii=False)
           
            
if __name__ == "__main__":
    main_body()
