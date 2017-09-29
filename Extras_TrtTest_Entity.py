# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
from pprint import pprint;
import os;
import json;

#python Extras_TrtTest_Entity.py ./datums_03_2016.json

def main_body():
    parser = argparse.ArgumentParser(prog='Extras with Treatment Test and Entities', usage='Extras_TrtTest_Entity.py <datums_JSON_file>', description='Script to get details about the Extras fields')
    parser.add_argument('datums_file', help='newest version of the datums file given by Mark')    
    args = parser.parse_args()

    file_handle = open(args.datums_file, "rb")
    myfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
    file_handle.close();
    
    elem_multi = []        
    for elem in myfile_data:    # elem represents a datum and is a dictionary
        if 'extras' in elem:
            if len(elem["extras"]) == 0:
                continue
            uniq_entities = []; thisdatum_reqs = False; thisdatum_dnreqs = False
            for ext in elem["extras"]:  # ext is a list element which is actually a dictionary
                if ext["type"] == "reqs":
                    if "entities" in ext:  
                        uniq_entities.extend(ext["entities"])                         
                    thisdatum_reqs = True;
                if ext["type"] == "does not req": 
                    if "entities" in ext:  
                        uniq_entities.extend(ext["entities"])                                   
                    thisdatum_dnreqs = True;

            if thisdatum_reqs or thisdatum_dnreqs:                
                uniq_entities = list(set(uniq_entities)); trtmnttest_list = []
                for ent in uniq_entities:
                    for ext in elem["extras"]:                        
                        if "entities" in ext and (ext["type"] == "does not req" or ext["type"] == "reqs") and ent in ext["entities"]:
                            if "test" in ext["treatment"]:  
                                trtmnttest_list.append(ext["treatment"]["test"])
                    trtmnttest_list = list(set(trtmnttest_list))
                    if len(trtmnttest_list) > 1:
                        elem_multi.append(elem)
                    break          

    print "The no of datums in Datum KB having multiple treatment tests for a single protein entity: ", len(elem_multi)
    json_fn = open("Datums_Extras_Multi_Trtments.json", 'w')
    json.dump(elem_multi, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()

if __name__ == "__main__":
    main_body()
