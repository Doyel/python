# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import re;
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import math;
from pprint import pprint;
import cPickle as pickle;
from random import shuffle
from Utilities import read_config_file

# python StratifyArticles.py ./PubMedIds_with_Extras.txt ./datums_03_2016.json ./Total_PubMedIds_John.txt ./protein_detect_output      

def main_body():
    global pmid_num_reqs_dnreqs; global openaccess
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='StratifyArticles', usage='StratifyArticles.py <PubMedfilelist> <datumKB_file> <John_PubMedfilelist> <proteins_detected_dir>', description='Script to stratify articles according extras type')    
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids that have relevant extras datums in the Datum KB')
    parser.add_argument('datums_file', help='newest version of the datums file given by Mark')
    parser.add_argument('John_PubMedfilelist', help='File listing the PubMed ids of all the files sent by John')
    parser.add_argument('protein_outputdir', help='Path where files containing the detected proteins are located')    
    args = parser.parse_args()

    openaccess = []
    with open(os.path.join(txt_location, "PubMedIDS_with_Extras_OpenAccess.txt"), 'rt') as f:
        for pmid in f:
            if not pmid.strip(): continue   
            
            fname_string1 = '%s*' % (os.path.join(parent_location, "JSON_SENTENCE_ASSOCS", pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                #print "The segmented sentences file for PubMedId: ", pmid.strip(), " was not sent by John"
                continue              
            openaccess.append(pmid.strip()) 
    print "No. of articles that are in OpenAccess and also in John's list: ", len(openaccess)
        
    mypmids = []    # PubmedIds from the Datum KB that have extras of type 'reqs' or 'do not reqs' and are present in John's list and have extracted sentences and proteins were detected in the sentences
    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f1:
        print "Processing file: ", args.PubMedfilelist
        for pmid in f1:
            if not pmid.strip(): continue          
            #if pmid.strip() in openaccess: continue

            fname_string1 = '%s*' % (os.path.join(parent_location, "JSON_SENTENCE_ASSOCS", pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                #print "The segmented sentences file for PubMedId: ", pmid.strip(), " was not sent by John"
                continue                   

            fname_string1 = '%s*' % (os.path.join(args.protein_outputdir, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                #print "No proteins were detected for PubMedId: ", pmid.strip()
                continue             

            mypmids.append(pmid.strip())

    mypmids_noextras = []   # PubmedIds sent by John that do not have extras of our interest, have extracted sentences and proteins were detected in the sentences
    with open(os.path.join(txt_location, args.John_PubMedfilelist), 'rt') as f2:
        print "Processing file: ", args.John_PubMedfilelist
        for pmid in f2:
            if not pmid.strip(): continue 
            if pmid.strip() in mypmids: continue
            #if pmid.strip() in openaccess: continue
            
            fname_string1 = '%s*' % (os.path.join(args.protein_outputdir, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                #print "No proteins were detected for PubMedId: ", pmid.strip()
                continue
    
            mypmids_noextras.append(pmid.strip())

    file_handle = open(os.path.join(json_location, args.datums_file), "rb")
    datumfile = json.load(file_handle)    # datumfile is a list of dictionaries
    file_handle.close();

    mydict = {}; pmid_num_reqs_dnreqs = {}
    for elem in datumfile:    # elem represents a datum and is a dictionary 
        curr_pmid = elem["source"]["pmid"].strip()       
        if 'extras' in elem and len(elem["extras"]) > 0:           
            for ext in elem["extras"]:  # elem["extras"] is a list and ext is a dictionary
                if ext["type"] == "reqs" or ext["type"] == "does not req":
                    
                    if curr_pmid not in pmid_num_reqs_dnreqs:
                        pmid_num_reqs_dnreqs[curr_pmid] = {}
                    if ext["type"] == "reqs":
                        if "reqs" in pmid_num_reqs_dnreqs[curr_pmid]:   pmid_num_reqs_dnreqs[curr_pmid]["reqs"]+=1                            
                        else:   pmid_num_reqs_dnreqs[curr_pmid]["reqs"] = 1                            
                    else: 
                        if "dnreqs" in pmid_num_reqs_dnreqs[curr_pmid]: pmid_num_reqs_dnreqs[curr_pmid]["dnreqs"]+=1                            
                        else:   pmid_num_reqs_dnreqs[curr_pmid]["dnreqs"] = 1
                            
                    if curr_pmid not in mydict:
                        mydict[curr_pmid] = {}  # There will be an entry for a PMID in mydict only if that article contains atleast one datum containing atleast one  extra element datum which is of type "reqs" or "does not req"
                    if "entities" in ext:
                        for ent in ext["entities"]:     # ext["entities"] is a list and may have more than one element in it
                            ent = ent.lower()
                            if isinstance(ext["treatment"]["test"], list):      # Assumption: If 'entities' exist, then 'treatment' will also exist in the same 'extra' element
                                print "test is a list! ", curr_pmid
                                exit()
                            if ent in mydict[curr_pmid]:                                                                    
                                if ext["type"] == "reqs":
                                    if ext["treatment"]["test"] not in mydict[curr_pmid][ent]["reqstest"]:   
                                        mydict[curr_pmid][ent]["reqstest"].append(ext["treatment"]["test"])
                                else:
                                    if ext["treatment"]["test"] not in mydict[curr_pmid][ent]["dnreqstest"]:   
                                        mydict[curr_pmid][ent]["dnreqstest"].append(ext["treatment"]["test"])
                            else:
                                mydict[curr_pmid][ent] = {}                                
                                if ext["type"] == "reqs":  
                                    mydict[curr_pmid][ent]["reqstest"] = [ext["treatment"]["test"]]    
                                    mydict[curr_pmid][ent]["dnreqstest"] = []                                                                     
                                else:   
                                    mydict[curr_pmid][ent]["reqstest"] = []
                                    mydict[curr_pmid][ent]["dnreqstest"] = [ext["treatment"]["test"]]                                                                        
                    else:   # 'entities' does not exist in ext
                        if "treatment" in ext:  # 'entities' does not exist BUT 'treatment' exists in ext -----> I doubt that control will ever come here!
                            print "AHEM! Control has come here!!!"  # [UPDATE]: I was right! Control never came here.
                            if isinstance(ext["treatment"]["test"], list):  
                                print "test is a list! ", curr_pmid
                                exit()
                            if "absent" in mydict[curr_pmid]:                                    
                                if ext["type"] == "reqs":   mydict[curr_pmid]["absent"]["reqstest"].append(ext["treatment"]["test"])
                                else:   mydict[curr_pmid]["absent"]["dnreqstest"].append(ext["treatment"]["test"])
                            else:
                                mydict[curr_pmid]["absent"] = {}
                                if ext["type"] == "reqs":  mydict[curr_pmid]["absent"]["reqstest"] = [ext["treatment"]["test"]]                                                                         
                                else:   mydict[curr_pmid]["absent"]["dnreqstest"] = [ext["treatment"]["test"]]                        
    for key in mydict.keys():
        if not mydict[key]:   # A datum that has one or more 'extra' datums of type 'reqs' or 'do not reqs' may not have any 'entities' entry in any of its 'extra' datums
            del mydict[key]

    protein_names = pickle.load( open( "list_UniProtId_name_synonym.p", "rb" ) )
    id_to_proteins = {}
    for itm in protein_names:
        if itm['id'].strip() not in id_to_proteins:
            id_to_proteins[itm['id'].strip()] = [itm['name'].strip()]
        else:
            id_to_proteins[itm['id'].strip()].append(itm['name'].strip())
        if 'synonym' in itm: id_to_proteins[itm['id'].strip()].extend(itm['synonym'])
    pickle.dump( id_to_proteins, open( "UniProtid_to_proteins.p", "wb" ) )

    print "Created the updated UniProtid_to_proteins.p"
    print "Exiting out before generating stratified train and test lists!"
    print "Please comment out the exit statement at line 154 to continue further"
    exit()

    pickle.dump(mydict, open("pubmedid_extras.p", "wb"))
    pickle.dump(pmid_num_reqs_dnreqs, open("pmid_num_reqs_dnreqs_extras.p", "wb"))

    print "The total no. of pmids that DO NOT have interesting extra elements in the DatumKB, are in John's set, have sentences extracted from them and proteins were detected in them: ", len(mypmids_noextras)
    training_set_noextras_size = math.ceil(0.8 * float(len(mypmids_noextras)))
    training_set_noextras_pmids, test_set_noextras_pmids = randomize_list(mypmids_noextras, training_set_noextras_size)    
    
    print "The total no. of pmids that have interesting extra elements in the DatumKB, are also in John's set, have sentences extracted from them and proteins were detected in them: ", len(mypmids)    
    print "Proportion of reqs and dnreqs extras over all PubMed articles that have interesting extra datums, are also in John's set, have sentences extracted from them and proteins were detected in them: "
    get_proportion(mypmids); print "\n"    

    mark_four_bins(mypmids, training_set_noextras_pmids, test_set_noextras_pmids)    
    exit() # Dont want to execute below statements!

    training_set_extras_size = math.ceil(0.8 * float(len(mypmids)))
    training_set_extras_pmids, test_set_extras_pmids = create_sets(mypmids, training_set_extras_size)    

    training_set = []; test_set = []
    training_set.extend(training_set_extras_pmids); training_set.extend(training_set_noextras_pmids)
    test_set.extend(test_set_extras_pmids); test_set.extend(test_set_noextras_pmids)

    pmid_fn1 = open(os.path.join(txt_location, "Training_Set_PubMedIds.txt"), 'w')
    for pm in training_set:
        pmid_fn1.write(pm+'\n')
    pmid_fn1.close()

    pmid_fn1 = open(os.path.join(txt_location, "Test_Set_PubMedIds.txt"), 'w')
    for pm in test_set:
        pmid_fn1.write(pm+'\n')
    pmid_fn1.close()
                    

def get_proportion(mypmids):
    num_reqs = 0.0; num_dnreqs = 0.0
    for pmid in mypmids:        
        if "reqs" in pmid_num_reqs_dnreqs[pmid]:    num_reqs += pmid_num_reqs_dnreqs[pmid]["reqs"]
        if "dnreqs" in pmid_num_reqs_dnreqs[pmid]:    num_dnreqs += pmid_num_reqs_dnreqs[pmid]["dnreqs"]
    print "The proportion of 'reqs' extras: ", float(100 * (num_reqs / (num_reqs + num_dnreqs)))
    print "The proportion of 'dnreqs' extras: ", float(100 * (num_dnreqs / (num_reqs + num_dnreqs)))    


def create_sets(mypmids, training_set_size):
    print "Inside Create Sets: "
    myresponse = "N"
    while myresponse == "N" or myresponse == "n":
        training_set_pmids, test_set_pmids = randomize_list(mypmids, training_set_size)
        print "Training set metric: "
        get_proportion(training_set_pmids)
        print "Test set metric: "
        get_proportion(test_set_pmids)
        myresponse = raw_input("Are you happy with the proportions?  [Y/N]:   ")
        myresponse = myresponse.strip()
        print "\n"
    return training_set_pmids, test_set_pmids    
 

def randomize_list(mypmids, training_set_size):
    shuffle(mypmids)
    training_set_pmids = mypmids[:int(training_set_size)]
    test_set_pmids = mypmids[int(training_set_size):]
    return training_set_pmids, test_set_pmids


def mark_four_bins(mypmids, training_set_noextras_pmids, test_set_noextras_pmids):
    print "Inside Mark\'s four bins method: "
    only_reqs_pmids = []; only_dnreqs_pmids = []; reqs_dnreqs_pmids = []
    for pmid in mypmids:
        if "reqs" not in pmid_num_reqs_dnreqs[pmid]: 
            only_dnreqs_pmids.append(pmid); continue
        if "dnreqs" not in pmid_num_reqs_dnreqs[pmid]: 
            only_reqs_pmids.append(pmid); continue
        reqs_dnreqs_pmids.append(pmid)

    shuffle(only_reqs_pmids); shuffle(only_dnreqs_pmids); shuffle(reqs_dnreqs_pmids)
    training_set_onlyreqs_size = math.ceil(0.8 * float(len(only_reqs_pmids)))
    training_set_onlydnreqs_size = math.ceil(0.8 * float(len(only_dnreqs_pmids)))
    training_set_reqsdnreqs_size = math.ceil(0.8 * float(len(reqs_dnreqs_pmids)))

    mark_training_set = []; mark_test_set = []
    mark_training_set.extend(only_reqs_pmids[:int(training_set_onlyreqs_size)])
    mark_test_set.extend(only_reqs_pmids[int(training_set_onlyreqs_size):])
    mark_training_set.extend(only_dnreqs_pmids[:int(training_set_onlydnreqs_size)])
    mark_test_set.extend(only_dnreqs_pmids[int(training_set_onlydnreqs_size):])
    mark_training_set.extend(reqs_dnreqs_pmids[:int(training_set_reqsdnreqs_size)])
    mark_test_set.extend(reqs_dnreqs_pmids[int(training_set_reqsdnreqs_size):])

    #mark_training_set.extend(openaccess)
    
    print "Mark\'s four bins - Training set metric: "
    get_proportion(mark_training_set)
    print "Mark\'s four bins - Test set metric: "
    get_proportion(mark_test_set)
    
    mark_training_set.extend(training_set_noextras_pmids)
    mark_test_set.extend(test_set_noextras_pmids)

    pmid_fn1 = open(os.path.join(txt_location, "Mark_Training_Set_PubMedIds.txt"), 'w')
    for pm in mark_training_set:
        pmid_fn1.write(pm+'\n')
    pmid_fn1.close()

    pmid_fn1 = open(os.path.join(txt_location, "Mark_Test_Set_PubMedIds.txt"), 'w')
    for pm in mark_test_set:
        pmid_fn1.write(pm+'\n')
    pmid_fn1.close()    
    print "\n\n"

    ctr=0; notincluded = []
    for pmid in openaccess:
        if pmid not in mark_training_set:
            ctr+=1; notincluded.append(pmid)
    print ctr, " articles not present in training set list"
    print notincluded


if __name__ == "__main__":
    main_body()
