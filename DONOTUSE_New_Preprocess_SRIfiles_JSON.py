# -*- coding: utf-8 -*-
__author__ = 'Doyel'


import json;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import re;
import copy;
import argparse;
import os;
import glob;
import cPickle as pickle;
from pprint import pprint;
import difflib;


#python Preprocess_SRIfiles_JSON.py ../../extracted-proto-datums/ SRI_json_Filelist.txt /ua/ml-group/big-mechanism-project/preprocessed_SRIextracted_datums

def __main__():
    parser = argparse.ArgumentParser(prog='Preprocess SRI files', usage='Preprocess_SRIfiles_JSON.py <dirpath> <filelist> <preprocessed_dir>', description='Script that pre-processes JSON files sent by SRI and adds new fields in each datum')
    parser.add_argument('dir_path', help='Path where the JSON files obtained from SRI are located')
    parser.add_argument('JSONfilelist', help='File listing the names of JSON files to be processed')
    parser.add_argument('preprocessed_dir_path', help='Path where the pre-processed JSON files are written')
    args = parser.parse_args()

    if os.path.isfile("evidenceid_evidsent_dict.p") and os.path.isfile("datumid_evidenceid_dict.p") and os.path.isfile("xmlpath_evidenceid_dict.p"):
        print "Found evidenceid_evidsent_dict.p, datumid_evidenceid_dict.p and xmlpath_evidenceid_dict.p in current directory. Using them!"
        evid_evidsent = pickle.load( open( "evidenceid_evidsent_dict.p", "rb" ) )
        datum_evidenceid = pickle.load( open( "datumid_evidenceid_dict.p", "rb" ) )
        xmlpath_evidenceid = pickle.load( open( "xmlpath_evidenceid_dict.p", "rb" ) )        
    else:
        print "Could not find evidenceid_evidsent_dict.p, datumid_evidenceid_dict.p and xmlpath_evidenceid_dict.p in current directory. Building them from scratch!"
        evid_evidsent = {}
        datum_evidenceid = {}
        xmlpath_evidenceid = {}        

    with open(args.JSONfilelist, 'rt') as f1:   ##################### Please make sure that JSONfilelist is a list of PMCIDs that have not been processed previously #########################
        for pmcid in f1:
            if not pmcid.strip(): continue
            fname_string = '%s' % (os.path.join(args.dir_path, pmcid.strip()))
            fname = glob.glob(fname_string)
            if len(fname) == 0:
                print "The specified JSON file was not found: ", pmcid.strip()
                continue
            elif len(fname) > 1:
                print "Two files having the same name exist: ", pmcid.strip()
                exit()     
            myfile = fname[0].strip()
            #print myfile
            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
            file_handle.close();

            if len(myfile_data) == 0:   # If the file has no datums in it then continue to the next json file
                print "Skipping file: ", myfile; continue

            pmid = os.path.splitext(os.path.basename(myfile))[0]
            basefilename = os.path.basename(myfile)
            print "Processing PMID: ", pmid

            if pmid not in xmlpath_evidenceid:      
                xmlpath_evidenceid[pmid] = {}
            else:
                pmcid_processed = True
            if pmid not in evid_evidsent:
                evid_evidsent[pmid] = {}
            else:
                pmcid_processed = True
            if pmid not in datum_evidenceid:
                datum_evidenceid[pmid] = {}
            else:
                pmcid_processed = True
            if pmcid_processed:
                print "***************************Re-processing PMCID: ", pmid, ", had already been processed earlier***************************"

            xpath_evid = {};
            for elem_dict in myfile_data:
                if len(elem_dict["evidence"]) > 1:  print "evidence has more than 1 element in its list"

                xmlpaths = {}
                if "subject" in elem_dict["map"]:
                    new_check_datumfield(elem_dict["map"]["subject"], xmlpaths, "subject", pmid, elem_dict["evidence"][0])    # xmlpaths: pass by reference
                                    
                if "assay" in elem_dict["map"]:
                    new_check_datumfield(elem_dict["map"]["assay"], xmlpaths, "assay", pmid, elem_dict["evidence"][0])    # xmlpaths: pass by reference
                    
                if "change" in elem_dict["map"]:
                    new_check_datumfield(elem_dict["map"]["change"], xmlpaths, "change", pmid, elem_dict["evidence"][0])    # xmlpaths: pass by reference
                    
                if "treatment" in elem_dict["map"]:
                    new_check_datumfield(elem_dict["map"]["treatment"], xmlpaths, "treatment", pmid, elem_dict["evidence"][0])    # xmlpaths: pass by reference
                    
                for field in xmlpaths:
                    if xmlpaths[field]["paragraphXpath"] not in xpath_evid:
                        xpath_evid[xmlpaths[field]["paragraphXpath"]] = [xmlpaths[field]["evidence"]]
                    else:
                        if xmlpaths[field]["evidence"] not in xpath_evid[xmlpaths[field]["paragraphXpath"]]:
                            xpath_evid[xmlpaths[field]["paragraphXpath"]].append(xmlpaths[field]["evidence"])

            pprint (xpath_evid)
            evid_dict = {}; ctr_evid = 1; temp_xmlpath_evidenceid = {}
            for xpath in xpath_evid:
                tlist = list(xpath_evid[xpath]); currlen = len(tlist); prevlen=0;                
                while currlen != prevlen:
                    mylist = []; has_overlap = False;
                    tlist.sort(lambda x,y: cmp(len(x), len(y)))                
                    for i in xrange(len(tlist)):
                        if has_overlap: 
                            has_overlap = False; continue;
                        is_substr = False;
                        for j in xrange(i+1, len(tlist)):
                            if tlist[i] in tlist[j]:    #Check whether tlist[i] is a substring of tlist[j]
                                is_substr = True; break
                            combinedstr, has_overlap = check_overlap(tlist[i], tlist[j])    # tlist[j] should end with tlist[i] or tlist[i] should start with tlist[j]
                            if has_overlap:
                                mylist.append(combinedstr);
                                break
                        if not is_substr and not has_overlap:   mylist.append(tlist[i])
                    prevlen = len(tlist); currlen = len(mylist); tlist = list(mylist)       # Assign the contents of mylist to tlist
                print tlist; temp_xmlpath_evidenceid[xpath] = [];                 
                sys.exit(); 
                for evid in tlist:
                    if evid in evid_dict:
                        print "The evidence sentence: ", evid, " is present in multiple XPaths within the article: ", pmid
                        sys.exit()
                    evid_dict[evid] = pmid + "_" + str(ctr_evid); ctr_evid+=1;                    
                    temp_xmlpath_evidenceid[xpath].append(evid_dict[evid])
                xpath_evid[xpath] = list(tlist)

            pprint (xpath_evid)
            xmlpath_evidenceid[pmid] = copy.deepcopy(temp_xmlpath_evidenceid) 
            evid_evidsent[pmid] = dict((v,k) for k,v in evid_dict.iteritems())                               
            sys.exit()   
              
            ctr_datum = 1; temp_datum_evidenceid = {}            
            for elem_dict in myfile_data:                                                               
                ####### Assumption: There is ONLY one evidence sentence (and hence ONE evidence_id) for each field of a datum (i.e. the datum represented by elem_dict)
                elem_dict["datum_id"] = pmid + "_" + str(ctr_datum); ctr_datum+=1;

                temp_datum_evidenceid[elem_dict["datum_id"]] = {}
                if "subject" in elem_dict["map"]:
                    new_check_datumfield(elem_dict["map"]["subject"], xmlpaths, "subject", pmid, elem_dict["evidence"][0])    # xmlpaths: pass by reference
                    for sent in xpath_evid[xmlpaths[field]["paragraphXpath"]]:
                        if xmlpaths[field]["evidence"] in sent:                            
                            temp_datum_evidenceid[elem_dict["datum_id"]]["subject"] = [evid_dict[sent]]
                            break
                    if 'subject' not in temp_datum_evidenceid[elem_dict["datum_id"]]:   print "Something went wrong! Subject key not present"
                else:   temp_datum_evidenceid[elem_dict["datum_id"]]["subject"] = []

                if "assay" in elem_dict["map"]:
                    new_check_datumfield(elem_dict["map"]["assay"], xmlpaths, "assay", pmid, elem_dict["evidence"][0])    # xmlpaths: pass by reference
                    for sent in xpath_evid[xmlpaths[field]["paragraphXpath"]]:
                        if xmlpaths[field]["evidence"] in sent:                            
                            temp_datum_evidenceid[elem_dict["datum_id"]]["assay"] = [evid_dict[sent]]
                            break
                    if 'assay' not in temp_datum_evidenceid[elem_dict["datum_id"]]:   print "Something went wrong! Assay key not present"
                else:   temp_datum_evidenceid[elem_dict["datum_id"]]["assay"] = []

                if "change" in elem_dict["map"]:
                    new_check_datumfield(elem_dict["map"]["change"], xmlpaths, "change", pmid, elem_dict["evidence"][0])    # xmlpaths: pass by reference
                    for sent in xpath_evid[xmlpaths[field]["paragraphXpath"]]:
                        if xmlpaths[field]["evidence"] in sent:                            
                            temp_datum_evidenceid[elem_dict["datum_id"]]["change"] = [evid_dict[sent]]
                            break
                    if 'change' not in temp_datum_evidenceid[elem_dict["datum_id"]]:   print "Something went wrong! Change key not present"
                else:   temp_datum_evidenceid[elem_dict["datum_id"]]["change"] = []

                if "treatment" in elem_dict["map"]:
                    new_check_datumfield(elem_dict["map"]["treatment"], xmlpaths, "treatment", pmid, elem_dict["evidence"][0])    # xmlpaths: pass by reference
                    for sent in xpath_evid[xmlpaths[field]["paragraphXpath"]]:
                        if xmlpaths[field]["evidence"] in sent:                            
                            temp_datum_evidenceid[elem_dict["datum_id"]]["treatment"] = [evid_dict[sent]]
                            break
                    if 'treatment' not in temp_datum_evidenceid[elem_dict["datum_id"]]:   print "Something went wrong! Treatment key not present"
                else:   temp_datum_evidenceid[elem_dict["datum_id"]]["treatment"] = []
                                                                                                                     
            outfile = os.path.join(args.preprocessed_dir_path, basefilename)
            json_fn = open(outfile, 'w')
            json.dump(myfile_data, json_fn, indent=4, ensure_ascii=False)
            datum_evidenceid[pmid] = copy.deepcopy(temp_datum_evidenceid)
            

    pickle.dump( evid_evidsent, open( "evidenceid_evidsent_dict.p", "wb" ) ) 
    pickle.dump( datum_evidenceid, open( "datumid_evidenceid_dict.p", "wb" ) )
    pickle.dump( xmlpath_evidenceid, open( "xmlpath_evidenceid_dict.p", "wb" ) )


def new_check_datumfield(mylist, xmlpaths, field, pmid, evidence):
    if len(mylist) > 1 and field <> "change":   # The datum field - change, always has multiple elements in elem_dict["map"]["change"]
        print field, " has more than 1 element for pmid: ", pmid
        pprint(mylist)
    mydict = mylist[0]; xmlpaths[field] = {}; 
    if "sectionXpath" in mydict["Provenance"]:  # sectionXpath may not be present in mydict["Provenance"] in some datums!!!
        if mydict["Provenance"]["paragraphXpath"].startswith(mydict["Provenance"]["sectionXpath"], 0):            
            #xmlpaths[field]["sectionXpath"] = mydict["Provenance"]["sectionXpath"]
            xmlpaths[field]["paragraphXpath"] = mydict["Provenance"]["paragraphXpath"]
            xmlpaths[field]["evidence"] = evidence      # Should be mydict["Provenance"]["evidence"] 
        else:
            print "paragraphXpath is different from sectionXpath"
            print mydict["Provenance"]["sectionXpath"]
            print mydict["Provenance"]["paragraphXpath"]
            exit()
    else:
        # Test for the key 'paragraphXpath' in mydict["Provenance"]. If it doesnt exist, then skip this datum altogether (return bad_datum = 1)
        #xmlpaths[field]["sectionXpath"] = ""
        xmlpaths[field]["paragraphXpath"] = mydict["Provenance"]["paragraphXpath"]
        xmlpaths[field]["evidence"] = evidence          # Should be mydict["Provenance"]["evidence"] 
            
            
def check_overlap(strI, strJ):
    strI = strI.strip(); strJ = strJ.strip()
    s = difflib.SequenceMatcher(None,strI, strJ)
    pos_a, pos_b, size = s.find_longest_match(0, len(strI), 0, len(strJ))
    if pos_a == 0 and pos_b == 0 and size > 0:
        print "The first str: ", strI, " is a substr of the second string: ", strJ
        print "Control should never come here!"
        return "", False
    if pos_a == 0 and pos_b == 0 and size == 0:
        print "No match between first str: ", strI, " and second string: ", strJ
        return "", False 
    if (pos_a == 0 and pos_b > 0) or (pos_a > 0 and pos_b == 0):
        if size >= 10:
            if pos_a == 0:
                return strJ + strI[pos_a+size:], True
            else:
                return strI + strJ[pos_b+size:], True
        else:
            print "Too less match between first str: ", strI, " and second string: ", strJ
            print "Treating as no overlap!"
            return "", False 
    print "Some phrase is matching in the middle of the strings!"
    print "The first str: ", strI, " The second string: ", strJ
    return "", False


__main__();
                
