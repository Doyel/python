# -*- coding: utf-8 -*-

__author__ = 'Doyel'



import os
import json;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import argparse;
import cPickle as pickle;
import glob
import math
from pprint import pprint
from Utilities import read_config_file, check_for_unicode


# check_protsWithUnicodeChars::     python Experiment_Miscellenous.py Total_PubMedIds_John_NoOpenAccess.txt ./../protein_detect_output
# check_protsWithUnicodeChars::     python Experiment_Miscellenous.py PubMedIDS_with_Extras_OpenAccess.txt ./../protein_detect_output_OpenAccess
# check_sentids_NonOAarticles::     python Experiment_Miscellenous.py Total_PubMedIds_John_NoOpenAccess.txt ./../JSON_SENTENCE_ASSOCS NonOA_graph_dict.json
# check_sentids_OAarticles::        python Experiment_Miscellenous.py PubMedIDS_with_Extras_OpenAccess.txt ./../JSON_SENTENCE_ASSOCS_OpenAccess OA_graph_dict.json
# check_3proteins_passageDict::     python Experiment_Miscellenous.py ./../datafiles/datafiles_WithWeights_BaseLine1_Aug15

def __main__():
    # check_protsWithUnicodeChars()
    # check_sentids_NonOAarticles()
    # check_sentids_OAarticles()
    #check_3proteins_passageDict()

    prob = 0
    for i in xrange(1, 7):
        print i
        temp = (1 - math.pow((i/6), 4)) * ((math.pow(i, 2) - math.pow((i-1), 2))/36)
        prob += temp
    print prob



def check_3proteins_passageDict():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    proteinsToCheck = {"erk", "vegf-a", "raf-1"}    # proteinsToCheck is a set!!!!
    passage_dicts = ["train_passage_dict.json", "test_passage_dict.json", "openaccess_passage_dict.json"]

    parser = argparse.ArgumentParser(prog='Experiments',
                                     usage='Experiment_Miscellenous.py <passagedict_loc>',
                                     description='Use this script for checking / testing miscellenous stuff')
    parser.add_argument('passagedict_loc', help='Path where the passage dict files are located')
    args = parser.parse_args()

    myarticles = []
    for fname in passage_dicts:
        passagedictfile = os.path.join(args.passagedict_loc, fname)
        if os.path.isfile(passagedictfile):
            print "Found - " + fname + "! Everything is perfect in the world!"
            sys.stdout.flush()
            file_handle = open(passagedictfile, "rb")
            passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
            file_handle.close()
        else:
            print "Couldn't locate the dictionary - " + fname + " in the given directory! Please give the correct path."
            exit()

        for pmid in passage_dict:
            all_proteins = []
            for label in passage_dict[pmid]:
                for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                    all_proteins.extend(passage_dict[pmid][label][myclass].keys())
            all_proteins_uniq = list(set(all_proteins))
            #if pmid == "19050761":
            #    print all_proteins_uniq
            if proteinsToCheck.issubset(set(all_proteins_uniq)):
                myarticles.append(pmid)
    print "Articles that have the proteins: ", ",".join(list(proteinsToCheck)), " are: ", myarticles


def check_sentids_NonOAarticles():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Experiments', usage='Experiment_Miscellenous.py <filelist> <extracted_sent_dir> <networkX_dict>', description='Use this script for checking / testing miscellenous stuff')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('sent_dir_path', help='Path where the protein match files are located')
    parser.add_argument('networkX_dict', help='Name of the dict file that stores token details and networkX graph')
    args = parser.parse_args()

    if os.path.isfile(os.path.join(json_location, args.networkX_dict)):
        print "Found - " + args.networkX_dict + "! Everything is perfect in this world!"
        file_handle = open(os.path.join(json_location, args.networkX_dict), "rb")
        networkX_data = json.load(file_handle)  # networkX_data is a dictionary
        file_handle.close()
    else:
        print "The dictionary file - " + args.networkX_dict + " was not found in the JSON directory. Exiting!"
        exit()

    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f1:
        for pmid in f1:
            pmid = pmid.strip()
            if not pmid.strip(): continue
            fname_string1 = '%s' % (os.path.join(args.sent_dir_path, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                continue
            myfile = fname[0].strip()
            print "PMID: ", pmid

            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)  # myfile_data is a list of dictionaries
            file_handle.close()

            if "sentenceList" not in myfile_data: continue;
            if "@items" not in myfile_data["sentenceList"] or len(myfile_data["sentenceList"]["@items"]) == 0:
                continue

            sentidList_CoreNLP = sorted([int(k) for k in networkX_data[pmid].keys()])

            for itm in myfile_data["sentenceList"]["@items"]:  # itm is a dictionary. It represent one line of text in the file
                if itm["index"] != sentidList_CoreNLP[itm["index"]]:
                    print "The sentence ids are not matching! Please check!"
                    print itm["index"], "\t", sentidList_CoreNLP[itm["index"]]
                    exit()


def check_sentids_OAarticles():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Experiments', usage='Experiment_Miscellenous.py <filelist> <extracted_sent_dir> <networkX_dict>', description='Use this script for checking / testing miscellenous stuff')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('sent_dir_path', help='Path where the protein match files are located')
    parser.add_argument('networkX_dict', help='Name of the dict file that stores token details and networkX graph')
    args = parser.parse_args()

    if os.path.isfile(os.path.join(json_location, args.networkX_dict)):
        print "Found - " + args.networkX_dict + "! Everything is perfect in this world!"
        file_handle = open(os.path.join(json_location, args.networkX_dict), "rb")
        networkX_data = json.load(file_handle)  # networkX_data is a dictionary
        file_handle.close()
    else:
        print "The dictionary file - " + args.networkX_dict + " was not found in the JSON directory. Exiting!"
        exit()

    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f1:
        for pmid in f1:
            pmid = pmid.strip()
            if not pmid.strip(): continue
            fname_string1 = '%s*' % (os.path.join(args.sent_dir_path, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                continue
            myfile = fname[0].strip()
            print "PMID: ", pmid

            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)  # myfile_data is a list of dictionaries
            file_handle.close()

            sentidList_CoreNLP = sorted([int(k) for k in networkX_data[pmid].keys()])

            if len(myfile_data["paragraphList"]) == 0: continue;
            for para in myfile_data["paragraphList"]:  # para is a dictionary
                if "sentenceList" not in para or "@items" not in para["sentenceList"] or len(para["sentenceList"]["@items"]) == 0:
                    continue

                for itm in para["sentenceList"]["@items"]:  # itm represents a sentence in the article
                    if itm["absoluteId"] != sentidList_CoreNLP[itm["absoluteId"]]:
                        print "The sentence ids are not matching! Please check!"
                        print itm["absoluteId"], "\t", sentidList_CoreNLP[itm["absoluteId"]]
                        exit()


def check_protsWithUnicodeChars():
    global DatumKBdict; global id_to_proteins
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Experiments', usage='Experiment_Miscellenous.py <filelist> <prot_matchfiles_dir>', description='Use this script for checking / testing miscellenous stuff')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('match_dir_path', help='Path where the protein match files are located')
    args = parser.parse_args()

    DatumKBdict = pickle.load(open("pubmedid_extras.p", "rb"))
    id_to_proteins = pickle.load(open("UniProtid_to_proteins.p", "rb"))
    allowed_labels = ['reqs', 'RNAi', 'KO', 'omission']

    mydict = {}
    for label in allowed_labels:
        print "Processing label: ", label
        mydict[label] = {}
        with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f1:
            for pmid in f1:
                if not pmid.strip(): continue
                pmid = pmid.strip()

                fname_string = '%s_Paragraphs*' % (os.path.join(args.match_dir_path, pmid.strip()))
                fname = glob.glob(fname_string)
                if len(fname) == 0:
                    #print "The Paragraphs file for PubMedId: ", pmid.strip(), "was not found!"  # This statement will never get printed
                    continue
                myparafile = fname[0].strip()

                file_handle = open(myparafile, "rb")
                myparafile_data = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
                file_handle.close()

                for elem in myparafile_data:
                    for prot in elem["matches"].keys(): # prot is a matched protein which is in lower case
                        found = False   # 'found' may be false in 2 cases - When the protein entity is not a key in the dict - DatumKBdict[pmid] or when the protein entity is present
                        # as a key but there is no extra datum in that article which contains the protein entity and the given label as a treatment test in the same datum.
                        if pmid in DatumKBdict:
                            if prot in DatumKBdict[pmid].keys():
                                if label == "reqs":
                                    if len(DatumKBdict[pmid][prot]['reqstest']) > 0: found = True    # We check the length because of lines 114 and 116 in StratifyArticles.py
                                elif label == "dnreqs":
                                    if len(DatumKBdict[pmid][prot]['dnreqstest']) > 0: found = True
                                else:
                                    for test in DatumKBdict[pmid][prot].keys():      # test can only be 'reqstest' or 'dnreqstest'
                                        if label in DatumKBdict[pmid][prot][test]:    found = True
                                foundprot = prot
                            else:
                                found,foundprot = find_protsyn(elem["matches"][prot]["UniProtId"], label, pmid.strip())

                        if found and check_for_unicode(prot):
                            if pmid in mydict[label]:
                                if prot not in mydict[label][pmid]:
                                    mydict[label][pmid].append(prot)
                            else:
                                mydict[label][pmid] = [prot]
    pprint(mydict)


def find_protsyn(uniprotidstr, label, pmid):
    found = False
    uniprotidlist = uniprotidstr.split(",")
    for upid in uniprotidlist:
        for syn in id_to_proteins[upid.strip()]:
            syn = syn.lower()
            if syn in DatumKBdict[pmid].keys():
                if label == "reqs":
                    if len(DatumKBdict[pmid][syn]['reqstest']) > 0: found = True     # We check the length because of lines 114 and 116 in StratifyArticles.py
                elif label == "dnreqs":
                    if len(DatumKBdict[pmid][syn]['dnreqstest']) > 0: found = True
                else:
                    for test in DatumKBdict[pmid][syn].keys():
                        if label in DatumKBdict[pmid][syn][test]:    found = True
                if found:
                    return found, syn
    return False, "NIL"


__main__()



