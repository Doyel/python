import argparse;
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import copy;
from pprint import pprint;
import cPickle as pickle;
from Utilities import read_config_file


# python Create_PosNegDict_OpenAccess.py RNAi PubMedIDS_with_Extras_OpenAccess.txt ./../protein_detect_output_OpenAccess ./../JSON_SENTENCE_ASSOCS_OpenAccess


def main_body():
    global DatumKBdict; global id_to_proteins; global absId_para_sentId; global myfile_data
    global parent_location; global txt_location; global json_location; global pickle_location
    parent_location, txt_location, json_location, pickle_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Create_PosNegDict', usage='Create_PosNegDict_OpenAccess.py <label> <PubMedfilelist> <matchfiles_dir_path> <segmented_sent_dir_path>', description='Script to create dictionaries for pos and neg proteins wrt the label')
    parser.add_argument('label', help='The field we are trying to predict')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('match_dir_path', help='Path where the protein match files are located')
    parser.add_argument('sentence_dir_path', help='Path where files containing the segmented sentences are located')

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()
    
    if os.path.isfile("openaccess_posprotdict.p") and os.path.isfile("openaccess_negprotdict.p"):
        print "Found - openaccess_posprotdict.p and openaccess_negprotdict.p! Only recreating the given label."
        posprotdict = pickle.load(open("openaccess_posprotdict.p", "rb"))
        negprotdict = pickle.load(open("openaccess_negprotdict.p", "rb"))
    else:
        print "The dictionaries - openaccess_posprotdict.p and openaccess_negprotdict.p was not found in current directory. Creating new ones!"
        posprotdict = {}; negprotdict = {}

    if os.path.isfile("DatumKBProteins_UniProtId.p"):       # This dict is different from the dict used in Protein Detector
        print "Found - DatumKBProteins_UniProtId.p! Only extending the dict."
        DatumKBProtein_UniprotID_dict = pickle.load(open("DatumKBProteins_UniProtId.p", "rb"))
    else:
        print "The dictionary - DatumKBProteins_UniProtId.p was not found in current directory. Creating a new one!"
        DatumKBProtein_UniprotID_dict = {}

    label = args.label
    DatumKBdict = pickle.load( open( "pubmedid_extras.p", "rb" ) )
    id_to_proteins = pickle.load( open( "UniProtid_to_proteins.p", "rb" ) )
    absId_para_sentId = pickle.load( open( "absId_para_sentId.p", "rb" ) )

    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f1:
        for pmid in f1:
            if not pmid.strip(): continue
            pmid = pmid.strip()
            
            fname_string = '%s*' % (os.path.join(args.sentence_dir_path, pmid.strip()))
            fname = glob.glob(fname_string)            
            if len(fname) == 0:
                print "The segmented sentences file for PubMedId: ", pmid.strip(), "was not found!"
                continue                 
            myfile = fname[0].strip()

            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)    # myfile_data is a dictionary
            file_handle.close();

            if len(myfile_data["paragraphList"]) == 0:       # myfile_data["paragraphList"] is a list of dicts
                print "No sentences extracted for PubMedId: ", pmid.strip()
                continue            

            fname1_string = '%s_Paragraphs*' % (os.path.join(args.match_dir_path, pmid.strip()))
            fname1 = glob.glob(fname1_string)
            if len(fname1) == 0:
                print "The Paragraphs file for PubMedId: ", pmid.strip(), "was not found!"   # This statement will never get printed
                continue                 
            myparafile = fname1[0].strip()

            file_handle = open(myparafile, "rb")
            myparafile_data = json.load(file_handle)    # mymatchfile_data is a list of dictionaries
            file_handle.close()
            
            if pmid not in posprotdict:
                posprotdict[pmid.strip()] = {}
            if pmid not in negprotdict:            
                negprotdict[pmid.strip()] = {}

            posprotdict[pmid.strip()][label] = {}; negprotdict[pmid.strip()][label] = {}
            for elem in myparafile_data:                
                for prot in elem["matches"].keys(): # prot is a matched protein which is in lower case                                                                                                                     
                    found = False   # 'found' may be false in 2 cases - When the protein entity is not a key in the dict - DatumKBdict[pmid] or when the protein entity is present    
                    # as a key but there is no extra datum in that article which contains the protein entity and the given label as a treatment test in the same datum.
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

                    if found:                        
                        populatedict(posprotdict, label, pmid.strip(), foundprot, elem, prot)   # foundprot can be either a synonym or the exact protein that was matched
                    else:
                        populatedict(negprotdict, label, pmid.strip(), prot, elem, prot)
                        foundprot = prot

                    if foundprot not in DatumKBProtein_UniprotID_dict:  # What if 2 different proteins have the same synonym which is a DatumKB prot
                        DatumKBProtein_UniprotID_dict[foundprot] = [[prot], elem["matches"][prot]["UniProtId"]]
                    else:
                        if prot not in DatumKBProtein_UniprotID_dict[foundprot][0]:
                            DatumKBProtein_UniprotID_dict[foundprot][0].append(prot)
                            if DatumKBProtein_UniprotID_dict[foundprot][1] != elem["matches"][prot]["UniProtId"]:
                                DatumKBProtein_UniprotID_dict[foundprot][1] += ", " + elem["matches"][prot]["UniProtId"]

            continue
            allposSentences = []
            for prot in posprotdict[pmid][label]:
                for absid in posprotdict[pmid][label][prot]:
                    allposSentences.append(absid)

            allposSentences.sort()
            for prot in negprotdict[pmid][label]:
                for absid in negprotdict[pmid][label][prot].keys():
                    if absid in allposSentences:                        
                        negprotdict[pmid][label][prot].pop(absid)
            print "Control is coming here!"     # This statement should not print
            for prot in negprotdict[pmid][label].keys():
                if not negprotdict[pmid][label][prot]:                    
                    negprotdict[pmid][label].pop(prot)

    num_pmids = 0
    for pmid in posprotdict:
        if posprotdict[pmid][label]:
            num_pmids += 1
    print "No of PMIDs that have extra datums with treatment type - ", label, " are: ", num_pmids

    pickle.dump(posprotdict, open("openaccess_posprotdict.p", "wb"))

    pickle.dump(negprotdict, open("openaccess_negprotdict.p", "wb"))

    pickle.dump(DatumKBProtein_UniprotID_dict, open("DatumKBProteins_UniProtId.p", "wb"))
       

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


def populatedict(myprotdict, label, pmid, myprot, mydict, matched_prot):
    if myprot not in myprotdict[pmid][label]:
        myprotdict[pmid][label][myprot] = {}
        myprotdict[pmid][label][myprot]["AbsoluteIds"] = {}
        myprotdict[pmid][label][myprot]["Mentions"] = []
        myprotdict[pmid][label][myprot]["MatchedProtein"] = [matched_prot]  # matched_prot is the actual protein phrase that appears in the article
        myprotdict[pmid][label][myprot]["UniProtId"] = mydict["matches"][matched_prot]["UniProtId"]
    
    for myid in mydict["Neighborhood"]:
        absId = int(myid.split("_")[1])
        paraId = absId_para_sentId[pmid][absId]["paraId"]
        sentId = absId_para_sentId[pmid][absId]["sentId"]
        if absId not in myprotdict[pmid][label][myprot]["AbsoluteIds"]:
            myprotdict[pmid][label][myprot]["AbsoluteIds"][absId] = [myfile_data["paragraphList"][paraId]["sentenceList"]["@items"][sentId]]
        
    for mention in mydict["matches"][matched_prot]["Mentions"]:     # myprot could be different from matched_prot (i.e. myprot could be some synonym for matched_prot
        if mention not in myprotdict[pmid][label][myprot]["Mentions"]:
            myprotdict[pmid][label][myprot]["Mentions"].append(mention)

    if matched_prot not in myprotdict[pmid][label][myprot]["MatchedProtein"]:
        myprotdict[pmid][label][myprot]["MatchedProtein"].append(matched_prot)


if __name__ == "__main__":
    main_body()
