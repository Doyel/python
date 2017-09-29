# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import collections;
import string;
from pprint import pprint;
import cPickle as pickle;
from random import shuffle
from nltk.tokenize import StanfordTokenizer
from nltk.tokenize import sent_tokenize, word_tokenize

# Mark_Training_Set_PubMedIds_Jul7.txt contains all 61 OpenAccess articles that have interesting 'extras' datums in the DatumKB
# python Create_Labelled_DataFile.py <label> ./Mark_Training_Set_PubMedIds_Jul7.txt ./protein_detect_output ./pubmedid_extras.p <filename>

# Mark_Training_Set_PubMedIds_Jul7_NoOpenAccess.txt does not contain the 61 OpenAccess articles because the *Paragraphs* files
# for the 61 articles are in a different format than the *Paragraphs* files for the rest non-OpenAccess articles
# python Create_Labelled_DataFile.py <label> ./Mark_Training_Set_PubMedIds_Jul7_NoOpenAccess.txt ./protein_detect_output ./pubmedid_extras.p <filename>

def main_body():
    global pmid_num_reqs_dnreqs; global mydict; global id_to_proteins   
    parser = argparse.ArgumentParser(prog='Create_Labelled_DataFile.py', usage='Create_Labelled_DataFile.py <label> <PMIDlist> <proteins_detected_dir> <DatumKB_dict> <datafilename>', description='Script to create labelled data file')    
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('PMIDlist', help='File listing the PubMed ids to be included in the training file')        
    parser.add_argument('protein_outputdir', help='Path where files containing the detected proteins are located')    
    parser.add_argument('DatumKB_dict', help='Dictionary created from the datum KB given by Mark')
    parser.add_argument('Out_filename', help='Name of the labelled data file')
    args = parser.parse_args()
    
    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()
    
    pmid_num_reqs_dnreqs = pickle.load( open( "pmid_num_reqs_dnreqs_extras.p", "rb" ) )
    mydict = pickle.load( open( args.DatumKB_dict, "rb" ) )
    id_to_proteins = pickle.load( open( "UniProtid_to_proteins.p", "rb" ) )
                
    mypmids = []
    with open(args.PMIDlist, 'rt') as f1:        
        for pmid in f1:
            if not pmid.strip(): continue                                  
            mypmids.append(pmid.strip())

    get_proportion(mypmids)
    filename = args.Out_filename    
    create_datafile(mypmids, args.protein_outputdir, filename, args.label)


def create_datafile(mypmids, protein_outputdir, filename, label): 
    myrecords = []       
    for pmid in mypmids:
        print "Processing: ", pmid
        fname_string = '%s*' % (os.path.join(protein_outputdir, pmid))
        fname = glob.glob(fname_string)
        if len(fname) != 2:
            print "Either one of \'ProteinMatches\' or \'Paragraphs\' file is not present for PubMedId: ", pmid
            continue   
        if "Paragraphs" in fname[0]:    parafile = fname[0].strip()                  
        else:   parafile = fname[1].strip()

        file_handle = open(parafile, "rb")
        myfile = json.load(file_handle)    # myfile is a list of dictionaries
        file_handle.close();
                    
        for elem in myfile:    # elem is a dictionary
            vw_tag = elem["id"]; vw_label = "-1"
            if pmid in mydict.keys():
                for prot in elem["matches"].keys(): # prot is a matched protein which is in lower case                                                                                                                     
                    found = False   # 'found' may be false in 2 cases - When the protein entity is not a key in the dict - mydict[pmid] or when the protein entity is present    
                    # as a key but there is no extra datum in that article which contains the protein entity and the given label as a treatment test in the same datum.
                    if prot in mydict[pmid].keys():
                        if label == "reqs":
                            if len(mydict[pmid][prot]['reqstest']) > 0: found = True    # We check the length because of lines 114 and 116 in StratifyArticles.py
                        elif label == "dnreqs":
                            if len(mydict[pmid][prot]['dnreqstest']) > 0: found = True
                        else:
                            for test in mydict[pmid][prot].keys():      # test can only be 'reqstest' or 'dnreqstest'
                                if label in mydict[pmid][prot][test]:    found = True
                        foundprot = prot                    
                    else:
                        found,foundprot = find_protsyn(elem["matches"][prot], label, pmid) 
                        #if found: print foundprot, elem["id"]
                    if found:
                        vw_tag = update_tag(vw_tag, foundprot)
                        vw_label = "1"    
            if vw_label == "-1":
                 vw_tag = update_tag(vw_tag, elem["matches"].keys()[0]) 
            lex_feat_str = create_lexical_features(elem["paragraph"])
            rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str
            myrecords.append(rec)

    data_fn1 = open(filename, 'w')
    for rc in myrecords:
        data_fn1.write(rc+'\n')
    data_fn1.close()


def find_protsyn(uniprotidstr, label, pmid):
    found = False
    uniprotidlist = uniprotidstr.split(",")
    for upid in uniprotidlist:
        for syn in id_to_proteins[upid.strip()]:
            syn = syn.lower()
            if syn in mydict[pmid].keys():
                if label == "reqs":
                    if len(mydict[pmid][syn]['reqstest']) > 0: found = True     # We check the length because of lines 114 and 116 in StratifyArticles.py
                elif label == "dnreqs":
                    if len(mydict[pmid][syn]['dnreqstest']) > 0: found = True
                else:
                    for test in mydict[pmid][syn].keys():
                        if label in mydict[pmid][syn][test]:    found = True
                if found:
                    return found, syn
    return False, "NIL"


def update_tag(vw_tag, foundprot):    
    foundprot = foundprot.replace (" ", "__")  
    if foundprot not in vw_tag:         
        vw_tag = vw_tag + "--" + foundprot
    return vw_tag


def create_lexical_features(paragraph):
    #tokens = StanfordTokenizer().tokenize(paragraph.strip())
    tokens = word_tokenize(paragraph.strip())
    #for t in sent_tokenize(paragraph):
        #tokens.extend(word_tokenize(t))
    countdict = collections.Counter(tokens)
    cleanse_punct(countdict)
    mystr = ' '.join(['%s:%s' % (k,v) for k,v in countdict.iteritems()])
    return mystr


def cleanse_punct(countdict):
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)
    for word in countdict.keys():
        if (":" in word) or (len(word.strip()) == 0):
            del countdict[word]
            continue
        out = word.translate(remove_punctuation_map)
        if len(out.strip()) == 0:
            del countdict[word]


def get_proportion(mypmids):
    num_reqs = 0.0; num_dnreqs = 0.0
    for pmid in mypmids:
        if pmid not in pmid_num_reqs_dnreqs:    continue
        if "reqs" in pmid_num_reqs_dnreqs[pmid]:    num_reqs += pmid_num_reqs_dnreqs[pmid]["reqs"]
        if "dnreqs" in pmid_num_reqs_dnreqs[pmid]:    num_dnreqs += pmid_num_reqs_dnreqs[pmid]["dnreqs"]
    print "The proportion of 'reqs' extras: ", float(100 * (num_reqs / (num_reqs + num_dnreqs)))
    print "The proportion of 'dnreqs' extras: ", float(100 * (num_dnreqs / (num_reqs + num_dnreqs))) 


if __name__ == "__main__":
    main_body()
