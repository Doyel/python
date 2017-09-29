# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
import sys
import glob
import json
reload(sys)
sys.setdefaultencoding("utf-8")
import os
import collections
from bs4 import BeautifulSoup
from bs4 import NavigableString
from bs4 import Tag
from pprint import pprint
import cPickle as pickle
import copy
import difflib
import string
from string import punctuation
import nltk.data
import re

# python Group_Passages.py ./PubMedIDS_with_Extras_OpenAccess.txt ./vw_RNAi_Training_File_Jul7.txt ./protein_detect_output ./downloaded_articles
# python Group_Passages.py ./Test_PMIDs.txt ./vw_RNAi_Training_File_Jul7.txt ./protein_detect_output ./downloaded_articles

def main_body():
    global sent_detector; global reg_expr; global spanctr; global problem_sentences; global manual_sentences_match; global ignore_sentences
    parser = argparse.ArgumentParser(prog='Group_Passages', usage='Group_Passages.py <PubMedfilelist> <training_file> <proteins_detected_dir> <downloaded_articles_dir>', description='Script to find and group adjacent passages')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids that are in OpenAccess and are also in John list')
    parser.add_argument('training_file', help='Training File')
    parser.add_argument('protein_outputdir', help='Path where files containing the detected proteins are located')
    parser.add_argument('downloaded_articles_dir', help='Path where openaccess articles were downloaded')
    args = parser.parse_args()

    sent_detector = nltk.data.load('tokenizers/punkt/english.pickle')
    reg_expr = re.compile('^[S]?[0-9][A-Za-z]?\W*')

    if os.path.isfile("prot_passages.p"):
        prot_passages = pickle.load(open("prot_passages.p", "rb"))
    else:
        print "The dictionary - prot_passages was not found in current directory. Creating it!"
        openaccess = []; id_passages = {}
        with open(args.PubMedfilelist, 'rt') as f:
            for pmid in f:
                if not pmid.strip() or pmid.startswith("#"): continue

                fname_string = '%s*Paragraphs*' % (os.path.join(args.protein_outputdir, pmid.strip()))
                fname = glob.glob(fname_string)
                if len(fname) == 0:
                    # print "No proteins were detected for PubMedId: ", pmid.strip()
                    continue

                file_handle = open(fname[0].strip(), "rb")
                myfile = json.load(file_handle)  # myfile is a list of dictionaries
                file_handle.close()

                id_passages[pmid.strip()] = {}
                for elem in myfile:
                    mystring = elem["paragraph"].replace('\n', ' ').replace('\r', ' ')
                    id_passages[pmid.strip()][elem["id"]] = " ".join(mystring.strip(' \t\n\r').split())

                openaccess.append(pmid.strip())
        print "In the Datum KB, no. of articles that are in OpenAccess and also in John's list: ", len(openaccess)

        print "Collecting all training instances having a label of 1"
        prot_passages = {}
        with open(args.training_file, 'rt') as f1:
            for line in f1:
                label_tag = line.split("|")[0]
                label = label_tag.split()[0].strip()
                id_prots = label_tag.split()[1].split("--")
                pmid = id_prots[0].strip().split("_")[0].strip()
                id = int(id_prots[0].strip().split("_")[1].strip())
                prots = list(set(id_prots[1:]))

                if pmid in openaccess and label == "1":
                    if pmid not in prot_passages:
                        prot_passages[pmid] = {}
                    for protein in prots:
                        if protein not in prot_passages[pmid]:
                            prot_passages[pmid][protein] = {}   # prot_passages[pmid][protein] is a dictionary
                        prot_passages[pmid][protein][id] = id_passages[pmid][pmid + "_" + str(id)]
        #pprint(prot_passages["21107320"]); exit()


        print "Grouping sentences into paragraphs"


if __name__ == "__main__":
    main_body()
