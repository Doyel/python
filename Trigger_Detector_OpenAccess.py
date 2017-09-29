# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import re;
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os
import string
import operator;
from nltk import ngrams
from nltk.tokenize import word_tokenize
from pprint import pprint;
import cPickle as pickle;
from Utilities import read_config_file, check_for_unicode

# PubMedIDS_with_Extras_OpenAccess.txt contains the 61 OpenAccess articles. The extracted sentences files for these articles have a different format.
# python Trigger_Detector_OpenAccess.py RNAi ./PubMedIDS_with_Extras_OpenAccess.txt ./JSON_SENTENCE_ASSOCS_OpenAccess ./trigger_words_dict.p  -c
# python Trigger_Detector_OpenAccess.py reqs ./PubMedIDS_with_Extras_OpenAccess.txt ./JSON_SENTENCE_ASSOCS_OpenAccess ./trigger_words_dict.p

def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Trigger_Detector_OpenAccess', usage='Trigger_Detector_OpenAccess.py <PubMedfilelist> <rootDir> <trigg_dict_loc> [-c]', description='Script to read the segmented sentences and recognize trigger words - ONLY for OpenAccess articles')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for trigger word detection')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('rootDir', help='Directory where the files containing individual sentences are located')
    parser.add_argument('trigg_dict_loc', help='Location of the word list that contains Mark\'s trigger words')
    parser.add_argument('-c', action='store_true')      # For a clean build of TriggerDetails dict key

    args = parser.parse_args()
    trigger_dict = pickle.load(open(args.trigg_dict_loc, "rb"))
    label = args.label

    if args.c:
        print "Performing a clean build of the TriggerDetails dict key for each sentence"

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if label not in allowed_labels:
        print "The given label is not allowed: ", label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f1:
        for pmid in f1:
            pmid = pmid.strip()
            print "PMID: ", pmid
            if not pmid.strip(): continue
            fname_string1 = '%s*' % (os.path.join(args.rootDir, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                continue
            myfile = fname[0].strip()

            filename = os.path.basename(myfile)
            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)  # myfile_data is a list of dictionaries
            file_handle.close()

            if len(myfile_data["paragraphList"]) == 0: continue;
            for para in myfile_data["paragraphList"]:  # para is a dictionary
                if "sentenceList" not in para or "@items" not in para["sentenceList"] or len(para["sentenceList"]["@items"]) == 0:
                    continue

                for itm in para["sentenceList"]["@items"]:  # itm represents a sentence in the article
                    if args.c and "TriggerDetails" in itm:
                        itm.pop("TriggerDetails", None)

                    mysent = itm["sentenceText"]
                    if check_for_unicode(mysent):
                        mysent = mysent.decode("utf-8")
                    mysent = mysent.lower()             # Converting the sentence to lowercase
                    mysent = mysent.replace('-', ' ')   # Please note that this substitution does not affect char offsets
                    tokens = word_tokenize(mysent)
                    tokens_spans = get_span_tokens(mysent, tokens)
                    clean_tokens_spans = cleanse_punct(tokens_spans, pmid)

                    for word in trigger_dict[label]:    # Trigger words are always in lowercase
                        temp = {}
                        temp["triggerWord"] = word
                        temp["label"] = label
                        found_idx = searchList(word, clean_tokens_spans)
                        if len(found_idx) > 0:
                            for idx in found_idx:
                                if "triggerCharOffset" not in temp:
                                    temp["triggerCharOffset"] = [(clean_tokens_spans[idx][1][0] + itm["charOffset"], clean_tokens_spans[idx][1][1] + itm["charOffset"], "rgb(239, 17, 255);")]
                                else:
                                    temp["triggerCharOffset"].append((clean_tokens_spans[idx][1][0] + itm["charOffset"], clean_tokens_spans[idx][1][1] + itm["charOffset"], "rgb(239, 17, 255);"))  # A sentence may contain multiple instances of the same trigger word

                            if "TriggerDetails" not in itm:
                                itm["TriggerDetails"] = [temp]
                            else:
                                itm["TriggerDetails"].append(temp)  # A sentence may contain more than one trigger words

            myfile_write_fname = os.path.join(args.rootDir, filename)
            os.remove(myfile_write_fname)
            json_fn = open(myfile_write_fname, 'w')
            json.dump(myfile_data, json_fn, indent=4, ensure_ascii=False)
            json_fn.close()


def cleanse_punct(tokens_spans, pmid):
    clean_tokens = []
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)
    for word_tuple in tokens_spans:
        out = word_tuple[0].translate(remove_punctuation_map)
        if len(out.strip()) > 0:
            clean_tokens.append(word_tuple)
    return clean_tokens


def get_span_tokens(sent, tokens):
    span_tokens = []; offset = 0
    for token in tokens:
        if token == "``" or token == "''":
            token = unicode("\"")
        offset = sent.find(token, offset)
        span_tokens.append((token, (offset, offset + len(token))))
        assert token == sent[offset:offset + len(token)]
        offset += len(token)
    return span_tokens


def searchList(word, clean_tokens_spans):
    return [i for i, x in enumerate(clean_tokens_spans) if x[0] == word]


if __name__ == "__main__":
    main_body()
