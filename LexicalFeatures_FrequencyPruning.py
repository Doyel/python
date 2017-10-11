# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
import json;
import ntpath
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import cPickle as pickle;
import string
from nltk.tokenize import word_tokenize
from collections import Counter;
from pprint import pprint;
from Utilities import read_config_file
from Frequency_Prune import isPunct, create_dict


# python LexicalFeatures_FrequencyPruning.py 10 ./../datafiles/datafiles_WithWeights_BaseLine1_Oct15_NoFrequencyPrune 50
# python LexicalFeatures_FrequencyPruning.py 20 ./../datafiles/datafiles_WithWeights_BaseLine1_Oct15_NoFrequencyPrune 50
# python LexicalFeatures_FrequencyPruning.py 40 ./../datafiles/datafiles_WithWeights_BaseLine1_Oct15_NoFrequencyPrune 50


def main_body():
    global stopwords; global parent_location; global txt_location; global json_location; global pickle_location
    parent_location, txt_location, json_location, pickle_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='LexicalFeatures_FrequencyPruning', usage='LexicalFeatures_FrequencyPruning.py <count> <passage_dict_loc>', description='Script to perform frequency analysis of lexical features')
    parser.add_argument('count', help='We will throw away words whose frequency is below count')
    parser.add_argument('psgdictloc', help='Location of the OpenAccess and Non-OA passage dict json files')
    parser.add_argument('numdisplay', help='No. of highest frequency words to display')

    args = parser.parse_args()
    allowed_labels = ['RNAi']     #allowed_labels = ['RNAi', 'KO', 'omission']

    stopwords = pickle.load(open(os.path.join(pickle_location, 'stop_words_dict.p'), 'rb'))
    pruned_filename = "pruned_words_dict_" + str(args.count) + ".p"
    allowed_filename = "allowed_words_dict_" + str(args.count) + ".p"
    allowed_alpha_filename = "allowed_words_alphabet_dict_" + str(args.count) + ".p"

    if os.path.isfile(os.path.join(pickle_location, "all_words_dict.p")):
        print "Found - all_words_dict.p!"
        myworddict = pickle.load(open(os.path.join(pickle_location, "all_words_dict.p"), "rb"))

        #create_pickled_files(myworddict, pruned_filename, allowed_filename, allowed_alpha_filename, int(args.count))
    else:
        print "The dictionary - all_words_dict.p was not found in the pickles directory. Creating a new file!"
        myworddict = dict()
        myworddict["ARTICLE"] = {}
        myworddict["INSTANCE"] = {}

        for fname in ["train_passage_dict.json", "openaccess_passage_dict.json"]:
            passagedictfile = os.path.join(args.psgdictloc, fname)
            if os.path.isfile(passagedictfile):
                print "Found - " + fname + "! Loading it into memory!"
                file_handle = open(passagedictfile, "rb")
                passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
                file_handle.close()
            else:
                print "Couldn't locate the dictionary - " + fname + " in the given location! Please give the correct path."
                exit()

            for label in allowed_labels:
                print "Label: ", label
                if label not in myworddict["ARTICLE"]:
                    myworddict["ARTICLE"][label] = Counter()
                if label not in myworddict["INSTANCE"]:
                    myworddict["INSTANCE"][label] = Counter()
                for pmid in passage_dict:
                    print "PMID: ", pmid
                    article_level_tokens = [] # Stores features extracted from all instances corresponding to all proteins in an article
                    for myclass in passage_dict[pmid][label]:
                        for prot in passage_dict[pmid][label][myclass]:  # prot is the DatumKB protein. It could be a synonym of the actual protein mentioned in the article
                            for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label]["Pos"][prot]["passageDetails"] is a list of dicts.
                                instance_level_tokens = create_lexical_tokens(contg_psg_dict["textOfInterest"])    # Stores features extracted from a single instance corresponding to a single protein
                                article_level_tokens.extend(instance_level_tokens)
                                instcountdict = Counter(list(set(instance_level_tokens)))   # All tokens will have a count of 1 within instcountdict
                                myworddict["INSTANCE"][label] = myworddict["INSTANCE"][label] + instcountdict
                    artcountdict = Counter(list(set(article_level_tokens)))     # All tokens will have a count of 1 within artcountdict
                    myworddict["ARTICLE"][label] = myworddict["ARTICLE"][label] + artcountdict
        pickle.dump(myworddict, open(os.path.join(pickle_location, "all_words_dict.p"), "wb"))

    create_pickled_files(myworddict, pruned_filename, allowed_filename, allowed_alpha_filename, int(args.count))


def create_lexical_tokens(sentList):
    tokens = []
    for sent in sentList:
        tokens.extend(word_tokenize(sent.strip()))
    return tokens


def create_pickled_files(myworddict, pruned_filename, allowed_filename, allowed_alpha_filename, count):
    pruned_words = dict(); allowed_words = dict(); allowed_words_dict = dict()
    for granularity in myworddict:
        pruned_words[granularity] = {}; allowed_words[granularity] = {}; allowed_words_dict[granularity] = {}
        print "Granularity: ", granularity
        for label in myworddict[granularity]:
            print "Label: ", label
            # pruned_words[granularity][label] = {}; allowed_words[granularity][label] = {}
            print "Most common words in our vocabulary, WITHOUT any stopword pruning and lower-limit frequency pruning"
            print myworddict[granularity][label].most_common(50)
            print "No. of unique words in all articles (OA and Non-OA): ", len(myworddict[granularity][label].keys())
            freq_pruned_words = {k: v for (k, v) in myworddict[granularity][label].items() if v <= count}
            stop_pruned_words = {k: v for (k, v) in myworddict[granularity][label].items() if isbadword(k.lower())}
            pruned_words[granularity][label] = merge_dicts(freq_pruned_words, stop_pruned_words)
            print "No. of pruned words: ", len(pruned_words[granularity][label].keys())

            temp_allowed_words = {k: v for (k, v) in myworddict[granularity][label].items() if v > count}
            allowed_words[granularity][label] = {k: v for (k, v) in temp_allowed_words.items() if not isbadword(k.lower())}
            print "No. of allowed words: ", len(allowed_words[granularity][label].keys())
            print "Most common 50 words after lower-limit frequency pruning, stopword and punctuation pruning"
            print allowed_words[granularity][label].most_common(50)
            allowed_words_dict[granularity][label] = {}
            create_dict(allowed_words[granularity][label], allowed_words_dict[granularity][label])
            print "No. of chars that are keys in the allowed words alphabet dict: ", len(allowed_words_dict[granularity][label].keys())

    pickle.dump(pruned_words, open(os.path.join(pickle_location, pruned_filename), "wb"))
    pickle.dump(allowed_words, open(os.path.join(pickle_location, allowed_filename), "wb"))
    pickle.dump(allowed_words_dict, open(os.path.join(pickle_location, allowed_alpha_filename), "wb"))


def isbadword(myword):
    if myword in ["fig", "figure", "a"]:  # Assumption: There is no positive protein in DatumKB "Fig"
        return True
    if len(myword) == 1 and isPunct(myword):
        return True
    if myword.isnumeric():
        return True
    start_char = myword[0]
    try:
        if myword in stopwords[start_char]:
            return True
    except KeyError:
        return False
    return False


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


if __name__ == "__main__":
    main_body()
