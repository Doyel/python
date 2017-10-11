# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
from nltk.tokenize import sent_tokenize, word_tokenize
from nltk.tokenize import StanfordTokenizer
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import string
from collections import Counter;
from pprint import pprint;
import cPickle as pickle;
from Utilities import read_config_file

# python Frequency_Prune.py Mark_Training_Set_PubMedIds_Jul7_NoOpenAccess.txt PubMedIDS_with_Extras_OpenAccess.txt 10

def main_body():
    global stopwords; global parent_location; global txt_location; global json_location; global pickle_location
    parent_location, txt_location, json_location, pickle_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Frequency_Prune', usage='Frequency_Prune.py <PubMedIdList NOT in OA> <PubMedIdList in OA> <count>', description='Script to create a dict having words that occur less than k times')
    parser.add_argument('PubMedfilelistNOTOA', help='File listing the training PubMed ids to be processed that are NOT in OpenAccess')
    parser.add_argument('PubMedfilelistOA', help='File listing the training PubMed ids to be processed that are in OpenAccess')
    parser.add_argument('count', help='Count of words we need to throw away')

    args = parser.parse_args()
    stopwords = pickle.load(open(os.path.join(pickle_location, 'stop_words_dict.p'), 'rb'))
    pruned_filename = "pruned_words_dict_" + str(args.count) + ".p"
    allowed_filename = "allowed_words_dict_" + str(args.count) + ".p"
    allowed_alpha_filename = "allowed_words_alphabet_dict_" + str(args.count) + ".p"

    if os.path.isfile(os.path.join(pickle_location, "all_words_dict.p")):
        print "Found - all_words_dict.p!"
        myworddict = pickle.load(open(os.path.join(pickle_location, "all_words_dict.p"), "rb"))

        create_pickled_files(myworddict, pruned_filename, allowed_filename, allowed_alpha_filename, int(args.count))
    else:
        print "The dictionary - all_words_dict.p was not found in the pickles directory. Creating a new file!"
        myworddict = Counter()

        sentDirOA = os.path.join(parent_location, 'JSON_SENTENCE_ASSOCS_OpenAccess')
        with open(os.path.join(txt_location, args.PubMedfilelistOA), 'rt') as f1:
            for pmid in f1:
                if not pmid.strip(): continue
                fname_string1 = '%s*' % (os.path.join(sentDirOA, pmid.strip()))
                fname = glob.glob(fname_string1)
                if len(fname) == 0:
                    print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                    continue
                myfile = fname[0].strip()

                file_handle = open(myfile, "rb")
                myfile_data = json.load(file_handle)  # myfile_data is a list of dictionaries
                file_handle.close()

                if len(myfile_data["paragraphList"]) == 0: continue;
                print "OpenAccess PMID: ", pmid.strip()
                sys.stdout.flush()
                for para in myfile_data["paragraphList"]:  # para is a dictionary
                    if "sentenceList" not in para:  continue;
                    if "@items" not in para["sentenceList"]: continue;
                    if len(para["sentenceList"]["@items"]) == 0: continue

                    for itm in para["sentenceList"]["@items"]:  # itm is a dictionary. It represent one sentence in the current paragraph
                        #words = StanfordTokenizer().tokenize(itm["cleansedText"].strip())
                        words = word_tokenize(itm["cleansedText"].strip())
                        countdict = Counter(words)
                        myworddict = myworddict + countdict

        print "No. of unique words in the OpenAccess articles: ", len(myworddict.keys())

        sentDirNOTOA = os.path.join(parent_location, 'JSON_SENTENCE_ASSOCS')
        with open(os.path.join(txt_location, args.PubMedfilelistNOTOA), 'rt') as f1:
            for pmid in f1:
                if not pmid.strip(): continue
                fname_string1 = '%s' % (os.path.join(sentDirNOTOA, pmid.strip()))
                fname = glob.glob(fname_string1)
                if len(fname) == 0:
                    print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                    continue
                myfile = fname[0].strip()

                file_handle = open(myfile, "rb")
                myfile_data = json.load(file_handle)  # myfile_data is a dictionary
                file_handle.close()

                if "sentenceList" not in myfile_data: continue;
                if "@items" not in myfile_data["sentenceList"]: continue;  # myfile_data["sentenceList"] is a dictionary
                print "PMID: ", pmid.strip()
                sys.stdout.flush()
                for itm in myfile_data["sentenceList"]["@items"]:  # itm is a dictionary. It represent one line of text in the file
                    #words = StanfordTokenizer().tokenize(itm["cleansedText"].strip())
                    words = word_tokenize(itm["cleansedText"].strip())
                    countdict = Counter(words)
                    myworddict = myworddict + countdict

        pickle.dump(myworddict, open(os.path.join(pickle_location, "all_words_dict.p"), "wb"))
        create_pickled_files(myworddict, pruned_filename, allowed_filename, allowed_alpha_filename, int(args.count))


def getNumOfUnicodeStr(pruned_words, value):
    num_uni = 0; num_ascii = 0
    for word in pruned_words:
        if isStrUnicode(word):
            num_uni += 1
        else:
            #print word
            num_ascii += 1
        if pruned_words[word] > value:
            print "pruned_words_dict.p has tokens that occur more than ", value, " times"
            exit()
    return (num_uni, num_ascii)


def isStrUnicode(mystring):
    try:
        mystring.decode('ascii')
    except UnicodeDecodeError:
        #print mystring
        return True
    else:
        return False


def create_dict(words_count_dict, words_alphabet_dict):
    for word in words_count_dict:
        start_char = word[0]
        if start_char in words_alphabet_dict:
            words_alphabet_dict[start_char].append(word)
        else:
            words_alphabet_dict[start_char] = [word]


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


def create_pickled_files(myworddict, pruned_filename, allowed_filename, allowed_alpha_filename, count):
    print "Most common words in our vocabulary, WITHOUT any stopword pruning and lower-limit frequency pruning"
    print myworddict.most_common(50)
    print "No. of unique words in all articles (OA and Non-OA): ", len(myworddict.keys())
    freq_pruned_words = Counter({k: v for (k, v) in myworddict.items() if v <= count})
    stop_pruned_words = Counter({k: v for (k, v) in myworddict.items() if isbadword(k.lower())})
    pruned_words = freq_pruned_words + stop_pruned_words
    print "No. of pruned words: ", len(pruned_words.keys())
    # (num_uniStr, num_ascii) = getNumOfUnicodeStr(pruned_words, count)
    # print "No. of pruned words that have unicode characters in them: ", num_uniStr
    # print "No. of pruned words that DO NOT have unicode characters in them: ", num_ascii

    temp_allowed_words = {k: v for (k, v) in myworddict.items() if v > count}
    allowed_words = Counter({k: v for (k, v) in temp_allowed_words.items() if not isbadword(k.lower())})
    print "No. of allowed words: ", len(allowed_words.keys())
    print "Most common 50 words after lower-limit frequency pruning, stopword and punctuation pruning"
    print allowed_words.most_common(50)
    allowed_words_dict = {}
    create_dict(allowed_words, allowed_words_dict)
    print "No. of chars that are keys in the allowed words alphabet dict: ", len(allowed_words_dict.keys())

    pickle.dump(pruned_words, open(os.path.join(pickle_location, pruned_filename), "wb"))
    pickle.dump(allowed_words, open(os.path.join(pickle_location, allowed_filename), "wb"))
    pickle.dump(allowed_words_dict, open(os.path.join(pickle_location, allowed_alpha_filename), "wb"))


def isPunct(mychar):
    for c in set(string.punctuation):           # prot is "IL" and sentence contains words like "IL-4" and "IL-5"
        if c == mychar:
            return 1                # [UPDATE] 6/23/2017: Lets consider "-" as a punctuation Eg: DDX41-shrna
    if mychar in [u"\u2013", u"\u2014", u"\u2212", u"\u2018", u"\u2019", u"\u201d", u"\u201c"]:
        return 1
    return 0


if __name__ == "__main__":
    main_body()
