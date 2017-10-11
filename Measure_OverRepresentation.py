# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
from nltk import ngrams, corpus
import string
import json
from pprint import pprint
import cPickle as pickle;
import nltk
from nltk.tokenize import word_tokenize
from nltk.collocations import *
from Utilities import read_config_file


# python Measure_OverRepresentation.py <label> <psg_dict_loc> <prune_freq> <num_grams> [-s]
# python Measure_OverRepresentation.py RNAi ./datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0 1 50


def main_body():
    global stopwords; global parent_location; global txt_location; global json_location; global pickle_location
    parent_location, txt_location, json_location, pickle_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Measure_OverRepresentation', usage='Measure_OverRepresentation.py <label> <psg_dict_loc> <prune_freq> <num_topngrams> [-s]', description='Script to measure overrepresentation')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for which to measure overrepresentation')
    parser.add_argument('psg_dict_loc', help='The location where the passage_dict json files are present')
    parser.add_argument('prune_freq', help='The positive ngram frequency below which ngrams are to be ignored')
    parser.add_argument('num_topngrams', help='The no. of highest scored ngrams to be printed')
    parser.add_argument('-s', action='store_true', help='Switch to control whether to ignore ngrams with stopwords')

    args = parser.parse_args()
    #stopwords = corpus.stopwords.words('english')
    stopwords = pickle.load(open(os.path.join(pickle_location, 'stop_words_dict.p'), 'rb'))

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    print "\n*******************************************************"
    if args.s:
        print "Ngrams containing stopwords will be ignored!"
        stopwordswitch = True
        stopword_str = "_WithStopwords"
    else:
        print "Ngrams containing stopwords will NOT be ignored!"
        stopwordswitch = False
        stopword_str = "_NoStopWords"
    print "*******************************************************"

    label = args.label
    unigram_dict = {}; bigram_dict = {}; trigram_dict = {}
    positive_passage_text = []

    print "\n\n"
    print "Processing for Label: ", label, " and Frequency Prune: ", args.prune_freq, "\n"
    sys.stdout.flush()

    for passagedictfile in ["openaccess_passage_dict.json", "train_passage_dict.json"]:

        fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
        if os.path.isfile(fullpath_passagedictfile):
            file_handle = open(fullpath_passagedictfile, "rb")
            passage_dict = json.load(file_handle)
            file_handle.close()
        else:
            print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
            exit()

        for pmid in passage_dict:
            for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label][myclass] is a dict
                    all_sentences_protein = []      # All sentences for current protein will be of the same class - myclass

                    for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                        all_sentences_protein.append(contg_psg_dict["textOfInterest"])      # sentences_protein is a list of lists
                        if myclass == "Pos":
                            positive_passage_text.extend(contg_psg_dict["textOfInterest"])

                    tokens = create_tokens(all_sentences_protein)   # tokens is a list of lists
                    for n in range(1, 3):
                        if n == 1:
                            create_ngrams(tokens, n, unigram_dict, myclass, stopwordswitch)
                        elif n == 2:
                            create_ngrams(tokens, n, bigram_dict, myclass, stopwordswitch)
                        #elif n == 3:
                            #create_ngrams(tokens, n, trigram_dict, myclass, stopwordswitch)

    print "Top uni-grams based on Mark's overrepresentation score:"
    print "------------------------------------------------------"
    sys.stdout.flush()
    fname = os.path.join(args.psg_dict_loc, label, "Unigram_Score_" + label + "_FrequencyPrune_" + args.prune_freq + stopword_str + ".txt")
    get_sorted_ngrams(unigram_dict, int(args.prune_freq), int(args.num_topngrams), fname)
    print "\n\n"
    print "Top bi-grams based on Mark's overrepresentation score:"
    print "------------------------------------------------------"
    sys.stdout.flush()
    fname = os.path.join(args.psg_dict_loc, label, "Bigram_Score_" + label + "_FrequencyPrune_" + args.prune_freq + stopword_str + ".txt")
    get_sorted_ngrams(bigram_dict, int(args.prune_freq), int(args.num_topngrams), fname)

    # print "\n\n"
    # print "Top tri-grams based on Mark's overrepresentation score:"
    # print "------------------------------------------------------"
    # sys.stdout.flush()
    # get_sorted_ngrams(trigram_dict, int(args.prune_freq), int(args.num_topngrams))

    # positive_tokens = create_tokens(positive_passage_text)      # positive_passage_text is a list of sentences that constitute the positive training instances
    #
    # bigram_measures = nltk.collocations.BigramAssocMeasures()
    # trigram_measures = nltk.collocations.TrigramAssocMeasures()
    #
    # print "\n\n"
    # print "Top collocated bi-grams in positive passages based on PMI score:"
    # print "---------------------------------------------------------------"
    # sys.stdout.flush()
    # finder = BigramCollocationFinder.from_words(positive_tokens)
    # finder.apply_freq_filter(int(args.prune_freq))
    # #finder.apply_word_filter(lambda w: w.lower() in stopwords)     # We remove all bigrams that contain stopwords from candidature
    # finder.nbest(bigram_measures.pmi, int(args.num_topngrams))
    # scored = finder.score_ngrams(bigram_measures.pmi)
    # sorted(scored, key=lambda x: x[1])
    # pprint(scored[:int(args.num_topngrams)])
    # sys.stdout.flush()
    #
    # print "\n\n"
    # print "Top collocated tri-grams in positive passages based on PMI score:"
    # print "----------------------------------------------------------------"
    # sys.stdout.flush()
    # finder = TrigramCollocationFinder.from_words(positive_tokens)
    # finder.apply_freq_filter(int(args.prune_freq))
    # #finder.apply_word_filter(lambda w: w.lower() in stopwords)
    # finder.nbest(trigram_measures.pmi, int(args.num_topngrams))
    # scored = finder.score_ngrams(trigram_measures.pmi)      # Assumption: scored is a list of tuples (ngram, score)
    # sorted(scored, key=lambda x: x[1])
    # pprint(scored[:int(args.num_topngrams)])
    # sys.stdout.flush()


def create_ngrams(tokens, n, mydict, myclass, stopwordswitch):    # tokens is a list of lists
    myngrams = []
    for tkns_list in tokens:    # tkns_list contains tokens from a single contiguous evidence passage
        myngrams.extend(list(ngrams(tkns_list, n)))
    myngrams = list(set(myngrams))      # Unique list of ngrams
    for each_ngram in myngrams:
        if stopwordswitch and contains_stopwords(each_ngram):  continue
        if each_ngram in mydict:
            if myclass == "Pos":
                mydict[each_ngram][0] += 1
            else:
                mydict[each_ngram][1] += 1
        else:
            if myclass == "Pos":
                mydict[each_ngram] = [1, 0]
            else:
                mydict[each_ngram] = [0, 1]


def cleanse_punct(tokens):
    clean_tokens = []
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)
    for word in tokens:
        out = word.translate(remove_punctuation_map)
        if len(out.strip()) > 0:
            clean_tokens.append(word)
    return clean_tokens


def contains_stopwords(ngram_tup):
    for word in ngram_tup:
        if isstopword(word.lower()):
            #print ngram_tup
            return True
    #print ngram_tup
    return False


def isstopword(myword):
    start_char = myword[0]
    try:
        if myword in stopwords[start_char]:
            return True
        else:
            return False
    except KeyError:
        return False


def get_sorted_ngrams(mydict, freq, n, fname):
    for each_ngram in mydict.keys():
        num_instances = float(mydict[each_ngram][0] + mydict[each_ngram][1])    # Total number of instances in which that ngram has occurred
        if mydict[each_ngram][0] <= freq:      # If an ngram occurs 1 time or less than 1 time in a positive instance then ignore that ngram
            #print "Saswati: ", mydict[each_ngram][0], "\t", num_instances
            del mydict[each_ngram]
            continue
        mydict[each_ngram].append(float(mydict[each_ngram][0]) / num_instances)

    ngrams_sorted = sorted(mydict, key=lambda k: mydict[k][2], reverse=True)
    for i in range(1, n):
        print ngrams_sorted[i], "\t", mydict[ngrams_sorted[i]]
        sys.stdout.flush()

    data_fn1 = open(fname, 'w')
    for k in ngrams_sorted:
        list_str = ", ".join(map(str, mydict[k]))
        ngram_str = ", ".join(k)
        data_fn1.write(ngram_str+"\t\t:-:\t\t"+list_str+"\n")
    data_fn1.close()


def create_tokens(contgparaList):
    tokens = []
    for sentList in contgparaList:  # sentList represents a list of sentences that constitute a contiguous passage
        tokens_contgpsg = []    # tokens_contgpsg contains tokens obtained from sentences that constitute a single contiguous passage
        for sent in sentList:
            sent = sent.lower()     # Trigger words should be case-insensitive!!!
            sent = sent.replace('-', ' ')
            tokens_contgpsg.extend(word_tokenize(sent.strip()))
        tokens_contgpsg = cleanse_punct(tokens_contgpsg)
        tokens.append(tokens_contgpsg)
    return tokens


if __name__ == "__main__":
    main_body()