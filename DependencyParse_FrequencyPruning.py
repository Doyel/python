# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
import json;
import ntpath
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import operator
import shutil
from pprint import pprint;
from Utilities import read_config_file


# python DependencyParse_FrequencyPruning.py RNAi Dependency_Features.json ./../datafiles/datafiles_WithWeights_BaseLine2_Aug31 DepParseFeatures_FrequencyCnts.json 50
# python DependencyParse_FrequencyPruning.py omission Dependency_Features.json ./../datafiles/datafiles_WithWeights_BaseLine2_Aug31 DepParseFeatures_FrequencyCnts.json 50
# python DependencyParse_FrequencyPruning.py KO Dependency_Features.json ./../datafiles/datafiles_WithWeights_BaseLine2_Aug31 DepParseFeatures_FrequencyCnts.json 50


def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='DependencyParse_FrequencyPruning', usage='DependencyParse_FrequencyPruning.py <depParse_features> <passage_dict_loc> <outputJSON> <numfeatures>', description='Script to perform frequency pruning of the dependency parse features')
    parser.add_argument('label', help='The predicate for which frequency analysis is to be done')
    parser.add_argument('depParseFeatJSON', help='JSON file that stores the dependency parse features for every sent')
    parser.add_argument('psgdictloc', help='Location of the OpenAccess and Non-OA passage dict json files')
    parser.add_argument('depParseFreqJSON', help='JSON file that stores the frequency of dependency parse features')
    parser.add_argument('numfeatures', help='No of highest frequency dependency parse features to print')


    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    label = args.label

    if os.path.isfile(os.path.join(json_location, args.depParseFreqJSON)):
        print "Found - " + args.depParseFreqJSON + "! Re-computing dep parse features frequency for the given predicate: ", label
        file_handle = open(os.path.join(json_location, args.depParseFreqJSON), "rb")
        freqdict = json.load(file_handle)
        file_handle.close()
    else:
        print "The dependency parse frequency counts file - " + args.depParseFreqJSON + " was not found in the JSON directory. Creating a new one!"
        freqdict = dict()
        freqdict["ARTICLE"] = {}
        freqdict["PROTEIN"] = {}
        freqdict["INSTANCE"] = {}

    freqdict["ARTICLE"][label] = {}
    freqdict["PROTEIN"][label] = {}
    freqdict["INSTANCE"][label] = {}

    if os.path.isfile(os.path.join(json_location, args.depParseFeatJSON)):
        print "Found - " + args.depParseFeatJSON + "! Loading the dependency parse features JSON file into memory!"
        file_handle = open(os.path.join(json_location, args.depParseFeatJSON), "rb")
        depfeat_data = json.load(file_handle)
        file_handle.close()
    else:
        print "The dictionary file - " + args.depParseFeatJSON + " was not found in the JSON directory. Exiting!"
        exit()

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

        for pmid in passage_dict:
            if not passage_dict[pmid][label]["Pos"]:
                continue
            print "PMID: ", pmid
            article_level_features = [] # Stores features extracted from all positive instances corresponding to all positive proteins in an article
            for prot in passage_dict[pmid][label]["Pos"]:  # prot is the DatumKB protein. It could be a synonym of the actual protein mentioned in the article
                #print prot
                protein_level_features = [] # Stores features extracted from all postive instances corresponding to a single positive protein
                matched_prots = passage_dict[pmid][label]["Pos"][prot]["MatchedProtein"]
                for contg_psg_dict in passage_dict[pmid][label]["Pos"][prot]["passageDetails"]:  # passage_dict[pmid][label]["Pos"][prot]["passageDetails"] is a list of dicts.
                    instance_level_features = []    # Stores features extracted from a single positive instance corresponding to a single positive protein
                    sent_ids = getSentenceIds(contg_psg_dict, fname)  # contg_psg_dict represents one training instance
                    for sentid in sent_ids:   # sent_ids contain the ids of sentences that make up a single positive training instance represented by the dict: contg_psg_dict
                        sentid = str(sentid); found = False
                        #print sentid
                        for mtchprot in matched_prots:  # There may be 2 or more matched proteins "mtchprot" for the same DatumKB protein "prot" in a single sentence
                            if mtchprot in depfeat_data[pmid][sentid]:
                                found = True    # Assumption: Each sentence will contain atleast one protein mention as "Neighborhood" is 0 in Prepare_Paragraphs
                                article_level_features.extend(depfeat_data[pmid][sentid][mtchprot].keys())    # I only consider features originating from the proteins which are positive wrt given predicate
                                protein_level_features.extend(depfeat_data[pmid][sentid][mtchprot].keys())
                                instance_level_features.extend(depfeat_data[pmid][sentid][mtchprot].keys())
                        if not found:
                            print "None of the matched proteins was found in the current sentence: " + sentid
                            exit()  # This shouldn't happen because sentid is part of a "Pos" training instance that corresponds to a positive protein wrt the given predicate

                    populate_freqcount_dict(list(set(instance_level_features)), freqdict["INSTANCE"][label])

                populate_freqcount_dict(list(set(protein_level_features)), freqdict["PROTEIN"][label])

            populate_freqcount_dict(list(set(article_level_features)), freqdict["ARTICLE"][label])

    output_dir = os.path.join(txt_location, "DepParse_FreqCnts", label)
    ensure_dir(output_dir)
    sort_features(output_dir, "ARTICLE", freqdict["ARTICLE"][label], args.numfeatures)
    sort_features(output_dir, "PROTEIN", freqdict["PROTEIN"][label], args.numfeatures)
    sort_features(output_dir, "INSTANCE", freqdict["INSTANCE"][label], args.numfeatures)

    json_fn = open(os.path.join(json_location, args.depParseFreqJSON), 'w')
    json.dump(freqdict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


def populate_freqcount_dict(uniq_features_list, mydict):
    for feat in uniq_features_list:
        if feat in mydict:
            mydict[feat] += 1
        else:
            mydict[feat] = 1


def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)


def sort_features(output_dir, mystr, mydict, numfeat):
    sorted_list = sorted(mydict.items(), key=operator.itemgetter(1), reverse=True)
    print "Printing highest frequency dependency parse features wrt " + mystr + " level"
    for i in xrange(int(numfeat)):
        print sorted_list[i]

    txtfilename = "FrequencyCounts_DepParseFeatures_" + mystr + ".txt"
    data_fn1 = open(os.path.join(output_dir, txtfilename), 'w')
    for elem in sorted_list:
        data_fn1.write(elem[0] + "\t\t:-:\t\t" + str(elem[1]) + "\n")
    data_fn1.close()


def getSentenceIds(train_instance_dict, fname):
    sent_ids = []
    if fname.startswith("openaccess"):
        for para_dict in train_instance_dict["details"]:    # train_instance_dict["details"] is a list of paragraphs
            for sent_dict in para_dict["sentenceDetails"]:  # para_dict["sentenceDetails"] is a list of positive sentences within a paragraph
                sent_ids.append(sent_dict["absoluteId"])
    else:
        for sent_dict in train_instance_dict["sentenceDetails"]:
            sent_ids.append(sent_dict["index"])
    return sent_ids


if __name__ == "__main__":
    main_body()
