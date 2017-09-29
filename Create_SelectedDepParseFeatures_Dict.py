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
from Utilities import read_config_file, silentremove


# python Create_SelectedDepParseFeatures_Dict.py RNAi Dependency_Features.json DepParseFeatures_FrequencyCnts.json ARTICLE 50 Selected_DepParseFeatures.json -c
# python Create_SelectedDepParseFeatures_Dict.py omission Dependency_Features.json DepParseFeatures_FrequencyCnts.json ARTICLE 50 Selected_DepParseFeatures.json
# python Create_SelectedDepParseFeatures_Dict.py KO Dependency_Features.json DepParseFeatures_FrequencyCnts.json ARTICLE 50 Selected_DepParseFeatures.json


def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Create_SelectedDepParseFeatures_Dict',
                                     usage='Create_SelectedDepParseFeatures_Dict.py <label> <depParse_features> <depParse_FreqCnt> <granularity> <minfreq> <outputJSON> [-c]',
                                     description='Script to built depparse features dict having user given frequency counts')
    parser.add_argument('label', help='The predicate for which frequency analysis is to be done')
    parser.add_argument('depParseFeatJSON', help='JSON file that stores the dependency parse features for every sent')
    parser.add_argument('depParseFreqJSON', help='JSON file that stores the frequency of dependency parse features')
    parser.add_argument('granularity', help='Granularity of the frequency counts to be used')
    parser.add_argument('minfreq', help='Min frequency count of dep parse features that are to be selected')
    #parser.add_argument('maxfreq', help='Max frequency count of dep parse features that are to be selected')
    parser.add_argument('outputJSON', help='JSON file that stores the selected dep parse features with provided freq')
    parser.add_argument('-c', action='store_true')  # For a clean build of the output json file

    args = parser.parse_args()

    if args.c:
        print "Performing a clean build of the selected dependency path features json file"
        silentremove(os.path.join(json_location, args.outputJSON))

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    label = args.label

    allowed_granularity = ['ARTICLE', 'PROTEIN', 'INSTANCE']
    if args.granularity not in allowed_granularity:
        print "The given granularity level is not allowed: ", args.granularity, ". It may only be any one of \'ARTICLE\', \'PROTEIN\', \'INSTANCE\'"
        print "Please try again!"
        exit()

    granularity = args.granularity
    mincnt = int(args.minfreq)
    #maxcnt = int(args.maxfreq)

    if os.path.isfile(os.path.join(json_location, args.depParseFeatJSON)):
        print "Found - " + args.depParseFeatJSON + "! The world is beautiful!"
        file_handle = open(os.path.join(json_location, args.depParseFeatJSON), "rb")
        featdict = json.load(file_handle)
        file_handle.close()
    else:
        print "The dependency parse features file - " + args.depParseFeatJSON + " was not found in the JSON directory. Exiting!"
        exit()

    if os.path.isfile(os.path.join(json_location, args.depParseFreqJSON)):
        print "Found - " + args.depParseFreqJSON + "! Everything is perfect in the world!"
        file_handle = open(os.path.join(json_location, args.depParseFreqJSON), "rb")
        freqdict = json.load(file_handle)
        file_handle.close()
    else:
        print "The dependency parse frequency counts file - " + args.depParseFreqJSON + " was not found in the JSON directory. Exiting!"
        exit()

    if os.path.isfile(os.path.join(json_location, args.outputJSON)):
        print "Found - " + args.outputJSON + "! Re-SELECTing the dep parse features for the given predicate: ", label
        file_handle = open(os.path.join(json_location, args.outputJSON), "rb")
        selecteddict = json.load(file_handle)
        file_handle.close()
    else:
        print "The dependency parse file containing the selected features - " + args.outputJSON + " was not found in the JSON directory. Creating a new one!"
        selecteddict = dict()

    selecteddict[label] = {}

    selected_features = set()
    for feat in freqdict[granularity][label]:
        if freqdict[granularity][label][feat] >= mincnt:
            #print feat, "\t", freqdict[granularity][label][feat]
            selected_features.add(feat)
    print "No. of features selected within the provided frequency range: ", len(selected_features)

    for pmid in featdict:
        # print "PMID: ", pmid
        for sentid in featdict[pmid]:
            for prot in featdict[pmid][sentid]:
                prot_features = set(featdict[pmid][sentid][prot].keys())
                # print "For SentId: ", sentid, " Features extracted from current protein: ", prot, "\n", prot_features
                common_features = selected_features.intersection(prot_features)
                # print "Common features: ", common_features
                if len(common_features) > 0:
                    if pmid not in selecteddict[label]:
                        selecteddict[label][pmid] = {}

                    if sentid not in selecteddict[label][pmid]:
                        selecteddict[label][pmid][sentid] = {}

                    selecteddict[label][pmid][sentid][prot] = list(common_features)
                #exit()

    json_fn = open(os.path.join(json_location, args.outputJSON), 'w')
    json.dump(selecteddict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


if __name__ == "__main__":
    main_body()
