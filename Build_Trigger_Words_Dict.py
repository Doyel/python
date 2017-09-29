# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
from pprint import pprint
import cPickle as pickle;


# python Build_Trigger_Words_Dict.py RNAi /ua/ml-group/big-mechanism-project/trigger_words/RNAi_triggers
# python Build_Trigger_Words_Dict.py omission /ua/ml-group/big-mechanism-project/trigger_words/omission_triggers
# python Build_Trigger_Words_Dict.py KO /ua/ml-group/big-mechanism-project/trigger_words/KO_triggers


def main_body():
    global my_trigger_dict
    parser = argparse.ArgumentParser(prog='Build_Trigger_Words_Dict', usage='Build_Trigger_Words_Dict.py <label> <file_loc>', description='Script to build a dict for trigger words corresponding to individual labels')
    parser.add_argument('label', help='The field for which we want to collect the trigger words')
    parser.add_argument('file_loc', help='The location where the trigger words selected by Mark are present')

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    triggerdictfile = "trigger_words_dict.p"
    if os.path.isfile(triggerdictfile):
        my_trigger_dict = pickle.load(open(triggerdictfile, "rb"))
        print "trigger_words_dict.p was found in current dir. Recreating the trigger words for the label: ", args.label
    else:
        print "Couldn't locate the dictionary - " + triggerdictfile + "! Creating the dict anew."
        my_trigger_dict = {}

    my_trigger_dict[args.label] = []
    create_trigger_dict_new(args.label, args.file_loc)
    my_trigger_dict[args.label] = list(set(my_trigger_dict[args.label]))
    #pprint(my_trigger_dict)
    print "No. of trigger words for label: ", args.label, " \t - ", len(my_trigger_dict[args.label])
    pickle.dump(my_trigger_dict, open(triggerdictfile, "wb"))


def create_trigger_dict(label, path_triggerfile):
    global my_trigger_dict
    if not os.path.isfile(path_triggerfile):
        print "Couldn't locate the trigger words file - " + path_triggerfile
        exit()

    with open(path_triggerfile, 'rt') as f:
        for line in f:
            if line.startswith("*"):
                word = line.split(":-:")[0].strip().lstrip("*")
                if len(word) > 0:
                    my_trigger_dict[label].append(word.lower())


def create_trigger_dict_new(label, path_triggerfile):
    global my_trigger_dict
    if not os.path.isfile(path_triggerfile):
        print "Couldn't locate the trigger words file - " + path_triggerfile
        exit()

    with open(path_triggerfile, 'rt') as f:
        for line in f:
            if len(line.strip()) > 0:
                my_trigger_dict[label].append(line.strip().lower())


if __name__ == "__main__":
    main_body()