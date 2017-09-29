# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
from nltk import ngrams, corpus
import string
import cPickle as pickle;
import json
from pprint import pprint


# python Compare_TestDatafiles.py ./datafiles_BaseLine1/RNAi/vw_RNAi_Test_File.txt ./datafiles_BaseLine2/RNAi/vw_RNAi_Test_File.txt
# python Compare_TestDatafiles.py ./datafiles_BaseLine1/reqs/vw_reqs_Test_File.txt ./datafiles_BaseLine2/reqs/vw_reqs_Test_File.txt


def main_body():
    parser = argparse.ArgumentParser(prog='Compare_TestDatafiles', usage='Compare_TestDatafiles.py <testfile_loc1> <testfile_loc2>', description='Script to compare 2 test data files')
    parser.add_argument('testfile_loc1', help='Location of the first test data file')
    parser.add_argument('testfile_loc2', help='Location of the second test data file')

    args = parser.parse_args()

    if not os.path.isfile(args.testfile_loc1) or not os.path.isfile(args.testfile_loc2):
        print "Either one of the provided test files was not found at the given location. Please Check!"
        exit()

    testfile_dict1 = {}
    with open(args.testfile_loc1, 'rt') as f:
        for line in f:
            line = line.strip()
            elements = line.split("|")
            tag = elements[0].split()[1].strip()
            if tag in testfile_dict1:
                print "Duplicate tag:", tag, " in test file: ", args.testfile_loc1
            testfile_dict1[tag] = elements[1].strip()

    testfile_dict2 = {}
    with open(args.testfile_loc2, 'rt') as f:
        for line in f:
            line = line.strip()
            elements = line.split("|")
            tag = elements[0].split()[1].strip()
            if tag in testfile_dict2:
                print "Duplicate tag:", tag, " in test file: ", args.testfile_loc2
            testfile_dict2[tag] = elements[1].strip()

    if len(testfile_dict1.keys()) == len(testfile_dict2.keys()):
        print "Check for equal no.of test instances passed!"
    else:
        print "Check for equal no.of test instances failed!"

    if set(testfile_dict1.keys()) == set(testfile_dict2.keys()):
        print "Check for the same tags present in both test datafiles passed!"
    else:
        print "Check for the same tags present in both test datafiles failed!"

    instance_diff = False
    for tag in testfile_dict1:
        if testfile_dict1[tag] != testfile_dict2[tag]:
            print testfile_dict1[tag]
            print "------------------------------------------------------"
            print testfile_dict2[tag]
            instance_diff = True
            break
    if instance_diff:
        print "Check for same test instances failed!"
    else:
        print "Check for same test instances passed!"


if __name__ == "__main__":
    main_body()
