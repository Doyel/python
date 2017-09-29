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


# python SaveMyAss.py <Test File>   

def main_body():     
    parser = argparse.ArgumentParser(prog='SaveMyAss.py', usage='SaveMyAss.py <Test File>', description='Script to save my ass from Mark')        
    parser.add_argument('test_filename', help='Name of the test data file')
    args = parser.parse_args()

    tfile = args.test_filename
    myrecords = []

    with open(tfile, 'r') as f1:        
        for line in f1:
            if line.strip() == "": continue
            first_part = line.split("|")[0]; second_part = line.split("|")[1]
            prot_str = first_part.strip().split("--", 1)[1]
            if " " in prot_str:
                new_prot_str = prot_str.replace (" ", "__")
            else:
                new_prot_str = prot_str
            label_index = first_part.strip().split("--", 1)[0]
            new_line = label_index + "--" + new_prot_str + "|" + second_part
            myrecords.append(new_line)

    data_fn1 = open(tfile, 'w');
    for rc in myrecords:
        data_fn1.write(rc)
    data_fn1.close()


if __name__ == "__main__":
    main_body()
            
            
