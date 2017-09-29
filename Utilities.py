# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import json;
import shutil
import subprocess
import string
import ntpath
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
from collections import Counter;
import errno, os
from nltk.tokenize import word_tokenize


def tokenizeOnHyphen(tokens):
    new_tokens = tokenizeDash(tokens, "/")
    new_tokens = tokenizeDash(new_tokens, "-")
    new_tokens = tokenizeDash(new_tokens, u"\u2212")
    new_tokens = tokenizeDash(new_tokens, u"\u2014")
    new_tokens = tokenizeDash(new_tokens, u"\u2013")
    return new_tokens


def tokenizeDash(tokens, ch_dash):
    new_tokens = []
    for tok in tokens:
        if len(tok) > 1 and ch_dash in tok:
            my_arr = tok.split(ch_dash)
            for i, elem in enumerate(my_arr):
                if len(elem) > 0:
                    new_tokens.append(elem)
                if i < len(my_arr) - 1:
                    new_tokens.append(ch_dash)
        else:
            new_tokens.append(tok)
    return new_tokens


def create_lexical_features(sentList):
    #tokens = StanfordTokenizer().tokenize(paragraph.strip())
    tokens = []
    for sent in sentList:
        tokens.extend(word_tokenize(sent.strip()))
    countdict = Counter(tokens)
    cleanse_punct(countdict)
    #keep_allowed_words_count(countdict)
    if not countdict:   # If all tokens are deleted due to frequency pruning
        return ""
    mystr = ' '.join(['%s:%s' % (k, v) for k, v in countdict.iteritems()])
    return mystr


def cleanse_punct(countdict):
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)
    for word in countdict.keys():
        if (":" in word) or ("|" in word) or (len(word.strip()) == 0):
            del countdict[word]
            continue
        out = word.translate(remove_punctuation_map)
        if len(out.strip()) == 0:
            del countdict[word]
            continue


def read_config_file(config_path):
    if not os.path.isfile(config_path):
        filename = ntpath.basename(config_path)
        print "Could not locate the configuration file - " + filename + "!"
        print "Please place the above file in the current directory and try again! Exitting!!!"
        exit()
    with open(config_path, 'rt') as f1:
        for line in f1:
            parts = line.split(":")
            if parts[0].strip() == "PARENT":
                parent_path = parts[1].strip()
            if parts[0].strip() == "JSON":
                json_folder = parts[1].strip()
                json_path = os.path.join(parent_path, json_folder)
            if parts[0].strip() == "TXT":
                txt_folder = parts[1].strip()
                txt_path = os.path.join(parent_path, txt_folder)
    return parent_path, txt_path, json_path


def check_for_unicode(mystring):
    try:
        mystring.decode('ascii')
    except UnicodeDecodeError:
        # Unicode and we need to decode
        #print mystring
        return True
    except UnicodeEncodeError:
        # Already decoded, hence no need to decode
        print "Control is coming here!"
        return False
    else:
        # String is not unicode
        return False


def writefile(filepath, records):
    data_fn = open(filepath, 'w')
    for rc in records:
        data_fn.write(rc + '\n')
    data_fn.close()


def readfile(records, filepath):
    with open(filepath, 'rt') as f1:
        for line in f1:
            records.append(line.strip())


def ensure_dir(file_path):
    if os.path.exists(file_path):
        shutil.rmtree(file_path)
    os.mkdir(file_path)


def run_command(cmd):
    #print "-----------------------------------------------------------------------------------------------------------"
    #print cmd
    #print "-----------------------------------------------------------------------------------------------------------"
    my_env = os.environ.copy()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

    for line in p.stderr.readlines():  # Do NOT comment this for statement
        #print line.strip()
        pass


def get_AUPR(cmd):
    #print "-----------------------------------------------------------------------------------------------------------"
    #print cmd
    #print "-----------------------------------------------------------------------------------------------------------"
    my_env = os.environ.copy()
    #print "Saswati: ", my_env["JAVA_JARS_HOME"]
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

    for line in p.stdout.readlines():  # Do NOT comment this for statement
        # print line.strip()
        if line.strip().startswith("Area Under the Curve for Precision - Recall"):
            auc = line.strip().split()[-1]
    return float(auc)


def silentremove(filepath):
    try:
        os.remove(filepath)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured
