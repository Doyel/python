# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno
import glob
import shutil, ntpath
import json
import subprocess
import cPickle as pickle;
from random import shuffle
from Create_MLFile import create_lexical_features, pickle_init

# This script was written due 2 an SRI request for LOOCV scores of extras predicates for evaluation purposes - 9/14/2017
# python Create_PredictionFile_LOOCV.py KO ./../SRI_workDir/extras_WorkDir/datafiles/LOOCV_workdir ./../SRI_workDir/extras_WorkDir/datafiles/train_passage_dict.json vw_subject_Training_AllPredictions_bfgs_L2-01.txt

def main_body():
    global passage_dict
    parser = argparse.ArgumentParser(prog='Create_PredictionFile_LOOCV.py', usage='Create_PredictionFile_LOOCV.py <label> <workdir> <passagedictPath> <output_fname> <allowedWords>', description='Script to create labelled ML file')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('workdir', help='The directory in which to store the ML data files')
    parser.add_argument('passagedictPath', help='The directory in which the passage dict JSON file is located')
    parser.add_argument('output_predictfile', help='The prediction file to be generated from the training passage dict')
    parser.add_argument('allowedWords', help='Pickle file that lists all words allowed to be lexical features')

    args = parser.parse_args()
    pickle_init(args.allowedWords)

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    label = args.label
    passagedictfile = args.passagedictPath
    fname = ntpath.basename(args.passagedictPath).strip()

    ensure_dir(os.path.join(args.workdir, label))
    curr_workdir = os.path.join(args.workdir, label)

    expected_fname = "train_passage_dict.json"
    if fname != expected_fname:
        print "Incorrect passage dict file was provided! Expected file is: ", expected_fname
        exit()

    if os.path.isfile(passagedictfile):
        print "Found - " + fname + "! Everything is perfect in the world!"
        sys.stdout.flush()
        file_handle = open(passagedictfile, "rb")
        passage_dict = json.load(file_handle)
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + fname + " in the given directory! Please give the correct path."
        exit()

    method = "bfgs"
    method_param = " --bfgs --mem 20"
    l2_val = 10
    l2_str = "_L2-" + getstr(l2_val)
    l2_cmd = " --l2 " + str(l2_val)

    all_predictions = []
    test_file = os.path.join(curr_workdir, "vw_" + label + "_Test_File.txt")
    train_file = os.path.join(curr_workdir, "vw_" + label + "_Training_File.txt")
    model_file = os.path.join(curr_workdir, label + "_predictor_model_" + method + l2_str + ".vw")
    predictions_file = os.path.join(curr_workdir, "vw_" + label + "_predictions_" + method + l2_str + ".txt")

    for test_pmid in passage_dict.keys():
        clean_files(curr_workdir, "*.cache")
        clean_files(curr_workdir, "*.vw")
        clean_files(curr_workdir, "*_predictions_*")
        test_records = generate_instances(test_pmid, label)
        print "PMID: ", test_pmid
        sys.stdout.flush()
        train_records = []
        for pmid in passage_dict:
            if pmid == test_pmid:
                continue
            train_records.extend(generate_instances(pmid, label))

        shuffle(test_records); shuffle(train_records)
        writefile(test_file, test_records)
        writefile(train_file, train_records)

        cmd = "vw -d " + train_file + " -c --passes 100 -f " + model_file + " --loss_function logistic" + method_param + l2_cmd
        run_command(cmd)

        cmd = "vw -d " + test_file + " -t -i " + model_file + " -p /dev/stdout | /ua/ml-group/big-mechanism-project/vowpal_wabbit/utl/logistic -0 > " + predictions_file
        run_command(cmd)

        readfile(all_predictions, predictions_file)

    writefile(os.path.join(curr_workdir, args.output_predictfile), all_predictions)


def readfile(records, filepath):
    with open(filepath, 'rt') as f1:
        for line in f1:
            records.append(line.strip())


def writefile(filepath, records):
    silentremove(filepath)
    data_fn = open(filepath, 'w')
    for rc in records:
        data_fn.write(rc + '\n')
    data_fn.close()


def clean_files(fold_dir, pattern):
    for cache_file in glob.glob(os.path.join(fold_dir, pattern)):
        #print "Removing file: ", cache_file
        os.remove(cache_file)


def ensure_dir(file_path):
    if os.path.exists(file_path):
        shutil.rmtree(file_path)
    os.mkdir(file_path)


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:     # errno.ENOENT = no such file or directory
            raise   # re-raise exception if a different error occured


def getstr(myval):
    mystr = str(myval)
    if "." in mystr:
        mystr = mystr.replace(".", "")
    return mystr


def run_command(cmd):
    my_env = os.environ.copy()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

    for line in p.stderr.readlines():  # Do NOT comment this for statement
        #print line.strip()
        pass


def generate_instances(pmid, label):
    myrecords = []
    for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
        for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label][myclass] is a dict
            for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                vw_label = contg_psg_dict["class"]  # contg_psg_dict represents one training instance
                vw_tag = contg_psg_dict["tag"]
                lex_feat_str = create_lexical_features(contg_psg_dict["textOfInterest"])
                if len(lex_feat_str) == 0:
                    continue
                if "weight" in contg_psg_dict:
                    rec = vw_label + " " + contg_psg_dict["weight"] + " " + vw_tag + "|LexicalFeatures " + lex_feat_str  # + " |ProteinFeature " + prot
                else:
                    rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str  # + " |ProteinFeature " + prot
                myrecords.append(rec)
    return myrecords


if __name__ == "__main__":
    main_body()