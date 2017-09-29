# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import os
import shutil
import subprocess
import glob
import string, ntpath
import json
import cPickle as pickle;
from Utilities import read_config_file


# python Generate_DepParses.py ./PubMedIDS_with_Extras_OpenAccess.txt ./JSON_SENTENCE_ASSOCS_OpenAccess -o
# python Generate_DepParses.py ./Total_PubMedIds_John_NoOpenAccess.txt ./JSON_SENTENCE_ASSOCS


def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Generate_DependencyParses', usage='Generate_DepParses.py <PubMedfilelist> <Extr_Sent_Dir> [-o]', description='Script to read the segmented sentences and generate dependency trees')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('Extr_Sent_Dir', help='Directory where the extracted sentences are located')
    parser.add_argument('-o', action='store_true', help='Switch to indicate whether OpenAccess articles are being processed')

    args = parser.parse_args()

    folderName = "Sentences_DepParse"
    fileListName = "PMIDList_DepParse.txt"
    outputFolder = os.path.join(args.Extr_Sent_Dir, folderName)
    ensure_dir(outputFolder)

    if args.o:
        collect_sentences_OpenAccess(os.path.join(txt_location, args.PubMedfilelist), args.Extr_Sent_Dir, outputFolder, fileListName)
    else:
        collect_sentences_NoOpenAccess(os.path.join(txt_location, args.PubMedfilelist), args.Extr_Sent_Dir, outputFolder, fileListName)

    print "\nGenerating the dependency parses for sentences extracted from every article\n"     # "sampleProps.properties" should be located in the current directory
    cmd = "java -cp $CLASSPATH -Xmx2g edu.stanford.nlp.pipeline.StanfordCoreNLP -props ./sampleProps.properties -filelist " + os.path.join(outputFolder, fileListName) + " -outputDirectory " + outputFolder
    print cmd
    my_env = os.environ.copy()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

    for line in p.stderr.readlines():  # Do NOT comment this for statement
        #print line.strip()
        pass


def collect_sentences_OpenAccess(PubMedfilelist, Extr_Sent_Dir, outputFolder, fileListName):
    fileList = []
    with open(PubMedfilelist, 'rt') as f1:
        for pmid in f1:
            pmid = pmid.strip()
            print "PMID: ", pmid
            if not pmid.strip(): continue
            fname_string1 = '%s*' % (os.path.join(Extr_Sent_Dir, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                continue
            myfile = fname[0].strip()

            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)  # myfile_data is a list of dictionaries
            file_handle.close()

            sentence_list = []
            if len(myfile_data["paragraphList"]) == 0: continue;
            for para in myfile_data["paragraphList"]:  # para is a dictionary
                if "sentenceList" not in para or "@items" not in para["sentenceList"] or len(para["sentenceList"]["@items"]) == 0:
                    continue

                for itm in para["sentenceList"]["@items"]:  # itm represents a sentence in the article
                    sentence_list.append(itm["cleansedText"])

            if len(sentence_list) > 0:
                filename = os.path.basename(myfile)
                outputfile = os.path.join(outputFolder, filename + "_SentencesOnly.txt")
                fileList.append(outputfile)

                data_fn1 = open(outputfile, 'w')
                for rc in sentence_list:
                    data_fn1.write(rc + '\n')
                data_fn1.close()

    data_fn = open(os.path.join(outputFolder, fileListName), 'w')
    for rc in fileList:
        data_fn.write(rc + '\n')
    data_fn.close()


def collect_sentences_NoOpenAccess(PubMedfilelist, Extr_Sent_Dir, outputFolder, fileListName):
    fileList = []
    with open(PubMedfilelist, 'rt') as f1:
        for pmid in f1:
            pmid = pmid.strip()
            print "PMID: ", pmid
            if not pmid.strip(): continue
            fname_string1 = '%s' % (os.path.join(Extr_Sent_Dir, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                continue
            myfile = fname[0].strip()

            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)  # myfile_data is a list of dictionaries
            file_handle.close()

            sentence_list = []
            if "sentenceList" not in myfile_data: continue;
            if "@items" not in myfile_data["sentenceList"] or len(myfile_data["sentenceList"]["@items"]) == 0:
                continue

            for itm in myfile_data["sentenceList"]["@items"]:  # itm is a dictionary. It represent one line of text in the file
                sentence_list.append(itm["cleansedText"])

            if len(sentence_list) > 0:
                filename = os.path.basename(myfile)
                outputfile = os.path.join(outputFolder, filename + "_SentencesOnly.txt")
                fileList.append(outputfile)

                data_fn1 = open(outputfile, 'w')
                for rc in sentence_list:
                    data_fn1.write(rc + '\n')
                data_fn1.close()

    data_fn = open(os.path.join(outputFolder, fileListName), 'w')
    for rc in fileList:
        data_fn.write(rc + '\n')
    data_fn.close()


def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.mkdir(file_path)


if __name__ == "__main__":
    main_body()