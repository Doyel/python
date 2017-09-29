# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
from os import listdir
from os.path import isfile, join
import cPickle as pickle;
from Utilities import read_config_file


# python Extend_DatumKBProteins_UniProtID_dict.py ./protein_detect_output_OpenAccess_OLDFORMAT ./DatumKBProteins_UniProtId.p
# python Extend_DatumKBProteins_UniProtID_dict.py ./protein_detect_output_OLDFORMAT ./DatumKBProteins_UniProtId.p


def main_body():
    #global parent_location; global txt_location; global json_location
    #parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Extend_DatumKBProteins_UniProtID_dict', usage='Extend_DatumKBProteins_UniProtID_dict.py <matchfiles_dir_path> <uniprotdict>', description='Script to extend a dictionary that lists DatumKB proteins to their Uniprotids')
    parser.add_argument('match_dir_path', help='Path where the protein match files are located')
    parser.add_argument('uniprotdict', help='Path where the DatumKB proteins to UniprotIDs dict is located')

    args = parser.parse_args()

    if os.path.isfile(args.uniprotdict):  # This dict is different from the dict used in Protein Detector: dictionary_nameToUniProtId.p
        print "Found - DatumKBProteins_UniProtId.p! Only extending the dict."
        DatumKBProtein_UniprotID_dict = pickle.load(open("DatumKBProteins_UniProtId.p", "rb"))
    else:
        print "The dictionary - DatumKBProteins_UniProtId.p was not found in current directory."
        print "This script can only extend an existing DatumKBProteins_UniProtId.p!"
        print "It is not meant to be executed as a stand-alone script."
        # Create_PosNegDict_PdfArticles.py and Create_PosNegDict_OpenAccess.py "create" DatumKBProteins_UniProtId.p
        exit()

    onlyfiles = [join(args.match_dir_path, f) for f in listdir(args.match_dir_path) if isfile(join(args.match_dir_path, f)) and isParaFile(join(args.match_dir_path, f))]

    for myparafile in onlyfiles:
        file_handle = open(myparafile, "rb")
        myparafile_data = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
        file_handle.close()

        for elem in myparafile_data:
            for prot in elem["matches"]:
                if prot not in DatumKBProtein_UniprotID_dict:
                    if isinstance(elem["matches"][prot], dict):  # for new format of paragraph files
                        DatumKBProtein_UniprotID_dict[prot] = [[prot], elem["matches"][prot]["UniProtId"]]
                    else:   # for old format of paragraph files
                        DatumKBProtein_UniprotID_dict[prot] = [[prot], elem["matches"][prot]]
                else:
                    if prot not in DatumKBProtein_UniprotID_dict[prot][0]:
                        DatumKBProtein_UniprotID_dict[prot][0].append(prot)
                        if isinstance(elem["matches"][prot], dict):
                            if DatumKBProtein_UniprotID_dict[prot][1] != elem["matches"][prot]["UniProtId"]:
                                DatumKBProtein_UniprotID_dict[prot][1] += ", " + elem["matches"][prot]["UniProtId"]
                        else:
                            if DatumKBProtein_UniprotID_dict[prot][1] != elem["matches"][prot]:
                                DatumKBProtein_UniprotID_dict[prot][1] += ", " + elem["matches"][prot]

    pickle.dump(DatumKBProtein_UniprotID_dict, open("DatumKBProteins_UniProtId.p", "wb"))


def isParaFile(filename):
    if filename.find("_Paragraphs") > -1:
        return True
    else:
        return False


if __name__ == "__main__":
    main_body()