# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
import glob;
import json;
import unicodedata
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import cPickle as pickle;
from Utilities import read_config_file, check_for_unicode


# Full path to the UnicodeData.txt: /ua/ml-group/big-mechanism-project/PLEIADES/UserGuides/UnicodeData.txt
# python Clean_InvalidChars.py ./Total_PubMedIds_John_NoOpenAccess.txt ./JSON_SENTENCE_ASSOCS ../UserGuides/UnicodeData.txt


def main_body():
    global category_chars
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Clean_InvalidChars', usage='Clean_InvalidChars.py <PubMedfilelist> <sent_dir> <unicode_csv>', description='Script to clean invalid chars from SRI supplied segmented sentences files')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('Extr_Sent_Dir', help='Directory where the extracted sentences for all articles are located')
    parser.add_argument('unicode_csv', help='Path to where the unicode csv file is located')

    args = parser.parse_args()

    #get_UnicodeCategory(args.unicode_csv, "Cc")
    #get_UnicodeCategory(args.unicode_csv, "Cf")
    #problem_chars_dict = {
    #    "10023673": [u"\x80"],
    #    "10092233": [u"\u00AB"]
    #}
    problem_chars_dict = pickle.load(open("NonOA-BadChars.p", 'rb'))
    quote_strings = [u"\xad", u"\x80", u"\x8a", u"\x8b", u"\x87", u"\u2019\u2019-\u2019~", u"\u2018\u201c", u"\u201c\u2019", u"\u2019\u201c", u"\u2018\u201d", u"\u201d\u2018", u"\u2019\u201d", u"\u201d\u201d", u"\u201d\u201c"]

    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f1:
        for pmid in f1:
            pmid = pmid.strip()         #pmid = "12356687"
            sentence_data = getSentencesFile(args.Extr_Sent_Dir, pmid)
            if not sentence_data:
                continue  # No sentences were extracted for this article!
            print "PMID: ", pmid
            for itm in sentence_data["sentenceList"]["@items"]:  # itm is a dictionary. It represent one line of text in the file
                mysent = itm["sentenceText"]

                if check_for_unicode(mysent):
                    mysent = mysent.decode("utf-8")
                    for qstr in quote_strings:
                        if qstr in mysent:
                            mysent = mysent.replace(qstr, "")

                    if pmid in problem_chars_dict:
                        for inv_char in problem_chars_dict[pmid]:
                            if inv_char in mysent:
                                print "Replacing: ", inv_char.encode("unicode_escape")
                                print "Sentence Id: ", itm["index"]
                                mysent = mysent.replace(inv_char, "")

                    #subst_sent = mysent.replace(u"\u00a4", "_")
                    #print_categories(mysent)
                    #detect_UnicodeChar(mysent)

                itm["sentenceText"] = mysent.encode("utf-8")
                itm["cleansedText"] = " ".join(itm["sentenceText"].strip(' \t\n\r').split())
                if "ProteinDetails" in itm:
                    itm.pop("ProteinDetails", None)
                if "TriggerDetails" in itm:
                    itm.pop("TriggerDetails", None)
                if "HasProtein" in itm:
                    itm.pop("HasProtein", None)

            myfile_write_fname = os.path.join(args.Extr_Sent_Dir, pmid.strip())
            os.remove(myfile_write_fname)
            json_fn = open(myfile_write_fname, 'w')
            json.dump(sentence_data, json_fn, indent=4, ensure_ascii=False)
            json_fn.close()


def subst_invalid_chars(data, invalid_chars):
    return ''.join(c if c not in invalid_chars else "_" for c in data )


def print_categories(sent):
    for c in sent:
        print c, "\t", unicodedata.category(c)


def detect_UnicodeChar(mysent):
    for i, ch in enumerate(category_chars):
        print i, "\t", ch.encode("unicode_escape")
        newsent = subst_unicode_char(mysent, ch)
        print newsent, "\n"


def subst_unicode_char(data, ch):
    return ''.join(c if c != ch else "_" for c in data )

#def strip_control_chars(data):
#    return ''.join(c if unicodedata.category(c) != 'Cc' else "_" for c in data)

def getSentencesFile(Extr_Sent_Dir, pmid):
    sent_file = get_file(Extr_Sent_Dir, pmid)
    file_handle = open(sent_file, "rb")
    sentfile_data = json.load(file_handle)  # sentfile_data is a dictionary
    file_handle.close()
    if "sentenceList" in sentfile_data and "@items" in sentfile_data["sentenceList"] and len(sentfile_data["sentenceList"]["@items"]) > 0:
        return sentfile_data    # There are extracted sentences for this article!
    else:
        return {}


def get_file(file_dir, pattern):
    myfilelist = glob.glob(os.path.join(file_dir, pattern))
    if len(myfilelist) > 1:
        print "Multiple files were matched by glob! \n", myfilelist
        exit()
    return myfilelist[0]


def get_UnicodeCategory(unicode_fileloc, category):
    global category_chars
    category_chars = []
    with open(unicode_fileloc, 'rt') as f1:
        for line in f1:
            line = line.strip()
            fields = line.split(";")
            fields[0] = "\u" + fields[0].strip()
            if fields[2] == category:
                category_chars.append(fields[0].decode('unicode_escape'))
    #print "The no. of chars in the given category: ", category, " is: ", len(category_chars)


if __name__ == "__main__":
    main_body()
