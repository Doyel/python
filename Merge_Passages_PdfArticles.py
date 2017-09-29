import glob;
import argparse;
import sys;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno;
from operator import itemgetter
from itertools import *
from pprint import pprint;
import cPickle as pickle;


# python Merge_Passages_PdfArticles.py RNAi
# python Merge_Passages_PdfArticles.py RNAi -t


def main_body():
    parser = argparse.ArgumentParser(prog='Merge_Passages', usage='Merge_Passages_PdfArticles.py <label> [-t]', description='Script to group sentences into contiguous passages wrt the given label')
    parser.add_argument('label', help='The field we are trying to predict')
    #parser.add_argument('csvLoc', help='The directory in which to store the csv files')
    parser.add_argument('-t', action='store_true')

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    label = args.label

    if args.t:
        if os.path.isfile("test_passage_dict.json"):
            print "Found - test_passage_dict.json! Only recreating the given label."
            file_handle = open("test_passage_dict.json", "rb")
            passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
            file_handle.close()
        else:
            print "The dictionary - test_passage_dict.json was not found in current directory. Creating a new one!"
            passage_dict = {}
        protdictlist = ["test_posprotdict.p", "test_negprotdict.p"]
        passage_dict_json = "test_passage_dict.json"
        #csv_file_name = label + "_test.csv"
    else:
        if os.path.isfile("train_passage_dict.json"):
            print "Found - train_passage_dict.json! Only recreating the given label."
            file_handle = open("train_passage_dict.json", "rb")
            passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
            file_handle.close()
        else:
            print "The dictionary - train_passage_dict.json was not found in current directory. Creating a new one!"
            passage_dict = {}
        protdictlist = ["train_posprotdict.p", "train_negprotdict.p"]
        passage_dict_json = "train_passage_dict.json"
        #csv_file_name = label + "_train.csv"

    #csv_full_path = os.path.join(args.csvLoc, csv_file_name)
    #silentremove(csv_full_path)
    #csv_fd = open(csv_full_path, 'w')

    for dictname in protdictlist:
        if dictname.endswith("posprotdict.p"):
            print "Creating Positive instances for label: ", label
            myclass = "Pos"
        else: 
            print "Creating Negative instances for label: ", label
            myclass = "Neg"

        myprotdict = pickle.load(open(dictname, "rb"))

        for pmid in myprotdict:
            if pmid not in passage_dict:
                passage_dict[pmid] = {}
            if label not in passage_dict[pmid]:
                passage_dict[pmid][label] = {}   
            passage_dict[pmid][label][myclass] = {}

            for prot in myprotdict[pmid][label]:                    
                passage_dict[pmid][label][myclass][prot] = {}
                passage_dict[pmid][label][myclass][prot]["Mentions"] = list(myprotdict[pmid][label][prot]["Mentions"])
                passage_dict[pmid][label][myclass][prot]["MatchedProtein"] = list(myprotdict[pmid][label][prot]["MatchedProtein"])
                passage_dict[pmid][label][myclass][prot]["passageDetails"] = []
                absIds = myprotdict[pmid][label][prot]["AbsoluteIds"].keys()
                grouped_absIds = group_absoluteIds(absIds)    # Returns a sorted list of contiguous group of absolute id lists
                grouped_partial_absIds = group_partial_absoluteIds(grouped_absIds, myprotdict[pmid][label][prot]["AbsoluteIds"])  # Returns a list of list of tuples
                psg_ctr = 0; foundprot = prot.replace(" ", "__")
                for contg_grp in grouped_partial_absIds:    # contg_grp is a list of tuples
                    tempdict = {}
                    passage_text = get_passagetext(contg_grp, myprotdict[pmid][label][prot]["AbsoluteIds"])
                    tempdict["textOfInterest"] = passage_text       # passage_text is a list of sentences
                    tempdict["tag"] = pmid + "_" + str(psg_ctr) + "--" + foundprot
                    tempdict["id"] = psg_ctr
                    if myclass == "Pos":
                        tempdict["class"] = "1"
                    else:
                        tempdict["class"] = "-1"
                    #csvstr = pmid.strip() + "|_|" + prot + "|_|" + tempdict["class"] + "|_|" + tempdict["tag"] + "|_|" + passage_text
                    #csv_fd.write(csvstr + '\n')
                    tempdict["sentenceDetails"] = []                    
                    for idx in contg_grp:     # idx is a tuple
                        tempdict["sentenceDetails"].append(myprotdict[pmid][label][prot]["AbsoluteIds"][idx[0]][idx[1]])
                    psg_ctr += 1
                    passage_dict[pmid][label][myclass][prot]["passageDetails"].append(tempdict)
    #csv_fd.close()

    # Check whether each negative instance has a protein reference in it
    for pmid in passage_dict:
        for prot in passage_dict[pmid][label]["Neg"].keys():
            templist = []
            for contg_psg_dict in passage_dict[pmid][label]["Neg"][prot]["passageDetails"]:   # passage_dict[pmid][label]["Neg"][prot]["passageDetails"] is a list of contiguous passages, i.e. a list of dicts
                for sent_dict in contg_psg_dict["sentenceDetails"]:                 # contg_psg_dict["details"] is a list of dicts. Each dict represents text from a single para                    
                    if "HasProtein" in sent_dict:
                        templist.append(dict(contg_psg_dict))
                        break
                    
            if len(templist) == 0:
                print "In PMID: ", pmid, "\tRemoving protein: ", prot
                passage_dict[pmid][label]["Neg"].pop(prot)
            else:
                passage_dict[pmid][label]["Neg"][prot]["passageDetails"] = templist

    #pickle.dump(passage_dict, open(passage_dict_file, "wb"))
    json_fn = open(passage_dict_json, 'w')
    json.dump(passage_dict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


def get_passagetext(absids, mydict):
    absids.sort(key=itemgetter(0, 1))       # absids = sorted(absids, key=itemgetter(0, 1))   absids is a list of tuples
    mytext = []
    for idx in absids:
        mytext.append(mydict[idx[0]][idx[1]]["cleansedText"])

    return mytext


def get_absIds_samepara(pmid, absids):
    absids.sort(key=itemgetter(0, 1))
    mylist = []; prev_paraId = -1
    for i,idx in enumerate(absids):     # idx is a tuple
        if prev_paraId != absId_para_sentId[pmid][idx[0]]["paraId"]:
            if i > 0:
                mylist.append(templist)
            templist = [idx]        # Control will always come here when i = 0
        else:
            templist.append(idx)
        if i == len(absids) - 1:
            mylist.append(templist)
        prev_paraId = absId_para_sentId[pmid][idx[0]]["paraId"]

    return mylist


def group_absoluteIds(absids):
    absids.sort(); mylist = []
    for k, g in groupby(enumerate(absids), lambda (i, x): i - x):
        mylist.append(map(itemgetter(1), g))
    return mylist


def group_partial_absoluteIds(grouped_absIds, absids_dict):
    mylist = []
    for idlist in grouped_absIds:
        for j, myid in enumerate(idlist):
            for i, sent in enumerate(absids_dict[myid]):
                if j == 0 and i == 0:
                    templist = [(myid, i)]
                else:
                    if templist[-1][0] == myid:
                        mylist.append(templist)
                        templist = [(myid, i)]
                    else:
                        prevsent = absids_dict[templist[-1][0]][templist[-1][1]]
                        if "partialSentence" in prevsent and "partialSentence" in sent:
                            if prevsent["charOffset"] + prevsent["sentenceLen"] == prevsent["partialOffset"] + prevsent["partialLength"] and sent["partialOffset"] == sent["charOffset"]:
                                templist.append((myid, i))
                            else:
                                mylist.append(templist)
                                templist = [(myid, i)]
                        elif "partialSentence" in prevsent and "partialSentence" not in sent:
                            if prevsent["charOffset"] + prevsent["sentenceLen"] == prevsent["partialOffset"] + prevsent["partialLength"]:
                                templist.append((myid, i))
                            else:
                                mylist.append(templist)
                                templist = [(myid, i)]
                        elif "partialSentence" not in prevsent and "partialSentence" in sent:
                            if sent["partialOffset"] == sent["charOffset"]:
                                templist.append((myid, i))
                            else:
                                mylist.append(templist)
                                templist = [(myid, i)]
                        elif "partialSentence" not in prevsent and "partialSentence" not in sent:
                            templist.append((myid, i))
        mylist.append(templist)

    return mylist


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occured


if __name__ == "__main__":
    main_body()
