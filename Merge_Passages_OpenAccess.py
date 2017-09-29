import glob;
import argparse;
import sys;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os, errno, ntpath
from operator import itemgetter
from itertools import *
from pprint import pprint;
import cPickle as pickle;


# python Merge_Passages_OpenAccess.py RNAi
# python Merge_Passages_OpenAccess.py RNAi --feedback ./user_feedback_JSON/Apr4_2017/openaccess_passage_dict_WithFeedback.json


def main_body():
    global absId_para_sentId
    parser = argparse.ArgumentParser(prog='Merge_Passages', usage='Merge_Passages_OpenAccess.py <label> [--feedback <Path to feedback passage dict>]', description='Script to group sentences into contiguous passages wrt the given label')
    parser.add_argument('label', help='The field we are trying to predict')
    parser.add_argument('-f', '--feedback', help='Path to the feedback passage dict json file')     # Optional arg

    args = parser.parse_args()
    label = args.label
    absId_para_sentId = pickle.load(open("absId_para_sentId.p", "rb"))

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    if args.feedback is None:
        passage_dict_json = "openaccess_passage_dict.json"
        protdictnames = ["openaccess_posprotdict.p", "openaccess_negprotdict.p"]
        if os.path.isfile(passage_dict_json):
            print "Found - " + passage_dict_json + "! Only recreating the given label."
            file_handle = open(passage_dict_json, "rb")
            passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
            file_handle.close()
        else:
            print "The dictionary - " + passage_dict_json + " was not found in current directory. Creating a new one!"
            passage_dict = {}
    else:
        passage_dict_json = ntpath.basename(args.feedback)
        drname = ntpath.dirname(args.feedback)
        protdictnames = ["openaccess_posprotdict_WithFeedback.p", "openaccess_negprotdict_WithFeedback.p"]
        if os.path.isfile(args.feedback):
            print "Found - " + passage_dict_json + "! Updating the passage details for the given label."
            file_handle = open(args.feedback, "rb")
            passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
            file_handle.close()
        else:
            print "The feedback passage dictionary - " + passage_dict_json + " was not found in given directory. Please check!"
            exit()

    for dictname in protdictnames:
        if dictname.startswith("openaccess_posprotdict"):
            print "Creating Positive instances for label: ", label
            myclass = "Pos"
        else: 
            print "Creating Negative instances for label: ", label
            myclass = "Neg"

        if args.feedback is None:
            myprotdict = pickle.load(open(dictname, "rb"))
        else:
            myprotdict = pickle.load(open(os.path.join(drname, dictname), "rb"))

        for pmid in myprotdict:
            if pmid not in passage_dict:
                passage_dict[pmid] = {}
            if label not in passage_dict[pmid]:
                passage_dict[pmid][label] = {}
            if args.feedback is None:
                passage_dict[pmid][label][myclass] = {}     # I do not want to lose the "HighlightObject" and "Timestamp" keys for all positive, unique PMID-label-protein triples

            for prot in myprotdict[pmid][label]:
                if args.feedback is None:
                    passage_dict[pmid][label][myclass][prot] = {}   # I do not want to lose the "HighlightObject" and "Timestamp" keys for each unique PMID-label-protein triple
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
                    if len(passage_text) == 1 and len(passage_text[0].strip()) == 0:
                        continue    # Ignore "single-phrase" whole passages that contain only whitespace. Mostly these are "Pos" passages
                    tempdict["textOfInterest"] = passage_text       # passage_text is a list of sentences
                    tempdict["tag"] = pmid + "_" + str(psg_ctr) + "--" + foundprot
                    tempdict["id"] = psg_ctr
                    if myclass == "Pos":
                        tempdict["class"] = "1"
                    else:
                        tempdict["class"] = "-1"
                    #csvstr = pmid.strip() + "|_|" + prot + "|_|" + tempdict["class"] + "|_|" + tempdict["tag"] + "|_|" + passage_text
                    #csv_fd.write(csvstr + '\n')
                    tempdict["details"] = []
                    sameparagrpids = get_absIds_samepara(pmid, contg_grp)  # Returns a sorted list of absolute id lists
                    for sameparagrp in sameparagrpids:  # sameparagrp is a list of tuples
                        tempdetdict = {}
                        same_passage_text = get_passagetext(sameparagrp, myprotdict[pmid][label][prot]["AbsoluteIds"])
                        tempdetdict["Text"] = same_passage_text
                        tempdetdict["htmlIdxPath"] = absId_para_sentId[pmid][sameparagrp[0][0]]["htmlIdxPath"]
                        tempdetdict["sentenceDetails"] = []
                        for idx in sameparagrp:     # idx is a tuple
                            tempdetdict["sentenceDetails"].append(myprotdict[pmid][label][prot]["AbsoluteIds"][idx[0]][idx[1]])
                        tempdict["details"].append(tempdetdict)
                    psg_ctr += 1
                    passage_dict[pmid][label][myclass][prot]["passageDetails"].append(tempdict)

    # Check whether each negative instance has a protein reference in it
    for pmid in passage_dict:
        for prot in passage_dict[pmid][label]["Neg"].keys():
            templist = []
            for contg_psg_dict in passage_dict[pmid][label]["Neg"][prot]["passageDetails"]:   # passage_dict[pmid][label]["Neg"][prot]["passageDetails"]
                protreffound = False                                                          # is a list of contiguous passages, i.e. a list of dicts
                for para_dict in contg_psg_dict["details"]:                 # contg_psg_dict["details"] is a list of dicts. Each dict represents text from a single para
                    for sent_dict in para_dict["sentenceDetails"]:          # para_dict["sentenceDetails"] is a list of dicts. Each dict represents details about a sentence
                        if "HasProtein" in sent_dict:
                            templist.append(dict(contg_psg_dict))
                            protreffound = True; break
                    if protreffound:   break
            if len(templist) == 0:
                print "In PMID: ", pmid, "\tRemoving protein: ", prot
                passage_dict[pmid][label]["Neg"].pop(prot)
            else:
                passage_dict[pmid][label]["Neg"][prot]["passageDetails"] = templist

    if args.feedback == None:
        json_fn = open(passage_dict_json, 'w')
        json.dump(passage_dict, json_fn, indent=4, ensure_ascii=False)
        json_fn.close()
    else:
        json_fn = open(args.feedback, 'w')
        json.dump(passage_dict, json_fn, indent=4, ensure_ascii=False)
        json_fn.close()


def get_passagetext(absids, mydict):
    absids.sort(key=itemgetter(0, 1))       # absids = sorted(absids, key=itemgetter(0, 1))   absids is a list of tuples
    mytext = []
    for idx in absids:         # idx is a tuple
        if "partialSentence" in mydict[idx[0]][idx[1]]:
            mytext.append(mydict[idx[0]][idx[1]]["partialcleansedText"])
        else:
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
