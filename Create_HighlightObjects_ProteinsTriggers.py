# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
import sys
import glob
import json
reload(sys)
sys.setdefaultencoding("utf-8")
import os
import copy
from bs4 import BeautifulSoup
from bs4 import NavigableString
from bs4 import Tag
from pprint import pprint
import cPickle as pickle
from operator import itemgetter
import time


# python Create_HighlightObjects_ProteinsTriggers.py RNAi ./datafiles_TriggerHighlights_BaseLine1 ./downloaded_articles ./JSON_SENTENCE_ASSOCS_OpenAccess
# python Create_HighlightObjects_ProteinsTriggers.py omission ./datafiles_TriggerHighlights_BaseLine1 ./downloaded_articles ./JSON_SENTENCE_ASSOCS_OpenAccess
# python Create_HighlightObjects_ProteinsTriggers.py KO ./datafiles_TriggerHighlights_BaseLine1 ./downloaded_articles ./JSON_SENTENCE_ASSOCS_OpenAccess


def main_body():
    global passage_dict; global spanctr; global currtime
    parser = argparse.ArgumentParser(prog='Create_HighlightObjects_ProteinsTriggers.py', usage='Create_HighlightObjects_ProteinsTriggers.py <label> <psg_dict_loc> <downloaded_articles_dir> <rootDir>', description='Script to create highlight objects containing only proteins and trigger words')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('psg_dict_loc', help='The location where the OpenAccess passage_dict json file is present')
    parser.add_argument('downloaded_articles_dir', help='Path where openaccess html articles were downloaded')
    parser.add_argument('rootDir', help='Directory where the files containing individual sentences are located')
    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    passagedictfile = "openaccess_passage_dict.json"
    fullpath_passagedictfile = os.path.join(args.psg_dict_loc, passagedictfile)
    if os.path.isfile(fullpath_passagedictfile):
        print "Found the dictionary - " + passagedictfile + " Everything is perfect in the world!"
        file_handle = open(fullpath_passagedictfile, "rb")
        passage_dict = json.load(file_handle)  # mymatchfile_data is a list of dictionaries
        file_handle.close()
    else:
        print "Couldn't locate the dictionary - " + passagedictfile + "! Please place it in - ", args.psg_dict_loc
        exit()

    label = args.label
    pm_pmc = pickle.load(open("pmid_pmcid_Extras.p", "rb"))
    #prot_uniprotID = pickle.load(open("dictionary_nameToUniProtId.p", "rb"))
    #currtime = int(round(time.time()))
    currtime = 1478131000   # The timestamp represents this date:- 11-02-2016 18:56:40

    # Whether 2 contiguous passages end and start, respectively, in the same html <p> / <h3> node
    idxPath_ids = {}
    for pmid in passage_dict:
        for prot in passage_dict[pmid][label]["Pos"]:  # prot represents a protein that occurs in the DatumKB for that article and label
            #print passage_dict[pmid][label]["Pos"].keys()
            if pmid not in idxPath_ids:
                idxPath_ids[pmid] = {}
            idxPath_ids[pmid][prot] = {}
            for contg_psg_dict in passage_dict[pmid][label]["Pos"][prot]["passageDetails"]:  # contg_psg_dict is a li elem (which is a dict) of the list - passage_dict[pmid][label]["Pos"][prot]["passageDetails"]
                for i, para_dict in enumerate(contg_psg_dict["details"]):  # para_dict is a li elem (which is a dict) of the list - contg_psg_dict["details"]
                    paraprotOffsets = []
                    for sentence in para_dict["sentenceDetails"]:
                        if "ProteinDetails" in sentence:
                            for protein in sentence["ProteinDetails"]:
                                if protein["protName"].strip() in passage_dict[pmid][label]["Pos"][prot]["MatchedProtein"]:
                                    paraprotOffsets.extend(protein["protcharOffset"])
                    if para_dict["htmlIdxPath"] not in idxPath_ids[pmid][prot]:
                        idxPath_ids[pmid][prot][para_dict["htmlIdxPath"]] = paraprotOffsets
                    else:
                        idxPath_ids[pmid][prot][para_dict["htmlIdxPath"]].extend(paraprotOffsets)   # The same sentence can be evidence for 2 different pos proteins for the same article

    # Populate idxPath_ids with the details of paragraphs that contain the trigger words
    for pmid in idxPath_ids:
        htmlPath_Offsets = {}
        if not pmid.strip(): continue
        fname_string1 = '%s*' % (os.path.join(args.rootDir, pmid.strip()))
        fname = glob.glob(fname_string1)
        if len(fname) == 0:
            print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
            continue
        myfile = fname[0].strip()

        file_handle = open(myfile, "rb")
        myfile_data = json.load(file_handle)  # myfile_data is a list of dictionaries
        file_handle.close()

        if len(myfile_data["paragraphList"]) == 0: continue;
        for para in myfile_data["paragraphList"]:  # para is a dictionary
            if "sentenceList" not in para or "@items" not in para["sentenceList"] or len(para["sentenceList"]["@items"]) == 0:
                continue
            htmlpath = para["htmlIdxPath"]
            for itm in para["sentenceList"]["@items"]:  # itm represents a sentence in the article
                if "TriggerDetails" in itm:
                    for trigword in itm["TriggerDetails"]:  # A sentence can have multiple trigger words
                        if trigword["label"] == label:
                            if htmlpath in htmlPath_Offsets:
                                htmlPath_Offsets[htmlpath].extend(trigword["triggerCharOffset"])    # htmlPath_Offsets[htmlpath] is a list of lists
                            else:
                                htmlPath_Offsets[htmlpath] = copy.deepcopy(trigword["triggerCharOffset"])
        for prot in idxPath_ids[pmid]:
            for path in htmlPath_Offsets:
                if path in idxPath_ids[pmid][prot]:
                    idxPath_ids[pmid][prot][path].extend(htmlPath_Offsets[path])
                else:
                    idxPath_ids[pmid][prot][path] = copy.deepcopy(htmlPath_Offsets[path])

    for pmid in idxPath_ids:
        for prot in idxPath_ids[pmid]:
            for nodePath in idxPath_ids[pmid][prot]:
                idxPath_ids[pmid][prot][nodePath] = copy.deepcopy(sorted(idxPath_ids[pmid][prot][nodePath], key=itemgetter(0)))
    #pprint(idxPath_ids); exit()

    for pmid in idxPath_ids:
        #pmid = "21573184"
        print "PMID: ", pmid, "\tPMCID: ", pm_pmc[pmid], "\tLabel: ", label
        soup = getBeautifulSoup(args.downloaded_articles_dir, pm_pmc[pmid])
        for prot in idxPath_ids[pmid]:
            highlightlist = []; spanctr = -1
            for nodePath in idxPath_ids[pmid][prot]:
                copy_soup = copy.copy(soup)
                targetnode = highlightnode(pmid, copy_soup, nodePath, copy.deepcopy(idxPath_ids[pmid][prot][nodePath]))
                for child in targetnode.descendants:
                    if isinstance(child, Tag) and child.name == "span" and child.has_attr('id') and child["id"].startswith(pmid):
                        childlen = len(child.string)
                        childindexliststr = getchildpathindex(child, copy_soup.body)
                        if isinstance(child.previous_sibling, NavigableString):
                            offset = len(child.previous_sibling)
                        else:
                            offset = 0
                        span_tag_str = "<span data-highlighted=\"" + child["data-highlighted"] + "\" data-timestamp=\"" + str(child["data-timestamp"]) + "\" class=\"" + child["class"] + "\" style=\"" + child["style"] + "\"></span>"
                        templist = [span_tag_str, child.string, childindexliststr, str(offset), str(childlen)]
                        highlightlist.append(templist)
                #pprint(highlightlist); exit()
            passage_dict[pmid][label]["Pos"][prot]["HighlightObject"] = highlightlist
            passage_dict[pmid][label]["Pos"][prot]["Timestamp"] = currtime

    #pickle.dump(passage_dict, open("openaccess_passage_dict.p", "wb"))
    json_fn = open(fullpath_passagedictfile, 'w')
    json.dump(passage_dict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


def highlightnode(pmid, soup, nodePath, offsetList):
    mynode = getnodefromidxpath(nodePath, soup.body)  # target node that needs to be highlighted
    #print "Node Path: ", nodePath
    nextOffset = 0; mylist = []
    for child in mynode.descendants:
        if isinstance(child, NavigableString):
            if len(offsetList) == 0:
                break
            childstr = child
            currOffset = nextOffset
            nextOffset = currOffset + len(childstr)
            child_parts = partitionChildStr(offsetList, currOffset, childstr)
            if len(child_parts) == 0:
                continue
            curr_parent = child.parent
            child_idx = getchildindex(curr_parent, child)
            parent_path = getchildpathindex(curr_parent, soup.body)
            mylist.append([parent_path, curr_parent, child, child_parts, child_idx])

    if len(offsetList) > 0:
        print "Not all offsets were matched in the node: ", nodePath; pprint(offsetList); 
        print mynode
        exit()

    mydict = {}
    for elem in mylist:
        if elem[0] in mydict:   # The same parent can have multiple children which are navigable strings
            mydict[elem[0]].append([elem[1], elem[2], elem[3], elem[4]])  # [curr_parent, child, found_psg, child_idx]
        else:
            mydict[elem[0]] = [[elem[1], elem[2], elem[3], elem[4]]]

    for delem in mydict:
        mydict[delem].sort(key=lambda x: int(x[3]), reverse=True)

    for ppath in mydict:
        for lelem in mydict[ppath]:  # mydict[ppath] is a list of lists. lelem is a 4-element list: [curr_parent 0, child 1, child_parts 2, child_idx 3]
            childstr = lelem[1]
            parts = lelem[2]   # parts is a list of tuples
            child_idx = lelem[3]
            for i, (strg, switch, color) in enumerate(parts):
                navg_strg = NavigableString(strg)
                if switch:
                    span_tag = soup.new_tag("span")
                    span_tag = populatespanattr(pmid, span_tag, color)
                    span_tag.append(navg_strg)
                    if i == 0:
                        childstr.replace_with(span_tag)
                    else:
                        lelem[0].insert(child_idx, span_tag)
                else:
                    if i == 0:
                        childstr.replace_with(navg_strg)
                    else:
                        lelem[0].insert(child_idx, navg_strg)
                child_idx += 1
    return mynode


def partitionChildStr(offsetList, childOffset, childstr):
    parts = []; nextChildOffset = childOffset + len(childstr)
    while len(childstr) > 0:
        #print offsetList, "\n", parts
        if len(offsetList) == 0:
            break
        if offsetList[0][0] >= childOffset and offsetList[0][0] < nextChildOffset:
            #print "Child String: ", childstr.encode("utf-8"); print parts;
            if offsetList[0][1] > nextChildOffset:
                mymatch = childstr[(offsetList[0][0] - childOffset):]
                parts.extend([(childstr[0:(offsetList[0][0] - childOffset)], 0, "NULL"), (mymatch, 1, offsetList[0][2])])
                offsetList[0][0] = nextChildOffset
                childstr = ""
            elif offsetList[0][1] < nextChildOffset:
                    mymatch = childstr[(offsetList[0][0] - childOffset):(offsetList[0][1] - childOffset)]
                    parts.extend([(childstr[0:(offsetList[0][0] - childOffset)], 0, "NULL"), (mymatch, 1, offsetList[0][2])])
                    childstr = childstr[(offsetList[0][1] - childOffset):]
                    childOffset = offsetList[0][1]
                    #print "New Child Str: ", childstr.encode("utf-8"), "New child offset: ", childOffset
                    offsetList.pop(0)
                    #if len(sentenceList) > 0:  print "New! New! New! Sentence Text: ", sentenceList[0]["sentenceText"].encode("utf-8")
            else:
                mymatch = childstr[(offsetList[0][0] - childOffset):]
                parts.extend([(childstr[0:(offsetList[0][0] - childOffset)], 0, "NULL"), (mymatch, 1, offsetList[0][2])])
                offsetList.pop(0)
                childstr = ""
        else:
            break
    if len(childstr) > 0 and len(parts) > 0:
        parts.append((childstr, 0, "NULL"))

    newparts = []
    for elem in parts:
        if len(elem[0]) > 0:
            newparts.append(elem)
    return newparts


def getnodefromidxpath(idxpath, referrence):
    nodesidx = idxpath.split(':')
    prevnode = referrence
    while len(nodesidx) > 0:
        currnodeidx = int(nodesidx.pop(0))
        for i, child in enumerate(prevnode.children):
            if i == currnodeidx:
                prevnode = child
                break
    return prevnode


def getBeautifulSoup(download_dir, pmcid):
    fname_string = '%s' % (os.path.join(download_dir, pmcid.strip(), pmcid.strip() + ".html"))
    fname = glob.glob(fname_string)
    if len(fname) == 0:
        print "The HTML file was not found for PMCID: ", pmcid.strip()
        exit()
    myfile = fname[0].strip()
    soup = BeautifulSoup(open(myfile), "lxml")
    return soup


def getchildpathindex(child, referrence):   # child cannot be the same node as referrence
    childindexlist = []
    currchild = child
    for parent in child.parents:
        for j, c in enumerate(parent.children):
            if c is currchild:
                childindexlist.insert(0, str(j))
                break
        currchild = parent
        if parent is referrence:
            break
    return ":".join(childindexlist)


def getchildindex(parent, child):
    for i, c in enumerate(parent.children):
        if c is child:
            return i


def getspanid(pmid):
    global spanctr
    spanctr += 1
    #print spanctr
    return pmid + "_" + str(spanctr)


def populatespanattr(pmid, span_tag, bgcolor):
    span_tag["id"] = getspanid(pmid)
    span_tag["data-highlighted"] = "true"
    span_tag["data-timestamp"] = currtime
    span_tag["class"] = "highlighted"
    span_tag["style"] = "background-color: " + bgcolor
    return span_tag


if __name__ == "__main__":
    main_body()
