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

# python Create_HighlightObjects.py RNAi ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0 ./downloaded_articles

def main_body():
    global passage_dict; global spanctr; global currtime
    parser = argparse.ArgumentParser(prog='Create_HighlightObjects.py', usage='Create_HighlightObjects.py <label> <psg_dict_loc> <downloaded_articles_dir>', description='Script to create labelled ML file')
    parser.add_argument('label', help='The class (it may be a extra test or a extra type) for current data file')
    parser.add_argument('psg_dict_loc', help='The location where the OpenAccess passage_dict json file is present')
    parser.add_argument('downloaded_articles_dir', help='Path where openaccess html articles were downloaded')
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
    #currtime = int(round(time.time()))
    currtime = 1478131000   # The timestamp represents this date:- 11-02-2016 18:56:40

    # Whether 2 contiguous passages end and start, respectively, in the same html <p> / <h3> node
    idxPath_ids = {}
    for pmid in passage_dict:
        idxPath_ids[pmid] = {}
        for prot in passage_dict[pmid][label]["Pos"]:   # prot represents a key of the dict passage_dict[pmid][label]["Pos"]
            idxPath_ids[pmid][prot] = {}
            for contg_psg_dict in passage_dict[pmid][label]["Pos"][prot]["passageDetails"]:   # contg_psg_dict is a li elem (which is a dict) of the list - passage_dict[pmid][label]["Pos"][prot]["passageDetails"]
                for i, para_dict in enumerate(contg_psg_dict["details"]):     # para_dict is a li elem (which is a dict) of the list - contg_psg_dict["details"]
                    myent = str(contg_psg_dict["id"]) + "_" + str(i)
                    if para_dict["htmlIdxPath"] not in idxPath_ids[pmid][prot]:
                        idxPath_ids[pmid][prot][para_dict["htmlIdxPath"]] = [myent]
                    else:
                        idxPath_ids[pmid][prot][para_dict["htmlIdxPath"]].append(myent)

    # The entire dict idxPath_ids corresponds to the user-supplied label
    for pmid in idxPath_ids:
        print "PMID: ", pmid, "\tPMCID: ", pm_pmc[pmid], "\tLabel: ", label
        soup = getBeautifulSoup(args.downloaded_articles_dir, pm_pmc[pmid])
        for prot in idxPath_ids[pmid]:
            highlightlist = []; spanctr = -1
            for nodePath in idxPath_ids[pmid][prot]:
                sentenceList = copy.deepcopy(getSentenceList(pmid, label, prot, nodePath, idxPath_ids[pmid][prot][nodePath]))
                copy_soup = copy.copy(soup)
                targetnode = highlightnode(pmid, copy_soup, nodePath, sentenceList)
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


def highlightnode(pmid, soup, nodePath, sentenceList):
    mynode = getnodefromidxpath(nodePath, soup.body)  # target node that needs to be highlighted
    nextOffset = 0; mylist = []
    for child in mynode.descendants:
        if isinstance(child, NavigableString):
            if len(sentenceList) == 0:
                break
            childstr = child
            currOffset = nextOffset
            nextOffset = currOffset + len(childstr)
            child_parts = partitionChildStr(sentenceList, currOffset, childstr)
            if len(child_parts) == 0:
                continue
            child_parts = mergeChildParts(child_parts, childstr)
            curr_parent = child.parent
            child_idx = getchildindex(curr_parent, child)
            parent_path = getchildpathindex(curr_parent, soup.body)
            mylist.append([parent_path, curr_parent, child, child_parts, child_idx])

    if len(sentenceList) > 0:
        print "Not all sentences were matched in the node: ", nodePath; pprint(sentenceList); exit()

    mydict = {}
    for elem in mylist:
        if elem[0] in mydict:
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
            for i, (strg, switch) in enumerate(parts):
                navg_strg = NavigableString(strg)
                if switch:
                    span_tag = soup.new_tag("span")
                    span_tag = populatespanattr(pmid, span_tag)
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


def partitionChildStr(sentenceList, childOffset, childstr):
    parts = []; nextChildOffset = childOffset + len(childstr)
    while len(childstr) > 0:
        if len(sentenceList) == 0:
            break
        if sentenceList[0]["charOffset"] >= childOffset and sentenceList[0]["charOffset"] < nextChildOffset:
            #print "Sentence Text: ", sentenceList[0]["sentenceText"].encode("utf-8"); print "Child String: ", childstr.encode("utf-8"); print parts;
            if sentenceList[0]["sentenceLen"] > (nextChildOffset - sentenceList[0]["charOffset"]):
                mymatch = childstr[(sentenceList[0]["charOffset"] - childOffset):]
                parts.extend([(childstr[0:(sentenceList[0]["charOffset"] - childOffset)], 0), (mymatch, 1)])
                sentenceList[0]["sentenceText"] = sentenceList[0]["sentenceText"][(nextChildOffset - sentenceList[0]["charOffset"]):]
                sentenceList[0]["charOffset"] = nextChildOffset
                sentenceList[0]["sentenceLen"] = len(sentenceList[0]["sentenceText"])
                childstr = ""
            elif sentenceList[0]["sentenceLen"] < (nextChildOffset - sentenceList[0]["charOffset"]):
                    mymatch = childstr[(sentenceList[0]["charOffset"] - childOffset):(sentenceList[0]["charOffset"] - childOffset)+sentenceList[0]["sentenceLen"]]
                    parts.extend([(childstr[0:(sentenceList[0]["charOffset"] - childOffset)], 0), (mymatch, 1)])
                    childstr = childstr[(sentenceList[0]["charOffset"] - childOffset)+sentenceList[0]["sentenceLen"]:]
                    childOffset = sentenceList[0]["charOffset"] + sentenceList[0]["sentenceLen"]
                    #print "New Child Str: ", childstr.encode("utf-8"), "New child offset: ", childOffset
                    sentenceList.pop(0)
                    #if len(sentenceList) > 0:  print "New! New! New! Sentence Text: ", sentenceList[0]["sentenceText"].encode("utf-8")
            else:
                mymatch = childstr[(sentenceList[0]["charOffset"] - childOffset):]
                parts.extend([(childstr[0:(sentenceList[0]["charOffset"] - childOffset)], 0), (mymatch, 1)])
                sentenceList.pop(0)
                childstr = ""
        else:
            break
    if len(childstr) > 0 and len(parts) > 0:
        parts.append((childstr, 0))
    return parts


def mergeChildParts(child_parts, childstr):
    merged_parts = []; origstr = childstr
    while len(childstr) > 0:
        if len(child_parts) == 0:
            break
        strLocate = child_parts[0][0]
        if strLocate == "":
            child_parts.pop(0); continue
        strIdx = childstr.find(strLocate)
        if strIdx == -1:
            print "Could not locate str, index is -1.\t", origstr.encode("utf-8"), "\t", strLocate.encode("utf-8"); exit()
        if strIdx > 0:
            extraStr = childstr[0:strIdx]
            if extraStr.isspace():
                merged_parts.append((extraStr+strLocate, child_parts[0][1]))
                childstr = childstr[strIdx+len(strLocate):]
                child_parts.pop(0)
            else:
                print "Non white space characters at the beginning of childstr.\t", origstr.encode("utf-8"), "\t", strLocate.encode("utf-8"); exit()
        elif strIdx == 0:
            merged_parts.append(child_parts[0])
            childstr = childstr[len(strLocate):]
            child_parts.pop(0)
    if len(childstr) > 0:
        if childstr.isspace():
            lastpart = merged_parts[-1]
            merged_parts[-1] = (lastpart[0]+childstr,lastpart[1])
        else:
            print "Non white space characters at the end of childstr.\t", origstr.encode("utf-8"), "\t", childstr.encode("utf-8"); exit()

    # At this point the entire childstr has been partitioned!
    prevlen = 0
    while prevlen != len(merged_parts):
        prevlen = len(merged_parts)
        if prevlen == 1: break
        new_merged_parts = []; skip_idx = []
        for i,elem in enumerate(merged_parts):
            if i in skip_idx:   continue;
            if i+1 < prevlen and elem[0].isspace() and elem[1] == 0:
                new_merged_parts.append((elem[0] + merged_parts[i + 1][0], merged_parts[i + 1][1]))
                skip_idx.append(i + 1); continue
            if i+1 < prevlen and elem[1] == merged_parts[i+1][1]:
                new_merged_parts.append((elem[0] + merged_parts[i + 1][0], merged_parts[i + 1][1]))
                skip_idx.append(i + 1); continue
            if i+1 < prevlen and merged_parts[i+1][0].isspace() and merged_parts[i+1][1] == 0:
                new_merged_parts.append((elem[0] + merged_parts[i + 1][0], elem[1]))
                skip_idx.append(i + 1); continue
            new_merged_parts.append(elem)
        merged_parts = list(new_merged_parts)

    return merged_parts
        

def getSentenceList(pmid, label, prot, nodePath, mylist):
    sentenceList = []
    for elem in mylist:
        contg_psg_dict_id = int(elem.split("_")[0])
        para_dict_id = int(elem.split("_")[1])
        para_dict = passage_dict[pmid][label]["Pos"][prot]["passageDetails"][contg_psg_dict_id]["details"][para_dict_id]
        if para_dict["htmlIdxPath"] != nodePath:
            print "Something has gone wrong!"
            exit()
        sentenceList.extend(para_dict["sentenceDetails"])
    newlist = sorted(sentenceList, key=itemgetter("charOffset"))
    return newlist


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


def populatespanattr(pmid, span_tag):
    span_tag["id"] = getspanid(pmid)
    span_tag["data-highlighted"] = "true"
    span_tag["data-timestamp"] = currtime
    span_tag["class"] = "highlighted"
    span_tag["style"] = "background-color: rgb(255, 255, 123);"
    return span_tag


if __name__ == "__main__":
    main_body()
