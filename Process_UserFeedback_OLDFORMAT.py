# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
import sys
import json
reload(sys)
sys.setdefaultencoding("utf-8")
import os
import re
import copy
import glob
import ntpath
from pprint import pprint
import cPickle as pickle
from bs4 import BeautifulSoup
from bs4 import NavigableString
from operator import itemgetter
from bs4 import Tag
import string
from New_Protein_Detector_OpenAccess import find_proteins, isstopword, fst_init, check_for_unicode


# python Process_UserFeedback_OLDFORMAT.py RNAi ./user_feedback_JSON/Apr4_2017/SavedHighlights_WithFeedback.json ./downloaded_articles ./JSON_SENTENCE_ASSOCS_OpenAccess ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0


def main_body():
    global spanctr
    parser = argparse.ArgumentParser(prog='Process_UserFeedback', usage='Process_UserFeedback.py <label> <feedback_file> <downloaded_articles_dir> <extracted_sentdir>', description='Script to process user feedback on a given label')
    parser.add_argument('label', help='The field we are trying to process feedback')
    parser.add_argument('hilightfile', help='Path to file that contains user feedback')
    parser.add_argument('downloaded_articles_dir', help='Path where openaccess html articles were downloaded')
    parser.add_argument('extracted_sentdir', help='Path where the extracted sentences files for the openaccess html articles are located')
    parser.add_argument('prot_dict_loc', help='The location where the posprot and negprot dict files are present')

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    drname = ntpath.dirname(args.hilightfile)
    file_handle = open(args.hilightfile, "rb")
    feedback_data = json.load(file_handle)
    file_handle.close()

    label = args.label
    pm_pmc = pickle.load(open("pmid_pmcid_Extras.p", "rb"))
    posprotdict = pickle.load(open(os.path.join(args.prot_dict_loc, "openaccess_posprotdict.p"), "rb"))
    negprotdict = pickle.load(open(os.path.join(args.prot_dict_loc, "openaccess_negprotdict.p"), "rb"))
    htmlIdx_paraId = pickle.load(open("htmlIdx_paraId.p", "rb"))
    fst_init("Minimized_Protein_Names_FST.fst", "SymbolTable_Protein_Names_FST.sym")

    # file_handle = open("Processed_Timestamp.json", "rb")
    # file_data = json.load(file_handle)
    # lasttimestamp = file_data["processedTimestamp"]
    # currtimestamp = int(round(time.time()))
    # file_data["processedTimestamp"] = currtimestamp
    # file_data["dateTime"] = datetime.datetime.fromtimestamp(currtimestamp).strftime('%Y-%m-%d %H:%M:%S')
    # file_handle.close()

    # We will process ALL user feedback that is provided to us in SavedHighlights_WithFeedback.json
    # This file is basically a dict that will contain every user feedback that was ever given, since the beginning of time!
    # ["20190815", "24106086", "21629263",      "22833096", "16009723", "18583988",     "14517278", "19783983", "18411307",     "20890305", "18504304", "19234442",     "22245064", "20141835", "21573184"]
    span_dict = {}; ignore_list = []        #["20547488", "22433566", "20026654"]
    print "Creating a deeply-nested dictionary containing the span details"
    for pmid in feedback_data:
        if label not in feedback_data[pmid]:    continue    # We process all of feedback data w.r.t the given label
        if pmid.strip() in ignore_list: continue
        print "PMID: ", pmid, "\tPMCID: ", pm_pmc[pmid], "\tLabel: ", label
        for prot in feedback_data[pmid][label]:
            #if pmid == "20026654" and prot in ["mdc1"]:    continue
            for span in feedback_data[pmid][label][prot]["HighlightObject"]:    # "HighlightObject" is a list of lists. 'span' is a list with 5 elements
                idxpathlist = [int(n) for n in span[2].split(":")]
                spanidx = idxpathlist[-1]; spandetails = {}
                populatespandict(span, spanidx, spandetails)

                if pmid not in span_dict:
                    span_dict[pmid] = {}
                if prot not in span_dict[pmid]:
                    span_dict[pmid][prot] = {}

                tempdict = span_dict[pmid][prot]
                for i, idx in enumerate(idxpathlist):
                    if i != len(idxpathlist)-1:
                        if idx not in tempdict:
                            tempdict[idx] = {}
                        tempdict = tempdict[idx]
                    else:
                        if idx in tempdict:
                            print "idx already present!"; exit()
                        tempdict[idx] = [spandetails]

    # Flatten the span_dict so that we have an "ordered" list of span elements
    span_list = {}
    for pmid in span_dict:
        span_list[pmid] = {}
        for prot in span_dict[pmid]:
            mylist = []
            flattenSpanDict(span_dict[pmid][prot], mylist, "")
            span_list[pmid][prot] = mylist

    print "Re-creating highlights within p and h3 nodes"
    para_childspan_details = {}; contg_highlights_dict = {}
    for pmid in span_list:
        soup = getBeautifulSoup(args.downloaded_articles_dir, pm_pmc[pmid])
        para_childspan_details[pmid] = {}; contg_highlights_dict[pmid] = {}
        article_data = getExtractedSentences(pmid, args.extracted_sentdir)
        print "PMID: ", pmid, "\tPMCID: ", pm_pmc[pmid]
        for prot in span_list[pmid]:
            copy_soup = copy.copy(soup); spanctr = -1
            para_childspan_details[pmid][prot] = {}; contg_highlights_dict[pmid][prot] = {}
            for span in span_list[pmid][prot]:      # span_list[pmid][prot] is a list of span dicts
                highlightnode(pmid, prot, copy_soup, span, para_childspan_details[pmid][prot])

            populate_negprotdict(negprotdict[pmid][label], pmid, prot, posprotdict[pmid][label][prot])
            posprotdict[pmid][label][prot] = {}  # The interface returns the WHOLE set of evidence sentences for any selected combination of - 'prot' and 'label'
            negprotdict[pmid][label][prot] = {}  # 'prot' wouldn't exist in negprotdict beforehand as it is a positive protein
            for paraidxpath in para_childspan_details[pmid][prot]:
                if paraidxpath not in htmlIdx_paraId[pmid]:
                    print "User has highlighted phrases that are not contained within \"targetted\" <p> and <h3> tags."
                    print "Ignoring current highlight!"
                    continue
                contg_highlights_dict[pmid][prot][paraidxpath] = {}
                childnode = getnodefromidxpath(paraidxpath, copy_soup.body)     # childnode may have both yellow and white highlights
                for color in ["Yellow", "White"]:
                    contg_passages = getContiguousPassages(para_childspan_details[pmid][prot][paraidxpath], childnode, color)
                    if len(contg_passages) > 0:
                        for each_psg in contg_passages:  # contg_passages is a list of dicts
                            sent_details = getSentenceList(article_data["paragraphList"][htmlIdx_paraId[pmid][paraidxpath]]["sentenceList"]["@items"], each_psg)
                            each_psg["sentenceDetails"] = sent_details
                            for each_sent in sent_details:  # sent_details is a list of dicts
                                if color == "Yellow":
                                    if each_sent["absoluteId"] in posprotdict[pmid][label][prot]:
                                        posprotdict[pmid][label][prot][each_sent["absoluteId"]].append(each_sent)   # There may be multiple contiguous phrases in a single sentence
                                    else:
                                        posprotdict[pmid][label][prot][each_sent["absoluteId"]] = [each_sent]
                                else:
                                    if "partialSentence" in each_sent:
                                        if "HasProtein" in each_sent:
                                            each_sent.pop("HasProtein", None)
                                        if not checkProtein(each_sent):  # Call checkProtein for only those partial sentences that were unhighlighted
                                            continue                                                        # by user for a protein that was given feedback on
                                        else:
                                            each_sent["HasProtein"] = "Yes"
                                    if each_sent["absoluteId"] in negprotdict[pmid][label][prot]:
                                        negprotdict[pmid][label][prot][each_sent["absoluteId"]].append(each_sent)
                                    else:
                                        negprotdict[pmid][label][prot][each_sent["absoluteId"]] = [each_sent]
                        contg_highlights_dict[pmid][prot][paraidxpath][color] = contg_passages

            if not negprotdict[pmid][label][prot]:
                negprotdict[pmid][label].pop(prot)
            else:
                for absid in negprotdict[pmid][label][prot]:
                    if len(negprotdict[pmid][label][prot][absid]) > 1:
                        newlist = sorted(negprotdict[pmid][label][prot][absid], key=itemgetter('partialOffset'))
                        negprotdict[pmid][label][prot][absid] = copy.deepcopy(newlist)
            if not posprotdict[pmid][label][prot]:
                posprotdict[pmid][label].pop(prot)
            else:
                for absid in posprotdict[pmid][label][prot]:  # prot is the protein for which the user has provided feedback
                    if len(posprotdict[pmid][label][prot][absid]) > 1:  # A protein may have either the full sentence OR multiple mutually exclusive parts of that sentence, but not both
                        newlist = sorted(posprotdict[pmid][label][prot][absid], key=itemgetter('partialOffset'))
                        posprotdict[pmid][label][prot][absid] = copy.deepcopy(newlist)

    # print "Update openaccess_posprotdict.p and openaccess_negprotdict.p"
    # for pmid in posprotdict:
    #     allposSentences = {}
    #     for prot in posprotdict[pmid][label]:
    #         for absid in posprotdict[pmid][label][prot]:
    #             if absid in allposSentences:
    #                 if prot not in allposSentences[absid]:
    #                     allposSentences[absid].append(prot)
    #             else:
    #                 allposSentences[absid] = [prot]
    #
    #     for prot in negprotdict[pmid][label]:
    #         for absid in negprotdict[pmid][label][prot].keys():
    #             if absid in allposSentences:
    #                 if checkfullsentence(absid, allposSentences[absid], posprotdict[pmid][label]):
    #                     negprotdict[pmid][label][prot].pop(absid)
    #                 else:   # All the proteins in allposSentences[absid], each have parts (one or multiple mutually exlcusive) of the sentence - absid
    #                     if len(negprotdict[pmid][label][prot][absid]) == 1 and "partialSentence" not in negprotdict[pmid][label][prot][absid][0]:
    #                         negprotdict[pmid][label][prot].pop(absid)
    #                     else:   # A protein may have either the full sentence OR multiple mutually exclusive parts of that sentence, but not both
    #                         newlist = []
    #                         for partsent in negprotdict[pmid][label][prot][absid]:
    #                             if not checkoverlap(partsent, absid, allposSentences[absid], posprotdict[pmid][label]):
    #                                 newlist.append(partsent)
    #                         if len(newlist) == 0:
    #                             negprotdict[pmid][label][prot].pop(absid)
    #                         else:
    #                             negprotdict[pmid][label][prot][absid] = newlist
    #
    #     for prot in negprotdict[pmid][label].keys():
    #         if not negprotdict[pmid][label][prot]:
    #             negprotdict[pmid][label].pop(prot)

    pickle.dump(posprotdict, open(os.path.join(drname, "openaccess_posprotdict_WithFeedback.p"), "wb"))
    pickle.dump(negprotdict, open(os.path.join(drname, "openaccess_negprotdict_WithFeedback.p"), "wb"))
    # json_fn = open(os.path.join(drname, "contiguous_highlights_details.json"), 'w')
    # json.dump(contg_highlights_dict, json_fn, indent=4, ensure_ascii=False)
    # json_fn.close()


def getBeautifulSoup(download_dir, pmcid):
    fname_string = '%s' % (os.path.join(download_dir, pmcid.strip(), pmcid.strip() + ".html"))
    fname = glob.glob(fname_string)
    if len(fname) == 0:
        print "The HTML file was not found for PMCID: ", pmcid.strip()
        exit()
    myfile = fname[0].strip()
    soup = BeautifulSoup(open(myfile), "lxml")
    return soup


def getExtractedSentences(pmid, rootDir):
    fname_string = '%s*' % (os.path.join(rootDir, pmid.strip()))
    fname = glob.glob(fname_string)
    if len(fname) == 0:
        print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
        exit()
    myfile = fname[0].strip()
    file_handle = open(myfile, "rb")
    myfile_data = json.load(file_handle)  # myfile_data is a list of dictionaries
    file_handle.close()
    return myfile_data


def populate_negprotdict(dest_dict, pmid, prot, source_dict):
    if prot in dest_dict:
        print "A positive protein already present in negprotdict: ", prot, "\t for PMID: ", pmid
        exit()
    dest_dict[prot] = copy.deepcopy(source_dict)


def populatespandict(span, spanidx, spandict):
    spandict["idx"] = spanidx
    spandict["highlightText"] = span[1]
    spandict["spanPath"] = span[2]
    spandict["highlightOffset"] = int(span[3])
    spandict["highlightLength"] = int(span[4])
    details = re.findall('"([^"]*)"', span[0])
    for det in details:
        if det.isdigit():
            spandict["data-timestamp"] = det
            continue
        if det.startswith("background-color"):
            spandict["style"] = det
            if "rgb(255, 255, 123)" in det:
                spandict["highlightColor"] = "Yellow"
            else:
                spandict["highlightColor"] = "White"
    spandict["data-highlighted"] = "true"
    spandict["class"] = "highlighted"


def getEnclosingP_H3(childpath, spanidx, copy_soup):
    childnode = getnodefromidxpath(childpath, copy_soup.body)
    if childnode.name == "p" or childnode.name == "h3":
        return childpath, str(spanidx)
    foundp = False; childpathidxlist = childpath.split(":")
    childrelpathidx = [childpathidxlist.pop()]  # childrelpathidx stores the index of childnode at this point
    for parent in childnode.parents:
        if parent.name == "body":
            break
        if parent.name == "p" or parent.name == "h3":
            foundp = True; break
        childrelpathidx.insert(0, childpathidxlist.pop())
    if foundp:
        childrelpath = ":".join(childrelpathidx)
        parentidxpath = childpath[0:childpath.rfind(childrelpath)-1]
        return parentidxpath, childrelpath+":"+str(spanidx)
    else:
        return "", ""       # User has highlighted text that is not enclosed in a <p> or a <h3> tag


def checkBadHighlights(pmid, prot, spandict):
    if pmid == "16009723" and prot == "mdc1" and spandict["spanPath"] == "3:3:13:7:13:2:6:0:4":     # Highlight is inside a <table>
        return True
    if pmid =="24106086" and prot == "jak3":
        if spandict["spanPath"].startswith("3:5:3:9:"):
            spandict["spanPath"] = "3:5:3:7:" + spandict["spanPath"][8:]
        if spandict["spanPath"].startswith("3:5:3:14:"):
            spandict["spanPath"] = "3:5:3:10:" + spandict["spanPath"][9:]
        return False


def highlightnode(pmid, prot, copy_soup, spandict, mydict):
    badhighlight = checkBadHighlights(pmid, prot, spandict)
    if badhighlight:    return ""

    childidxpath = ":".join(spandict["spanPath"].split(":")[:-1])
    (parentidxpath, childrelpath) = getEnclosingP_H3(childidxpath, spandict["idx"], copy_soup)
    if parentidxpath == "":
        return ""
    else:
        if parentidxpath not in mydict:
            mydict[parentidxpath] = {}
        if childrelpath in mydict[parentidxpath]:
            print "The child relative path already exists for this paragraph! How can that be?!"
            exit()
        if parentidxpath+":"+childrelpath != spandict["spanPath"]:
            print "I have not been able to split the paths correctly! Shame on me :("
            exit()
        mydict[parentidxpath][childrelpath] = spandict

    childnode = getnodefromidxpath(childidxpath, copy_soup.body)
    elIndex = spandict["idx"]
    if elIndex > 0 and isinstance(getnodefromidxpath(str(elIndex-1), childnode), NavigableString):
        elIndex -= 1
    navStrNode = getnodefromidxpath(str(elIndex), childnode)
    if not isinstance(navStrNode, NavigableString):
        print "The node is not a NavigableString"
        print childnode
        exit()

    strleft = navStrNode[0:spandict["highlightOffset"]]
    mystr = navStrNode[spandict["highlightOffset"]:spandict["highlightOffset"]+spandict["highlightLength"]]
    assert mystr == spandict["highlightText"]
    strright = navStrNode[spandict["highlightOffset"]+spandict["highlightLength"]:]

    span_tag = copy_soup.new_tag("span")
    span_tag = populatespanattr(pmid, span_tag, spandict)
    span_tag.append(NavigableString(mystr))

    if len(strleft) != 0:
        navStrNode.replace_with(NavigableString(strleft))
        childnode.insert(elIndex+1, span_tag)
        elIndex += 1
    else:
        navStrNode.replace_with(span_tag)

    if len(strright) != 0:
        childnode.insert(elIndex+1, NavigableString(strright))

    return childnode


def getContiguousPassages(paraDetails, hilitednode, color):
    mylist = []
    nextOffset = 0; startofpsg = False
    for child in hilitednode.descendants:   # We scan through the ENTIRE <p> or <h3> node
        if isinstance(child, NavigableString):  # We will only concentrate on strings whose parents are span tags
            childstr = child
            currOffset = nextOffset
            nextOffset = currOffset + len(childstr)
            if child.parent is hilitednode:     # All highlighted phrases will be enclosed in span tags, thus the parent of any highlighted phrase cannot be a <p> or a <h3> tag
                if startofpsg:
                    addtoList(mylist, psgText, psgOffset, psgLength)
                startofpsg = False      # If child.parent is either <p> or <h3>. This indicates a break in contiguity
                continue
            spanrelidxpath = getchildpathindex(child.parent, hilitednode)   # child.parent may be a highlighted span tag or may be any other tag like "<sup>", "<b>", "<a ref>"
            if spanrelidxpath in paraDetails and paraDetails[spanrelidxpath]["highlightColor"] == color:
                if not startofpsg:  # Detected the start of a new contiguous passage
                    startofpsg = True
                    psgOffset = currOffset
                    psgLength = len(childstr)
                    psgText = childstr
                else:   # Continuation of a contiguous passage
                    psgLength += len(childstr)
                    psgText += childstr
            else:       # If child.parent is not a yellow colored highlighted span node, i.e. if its any other node
                if startofpsg:
                    addtoList(mylist, psgText, psgOffset, psgLength)    # This indicates a break in contiguity
                startofpsg = False
    if startofpsg:      # If the last sentence of a passage is highlighted
        addtoList(mylist, psgText, psgOffset, psgLength)

    return mylist


def addtoList (mylist, text, offset, length):
    tempdict = {}
    tempdict["psgText"] = text
    tempdict["psgOffset"] = offset
    tempdict["psgLength"] = length
    mylist.append(tempdict)


def getSentenceList(origList, psg):
    psgtext = psg["psgText"]; psgoffset = psg["psgOffset"]
    psglength = psg["psgLength"]; sentenceList = []
    (psgtext, psgoffset, psglength) = skipleadingspace(psgtext, psgoffset, psglength)
    for sentence in origList:
        if len(psgtext.strip()) == 0:
            break
        if psgoffset > sentence["charOffset"]:
            if psgoffset > sentence["charOffset"] + sentence["sentenceLen"]:
                continue
            elif psgoffset <= sentence["charOffset"] + sentence["sentenceLen"]:
                if psglength > (sentence["charOffset"] + sentence["sentenceLen"]) - psgoffset:     # psg spills over to the next adjacent sentence
                    tempdict = copy.deepcopy(sentence)
                    tempdict["partialSentence"] = True
                    tempdict["partialOffset"] = psgoffset
                    partialText = sentence["sentenceText"][psgoffset - sentence["charOffset"]:]
                    tempdict["partialText"] = checkBoundaries(partialText, sentence, psgoffset - sentence["charOffset"])
                    tempdict["partialcleansedText"] = " ".join(tempdict["partialText"].strip(' \t\n\r').split())
                    tempdict["partialLength"] = len(tempdict["partialText"])
                    assert partialText == psgtext[:len(partialText)]
                    psgtext = psgtext[tempdict["partialLength"]:]
                    psgoffset += tempdict["partialLength"]
                    psglength -= tempdict["partialLength"]
                    (psgtext, psgoffset, psglength) = skipleadingspace(psgtext, psgoffset, psglength)
                    if tempdict["partialcleansedText"] != sentence["cleansedText"]:
                        sentenceList.append(tempdict)
                    else:
                        sentenceList.append(sentence)
                elif psglength <= (sentence["charOffset"] + sentence["sentenceLen"]) - psgoffset:      # The equality condition can NEVER happen because psgoffset > sentence["charOffset"]
                    tempdict = copy.deepcopy(sentence)
                    tempdict["partialSentence"] = True
                    tempdict["partialOffset"] = psgoffset
                    partialText = sentence["sentenceText"][psgoffset - sentence["charOffset"]:(psgoffset - sentence["charOffset"]) + psglength]
                    tempdict["partialText"] = checkBoundaries(partialText, sentence, psgoffset - sentence["charOffset"])
                    tempdict["partialcleansedText"] = " ".join(tempdict["partialText"].strip(' \t\n\r').split())
                    tempdict["partialLength"] = len(tempdict["partialText"])
                    assert partialText == psgtext
                    psgtext = psgtext[tempdict["partialLength"]:]    # This should be an empty string now
                    if tempdict["partialcleansedText"] != sentence["cleansedText"]:
                        sentenceList.append(tempdict)
                    else:
                        sentenceList.append(sentence)
            #elif psgoffset == sentence["charOffset"] + sentence["sentenceLen"]:     # psgoffset points to whitespace between sentence boundaries
            #    print "Something went wrong during partitioning of contiguous passages into sentences"
            #    exit()
        elif psgoffset == sentence["charOffset"]:
            if psglength > sentence["sentenceLen"]:
                sentenceList.append(sentence)
                psgtext = psgtext[sentence["sentenceLen"]:]
                psgoffset += sentence["sentenceLen"]
                psglength -= sentence["sentenceLen"]
                (psgtext, psgoffset, psglength) = skipleadingspace(psgtext, psgoffset, psglength)
            elif psglength == sentence["sentenceLen"]:
                sentenceList.append(sentence)
                psgtext = psgtext[sentence["sentenceLen"]:]         # This should be an empty string now
            elif psglength < sentence["sentenceLen"]:
                tempdict = copy.deepcopy(sentence)
                tempdict["partialSentence"] = True
                tempdict["partialOffset"] = psgoffset
                partialText = sentence["sentenceText"][0:psglength]
                tempdict["partialText"] = checkBoundaries(partialText, sentence, 0)
                tempdict["partialcleansedText"] = " ".join(tempdict["partialText"].strip(' \t\n\r').split())
                tempdict["partialLength"] = len(tempdict["partialText"])
                assert partialText == psgtext
                psgtext = psgtext[tempdict["partialLength"]:]  # This should be an empty string now
                if tempdict["partialcleansedText"] != sentence["cleansedText"]:
                    sentenceList.append(tempdict)
                else:
                    sentenceList.append(sentence)
        elif psgoffset < sentence["charOffset"]:    # Should never happen!!! Because I am going through all the sentences in a paragraph in order
            print "Something went wrong during partitioning of contiguous passages into sentences"
            print "The passage offset is less than the sentence offset"
            exit()

    return sentenceList


def skipleadingspace(psgtext, psgoffset, psglength):
    numSpace = len(psgtext) - len(psgtext.lstrip())
    if numSpace > 0:
        psgtext = psgtext[numSpace:]
        psgoffset += numSpace
        psglength -= numSpace
    return (psgtext, psgoffset, psglength)


def checkProtein(sentdict):
    cleansed_text = " ".join(sentdict["partialText"].strip(' \t\n\r').split())
    if check_for_unicode(cleansed_text):
        cleansed_text = cleansed_text.decode("utf-8")
    new_text = cleansed_text.replace("\n", " ").replace("\r", "")
    matches = find_proteins(new_text.lower(), 19)

    # matches = []; ctr = 0; old_text = ""
    # while len(matches) == 0 and ctr <= 3:
    #     if new_text != old_text:
    #         matches = find_proteins(new_text.lower(), 19)
    #     if len(matches) > 0:    break
    #     ctr += 1; old_text = new_text
    #     if ctr == 1:
    #         if "-" in cleansed_text:
    #             new_text = cleansed_text.replace("-", " ")
    #     if ctr == 2:
    #         if "-" in cleansed_text:
    #             new_text = cleansed_text.replace("-", "")
    #     if ctr == 3:
    #         if "/" in cleansed_text:
    #             new_text = cleansed_text.replace("/", " ")

    if len(matches) > 0:
        for m in matches:
            match_name = m[0].strip()
            exact_mention = cleansed_text[m[1][0]:m[1][1]]
            if isstopword(match_name, exact_mention):
                continue
            else:
                return True
    return False


def checkBoundaries(partialText, sentence, reloffset):
    origPartialText = partialText
    mysent = sentence["sentenceText"]
    mysent_len = sentence["sentenceLen"]
    if mysent.startswith(partialText):                                  # Need to check the end sentence boundary
        end_idx = len(partialText)
        if not mysent[end_idx - 1].isspace() and not isPunct(mysent[end_idx - 1]):
            while end_idx < mysent_len and not mysent[end_idx].isspace() and not isPunct(mysent[end_idx]):
                partialText = partialText + mysent[end_idx]
                end_idx += 1
    elif mysent.endswith(partialText):                                  # Need to check the starting sentence boundary
        start_idx = reloffset - 1
        if not mysent[reloffset].isspace() and not isPunct(mysent[reloffset]):  # Its not a space AND its not a punct
            while start_idx > -1 and not mysent[start_idx].isspace() and not isPunct(mysent[start_idx]):
                partialText = mysent[start_idx] + partialText
                start_idx -= 1
    else:                                                               # Check both sentence boundaries
        start_idx = reloffset - 1
        end_idx = reloffset + len(partialText)
        if not mysent[end_idx - 1].isspace() and not isPunct(mysent[end_idx - 1]):
            while end_idx < mysent_len and not mysent[end_idx].isspace() and not isPunct(mysent[end_idx]):
                partialText = partialText + mysent[end_idx]
                end_idx += 1
        if not mysent[reloffset].isspace() and not isPunct(mysent[reloffset]):
            while start_idx > -1 and not mysent[start_idx].isspace() and not isPunct(mysent[start_idx]):
                partialText = mysent[start_idx] + partialText
                start_idx -= 1

    # if origPartialText != partialText:
    #     print "Original Partial Text: ", origPartialText
    #     print "New Partial Text: ", partialText
    return partialText


def isPunct(mychar):
    for c in set(string.punctuation):           # "DDX41-shrna", "IL-4" and "IL-5"
        if c != "-" and c == mychar: return 1   # Sentence does contain the word "IL" standalone, ofcourse!
    return 0


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


def getspanid(pmid):
    global spanctr
    spanctr += 1
    return pmid + "_" + str(spanctr)


def populatespanattr(pmid, span_tag, spandict):
    span_tag["id"] = getspanid(pmid)
    span_tag["data-highlighted"] = "true"
    span_tag["data-timestamp"] = spandict["data-timestamp"]
    span_tag["class"] = "highlighted"
    span_tag["style"] = spandict["style"]
    return span_tag


def flattenSpanDict(DictIn, mylist, parent_key):
    for key, value in sorted(DictIn.items()):
        if isinstance(value, dict):
            flattenSpanDict(value, mylist, concatenate_keys(parent_key, key))
        elif isinstance(value, list):
            if value[0]["spanPath"] != concatenate_keys(parent_key, key):
                print "Something went wrong trying to flatten the dict"
                print value[0]["spanPath"]; print concatenate_keys(parent_key, key); pprint(value)
                exit()
            mylist.append(value[0])


def concatenate_keys(parent_key, key):
    if parent_key == "":
        return str(key)
    else:
        return parent_key + ":" + str(key)


def checkfullsentence(absid, proteinlist, posSentences):
    for prot in proteinlist:
        for sent in posSentences[prot][absid]:
            if len(posSentences[prot][absid]) == 1 and "partialSentence" not in sent:
                return True
    return False


def checkoverlap(partsent, absid, proteinlist, posSentences):
    for prot in proteinlist:
        for sent in posSentences[prot][absid]:
            if partsent["partialOffset"] > sent["partialOffset"] and partsent["partialOffset"] < sent["partialOffset"] + sent["partialLength"]:
                return True
    return False


if __name__ == "__main__":
    main_body()
