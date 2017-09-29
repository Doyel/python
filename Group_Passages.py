# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
import sys
import glob
import json
reload(sys)
sys.setdefaultencoding("utf-8")
import os
import collections
from bs4 import BeautifulSoup
from bs4 import NavigableString
from bs4 import Tag
from pprint import pprint
import cPickle as pickle
import copy
import difflib
import string
from string import punctuation
import nltk.data
import re

# python Group_Passages.py ./PubMedIDS_with_Extras_OpenAccess.txt ./vw_RNAi_Training_File_Jul7.txt ./protein_detect_output ./downloaded_articles
# python Group_Passages.py ./Test_PMIDs.txt ./vw_RNAi_Training_File_Jul7.txt ./protein_detect_output ./downloaded_articles

def main_body():
    global sent_detector; global reg_expr; global spanctr; global problem_sentences; global manual_sentences_match; global ignore_sentences
    parser = argparse.ArgumentParser(prog='Group_Passages', usage='Group_Passages.py <PubMedfilelist> <training_file> <proteins_detected_dir> <downloaded_articles_dir>', description='Script to find and group adjacent passages')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids that are in OpenAccess and are also in John list')
    parser.add_argument('training_file', help='Training File')
    parser.add_argument('protein_outputdir', help='Path where files containing the detected proteins are located')
    parser.add_argument('downloaded_articles_dir', help='Path where openaccess articles were downloaded')
    args = parser.parse_args()

    sent_detector = nltk.data.load('tokenizers/punkt/english.pickle')
    reg_expr = re.compile('^[S]?[0-9][A-Za-z]?\W*')
    problem_sentences = {}
    with open("Problem_Sentence_Matches.txt", 'rt') as f:
        for line in f:
            if line.strip() == "": continue
            line = line.strip()
            details = line.split("|_|")
            details[1] = unicode(details[1], "utf-8")
            details[2] = unicode(details[2], "utf-8")
            if details[0] in problem_sentences:
                if details[1] in problem_sentences[details[0]]:
                    problem_sentences[details[0]][details[1]].append(details[2])
                else:
                    problem_sentences[details[0]][details[1]] = [details[2]]
            else:
                problem_sentences[details[0]] = {}
                problem_sentences[details[0]][details[1]] = [details[2]]
    #pprint(problem_sentences); exit()

    manual_sentences_match = {}
    with open("Manual_Sentence_Matches.txt", 'rt') as f:
        for line in f:
            if line.strip() == "": continue
            line = line.strip()
            details = line.split("|_|")
            details[1] = unicode(details[1], "utf-8")
            details[2] = unicode(details[2], "utf-8")
            if details[0] in manual_sentences_match:
                if details[1] in manual_sentences_match[details[0]]:
                    manual_sentences_match[details[0]][details[1]][details[2]] = details[3]
                else:
                    manual_sentences_match[details[0]][details[1]] = {}
                    manual_sentences_match[details[0]][details[1]][details[2]] = details[3]
            else:
                manual_sentences_match[details[0]] = {}
                manual_sentences_match[details[0]][details[1]] = {}
                manual_sentences_match[details[0]][details[1]][details[2]] = details[3]
    #pprint(manual_sentences_match); exit()

    ignore_sentences = {}
    with open("Ignore_Sentences.txt", 'rt') as f:
        for line in f:
            if line.strip() == "": continue
            line = line.strip()
            details = line.split("|_|")
            details[1] = unicode(details[1], "utf-8")
            if details[0] in ignore_sentences:
                ignore_sentences[details[0]].append(details[1])
            else:
                ignore_sentences[details[0]] = [details[1]]

    if os.path.isfile("prot_passages.p"):
        prot_passages = pickle.load(open("prot_passages.p", "rb"))
    else:
        print "The dictionary - prot_passages was not found in current directory. Creating it!"
        openaccess = []; id_passages = {}
        with open(args.PubMedfilelist, 'rt') as f:
            for pmid in f:
                if not pmid.strip() or pmid.startswith("#"): continue

                fname_string = '%s*Paragraphs*' % (os.path.join(args.protein_outputdir, pmid.strip()))
                fname = glob.glob(fname_string)
                if len(fname) == 0:
                    # print "No proteins were detected for PubMedId: ", pmid.strip()
                    continue

                file_handle = open(fname[0].strip(), "rb")
                myfile = json.load(file_handle)  # myfile is a list of dictionaries
                file_handle.close()

                id_passages[pmid.strip()] = {}
                for elem in myfile:
                    mystring = elem["paragraph"].replace('\n', ' ').replace('\r', ' ')
                    id_passages[pmid.strip()][elem["id"]] = " ".join(mystring.strip(' \t\n\r').split())

                openaccess.append(pmid.strip())
        print "In the Datum KB, no. of articles that are in OpenAccess and also in John's list: ", len(openaccess)

        print "Collecting all training instances having a label of 1"
        prot_passages = {}
        with open(args.training_file, 'rt') as f1:
            for line in f1:
                label_tag = line.split("|")[0]
                label = label_tag.split()[0].strip()
                id_prots = label_tag.split()[1].split("--")
                pmid = id_prots[0].strip().split("_")[0].strip()
                id = int(id_prots[0].strip().split("_")[1].strip())
                prots = list(set(id_prots[1:]))

                if pmid in openaccess and label == "1":
                    if pmid not in prot_passages:
                        prot_passages[pmid] = {}
                    for protein in prots:
                        if protein not in prot_passages[pmid]:
                            prot_passages[pmid][protein] = {}   # prot_passages[pmid][protein] is a dictionary
                        prot_passages[pmid][protein][id] = id_passages[pmid][pmid + "_" + str(id)]
        #pprint(prot_passages["21107320"]); exit()


        print "Grouping sentences into paragraphs"
        for pmid in prot_passages.keys():
            print "PMID: ", pmid
            for protein in prot_passages[pmid]:
                od = collections.OrderedDict(sorted(prot_passages[pmid][protein].items()))
                #print "Divya: ", json.dumps(od, indent=4)
                mygroup = ""; templist = []
                for id, psg in od.iteritems():
                    #print "id: ", id
                    #if id == 44:
                        #print mygroup.encode('utf-8'); print "\n"; print psg.encode('utf-8')
                        #exit()
                    if mygroup == "":
                        mygroup = psg; myidlist = [id]
                        if next(reversed(od)) == id:  # current id is the last entry of the ordered dict, then wrap up!
                            tempdict = {}
                            tempdict["id_list"] = myidlist
                            tempdict["grouped_evid"] = mygroup
                            templist.append(tempdict)
                    else:
                        common = find_hooks(mygroup, psg)
                        #print "Saswati Common: ", common.encode('utf-8')
                        #if id == 44:
                            #exit()
                        if common:
                            mygroup = common
                            myidlist.append(id)
                            if next(reversed(od)) == id:    # current id is the last entry of the ordered dict, then wrap up!
                                tempdict = {}
                                tempdict["id_list"] = myidlist
                                tempdict["grouped_evid"] = mygroup
                                templist.append(tempdict)
                        else:
                            tempdict = {}
                            tempdict["id_list"] = myidlist
                            tempdict["grouped_evid"] = mygroup   # tempdict["groupd_evid"] is a string currently
                            templist.append(tempdict)
                            mygroup = psg; myidlist = [id]
                prot_passages[pmid][protein] = templist     # prot_passages[pmid][protein] is a list
            #pprint(prot_passages[pmid])
            #exit()

        print "Locating paragraphs in HTML documents"
        pm_pmc = pickle.load(open("pmid_pmcid_Extras.p", "rb"))
        for pmid in prot_passages.keys():
            print "PMID: ", pmid, "\tPMCID: ", pm_pmc[pmid]
            soup = getBeautifulSoup(args.downloaded_articles_dir, pm_pmc[pmid])
            for protein in prot_passages[pmid]:
                soup_copy = copy.copy(soup); prevhtmlpathidx = -1
                for mydict in prot_passages[pmid][protein]:     # prot_passages[pmid][protein] is a list of dictionaries
                    evidence = mydict["grouped_evid"]    # I need to find the contents of evidence in the html article
                    mydict["located_evid"] = list(parse_html(evidence, soup_copy, prevhtmlpathidx))    # mydict["located_evid"] is a list of dictionaries
                    prevhtmlpathidx = mydict["located_evid"][-1]["htmlpathindex"]
                    #pprint(mydict)

            pprint(prot_passages[pmid])
        pickle.dump(prot_passages, open("prot_passages.p", "wb"))
        exit()

    print "Creating highlight objects"
    pm_pmc = pickle.load(open("pmid_pmcid_Extras.p", "rb"))
    located_prot_passages = {}
    for pmid in prot_passages.keys():
        print "PMID: ", pmid, "\tPMCID: ", pm_pmc[pmid]
        orig_soup = getBeautifulSoup(args.downloaded_articles_dir, pm_pmc[pmid])
        spanctr = -1
        located_prot_passages[pmid] = {}
        for protein in prot_passages[pmid]:
            soup = copy.copy(orig_soup); body = soup.body   # referrence node
            tempdict = {}; located_prot_passages[pmid][protein] = {}
            for mydict in prot_passages[pmid][protein]:  # prot_passages[pmid][protein] is a list of dictionaries
                for evid in mydict["located_evid"]:     # mydict["located_evid"] is a list of dictionaries
                    if evid["htmlpath"] in tempdict:
                        tempdict[evid["htmlpath"]].append(evid["evidence"])
                    else:
                        tempdict[evid["htmlpath"]] = [evid["evidence"]]
            located_prot_passages[pmid][protein]["grouped_htmlpaths"] = tempdict

            highlight_object = []
            for htmlpath in located_prot_passages[pmid][protein]["grouped_htmlpaths"]:
                if htmlpath.strip() == "":  continue
                #print htmlpath
                targetnode = getnodefromxpath(htmlpath, body)   # target node that needs to be highlighted
                #print targetnode
                #print located_prot_passages[pmid][protein]["grouped_htmlpaths"][htmlpath]
                #soup_copy = copy.copy(targetnode)
                for evid_sent in located_prot_passages[pmid][protein]["grouped_htmlpaths"][htmlpath]:
                    targetnode = updatenode(pmid, soup, targetnode, evid_sent)

                    #newnode = getnewnode(soup, targetnode, evid_sent, pmid)
                    #print newnode
                    #targetnode.replace_with(newnode)
                    #targetnode = getnodefromxpath(htmlpath, body)
                #pprint(located_prot_passages[pmid][protein]["grouped_htmlpaths"][htmlpath])
                print targetnode
                for i, child in enumerate(targetnode.descendants):
                    #print targetnode
                    if isinstance(child, Tag) and child.name == "span" and child.has_attr('id') and child["id"].startswith(pmid):
                        childlen = len(child.string)
                        childindexlist = getchildpathindex(child, body)
                        if isinstance(child.previous_sibling, NavigableString):
                            offset = len(child.previous_sibling)
                        else:
                            offset = 0
                        span_tag_str = "<span data-highlighted=\"" + child["data-highlighted"] + "\" class=\"" + child["class"] + "\" style=\"" + child["style"] + "\"></span>"
                        #print span_tag_str
                        templist = [span_tag_str, child.string, ":".join(childindexlist), str(offset), str(childlen)]
                        highlight_object.append(templist)
                #pprint(highlight_object); exit()
            located_prot_passages[pmid][protein]["highlights"] = highlight_object


def getmatcheslist(pmid, mynode, evidence, mylist):
    if not isinstance(evidence, unicode):
        evidence = unicode(evidence, "utf-8")

    for child in mynode.descendants:
        if len(evidence.strip()) == 0 or evidence.strip() == ".": break

        if isinstance(child, NavigableString):
            childstr = child
            if not isinstance(childstr, unicode):
                childstr = unicode(childstr, "utf-8")
            print "Child Str: ", childstr.encode("utf-8");
            print "Evidence Sentence: ", evidence.encode("utf-8")
            found_psg = findmatch(pmid, evidence, childstr)  # found_psg is a clean string
            print "Found Passage: ", found_psg.encode("utf-8")
            if found_psg != "":

                if (evidence.find(found_psg) + len(found_psg) + 1) < len(evidence) and evidence[evidence.find(found_psg) + len(found_psg)].isspace():
                    print "Mamoni: ", evidence[evidence.find(found_psg): len(found_psg)].encode("utf-8"), "***", evidence[evidence.find(found_psg) + len(found_psg)]
                    print "Divya: ", evidence[evidence.find(found_psg) + len(found_psg) + 1:].encode("utf-8")
                    evidence = evidence[evidence.find(found_psg) + len(found_psg) + 1:]  # found_psg will always be at the beginning of evidence
                else:
                    print "Saswati: ", evidence[evidence.find(found_psg) + len(found_psg):].encode("utf-8")
                    evidence = evidence[evidence.find(found_psg) + len(found_psg):]
                curr_parent = child.parent
                child_idx = getchildindex(curr_parent, child)
                parent_path = get_htmlpath(curr_parent)
                mylist.append([parent_path, curr_parent, child, found_psg, child_idx])
    return (evidence, mylist)


def updatenode(pmid, soup, mynode, evidence):
    mylist = []
    evidence, mylist = getmatcheslist(pmid, mynode, evidence, mylist)

    if len(evidence.strip()) > 0 and evidence.strip() != ".":
        if pmid in ignore_sentences:
            for evid_sent in ignore_sentences[pmid]:
                if checkmatch(evid_sent.strip(), evidence.strip(), 0.9) and checkmatch(evidence.strip(),
                                                                                       evid_sent.strip(), 0.9):
                    evidence = ""
        else:
            num_words = [x.strip() for x in evidence.split() if len(x.strip()) > 0]
            if not(len(num_words) == 1 and num_words[0].endswith(".")):
                print pmid
                print "Could not find part of the evidence sentence: ", evidence.encode("utf-8")
                print mynode
                exit()

    mydict = {}
    for elem in mylist:
        if elem[0] in mydict:
            mydict[elem[0]].append([elem[1], elem[2], elem[3], elem[4]])  # [curr_parent, child, found_psg, child_idx]
        else:
            mydict[elem[0]] = [[elem[1], elem[2], elem[3], elem[4]]]

    for delem in mydict:
        mydict[delem].sort(key=lambda x: int(x[3]), reverse=True)

    for ppath in mydict:
        for lelem in mydict[ppath]:  # mydict[ppath] is a list of lists. lelem is a 4-element list: [curr_parent 0, child 1, found_psg 2, child_idx 3]
            childstr = lelem[1];
            found_psg = lelem[2];
            child_idx = lelem[3]
            # print lelem[0]; print childstr; print child_idx
            parts = partitionChild(pmid, childstr, found_psg)
            # print parts; #exit()
            for i, (strg, switch) in enumerate(parts):
                navg_strg = NavigableString(strg)
                if switch:
                    span_tag = soup.new_tag("span")
                    span_tag = populatespanattr(pmid, span_tag)
                    span_tag.append(navg_strg)
                    if i == 0:
                        childstr.replace_with(span_tag)
                        # print lelem[0]; exit()
                    else:
                        lelem[0].insert(child_idx, span_tag)
                else:
                    if i == 0:
                        childstr.replace_with(navg_strg)
                    else:
                        lelem[0].insert(child_idx, navg_strg)
                child_idx += 1
    return mynode


def getnodefromxpath(xpath, prevnode):
    nodes = xpath.split(' > ')
    #print nodes
    while len(nodes) > 0:
        currnodestr = nodes.pop(0)
        if ':' in currnodestr:
            childno = currnodestr.split(':')[1]
            currnodestr = currnodestr.split(':')[0]
        else:
            childno = '1'
        currno = 1
        for child in prevnode.children:
            if child.name == currnodestr:
                if currno == int(childno):
                    prevnode = child
                    break
                else:
                    currno += 1
    return prevnode


def getnewnode(soup, mynode, evidence, pmid):
    my_new_tag = soup.new_tag(mynode.name)
    for i, child in enumerate(mynode.children):

        if len(evidence.strip()) == 0 or evidence.strip() == ".": break

        if isinstance(child, NavigableString):
            childstr = child                    # childstr may have multiple newlines and spaces in it
            found_psg = findmatch(evidence, childstr)  # found_psg is a clean, well formatted string
            if found_psg == "":
                navg_strg = NavigableString(childstr)
                my_new_tag.append(navg_strg)
            else:
                if (evidence.find(found_psg) + len(found_psg) + 1) < len(evidence) and evidence[evidence.find(found_psg) + len(found_psg) + 1].isspace():
                    evidence = evidence[evidence.find(found_psg) + len(found_psg) + 1:]  # found_psg will always be at the beginning of evidence
                else:
                    evidence = evidence[evidence.find(found_psg) + len(found_psg):]
                parts = partitionChild(childstr, found_psg)
                for strg, switch in parts:
                    navg_strg = NavigableString(strg)
                    if switch:
                        span_tag = soup.new_tag("span")
                        span_tag = populatespanattr(pmid, span_tag)
                        span_tag.append(navg_strg)
                        my_new_tag.append(span_tag)
                    else:
                        my_new_tag.append(navg_strg)

        else:
            new_child_tag = copy.copy(child)
            for greatchild in child.descendants:    # child is a tag. It may be a nested tag
                if isinstance(greatchild, NavigableString):
                    curr_parent = greatchild.parent
                    childstr = greatchild
                    #child_idx = getchildindex(curr_parent, greatchild)

                    found_psg = findmatch(evidence, childstr)  # found_psg is a clean string
                    if found_psg != "":
                        child_clean = " ".join(child.strip(' \t\n\r').split())
                        if len(found_psg) < 0.9 * len(child_clean):
                            print childstr; print found_psg; print evidence; exit()

                        if (evidence.find(found_psg) + len(found_psg) + 1) < len(evidence) and evidence[evidence.find(found_psg) + len(found_psg) + 1].isspace():
                            evidence = evidence[evidence.find(found_psg) + len(found_psg) + 1:]  # found_psg will always be at the beginning of evidence
                        else:
                            evidence = evidence[evidence.find(found_psg) + len(found_psg):]
                        span_tag = soup.new_tag("span")
                        span_tag = populatespanattr(pmid, span_tag)
                        greatchild.wrap(span_tag)

            my_new_tag.append(new_child_tag)

    if len(evidence) > 0 and evidence.strip() != ".":
        print "Could not find part of the evidence sentence: ", evidence.encode("utf-8")
        exit()
    return my_new_tag


def findmatch(pmid, evidence, child):
    # child will have lots of whitespace characters within it. evidence will usually be a well-formed, well-formatted string
    evidence_clean = " ".join(evidence.strip(' \t\n\r').split())
    print "Clean Evidence: ", evidence_clean.encode("utf-8")
    child_clean = " ".join(child.strip(' \t\n\r').split())
    print "Clean Child: ", child_clean.encode("utf-8")

    if not isinstance(evidence_clean, unicode):
        evidence_clean = unicode(evidence_clean, "utf-8")
    if not isinstance(child_clean, unicode):
        child_clean = unicode(child_clean, "utf-8")

    if pmid in problem_sentences:
        for child_sent in problem_sentences[pmid]:
            if checkmatch(child_clean.strip(), child_sent.strip(), 0.9) and checkmatch(child_sent.strip(),
                                                                                       child_clean.strip(), 0.9):
                print "***************** Coming here! ******************"
                for evid_sent in problem_sentences[pmid][child_sent.strip()]:
                    print "Dictionary Key: ", child_sent.encode("utf-8")
                    print "Clean Evidence: ", evidence_clean.encode("utf-8")
                    print "Dictionary Evidence: ", evid_sent.encode("utf-8")
                    s = difflib.SequenceMatcher(None, evidence_clean.strip(), evid_sent.strip())
                    pos_a1, pos_b1, size1 = s.find_longest_match(0, len(evidence_clean.strip()), 0,
                                                                 len(evid_sent.strip()))
                    # print size1 #
                    if checkmatch(evidence_clean.strip(), evid_sent.strip(), 0.9):
                        return ""

    lcs_match_pmid = ["22245064", "20141835"]
    if pmid in lcs_match_pmid:
        common_str = longest_common_substring(evidence_clean, child_clean)
        if len(common_str.strip()) > 1 and evidence_clean.startswith(common_str):
            print "LCS"
            return common_str

    if pmid in manual_sentences_match:
        for child_sent in manual_sentences_match[pmid]:
            if checkmatch(child_clean.strip(), child_sent.strip(), 0.9) and checkmatch(child_sent.strip(),
                                                                                       child_clean.strip(), 0.9):
                for evid_sent in manual_sentences_match[pmid][child_sent.strip()]:
                    if checkmatch(evidence_clean.strip(), evid_sent.strip(), 0.9) and checkmatch(evid_sent.strip(), evidence_clean.strip(), 0.9):
                        mysize = int(manual_sentences_match[pmid][child_sent.strip()][evid_sent.strip()])
                        return evidence_clean[:mysize]


        #evidence_tokens = [x.strip() for x in evidence_clean.split() if len(x.strip()) > 0]
        #child_tokens = [x.strip() for x in child_clean.split() if len(x.strip()) > 0]

        #first_token_idx = getIndex(child_tokens, evidence_tokens[0])

    s = difflib.SequenceMatcher(None, evidence_clean, child_clean)
    pos_a, pos_b, size = s.find_longest_match(0, len(evidence_clean), 0, len(child_clean))
    if pos_a == 0 and size > 0:

        period_space_list = ["19274086"]
        if pmid in period_space_list and ". )" in evidence_clean:
            matches = s.get_matching_blocks()
            for midx in range(len(matches)-1):
                if midx == 0 or abs(matches[midx][0] - matches[midx][1]) == 1:   continue
                else:   return evidence_clean[:size]
            return evidence_clean[:matches[-2][0]+matches[-2][2]]

        if not (size == 1 and pos_a+1 < len(evidence_clean) and pos_b+1 < len(child_clean) and
                child_clean[pos_b+1].isalpha() and child_clean[pos_b+1] != evidence_clean[pos_a+1]):

            if pmid in problem_sentences:
                for child_sent in problem_sentences[pmid]:
                    if checkmatch(child_clean.strip(), child_sent.strip(), 0.9) and checkmatch(child_sent.strip(), child_clean.strip(), 0.9):
                        print "***************** Coming here! ******************"
                        for evid_sent in problem_sentences[pmid][child_sent.strip()]:
                            print "Dictionary Key: ", child_sent.encode("utf-8")
                            print "Clean Evidence: ", evidence_clean.encode("utf-8")
                            print "Dictionary Evidence: ", evid_sent.encode("utf-8")
                            s = difflib.SequenceMatcher(None, evidence_clean.strip(), evid_sent.strip())
                            pos_a1, pos_b1, size1 = s.find_longest_match(0, len(evidence_clean.strip()), 0, len(evid_sent.strip()))
                            #print size1 #
                            if checkmatch(evidence_clean.strip(), evid_sent.strip(), 0.9):
                                return ""

            print "Difflib"
            return evidence_clean[:size]
    punct_PMID_list = ["18504304", "19274086", "20190815"]
    if pmid not in punct_PMID_list and len(child_clean.strip()) > 0 and child_clean[-1] in punctuation and child_clean[-1] == evidence_clean[0]:
        print "Punct"
        return child_clean[-1]
    common_str = longest_common_substring(evidence_clean, child_clean)
    if len(common_str.strip()) > 1 and evidence_clean.startswith(common_str):
        print "LCS"
        return common_str
    return ""


def findmatch_NOSPACE(evidence, child):
    # child will have lots of whitespace characters within it. evidence will usually be a well-formed, well-formatted string
    evidence_clean = " ".join(evidence.strip(' \t\n\r').split())
    #print "Clean Evidence: ", evidence_clean
    evidence_nospace = evidence_clean.replace(" ", "")
    child_clean = " ".join(child.strip(' \t\n\r').split())
    child_nospace = child_clean.replace(" ", "")

    s = difflib.SequenceMatcher(None, evidence_nospace, child_nospace)
    pos_a, pos_b, size = s.find_longest_match(0, len(evidence_nospace), 0, len(child_nospace))
    if pos_a == 0 and size > 0:
        mymatch_nospace = evidence_nospace[:size]
        #print "Match No Space: ", mymatch_nospace
        mymatch = ""; j = 0
        for i in range(len(mymatch_nospace)):
            while evidence_clean[j].isspace():
                mymatch += evidence_clean[j]
                j += 1
            if evidence_clean[j] == mymatch_nospace[i]:
                mymatch += evidence_clean[j]
                j += 1
        #print "My shitty match: ", mymatch
        return mymatch
    return ""


def getspanid(pmid):
    global spanctr
    spanctr += 1
    #print spanctr
    return pmid + "_" + str(spanctr)


def populatespanattr(pmid, span_tag):
    span_tag["id"] = getspanid(pmid)
    span_tag["data-highlighted"] = "true"
    # span_tag["data-timestamp"] = "1468083232995"
    span_tag["class"] = "highlighted"
    span_tag["style"] = "background-color: rgb(255, 255, 123);"
    return span_tag


def parse_html(origpsg, soup, prevhtmlpathidx):
    output_list = []; psg = origpsg
    psg_list = get_sentences(psg)

    if len(psg_list[0].split()) < 3 and psg_list[0].endswith("."):
        psg_list = psg_list[1:]

    found_list = [0] * len(psg_list)
    body = soup.body
    ptag_list = body.find_all(p_or_h3) # All p tags in the article that have ids
    #print ptag_list; exit()

    while len(psg_list) > 0:
        for idx, p_tag in enumerate(ptag_list):
            if prevhtmlpathidx != -1:
                if idx < prevhtmlpathidx: continue

            if len(psg_list) == 0:  break

            orig_ptag_sentlist = get_sentences(p_tag.get_text())
            ptag_sentlist = copy.copy(orig_ptag_sentlist)
            accounted_taglist = [0] * len(ptag_sentlist)
            #print "Original Passage: ", origpsg
            #print "Trying to find: ", psg_list
            for i, psent in enumerate(psg_list):
                for j, tagsent in enumerate(ptag_sentlist):
                    #print "Saswati_Divya: "; print psent; print tagsent; print checkmatch(psent, tagsent, 0.9)
                    if checkmatch(psent, tagsent, 0.9) and not accounted_taglist[j]:
                        found_list[i] = 1
                        accounted_taglist[j] = 1
                        break
                if found_list[i] == 1:
                    if j < len(ptag_sentlist)-1:
                        ptag_sentlist = ptag_sentlist[j+1:]
                        accounted_taglist = [0] * len(ptag_sentlist)
                    else:
                        ptag_sentlist = []; accounted_taglist = []
                if i == 0 and found_list[i] == 0:
                    break

            if 1 in found_list:
                if not check_foundlist(found_list):
                    print "Some sentences could not be found in the middle of psg"
                    print psg_list; print found_list; print orig_ptag_sentlist; exit()
                index_0 = getIndex(found_list, 0)   # index_0 contains all indices that are 0 in found_list
                first_sent_found = psg_list[0]
                if index_0 == -1:
                    last_sent_found = psg_list[-1]
                    psg_list = []; found_list = []
                else:
                    last_sent_found = psg_list[index_0[0]-1]    # index_0 is a list
                    psg_list = psg_list[index_0[0]:]
                    found_list = [0] * len(psg_list)
                #print "First Sent: ", first_sent_found
                #print "Last Sent: ", last_sent_found
                #print "Current Passage: ", psg
                found_psg = extract_found_psg(psg, first_sent_found, last_sent_found)
                psg = psg[psg.find(last_sent_found) + len(last_sent_found):]
                #print "Found Passage: ", found_psg
                temp_dict = {}
                temp_dict["htmlpath"] = get_htmlpath(p_tag)
                temp_dict["htmlpathindex"] = idx
                temp_dict["evidence"] = found_psg
                #temp_dict["highlight_details"] = get_details(p_tag, found_psg)
                output_list.append(temp_dict)
                #pprint(temp_dict)

            #targetnode = body.select(mypath)    # Returns a list of nodes. Eg Usage: display_strings(targetnode[0])
        #pprint(output_list)
        #exit()
        if len(found_list) > 0 and found_list[0] == 0:  # I cannot locate the first sentence of psg_list even after going over all the tags in the article
            temp_dict = {}
            temp_dict["htmlpath"] = ""
            temp_dict["htmlpathindex"] = -1
            temp_dict["evidence"] = psg_list[0]
            #print "Cant find Passage: ", psg_list[0]
            output_list.append(temp_dict)
            if len(psg_list) > 1:
                psg = psg[psg.find(psg_list[0]) + len(psg_list[0]):]
                psg_list = psg_list[1:]
                found_list = [0] * len(psg_list)
                #print "Current Passage: ", psg
            else:
                psg_list = []
                found_list = []

    return output_list


def p_or_h3(tag):
    return (tag.name == "p" and tag.has_attr('id')) or (tag.name == "h3" and tag["class"][0] == "section-title")


def get_sentences(mytext):
    mytext = " ".join(mytext.strip(' \t\n\r').split())
    sentences = sent_detector.tokenize(mytext.strip())
    #print "Inside get_sentences: ", sentences
    skip_idx = []; combined_sentences = []
    for i, sent in enumerate(sentences):
        if i in skip_idx: continue
        tempsent = sent; tempidx = i
        inside_while = False
        while tempsent.endswith("Fig."):
            inside_while = True
            if tempidx < (len(sentences)-1) and reg_expr.match(sentences[tempidx + 1].strip()):
                skip_idx.append(tempidx + 1)
                tempsent += " " + sentences[tempidx + 1]
                tempidx += 1
            else:
                break
        if inside_while:
            combined_sentences.append(tempsent)
        else:
            combined_sentences.append(sent)
    #out_sentences = cleanse_punct_shortsent(combined_sentences)
    out_sentences = combined_sentences
    return out_sentences


def cleanse_punct_shortsent(mylist):
    out_list = []
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)
    for sent in mylist:
        num_words = 0
        if len(sent.strip()) == 0: continue
        for word in sent.split():
            out = word.translate(remove_punctuation_map)
            if len(out.strip()) > 1:    # If the length of a word is atleast 2 characters then we consider it as a valid word
                num_words += 1
        if num_words > 3:       # If the sentence has atleast 4 words then we consider it as a valid sentence. False positives for strings like - "Luke et al"
            out_list.append(sent)
    return out_list


def extract_found_psg(psg, first_sent, last_sent):
    if psg.find(first_sent) != -1 and psg.find(last_sent) != -1:
        return psg[psg.find(first_sent):psg.find(last_sent) + len(last_sent)]
    else:
        print "Passage: ", psg.encode("utf-8")
        print "First sentence: ", first_sent.encode("utf-8"), "------------>", psg.find(first_sent)
        print "Last sentence: ", last_sent.encode("utf-8"), "------------>", psg.find(last_sent)
        exit()


def getBeautifulSoup(download_dir, pmcid):
    fname_string = '%s' % (os.path.join(download_dir, pmcid.strip(), pmcid.strip() + ".html"))
    fname = glob.glob(fname_string)
    if len(fname) == 0:
        print "The HTML file was not found for PMCID: ", pmcid.strip()
        exit()
    myfile = fname[0].strip()
    soup = BeautifulSoup(open(myfile), "lxml")
    return soup


def getIndex(mylist, myvalue):
    otpt = []
    for i,val in enumerate(mylist):
        if val == myvalue:
            otpt.append(i)
    if len(otpt) == 0: return -1
    else:   return otpt


def check_foundlist(mylist):
    list_1 = getIndex(mylist, 1)
    for i, elem in enumerate(list_1):
        if i > 0:
            if elem - prev > 1: return False
        prev = elem
    return True


def get_element(node):
    # for XPATH we have to count only for nodes with same type!
    length=0
    for sibling in node.previous_siblings:
        if sibling.name == node.name:
            length+=1
    length+=1
    if (length) > 1:
        return '%s:%s' % (node.name, length)
    else:
        return node.name


def get_htmlpath(node):
    path = [get_element(node)]
    for parent in node.parents:
        if parent.name == 'body':
            break
        path.insert(0, get_element(parent))
    return ' > '.join(path)


def display_strings(node):
    for para_str in node.stripped_strings:
        for sent in para_str.split("."):
            print(repr(sent))


def longest_common_substring(s1, s2):
    m = [[0] * (1 + len(s2)) for i in xrange(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in xrange(1, 1 + len(s1)):
        for y in xrange(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]


def checkmatch(strI, strJ, percent):
    strI = strI.strip(); strJ = strJ.strip()
    s = difflib.SequenceMatcher(None,strI, strJ)
    pos_a, pos_b, size = s.find_longest_match(0, len(strI), 0, len(strJ))
    # print size, len(strI); print strI.encode("utf-8"); print strJ.encode("utf-8")
    if pos_a == 0 and pos_b == 0 and size == 0: return False
    if float(size) >= percent * float(len(strI)):   return True
    else:
        strI_r = convert_str(strI); strJ_r = convert_str(strJ)
        s1 = difflib.SequenceMatcher(None,strI_r, strJ_r)
        pos_a, pos_b, size = s1.find_longest_match(0, len(strI_r), 0, len(strJ_r))
        if pos_a == 0 and pos_b == 0 and size == 0: return False
        if float(size) >= percent * float(len(strI_r)): return True
        else:
            num_diff = 0
            for i, s in enumerate(difflib.ndiff(strI_r, strJ_r)):
                if s[0] == ' ':
                    continue
                elif s[0] == '-' or s[0] == '+':
                    num_diff += 1
            if float(num_diff) <= (1.0-percent) * float(len(strI_r)):   return True
            else:   return False


def convert_str(mystr):
    mystr = mystr.replace(" ", "")  # Remove blank spaces from string
    # mystr = mystr.encode('ascii', 'ignore')     # Remove unicode characters from string
    # mystr = mystr.translate(string.maketrans("", ""), string.punctuation)   # Remove punctuation from string
    mystr = mystr.lower()   # Convert string to lower case
    return mystr


def find_hooks(mygroup, psg):
    group_list = get_sentences(mygroup)
    #group_list = mygroup.split(".")
    #group_list = [x for x in group_list if len(x.strip()) > 0]
    #print "Group List: ", group_list, "\n"
    psg_list = get_sentences(psg)
    #psg_list = psg.split(".")
    #psg_list = [x for x in psg_list if len(x.strip()) > 0]
    #print "Original Passage: ", psg.encode("utf-8")
    #print "Psg List: ", psg_list, "\n"

    group_hooks = []
    for k, gsent in enumerate(group_list):
        if psg_list[0].strip() == gsent.strip():    # The first sentence of psg may occur in its entirety at multiple locations within mygroup
            group_hooks.append(k)
    #print "Group Hooks: ", group_hooks
    if len(group_hooks) == 0:
        common = find_common(mygroup, psg, group_list, psg_list)
    else:
        while len(group_hooks) > 0:
            temp_group_list = group_list[group_hooks[0]:]
            common = find_common(mygroup, psg, temp_group_list, psg_list)
            #print "Current hook: ", group_hooks[0]
            if common.strip() == "" and len(group_hooks) > 0:
                group_hooks.pop(0)
            else:
                break
    return common


def find_common(mygroup, psg, group_list, psg_list):
    mymap = []; partial_map = []
    for j, psent in enumerate(psg_list):
        if len(mymap) == 0:
            if j > 0:   return ""
            for i, gsent in enumerate(group_list):
                if psent.strip() == gsent.strip():
                    mymap.append((i, j))
                    break
                else:
                    commstr = longest_common_substring(gsent.strip(), psent.strip())
                    if commstr.strip() != "" and gsent.endswith(commstr) and psent.startswith(commstr):
                        mymap.append((i, j)); partial_map.append((i,j))
                        break
        else:
            (g,p) = mymap[-1]
            if g == len(group_list)-1:   # g is the last sentence of mygroup
                break
            if psent.strip() == group_list[g+1].strip():
                mymap.append((g+1, j))
            else:
                commstr = longest_common_substring(group_list[g+1].strip(), psent.strip())
                if group_list[g+1].endswith(commstr) and psent.startswith(commstr) and len(mymap) == j:
                    mymap.append((g+1, j)); partial_map.append((g+1,j))
                else:
                    break
    #print mymap
    (g0, p0) = mymap[0]; (gE, pE) = mymap[-1]
    if p0 == 0 and gE == len(group_list)-1:
        if pE == len(psg_list) - 1:
            return mygroup
        for g,p in mymap:
            if len(partial_map) > 0 and (g,p) == partial_map[-1] and (g,p) == mymap[-1]:
                end = psg.find(group_list[g]) + len(group_list[g])
                part_sent = " " + psg_list[p][end:]
                psg_list.insert(pE+1, part_sent)
            else:
                start = psg.find(psg_list[p])   # start may not be 0 always, because of . chars in psg
                psg = psg[start + len(psg_list[p]):]
        return mygroup + " " + psg[psg.find(psg_list[pE+1]):]     # TODO: I might be missing a space between the two strings being merged
    else:
        return ""


def getchildpathindex(child, referrence):
    childindexlist = []
    currchild = child
    for parent in child.parents:
        for j, c in enumerate(parent.children):
            if c == currchild:
                childindexlist.insert(0, str(j))
        currchild = parent
        if parent == referrence:
            break
    return childindexlist


def getchildindex(parent, child):
    for i, c in enumerate(parent.children):
        if c == child:
            return i


def partitionChild(pmid, childstr, found_psg):
    # Assumption: found_psg is a contiguous block of string
    found_tokens = [x.strip() for x in found_psg.split() if len(x.strip()) > 0]
    child_tokens = [x.strip() for x in childstr.split() if len(x.strip()) > 0]

    problem_pmids = ["19274086"]
    if pmid in problem_pmids and found_tokens[0].strip() == "),":
        found_tokens = found_tokens[1:]

    first_token_idx = getIndex(child_tokens, found_tokens[0])
    if first_token_idx == -1:
        print "Something went wrong! Cant find first token of found passage: ", found_psg.encode("utf-8"), " in childstr: ", childstr.encode("utf-8")
        exit()

    childstr_copy = childstr
    for i in range(len(first_token_idx)):
        cidx = childstr_copy.find(found_tokens[0])
        mismatch = False; mymatch = ""
        for j in range(len(found_tokens)):
            for k in range(len(found_tokens[j])):
                while childstr_copy[cidx].isspace():
                    mymatch += childstr_copy[cidx]
                    cidx += 1
                if childstr_copy[cidx] == found_tokens[j][k]:
                    mymatch += childstr_copy[cidx]
                    cidx += 1
                else:
                    mismatch = True
                    break
            if mismatch:    break
        if mismatch:
            childstr_copy = childstr_copy[cidx:]
        else:
            break

    olist = []
    aftermatch = childstr[childstr.find(mymatch) + len(mymatch):]
    if len(aftermatch.strip()) == 0:
        mymatch += aftermatch
        aftermatch = ""
    beforematch = childstr[:childstr.find(mymatch)]
    if len(beforematch.strip()) == 0:
        mymatch = beforematch + mymatch
        beforematch = ""
    mylist = [beforematch, mymatch, aftermatch]
    for elem in mylist:
        if len(elem) == 0:
            continue
        else:
            if elem == mymatch:
                olist.append((elem, 1))
            else:
                olist.append((elem, 0))
    return olist


if __name__ == "__main__":
    main_body()
