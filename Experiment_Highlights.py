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
from nltk.tokenize import StanfordTokenizer
import re


def main_body():
    markup = '<a href="http://example.com/">I linked to <i>example.com</i></a>'
    soup = BeautifulSoup(markup, "lxml")
    a_tag = soup.a

    new_tag = soup.new_tag("b")
    new_tag.string = "example.net"
    mystr = a_tag.i.string
    mystr.replace_with("Saswati")
    print a_tag
    #exit()

    evidence = "Cras congue id est sit amet mattis. Sed in metus et orci eleifend commodo. Phasellus at odio imperdiet, 21, 22 efficitur augue in, pulvinar sapien. Pellentesque leo nulla, porta non lectus eu, ullamcorper semper est. Nunc convallis risus vel mauris accumsan, in rutrum odio sodales. Vestibulum ante ipsum primis in faucibus orci luctus et ultrices posuere cubilia Curae"
    #evidence = "Nam felis velit, ullamcorper eu turpis ut, hendrerit accumsan augue. Nulla et purus sem. Ut at hendrerit purus. Phasellus mollis commodo ante eu mollis."
    soup = BeautifulSoup(open("mir3z.html"), "lxml")
    sandbox = soup.find(id="sandbox")   # referrence node

    first_p = sandbox.p     # target node

    for i, child in enumerate(sandbox.children):
        if child == first_p: break

    print first_p
    first_p = updatenode(soup, first_p, evidence)
    print first_p
            #pprint(lelem[0])
    #print first_p

    print "Saswati <3 Divya"
    for i, child in enumerate(first_p.descendants):
        if isinstance(child, Tag) and child.name == "span":
            childlen = len(child.string)
            childindexlist = getchildpathindex(child, sandbox)
            if isinstance(child.previous_sibling, NavigableString):
                offset = len(child.previous_sibling)
            else:
                offset = 0

            print child, childlen, offset, ":".join(childindexlist)

    exit()



    my_new_tag = soup.new_tag(first_p.name)
    for i, child in enumerate(first_p.children):
        if len(evidence.strip()) == 0: break
        if isinstance(child, NavigableString):
            childstr = child
        else:
            childstr = child.string     # child is a tag

        found_psg = findmatch(evidence, childstr)   # found_psg is a clean string
        if found_psg == "":
            if isinstance(child, NavigableString):
                navg_strg = NavigableString(childstr)
                my_new_tag.append(navg_strg)
            else:
                #new_child_tag = soup.new_tag(child.name)
                #navg_strg = NavigableString(childstr)
                #new_child_tag.append(navg_strg)
                new_child_tag = copy.copy(child)
                my_new_tag.append(new_child_tag)
        else:
            evidence = evidence[evidence.find(found_psg) + len(found_psg) + 1:]    # found_psg will always be at the beginning of evidence
            parts = partitionChild(childstr, found_psg)
            if isinstance(child, NavigableString):
                for strg, switch in parts:
                    navg_strg = NavigableString(strg)
                    if switch:
                        span_tag = soup.new_tag("span")
                        span_tag.append(navg_strg)
                        my_new_tag.append(span_tag)
                    else:
                        my_new_tag.append(navg_strg)
            else:
                new_child_tag = soup.new_tag(child.name)
                for strg, switch in parts:
                    navg_strg = NavigableString(strg)
                    if switch:
                        span_tag = soup.new_tag("span")
                        span_tag.append(navg_strg)
                        new_child_tag.append(span_tag)
                    else:
                        new_child_tag.append(navg_strg)
                my_new_tag.append(new_child_tag)

    first_p.replace_with(my_new_tag)

    first_p = sandbox.p

    for i, child in enumerate(first_p.descendants):
        print child

    print "Saswati <3 Divya"
    for i, child in enumerate(first_p.descendants):
        if isinstance(child, Tag) and child.name == "span":
            childlen = len(child.string)
            childindexlist = getchildpathindex(child, sandbox)
            if isinstance(child.previous_sibling, NavigableString):
                offset = len(child.previous_sibling)
            else:
                offset = 0

            print child, childlen, offset, ":".join(childindexlist)


def updatenode(soup, mynode, evidence):
    mylist = []
    for child in mynode.descendants:

        if len(evidence.strip()) == 0 or evidence.strip() == ".": break

        if isinstance(child, NavigableString):
            childstr = child
            found_psg = findmatch(evidence, childstr)  # found_psg is a clean string
            if found_psg != "":
                if (evidence.find(found_psg) + len(found_psg) + 1) < len(evidence) and evidence[
                                    evidence.find(found_psg) + len(found_psg) + 1].isspace():
                    evidence = evidence[evidence.find(found_psg) + len(
                        found_psg) + 1:]  # found_psg will always be at the beginning of evidence
                else:
                    evidence = evidence[evidence.find(found_psg) + len(found_psg):]
                curr_parent = child.parent
                child_idx = getchildindex(curr_parent, child)
                parent_path = get_htmlpath(curr_parent)
                mylist.append([parent_path, curr_parent, child, found_psg, child_idx])

    if len(evidence.strip()) > 0 and evidence.strip() != ".":
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
    # print first_p
    for ppath in mydict:
        for lelem in mydict[ppath]:  # mydict[ppath] is a list of lists. lelem is a 4-element list: [curr_parent 0, child 1, found_psg 2, child_idx 3]
            childstr = lelem[1];
            found_psg = lelem[2];
            child_idx = lelem[3]
            # print lelem[0]; print childstr; print child_idx
            parts = partitionChild(childstr, found_psg)
            # print parts; #exit()
            for i, (strg, switch) in enumerate(parts):
                navg_strg = NavigableString(strg)
                if switch:
                    span_tag = soup.new_tag("span")
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


def getchildfromindex(parent, child_idx):
    for i, c in enumerate(parent.children):
        if i == child_idx:
            return c


def locate_parent_in_copy(curr_parent, copy_first_p):
    if curr_parent == copy_first_p:
        return copy_first_p
    else:
        for child in copy_first_p.descendants:
            if child == curr_parent:
                return child


def findmatch(evidence, child):
    # child will have lots of whitespace characters within it. evidence will usually be a well-formed, well-formatted string
    evidence_clean = " ".join(evidence.strip(' \t\n\r').split())
    print "Clean Evidence: ", evidence_clean
    evidence_nospace = evidence_clean.replace(" ", "")
    child_clean = " ".join(child.strip(' \t\n\r').split())
    child_nospace = child_clean.replace(" ", "")

    s = difflib.SequenceMatcher(None, evidence_nospace, child_nospace)
    pos_a, pos_b, size = s.find_longest_match(0, len(evidence_nospace), 0, len(child_nospace))
    if pos_a == 0 and size > 0:
        mymatch_nospace = evidence_nospace[:size]
        print "Match No Space: ", mymatch_nospace
        mymatch = ""; j = 0
        for i in range(len(mymatch_nospace)):
            while evidence_clean[j].isspace():
                mymatch += evidence_clean[j]
                j += 1
            if evidence_clean[j] == mymatch_nospace[i]:
                mymatch += evidence_clean[j]
                j += 1

        print "My shitty match: ", mymatch
        return mymatch
    return ""


def findmatch_old(evidence, child):
    # child will have lots of whitespace characters within it. evidence will usually be a well-formed, well-formatted string
    evidence_clean = " ".join(evidence.strip(' \t\n\r').split())
    child_clean = " ".join(child.strip(' \t\n\r').split())
    s = difflib.SequenceMatcher(None, evidence_clean, child_clean)
    pos_a, pos_b, size = s.find_longest_match(0, len(evidence_clean), 0, len(child_clean))
    if pos_a == 0 and size > 0:
        return evidence_clean[:size]
    return ""


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


def get_element(node):
    # for XPATH we have to count only for nodes with same type!
    length = 0
    for sibling in node.previous_siblings:
        if sibling.name == node.name:
            length += 1
    length += 1
    if (length) > 1:
        return '%s:%s' % (node.name, length)
    else:
        return node.name


def get_htmlpath(node):
    path = [get_element(node)]
    for parent in node.parents:
        if parent.has_attr('id') and parent["id"] == 'sandbox':
            break
        path.insert(0, get_element(parent))
    return ' > '.join(path)


def getIndex(mylist, myvalue):
    otpt = []
    for i,val in enumerate(mylist):
        if val == myvalue:
            otpt.append(i)
    if len(otpt) == 0: return -1
    else:   return otpt


def partitionChild(childstr, found_psg):
    # Assumption: found_psg is a contiguous block of string
    found_tokens = [x.strip() for x in found_psg.split() if len(x.strip()) > 0]
    child_tokens = [x.strip() for x in childstr.split() if len(x.strip()) > 0]

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
