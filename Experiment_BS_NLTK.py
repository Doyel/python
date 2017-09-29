# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
from bs4 import BeautifulSoup;
import nltk.data;
from pprint import pprint;
import os;
import glob;
import copy
import cPickle as pickle;

# python Experiment_BS_NLTK.py ./static/articles ./Test_PMCID.txt
# python Experiment_BS_NLTK.py ./downloaded_articles ./Test_PMCID.txt

def main_body():
    global sent_detector
    parser = argparse.ArgumentParser(prog='Use BeautifulSoup NLTK', usage='Experiment_BS_NLTK.py <downloaded_articles_dir> <filelist>', description='Script to experiment with HTML files')    
    parser.add_argument('download_dir', help='Dir where the downloaded articles are located')
    parser.add_argument('JSONfilelist', help='File listing the PMC ids of articles')
    args = parser.parse_args()

    sent_detector = nltk.data.load('tokenizers/punkt/english.pickle')

    with open(args.JSONfilelist, 'rt') as f1:
        for pmcid in f1:
            if not pmcid.strip(): continue
            fname_string = '%s' % (os.path.join(args.download_dir, "PMC"+pmcid.strip(), "PMC"+pmcid.strip()+".html"))
            fname = glob.glob(fname_string)
            if len(fname) == 0:
                print "The HTML file was not found for PMCID: ", pmcid.strip()
                continue    
            myfile = fname[0].strip()

            soup1 = BeautifulSoup(open(myfile), "lxml")

            #tag_list1 = soup1.body.find_all(p_or_h3)

            # for p_tag in body.find_all("p", id=True):
            #     mypath = get_xpath(p_tag); print mypath
            #     display_strings(p_tag)
            #
            #     targetnode = body.select(mypath) #getnodefromxpath(mypath, copy.copy(body))
            #     display_strings(targetnode[0])

            p_path = "3:1:13"
            p_node = getnodefromidxpath(p_path, soup1.body)
            print p_node
            for i, child in enumerate(p_node.children):
                print "Index No: ", i
                print child


def p_or_h3(tag):
    return (tag.name == "p" and tag.has_attr('id')) or (tag.name == "h3" and tag["class"][0] == "section-title")


def get_element(node):
    # for XPATH we have to count only for nodes with same type!
    length=0
    for sibling in node.previous_siblings:
        if sibling.name == node.name:
            length+=1
    length+=1
    if (length) > 1:
        return '%s:nth-of-type(%s)' % (node.name, length)
    else:
        return node.name


def get_xpath(node):
    path = [get_element(node)]
    for parent in node.parents:        
        if parent.name == 'body':
            break
        path.insert(0, get_element(parent))
    return ' > '.join(path)


def getnodefromxpath(xpath, prevnode):
    nodes = xpath.split(' > ')
    while len(nodes) > 0:
        currnodestr = nodes.pop(0)
        if ':' in currnodestr:                        
            childno = currnodestr.split(':')[1]
            currnodestr = currnodestr.split(':')[0]
        else:
            childno = '1'
        currno=1
        for child in prevnode.children:
            if child.name == currnodestr:
                if currno == int(childno):
                    prevnode = copy.copy(child)
                    break
                else:
                    currno+=1
    return prevnode


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


def display_strings(node):    
    for para_str in node.strings:
        for sent in sent_detector.tokenize(para_str.strip()):
            print(repr(sent))


if __name__ == "__main__":
    main_body()
