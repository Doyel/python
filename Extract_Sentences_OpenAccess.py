# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
from bs4 import BeautifulSoup;
from nltk.tokenize.punkt import PunktSentenceTokenizer, PunktParameters
from pprint import pprint
import sys
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import glob;
import json
import cPickle as pickle;
import re
from Utilities import read_config_file

# python Extract_Sentences_OpenAccess.py ./PubMedIDS_with_Extras_OpenAccess.txt ./downloaded_articles ./JSON_SENTENCE_ASSOCS_OpenAccess

def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Extract_Sentences', usage='Extract_Sentences.py <OpenAccess_RNAifilelist> <downloaded_articles_dir> <write_dir>', description='Script to extract sentences from HTML articles')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids that are in OpenAccess and are also in John list')
    parser.add_argument('downloaded_articles_dir', help='Path where openaccess html articles were downloaded')
    parser.add_argument('write_dir', help='Path where the file containing the extracted sentences are written into')
    args = parser.parse_args()

    #sent_detector = nltk.data.load('tokenizers/punkt/english.pickle')
    punkt_param = PunktParameters()
    abbreviation = ['fig']
    punkt_param.abbrev_types = set(abbreviation)
    tokenizer = PunktSentenceTokenizer(punkt_param)

    pm_pmc = pickle.load(open("pmid_pmcid_Extras.p", "rb"))
    htmlIdx_paraId = {}; absId_para_sentId = {}

    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f:
        for pmid in f:
            if not pmid.strip() or pmid.startswith("#"): continue
            pmid = pmid.strip()
            print "PMID: ", pmid, "\tPMCID: ", pm_pmc[pmid]

            htmlIdx_paraId[pmid] = {}; absId_para_sentId[pmid] = {}
            soup = getBeautifulSoup(args.downloaded_articles_dir, pm_pmc[pmid])
            tag_list = soup.body.find_all(p_or_h3)
            sentences = {}; para_ctr = 0; abs_sent_ctr = 0
            sentences["pmid"] = pmid.strip()
            sentences["paragraphList"] = []
            # p = re.compile('^[S]?[0-9][A-Za-z]?\W*')

            for tag in tag_list:
                tempdict = {}
                tag_text = tag.get_text()
                cleansed_tag_text = " ".join(tag_text.strip(' \t\n\r').split())
                if not cleansed_tag_text:
                    continue
                tempdict["paragraphId"] = para_ctr
                tempdict["rawText"] = tag_text
                tempdict["paragraphLen"] = len(tag_text)
                tempdict["htmlIdxPath"] = getchildpathindex(tag, soup.body)
                tempdict["sentenceList"] = {}
                tempdict["sentenceList"]["@items"] = []

                checknode = getnodefromidxpath(tempdict["htmlIdxPath"], soup.body)
                if checknode.get_text() != tag_text:
                    print "Idx path points to wrong node. Please check! \n", checknode.get_text().encode("utf-8"), "\n\n\n", tag_text.encode("utf-8")
                    exit()
                if tempdict["htmlIdxPath"] in htmlIdx_paraId[pmid]:
                    print "Two distinct nodes have the html idx path!"
                    exit()
                htmlIdx_paraId[pmid][tempdict["htmlIdxPath"]] = tempdict["paragraphId"]

                span_sentences = tokenizer.span_tokenize(tag_text)

                sent_ctr = 0
                for (start, end) in span_sentences:
                    temp = dict()
                    temp["sentenceId"] = sent_ctr
                    temp["absoluteId"] = abs_sent_ctr
                    temp["sentenceText"] = tag_text[start:end]
                    temp["cleansedText"] = " ".join(temp["sentenceText"].strip(' \t\n\r').split())
                    if not temp["cleansedText"]:
                        print "Control should not come here!"   # I was correct! Control never came here :)
                        continue
                    temp["charOffset"] = start
                    temp["sentenceLen"] = end - start
                    sent_ctr += 1; abs_sent_ctr += 1
                    tempdict["sentenceList"]["@items"].append(temp)

                    if temp["absoluteId"] in absId_para_sentId[pmid]:
                        print "Current absolute id already present in dict. Please check!", "\t", temp["absoluteId"]
                        exit()
                    absId_para_sentId[pmid][temp["absoluteId"]] = {}
                    absId_para_sentId[pmid][temp["absoluteId"]]["sentId"] = temp["sentenceId"]
                    absId_para_sentId[pmid][temp["absoluteId"]]["paraId"] = tempdict["paragraphId"]
                    absId_para_sentId[pmid][temp["absoluteId"]]["htmlIdxPath"] = tempdict["htmlIdxPath"]

                para_ctr += 1
                sentences["paragraphList"].append(tempdict)

            fullname_write_fname = os.path.join(args.write_dir, pmid.strip()+"_OA")
            json_fn = open(fullname_write_fname, 'w')
            json.dump(sentences, json_fn, indent=4, ensure_ascii=False)
            json_fn.close()

    pickle.dump(htmlIdx_paraId, open("htmlIdx_paraId.p", "wb"))
    pickle.dump(absId_para_sentId, open("absId_para_sentId.p", "wb"))


def getBeautifulSoup(download_dir, pmcid):
    fname_string = '%s' % (os.path.join(download_dir, pmcid.strip(), pmcid.strip() + ".html"))
    fname = glob.glob(fname_string)
    if len(fname) == 0:
        print "The HTML file was not found for PMCID: ", pmcid.strip()
        exit()
    myfile = fname[0].strip()
    soup = BeautifulSoup(open(myfile), "lxml")
    return soup


def p_or_h3(tag):
    return (tag.name == "p" and tag.has_attr('id')) or (tag.name == "h3" and tag["class"][0] == "section-title")


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


def getchildpathindex(child, referrence):
    childindexlist = []
    currchild = child
    for parent in child.parents:
        for j, c in enumerate(parent.children):
            if c == currchild:
                childindexlist.insert(0, str(j))
                break
        currchild = parent
        if parent == referrence:
            break
    return ":".join(childindexlist)


if __name__ == "__main__":
    main_body()
