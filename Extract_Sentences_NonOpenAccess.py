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

# Save this script at the location: /ua/ml-group/big-mechanism-project/PLEIADES/workDir/Eric_SRI_Perl
# python Extract_Sentences_NonOpenAccess.py pubmedids_EricSRI_WWWIndicator.txt ./downloaded_articles ./JSON_SENTENCE_ASSOCS_SRI

def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Extract_Sentences',
                                     usage='Extract_Sentences_NonOpenAccess.py <PubMedfilelist> <downloaded_articles_dir> <write_dir>',
                                     description='Script to extract sentences from HTML articles and write the output into a JSON file having the same format as NonOpenAccess articles')
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
    # txt_location = "./"  # Temporary: Please comment out after Eric SRI deliverables are sent

    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f:
        for line in f:
            pmid = line.split("|")[0].strip()
            if len(line.split("|")) == 2:
                download_web = bool(line.split("|")[1].strip())
            else:   # Default: Article is OpenAccess and was downloaded in NXML format before being converted to html
                print "Treating current article as OpenAccess since the identifier for indicating whether it was downloaded from web is absent!"
                download_web = False

            if not pmid or pmid.startswith("#"): continue

            if pm_pmc[pmid] != "Does Not Exist":
                print "PMID: ", pmid, "\tPMCID: ", pm_pmc[pmid]
                soup = getBeautifulSoup(args.downloaded_articles_dir, pm_pmc[pmid])
            else:
                print "PMID: ", pmid, "\tDoes not have PMCID"
                soup = getBeautifulSoup(args.downloaded_articles_dir, "PubMed"+pmid)

            if download_web:
                tag_list = soup.body.find_all(p_or_h3_web)
            else:
                tag_list = soup.body.find_all(p_or_h3)

            sentences = {}; abs_sent_ctr = 0
            sentences["@type"] = "com.sri.projects.bigmech.jjnpipeline.BigmechPipeline$DocumentInfo"
            sentences["pmid"] = pmid.strip()
            sentences["sentenceList"] = {}
            sentences["sentenceList"]["@type"] = "java.util.ArrayList"
            sentences["sentenceList"]["@items"] = []
            # p = re.compile('^[S]?[0-9][A-Za-z]?\W*')

            for tag in tag_list:
                tag_text = tag.get_text()
                cleansed_tag_text = " ".join(tag_text.strip(' \t\n\r').split())
                if not cleansed_tag_text:
                    continue

                htmlIdxPath = getchildpathindex(tag, soup.body)
                checknode = getnodefromidxpath(htmlIdxPath, soup.body)
                if checknode.get_text() != tag_text:
                    print "Idx path points to wrong node. Please check! \n", checknode.get_text().encode("utf-8"), "\n\n\n", tag_text.encode("utf-8")
                    exit()

                span_sentences = tokenizer.span_tokenize(tag_text)

                for (start, end) in span_sentences:
                    temp = dict()
                    temp["@type"] = "com.sri.projects.bigmech.jjnpipeline.BigmechPipeline$SentenceInfo"
                    temp["index"] = abs_sent_ctr
                    temp["sentenceText"] = tag_text[start:end]
                    temp["cleansedText"] = " ".join(temp["sentenceText"].strip(' \t\n\r').split())
                    if not temp["cleansedText"]:
                        print "Control should not come here!"   # I was correct! Control never came here :)
                        continue
                    temp["startCharPos"] = start
                    temp["sentenceLen"] = end - start
                    temp["endCharPos"] = end
                    abs_sent_ctr += 1
                    sentences["sentenceList"]["@items"].append(temp)

            fullname_write_fname = os.path.join(args.write_dir, pmid.strip())
            json_fn = open(fullname_write_fname, 'w')
            json.dump(sentences, json_fn, indent=4, ensure_ascii=False)
            json_fn.close()


def getBeautifulSoup(download_dir, pmcid):
    fname_string = '%s*' % (os.path.join(download_dir, pmcid.strip(), pmcid.strip() + ".htm"))
    fname = glob.glob(fname_string)
    if len(fname) == 0:
        print "The HTML file was not found for PMCID: ", pmcid.strip()
        exit()
    myfile = fname[0].strip()
    soup = BeautifulSoup(open(myfile), "lxml")
    return soup


def p_or_h3(tag):
    return (tag.name == "p" and tag.has_attr('id')) or (tag.name == "h3" and tag["class"][0] == "section-title")


def p_or_h3_web(tag):
    return (tag.name == "p") or (tag.name == "h3")


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
