# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import copy;
import xml.etree.ElementTree as ET
import networkx as nx
from networkx.readwrite import json_graph
import string
from pprint import pprint;
import cPickle as pickle;


# python Generate_NetworkX_Dict_OpenAccess.py ./PubMedIDS_with_Extras_OpenAccess.txt ./JSON_SENTENCE_ASSOCS_OpenAccess ./OA_graph_dict.json


def main_body():
    global absId_para_sentId
    parser = argparse.ArgumentParser(prog='Generate_NetworkX_Dict_OpenAccess', usage='Generate_NetworkX_Dict_OpenAccess.py <PubMedfilelist> <sent_dir> <outputFile>', description='Script to generate NetworkX dict')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('Extr_Sent_Dir', help='Directory where the extracted sentences for all articles are located')
    parser.add_argument('outputJSON', help='Output JSON filename that will store the graph, prot and trigger details')

    args = parser.parse_args()
    depParse_folder = "Sentences_DepParse"
    absId_para_sentId = pickle.load(open("absId_para_sentId.p", "rb"))
    mydict = {}
    with open(args.PubMedfilelist, 'rt') as f1:
        for pmid in f1:
            pmid = pmid.strip()
            if not pmid: continue
            sentence_data = getSentencesFile(args.Extr_Sent_Dir, pmid)
            if not sentence_data:   continue    # No sentences were extracted for this article!

            print "PMID: ", pmid
            mydict[pmid] = {}
            root = getXMLTreeRoot(os.path.join(args.Extr_Sent_Dir, depParse_folder), pmid + "_*xml")
            sentences = root.find('./document/sentences')
            for sent in sentences:
                sentid = checkforMultipleSent(sent)
                extracted_sent = getSentence(pmid, sentence_data, sentid)
                mydict[pmid][sentid] = {}
                for child in sent:
                    if child.tag == "tokens":
                        mydict[pmid][sentid]["Tokens"] = getTokensDict(child)
                        tokens = getTokens(mydict[pmid][sentid]["Tokens"])
                        copied_sent = copy.deepcopy(extracted_sent)
                        if "ProteinDetails" in extracted_sent:
                            mydict[pmid][sentid]["Proteins"] = getProteinDict(copied_sent, tokens)
                        if "TriggerDetails" in extracted_sent:
                            mydict[pmid][sentid]["Triggers"] = getTriggerDict(copied_sent, tokens)

                    if child.get("type") == "enhanced-plus-plus-dependencies":
                        G = get_NetworkX_DiGraph(child, mydict[pmid][sentid]["Tokens"])
                        data = json_graph.node_link_data(G)
                        mydict[pmid][sentid]["DiGraph"] = json.dumps(data)
                        if "Triggers" in mydict[pmid][sentid]:
                            U = get_NetworkX_UnGraph(child)
                            un_data = json_graph.node_link_data(U)
                            mydict[pmid][sentid]["UnGraph"] = json.dumps(un_data)
    json_fn = open(args.outputJSON, 'w')
    json.dump(mydict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


def getProteinDict(extracted_sent, stanford_tokens):
    mydict = dict()
    for prot_dict in extracted_sent["ProteinDetails"]:  # There may be multiple proteins detected in a single sentence
        if prot_dict["protName"] in mydict:
            print "For sentence id: ", extracted_sent["absoluteId"]
            print "The same protein already exists in the dictionary: ", prot_dict["protName"]
            exit()
        mydict[prot_dict["protName"]] = getProtTokenids(prot_dict, extracted_sent, stanford_tokens)
    return mydict


def getProtTokenids(prot_dict, sent_dict, tokens):
    mylist = list()
    mysent = cleanse_quotes(sent_dict["sentenceText"])
    if check_for_unicode(mysent):
        mysent = mysent.decode("utf-8")
    tokens_spans = get_span_tokens(mysent, tokens, sent_dict["charOffset"])
    # The same protein may be detected multiple times in the same sentence, hence the for loop below
    for offsets in prot_dict["protcharOffset"]: # offsets is a list with 3 elements: [startpos, endpos, highlightColor]
        tokenids = []   # The protein spanned by 'offsets' could be a multi-token protein, hence 'tokenids' is a list!
        for i, tok in enumerate(tokens_spans):    # tok is a tuple having the format: (word, (startpos, endpos))
            if offsets[0] >= tok[1][0]:
                if offsets[0] > tok[1][1]:
                    continue
                elif offsets[0] < tok[1][1]:
                    tokenids.append(i + 1)
                    if offsets[1] > tok[1][1]:
                        offsets[0] = tokens_spans[i + 1][1][0]
            elif offsets[0] < tok[1][0]:
                break
        if not tokenids:
            print "We couldnt find token ids for the protein: ", prot_dict["protName"]
            print "Sentence id: ", sent_dict["absoluteId"]
            print "Sentence: ", sent_dict["sentenceText"]
            exit()
        mylist.append(tokenids)
    return mylist


def getTriggerDict(extracted_sent, stanford_tokens):
    mydict = dict()
    for trig_dict in extracted_sent["TriggerDetails"]:  # There may be multiple proteins detected in a single sentence
        if trig_dict["label"] not in mydict:
            mydict[trig_dict["label"]] = {}
        if trig_dict["triggerWord"] in mydict[trig_dict["label"]]:
            print "For sentence id: ", extracted_sent["absoluteId"], " and label: ", trig_dict["label"]
            print "The same trigger word already exists in the dictionary: ", trig_dict["triggerWord"]
            exit()
        mydict[trig_dict["label"]][trig_dict["triggerWord"]] = getTriggerTokenids(trig_dict, extracted_sent, stanford_tokens)
    return mydict


def getTriggerTokenids(trig_dict, sent_dict, tokens):
    mylist = list()
    mysent = cleanse_quotes(sent_dict["sentenceText"])
    if check_for_unicode(mysent):
        mysent = mysent.decode("utf-8")
    tokens_spans = get_span_tokens(mysent, tokens, sent_dict["charOffset"])
    # The same trigger may be detected multiple times in the same sentence for a given label, hence the for loop below
    for offsets in trig_dict["triggerCharOffset"]: # offsets is a list with 2 elements: [startpos, endpos]
        tokenids = []   # Current trigger words are only single-token triggers
        for i, tok in enumerate(tokens_spans):    # tok is a tuple having the format: (word, (startpos, endpos))
            if offsets[0] >= tok[1][0]:
                if offsets[0] > tok[1][1]:
                    continue
                elif offsets[0] < tok[1][1]:
                    tokenids.append(i + 1)
                    if offsets[1] > tok[1][1]:
                        offsets[0] = tokens_spans[i + 1][1][0]
            elif offsets[0] < tok[1][0]:
                break
        if not tokenids:
            print "We couldnt find token ids for the trigger word: ", trig_dict["triggerWord"]
            print "Sentence id: ", sent_dict["absoluteId"]
            print "Sentence: ", sent_dict["sentenceText"]
            exit()
        mylist.append(tokenids)
    return mylist


def get_span_tokens(sent, tokens, startChar):
    span_tokens = []; offset = 0
    for token in tokens:
        if token == "``" or token == "''":
            token = "\""
        offset = sent.find(token, offset)
        span_tokens.append((token, (startChar + offset, startChar + offset + len(token))))
        assert token == sent[offset:offset + len(token)]
        offset += len(token)
    return span_tokens


def cleanse_quotes(mystr):
    if "``" in mystr:
        mystr = mystr.replace("``", "\"")
    if "''" in mystr:
        mystr = mystr.replace("''", "\"")
    return mystr


def check_for_unicode(mystring):
    try:
        mystring.decode('ascii')
    except UnicodeDecodeError:
        #print mystring
        return True
    else:
        return False


def getSentencesFile(Extr_Sent_Dir, pmid):
    sent_file = get_file(Extr_Sent_Dir, pmid+"_*")
    file_handle = open(sent_file, "rb")
    sentfile_data = json.load(file_handle)  # sentfile_data is a dictionary
    file_handle.close()
    if "paragraphList" in sentfile_data and len(sentfile_data["paragraphList"]) > 0:
        return sentfile_data    # There are extracted sentences for this article!
    else:
        return {}


def getSentence(pmid, sentence_data, absId):
    paraId = absId_para_sentId[pmid][absId]["paraId"]
    sentId = absId_para_sentId[pmid][absId]["sentId"]
    mysent = sentence_data["paragraphList"][paraId]["sentenceList"]["@items"][sentId]
    if mysent["absoluteId"] != absId:
        print "We got the wrong sentence: \n", mysent["cleansedText"]
        print "Sentence Id obtained from the xml file: ", absId + 1
        exit()
    return mysent


def getXMLTreeRoot(filepath, pattern):
    xml_file = get_file(filepath, pattern)
    xml_tree = ET.parse(xml_file)
    root = xml_tree.getroot()
    return root


def compare_tokens(edge_tokens_dict, xml_tokens_dict):
    edge_num_tokens = edge_tokens_dict.keys()
    xml_num_tokens = xml_tokens_dict.keys()
    if set(edge_num_tokens) != set(xml_num_tokens):
        print "The token ids obtained from the dependency nodes are different from the token ids obtained from the \"tokens\" xml element node"
        print "Token ids obtained from the dependency node elements: "; pprint(edge_tokens_dict)
        print "Token ids obtained from the \"tokens\" node element: "; pprint(xml_tokens_dict)
        exit()
    for tokenid in edge_num_tokens:
        if edge_tokens_dict[tokenid] != xml_tokens_dict[tokenid]["word"]:
            print "The token ids do not point to the same word: ", edge_tokens_dict[tokenid], "\t", xml_tokens_dict[tokenid]["word"]
            exit()


def getTokens(tokens_dict):
    mylist = []
    for tk_id in tokens_dict:
        mylist.append(tokens_dict[tk_id]["word"])
    return mylist


def getTokensDict(xml_element):
    tokens_dict = dict()
    for token in xml_element:
        tokenid = int(token.get("id"))
        tokens_dict[tokenid] = {}
        tokens_dict[tokenid]["word"] = token.find("word").text
        tokens_dict[tokenid]["pos"] = token.find("POS").text
    return tokens_dict


def get_NetworkX_DiGraph(xml_element, xml_tokens_dict):
    tokens_dict = dict(); G = nx.DiGraph()
    for dep in xml_element:
        edge_name = dep.get("type")
        governor_word = dep.find("governor").text
        governor_tokenid = int(dep.find("governor").get("idx"))
        dependent_word = dep.find("dependent").text
        dependent_tokenid = int(dep.find("dependent").get("idx"))
        if governor_tokenid > 0 and governor_tokenid not in tokens_dict:
            tokens_dict[governor_tokenid] = governor_word
        if dependent_tokenid > 0 and dependent_tokenid not in tokens_dict:
            tokens_dict[dependent_tokenid] = dependent_word
        G.add_edge(governor_tokenid, dependent_tokenid, name=edge_name, weight=1)
        G.node[governor_tokenid]['word'] = governor_word
        G.node[dependent_tokenid]['word'] = dependent_word
    compare_tokens(tokens_dict, xml_tokens_dict)
    return G


def get_NetworkX_UnGraph(xml_element):
    G = nx.Graph()
    for dep in xml_element:
        edge_name = dep.get("type")
        governor_word = dep.find("governor").text
        governor_tokenid = int(dep.find("governor").get("idx"))
        dependent_word = dep.find("dependent").text
        dependent_tokenid = int(dep.find("dependent").get("idx"))
        G.add_edge(governor_tokenid, dependent_tokenid, name=edge_name, weight=1)
        G.node[governor_tokenid]['word'] = governor_word
        G.node[dependent_tokenid]['word'] = dependent_word
    return G


def get_file(file_dir, pattern):
    myfilelist = glob.glob(os.path.join(file_dir, pattern))
    if len(myfilelist) > 1:
        print "Multiple files were matched by glob! \n", myfilelist
        exit()
    return myfilelist[0]


def checkforMultipleSent(sent):
    sentid = int(sent.get("id"))
    lineid = int(sent.get("line"))
    if lineid != sentid:
        print "Multiple lines were detected in a single sentence: ", sentid
        exit()
    return sentid-1 # Sent ids in the extracted sentences file start from 0 whereas they start from 1 in the xml files

if __name__ == "__main__":
    main_body()