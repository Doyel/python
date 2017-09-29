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
import unicodedata


# Total_PubMedIds_John_NoOpenAccess.txt does not contain the 61 OpenAccess articles because the format of the extracted sentences files for the
# OpenAccess articles are different from the format of the extracted sentences files for the other articles sent by John
# python Generate_NetworkX_Dict_PdfArticles.py ./Total_PubMedIds_John_NoOpenAccess.txt ./JSON_SENTENCE_ASSOCS ./NonOA_graph_dict.json


def main_body():
    global sentid; global pmid; global problem_pmids
    parser = argparse.ArgumentParser(prog='Generate_NetworkX_Dict_PdfArticles', usage='Generate_NetworkX_Dict_PdfArticles.py <PubMedfilelist> <sent_dir> <outputFile>', description='Script to generate NetworkX dict')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('Extr_Sent_Dir', help='Directory where the extracted sentences for all articles are located')
    parser.add_argument('outputJSON', help='Output JSON filename that will store the graph, prot and trigger details')

    args = parser.parse_args()
    depParse_folder = "Sentences_DepParse"

    problem_pmids = []  #["10023673"]

    mydict = {}
    with open(args.PubMedfilelist, 'rt') as f1:
        for pmid in f1:
            pmid = pmid.strip()
            #pmid = "8916957"
            if not pmid or pmid in problem_pmids: continue
            sentence_data = getSentencesFile(args.Extr_Sent_Dir, pmid)
            if not sentence_data:   continue    # No sentences were extracted for this article!

            print "PMID: ", pmid
            mydict[pmid] = {}
            root = getXMLTreeRoot(os.path.join(args.Extr_Sent_Dir, depParse_folder), pmid + "_*xml")
            sentences = root.find('./document/sentences')
            for sent in sentences:
                sentid = checkforMultipleSent(sent)
                if pmid == "1460037" and sentid == 9:   # intermediate ( K d  -lo-’ M)
                    continue
                if pmid == "8916957" and sentid == 183:  # s~ggested,’’-’~
                    continue
                extracted_sent = getSentence(sentence_data, sentid)
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
            #exit()
    json_fn = open(args.outputJSON, 'w')
    json.dump(mydict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()
    problem_pmids = list(set(problem_pmids))
    print "Problem PMIDS: ", problem_pmids
    save_problem_pmids()


def getProteinDict(extracted_sent, stanford_tokens):
    mydict = dict()
    for prot_dict in extracted_sent["ProteinDetails"]:  # There may be multiple proteins detected in a single sentence
        if prot_dict["protName"] in mydict:
            print "For sentence id: ", extracted_sent["index"]
            print "The same protein already exists in the dictionary: ", prot_dict["protName"]
            exit()
        mydict[prot_dict["protName"]] = getProtTokenids(prot_dict, copy.deepcopy(extracted_sent), copy.deepcopy(stanford_tokens))
    return mydict


def getProtTokenids(prot_dict, sent_dict, tokens):
    mylist = list()
    mysent = sent_dict["sentenceText"]
    if check_for_unicode(mysent):
        mysent = mysent.decode("utf-8")
    tokens, mysent = handle_quotes(mysent, tokens)
    tokens_spans = get_span_tokens(mysent, tokens, sent_dict["startCharPos"])
    # The same protein may be detected multiple times in the same sentence, hence the for loop below
    for offsets in prot_dict["protcharOffset"]: # offsets is a list with 2 elements: [startpos, endpos]
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
            print "Sentence id: ", sent_dict["index"]
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
            print "For sentence id: ", extracted_sent["index"], " and label: ", trig_dict["label"]
            print "The same trigger word already exists in the dictionary: ", trig_dict["triggerWord"]
            exit()
        mydict[trig_dict["label"]][trig_dict["triggerWord"]] = getTriggerTokenids(trig_dict, copy.deepcopy(extracted_sent), copy.deepcopy(stanford_tokens))
    return mydict


def getTriggerTokenids(trig_dict, sent_dict, tokens):
    mylist = list()
    mysent = sent_dict["sentenceText"]
    if check_for_unicode(mysent):
        mysent = mysent.decode("utf-8")
    tokens, mysent = handle_quotes(mysent, tokens)
    tokens_spans = get_span_tokens(mysent, tokens, sent_dict["startCharPos"])
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
            print "Sentence id: ", sent_dict["index"]
            print "Sentence: ", sent_dict["sentenceText"]
            exit()
        mylist.append(tokenids)
    return mylist


def get_span_tokens(sent, tokens, startChar):
    span_tokens = []; offset = 0
    for token in tokens:
        offset = sent.find(token, offset)
        span_tokens.append((token, (startChar + offset, startChar + offset + len(token))))
        try:
            assert token == sent[offset:offset + len(token)], "token span does not point to the same string: %r \nSentence: %r" % (token, sent)
        except AssertionError:
            print "Offending token: ", token
            print "Sentence Id: ", sentid, "\n", sent
            print sent.encode('raw_unicode_escape')
            print tokens
            exit()
        offset += len(token)
    return span_tokens


def handle_quotes(mystr, tokens):
    global problem_pmids
    mystr = mystr.replace("\n", " ")
    if pmid == "18675783" and u"\xad" in mystr:
        mystr = mystr.replace(u"\xad", "")
        problem_pmids.append(pmid+":"+u"\xad")

    #print tokens
    for i, tok in enumerate(tokens):
        if tok == "``":
            if u"\u201c" in mystr:
                tokens[i] = u"\u201c"
            elif u"\u00AB" in mystr:
                tokens[i] = u"\u00AB"
                problem_pmids.append(pmid + ":" + u"\u00AB")
            elif mystr.count(u"\u2018") >= 2:
                tokens[i] = u"\u2018\u2018"
            else:
                tokens[i] = "\""
        if tok == "''":
            if u"\u201d" in mystr and "\"" not in mystr:
                tokens[i] = u"\u201d"
            elif u"\u201d" in mystr and "\"" in mystr:
                tokens[i] = u"\""
                mystr = mystr.replace(u"\u201d", "\"")
            elif mystr.count(u"\u2019") >= 2 and "\"" not in mystr:
                tokens[i] = u"\u2019\u2019"
            elif mystr.count(u"\u2019") >= 2 and "\"" in mystr:
                tokens[i] = "\""
                mystr = mystr.replace(u"\u2019", "\"")
            else:
                tokens[i] = "\""
        if tok == "\'s":
            if u"\u2019s" in mystr:
                tokens[i] = u"\u2019s"
            elif "\"s" in mystr:
                tokens[i] = "\"s"
        if tok == "--" and (u"\u2014" in mystr or u"\u2013" in mystr):
            tokens[i] = u"\u2013"
            mystr = mystr.replace(u"\u2014", u"\u2013")
        if tok == "`":
            if u"\u2018" in mystr:
                tokens[i] = u"\u2018"
            elif u"\u2039" in mystr:
                tokens[i] = u"\u2039"
                problem_pmids.append(pmid + ":" + u"\u2039")
            elif "`" not in mystr:
                tokens[i] = "'"
        if tok == "'" or ("'" in tok and len(tok) > 1 and tok != "''" and tok != "``'" and tok != "'``" and tok != "'''"):
            if u"\u2019" in mystr and u"\u2019\u2019" not in mystr:
                tokens[i] = tok.replace("'", u"\u2019")
            elif u"\u203a" in mystr:
                tokens[i] = u"\u203a"
        if u'\xa0' in tok and u'\xa0' not in mystr:
            tokens[i] = tok.replace(u'\xa0', u' ')
        if tok == "cents" and u"\u00A2" in mystr:
            tokens[i] = u"\u00A2"
            problem_pmids.append(pmid+":"+u"\u00A2")
        if tok == "#" and u"\u00A3" in mystr:
            tokens[i] = u"\u00A3"
            problem_pmids.append(pmid+":"+u"\u00A3")
        if tok == "$":
            if u"\x80" in mystr:
                tokens[i] = u"\x80"
                problem_pmids.append(pmid+":"+u"\x80")
            elif u"\u20ac" in mystr:
                tokens[i] = u"\u20ac"
                problem_pmids.append(pmid+":"+u"\u20ac")
            elif u"\u00a4" in mystr:
                tokens[i] = u"\u00a4"
                problem_pmids.append(pmid+":"+u"\u00a4")
        if tok == "-" and u"\xad" in mystr:
            tokens[i] = u"\xad"
            problem_pmids.append(pmid+":"+u"\xad")
        if tok == "..." and ". . ." in mystr:
            tokens[i] = ". . ."
        if tok == "```" and u"\u2018\u201c" in mystr:
            tokens[i] = u"\u2018\u201c"
        if tok == "``'" and u"\u201c\u2019" in mystr:
            tokens[i] = u"\u201c\u2019"
        if tok == "'``" and u"\u2019\u201c" in mystr:
            tokens[i] = u"\u2019\u201c"
        if tok == "`''" and u"\u2018\u201d" in mystr:
            tokens[i] = u"\u2018\u201d"
        if tok == "''`" and u"\u201d\u2018" in mystr:
            tokens[i] = u"\u201d\u2018"
        if tok == "'''" and u"\u2019\u201d" in mystr:
            tokens[i] = u"\u2019\u201d"
        if tok == "''''" and u"\u201d\u201d" in mystr:
            tokens[i] = u"\u201d\u201d"
        if tok == "''``" and u"\u201d\u201c" in mystr:
            tokens[i] = u"\u201d\u201c"
    return tokens, mystr


def check_for_unicode(mystring):
    try:
        mystring.decode('ascii')
    except UnicodeDecodeError:
        #print mystring
        return True
    else:
        return False


def save_problem_pmids():
    mydict = {}
    for elem in problem_pmids:
        artid = elem.split(":")[0].strip()
        badchar = elem.split(":")[1].strip()
        if artid in mydict:
            if badchar not in mydict[artid]:
                mydict[artid].append(badchar)
        else:
            mydict[artid] = [badchar]
    pickle.dump(mydict, open("NonOA-BadChars.p", "wb"))


def getSentencesFile(Extr_Sent_Dir, pmid):
    sent_file = get_file(Extr_Sent_Dir, pmid)
    file_handle = open(sent_file, "rb")
    sentfile_data = json.load(file_handle)  # sentfile_data is a dictionary
    file_handle.close()
    if "sentenceList" in sentfile_data and "@items" in sentfile_data["sentenceList"] and len(sentfile_data["sentenceList"]["@items"]) > 0:
        return sentfile_data    # There are extracted sentences for this article!
    else:
        return {}


def getSentence(sentence_data, sentid):
    mysent = sentence_data["sentenceList"]["@items"][sentid]
    if mysent["index"] != sentid:
        print "We got the wrong sentence: \n", mysent["cleansedText"]
        print "Sentence Id obtained from the xml file: ", sentid + 1
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

    #abbrev = ["AL.", "Inc.", "RT."]
    for i, tok in enumerate(mylist):
        if tok == "-LRB-":
            mylist[i] = "("
        elif tok == "-RRB-":
            mylist[i] = ")"
        elif tok == "-LSB-":
            mylist[i] = "["
        elif tok == "-RSB-":
            mylist[i] = "]"
        elif tok == "-LCB-":
            mylist[i] = "{"
        elif tok == "-RCB-":
            mylist[i] = "}"
        elif tok.endswith(".") and len(tok) > 1 and i < len(mylist) - 1 and mylist[i + 1] == "." and tok != "...":
            mylist[i] = tok[:-1]
        elif "-LRB-" in tok and "-RRB-" in tok:
            mylist[i] = tok.replace("-LRB-", "(")
            mylist[i] = mylist[i].replace("-RRB-", ")")
        elif "-RRB-" in tok:
            mylist[i] = tok.replace("-RRB-", ")")
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


def ignoreSent(pmid, sentid):
    if pmid == "10022866" and sentid in [139, 162]:    # ['1786', 'MAROUN', 'ET', 'AL.', '.', 'MOL', '.'] Extra period token
        return True
    return False


if __name__ == "__main__":
    main_body()
