# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os
import re
import string
import operator;
from pprint import pprint;
import cPickle as pickle;
from nltk import ngrams
from nltk.tokenize import word_tokenize
import fst;
from Utilities import read_config_file, check_for_unicode, tokenizeOnHyphen

# PubMedIDS_with_Extras_OpenAccess.txt contains the 61 OpenAccess articles. The extracted sentences files for these articles have a different format.
# python New_Protein_Detector_OpenAccess.py ./Minimized_Protein_Names_FST.fst ./SymbolTable_Protein_Names_FST.sym 19 ./PubMedIDS_with_Extras_OpenAccess.txt

def main_body():
    parser = argparse.ArgumentParser(prog='New_Protein_Detector_OpenAccess', usage='New_Protein_Detector_OpenAccess.py <FST file> <SymbolTable file> <Maximum Length of a protein name> <PubMedfilelist> [--file <Optional Filename to be parsed>]', description='Script to read the segmented sentences and recognize protein names')
    parser.add_argument('fst_file', help='Name of the binary FST file')
    parser.add_argument('symboltable_file', help='Name of the binary symbol table file')
    parser.add_argument('max_prot_len', help='Maximum Length of a protein name')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('--file', help='Name of file to be parsed')

    args = parser.parse_args()
    #print args.fst_file, args.symboltable_file, args.max_prot_len, args.file
    if args.file==None:
        fst_init(args.fst_file, args.symboltable_file)
        read_dirs(int(args.max_prot_len), os.path.join(txt_location, args.PubMedfilelist))
    else:
        pass  #read_file(args.fst_file, args.symboltable_file, int(args.max_prot_len), args.file)


def fst_init(fst_file, symboltable_file):
    global my_fst; global syms; global stopwords; global common_english_words
    global parent_location; global txt_location; global json_location; global pickle_location
    parent_location, txt_location, json_location, pickle_location = read_config_file("./config.cfg")

    my_fst = fst.read_std(os.path.join(parent_location, "FST_files", fst_file))
    syms = fst.read_symbols(os.path.join(parent_location, "FST_files", symboltable_file))
    stopwords = pickle.load(open(os.path.join(pickle_location, 'stop_words_dict.p'), 'rb'))
    common_english_words = pickle.load(open('common_english_words_proteins_dict.p', 'rb'))


def read_dirs(max_prot_len, PubMedfilelist):
    #fst_init(fst_file, symboltable_file)
    dict_nameToid = pickle.load(open('dictionary_nameToUniProtId.p', 'rb')); max_ngram_len = 0
    rootDir = '/ua/ml-group/big-mechanism-project/PLEIADES/workDir/JSON_SENTENCE_ASSOCS_OpenAccess'
    writeDir = '/ua/ml-group/big-mechanism-project/PLEIADES/workDir/protein_detect_output_OpenAccess'
    problem_proteins = ["plcε", "ckiε", "ck1ε", "fcεri", "14-3-3ε"]

    with open(PubMedfilelist, 'rt') as f1:
        for pmid in f1:
            if not pmid.strip(): continue                
            fname_string1 = '%s*' % (os.path.join(rootDir, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                continue                 
            myfile = fname[0].strip()

            write_fname = pmid.strip() + '_ProteinMatches.json'            
            fullname_write_fname = os.path.join(writeDir, write_fname) 
            if os.path.exists(fullname_write_fname):
                print 'Skipping file: ', myfile, ' as it has already been processed!'
                continue
            print 'Processing File: ', pmid.strip()
            sys.stdout.flush()
            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
            file_handle.close()
            
            myjson=[]
            if len(myfile_data["paragraphList"]) == 0: continue;
            for para in myfile_data["paragraphList"]:   # para is a dictionary
                if "sentenceList" not in para:  continue;
                if "@items" not in para["sentenceList"]: continue;
                if len(para["sentenceList"]["@items"]) == 0: continue

                for itm in para["sentenceList"]["@items"]:  # itm is a dictionary. It represent one sentence in the current paragraph

                    if "ProteinDetails" in itm:
                        itm.pop("ProteinDetails", None)
                    if "TriggerDetails" in itm:
                        itm.pop("TriggerDetails", None)
                    if "HasProtein" in itm:
                        itm.pop("HasProtein", None)
                    mysent = itm["cleansedText"]
                    if check_for_unicode(mysent):
                        mysent = mysent.decode("utf-8")
                    matches = find_proteins(mysent.lower(), max_prot_len)

                    temp = {}
                    if len(matches) > 0:
                        temp["id"] = pmid.strip() + "_" + str(itm["absoluteId"])
                        temp["text"] = mysent
                        temp["matches"] = {}
                        for m in matches:   # matches is a list of tuples. So m looks like (match_strg, absindex)
                            match_name = m[0].strip()
                            exact_mention = temp["text"][m[1][0]:m[1][1]]
                            assert exact_mention.lower() == match_name
                            if isstopword(match_name, exact_mention):  continue;
                            if match_name in problem_proteins:
                                matched_ids = getid_probprot(match_name)
                            else:
                                matched_ids = dict_nameToid[match_name]
                            if match_name == "cab-":  # Temporary Hack
                                if mysent.lower().find("cab-1") > -1:
                                    match_name = "cab-1"
                                else:
                                    continue
                            matched_ids_string = ', '.join(matched_ids)
                            if match_name in temp["matches"]:
                                temp["matches"][match_name]["Index"].append(m[1])
                                if exact_mention not in temp["matches"][match_name]["Mentions"]:
                                    temp["matches"][match_name]["Mentions"].append(exact_mention)
                            else:
                                temp["matches"][match_name] = {}
                                temp["matches"][match_name]["UniProtId"] = matched_ids_string
                                temp["matches"][match_name]["Index"] = [m[1]]
                                temp["matches"][match_name]["Mentions"] = [exact_mention]
                        if temp["matches"]: 
                            myjson.append(temp)
                            itm["HasProtein"] = "Yes"
                            populateproteinkey(itm, temp)

            if len(myjson) > 0:
                json_fn = open(fullname_write_fname, 'w')
                json.dump(myjson, json_fn, indent=4, ensure_ascii=False)
                json_fn.close()

            myfile_write_fname = os.path.join(rootDir, pmid.strip()+"_OA")
            os.remove(myfile_write_fname)
            json_fn = open(myfile_write_fname, 'w')
            json.dump(myfile_data, json_fn, indent=4, ensure_ascii=False)
            json_fn.close()
            
    
def find_proteins(sentence, max_prot_len):
    numTokens = getNumberOfTokens(sentence)
    if numTokens >= max_prot_len:
        max_ngram_len = max_prot_len
    else:
        max_ngram_len = numTokens
    sent_queue = [sentence+"|_|0|_|"+str(max_ngram_len)]; matches = []
    while len(sent_queue) > 0:
        front_elem = sent_queue.pop()
        sent, abs_start_idx, curr_ngram_len = front_elem.split("|_|")
        for size in range(int(curr_ngram_len), 0, -1):
            ngram_list = getNgrams(sent, int(abs_start_idx), size)
            match = False; match_strg = ""
            for strg, absindex in ngram_list:
                str_fst = fst.linear_chain(strg,syms=syms)
                str_fst.arc_sort_input()
                intersect_fst = my_fst & str_fst
                try:
                    for path in intersect_fst.paths():  # if there is a match, intersect_fst should only have one path
                       match_strg = ''.join(syms.find(arc.ilabel) for arc in path)                    
                    if match_strg != strg:
                        print "Something is wrong! Matched FST String: " + match_strg + "\t Ngram from text: " +  strg
                        exit()
                    match = True; break
                except KeyError:
                    match = False
            if match:
                matches.append((match_strg, absindex))
                relindex = getRelativeIndices(int(abs_start_idx), absindex)
                try:
                    assert match_strg == sent[relindex[0]:relindex[1]]
                except AssertionError:
                    print "Matched String: ", match_strg
                    print "Sentence: ", sent
                    print "Protein found in above sentence: ", sent[relindex[0]:relindex[1]]
                    exit()
                left_partition, lnumspaces = getStrippedSent(sent[:relindex[0]], "L")
                right_partition, rnumspaces = getStrippedSent(sent[relindex[1]:], "R")
                if len(left_partition) > 0:
                    sent_queue.append(left_partition+"|_|"+abs_start_idx+"|_|"+str(size-1))
                if len(right_partition) > 0:
                    sent_queue.append(right_partition + "|_|" + str(absindex[1]+rnumspaces) + "|_|" + str(size))
                break  

    return matches  


def isstopword(myword, exactmention):
    if myword == "fig" and exactmention == "Fig":   # Assumption: There is no positive protein in DatumKB "Fig"
        return True
    start_char = myword[0]
    try:
        if myword in stopwords[start_char]:
            return True
        else:
            if start_char in common_english_words and myword in common_english_words[start_char] and exactmention.islower():
                return True
            return False
    except KeyError:
        if start_char in common_english_words and myword in common_english_words[start_char] and exactmention.islower():
            return True
        return False


def myngrams(inpt, n):
  inpt = inpt.split(' ')
  output = []
  for i in range(len(inpt)-n+1):
    output.append(inpt[i:i+n])
  return output


def getid_probprot(match_name):
    if match_name == "plcε":
        return ["Q9P212"]
    elif match_name == "ckiε":
        return ["P49674"]
    elif match_name == "ck1ε" or match_name == "ck 1ε":
        return ["P48729"]
    elif match_name == "fcεri":
        return ["P12319"]
    elif match_name == "14-3-3ε":
        return ["P62258"]


def getNumberOfTokens(sentence):
    tokens = word_tokenize(sentence)
    tokens = tokenizeOnHyphen(tokens)
    return len(tokens)


def getNgrams(sent, abs_start_idx, size):
    # sent may actually be a part of the original sentence
    tokens = word_tokenize(sent)
    tokens = tokenizeOnHyphen(tokens)
    tokens_spans = get_span_tokens(sent, tokens)
    tokens = [t[0] for t in tokens_spans]
    ngram_spans = []
    for i, each_ngram in enumerate(list(ngrams(tokens, size))):
        start_idx = tokens_spans[i][1][0]              # An element of tokens_spans looks like: (u'guanine', (112, 119))
        end_idx = tokens_spans[i+(size-1)][1][1]
        ngram_string = sent[start_idx:end_idx]
        ngram_spans.append((ngram_string,(abs_start_idx+start_idx, abs_start_idx+end_idx)))
    return ngram_spans


def get_span_tokens(sent, tokens):
    span_tokens = []; offset = 0
    for token in tokens:
        if token == "``" or token == "''":
            token = "\""
        offset = sent.find(token, offset)
        span_tokens.append((token, (offset, offset + len(token))))
        try:
            assert token == sent[offset:offset + len(token)]
        except AssertionError:
            print "Offending token: ", token
            print "Token found in sentence: ", sent[offset:offset + len(token)]
            print "Sentence: ", sent
            print "All Tokens: ", tokens
            exit()
        offset += len(token)
    return span_tokens


def getRelativeIndices(abs_start_idx, absindex):
    start_idx = absindex[0] - abs_start_idx
    end_idx = absindex[1] - abs_start_idx
    return (start_idx, end_idx)


def getStrippedSent(mysentpart, LorR):
    if len(mysentpart.strip()) == 0:
        return "", 0
    orig_len = len(mysentpart)
    if LorR == "L":
        sentpart = mysentpart.rstrip()
        new_len = len(sentpart)
    else:
        sentpart = mysentpart.lstrip()
        new_len = len(sentpart)
    return sentpart, orig_len - new_len


def populateproteinkey(sentDetails, protMatch):
    sent_lower = sentDetails["sentenceText"].lower()
    prot_all_indices = []
    sentDetails["ProteinDetails"] = []
    sorted_matches = sorted(protMatch["matches"], key=len, reverse=True)
    for prot in sorted_matches:
        tokens = prot.lower().split()
        pattern = re.compile("\s+".join(tokens))
        start_idx = 0; temp = {}
        temp["protName"] = prot
        temp["protMentions"] = protMatch["matches"][prot]["Mentions"]
        for i in range(len(protMatch["matches"][prot]["Index"])):
            protfound = False
            while not protfound and start_idx < len(sent_lower):
                m = pattern.search(sent_lower, start_idx)
                if m:
                    start_idx = m.end()
                    assert prot == " ".join(sent_lower[m.start():m.end()].strip(' \t\n\r').split())
                    if check_protein_wordbounds(m.start(), m.end(), sent_lower):
                        if check_protein_overlap(sentDetails["charOffset"] + m.start(), sentDetails["charOffset"] + m.end(), prot_all_indices):
                            protfound = False   # There should not be any overlap between previously detected proteins in the same sentence
                            continue
                        else:
                            protfound = True
                    else:
                        protfound = False
                        continue
                    prot_all_indices.append((sentDetails["charOffset"] + m.start(), sentDetails["charOffset"] + m.end()))
                    if i == 0:
                        temp["protcharOffset"] = [(sentDetails["charOffset"] + m.start(), sentDetails["charOffset"] + m.end(), "rgb(33, 255, 17);")]
                    else:
                        temp["protcharOffset"].append((sentDetails["charOffset"] + m.start(), sentDetails["charOffset"] + m.end(), "rgb(33, 255, 17);"))
                else:
                    start_idx = len(sent_lower)
            if not protfound:
                print "Something went wrong! The protein: ", prot, " was not found in the sentence: \n", sentDetails["sentenceText"]
                print sentDetails["sentenceText"].encode('raw_unicode_escape')
                print "The matched protein is part of a bigger word!"
                print "UniProtID: ", protMatch["matches"][prot]["UniProtId"]
                exit()

        if len(protMatch["matches"][prot]["Index"]) != len(temp["protcharOffset"]):
            print "The no. of proteins found in sentence do not match the no. of proteins detected by FSA"
            print sentDetails["charOffset"]; print protMatch["matches"][prot]["Index"]; print temp["protcharOffset"]
            print prot; print protMatch["matches"][prot]["UniProtId"]; print sentDetails["sentenceText"]
            exit()

        sentDetails["ProteinDetails"].append(temp)


def check_protein_wordbounds(start, end, sentence):
    if start > 0 and end < len(sentence):
        if (not sentence[start - 1].isspace() and not isPunct(sentence[start - 1])) or \
                (not sentence[end].isspace() and not isPunct(sentence[end])):  # match at the middle of the sentence
            return False
    if start == 0 and end < len(sentence):  # match at the beginning of the sentence
        if not sentence[end].isspace() and not isPunct(sentence[end]):
            return False
    if start > 0 and end == len(sentence):  # match at the end of the sentence
        if not sentence[start - 1].isspace() and not isPunct(sentence[start - 1]):
            print "Protein found is part of another bigger word!", sentence[start:end], "\n", sentence
            exit()
    return True


def isPunct(mychar):
    for c in set(string.punctuation):           # prot is "IL" and sentence contains words like "IL-4" and "IL-5"
        #if c != "-" and c == mychar: return 1   # Sentence does contain the word "IL" standalone, ofcourse!
        if c == mychar:
            return 1                # [UPDATE] 6/23/2017: Lets consider "-" as a punctuation Eg: DDX41-shrna
    if mychar == u"\u2013" or mychar == u"\u2014" or mychar == u"\u2212":
        return 1
    return 0                                    # Do not consider "-" as a punctuation because word_tokenize doesn't!


def check_protein_overlap(start_idx, end_idx, prot_all_indices):
    #print start_idx, end_idx
    #print prot_all_indices
    for indx in prot_all_indices:
        if start_idx >= indx[0] and end_idx <= indx[1]:
            return True
    return False


if __name__ == "__main__":
    main_body()
