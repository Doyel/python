# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import re;
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import operator;
from pprint import pprint;
import cPickle as pickle;
import fst;

# PubMedIDS_with_Extras_OpenAccess.txt contains the 61 OpenAccess articles. The extracted sentences files for these articles have a different format.
# python Protein_Detector_OpenAccess.py ./Minimized_Protein_Names_FST.fst ./SymbolTable_Protein_Names_FST.sym 19 ./PubMedIDS_with_Extras_OpenAccess.txt

def main_body():
    parser = argparse.ArgumentParser(prog='Protein_Detector', usage='Protein_Detector.py <FST file> <SymbolTable file> <Maximum Length of a protein name> <PubMedfilelist> [--file <Optional Filename to be parsed>]', description='Script to read the segmented sentences and recognize protein names')
    parser.add_argument('fst_file', help='Name of the binary FST file')
    parser.add_argument('symboltable_file', help='Name of the binary symbol table file')
    parser.add_argument('max_prot_len', help='Maximum Length of a protein name')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('--file', help='Name of file to be parsed')

    args = parser.parse_args()
    #print args.fst_file, args.symboltable_file, args.max_prot_len, args.file
    if args.file==None:
        read_dirs(args.fst_file, args.symboltable_file, int(args.max_prot_len), args.PubMedfilelist)
    else:
        pass  #read_file(args.fst_file, args.symboltable_file, int(args.max_prot_len), args.file)


def fst_init(fst_file, symboltable_file):
    global my_fst
    global syms
    global badfile_bit
    global stopwords
    my_fst = fst.read_std(fst_file)
    syms = fst.read_symbols(symboltable_file)
    stopwords = pickle.load(open('stop_words_dict.p', 'rb'))


def read_dirs(fst_file, symboltable_file, max_prot_len, PubMedfilelist):
    fst_init(fst_file, symboltable_file)
    dict_nameToid = pickle.load(open('dictionary_nameToUniProtId.p', 'rb')); max_ngram_len = 0
    rootDir = '/ua/ml-group/big-mechanism-project/PLEIADES/Sep2016_SMA/JSON_SENTENCE_ASSOCS_OpenAccess'
    writeDir = '/ua/ml-group/big-mechanism-project/PLEIADES/Sep2016_SMA/protein_detect_output_OpenAccess'
    problem_proteins = ["plcε", "ckiε", "ck1ε", "fcεri", "14-3-3ε"]
    badfilename = writeDir + '/BadFiles.txt'
    bf = open(badfilename, 'a'); 

    with open(PubMedfilelist, 'rt') as f1:
        for pmid in f1:
            if not pmid.strip(): continue                
            fname_string1 = '%s*' % (os.path.join(rootDir, pmid.strip()))
            fname = glob.glob(fname_string1)
            if len(fname) == 0:
                print "Could not locate extracted sentences file for PubMedID: ", pmid.strip()
                continue                 
            myfile = fname[0].strip()
            badfile_bit = False
            
            filename = os.path.basename(myfile)
            write_fname = pmid.strip() + '_ProteinMatches.json'            
            fullname_write_fname = os.path.join(writeDir, write_fname) 
            if os.path.exists(fullname_write_fname):
                print 'Skipping file: ', myfile, ' as it has already been processed!'
                continue
            print 'Processing File: ', pmid.strip()
            
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

                    matches = []; ctr = 0; old_text = ""; new_text = itm["cleansedText"].replace("\n", " ").replace("\r", "")
                    while len(matches) == 0 and ctr <= 3:
                        if new_text != old_text:
                            matches = find_proteins(new_text.lower(), max_prot_len)
                        if len(matches) > 0:    break
                        ctr += 1; old_text = new_text
                        if ctr == 1:
                            if "-" in itm["cleansedText"]:
                                new_text = itm["cleansedText"].replace("-", " ")
                        if ctr == 2:
                            if "-" in itm["cleansedText"]:
                                new_text = itm["cleansedText"].replace("-", "")
                        if ctr == 3:
                            if "/" in itm["cleansedText"]:
                                new_text = itm["cleansedText"].replace("/", " ")

                    temp = {}
                    if len(matches) > 0:
                        temp["id"] = pmid.strip() + "_" + str(itm["absoluteId"])
                        temp["text"] = itm["cleansedText"]
                        temp["matches"] = {}
                        for m in matches:
                            match_name = m.split("|_|")[0].strip()
                            if isstopword(match_name):  continue;
                            if match_name in problem_proteins:
                                matched_ids = getid_probprot(match_name)
                            else:
                                matched_ids = dict_nameToid[match_name]
                            matched_ids_string = ', '.join(matched_ids)
                            temp["matches"][match_name] = matched_ids_string
                        if temp["matches"]: 
                            myjson.append(temp)
                            itm["HasProtein"] = "Yes"

            if len(myjson) > 0:
                json_fn = open(fullname_write_fname, 'w')
                json.dump(myjson, json_fn, indent=4, ensure_ascii=False)
                json_fn.close()

            if badfile_bit:
                bf.write(myfile+'\n')

            myfile_write_fname = os.path.join(rootDir, pmid.strip()+"_OA")
            os.remove(myfile_write_fname)
            json_fn = open(myfile_write_fname, 'w')
            json.dump(myfile_data, json_fn, indent=4, ensure_ascii=False)
            json_fn.close()

    bf.close()
            
    
def find_proteins(sentence, max_prot_len):
    global badfile_bit;
    if len(sentence.split(' ')) >= max_prot_len:
        max_ngram_len = max_prot_len
    else:
        max_ngram_len = len(sentence.split(' '))
    sent_queue = [sentence+"|_|"+str(max_ngram_len)]; matches = []
    while len(sent_queue) > 0:
        front_elem = sent_queue.pop(); 
        sent, max_ngram_len  = front_elem.split("|_|");
        max_ngram_len = int(max_ngram_len); #print sent
        for size in range(max_ngram_len, 0, -1):                            
            ngram_list = [' '.join(x) for x in ngrams(sent, size)]
            match = False; match_strg = ""
            for index, strg in enumerate(ngram_list):
                if len(strg.strip())==0: continue
                str_fst = fst.linear_chain(strg,syms=syms)
                str_fst.arc_sort_input()
                intersect_fst = my_fst & str_fst
                try:
                    for path in intersect_fst.paths():  # if there is a match, intersect_fst should only have one path
                       match_strg = ''.join(syms.find(arc.ilabel) for arc in path)                    
                    if match_strg != strg:
                        print "Something is wrong! Matched FST String: " + match_strg + "\t Ngram from text: " +  strg                        
                    match = True; break
                except KeyError:
                    match = False;
            if match:
                try:
                    partitioned_strg = splice_sentence(match_strg, sent, size, index, len(ngram_list)-1)
                    sent_queue.extend(partitioned_strg)
                    # The below 4 lines might be buggy. Praying I dont have to come back to this shit again
                    full_sent = sentence.decode("utf-8"); sent = sent.decode("utf-8")
                    sub_sent_start_idx = full_sent.find(sent)   # match_strg is a substring of sent
                    if index == 0:    # Match occurred at the beginning of sent
                        start_idx = sub_sent_start_idx + sent.find(match_strg+" ")
                        end_idx = start_idx + len(match_strg)
                    elif index == len(ngram_list)-1:   # Match occurred at the end of sent
                        start_idx = sub_sent_start_idx + sent.find(" "+match_strg)
                        end_idx = (start_idx + 1) + len(match_strg); start_idx += 1
                    else:   # Match occurred at the middle of sent
                        start_idx = sub_sent_start_idx + sent.find(" "+match_strg+" ")
                        end_idx = (start_idx + 1) + len(match_strg); start_idx += 1
                    matches.append(match_strg+"|_|"+str(start_idx)+","+str(end_idx))
                except UnicodeDecodeError:
                    badfile_bit = True;
                    pass
                break  

    return matches  


def isstopword(myword):
    lower_var = myword.strip().lower()
    start_char = lower_var[0]
    try:
        if myword.lower() in stopwords[start_char]:
            return True
        else:
            return False
    except KeyError:
        return False


def ngrams(inpt, n):
  inpt = inpt.split(' ')
  output = []
  for i in range(len(inpt)-n+1):
    output.append(inpt[i:i+n])
  return output


def splice_sentence(matched_string, full_sentence, ngram_size, index, ngram_len):
    #print full_sentence; print matched_string
    #print "Splice_sentence: ",full_sentence 
    full_sentence = full_sentence.decode("utf-8")
    if matched_string == full_sentence:
        return []
    if index == 0:    # Match occurred at the beginning of the sentence
        start_index = full_sentence.find(matched_string+" ")
        end_index = start_index + len(matched_string)
        #print start_index, end_index
    elif index == ngram_len:   # Match occurred at the end of the sentence
        start_index = full_sentence.find(" "+matched_string)
        end_index = (start_index + 1) + len(matched_string)
    else:   # Match occurred at the middle of the sentence
        start_index = full_sentence.find(" "+matched_string+" ")
        end_index = (start_index + 1) + len(matched_string)
    partition1 = full_sentence[:start_index].strip().encode("utf-8"); #print "Partition 1: ", partition1
    partition2 = full_sentence[end_index:].strip().encode("utf-8"); #print "Partition 2: ", partition2
    if len(partition1) == 0 and len(partition2) > 0:
        return [partition2+"|_|"+str(ngram_size)]
    elif len(partition1) > 0 and len(partition2) == 0:
        if ngram_size == 1:
            return []
        else:
            return [partition1+"|_|"+str(ngram_size-1)]
    elif len(partition1) == 0 and len(partition2) == 0:
        return []
    else:
        if ngram_size == 1:
            #print "Partition 1: ", partition1
            #print "Partition 2: ", partition2
            return [partition2+"|_|"+str(ngram_size)]
        else:
            return [partition1+"|_|"+str(ngram_size-1), partition2+"|_|"+str(ngram_size)]


def check_for_unicode(mystring):
    try:
        mystring.decode('ascii')
    except UnicodeDecodeError:
        print mystring
        return True
    else:
        return False


def getid_probprot(match_name):
    if match_name == "plcε":
        return list("Q9P212")
    elif match_name == "ckiε":
        return list("P49674")
    elif match_name == "ck1ε":
        return list("P48729")
    elif match_name == "fcεri":
        return list("P12319")
    elif match_name == "14-3-3ε":
        return list("P62258")


if __name__ == "__main__":
    main_body()
