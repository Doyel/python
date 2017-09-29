# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import fst;
import json;
import csv;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
import re, os
import cPickle as pickle;
from pprint import pprint;
from nltk.tokenize import word_tokenize
from Utilities import read_config_file, check_for_unicode, tokenizeOnHyphen

# python Protein_Extractor_UniProtID.py pro.obo.txt

def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Protein_Extractor_UniProtID', usage='Protein_Extractor_UniProtID.py <ontology_file>', description='Script to extract protein names from an ontology')
    parser.add_argument('ProtOntofile', help='Ontology file to be parsed')

    args = parser.parse_args()

    ontologyDict(args.ProtOntofile)


def ontologyDict(ontolFile):
    l=[]
    #print(ontolFile)
    with open(os.path.join(txt_location, ontolFile), 'rt') as f:
        for line in f:
            if line.strip() and line.strip().split()[0] == "[Term]":
                d = {}; listSynonym = []; hasSynm = 0;
                for item in f:
                    if not item.strip():
                        break
                    (key, value) = item.split(':',1)
                    if key.strip() == "synonym":
                        hasSynm = 1
                        listSynonym.append(value.strip().lower())
                        continue
                    d[key.strip()]=value.strip()
                if hasSynm:
                    d['synonym'] = listSynonym
                l.append(d)
               
    print "Total no of terms in ontology: ", len(l)  # l is a list of dictionaries
    extractProNames(l)


def extractProNames(l):
    HGNC_dict_BasicOps_lower = {}; HGNC_hugosym_spno_lower = {}
    # HGNC_dict_BasicOps will have all hugosyms as keys even if a hugosym does not have a synonym. In that case, the value of that dictionary entry will be an empty list.
    parse_BasicOps(os.path.join(json_location, 'basicops.json'), HGNC_dict_BasicOps_lower, HGNC_hugosym_spno_lower)
    #pprint(HGNC_dict_BasicOps_lower); exit()

    HGNC_gene_proteins_lower = {}; HGNC_symbol_id_lower = {}
    parse_HGNC_gene_proteins(HGNC_gene_proteins_lower, HGNC_symbol_id_lower)
    #pprint(HGNC_gene_proteins_lower); exit()

    HGNC_hgncid_uniprotid = {}; HGNC_symbol_uniprot_lower = {}
    parse_HGNC_gene_uniprot(os.path.join(txt_location, 'HGNC_gene_with_UNIPROT_IDs.txt'), HGNC_hgncid_uniprotid, HGNC_symbol_uniprot_lower)

    BasicOps_proteins = HGNC_dict_BasicOps_lower.keys()
    HGNC_proteins = HGNC_gene_proteins_lower.keys()
    Only_BasicOps = list(set(BasicOps_proteins).difference(HGNC_proteins)) # Proteins that are only present in BasicOps.json but not in HGNC_gene_with_protein_product.txt
    #print Only_BasicOps; exit()
  
    protein_names=[]; max_prot_len = 0; dict_name_id = {}; tot_items = 0; 
    for li in l: # li is a dictionary
        d = {};

        if li['id'].split(':')[0].strip() == "HGNC":
            prot_name = chkParens(li['name'])

        # Get all the proteins that are not present in HGNC_gene_proteins.txt but present in pro.obo.txt
        if (li['id'].split(':')[0].strip() == "HGNC" and prot_name.lower() not in HGNC_gene_proteins_lower) or (li.has_key('is_a') and (li['is_a'].strip() == "PR:000003292 ! NF-kappaB inhibitor" or li['is_a'].strip() == "PR:000001774 ! inhibitor of nuclear factor kappa-B kinase complex alpha/beta subunits")):
            if li['id'].split(':')[0].strip() == "HGNC":
                d['id'] = HGNC_hgncid_uniprotid[li['id'].strip()]
            else:
                d['id'] = li['id'].strip(); #print d['id']
            #match = re.search(r'.*, .*', li['name'])
            #d['name'] may contain names of the form - abcdef, mitochondria
            prot_name = chkParens(li['name']); #print prot_name
            if prot_name[-1] == '.':
                print prot_name
            if len(prot_name.split()) > max_prot_len:
                max_prot_len = len(prot_name.split())
            d['name'] = prot_name.lower(); tot_items +=1

            if d['name'] in dict_name_id:            
                if d['id'] not in dict_name_id[d['name']]:
                    dict_name_id[d['name']].append(d['id'])
            else: dict_name_id[d['name']] = [d['id']]

            protname_var = generate_variants(d['name'])     # Get all variants of the protein name

            listSynonym = []; has_synonyms = False
            if li.has_key('synonym'):   # Both PR and HGNC terms may have synonyms
                has_synonyms = True         
                for syn in li['synonym']:
                    match = re.match(r'^\".*?\"',syn)
                    extr_syn = syn[match.start()+1:match.end()-1]
                    listSynonym.append(chkParens(extr_syn.lower())); 
                    if len(listSynonym[-1].split()) > max_prot_len:
                        max_prot_len = len(listSynonym[-1].split())
                    listSynonym.extend(generate_variants(listSynonym[-1]))      # Get all variants of syn, add to listSynonym

            if (li['id'].split(':')[0].strip() == "HGNC") and (d['name'] in HGNC_dict_BasicOps_lower):
                basicops_syns = HGNC_dict_BasicOps_lower[d['name']]     # Get synonyms from BasicOps for HGNC terms that are present in pro.obo.txt but not in HGNC_gene_with_protein_product.txt
                if len(basicops_syns) > 0:
                    has_synonyms = True; 
                    for s in basicops_syns:
                        listSynonym.append(s);
                        if len(s.split()) > max_prot_len:
                            max_prot_len = len(s.split())
                        listSynonym.extend(generate_variants(s))                    

            if protname_var:
                listSynonym.extend(protname_var)    # Append variants of d['name'] to listSynonym, set has_synonyms to True  
                has_synonyms = True                    
                                                           
            if has_synonyms:
                listSynonym_uniq = list(set(listSynonym)) 
                d['synonym'] = listSynonym_uniq
                tot_items += len(listSynonym_uniq)
                for sym in listSynonym_uniq:
                    if sym in dict_name_id:            
                        if d['id'] not in dict_name_id[sym]:
                            dict_name_id[sym].append(d['id'])
                    else: dict_name_id[sym] = [d['id']]

            protein_names.append(d)

    for symbol in HGNC_gene_proteins_lower.keys():
        d = {}
        d['id'] = HGNC_hgncid_uniprotid[HGNC_symbol_id_lower[symbol]]
        d['name'] = symbol; tot_items +=1
        if len(symbol.split()) > max_prot_len:
            max_prot_len = len(symbol.split())
        if d['name'] in dict_name_id:            
            if d['id'] not in dict_name_id[d['name']]:
                dict_name_id[d['name']].append(d['id'])
        else: dict_name_id[d['name']] = [d['id']]
        protname_var = generate_variants(d['name'])     # Get all variants of the protein name. protname_var may be empty
        
        # Hacks suggested by Mark
        if symbol == "hras":
            protname_var.extend(["h-ras", "h ras"])
        if symbol == "mras":
            protname_var.extend(["m-ras", "m ras"])
        if symbol == "nras":
            protname_var.extend(["n-ras", "n ras"])

        HGNC_syns = HGNC_gene_proteins_lower[symbol]
        if symbol in HGNC_dict_BasicOps_lower:
            BasicOps_syns = HGNC_dict_BasicOps_lower[symbol]
        else: BasicOps_syns = []
        synonyms_union = HGNC_syns + BasicOps_syns    # synonyms_union may be an [] list
        synonyms = list(set(synonyms_union))

        all_synonyms = []
        if len(synonyms) > 0:            
            for syn in synonyms:
                all_synonyms.append(syn);
                if len(syn.split()) > max_prot_len:
                    max_prot_len = len(syn.split())
                all_synonyms.extend(generate_variants(syn))
            
        if protname_var:
                all_synonyms.extend(protname_var)
        all_synonyms_uniq = list(set(all_synonyms))

        if len(all_synonyms_uniq) > 0:
            d['synonym'] = all_synonyms_uniq
            tot_items += len(all_synonyms_uniq)
            for sym in all_synonyms_uniq:
                if sym in dict_name_id:            
                    if d['id'] not in dict_name_id[sym]:
                        dict_name_id[sym].append(d['id'])
                else: dict_name_id[sym] = [d['id']]

        protein_names.append(d)

    # Include all the proteins that are only present in BasicOps.json but not in HGNC_gene_proteins.txt
    for hugosym in Only_BasicOps:
        d = {}
        if hugosym == "tcrb1":
            d['id'] = "P01850"
        elif hugosym == "trb@":
            d['id'] = "P04435"
        else:
            if hugosym in HGNC_symbol_uniprot_lower:
                d['id'] = HGNC_symbol_uniprot_lower[hugosym]
            else:
                d['id'] = HGNC_hugosym_spno_lower[hugosym]
        d['name'] = hugosym; tot_items +=1
        if d['name'] in dict_name_id:            
            if d['id'] not in dict_name_id[d['name']]:
                dict_name_id[d['name']].append(d['id'])
        else: dict_name_id[d['name']] = [d['id']]
        protname_var = generate_variants(hugosym)   # Get all variants of the protein name. protname_var may be empty
        bsyns = HGNC_dict_BasicOps_lower[hugosym]   # bsyns may be a [] list

        all_bsyns = []
        if len(bsyns) > 0:            
            for syn in bsyns:
                all_bsyns.append(syn);
                if len(syn.split()) > max_prot_len:
                    max_prot_len = len(syn.split())
                all_bsyns.extend(generate_variants(syn))
            
        if protname_var:
                all_bsyns.extend(protname_var)
        all_bsyns_uniq = list(set(all_bsyns))

        if len(all_bsyns_uniq) > 0:
            d['synonym'] = all_bsyns_uniq
            tot_items += len(all_bsyns_uniq)
            for sym in all_bsyns_uniq:
                if sym in dict_name_id:            
                    if d['id'] not in dict_name_id[sym]:
                        dict_name_id[sym].append(d['id'])
                else: dict_name_id[sym] = [d['id']]

        protein_names.append(d)

    print "Total no of protein names and synonyms extracted: ", tot_items
    print "Max length of protein name string: ", max_prot_len
    pickle.dump(dict_name_id, open("dictionary_nameToUniProtId.p", "wb"))
    pickle.dump(protein_names, open("list_UniProtId_name_synonym.p", "wb"))

    id_to_name_synonym = {}    
    for itm in protein_names:
        id_to_name_synonym[itm['id'].strip()] = {}
        id_to_name_synonym[itm['id'].strip()]['name'] = itm['name'].strip()
        if 'synonym' in itm: id_to_name_synonym[itm['id'].strip()]['synonym'] = itm['synonym'] 
    pickle.dump( id_to_name_synonym, open( "UniProtId_to_name_synonym.p", "wb" ) )

    # print "Creating Protein name to UniPROT Id dict and protein synonyms list. Now Exiting!!!"
    # exit()
    create_fst_BasicOps(protein_names)


def chkParens(string):
    match=re.search(r' \(.*?\)$', string)
    if match:                
        name = re.sub(r' \(.*?\)$',"",string).strip()                
    else:
        name = string.strip()
    return name


def create_fst_BasicOps(protein_names):
    stop_words_dict = pickle.load( open( "stop_words_dict.p", "rb" ) )
    a = fst.Acceptor(); syms = fst.SymbolTable(); ctr = 0
    for index, item in enumerate(protein_names):
        if len(item['name'].strip())==0: continue
        
        if check_for_unicode(item['name'].strip()):
            #print item['name'].strip()
            item['name'] = item['name'].decode("utf-8")

        if index == 0:
            a = fst.linear_chain(item['name'],syms=syms)
        else:
            b = fst.linear_chain(item['name'],syms=syms)
            a = a|b

        a.remove_epsilon(); ctr+=1

        if item.has_key('synonym'):
            for syn in item['synonym']:
                if check_for_unicode(syn.strip()): 
                    #print syn.strip()
                    syn = syn.decode("utf-8")
                if len(syn.strip().split()) > 1: is_stopword = False
                else:   is_stopword = check_stop_word(syn.strip(), stop_words_dict)          
                if len(syn.strip())==0 or is_stopword: continue
                c = fst.linear_chain(syn,syms=syms)
                a = a|c
                a.remove_epsilon()
                ctr+=1
   
        if index % 100 == 0:
            print "Index of protein name being processed: ", ctr

    print "Total no of protein names and synonyms processed: ", ctr     # ctr should be equal to tot_items
    det_a = a.determinize()
    det_a.minimize()
    det_a.arc_sort_input()
    det_a.write(os.path.join(parent_location, 'Minimized_Protein_Names_FST.fst'))
    syms.write(os.path.join(parent_location, 'SymbolTable_Protein_Names_FST.sym'))


def generate_variants(mystrg):
    ret_val = generate_unicodechars(mystrg)
    if ret_val:
        all_strg = [mystrg, ret_val[0]]   # strg is the original protein string, ret_val[0] contains unicode characters
    else:
        all_strg = [mystrg]

    for strg in all_strg:
        if "/" in strg:
            new_str = strg.replace("/", " ")
            ret_val.append(new_str)
            new_str = strg.replace("/", "")
            ret_val.append(new_str)

        m = re.search(r'\d+$', strg)
        # if the string ENDS in digits (or the entire string is made up of digits), m will be a Match object, or None otherwise.
        if m is not None:
            idx = strg.find(m.group())
            prec_char = strg[idx-1]
            if prec_char == " ":
                new_str = strg[:idx-1] + m.group()
                ret_val.append(new_str)
                new_str = strg[:idx-1] + '-' + m.group()
                ret_val.append(new_str)
            elif prec_char == "-":      # There will be duplicate variants for strings like: "IFN γ-2"
                new_str = strg[:idx-1] + m.group()
                ret_val.append(new_str)
                new_str = strg[:idx-1] + ' ' + m.group()
                ret_val.append(new_str)
            elif prec_char.isalpha():
                new_str = strg[:idx] + '-' + m.group()
                ret_val.append(new_str)
                new_str = strg[:idx] + ' ' + m.group()
                ret_val.append(new_str)

    ret_val_uniq = list(set(ret_val))
    return ret_val_uniq


def generate_unicodechars(strg):
    ret_val = []
    greek = False
    unicode_chars = [u"\u03b1", u"\u03b2", u"\u03b4", u"\u03b3"]
    tokens = word_tokenize(strg)
    tokens = tokenizeOnHyphen(tokens)

    if "|" in strg:
        print strg
        print "The character '|' is present in the above protein name. Please use a different character!"
        exit()

    for tok in tokens:
        if tok in ["alpha", "beta", "delta", "kappa", "gamma"]:
            greek = True
            if tok == "alpha":
                strg = replace_tokenstr(strg, "alpha", u"\u03b1")
            elif tok == "beta":
                strg = replace_tokenstr(strg, "beta", u"\u03b2")
            elif tok == "delta":
                strg = replace_tokenstr(strg, "delta", u"\u03b4")
            elif tok == "kappa":
                strg = replace_tokenstr(strg, "kappa", u"\u03ba")
            elif tok == "gamma":
                strg = replace_tokenstr(strg, "gamma", u"\u03b3")
    if greek:
        ret_val.append(strg)    # strg has unicode characters!
        orig_strg = strg
        if "-" in strg:         # Generate variants for "IFN-α" or "IFN γ-2 type II" or "IFN γ-2"
            strg_copy = strg
            hyphen_idx = findOccurences(strg, "-")
            for idx in hyphen_idx:
                if (idx > 0 and strg[idx - 1] in unicode_chars) or (idx < len(strg)-1 and strg[idx + 1] in unicode_chars):
                    strg = strg[:idx] + " " + strg[idx + 1:]
                    strg_copy = strg_copy[:idx] + "|" + strg_copy[idx + 1:]
            if strg != orig_strg:
                ret_val.append(strg)    # strg contains a blank space in place of "-"
                ret_val.append(strg_copy.replace("|", ""))  # strg_copy deletes the "-" before or after unicode chars

    return list(set(ret_val))


def findOccurences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


def replace_tokenstr(strg, tok, greek):
    indices = [(m.start(), m.end()) for m in re.finditer(tok, strg)]
    if len(indices) == 1:
        return strg.replace(tok, greek)
    else:   # The token exists multiple times in the string - strg
        for (start, end) in indices:
            if (start > 0 and strg[start - 1].isalnum()) or (end < len(strg) and strg[end].isalnum()):    # check whether the token is part of a word
                continue
            else:
                return strg[:start] + greek + strg[end:]    # replace the first occurrence of token as a separate word
        print "Control should never come here!"; exit()


def parse_HGNC_gene_proteins(HGNC_gene_proteins, HGNC_symbol_id):
    with open(os.path.join(txt_location, 'HGNC_gene_with_protein_product.txt'), "rb") as f:
        headings = f.readline(); filewidth = len(headings.split('\t')); 
        #next(f) # skip headings       
        csv_file = csv.reader(f,delimiter='\t'); ctr = 0
        for line in csv_file:  # line is a list of fields
            if len(line) != filewidth: 
                print "Mismatch in file width" 
                exit()                   
            hgnc_id = line[0]; symbol = line[1].lower(); name = line[2].lower()
            alias_sym = line[8].lower(); alias_name = line[9].lower(); #ctr += 1
            #print symbol, "\t\t Name: ", name, "\t\t Alias Symb: ", alias_sym, "\t\t Alias Name: ", alias_name
            #if ctr == 10:    exit(); 
            if "|" in name: print "Symbol: ", symbol, " contains a list in the name field: ", name
            if name.strip() == "": name_list = []
            else: name_list = [name.strip().lower()]
            if alias_sym.strip() != "":
                alias_sym_list = alias_sym.split('|')
                alias_sym_list = [syn.strip().lower() for syn in alias_sym_list]
            else:   alias_sym_list = []
            if alias_name.strip() != "":
                alias_name_list = alias_name.split('|')
                alias_name_list = [syn.strip().lower() for syn in alias_name_list]
            else:   alias_name_list = []
            synonyms = alias_sym_list + alias_name_list + name_list
            synonyms_uniq = list(set(synonyms))

            if symbol.strip().lower() in HGNC_gene_proteins:
                print "Inside parse_HGNC_gene_proteins"
                print symbol.strip()
                temp_list = HGNC_gene_proteins[symbol.strip().lower()]
                temp_list.extend(synonyms_uniq)                        
                temp_list_uniq = list(set(temp_list))
                HGNC_gene_proteins[symbol.strip().lower()] = temp_list_uniq                    
            else:
                HGNC_gene_proteins[symbol.strip().lower()] = synonyms_uniq

            if symbol.strip().lower() in HGNC_symbol_id:
                print "Inside parse_HGNC_gene_proteins: symbol to id"
                print symbol.strip()
            else:
                HGNC_symbol_id[symbol.strip().lower()] = hgnc_id.strip()    #hgnc_id should not be lowercase


def parse_BasicOps(filepath, HGNC_dict_BasicOps, HGNC_hugosym_spno):
    file_handle = open(filepath, "rb")
    basicops_file = json.load(file_handle)
    file_handle.close();
    
    for prot in basicops_file.keys():
        if (basicops_file[prot]['modname'] == "PROTEINOPS") and ('metadata' in basicops_file[prot]) and ('hugosym' in basicops_file[prot]['metadata']):            
            if 'synonyms' in basicops_file[prot]['metadata']: 
                if isinstance(basicops_file[prot]['metadata']['synonyms'], list):
                    synonyms = [syn.strip().lower() for syn in basicops_file[prot]['metadata']['synonyms']]                    
                else: 
                    synonyms = [basicops_file[prot]['metadata']['synonyms'].strip().lower()] # Only one synonym which is a string                                    
                
                synonyms = [syn.replace("splice variants:", "").strip() if "splice variants:" in syn else syn for syn in synonyms]
                synonyms_uniq = list(set(synonyms))
            else:   synonyms_uniq = []

            hugosym = basicops_file[prot]['metadata']['hugosym'] # This may be a list            
            if isinstance(hugosym, list):   # I checked, there is no "none" within a list
                #print hugosym, synonyms_uniq

                for each_hugosym in hugosym: 
                    # Get all variants of each_hugosym, append them to synonyms_uniq, then list(set(synonyms_uniq))               
                    if each_hugosym.strip().lower() in HGNC_dict_BasicOps:
                        temp_list = HGNC_dict_BasicOps[each_hugosym.strip().lower()]
                        temp_list.extend(synonyms_uniq)
                        temp_list_uniq = list(set(temp_list))
                        HGNC_dict_BasicOps[each_hugosym.strip().lower()] = temp_list_uniq                    
                    else:
                        HGNC_dict_BasicOps[each_hugosym.strip().lower()] = synonyms_uniq
                    HGNC_hugosym_spno[each_hugosym.strip().lower()] = basicops_file[prot]['metadata']['spnumber']                                       
            else:  
                if hugosym == 'none':
                    hgnc_key = basicops_file[prot]["name"].strip().lower()
                else:
                    hgnc_key = hugosym.strip().lower()

                # Get all variants of hgnc_key, append them to synonyms_uniq, then list(set(synonyms_uniq))
                if hgnc_key.lower() in HGNC_dict_BasicOps:
                    #print hgnc_key.strip()
                    temp_list = HGNC_dict_BasicOps[hgnc_key.lower()]
                    temp_list.extend(synonyms_uniq)                        
                    temp_list_uniq = list(set(temp_list))
                    HGNC_dict_BasicOps[hgnc_key.lower()] = temp_list_uniq                    
                else:
                    HGNC_dict_BasicOps[hgnc_key.lower()] = synonyms_uniq
                HGNC_hugosym_spno[hgnc_key.lower()] = basicops_file[prot]['metadata']['spnumber']


def parse_HGNC_gene_uniprot(filepath, HGNC_id_uniprot, HGNC_symbol_uniprot):
    with open(filepath, "rb") as f:
        headings = f.readline(); filewidth = len(headings.split('\t')); 
        #next(f) # skip headings       
        csv_file = csv.reader(f,delimiter='\t');
        for line in csv_file:  # line is a list of fields
            if len(line) != filewidth: 
                print "Mismatch in file width" 
                exit()
            uniprotid = line[-1].strip(); hgnc_id = line[0].strip(); symbol = line[1].lower().strip()
            
            if hgnc_id in HGNC_id_uniprot:
                print hgnc_id
            else:
                HGNC_id_uniprot[hgnc_id] = uniprotid
            if symbol in HGNC_symbol_uniprot:
                print hgnc_id, symbol
            else:
                HGNC_symbol_uniprot[symbol] = uniprotid


def check_stop_word(var, stop_words_dict):
    if var.strip() == "": return True
    lower_var = var.strip().lower()
    start_char = lower_var[0]        
    try:        
        if lower_var in stop_words_dict[start_char]:
            return True
        else:
            return False
    except KeyError:
        return False


if __name__ == "__main__":
    main_body()
