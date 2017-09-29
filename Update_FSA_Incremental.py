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
import datetime;
import string
from pprint import pprint;
import cPickle as pickle;
import ntpath
import fst;
from Utilities import read_config_file, check_for_unicode
from Protein_Extractor_UniProtID import parse_BasicOps, parse_HGNC_gene_uniprot


# python Update_FSA_Incremental.py Minimized_Protein_Names_FST.fst SymbolTable_Protein_Names_FST.sym  incremental_basicops.json


def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Update_FSA_Incremental', usage='Update_FSA_Incremental.py <FST file> <SymbolTable file> <incremental_basicOps>', description='Script to include synonyms for proteins in our FSA on a one-off basis')
    parser.add_argument('fst_file', help='Name of the binary FST file')
    parser.add_argument('symboltable_file', help='Name of the binary symbol table file')
    parser.add_argument('newProt_file', help='Proteins that need to be included in the FSA')

    args = parser.parse_args()
    fst_folder = os.path.join(parent_location, "FST_files")

    dict_name_id = pickle.load(open("dictionary_nameToUniProtId.p", 'rb'))
    protein_names = pickle.load(open("list_UniProtId_name_synonym.p", 'rb'))

    HGNC_dict_BasicOps_lower = {}; HGNC_hugosym_spno_lower = {}
    # HGNC_dict_BasicOps will have all hugosyms as keys
    parse_BasicOps(os.path.join(json_location, args.newProt_file), HGNC_dict_BasicOps_lower, HGNC_hugosym_spno_lower)

    HGNC_hgncid_uniprotid = {}; HGNC_symbol_uniprot_lower = {}
    parse_HGNC_gene_uniprot(os.path.join(txt_location, 'HGNC_gene_with_UNIPROT_IDs.txt'), HGNC_hgncid_uniprotid, HGNC_symbol_uniprot_lower)

    new_protein_names = []
    for hugosym in HGNC_dict_BasicOps_lower:
        d = {}
        if hugosym in HGNC_symbol_uniprot_lower:
            d['id'] = HGNC_symbol_uniprot_lower[hugosym]
        else:
            d['id'] = HGNC_hugosym_spno_lower[hugosym]

        d['name'] = hugosym
        if d['name'] in dict_name_id:
            if d['id'] not in dict_name_id[d['name']]:
                print "Adding UniprotID: " + d['id'] + " for existing protein: " + d['name']
                dict_name_id[d['name']].append(d['id'])
        else:
            print "Creating new entry with key as protein: " + d['name'] + " having UniprotID: " + d['id']
            dict_name_id[d['name']] = [d['id']]

        bsyns = HGNC_dict_BasicOps_lower[hugosym]  # bsyns may be a [] list
        if len(bsyns) > 0:
            bsyns_uniq = list(set(bsyns))
            d['synonym'] = bsyns_uniq
            for sym in bsyns_uniq:
                if sym in dict_name_id:
                    if d['id'] not in dict_name_id[sym]:
                        print "Adding UniprotID: " + d['id'] + " for existing synonym: " + sym
                        dict_name_id[sym].append(d['id'])
                else:
                    print "Creating new entry with key as synonym: " + sym + " having UniprotID: " + d['id']
                    dict_name_id[sym] = [d['id']]

        new_protein_names.append(d)

    fsa_names_include = []
    for itm in protein_names:
        for newitm in new_protein_names:
            if itm['id'].strip() == newitm['id'] and itm['name'] == newitm['name']:
                print "Found entry for protein: ", newitm['name']
                new_syns = set(newitm['synonym']) - set(itm['synonym'])
                itm['synonym'].extend(new_syns)
                fsa_names_include.extend(list(new_syns))
                print "Updated entry for protein: ", itm['name']
                pprint(itm)
    print "Strings that need to be included in the FSA: ", ", ".join(fsa_names_include)

    id_to_name_synonym = {}
    for itm in protein_names:
        id_to_name_synonym[itm['id'].strip()] = {}
        id_to_name_synonym[itm['id'].strip()]['name'] = itm['name'].strip()
        if 'synonym' in itm: id_to_name_synonym[itm['id'].strip()]['synonym'] = itm['synonym']
    pickle.dump(id_to_name_synonym, open("UniProtId_to_name_synonym.p", "wb"))

    pickle.dump(dict_name_id, open("dictionary_nameToUniProtId.p", "wb"))
    pickle.dump(protein_names, open("list_UniProtId_name_synonym.p", "wb"))

    update_fst(fsa_names_include, fst_folder, args.fst_file, args.symboltable_file)


def update_fst(fsa_names_include, fst_folder, fst_file, symboltable_file):
    my_fst = fst.read_std(os.path.join(fst_folder, fst_file))
    syms = fst.read_symbols(os.path.join(fst_folder, symboltable_file))

    for prot in fsa_names_include:
        if check_for_unicode(prot.strip()):
            prot = prot.decode("utf-8")
        b = fst.linear_chain(prot, syms=syms)
        my_fst = my_fst | b

        my_fst.remove_epsilon()

    det_a = my_fst.determinize()
    det_a.minimize()
    det_a.arc_sort_input()
    save_fst(det_a, os.path.join(parent_location, fst_folder, 'Minimized_Protein_Names_FST.fst'))
    save_fst(syms, os.path.join(parent_location, fst_folder, 'SymbolTable_Protein_Names_FST.sym'))
    # det_a.write(os.path.join(parent_location, fst_folder, 'Minimized_Protein_Names_FST.fst'))
    # syms.write(os.path.join(parent_location, fst_folder, 'SymbolTable_Protein_Names_FST.sym'))


def save_fst(structure, fileloc):
    fname = ntpath.basename(fileloc).strip()
    dirname = ntpath.dirname(fileloc)
    datestr = "_" + datetime.datetime.today().strftime('%Y%m%d') + "."
    if os.path.isfile(fileloc):
        print "A file having the same name as: " + fname + " already exists at the path: ", dirname
        print "Making a copy of the original file: ", fname
        fullpath_renamed_passagedictfile = os.path.join(dirname, fname.split(".")[0] + datestr + fname.split(".")[1])
        os.rename(fileloc, fullpath_renamed_passagedictfile)
        print "Saving the file: ", fname
        structure.write(fileloc)


if __name__ == "__main__":
    main_body()
