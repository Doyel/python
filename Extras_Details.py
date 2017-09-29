# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
from pprint import pprint;
import os;
import json;

#python Extras_Details.py ./datums_03_2016.json

def main_body():
    parser = argparse.ArgumentParser(prog='Extras Details', usage='Extras_Details.py <datums_JSON_file>', description='Script to get details about the Extras fields')
    parser.add_argument('datums_file', help='newest version of the datums file given by Mark')    
    args = parser.parse_args()

    file_handle = open(args.datums_file, "rb")
    myfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
    file_handle.close()
    
    tot_datums = 0; tot_datums_extras = 0; extras_reqs_or_dnreqs = 0; all_pmids = []
    extras_reqs = 0; extras_dnreqs = 0; pmids = []; extras_reqs_dnreqs = 0
    subject_extra_entity_cnt = 0; subject_prot = ""
    for elem in myfile_data:    # elem represents a datum and is a dictionary
        if "subject" in elem and "entity" in elem["subject"]:
            subject_prot = elem["subject"]["entity"].lower()
        tot_datums += 1; all_pmids.append(elem["source"]["pmid"].strip())

        if 'extras' in elem and len(elem["extras"]) > 0:
            tot_datums_extras += 1; same_subject_extra_entity = True
            thisdatum_reqs = False; thisdatum_dnreqs = False
            for ext in elem["extras"]:  # elem["extras"] is a list of dictionaries
                if ext["type"] == "reqs":
                    thisdatum_reqs = True
                if ext["type"] == "does not req":              
                    thisdatum_dnreqs = True

                if len(subject_prot) > 0 and "entities" in ext and (thisdatum_reqs or thisdatum_dnreqs):
                    for ent in ext["entities"]:  # ext["entities"] is a list and may have more than one element in it
                        ent = ent.lower()
                        if ent != subject_prot:
                            same_subject_extra_entity = False
                        else:
                            print "PMID: ", elem["source"]["pmid"].strip(), "\t", subject_prot, "\t", ent

            if thisdatum_reqs or thisdatum_dnreqs:
                pmids.append(elem["source"]["pmid"].strip()); extras_reqs_or_dnreqs += 1
            if thisdatum_reqs and thisdatum_dnreqs:
                extras_reqs_dnreqs += 1
            if thisdatum_reqs and not(thisdatum_dnreqs):
                extras_reqs += 1
            if not(thisdatum_reqs) and thisdatum_dnreqs:
                extras_dnreqs += 1

            # There may be some datums containing multiple extras datums in which some of the extras datums could have same subject protein and extras entities protein
            if (thisdatum_reqs or thisdatum_dnreqs) and same_subject_extra_entity:
                subject_extra_entity_cnt += 1


    print "The no. of total datums in JSON file: ", tot_datums
    print "The total no. of datums that have the extras field: ", tot_datums_extras
    print "The total no. of datums having extras field of type reqs OR not reqs: ", extras_reqs_or_dnreqs
    print "The total no. of datums having extras field of type reqs AND not reqs: ", extras_reqs_dnreqs
    print "The total no. of datums having extras field of type - requires, ONLY: ", extras_reqs
    print "The total no. of datums having extras field of type - does not require, ONLY: ", extras_dnreqs

    print "The total no. of datums in which the subject entity is the same as all extras entities: ", subject_extra_entity_cnt

    pmids = list(set(pmids))
    print "The total no. of articles having extras field of type reqs OR not reqs: ", len(pmids)

    all_pmids = list(set(all_pmids))
    print "The total no. of articles present in the Datum KB: ", len(all_pmids)

    # pmid_fn1 = open("PubMedIds_with_Extras.txt", 'w');
    # for pm in pmids:
    #     pmid_fn1.write(pm+'\n')
    # pmid_fn1.close()
    #
    # pmid_fn2 = open("Total_PubMedIds_DatumKB.txt", 'w');
    # for pm in all_pmids:
    #     pmid_fn2.write(pm+'\n')
    # pmid_fn2.close()


if __name__ == "__main__":
    main_body()
