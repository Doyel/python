__author__ = 'Doyel'

import argparse
from xml.etree import ElementTree
from pprint import pprint;
import cPickle as pickle;
import os;
import glob

# python eSummary_Parser_OpenAccess.py ./eSummary_extras ./PubMedIDS_with_Extras_OpenAccess.txt

parser = argparse.ArgumentParser(prog='eSummary XML Parser', usage='eSummary_Parser_OpenAccess.py <inputfile>', description='Script to parse the eSummary file for an article and extract out its pmid if it exists')
parser.add_argument('dir_path', help='Path where the eSummary files are located')
parser.add_argument('PMIDfilelist', help='File listing the PMIDs to be processed')

args = parser.parse_args()
#print args.XMLinputfile


summary_dict= {}; pmid_pmcid = {}; num = 0; final_pmcid_list = []
with open(args.PMIDfilelist, 'rt') as f1:
    for pmid in f1:
        if not pmid.strip(): continue
        curr_pmid = ""
        eSumfname = "eSummary_Results_" + pmid.strip()
        fname_string = '%s*' % (os.path.join(args.dir_path, eSumfname))
        fname = glob.glob(fname_string)
        if len(fname) == 0:
            print "The specified JSON file was not found: ", pmid.strip()
            continue                 
        myfile = fname[0].strip()
        xml_tree = ElementTree.parse(myfile)
        root = xml_tree.getroot()
        for each_doc_sum in root:
            for eachitem in each_doc_sum:            
                if eachitem.tag == "Id":
                    curr_pmid = eachitem.text 
                    if curr_pmid in summary_dict:
                        print "pmid already present in summary_dict: ", curr_pmid
                        print "Please check for duplicate pmids in: ", args.PMIDfilelist
                        exit()               
                    summary_dict[curr_pmid] = {}
                if eachitem.get('Name') == "PubDate":   summary_dict[curr_pmid]["PubDate"] = eachitem.text.split('/')[0]            
                if eachitem.get('Name') == "FullJournalName":   summary_dict[curr_pmid]["FullJournalName"] = eachitem.text            
                if eachitem.get('Name') == "Title":   summary_dict[curr_pmid]["Title"] = eachitem.text            
                if eachitem.get('Name') == "Volume":   summary_dict[curr_pmid]["Volume"] = eachitem.text            
                if eachitem.get('Name') == "Issue":   summary_dict[curr_pmid]["Issue"] = eachitem.text            
                if eachitem.get('Name') == "Pages":   summary_dict[curr_pmid]["Pages"] = eachitem.text            
                if eachitem.get('Name') == "AuthorList":
                    authors = ""
                    for eachauthor in eachitem:
                         if eachauthor.get('Name') == "Author": authors += eachauthor.text + ", "
                    summary_dict[curr_pmid]["Authors"] = authors[:-2]        
                if eachitem.get('Name') == "ArticleIds":
                    mypmc = ""
                    for eachid in eachitem:
                        if eachid.get('Name') == "pubmed": mypmid = eachid.text
                        if eachid.get('Name') == "pmc": mypmc = eachid.text
                    if mypmc == "": mypmc = "Does Not Exist";   
                    else: num += 1; final_pmcid_list.append(mypmc.strip())                                                          
                    if curr_pmid == mypmid:
                        summary_dict[curr_pmid]["PMCID"] = mypmc
                        summary_dict[curr_pmid]["PMID"] = mypmid
                        pmid_pmcid[curr_pmid] = mypmc
                    else:
                        print "There is something wrong in eSummary output for pmid:", curr_pmid
                        exit()
                    
                #print eachitem.get('Name')
            temp_dict = summary_dict[curr_pmid] 
            if "PubDate" not in temp_dict: print "For pmid: ", curr_pmid, " there is no PubDate attribute"   
            if "FullJournalName" not in temp_dict: print "For pmid: ", curr_pmid, " there is no FullJournalName attribute"
            if "Title" not in temp_dict: print "For pmid: ", curr_pmid, " there is no Title attribute"
            if "Volume" not in temp_dict: print "For pmid: ", curr_pmid, " there is no Volume attribute"
            if "Issue" not in temp_dict: print "For pmid: ", curr_pmid, " there is no Issue attribute"
            if "Pages" not in temp_dict: print "For pmid: ", curr_pmid, " there is no Pages attribute"
            if "Authors" not in temp_dict: print "For pmid: ", curr_pmid, " there is no Authors attribute"
            if "PMCID" not in temp_dict: print "For pmid: ", curr_pmid, " there is no PMCID attribute"
            if "PMID" not in temp_dict: print "For pmid: ", curr_pmid, " there is no PMID attribute"

#pprint(summary_dict)
pickle.dump( summary_dict, open( "summary_dict_Extras_OpenAccess.p", "wb" ) )
print "No. of articles in summary_dict: ", len(summary_dict.keys())
 

