import argparse;
import sys;
import glob;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import copy;
from pprint import pprint;
import cPickle as pickle;
from Utilities import read_config_file


# python Prepare_Paragraphs_PdfArticles.py ./Total_PubMedIds_John_NoOpenAccess.txt ./protein_detect_output ./JSON_SENTENCE_ASSOCS 0
# python Prepare_Paragraphs_PdfArticles.py ./Total_PubMedIds_John_NoOpenAccess.txt ./protein_detect_output_Neighborhood1 ./JSON_SENTENCE_ASSOCS 1


def __main__():
    global myfile_data
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Prepare_Paragraphs', usage='Prepare_Paragraphs_PdfArticles.py <PubMedfilelist> <matchfiles_dir_path> <segmented_sent_dir_path> <windowSize>', description='Script to create paragraphs around sentences with protein matches')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('match_dir_path', help='Path where the protein match files are located')
    parser.add_argument('sentence_dir_path', help='Path where files containing the segmented sentences are located')
    parser.add_argument('windowSize', type=int, help='No of sentences to be included in the context')

    args = parser.parse_args()

    with open(os.path.join(txt_location, args.PubMedfilelist), 'rt') as f1:
        for pmid in f1:
            if not pmid.strip(): continue
            
            fname_string = '%s*' % (os.path.join(args.sentence_dir_path, pmid.strip()))
            fname = glob.glob(fname_string)            
            if len(fname) == 0:
                print "The segmented sentences file for PubMedId: ", pmid.strip(), "was not found!"
                continue                 
            myfile = fname[0].strip()

            file_handle = open(myfile, "rb")
            myfile_data = json.load(file_handle)    # myfile_data is a dictionary
            file_handle.close();
            
            if "sentenceList" not in myfile_data or "@items" not in myfile_data["sentenceList"]: 
                print "No sentences extracted for PubMedId: ", pmid.strip()
                continue;  
          
            last_sent_absidx = len(myfile_data["sentenceList"]["@items"]) - 1
            if myfile_data["sentenceList"]["@items"][-1]["index"] != last_sent_absidx:
                print "The last sentence index does not tally: ", last_sent_absidx, myfile_data["sentenceList"]["@items"][-1]["index"]
                exit()

            fname1_string = '%s*' % (os.path.join(args.match_dir_path, pmid.strip()))
            fname1 = glob.glob(fname1_string)
            if len(fname1) == 0:
                print "The match file for PubMedId: ", pmid.strip(), "was not found!"   # This statement will never get printed
                continue                 
            mymatchfile = fname1[0].strip()

            file_handle = open(mymatchfile, "rb")
            mymatchfile_data = json.load(file_handle)    # mymatchfile_data is a list of dictionaries
            file_handle.close()
            
            write_fname = pmid.strip() + '_Paragraphs.json'        
            fullname_write_fname = os.path.join(args.match_dir_path, write_fname)
            if os.path.exists(fullname_write_fname):
                print 'Skipping file: ', mymatchfile, ' as it has already been processed!'
                continue
            
            print 'Processing File: ', pmid.strip()            
            wSize = args.windowSize

            for match in mymatchfile_data:
                match_absId = int(match["id"].split("_")[1])
                before = getpreviousids(match_absId, wSize, pmid.strip())
                after = getnextids(match_absId, last_sent_absidx, wSize, pmid.strip())
                match["Neighborhood"] = before + [match["id"]] + after
                                                        
            json_fn = open(fullname_write_fname, 'w')
            json.dump(mymatchfile_data, json_fn, indent=4, ensure_ascii=False)
            json_fn.close()     
            #exit()   


def getpreviousids(match_absId, wSize, pmid):
    mylist = []
    for i in xrange(wSize, 0, -1):
        idx = match_absId - i
        if idx >= 0:
            mylist.append(pmid+"_"+str(idx))
    return mylist


def getnextids(match_absId, max_absidx, wSize, pmid):
    mylist = []
    for i in xrange(1, wSize+1, 1):
        idx = match_absId + i
        if idx <= max_absidx:
            mylist.insert(0, pmid+"_"+str(idx))
    return mylist      


__main__();
