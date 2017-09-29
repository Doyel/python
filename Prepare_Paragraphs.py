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


# python Prepare_Paragraphs.py ./Total_PubMedIds_John_NoOpenAccess.txt ./protein_detect_output ./JSON_SENTENCE_ASSOCS
# Take +2 and -2 around the current sentence


def __main__():
    global myfile_data
    parser = argparse.ArgumentParser(prog='Prepare_Paragraphs', usage='Prepare_Paragraphs.py <PubMedfilelist> <matchfiles_dir_path> <segmented_sent_dir_path>', description='Script to create paragraphs around sentences with protein matches')
    parser.add_argument('PubMedfilelist', help='File listing the PubMed ids to be processed')
    parser.add_argument('match_dir_path', help='Path where the protein match files are located')
    parser.add_argument('sentence_dir_path', help='Path where files containing the segmented sentences are located')

    args = parser.parse_args()

    with open(args.PubMedfilelist, 'rt') as f1:
        for pmid in f1:
            if not pmid.strip(): continue
            
            fname_string = '%s' % (os.path.join(args.sentence_dir_path, pmid.strip()))
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
            end_sent_idx = len(myfile_data["sentenceList"]["@items"]) - 1
            
            fname1_string = '%s*' % (os.path.join(args.match_dir_path, pmid.strip()))
            fname1 = glob.glob(fname1_string)
            if len(fname1) == 0:
                print "The match file for PubMedId: ", pmid.strip(), "was not found!"   # This statement will never get printed
                continue                 
            mymatchfile = fname1[0].strip()

            file_handle = open(mymatchfile, "rb")
            mymatchfile_data = json.load(file_handle)    # myfile_data is a list of dictionaries
            file_handle.close();
            
            write_fname = pmid.strip() + '_Paragraphs.json'        
            fullname_write_fname = os.path.join(args.match_dir_path, write_fname)
            if os.path.exists(fullname_write_fname):
                print 'Skipping file: ', mymatchfile, ' as it has already been processed!'
                continue
            print 'Processing File: ', pmid.strip()

            myjson = []; plus_q = []; minus_q = []                        
            for itm in myfile_data["sentenceList"]["@items"]:
                if len(mymatchfile_data) == 0:  break;
                curr_idx = itm["index"]                
                curr_match = int(mymatchfile_data[0]["id"].split("_")[1])
                update_q(plus_q, minus_q, curr_idx, end_sent_idx)
                if curr_idx == curr_match:
                    temp_dict = copy.deepcopy(mymatchfile_data[0])
                    update_dict(temp_dict, plus_q, minus_q, pmid.strip())
                    myjson.append(temp_dict); mymatchfile_data.pop(0)
                    
            json_fn = open(fullname_write_fname, 'w')
            json.dump(myjson, json_fn, indent=4, ensure_ascii=False)
            json_fn.close()     
            #exit()   


def update_q(plus_q, minus_q, curr_idx, end_idx):
    if curr_idx == 0:
        for i in xrange(1, 3, 1):
            idx = curr_idx + i
            if idx <= end_idx:                  
                plus_q.append(copy.deepcopy(myfile_data["sentenceList"]["@items"][idx]))
    else:
        pidx = curr_idx + 2
        if pidx <= end_idx: 
            if len(plus_q) == 2:
                plus_q.pop(0)
            plus_q.append(copy.deepcopy(myfile_data["sentenceList"]["@items"][pidx]))
        else:
            if len(plus_q) > 0: plus_q.pop(0)                
        
    nidx = curr_idx - 1
    if nidx >= 0: 
        if len(minus_q) == 2:
            minus_q.pop(0)
        minus_q.append(copy.deepcopy(myfile_data["sentenceList"]["@items"][nidx]))
             

def update_dict(temp_dict, plus_q, minus_q, pmid):
    all_sents = []; all_idx = []
    for elem in minus_q:
        all_sents.append(elem["sentenceText"])
        all_idx.append(pmid + "_" + str(elem["index"]))
    
    all_sents.append(temp_dict["text"])
    all_idx.append(temp_dict["id"])
    
    for elem in plus_q:
        all_sents.append(elem["sentenceText"])
        all_idx.append(pmid + "_" + str(elem["index"]))   

    temp_dict["paragraph"] = " ".join(all_sents)
    temp_dict["para_ids"] = list(all_idx)


__main__();
