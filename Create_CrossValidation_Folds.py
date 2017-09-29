import argparse;
import os
import shutil
import subprocess
import cPickle as pickle;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
from random import shuffle


# python Create_CrossValidation_Folds.py vw_RNAi_Weights_1_1_Training_File.txt ./datafiles_TriggerHighlights_BaseLine2/RNAi/ 10 folds_1_1
# python Create_CrossValidation_Folds.py vw_RNAi_WithFeedback_Training_File_1_1.txt ./user_feedback_JSON/Apr4_2017/RNAi/ 10 folds_1_1


def main_body():
    parser = argparse.ArgumentParser(prog='Create_CrossValidation_Folds', usage='Create_CrossValidation_Folds.py <Training_File> <workdir> <num_folds> <folds_dirname>', description='Script to divide a training file into a certain no.of folds')
    parser.add_argument('trainFile', help='Training File containing training instances from all articles')
    parser.add_argument('workdir', help='Working dir where all the data file folds will be created')
    parser.add_argument('num_folds', help='Number of folds we want to create')
    parser.add_argument('folds_dirname', help='Name of folder that will contain the internal CV folds of training file')

    args = parser.parse_args()

    posrecords = []; negrecords = []
    extract_records(posrecords, negrecords, os.path.join(args.workdir, args.trainFile))
    total_recs = len(posrecords) + len(negrecords)
    shuffle(posrecords); shuffle(negrecords)
    ensure_dir(os.path.join(args.workdir, args.folds_dirname))
    fold_dir = os.path.join(args.workdir, args.folds_dirname)

    print "Separating the dataset into ", args.num_folds, " cross validation folds"
    folds = []
    create_folds(posrecords, folds, int(args.num_folds))
    create_folds(negrecords, folds, int(args.num_folds))
    check_foldsize(folds, total_recs)

    print "Creating train and test files for each cross validation fold"
    for i, elem in enumerate(folds):
        shuffle(elem)
        testfoldfile = args.trainFile.split(".")[0] + "_TestFold-" + str(i) + "." + args.trainFile.split(".")[1]
        writefile(os.path.join(fold_dir, testfoldfile), elem)
        trainfold_recs = []
        for j, myelem in enumerate(folds):
            if j == i:
                continue
            trainfold_recs.extend(myelem)
        shuffle(trainfold_recs)
        if len(elem) + len(trainfold_recs) != total_recs:
            print "Something went wrong while creating the train and test files for each fold"
            exit()
        trainfoldfile = args.trainFile.split(".")[0] + "_TrainFold-" + str(i) + "." + args.trainFile.split(".")[1]
        writefile(os.path.join(fold_dir, trainfoldfile), trainfold_recs)


def extract_records(posrecords, negrecords, trainfilename):
    with open(trainfilename, 'rt') as f:
        print "Processing file: ", os.path.basename(trainfilename)
        for train_instance in f:
            train_instance = train_instance.strip()
            currlabel = train_instance.split()[0].strip()
            if currlabel == "1":
                posrecords.append(train_instance)
            elif currlabel == "-1":
                negrecords.append(train_instance)


def writefile(filepath, records):
    data_fn = open(filepath, 'w')
    for rc in records:
        data_fn.write(rc + '\n')
    data_fn.close()


def ensure_dir(file_path):
    if os.path.exists(file_path):
        shutil.rmtree(file_path)
    os.mkdir(file_path)


def create_folds(myrecords, folds, num_folds):
    fold_size = len(myrecords) / num_folds
    rem = len(myrecords) % num_folds
    if len(folds) == num_folds:
        for i in xrange(num_folds):
            folds[i].extend(myrecords[0:fold_size])
            myrecords = myrecords[fold_size:]
    else:
        for i in xrange(num_folds):
            folds.append(myrecords[0:fold_size])
            myrecords = myrecords[fold_size:]

    if len(myrecords) != rem:
        print "Something went wrong while trying to separate out the folds"
        exit()
    for i in xrange(rem):
        folds[i].append(myrecords[i])


def check_foldsize(folds, total_recs):
    mysum = 0
    for fold in folds:
        mysum += len(fold)
    if mysum != total_recs:
        print "Fold sizes do not tally up. The folds couldnt get created correctly!"
        exit()


if __name__ == "__main__":
    main_body()