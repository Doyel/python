import argparse;
import os
import shutil
import subprocess
import cPickle as pickle;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
from random import shuffle


# python CrossValidation_ParamTuning.py RNAi vw_RNAi_Training_File_1_1.txt ./datafiles_TriggerHighlights_BaseLine2/RNAi 10 bfgs


def main_body():
    parser = argparse.ArgumentParser(prog='CrossValidation_ParamTuning', usage='CrossValidation_ParamTuning.py <label> <Training_File> <workdir> <num_folds> <opt_method>', description='Script to tune parameters and return the best values for the params')
    parser.add_argument('label', help='Predicate being analyzed')
    parser.add_argument('trainFile', help='Training File containing training instances from all articles')
    parser.add_argument('workdir', help='Working dir where all the data file folds will be created')
    parser.add_argument('num_folds', help='Number of folds we want to create')
    parser.add_argument('optim_method', help='The optimization method for VW')

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()
    label = args.label

    DatumKB_dict = pickle.load(open("pubmedid_extras.p", "rb"))

    if args.optim_method == "bfgs":
        method = " --bfgs --mem 20"
    elif args.optim_method == "ftrl":
        method = " --ftrl --holdout_off"
    else:
        print "Cannot interpret the provided optimization method: ", args.optim_method
        print "Only values allowed currently are \'bfgs\' and \'ftrl\'!"
        exit()

    posrecords = []; negrecords = []
    extract_records(posrecords, negrecords, os.path.join(args.workdir, args.trainFile))
    total_recs = len(posrecords) + len(negrecords)
    shuffle(posrecords); shuffle(negrecords)
    ensure_dir(os.path.join(args.workdir, "folds"))
    fold_dir = os.path.join(args.workdir, "folds")

    print "Separating the dataset into ", args.num_folds, " cross validation folds"
    folds = []
    create_folds(posrecords, folds, int(args.num_folds))
    create_folds(negrecords, folds, int(args.num_folds))
    check_foldsize(folds, total_recs)

    print "Creating train and test files for each cross validation fold"
    for i, elem in enumerate(folds):
        shuffle(elem)
        testfoldfile = args.trainFile.split(".")[0] + "_TestFold" + str(i) + "." + args.trainFile.split(".")[1]
        writefile(os.path.join(fold_dir, testfoldfile), elem)
        trainfold_recs = []
        for j, myelem in enumerate(folds):
            if j == i:
                continue
            trainfold_recs.extend(myelem)
        shuffle(trainfold_recs)
        trainfoldfile = args.trainFile.split(".")[0] + "_TrainFold" + str(i) + "." + args.trainFile.split(".")[1]
        writefile(os.path.join(fold_dir, trainfoldfile), trainfold_recs)

    L1_lambda_vals = [0, 0.1, 0.25, 0.5, 0.75, 0.01, 0.05, 0.001, 0.005, 2, 5, 10, 12, 15, 20, 50, 75, 100, 150, 250, 500]
    L2_lambda_vals = [0, 0.1, 0.25, 0.5, 0.75, 0.01, 0.05, 0.001, 0.005, 2, 5, 10, 12, 15, 20, 50, 75, 100, 150, 250, 500]
    ftrl_alpha_vals = [0, 0.1, 0.25, 0.5, 0.75, 0.01, 0.05, 0.001, 0.005, 2, 5, 10, 12, 15, 20, 50, 75, 100, 150, 250, 500]
    ftrl_beta_vals = [0, 0.1, 0.25, 0.5, 0.75, 0.01, 0.05, 0.001, 0.005, 2, 5, 10, 12, 15, 20, 50, 75, 100, 150, 250, 500]

    only_L1 = False; only_L2 = False; L1_L2 = False; ftrl_params = False; all_four = False


    ensure_dir(os.path.join(fold_dir, "models_folds"))
    models_dir = os.path.join(fold_dir, "models_folds")
    for l1_val in L1_lambda_vals:

        l1_str = "_L1-" + getstr(l1_val); l1_cmd = " --l1 " + str(l1_val)
        combined_predictions_file = os.path.join(models_dir, label + "_AllPredictions" + l1_str + ".txt")
        combined_confidences_file = os.path.join(models_dir, label + "_AllConfidences" + l1_str + ".txt")
        PRpoints_file = os.path.join(models_dir, label + "_PRPoints" + l1_str + ".csv")
        all_predictions = []; all_confidences = []

        for i in xrange(int(args.num_folds)):
            testfoldfile = args.trainFile.split(".")[0] + "_TestFold" + str(i) + "." + args.trainFile.split(".")[1]
            trainfoldfile = args.trainFile.split(".")[0] + "_TrainFold" + str(i) + "." + args.trainFile.split(".")[1]
            model_file = os.path.join(models_dir, label + "_Fold" + str(i) + "_VWModel" + l1_str + ".vw")
            predictions_file = os.path.join(models_dir, label + "_Fold" + str(i) + "_Predictions" + l1_str + ".txt")
            confidences_file = os.path.join(models_dir, label + "_Fold" + str(i) + "_Confidences" + l1_str + ".txt")

            cmd = "vw -d " + os.path.join(fold_dir, trainfoldfile) + " -c --passes 100 -f " + model_file + " --loss_function logistic" + method + l1_cmd
            run_command(cmd)

            cmd = "vw -d " + os.path.join(fold_dir, testfoldfile) + " -t -i " + model_file + " -p /dev/stdout | /ua/ml-group/big-mechanism-project/vowpal_wabbit/utl/logistic -0 > " + predictions_file
            run_command(cmd)

            cmd = "python Generate_Confidences.py " + label + " " + predictions_file + " " + confidences_file + " ./pubmedid_extras.p"
            run_command(cmd)

            readfile(all_predictions, predictions_file)
            readfile(all_confidences, confidences_file)

        writefile(combined_predictions_file, all_predictions)
        writefile(combined_confidences_file, all_confidences)
        label_recall = get_recall_denominator(combined_predictions_file, DatumKB_dict, label)

        cmd = "java -cp $PR_CURVE_JAR LRTest_CLArg " + combined_confidences_file + " " + PRpoints_file + " " + label_recall
        run_command(cmd)
        exit()


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


def getstr(myval):
    mystr = str(myval)
    if "." in mystr:
        mystr = mystr.replace(".", "")
    return mystr


def run_command(cmd):
    print "-----------------------------------------------------------------------------------------------------------"
    print cmd
    print "-----------------------------------------------------------------------------------------------------------"
    my_env = os.environ.copy()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

    for line in p.stderr.readlines():  # Do NOT comment this for statement
        print line.strip()
        # pass


def writefile(filepath, records):
    data_fn = open(filepath, 'w')
    for rc in records:
        data_fn.write(rc + '\n')
    data_fn.close()


def readfile(records, filepath):
    with open(filepath, 'rt') as f1:
        for line in f1:
            records.append(line.strip())


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


def get_recall_denominator(pfile, DatumKB_dict, label):
    test_pmids = []
    with open(pfile, 'rt') as f1:
        for line in f1:
            line_list = line.split()
            pmid = line_list[1].strip().split("_", 1)[0]
            test_pmids.append(pmid.strip())
    test_pmids = list(set(test_pmids))

    temp_dict = {}  # temp_dict is a flattened version of mydict. It only contains those PMIDs that have extras datums of interest and are included in the Test file list.
    create_act_dict(test_pmids, temp_dict, DatumKB_dict)
    ctr = 0
    for pmid in temp_dict.keys():
        for prot in temp_dict[pmid]:
            if label.strip() in temp_dict[pmid][prot]:
                ctr += 1
    # Use this value as the denominator for Recall, not the positives present in the test file
    print "Actual positives in the training file for label", label.strip(), "with respect to the Datum KB: ", ctr
    return str(ctr)


def create_act_dict(mypmids, tdict, mydict):
    print "Total no of pmids in the training file: ", len(mypmids)
    for pmid in mypmids:
        if pmid in mydict:
            tdict[pmid] = {}          # tdict contains only those PMIDs that exist both in DatumKB and the test file
            for ent in mydict[pmid]:
                mylist = []
                if len(mydict[pmid][ent]['reqstest']) > 0:
                    mylist.append("reqs"); mylist.extend(mydict[pmid][ent]["reqstest"])
                if len(mydict[pmid][ent]['dnreqstest']) > 0:
                    mylist.append("dnreqs"); mylist.extend(mydict[pmid][ent]["dnreqstest"])
                tdict[pmid][ent] = list(set(mylist))


if __name__ == "__main__":
    main_body()