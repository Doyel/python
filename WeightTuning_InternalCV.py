# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import os
import ntpath
import json
import shutil
import cPickle as pickle
import sys
reload(sys)
sys.setdefaultencoding("utf-8")
from random import shuffle
from Create_MLFile import create_lexical_features
from Create_CrossValidation_Folds import create_folds, check_foldsize
from ParamTuning_InternalCV import getstr, getFile, run_command, get_actual_pos_neg, get_AUPR, readfile, writefile, clean_files


# First tune parameters, then tune feedback weights!!! Do NOT tune Baseline weights!
# python WeightTuning_InternalCV.py RNAi ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/datafiles_Baseline2 10 bfgs
# python WeightTuning_InternalCV.py RNAi ./datafiles/datafiles_NoMutualExc_NoProtFeat_NoFreqPrune_Neighborhood0/datafiles_Baseline2 10 bfgs --feedback ./user_feedback_JSON/Apr4_2017/
# python WeightTuning_InternalCV.py KO /ua/ml-group/big-mechanism-project/PLEIADES/Sep2016_SMA/datafiles_TriggerHighlights_BaseLine2 10 bfgs
# python WeightTuning_InternalCV.py omission /ua/ml-group/big-mechanism-project/PLEIADES/Sep2016_SMA/datafiles_TriggerHighlights_BaseLine2 10 bfgs


def main_body():
    parser = argparse.ArgumentParser(prog='WeightTuning_InternalCV', usage='WeightTuning_InternalCV.py <label> <passageDict_dir> <num_folds> <optim_method> [--feedback <Path to feedback passage dict>]', description='Script to tune instance weights')
    parser.add_argument('label', help='Predicate being analyzed')
    parser.add_argument('passageDict_dir', help='Dir where all the original passage dicts are located')
    parser.add_argument('num_folds', help='Number of folds we want to create')
    parser.add_argument('optim_method', help='The optimization method for VW')
    parser.add_argument('-f', '--feedback', help='Path to the feedback passage dict json file')  # Optional arg

    args = parser.parse_args()
    label = args.label
    Instance_Weights = [1, 5, 10, 20, 30, 50, 75, 100, 125, 150]
    listpassagedicts = [os.path.join(args.passageDict_dir, "train_passage_dict.json")]

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    DatumKB_dict = pickle.load(open("pubmedid_extras.p", "rb"))

    if args.optim_method == "bfgs":
        method = " --bfgs --mem 20" # BFGS with holdout_off is a joke! Don't ever try the --holdout_off switch with BFGS
    elif args.optim_method == "ftrl":
        method = " --ftrl --holdout_off"
    elif args.optim_method == "default":
        method = ""
    else:
        print "Cannot interpret the provided optimization method: ", args.optim_method
        print "Values allowed currently are \'bfgs\' and \'ftrl\'!"
        exit()

    if args.feedback is None:       # Should NOT tune BaseLine Weights.
        OA_Multiples = [1]          # Acc. to my experiments, the AUC for various absolute weights is more or less similar for either optimization method!
        listpassagedicts.append(os.path.join(args.passageDict_dir, "openaccess_passage_dict.json"))
        folder_prefix = "Weights_"
        work_dir = create_workDir(args.passageDict_dir, label, args.optim_method)
    else:
        OA_Multiples = [1, 2, 3, 5, 10, 15, 20]
        listpassagedicts.append(os.path.join(args.feedback, "openaccess_passage_dict_WithFeedback.json"))
        folder_prefix = "Feedback_Weights_"
        work_dir = create_workDir(args.feedback, label, args.optim_method)

    for i, elem in enumerate(listpassagedicts):
        if os.path.isfile(elem):
            print "Found - " + ntpath.basename(elem)
            file_handle = open(elem, "rb")
            if i == 0:
                train_passage_dict = json.load(file_handle)
            else:
                OA_passage_dict = json.load(file_handle)
            file_handle.close()
        else:
            print "The passage dict file - " + ntpath.basename(elem) + " was not found in " + ntpath.dirname(elem)
            exit()

    best_weights_auc = (0.0, 0.0, 0.0)
    for wt in Instance_Weights:
        print "Non OpenAccess weight: ", wt
        update_passagedict(label, train_passage_dict, wt)
        train_instances = create_instances(label, train_passage_dict)
        folds = get_folds(train_instances, args.num_folds)
        for mult in OA_Multiples:
            print "OpenAccess weight: ", wt*mult
            update_passagedict(label, OA_passage_dict, wt*mult)
            OAtrain_instances = create_instances(label, OA_passage_dict)
            weight_str = getcomb_wgtstr(wt, wt*mult)
            ensure_dir(os.path.join(work_dir, folder_prefix + weight_str))
            fold_dir = os.path.join(work_dir, folder_prefix + weight_str)

            create_foldfiles(label, folds, len(train_instances) + len(OAtrain_instances), OAtrain_instances, fold_dir, weight_str)
            ensure_dir(os.path.join(fold_dir, args.optim_method))
            models_dir = os.path.join(fold_dir, args.optim_method)

            l2_str = "_L2-" + getstr(0.1); l2_cmd = " --l2 " + str(0.1)
            #l2_str = ""; l2_cmd = ""
            combined_predictions_file = os.path.join(models_dir, label + l2_str + "_AllPredictions.txt")
            combined_confidences_file = os.path.join(models_dir, label + l2_str + "_AllConfidences.txt")
            PRpoints_file = os.path.join(models_dir, label + l2_str + "_PRPointsTabDelim.txt")
            all_predictions = []

            for i in xrange(int(args.num_folds)):
                testfoldfile = getFile("Test", fold_dir, i)
                trainfoldfile = getFile("Train", fold_dir, i)
                model_file = os.path.join(models_dir, label + "_Fold-" + str(i) + l2_str + "_VWModel.vw")
                predictions_file = os.path.join(models_dir, label + "_Fold-" + str(i) + l2_str + "_Predictions.txt")

                cmd = "vw -d " + trainfoldfile + " -c --passes 100 -f " + model_file + " --loss_function logistic" + method + l2_cmd
                run_command(cmd)

                cmd = "vw -d " + testfoldfile + " -t -i " + model_file + " -p /dev/stdout | /ua/ml-group/big-mechanism-project/vowpal_wabbit/utl/logistic -0 > " + predictions_file
                run_command(cmd)

                readfile(all_predictions, predictions_file)

            writefile(combined_predictions_file, all_predictions)
            cmd = "python Generate_Confidences.py " + label + " " + combined_predictions_file + " " + combined_confidences_file + " ./pubmedid_extras.p"
            run_command(cmd)

            (actual_pos, actual_neg) = get_actual_pos_neg(combined_predictions_file, DatumKB_dict, label)

            cmd = "java -cp $PR_CURVE_JAR LRTest_CLArg_TabDelim " + combined_confidences_file + " " + PRpoints_file + " " + actual_pos
            run_command(cmd)
            clean_files("./", "Sorted_*txt")

            cmd = "java -jar $JAVA_JARS_HOME/auc.jar " + PRpoints_file + " PR " + actual_pos + " " + actual_neg
            auc = get_AUPR(cmd)
            print "The area under the PR Curve is: ", auc
            if (auc - best_weights_auc[0]) > 0.01:
                best_weights_auc = (auc, float(wt), float(wt*mult))
    print "************************************************************************************************************************************"
    print "The best AUC is: ", best_weights_auc[0]
    print "The best Non-OpenAccess weights are: ", best_weights_auc[1], "\t The best OpenAccess weights are: ", best_weights_auc[2]
    print "************************************************************************************************************************************"

    print "Creating optimally weighted Training file"
    optimal_NoOA_wt = best_weights_auc[1]
    optimal_OA_wt = best_weights_auc[2]
    update_passagedict(label, train_passage_dict, optimal_NoOA_wt)
    train_instances = create_instances(label, train_passage_dict)
    update_passagedict(label, OA_passage_dict, optimal_OA_wt)
    OAtrain_instances = create_instances(label, OA_passage_dict)
    weight_str = getcomb_wgtstr(int(optimal_NoOA_wt), int(optimal_OA_wt))

    train_records = train_instances + OAtrain_instances
    shuffle(train_records)
    if args.feedback is None:
        trainfilename = "vw_" + label + "_Weights_" + weight_str + "_Training_File.txt"
    else:
        trainfilename = "vw_" + label + "_Weights_" + weight_str + "_WithFeedback_Training_File.txt"
    writefile(os.path.join(work_dir, trainfilename), train_records)
    print "Training File created at: ", work_dir


def getcomb_wgtstr(weight1, weight2):
    mystr1 = str(weight1)
    mystr2 = str(weight2)
    return mystr1 + "_" + mystr2


def get_folds(train_records, num_folds):
    posrecords = []; negrecords = []
    for train_instance in train_records:
        train_instance = train_instance.strip()
        currlabel = train_instance.split()[0].strip()
        if currlabel == "1":
            posrecords.append(train_instance)
        elif currlabel == "-1":
            negrecords.append(train_instance)
    shuffle(posrecords); shuffle(negrecords)
    total_recs = len(posrecords) + len(negrecords)
    folds = []
    create_folds(posrecords, folds, int(num_folds))
    create_folds(negrecords, folds, int(num_folds))
    check_foldsize(folds, total_recs)
    return folds


def create_workDir(passagedict_dir, label, optim_method):
    if os.path.exists(os.path.join(passagedict_dir, label)):
        if os.path.exists(os.path.join(passagedict_dir, label, "Weights_" + optim_method)):
            shutil.rmtree(os.path.join(passagedict_dir, label, "Weights_" + optim_method))
        os.mkdir(os.path.join(passagedict_dir, label, "Weights_" + optim_method))
    else:
        os.makedirs(os.path.join(passagedict_dir, label, "Weights_" + optim_method))
    return os.path.join(passagedict_dir, label, "Weights_" + optim_method)


def ensure_dir(file_path):
    if os.path.exists(file_path):
        shutil.rmtree(file_path)
    os.mkdir(file_path)


def create_foldfiles(label, folds, total_recs, OAtrain_instances, fold_dir, weight_str):
    #print "Creating train and test files for each cross validation fold"
    for i, elem in enumerate(folds):
        shuffle(elem)
        testfoldfile = "vw_" + label + "_Weights_" + weight_str + "_TestFold-" + str(i) + ".txt"
        writefile(os.path.join(fold_dir, testfoldfile), elem)
        trainfold_recs = []
        for j, myelem in enumerate(folds):
            if j == i:
                continue
            trainfold_recs.extend(myelem)
        trainfold_recs.extend(OAtrain_instances)
        shuffle(trainfold_recs)
        if len(elem) + len(trainfold_recs) != total_recs:
            print "Something went wrong while creating the train and test files for each fold"
            exit()
        trainfoldfile = "vw_" + label + "_Weights_" + weight_str + "_TrainFold-" + str(i) + ".txt"
        writefile(os.path.join(fold_dir, trainfoldfile), trainfold_recs)


def update_passagedict(label, passage_dict, total_wt):
    temp_dict = {}
    for pmid in passage_dict:
        if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
            continue
        temp_dict[pmid] = {}
        for myclass in passage_dict[pmid][label]:
            for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
                if prot not in temp_dict[pmid]:  # The same protein can have pos as well as neg instances after user gives feedback
                    temp_dict[pmid][prot] = {}
                    temp_dict[pmid][prot]["numInstances"] = 0

                temp_dict[pmid][prot]["numInstances"] += len(passage_dict[pmid][label][myclass][prot]["passageDetails"])
                if myclass == "Pos":
                    temp_dict[pmid][prot]["Pos"] = len(passage_dict[pmid][label][myclass][prot]["passageDetails"])
                else:
                    temp_dict[pmid][prot]["Neg"] = len(passage_dict[pmid][label][myclass][prot]["passageDetails"])

    for pmid in temp_dict:
        for prot in temp_dict[pmid]:
            temp_dict[pmid][prot]["weight"] = format(total_wt / float(temp_dict[pmid][prot]["numInstances"]), '.2f')

    for pmid in passage_dict:
        if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
            continue
        for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label][myclass] is a dict
                for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                    contg_psg_dict["weight"] = temp_dict[pmid][prot]["weight"]


def create_instances(label, passage_dict):
    myrecords = []
    for pmid in passage_dict:
        if label not in passage_dict[pmid] or not passage_dict[pmid][label]:
            continue
        for myclass in passage_dict[pmid][label]:  # passage_dict[pmid][label] is a dict containing only 2 keys - "Pos" and "Neg"
            for prot in passage_dict[pmid][label][myclass]:  # passage_dict[pmid][label][myclass] is a dict
                for contg_psg_dict in passage_dict[pmid][label][myclass][prot]["passageDetails"]:  # passage_dict[pmid][label][myclass][prot] is a list of dicts.
                    vw_label = contg_psg_dict["class"]  # contg_psg_dict represents one training instance
                    vw_tag = contg_psg_dict["tag"]
                    lex_feat_str = create_lexical_features(contg_psg_dict["textOfInterest"])
                    if "weight" in contg_psg_dict:
                        rec = vw_label + " " + contg_psg_dict["weight"] + " " + vw_tag + "|LexicalFeatures " + lex_feat_str  # + " |ProteinFeature " + prot
                    else:
                        rec = vw_label + " " + vw_tag + "|LexicalFeatures " + lex_feat_str  # + " |ProteinFeature " + prot
                    myrecords.append(rec)

    return myrecords


if __name__ == "__main__":
    main_body()