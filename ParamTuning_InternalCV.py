# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import os
import glob
import copy
import shutil
import subprocess
import cPickle as pickle;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")


# python ParamTuning_InternalCV.py RNAi ./datafiles_TriggerHighlights_BaseLine2/RNAi/folds_1_1 10 bfgs


def main_body():
    parser = argparse.ArgumentParser(prog='ParamTuning_InternalCV', usage='ParamTuning_InternalCV.py <label> <foldDir> <num_folds> <opt_method>', description='Script to tune parameters and return the best values for the params')
    parser.add_argument('label', help='Predicate being analyzed')
    parser.add_argument('foldDir', help='Dir where the datafiles for all folds are located')
    parser.add_argument('num_folds', help='Number of folds')
    parser.add_argument('optim_method', help='The optimization method for VW')

    args = parser.parse_args()

    allowed_labels = ['reqs', 'dnreqs', 'RNAi', 'KO', 'omission', 'DKO', 'DRNAi']
    if args.label not in allowed_labels:
        print "The given label is not allowed: ", args.label, ". Label may only be any one of \'reqs\', \'dnreqs\', \'RNAi\', \'KO\', \'omission\', \'DKO\' \'DRNAi\'"
        print "Please try again!"
        exit()

    fold_dir = args.foldDir
    label = args.label

    DatumKB_dict = pickle.load(open("pubmedid_extras.p", "rb"))

    if args.optim_method == "bfgs":
        method = " --bfgs --mem 20" # BFGS with holdout_off is a joke! Don't ever try the --holdout_off switch with BFGS
    elif args.optim_method == "ftrl":
        method = " --ftrl --holdout_off"
    else:
        print "Cannot interpret the provided optimization method: ", args.optim_method
        print "Values allowed currently are \'bfgs\' and \'ftrl\'!"
        exit()

    clean_files(fold_dir, "*.cache")
    ensure_dir(os.path.join(fold_dir, args.optim_method))
    models_dir = os.path.join(fold_dir, args.optim_method)

    L1_lambda_vals = [0, 0.1, 0.25, 0.5, 0.75, 0.01, 0.05, 0.001, 0.005, 2, 5, 10, 12, 15, 20, 50, 75, 100, 150, 250, 500]
    L2_lambda_vals = [0, 0.1, 0.25, 0.5, 0.75, 0.01, 0.05, 0.001, 0.005, 2, 5, 10, 12, 15, 20, 50, 75, 100, 150, 250, 500]
    ftrl_alpha_vals = [0, 0.1, 0.25, 0.5, 0.75, 0.01, 0.05, 0.001, 0.005, 2, 5, 10, 12, 15, 20, 50, 75, 100, 150, 250, 500]
    ftrl_beta_vals = [0, 0.1, 0.25, 0.5, 0.75, 0.01, 0.05, 0.001, 0.005, 2, 5, 10, 12, 15, 20, 50, 75, 100, 150, 250, 500]

    print "L2 Regularization"; sys.stdout.flush()
    best_l2_auc = (0.0, 0.0)
    for l2_val in L2_lambda_vals:
        print "L2 Regularization lambda value: ", l2_val; sys.stdout.flush()
        l2_str = "_L2-" + getstr(l2_val)
        l2_cmd = " --l2 " + str(l2_val)
        combined_predictions_file = os.path.join(models_dir, label + "_AllPredictions" + l2_str + ".txt")
        combined_confidences_file = os.path.join(models_dir, label + "_AllConfidences" + l2_str + ".txt")
        PRpoints_file = os.path.join(models_dir, label + "_PRPointsTabDelim" + l2_str + ".txt")
        all_predictions = []

        for i in xrange(int(args.num_folds)):
            testfoldfile = getFile("Test", fold_dir, i)
            trainfoldfile = getFile("Train", fold_dir, i)
            model_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_VWModel" + l2_str + ".vw")
            predictions_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_Predictions" + l2_str + ".txt")

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
        print "The area under the PR Curve is: ", auc; sys.stdout.flush()
        #if l2_val == 10:    exit()
        if (auc - best_l2_auc[0]) > 0.01:
            best_l2_auc = (auc, l2_val)
    print "The best AUC is: ", best_l2_auc[0]; sys.stdout.flush()
    print "The best L2 lambda value is: ", best_l2_auc[1]; sys.stdout.flush()

    if args.optim_method == "ftrl":     # BFGS is insensitive to L1 Regularization values!

        print "L1 Regularization"
        sys.stdout.flush()
        best_l1_auc = (0.0, 0.0)
        for l1_val in L1_lambda_vals:
            print "L1 Regularization lambda value: ", l1_val
            sys.stdout.flush()
            l1_str = "_L1-" + getstr(l1_val);
            l1_cmd = " --l1 " + str(l1_val)
            combined_predictions_file = os.path.join(models_dir, label + "_AllPredictions" + l1_str + ".txt")
            combined_confidences_file = os.path.join(models_dir, label + "_AllConfidences" + l1_str + ".txt")
            PRpoints_file = os.path.join(models_dir, label + "_PRPointsTabDelim" + l1_str + ".txt")
            all_predictions = []

            for i in xrange(int(args.num_folds)):
                testfoldfile = getFile("Test", fold_dir, i)
                trainfoldfile = getFile("Train", fold_dir, i)
                model_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_VWModel" + l1_str + ".vw")
                predictions_file = os.path.join(models_dir,
                                                label + "_Fold-" + str(i) + "_Predictions" + l1_str + ".txt")

                cmd = "vw -d " + trainfoldfile + " -c --passes 100 -f " + model_file + " --loss_function logistic" + method + l1_cmd
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
            sys.stdout.flush()
            if (auc - best_l1_auc[0]) > 0.01:
                best_l1_auc = (auc, l1_val)
        print "The best AUC is: ", best_l1_auc[0]
        sys.stdout.flush()
        print "The best L1 lambda value is: ", best_l1_auc[1]
        sys.stdout.flush()

        print "L1 and L2 Regularization - Together"; sys.stdout.flush()
        best_l1_l2_auc = (0.0, 0.0, 0.0)
        for l1_val in L1_lambda_vals:
            print "L1 Regularization lambda value: ", l1_val; sys.stdout.flush()
            l1_str = "_L1-" + getstr(l1_val)
            l1_cmd = " --l1 " + str(l1_val)
            for l2_val in L2_lambda_vals:
                print "L2 Regularization lambda value: ", l2_val; sys.stdout.flush()
                l2_str = "_L2-" + getstr(l2_val)
                l2_cmd = " --l2 " + str(l2_val)
                combined_predictions_file = os.path.join(models_dir, label + "_AllPredictions" + l1_str + l2_str + ".txt")
                combined_confidences_file = os.path.join(models_dir, label + "_AllConfidences" + l1_str + l2_str + ".txt")
                PRpoints_file = os.path.join(models_dir, label + "_PRPointsTabDelim" + l1_str + l2_str + ".txt")
                all_predictions = []

                for i in xrange(int(args.num_folds)):
                    testfoldfile = getFile("Test", fold_dir, i)
                    trainfoldfile = getFile("Train", fold_dir, i)
                    model_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_VWModel" + l1_str + l2_str + ".vw")
                    predictions_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_Predictions" + l1_str + l2_str + ".txt")

                    cmd = "vw -d " + trainfoldfile + " -c --passes 100 -f " + model_file + " --loss_function logistic" + method + l1_cmd + l2_cmd
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
                print "The area under the PR Curve is: ", auc; sys.stdout.flush()
                if (auc - best_l1_l2_auc[0]) > 0.01:
                    best_l1_l2_auc = (auc, l1_val, l2_val)
        print "The best value for AUC is: ", best_l1_l2_auc[0]; sys.stdout.flush()
        print "The best L1 and L2 lambda values together are: L1 - ", best_l1_l2_auc[1], ", L2 - ", best_l1_l2_auc[2]; sys.stdout.flush()


        print "FTRL Alpha and Beta Regularization - Together"
        best_alpha_beta_auc = (0.0, 0.0, 0.0)
        for al_val in ftrl_alpha_vals:
            print "FTRL Alpha Regularization lambda value: ", al_val
            alpha_str = "_FTRL_alpha-" + getstr(al_val)
            alpha_cmd = " --ftrl_alpha " + str(al_val)
            for be_val in ftrl_beta_vals:
                #print "FTRL Beta Regularization lambda value: ", be_val
                beta_str = "_FTRL_beta-" + getstr(be_val)
                beta_cmd = " --ftrl_beta " + str(be_val)
                combined_predictions_file = os.path.join(models_dir, label + "_AllPredictions" + alpha_str + beta_str + ".txt")
                combined_confidences_file = os.path.join(models_dir, label + "_AllConfidences" + alpha_str + beta_str + ".txt")
                PRpoints_file = os.path.join(models_dir, label + "_PRPointsTabDelim" + alpha_str + beta_str + ".txt")
                all_predictions = []

                for i in xrange(int(args.num_folds)):
                    testfoldfile = getFile("Test", fold_dir, i)
                    trainfoldfile = getFile("Train", fold_dir, i)
                    model_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_VWModel" + alpha_str + beta_str + ".vw")
                    predictions_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_Predictions" + alpha_str + beta_str + ".txt")

                    cmd = "vw -d " + trainfoldfile + " -c --passes 100 -f " + model_file + " --loss_function logistic" + method + alpha_cmd + beta_cmd
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
                # print "The area under the PR Curve is: ", auc
                if (auc - best_alpha_beta_auc[0]) > 0.01:
                    best_alpha_beta_auc = (auc, al_val, be_val)
        print "The best value for AUC is: ", best_alpha_beta_auc[0]
        print "The best FTRL alpha and beta lambda values together are: Alpha - ", best_alpha_beta_auc[1], ", Beta - ", best_alpha_beta_auc[2]

        # print "L1 Reg, L2 Reg, FTRL Alpha and Beta Regularization - All 4 Together"
        # best_l1_l2_alpha_beta_auc = (0.0, 0.0, 0.0, 0.0, 0.0)
        # for l1_val in L1_lambda_vals:
        #     print "L1 Regularization lambda value: ", l1_val
        #     l1_str = "_L1-" + getstr(l1_val)
        #     l1_cmd = " --l1 " + str(l1_val)
        #     for l2_val in L2_lambda_vals:
        #         #print "L2 Regularization lambda value: ", l2_val
        #         l2_str = "_L2-" + getstr(l2_val)
        #         l2_cmd = " --l2 " + str(l2_val)
        #         for al_val in ftrl_alpha_vals:
        #             #print "FTRL Alpha Regularization lambda value: ", al_val
        #             alpha_str = "_FTRL_alpha-" + getstr(al_val)
        #             alpha_cmd = " --ftrl_alpha " + str(al_val)
        #             for be_val in ftrl_beta_vals:
        #                 #print "FTRL Beta Regularization lambda value: ", be_val
        #                 beta_str = "_FTRL_beta-" + getstr(be_val)
        #                 beta_cmd = " --ftrl_beta " + str(be_val)
        #                 combined_predictions_file = os.path.join(models_dir, label + "_AllPredictions" + l1_str + l2_str + alpha_str + beta_str + ".txt")
        #                 combined_confidences_file = os.path.join(models_dir, label + "_AllConfidences" + l1_str + l2_str + alpha_str + beta_str + ".txt")
        #                 PRpoints_file = os.path.join(models_dir, label + "_PRPointsTabDelim" + l1_str + l2_str + alpha_str + beta_str + ".txt")
        #                 all_predictions = []
        #
        #                 for i in xrange(int(args.num_folds)):
        #                     testfoldfile = getFile("Test", fold_dir, i)
        #                     trainfoldfile = getFile("Train", fold_dir, i)
        #                     model_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_VWModel" + l1_str + l2_str + alpha_str + beta_str + ".vw")
        #                     predictions_file = os.path.join(models_dir, label + "_Fold-" + str(i) + "_Predictions" + l1_str + l2_str + alpha_str + beta_str + ".txt")
        #
        #                     cmd = "vw -d " + trainfoldfile + " -c --passes 100 -f " + model_file + " --loss_function logistic" + method + l1_cmd + l2_cmd + alpha_cmd + beta_cmd
        #                     run_command(cmd)
        #
        #                     cmd = "vw -d " + testfoldfile + " -t -i " + model_file + " -p /dev/stdout | /ua/ml-group/big-mechanism-project/vowpal_wabbit/utl/logistic -0 > " + predictions_file
        #                     run_command(cmd)
        #
        #                     readfile(all_predictions, predictions_file)
        #
        #                 writefile(combined_predictions_file, all_predictions)
        #                 cmd = "python Generate_Confidences.py " + label + " " + combined_predictions_file + " " + combined_confidences_file + " ./pubmedid_extras.p"
        #                 run_command(cmd)
        #
        #                 (actual_pos, actual_neg) = get_actual_pos_neg(combined_predictions_file, DatumKB_dict, label)
        #
        #                 cmd = "java -cp $PR_CURVE_JAR LRTest_CLArg_TabDelim " + combined_confidences_file + " " + PRpoints_file + " " + actual_pos
        #                 run_command(cmd)
        #                 clean_files("./", "Sorted_*txt")
        #
        #                 cmd = "java -jar $JAVA_JARS_HOME/auc.jar " + PRpoints_file + " PR " + actual_pos + " " + actual_neg
        #                 auc = get_AUPR(cmd)
        #                 # print "The area under the PR Curve is: ", auc
        #                 if (auc - best_l1_l2_alpha_beta_auc[0]) > 0.01:
        #                     best_l1_l2_alpha_beta_auc = (auc, l1_val, l2_val, al_val, be_val)
        # print "The best value for AUC is: ", best_l1_l2_alpha_beta_auc[0]
        # print "The below 4 parameters were all tested together!"
        # print "The L1 and L2 lambda values are: L1 - ", best_l1_l2_alpha_beta_auc[1], ", L2 - ", best_l1_l2_alpha_beta_auc[2]
        # print "The FTRL alpha and beta lambda values are: Alpha - ", best_l1_l2_alpha_beta_auc[3], ", Beta - ", best_l1_l2_alpha_beta_auc[4]

        all_areas = [best_l1_auc[0], best_l2_auc[0], best_l1_l2_auc[0], best_alpha_beta_auc[0]] #, best_l1_l2_alpha_beta_auc[0]]
        max_index = all_areas.index(max(all_areas))
        if max_index == 0:
            print "L1 Regularization wins!"
            print best_l1_auc
        elif max_index == 1:
            print "L2 Regularization wins!"
            print best_l2_auc
        elif max_index == 2:
            print "A combination of L1 and L2 Regularization wins!"
            print best_l1_l2_auc
            print "L1 Lambda value is: ", best_l1_l2_auc[1], "\t L2 Lambda value is: ", best_l1_l2_auc[2]
        elif max_index == 3:
            print "A combination of FTRL Alpha and FTRL Beta Regularization wins!"
            print best_alpha_beta_auc
            print "FTRL Alpha value is: ", best_alpha_beta_auc[1], "\t FTRL Beta value is: ", best_alpha_beta_auc[2]
        # elif max_index == 4:
        #     print "A combination of L1, L2, FTRL Alpha and FTRL Beta Regularization wins!"
        #     print best_l1_l2_alpha_beta_auc
        #     print "L1 Lambda value is: ", best_l1_l2_alpha_beta_auc[1], "\t L2 Lambda value is: ", best_l1_l2_alpha_beta_auc[2], "\t FTRL Alpha value is: ", best_l1_l2_alpha_beta_auc[3], "\t FTRL Beta value is: ", best_l1_l2_alpha_beta_auc[4]

    # elif args.optim_method == "bfgs":
    #     all_areas = [best_l1_auc[0], best_l2_auc[0], best_l1_l2_auc[0]]
    #     max_index = all_areas.index(max(all_areas))
    #     if max_index == 0:
    #         print "L1 Regularization wins!"
    #         print best_l1_auc
    #     elif max_index == 1:
    #         print "L2 Regularization wins!"
    #         print best_l2_auc
    #     elif max_index == 2:
    #         print "A combination of L1 and L2 Regularization wins!"
    #         print best_l1_l2_auc
    #         print "L1 Lambda value is: ", best_l1_l2_auc[1], "\t L2 Lambda value is: ", best_l1_l2_auc[2]


def getstr(myval):
    mystr = str(myval)
    if "." in mystr:
        mystr = mystr.replace(".", "")
    return mystr


def getFile(fileType, fold_dir, i):
    search_str = "*_" + fileType + "Fold-" + str(i) + "*txt"
    myfilelist = glob.glob(os.path.join(fold_dir, search_str))
    if len(myfilelist) > 1:
        print "Multiple files were matched by glob! \n", myfilelist
        exit()
    return myfilelist[0]


def run_command(cmd):
    #print "-----------------------------------------------------------------------------------------------------------"
    #print cmd
    #print "-----------------------------------------------------------------------------------------------------------"
    my_env = os.environ.copy()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

    for line in p.stderr.readlines():  # Do NOT comment this for statement
        #print line.strip()
        pass


def get_AUPR(cmd):
    #print "-----------------------------------------------------------------------------------------------------------"
    #print cmd
    #print "-----------------------------------------------------------------------------------------------------------"
    my_env = os.environ.copy()
    #print "Saswati: ", my_env["JAVA_JARS_HOME"]
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=my_env)

    for line in p.stdout.readlines():  # Do NOT comment this for statement
        # print line.strip()
        if line.strip().startswith("Area Under the Curve for Precision - Recall"):
            auc = line.strip().split()[-1]
    return float(auc)


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


def clean_files(fold_dir, pattern):
    for cache_file in glob.glob(os.path.join(fold_dir, pattern)):
        #print "Removing file: ", cache_file
        os.remove(cache_file)


def get_actual_pos_neg(pfile, DatumKB_dict, label):
    pred_dict = {}     # This is basically all the PMIDs in the training file
    generate_prediction_confs(pfile, pred_dict)

    temp_dict = {}  # temp_dict is a flattened version of DatumKB_dict. It only contains those PMIDs that have extras datums of interest and are included in the Test file list.
    create_act_dict(pred_dict.keys(), temp_dict, DatumKB_dict)
    pos = 0
    for pmid in temp_dict.keys():
        for prot in temp_dict[pmid]:    # For a particular PMID, there may be some proteins in temp_dict which are not present in pred_dict
            if label.strip() in temp_dict[pmid][prot]:
                pos += 1
    # Use this value as the denominator for Recall, not the positives present in the test file
    # print "Actual positives in the training file for label ", label.strip(), " with respect to the Datum KB: ", pos

    neg = 0; total_instances = 0; pos_trainfile = 0;
    for pmid in pred_dict.keys():
        if pmid not in temp_dict:  # Article does not have any extras datum that are of interest (i.e. it does not have any 'reqs' and 'do not reqs' extras datums)
            for prot in pred_dict[pmid]:
                neg += 1; total_instances += 1
            continue
        for prot in pred_dict[pmid]:  # Article does has extras datum that are of interest to us (i.e. it has 'reqs' and 'do not reqs' extras datums)
            total_instances += 1
            if prot in temp_dict[pmid]:
                if label.strip() not in temp_dict[pmid][prot]:
                    neg += 1
                else:
                    pos_trainfile += 1
            else:  # Not all proteins detected in the article are entities of 'reqs' or 'dnreqs' extras datums
                neg += 1
    # print "Total no. of instances in the training file: ", total_instances
    # print "Actual negatives in the training file for label ", label.strip(), " is: ", neg
    # print "Actual positives in the training file for label ", label.strip(), " with respect to the training file: ", pos_trainfile
    return str(pos), str(neg)


def create_act_dict(mypmids, tdict, mydict):
    #print "Total no of pmids in the training file: ", len(mypmids)
    for pmid in mypmids:
        if pmid in mydict:  # If the PMID is present in the DatumKB
            tdict[pmid] = {}          # tdict contains only those PMIDs that exist both in DatumKB and the train file
            for ent in mydict[pmid]:
                mylist = []
                if len(mydict[pmid][ent]['reqstest']) > 0:
                    mylist.append("reqs"); mylist.extend(mydict[pmid][ent]["reqstest"])
                if len(mydict[pmid][ent]['dnreqstest']) > 0:
                    mylist.append("dnreqs"); mylist.extend(mydict[pmid][ent]["dnreqstest"])
                tdict[pmid][ent] = list(set(mylist))


def generate_prediction_confs(predict_file, pred_dict):
    temp_dict = {}
    create_pred_dict(predict_file, temp_dict)
    for pmid in temp_dict.keys():
        pred_dict[pmid] = {}
        for prot in temp_dict[pmid]:
            if "__" in prot:
                myprot = prot.replace("__", " ")
            else:
                myprot = prot
            pred_dict[pmid][myprot] = copy.deepcopy(temp_dict[pmid][prot])


def create_pred_dict(pfile, tdict):
    with open(pfile, 'rt') as f1:
        for line in f1:
            line_list = line.split()
            pmid = line_list[1].strip().split("_", 1)[0]
            prots = line_list[1].strip().split("--")[1:]
            prots = list(set(prots))    # prots is a distinct list of proteins extracted from the tag of an instance
            if pmid.strip() in tdict:
                for prt in prots:   # A tag will have multiple proteins only if they are positive
                    if prt.strip() in tdict[pmid.strip()]:
                        tdict[pmid.strip()][prt.strip()].append(float(line_list[0].strip()))
                    else:
                        tdict[pmid.strip()][prt.strip()] = [float(line_list[0].strip())]
            else:
                tdict[pmid.strip()] = {}
                for prt in prots:
                    tdict[pmid.strip()][prt.strip()] = [float(line_list[0].strip())]


if __name__ == "__main__":
    main_body()