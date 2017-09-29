# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse
from pprint import pprint
import sys
reload(sys)
sys.setdefaultencoding("utf-8")


# python Parse_ParamTuningLogs.py ./RNAi_Baseline_Weights_1_1_L1L2Together.txt


def main_body():
    parser = argparse.ArgumentParser(prog='Parse_ParamTuningLogs', usage='Parse_ParamTuningLogs.py <logfile>', description='Script to parse log file and return the best AUC value')
    parser.add_argument('logfile', help='Log file to parse')

    args = parser.parse_args()
    mydict = {}
    with open(args.logfile, 'rt') as f1:
        for line in f1:
            line = line.strip()
            if not line: continue
            if line.startswith("L1 Regularization lambda value"):
                L1_lambda = line.split(":")[1].strip()
                if L1_lambda not in mydict:
                    mydict[L1_lambda] = {}
                else:
                    print "Something went wrong! L1 lamda value is already present in dict: ", L1_lambda
                    exit()
            if line.startswith("L2 Regularization lambda value"):
                L2_lambda = line.split(":")[1].strip()
            if line.startswith("The area under the PR Curve is"):
                AUC_val = float(line.split(":")[1].strip())
                if L2_lambda not in mydict[L1_lambda]:
                    mydict[L1_lambda][L2_lambda] = AUC_val
                else:
                    print "Something went wrong! L2 lambda value: ", L2_lambda, " is already present within the dict for L1 lambda value: ", L1_lambda
                    exit()
    mylist = []
    for L1_val in mydict:
        for L2_val in mydict[L1_val]:
            mylist.append((mydict[L1_val][L2_val], L1_val, L2_val))

    mylist.sort(key=lambda x: x[0], reverse=True)
    pprint(mylist)
    print "Results for log file: ", args.logfile
    print "The largest AUC: ", mylist[0][0], "\nL1 value: ", mylist[0][1], "\t L2 value: ", mylist[0][2]


if __name__ == "__main__":
    main_body()