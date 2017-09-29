import argparse;
import sys;
reload(sys)
sys.setdefaultencoding("utf-8")
from random import shuffle


# python Concatenate_TrainingFiles.py <Training_File> <OpenAccess_Training_File> <Output_Training_File>


def main_body():
    parser = argparse.ArgumentParser(prog='Concatenate_TrainingFiles', usage='Concatenate_TrainingFiles.py <NoOA_Training_File> <OA_Training_File> <Output_Training_File>', description='Script to concatenate training files')
    parser.add_argument('trainFile', help='Training File which does not have any OpenAccess articles')
    parser.add_argument('OAtrainFile', help='Training File containing 61 OpenAccess articles')
    parser.add_argument('outfilename', help='Output Training File containing all 61 OpenAccess along with the other articles')

    args = parser.parse_args()

    myrecords = []
    with open(args.trainFile, 'rt') as f:
        for train_instance in f:
            train_instance = train_instance.strip()
            myrecords.append(train_instance)

    with open(args.OAtrainFile, 'rt') as f:
        for OAtrain_instance in f:
            OAtrain_instance = OAtrain_instance.strip()
            myrecords.append(OAtrain_instance)

    shuffle(myrecords)
    data_fn1 = open(args.outfilename, 'w')
    for rc in myrecords:
        data_fn1.write(rc + '\n')
    data_fn1.close()


if __name__ == "__main__":
    main_body()