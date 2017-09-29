import argparse;
import sys;
reload(sys)
import os
sys.setdefaultencoding("utf-8")


# python Insert_InstanceWeights.py <Training_File> <Positive Real num> <outputdir>


def main_body():
    parser = argparse.ArgumentParser(prog='Insert_InstanceWeights', usage='Insert_InstanceWeights.py <Training_File> <Instance Weight> <csv_loc>', description='Script to insert importance weights for every training instance')
    parser.add_argument('trainFile', help='Training File which does not have any importance weights')
    parser.add_argument('weight', help='Importance weight which must be a real number')
    parser.add_argument('csvLoc', help='The directory in which to store the weighted training file')

    args = parser.parse_args()
    weight = args.weight

    if not is_number(weight):
        print "Please provide a positive real number. No strings allowed!"
        exit()

    if float(weight) < 0.0:
        print "The importance weight should be a positive real number"
        exit()

    weight = float(weight)
    filetag = getfiletag(str(weight))

    print weight, "\t", filetag

    trainFile_full_path = os.path.join(args.csvLoc, args.trainFile)

    myrecords = []
    with open(trainFile_full_path, 'rt') as f:
        for train_instance in f:
            train_instance = train_instance.strip()
            parts = train_instance.split(" ", 1)
            parts.insert(1, str(weight))
            new_train_instance = " ".join(parts)
            #print new_train_instance; exit()
            myrecords.append(new_train_instance)

    outfilename = args.trainFile.split(".")[0] + "_" + filetag + "." + args.trainFile.split(".")[1]
    print outfilename
    outfile_full_path = os.path.join(args.csvLoc, outfilename)
    data_fn1 = open(outfile_full_path, 'w')
    for rc in myrecords:
        data_fn1.write(rc + '\n')
    data_fn1.close()


def getfiletag(ft):
    if "." in ft:
        ft = ft.replace(".", "_")
    return ft


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


if __name__ == "__main__":
    main_body()