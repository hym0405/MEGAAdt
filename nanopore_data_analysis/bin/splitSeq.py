#!/usr/bin/env python
import os
import sys
def main():
    In = sys.argv[1]
    ref = sys.argv[2]
    Out = sys.argv[3]
    f = open(ref, "rb")
    data = f.readlines()
    f.close()
    Pool = {}
    for i in range(len(data) / 4):
        label = data[4 * i][1:-1]
        seq = data[4 * i + 1][:-1]
        dum = data[4 * i + 2][:-1]
        quality = data[4 * i + 3][:-1]
        Pool[label] = [seq, dum, quality]

    f = open(In, "rb")
    data = f.readlines()
    f.close()
    for each in data:
        each = each[:-1]
        tmp = each.split("\t")
        label = tmp[0]
        sample = tmp[5]
        if tmp[5] != "NA":
            f = open(Out + "/" + tmp[5] + ".fastq", "a")
            f.writelines([">" + label + os.linesep])
            f.writelines([Pool[label][0] + os.linesep])
            f.writelines([Pool[label][1] + os.linesep])
            f.writelines([Pool[label][2] + os.linesep])
            f.close()
if __name__ == "__main__":
    main()
