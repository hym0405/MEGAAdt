#!/usr/bin/env python
import os
import sys
def main():
    In = sys.argv[1]
    Out = sys.argv[2]
    fileNum = int(sys.argv[3])
    f = open(In + ".0.tsv", "rb")
    data = f.readlines()
    f.close()
    fw = open(Out, "w")
    fw.writelines([data[0]])
    for i in range(fileNum):
        f = open(In + "." + str(i) + ".tsv", "rb")
        data = f.readlines()
        f.close()
        for each in data[1:]:
            fw.writelines([each])
    fw.close()
if __name__ == "__main__":
    main()
