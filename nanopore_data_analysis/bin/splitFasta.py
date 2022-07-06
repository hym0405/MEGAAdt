#!/usr/bin/env python
import os
import sys
def main():
    In = sys.argv[1]
    Out = sys.argv[2]
    fileNum = int(sys.argv[3])
    f = open(In, "rb")
    data = f.readlines()
    f.close()
    lineNum = len(data) / 2

    lineCount = lineNum / fileNum

    for i in range(fileNum - 1):
        f = open(Out + "." + str(i) + ".fasta", "w")
        for j in range(i * lineCount, (i + 1) * lineCount):
            f.writelines([data[2 * j]])
            f.writelines([data[2 * j + 1]])
        f.close()
    f = open(Out + "." + str(fileNum - 1) + ".fasta", "w")
    for j in range((fileNum - 1) * lineCount, lineNum):
        f.writelines([data[2 * j]])
        f.writelines([data[2 * j + 1]])
    f.close()
if __name__ == "__main__":
    main()

