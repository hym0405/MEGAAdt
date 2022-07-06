#!/usr/bin/env python
import os
import sys
def main():
    In = sys.argv[1]
    Out = sys.argv[2]
    f = open(In, "rb")
    data = f.readlines()
    f.close()
    f = open(Out, "w")
    for i in range(len(data) / 4):
        label = data[4 * i][1:-1]
        seq = data[4 * i + 1][:-1]
        f.writelines([">" + label + os.linesep])
        f.writelines([seq + os.linesep])
    f.close()
if __name__ == "__main__":
    main()
