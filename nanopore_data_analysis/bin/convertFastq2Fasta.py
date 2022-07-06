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
    count = 1
    for i in range(len(data) / 4):
        label = ">seq" + str(count)
        seq = data[4 * i + 1][:-2]
        f.writelines([label + "_" + str(len(seq)) + os.linesep])
        f.writelines([seq + os.linesep])
        count += 1
    f.close()
if __name__ == "__main__":
    main()
