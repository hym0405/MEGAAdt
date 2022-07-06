#!/usr/bin/env python
import os
import sys
def main():
    In = sys.argv[1]
    Out = sys.argv[2]
    cutoff = int(sys.argv[3])
    f = open(In, "rb")
    data = f.readlines()
    f.close()
    f = open(Out, "w")
    for i in range(len(data) / 4):
        if len(data[4 * i + 1][:-1]) >= cutoff:
            f.writelines([data[4 * i]])
            f.writelines([data[4 * i + 1]])
            f.writelines([data[4 * i + 2]])
            f.writelines([data[4 * i + 3]])
    f.close()
if __name__ == "__main__":
    main()
