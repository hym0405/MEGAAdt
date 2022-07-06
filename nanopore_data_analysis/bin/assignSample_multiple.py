#!/usr/bin/env python
import os
import sys
def checkSeq(seq1,seq2):
    if len(seq2) == 0:
        return 0
    if len(seq1) != len(seq2):
        return 9999
    else:
        count = 0
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                count += 1
        return count

def checkBC(BC1, BC2, BC_pool, cutoff):
    flag = 0
    out_pool = []
    for e in BC_pool:
        tmpLabel, tmpBC1, tmpBC2 = e
        delta1 = checkSeq(BC1,tmpBC1)
        delta2 = checkSeq(BC2,tmpBC2)
        if delta1 <= cutoff and delta2 <= cutoff:
            flag = 1
            out_pool.append([delta1 + delta2, tmpLabel])
    if flag == 0:
        return "NA"
    else:
        out_pool.sort()
        return out_pool[0][1]

def main():
    In = sys.argv[1]
    ref = sys.argv[2]
    Out = sys.argv[3]
    cutoff = int(sys.argv[4])
    sample_L = sys.argv[5]

    f = open(ref,"rb")
    data = f.readlines()
    f.close()
    Pool = []
    for each in data[1:]:
        each = each[:-1]
        tmp = each.split("\t")
        if tmp[0] == sample_L:
            Pool.append([tmp[6], tmp[4], tmp[5]])

    f = open(In, "rb")
    data = f.readlines()
    f.close()
    f = open(Out, "w")
    stat = {}
    for e in Pool:
        stat[e[0]] = 0
    stat["NA"] = 0
    for each in data:
        each = each[:-1]
        tmp = each.split("\t")
        tmpLabel = checkBC(tmp[2], tmp[3], Pool, cutoff)
        f.writelines([each + "\t" + tmpLabel + os.linesep])
        stat[tmpLabel] += 1
    f.close()
    for e in Pool:
        print e[0], stat[e[0]]
    print "NA", stat["NA"]
if __name__ == "__main__":
    main()



    
