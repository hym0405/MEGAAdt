#!/usr/bin/env python
import os
import sys

def rc(seq):
    pool = {"A": "T", "T": "A", "C": "G", "G": "C"}
    out = ""
    for e in seq:
        out = pool[e] + out
    return out

def rev(seq):
    out = ""
    for e in seq:
        out = e + out
    return out

def main():
    In = sys.argv[1]
    anchorRef = sys.argv[2]
    ref = sys.argv[3]
    Out = sys.argv[4]
    Out2 = sys.argv[5]
    sampleL = sys.argv[6]
    cutoff_mis = int(sys.argv[7])
    cutoff_gap = int(sys.argv[8])
    cutoff_length_sd = int(sys.argv[9])

    f = open(anchorRef, "rb")
    data = f.readlines()
    f.close()
    ref_pool_sample = {}
    for each in data[1:]:
        each = each[:-1]
        tmp = each.split("\t")
        ref_pool_sample[tmp[0]] = int(tmp[3])
    reference_length = ref_pool_sample[sampleL]
    cutoff_length_low = reference_length - cutoff_length_sd
    cutoff_length_high = reference_length + cutoff_length_sd

    f = open(ref,"rb")
    data = f.readlines()
    f.close()
    seq_pool = {}
    for i in range(len(data) / 4):
        tmpLabel = data[4 * i][1:-1]
        tmpSeq = data[4 * i + 1][:-1]
        tmpDum = data[4 * i + 2][:-1]
        tmpQuality = data[4 * i + 3][:-1]
        seq_pool[tmpLabel] = [tmpSeq, tmpDum, tmpQuality]

    f = open(In, "rb")
    data = f.readlines()
    f.close()
    f = open(Out, "w")
    f2 = open(Out2, "w")
    for each in data[1:]:
        each = each[:-1]
        tmp = each.split("\t")
        label, fstart, fend, fmis, fgap, rstart, rend, rmis, rgap, fstartRC, fendRC, fmisRC, fgapRC, rstartRC, rendRC, rmisRC, rgapRC = tmp
        fstart = int(fstart)
        fend = int(fend)
        rstart = int(rstart)
        rend = int(rend)
        fstartRC = int(fstartRC)
        fendRC = int(fendRC)
        rstartRC = int(rstartRC)
        rendRC = int(rendRC)
        flag = 0
        flag_RC = 0
        if int(fmis) <= cutoff_mis and int(fgap) <= cutoff_gap and int(rmis) <= cutoff_mis and int(rgap) <= cutoff_gap and fstart < rstart:
            flag = 1
        if int(fmisRC) <= cutoff_mis and int(fgapRC) <= cutoff_gap and int(rmisRC) <= cutoff_mis and int(rgapRC) <= cutoff_gap and fstartRC < rstartRC:
            flag_RC = 1

        if flag == 1 and flag_RC == 0:
            strand = "+"
            seqStart = fstart
            seqEnd = rend
            tmpSeq, tmpDum, tmpQuality = seq_pool[label]
            finalSeq = tmpSeq[(seqStart - 1): seqEnd]
            finalQuality = tmpQuality[(seqStart - 1): seqEnd]
            p5bc = tmpSeq[(seqStart - 13): (seqStart - 1)]
            if (seqEnd + 12) <= len(tmpSeq):
                p3bc = tmpSeq[seqEnd:(seqEnd + 12)]
            else:
                p3bc = ""

            length = len(finalSeq)
            if length >= cutoff_length_low and length <= cutoff_length_high:
                toWrite = [label, strand, p5bc, p3bc, length]
                f.writelines(["\t".join([str(e) for e in toWrite]) + os.linesep])
                f2.writelines([">" + label + os.linesep])
                f2.writelines([finalSeq + os.linesep])
                f2.writelines([tmpDum + os.linesep])
                f2.writelines([finalQuality + os.linesep])
        if flag == 0 and flag_RC == 1:
            strand = "-"
            seqStart = fstartRC
            seqEnd = rendRC
            tmpSeq, tmpDum, tmpQuality = seq_pool[label]
            finalSeq = rc(tmpSeq)[(seqStart - 1): seqEnd]
            finalQuality = rev(tmpQuality)[(seqStart - 1): seqEnd]
            p5bc = rc(tmpSeq)[(seqStart - 13): (seqStart - 1)]

            if (seqEnd + 12) <= len(tmpSeq):
                p3bc = rc(tmpSeq)[seqEnd:(seqEnd + 12)]
            else:
                p3bc = ""

            length = len(finalSeq)
            if length >= cutoff_length_low and length <= cutoff_length_high:
                toWrite = [label, strand, p5bc, p3bc, length]
                f.writelines(["\t".join([str(e) for e in toWrite]) + os.linesep])
                f2.writelines([">" + label + os.linesep])
                f2.writelines([finalSeq + os.linesep])
                f2.writelines([tmpDum + os.linesep])
                f2.writelines([finalQuality + os.linesep])
    f.close()
    f2.close()

if __name__ == "__main__":
    main()

