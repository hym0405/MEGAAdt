#!/usr/bin/env python
import os
import sys
from tqdm import tqdm
def statMutants(seq, refSeq):
    refCount = 0
    mutantsPool = []
    for i in range(len(seq)):
        tmpRefBase = refSeq[i]
        tmpSeqBase = seq[i]
        if tmpRefBase != "-":
            refCount += 1
            if tmpSeqBase != "-":
                if tmpSeqBase == tmpRefBase:
                    continue
                else:
                    mutantsPool.append("SNP_" + tmpRefBase + str(refCount) + tmpSeqBase)
                    continue
            else:
                mutantsPool.append("DEL_" + tmpRefBase + str(refCount))
                continue
        else:
            if tmpSeqBase != "-":
                mutantsPool.append("INS_" + str(refCount) + tmpSeqBase)
                continue
            else:
                continue
    return mutantsPool
def main():
    In = sys.argv[1]
    refMutants = sys.argv[2]
    Out = sys.argv[3]
    label = sys.argv[4]

    f = open(refMutants,"rb")
    data = f.readlines()
    f.close()

    ref_pool = {}
    for each in data[1:]:
        each = each[:-1]
        tmp = each.split("\t")
        if tmp[6] == label:
            ref_seq = tmp[7]
        ref_pool[tmp[6]] = tmp[8]


    f = open(In, "rb")
    data_raw = f.readlines()
    f.close()
    data_raw = [">" + label + os.linesep, ref_pool[label] + os.linesep] + data_raw
    fw = open(Out, "w")
    fw.writelines(["seqID\tmutations" + os.linesep])

    for i in tqdm(range(len(data_raw) / 2)):
        label = data_raw[2 * i][1:-1]
        seq = data_raw[2 * i + 1][:-1]
        f = open(Out + ".tmp.fasta", "w")
        f.writelines([">reference" + os.linesep])
        f.writelines([ref_seq + os.linesep])
        f.writelines([">reads" + os.linesep])
        f.writelines([seq + os.linesep])
        f.close()
        os.system("muscle -in " + Out + ".tmp.fasta -out " + Out + ".tmp.MSA.fasta &> /dev/null")
        f = open(Out + ".tmp.MSA.fasta", "rb")
        data = f.readlines()
        f.close()

        Pool = {}
        tmpLabel = ""
        tmpSeq = ""
        flag = 0
        for each in data:
            if each[0] == ">":
                if flag == 1:
                    Pool[tmpLabel] = tmpSeq
                tmpLabel = each[1:-1]
                tmpSeq = ""
            else:
                tmpSeq += each[:-1]
            flag = 1
        Pool[tmpLabel] = tmpSeq
        outSeqMSA = Pool["reads"]
        refSeqMSA = Pool["reference"]
        mutants = statMutants(outSeqMSA, refSeqMSA)
        fw.writelines([label + "\t" + ",".join(mutants) + os.linesep])
        os.system("rm -f " + Out + ".tmp.fasta")
        os.system("rm -f " + Out + ".tmp.MSA.fasta")

    fw.close()
if __name__ == "__main__":
    main()
