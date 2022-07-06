#!/usr/bin/env python
import os
import sys
from tqdm import tqdm
def rc(seq):
    pool = {"A": "T", "T": "A", "C": "G", "G": "C"}
    out = ""
    for e in seq:
        out = pool[e] + out
    return out

def statSeq(fileP, anchor_length):
    f = open(fileP, "rb")
    data = f.readlines()
    f.close()
    tmpSeq = ""
    tmpLabel = ""
    flag = 0
    seqPool = {}
    for each in data:
        if each[0] == ">":
            if flag == 1:
                seqPool[tmpLabel] = tmpSeq
            tmpLabel = each[1:-1]
            tmpSeq = ""
        else:
            tmpSeq += each[:-1]
        flag = 1
    seqPool[tmpLabel] = tmpSeq
    seq_anchor = seqPool["anchor"]
    seq_reads = seqPool["reads"]

    gaps_count = 0
    mis_count = 0
    reads_start = 0
    reads_end = 0
    reads_current = 0
    anchor_tmp_len = 0
    flag = 0
    for i in range(len(seq_anchor)):
        if seq_reads[i] != "-":
            reads_current += 1

        if seq_anchor[i] == "-" and flag == 0:
            continue

        if seq_anchor[i] != "-" and flag == 0:
            flag = 1
            reads_start = reads_current
            anchor_tmp_len = 1
            continue

        if seq_anchor[i] != "-" and flag == 1:
            anchor_tmp_len += 1

            if seq_anchor[i] != seq_reads[i]:
                mis_count += 1

            if anchor_tmp_len == anchor_length:
                flag = 2
                reads_end = reads_current
                continue

        if seq_anchor[i] == "-" and flag == 1:
            gaps_count += 1
            continue

        if seq_anchor[i] != "-" and flag == 2:
            print "Error"
            continue

    return reads_start, reads_end, mis_count, gaps_count

def main():
    In = sys.argv[1]
    ref = sys.argv[2]
    Out = sys.argv[3]
    sampleL = sys.argv[4]
    tempDir = sys.argv[5]
    f = open(ref, "rb")
    data = f.readlines()
    f.close()
    ref_pool_sample = {}
    for each in data[1:]:
        each = each[:-1]
        tmp = each.split("\t")
        ref_pool_sample[tmp[0]] = [tmp[1], tmp[2]]
    ref_pool = {"start": ref_pool_sample[sampleL][0], \
                "end": ref_pool_sample[sampleL][1]}

    f = open(In, "rb")
    data = f.readlines()
    f.close()
    fw = open(Out, "w")
    fw.writelines(["reads\tanchor_F_start\tanchor_F_end\tanchor_F_mismatch\tanchor_F_gap\t" + \
                "anchor_R_start\tanchor_R_end\tanchor_R_mismatch\tanchor_R_gap\t" + \
                "anchor_F_start_RC\tanchor_F_end_RC\tanchor_F_mismatch_RC\tanchor_F_gap_RC\t" + \
                "anchor_R_start_RC\tanchor_R_end_RC\tanchor_R_mismatch_RC\tanchor_R_gap_RC" + os.linesep])

    for i in tqdm(range(len(data) / 2)):
        label = data[2 * i][1:-1]
        seq = data[2 * i + 1][:-1]

        f = open(tempDir + ".tmp.fasta", "w")
        f.writelines([">anchor" + os.linesep])
        f.writelines([ref_pool["start"] + os.linesep])
        f.writelines([">reads" + os.linesep])
        f.writelines([seq + os.linesep])
        f.close()
        os.system("muscle -in " + tempDir + ".tmp.fasta -out " + tempDir + ".tmp.MSA.fasta &> /dev/null")
        toWrite = [label,] + [e for e in statSeq(tempDir + ".tmp.MSA.fasta", len(ref_pool["start"]))]
        os.system("rm -f " + tempDir + ".tmp.fasta")
        os.system("rm -f " + tempDir + ".tmp.MSA.fasta")

        f = open(tempDir + ".tmp.fasta", "w")
        f.writelines([">anchor" + os.linesep])
        f.writelines([ref_pool["end"] + os.linesep])
        f.writelines([">reads" + os.linesep])
        f.writelines([seq + os.linesep])
        f.close()
        os.system("muscle -in " + tempDir + ".tmp.fasta -out " + tempDir + ".tmp.MSA.fasta &> /dev/null")
        toWrite += [e for e in statSeq(tempDir + ".tmp.MSA.fasta", len(ref_pool["end"]))]
        os.system("rm -f " + tempDir + ".tmp.fasta")
        os.system("rm -f " + tempDir + ".tmp.MSA.fasta")


        f = open(tempDir + ".tmp.fasta", "w")
        f.writelines([">anchor" + os.linesep])
        f.writelines([ref_pool["start"] + os.linesep])
        f.writelines([">reads" + os.linesep])
        f.writelines([rc(seq) + os.linesep])
        f.close()
        os.system("muscle -in " + tempDir + ".tmp.fasta -out " + tempDir + ".tmp.MSA.fasta &> /dev/null")
        toWrite += [e for e in statSeq(tempDir + ".tmp.MSA.fasta", len(ref_pool["start"]))]
        os.system("rm -f " + tempDir + ".tmp.fasta")
        os.system("rm -f " + tempDir + ".tmp.MSA.fasta")

        f = open(tempDir + ".tmp.fasta", "w")
        f.writelines([">anchor" + os.linesep])
        f.writelines([ref_pool["end"] + os.linesep])
        f.writelines([">reads" + os.linesep])
        f.writelines([rc(seq) + os.linesep])
        f.close()
        os.system("muscle -in " + tempDir + ".tmp.fasta -out " + tempDir + ".tmp.MSA.fasta &> /dev/null")
        toWrite += [e for e in statSeq(tempDir + ".tmp.MSA.fasta", len(ref_pool["end"]))]
        os.system("rm -f " + tempDir + ".tmp.fasta")
        os.system("rm -f " + tempDir + ".tmp.MSA.fasta")

        fw.writelines(["\t".join([str(e) for e in toWrite]) + os.linesep])

    fw.close()
if __name__ == "__main__":
    main()
