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
    f.writelines(["seqID\tmutations" + os.linesep])
    for each in data[1:]:
        each = each[:-1]
        tmp = each.split("\t")
        label = tmp[0]
        mutations = tmp[1].split(",")
        mutation_pool = {}
        mutation_list = []
        for e in mutations:
            tmp2 = e.split("_")
            if tmp2[0] == "SNP":
                pos = tmp2[1][1:-1]
                mutation_list.append(e + "_" + pos)
            if tmp2[0] == "DEL":
                pos = tmp2[1][1:]
                mutation_list.append(e + "_" + pos)
            elif tmp2[0] == "INS":
                loc = tmp2[1][:-1]
                if loc in mutation_pool.keys():
                    mutation_pool[loc] += tmp2[1][-1]
                else:
                    mutation_pool[loc] = tmp2[1][-1]
        for e in mutation_pool.keys():
            mutation_list.append("INS_" + e + mutation_pool[e] + "_" + e)
        f.writelines([label + "\t" + ",".join(mutation_list) + os.linesep])
    f.close()
if __name__ == "__main__":
    main()

