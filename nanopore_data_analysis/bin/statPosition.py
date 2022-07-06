#!/usr/bin/env python
import os
import sys
def main():
    In = sys.argv[1]
    Out = sys.argv[2]
    f = open(In, "rb")
    data = f.readlines()
    f.close()
    Pool = {}
    Pool_star = []

    tmp = data[1][:-1].split("\t")
    Pool_star = tmp[1].split(",")

    for each in data[1:]:
        each = each[:-1]
        tmp = each.split("\t")
        tmp2 = tmp[1].split(",")
        for e in tmp2:
            Pool[e] = 0

    for each in data[2:]:
        each = each[:-1]
        tmp = each.split("\t")
        tmp2 = tmp[1].split(",")
        for e in tmp2:
            Pool[e] += 1

    Pool_list = []
    for e in Pool.keys():
        Pool_list.append([Pool[e], e])
    Pool_list.sort()
    Pool_list.reverse()

    f = open(Out, "w")
    f.writelines(["Position\tMutation\tcount\ttotalCount\tifTargets" + os.linesep])
    for each in Pool_list:
        count, mutation = each
        totalCount = len(data[1:]) - 1
        if mutation in Pool_star:
            iftarget = "Yes"
        else:
            iftarget = ""
        if count > 1 or iftarget == "Yes":
            toWrite = [mutation, str(count), str(totalCount), iftarget]
            tmp2 = mutation.split("_")
            if len(tmp2) > 1:
                f.writelines([tmp2[2] + "\t" + "\t".join(toWrite) + os.linesep])
    f.close()


if __name__ == "__main__":
    main()


