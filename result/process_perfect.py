#!/usr/bin/python
import sys
import glob
import re

# get number of solutions
solCount = {}
for filename in glob.glob("sims_r*_m*_n*.res"):
    res = re.match(r'sims_r([0-9]*)_m([0-9]*)_n([0-9]*).res', filename)
    r = int(res.group(1))
    m = int(res.group(2))
    n = int(res.group(3))

    with open(filename) as f:
        string = f.readline()
        if string == "":
            continue
        nSols = int(string.split()[0])
        if r not in solCount:
            solCount[r] = {}
        if m not in solCount[r]:
            solCount[r][m] = {}
        solCount[r][m][n] = nSols

recallMostRep = {}
for filename in glob.glob("sims_r*_m*_n*.recall"):
    res = re.match(r'sims_r([0-9]*)_m([0-9]*)_n([0-9]*).recall', filename)
    if not(res):
        continue
    r = int(res.group(1))
    m = int(res.group(2))
    n = int(res.group(3))

    with open(filename) as f:
        string = f.readline()
        if string == "":
            continue
        concordance = float(string.split()[0])
        if r not in recallMostRep:
            recallMostRep[r] = {}
        if m not in recallMostRep[r]:
            recallMostRep[r][m] = {}
        recallMostRep[r][m][n] = concordance

recallMedian = {}
for filename in glob.glob("sims_r*_m*_n*.median.recall"):
    res = re.match(r'sims_r([0-9]*)_m([0-9]*)_n([0-9]*).median.recall', filename)
    if not(res):
        continue
    r = int(res.group(1))
    m = int(res.group(2))
    n = int(res.group(3))

    with open(filename) as f:
        string = f.readline()
        if string == "":
            continue
        concordance = float(string.split()[0])
        if r not in recallMedian:
            recallMedian[r] = {}
        if m not in recallMedian[r]:
            recallMedian[r][m] = {}
        recallMedian[r][m][n] = concordance

print '"r","$m$","$n$","coverage","\#solutions","most representative recall","median recall"'
for r in solCount:
    for m in solCount[r]:
        for n in solCount[r][m]:
            print ",".join(map(str, [r, m, n, "inf", solCount[r][m][n], recallMostRep[r][m][n], recallMedian[r][m][n]]))

