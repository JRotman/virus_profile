# Note :
# 1) make sure the same read is not counted multiple times to the same genome


import sys
import random
import csv
import os
import numpy
import matplotlib
# matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt


import argparse
from collections import Counter


ap = argparse.ArgumentParser()
ap.add_argument('virusMegablastSample', help='virusMegablastSample')
ap.add_argument('readRength', help='')

args = ap.parse_args()


rLength = int(args.readRength)

# NR - number of reads
NR = 20
SP_DIFF_MAX = rLength - 20
COVERAGE_MIN = int(1.5 * rLength)

base = os.path.basename(args.virusMegablastSample)
prefix = os.path.splitext(base)[0]

print "Read ", args.virusMegablastSample

out = prefix + "_filtered.csv"

# print "Save to",out


dict = {}
genomesList = []
genomesList2 = set()
genomeReads = {}
genomeReads2 = {}


genomes = {}

genomeL = {}
genomeR = {}

genomeCoverage = {}
genomeStartPosition = {}
genomeNumberReads = {}

k = 0

# HWI-ST1148:179:C4BAKACXX:5:1101:11944:1949/1    gi|8486122|ref|NC_002016.1|
# 94.68   94      5       0     7   100     750     843     7e-35    147


# get list with genomes
check = 1
with open(args.virusMegablastSample, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for line in reader:
        genome = line[1]
        genomesList.append(genome)

for g in genomesList:
    genomeReads[g] = set()


# get reads from each genome
with open(args.virusMegablastSample, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for line in reader:
        genome = line[1]
        genomeReads[genome].add(line[0])

genomesList = list(set(genomesList))


# if there are more then 10 reads assigned onto the genome prepare the coordinates of LEFT and RIGHT
# coordinates of the covered region
for g in genomesList:

        k += 1
        genomeL[g] = []
        genomeR[g] = []


# dict = {}     ???????


for g in genomesList:
    for r in genomeReads[g]:
        dict[r] = []


# make dict for each read
with open(args.virusMegablastSample, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for line in reader:
        read = line[0]
        genome = line[1]
        dict[read].append(line)


# get LEFT and RIGHT coordinates of the covered region

dict2 = {}
single_map_dict = {}
for g in genomesList:
    dict2[g] = []
    single_map_dict[g] = []

nMultimapped = 0
random.seed()
for key, value in dict.iteritems():
    if len(value) > 1:

        gList = set()
        for v in value:
            if v[1] not in gList:
                dict2[v[1]].append(v)

        randIndex = random.randint(0, len(value) - 1)
        dict2[value[randIndex][1]].append(value[randIndex])

        nMultimapped += 1

    else:
        dict2[value[0][1]].append(value[0])

print "Number of multimapped reads: ", nMultimapped

checkList = set()
mMap = 0

dictGenome = {}
for key, value in dict2.items():
    dictGenome[key] = [0, 0, 0, [], [], ""]


# ['HWI-ST1148:179:C4BAKACXX:5:1101:17968:1893/1', 'gi|752901102|ref|NC_026427.1|',
#  '95.65', '92', '4', '0', '7', '98', '798', '889', '2e-35', ' 148']

for key, value in dict2.items():

    lGlobal = []
    rGlobal = []

    if len(value) > NR:

        for v in value:

            l = int(v[8])
            r = int(v[9])

            if l > r:
                l, r = r, l

            lGlobal.append(int(l))
            rGlobal.append(int(r))

        # dictGenome:
        # [ length, leftmost point, rightmost point, starting point array, coverage array, key ]
        dictGenome[key][0] = len(value)
        dictGenome[key][5] = line[1]
        dictGenome[key][1] = min(lGlobal)
        dictGenome[key][2] = max(rGlobal)
        dictGenome[key][3] = [0]*(dictGenome[key][2]+1)
        dictGenome[key][4] = [0]*(dictGenome[key][2]+1)

figure_number = 1

pathway = "eupathPlots/ETAG007/"

# iterate through each virus
for key, value in dict2.items():
    startsFull = []
    endsFull = []
    startIdx = []
    endIdx = []
    spPlots = []
    cPlots = []
    accept = 0
    overall_accept = 0

    # check that there are at least NR (20) reads
    if len(value) > NR:
        print str(key) + " : " + str(len(value))
        for v in value:
            l = int(v[8])
            r = int(v[9])
            if l > r:
                l, r = r, l

            dictGenome[key][3][l] += 1

            for i in range(l, r+1):
                dictGenome[key][4][i] += 1

        startIndex = 0
        endIndex = 0
        zeroCounter = 0
        lastEnd = 0
        starts = []
        ends = []
        first = 0
        for n in range(len(dictGenome[key][4])):
            if dictGenome[key][4][n] == 0:
                if startIndex > lastEnd and zeroCounter > 0:
                    sIdx = []
                    eIdx = []
                    if first == 0:
                        spPlots.append(dictGenome[key][3][startIndex - 1:endIndex + 1])
                        cPlots.append(dictGenome[key][4][startIndex - 1:endIndex + 1])
                        for x in range(startIndex - 1, endIndex + 1):
                            sIdx.append(x)
                        first = 1
                    else:
                        spPlots.append(dictGenome[key][3][startIndex:endIndex + 1])
                        cPlots.append(dictGenome[key][4][startIndex:endIndex + 1])
                        for x in range(startIndex, endIndex + 1):
                            sIdx.append(x)
                    startIdx.append(sIdx)
                    starts.append(startIndex)
                    ends.append(endIndex + 1)
                    startsFull.append(starts)
                    endsFull.append(ends)
                    lastEnd = endIndex
                    startIndex = 0
                    zeroCounter = 0
                else:
                    zeroCounter += 1
            elif startIndex == 0 and endIndex != 0:
                startIndex = n
                endIndex = n
            else:
                endIndex = n

        if dictGenome[key][4][n] != 0:
            sIdx = []
            eIdx = []
            if first == 0:
                spPlots.append(dictGenome[key][3][startIndex - 1:endIndex + 1])
                cPlots.append(dictGenome[key][4][startIndex - 1:endIndex + 1])
                for x in range(startIndex - 1, endIndex + 1):
                    sIdx.append(x)
                first = 1
            else:
                spPlots.append(dictGenome[key][3][startIndex:endIndex + 1])
                cPlots.append(dictGenome[key][4][startIndex:endIndex + 1])
                for x in range(startIndex, endIndex + 1):
                    sIdx.append(x)
            startIdx.append(sIdx)

        maxC = 0
        index = 0

        for pl in range(len(spPlots)):
            print key

            if len(cPlots[pl]) >= COVERAGE_MIN:
                print "Length Covered: " + str(len(cPlots[pl]))
                spDiff = []
                currIdx = 0
                lastIdx = 0
                for x in range(len(spPlots[pl])):
                    if spPlots[pl][lastIdx] == 0:
                        if spPlots[pl][x] > 0:
                            lastIdx = x
                    elif spPlots[pl][x] > 0:
                        # print "difference: " + str(x) + str(lastIdx)
                        if lastIdx != x:
                            spDiff.append(startIdx[pl][x] - startIdx[pl][lastIdx])
                            lastIdx = x
                spAverage = sum(spDiff) / len(spDiff)
                # print spDiff
                print "Average SP Difference: " + str(spAverage)
                if spAverage <= SP_DIFF_MAX:
                    print "ACCEPTED"
                    accept = 1
                    overall_accept = 1
                    if len(cPlots[pl]) > maxC:
                        maxC = len(cPlots[pl])
                        index = pl

            fig = plt.figure(figure_number, figsize=(10.24, 10.24))

            if accept == 1:
                fig.suptitle('%s_ACCEPT' % key, fontsize=16, fontweight='bold')
            else:
                fig.suptitle('%s_REJECT' % key, fontsize=16, fontweight='bold')

            ax = fig.add_subplot(211)
            ax.set_title('Starting Point Distribution', fontweight='bold')
            ax.set_xlabel('genome index')
            ax.set_ylabel('number of starting points')
            ax.plot(startIdx[pl], spPlots[pl])

            ax = fig.add_subplot(212)
            ax.set_title('Coverage', fontweight='bold')
            ax.set_xlabel('genome index')
            ax.set_ylabel('number of covering reads')
            ax.plot(startIdx[pl], cPlots[pl])

            if accept == 1:
                plt.savefig(pathway + "smallPlots/%s_%d.png" % (key, pl))
            # else:
            #    plt.savefig(pathway + "smallPlots/%s_%d.png" % (key, pl))
            plt.close(figure_number)
            figure_number += 1

    fig = plt.figure(figure_number, figsize=(10.24, 10.24))

    if overall_accept == 1:
        fig.suptitle('%s_ACCEPT' % key, fontsize=16, fontweight='bold')
    else:
        fig.suptitle('%s_REJECT' % key, fontsize=16, fontweight='bold')

    ax = fig.add_subplot(211)
    ax.set_title('Starting Point Distribution', fontweight='bold')
    ax.set_xlabel('genome index')
    ax.set_ylabel('number of starting points')
    ax.plot(dictGenome[key][3])

    ax = fig.add_subplot(212)
    ax.set_title('Coverage', fontweight='bold')
    ax.set_xlabel('genome index')
    ax.set_ylabel('number of covering reads')
    ax.plot(dictGenome[key][4])

    if overall_accept == 1:
        plt.savefig(pathway + "largePlots/%s.png" % key)
    # else:
    #     plt.savefig(pathway + "largePlots/%s.png" % key)
    plt.close(figure_number)
    figure_number += 1
