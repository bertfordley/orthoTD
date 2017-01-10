'''

tool to call python functions in the terminal

Usage:
for getting spearman, p-values, and time delay
    - python MultiTool.py <InDirectory> ccorr <s>

for getting manuscript gene list
    - python MultiTool.py <InDirectory> listgen <genelist>

for getting diel shifts from Spearman/TD tables
    - python MultiTool.py <Directory> findshifts

for prepending line in file
    - python MultiTool.py <Directory> lineprepend
    * edit line in code below if needed

for getting summary data
    - python MultiTool.py <Directory> getsummary <InDirectory>

for getting info on Kalanchoe genes in Kafe-genes-list.txt file from Kalanchoe masterfile
    - python MultiTool.py <Directory> addkafeinfo <InDirectory>

for getting all data (i.e. FPKM expression data, z-score data, shift tables) for genes in gene list file
    - python MultiTool.py <Directory> adddata <gene_list>

Input Arguments:

:InDirectory    directory with files to be worked on with respective function
:function       either ccorr, listgen, getshifts, or lineprepend
:s              spline number
:gene_list       text doc with gene list of orthologs - each species is a column


'''

import os, sys, glob, ntpath, re

InDirectory = sys.argv[1]
useFunction = sys.argv[2]
functions = ["ccorr", "listgen", "findshifts", "lineprepend", "getsummary", "addkafeinfo", "adddata"]

if useFunction.lower() not in functions:
    print "Please indicate which function to use:"
    print [x for x in functions]
    print "\n"
else:

    if useFunction.lower() == "ccorr":

        sys.path.insert(0, '')
        from Ccorr import Ccorr
        s = sys.argv[3]

        OutDirectory = InDirectory + "_Tables"
        if not os.path.exists(OutDirectory):
            os.makedirs(OutDirectory)

        if not os.path.exists(InDirectory + "_zscores"):
            os.makedirs(InDirectory + "_zscores")

        for file in glob.glob(os.path.join(InDirectory, '*.txt')):
            fileName = ntpath.basename(file)
            print "Working on file: " + fileName
            Ccorr(file, s, OutDirectory, InDirectory)

    elif useFunction.lower() == "listgen":
        genelist = sys.argv[3]
        sys.path.insert(0, '')
        from listGen import listGen
        manulist = []
        with open(genelist) as f:
            manu_list = f.read().splitlines()
            for line in manu_list:
                trips = re.split(r'\t', line.rstrip())
                pairs = [[trips[0], trips[2]],
                         [trips[1], trips[2]],
                         [trips[0], trips[1]]]
                manulist.append(pairs)

        for file in glob.glob(os.path.join(InDirectory, '*.txt')):
            listGen(file, manulist)

    elif useFunction.lower() == "findshifts":

        sys.path.insert(0, '')
        from findShifts import findShifts

        OutDirectory = InDirectory + "_KeyShifts"
        if not os.path.exists(OutDirectory):
            os.makedirs(OutDirectory)
            os.makedirs(OutDirectory + "/Kafe-Anco-12hr-Shift")
            os.makedirs(OutDirectory + "/Kafe-Anco-6hr-Shift")
            os.makedirs(OutDirectory + "/Anco-12hr-Shift")
            os.makedirs(OutDirectory + "/Anco-6hr-Shift")
            os.makedirs(OutDirectory + "/Kafe-12hr-Shift")
            os.makedirs(OutDirectory + "/Kafe-6hr-Shift")
        print OutDirectory
        for file in glob.glob(os.path.join(InDirectory, '*.txt')):
            fileName = ntpath.basename(file)
            print "Working on file: " + fileName
            findShifts(file, OutDirectory, 0.7)

    elif useFunction.lower() == 'lineprepend':
        # incase you need to add column headers to orthogroup files
        sys.path.insert(0, '')
        from lineprepend import linePrepend

        line = "Transcript\t2h\t4h\t6h\t8h\t10h\t12h\t14h\t16h\t18h\t20h\t22h\t24h"

        for file in glob.glob(os.path.join(InDirectory, '*.txt')):
            fileName = ntpath.basename(file)
            linePrepend(file,line)

    elif useFunction.lower() == "addkafeinfo":
        # makes new file for Kafe genes and adds TF and description info
        sys.path.insert(0, '')
        from addKafeInfo import addKafeInfo

        addKafeInfo(InDirectory)

    elif useFunction.lower() == "getsummary":
        # gets total genes per species, number of orthogroups, and number of Kafe TFs for each shift group
        sys.path.insert(0, '')
        from getSummary import getSummary

        getSummary(InDirectory)

    elif useFunction.lower() == "adddata":
        sys.path.insert(0, '')
        from addData import addData

        OutDirectory = InDirectory + "_geneLists"
        if not os.path.exists(OutDirectory):
            os.makedirs(OutDirectory)

        addData(InDirectory, OutDirectory, sys.argv[3])

