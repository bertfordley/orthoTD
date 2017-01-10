def addKafeInfo(InDirectory):
    '''
    pulls data from my masterfile for Kalanchoe and creates new file with only specified
    kalanchoe genes and their respective info (e.g. GO)

    '''
    import os, glob, re, ntpath

    folders = [name for name in os.listdir(InDirectory) if os.path.isdir(InDirectory + "/" + name)]

    masterDict = {}
    masterFile = open("../Kaladp_masterTable_20160314_all.txt", "r")
    masterFileList = masterFile.read().splitlines()
    masterFile.close()
    for line in masterFileList:
        geneLine = re.split(r'\t', line.rstrip())
        masterDict[geneLine[0]] = geneLine[1:]

    for folder in folders:
        print folder
        kafeGeneFile =  glob.glob(InDirectory + "/" + folder + "/*Kafe-Genes-List.txt")
        kafeFileName = os.path.splitext(ntpath.basename(kafeGeneFile[0]))[0]

        with open(kafeGeneFile[0], "r") as kafeFile:
            outfileKafe = open(InDirectory + "/" + folder + "/" + kafeFileName + "-Info.txt", 'w')
            kafeGeneList = kafeFile.read().splitlines()
            for kafeGene in kafeGeneList:
                geneInfo =  [[key,value] for key, value in masterDict.items() if kafeGene in key]
                geneID = geneInfo[0][0]
                TF = geneInfo[0][1][0].split(":")[0]
                geneDesc = geneInfo[0][1][1].split(":")[0]
                outfileKafe.write("%s\t%s\t%s\n" % (geneID, TF, geneDesc))
            outfileKafe.close()
