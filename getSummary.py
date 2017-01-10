
def getSummary(InDirectory):
    '''
    creates summary file from outputs of findShifts function
    summaries include number of genes per species per shift group
    number of Kalanchoe TFs per shift group
    number of orthogroups per shift group
    '''
    import os, glob, re, ntpath

    outfile = open(InDirectory + "/summary.txt", "w")

    folders = [name for name in os.listdir(InDirectory) if os.path.isdir(InDirectory + "/" + name)]

    for folder in folders:
        print folder
        geneFile =  glob.glob(InDirectory + "/" + folder + "/*Redundant-Genes.txt")
        shiftFile = glob.glob(InDirectory + "/" + folder + "/*-Shift.txt")
        infoFile = glob.glob(InDirectory + "/" + folder + "/*-Info.txt")
        fileName = os.path.splitext(ntpath.basename(geneFile[0]))[0]
        numbTF = 0

        with open(infoFile[0], "r") as infoInfo:
            infoList = infoInfo.read().splitlines()
            geneInfo = []
            for infoLine in infoList:
                info = re.split(r"\t", infoLine.rstrip())
                geneInfo.append(info)
            for line in geneInfo:
                if line[1] != "notTF":
                    numbTF += 1

        with open(geneFile[0], "r") as file1:
            kafeSingles = []
            ancoSingles = []
            arthSingles = []
            geneList = file1.read().splitlines()
            for pair in geneList:
                genes = re.split(r'\t', pair.rstrip())
                kafeSingles.append(genes[0])
                ancoSingles.append(genes[1])
                arthSingles.append(genes[2])
            kafeSingles = list(set(kafeSingles))
            ancoSingles = list(set(ancoSingles))
            arthSingles = list(set(arthSingles))
            outfileKafe = open(InDirectory + "/" + folder + "/Nonredundant-Kafe-Genes-List.txt", 'w')
            outfileAnco = open(InDirectory + "/" + folder + "/Nonredundant-Anco-Genes-List.txt", 'w')
            outfileArth = open(InDirectory + "/" + folder + "/Nonredundant-Arth-Genes-List.txt", 'w')
            for k in kafeSingles:
                outfileKafe.write('%s\n' % k)
            outfileKafe.close()
            for an in ancoSingles:
                outfileAnco.write('%s\n' % an)
            outfileAnco.close()
            for ar in arthSingles:
                outfileArth.write('%s\n' % ar)
            outfileKafe.close()
            outfile.write('\n%s\nKafe genes\t%s\n# of Kafe TF\t%d\nAnco genes\t%s\nArth genes\t%s\n' % (fileName,
                                                                                                        len(kafeSingles),
                                                                                                        numbTF,
                                                                                                        len(ancoSingles),
                                                                                                        len(arthSingles)))
        with open(shiftFile[0], 'r') as file2:
            orthoTest = False
            orthocount = 0
            orthos = file2.read().splitlines()
            for line in orthos:
                group = re.split(r'\t', line.rstrip())
                if len(group) > 1:
                    if orthoTest == False:
                        orthocount += 1
                        orthoTest = True
                else:
                    orthoTest = False
            outfile.write('Orthogroups\t%s\n' % orthocount)

    outfile.close()
