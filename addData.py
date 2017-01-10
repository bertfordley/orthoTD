def addData(InDirectory, OutDirectory, genelist_file):
    '''
    for getting all data associated with genes in the gene list
    :param genelist:
    :return:
    '''
    import os, ntpath, glob, re

    with open(genelist_file, 'r') as genes:
        genelist = genes.read().splitlines()

    fileName = os.path.splitext(ntpath.basename(genelist_file))[0]
    outfile = open(OutDirectory + "/" + fileName + "_allData.txt", 'w')

    for gene in genelist:
        idGene = re.split(r'\t', gene.rstrip())
        print idGene[0], idGene[1]
        outfile.write("\n%s\t%s\n" % (idGene[0], idGene[1]))
        header = ["Gene 1", "Gene 2", "Spearman Coeff", "p-val", "Time Delay"]
        orthogroup = None
        allSpeciesGenes = []
        # get 3 species group for each gene from _Tables folder
        # also get gene IDs for orthologues
        tables_folder = InDirectory + "_Tables_0.7_Threshold"
        for file in glob.glob(os.path.join(tables_folder, "*.txt")):
            orthogroup = os.path.splitext(ntpath.basename(file))[0]
            with open(file, 'r') as f:
                tables = f.read().splitlines()
                if any(idGene[1].upper() in line for line in [x.upper() for x in tables]):
                    outfile.write(orthogroup[:-6] + "\n")
                    outfile.write("\t".join(header) + "\n")
                    for line in tables:
                        group = re.split(r'\t', line.rstrip())
                        if idGene[1].upper() in [x.upper() for x in group]:
                            # outfile.write("\t".join(group) + "\n")
                            [allSpeciesGenes.append(geneName) for geneName in group[:2] if geneName not in
                             allSpeciesGenes]
        for file in glob.glob(os.path.join(tables_folder, "*.txt")):
            with open(file, 'r')  as f:
                tables = f.read().splitlines()
                for line in tables:
                    splitline = re.split('\t', line.rstrip())
                    for gene1 in allSpeciesGenes:
                        for gene2 in allSpeciesGenes:
                            if gene1 != gene2 and gene1 == splitline[0] and gene2 == splitline[1]:
                                outfile.write(line + "\n")

        # get fpkm for each gene
        outfile.write("FPKM expression" + "\n")
        fpkmTimePoints = ['Transcript', '8am', '10am', '12am', '2pm', '4pm', '6pm', '8pm', '10pm', '12pm','2am',
                          '4am', '6am']
        outfile.write("\t".join(fpkmTimePoints) + "\n")
        for file in glob.glob(os.path.join(InDirectory, "*.txt")):
            if not orthogroup:
                orthogroup = os.path.splitext(ntpath.basename(file))[0]
                outfile.write(orthogroup[:-6] + "\n")
            with open(file, 'r') as f:
                fpkm = f.read().splitlines()
                for line in fpkm:
                    for ortholog in allSpeciesGenes:
                        if any(x in line for x in [ortholog, ortholog.upper(), ortholog.lower()]):
                            outfile.write(line + "\n")

        # get z-score for each gene
        outfile.write("Z-Score expression" + "\n")
        zTimePoints = ['Transcript', '7am', '8am', '9am', '10am', '11am', '12am','1pm', '2pm', '3pm', '4pm', '5pm',
                      '6pm', '7pm', '8pm', '9pm', '10pm', '11pm', '12pm', '1am', '2am', '3am', '4am', '5am', '6am']
        outfile.write("\t".join(zTimePoints) + "\n")
        for file in glob.glob(os.path.join(InDirectory + "_zscores", "*.txt")):
            with open(file, 'r') as f:
                zscore = f.read().splitlines()
                for line in zscore:
                    for ortholog in allSpeciesGenes:
                        if any(x in line for x in [ortholog, ortholog.upper(), ortholog.lower()]):
                            outfile.write(line + "\n")

    outfile.close()