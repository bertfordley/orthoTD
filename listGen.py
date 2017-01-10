def listGen(input_file, genelist):

    '''
    creates a shift table file for genes in the gene list file

    Kafe -> Arth     Anco -> Arth     Kafe -> Anco

    :param input_file: three column list
    :return: list of gene pairs

    '''

    import re

    outfile = open("New_List.txt",'a')

    with open(input_file) as f:
        filedata = f.read().splitlines()
        genedata = []
        for line in filedata[3:]:
            splitline = re.split(r'\t', line.rstrip())
            for trips in genelist:
                for pair in trips:
                    if splitline[0][5:] == pair[0] and splitline[1][5:] == pair[1]:
                        genedata.append(splitline)
        if genedata:
            KafeArth = None
            AncoArth = None
            KafeAnco = None
            for pairs in genedata:
                if pairs[0][:4] == "Kala" and pairs[1][:4] == "Arth":
                    KafeArth = pairs
                elif pairs[0][:4] == "Anco" and pairs[1][:4] == "Arth":
                    AncoArth = pairs
                elif pairs[0][:4] == "Kala" and pairs[1][:4] == "Anco":
                    KafeAnco = pairs
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (KafeArth[0][5:],AncoArth[0][5:],
                                                                                KafeArth[1][5:],
                                                                                KafeArth[4],KafeArth[2],KafeArth[3],
                                                                                AncoArth[4],AncoArth[2],AncoArth[3],
                                                                                KafeAnco[4],KafeAnco[2],KafeAnco[3]))
