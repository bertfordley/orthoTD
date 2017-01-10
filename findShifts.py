def findShifts(input_file, OutDirectory, threshold):

    '''
    finds 12 and 6 hour (+/- 2 hours) shifts between three species

    Kafe -> Arth     Anco -> Arth     Kafe -> Anco

    :param input_file: Spearman and time delay tables for ortholog groups
    :return: list of gene pairs

    '''

    # zscoreMatrix.to_csv(InDirectory + "_zscores/" + fileName + "_zScores.txt", sep='\t')

    import re, os, ntpath

    fileName = os.path.splitext(ntpath.basename(input_file))[0]
    outfile1 = open(OutDirectory + "/Kafe-Anco-12hr-Shift/Kafe-Anco-12hr-Shift.txt", 'a')
    outfile1.write('%s\n' % fileName)
    outfile1genes = open(OutDirectory + "/Kafe-Anco-12hr-Shift/Kafe-Anco-12hr-Shift-Redundant-Genes.txt", 'a')

    outfile2 = open(OutDirectory + "/Kafe-Anco-6hr-Shift/Kafe-Anco-6hr-Shift.txt", 'a')
    outfile2.write('%s\n' % fileName)
    outfile2genes = open(OutDirectory + "/Kafe-Anco-6hr-Shift/Kafe-Anco-6hr-Shift-Redundant-Genes.txt", 'a')

    outfile3 = open(OutDirectory + "/Anco-12hr-Shift/Anco-12hr-Shift.txt", 'a')
    outfile3.write('%s\n' % fileName)
    outfile3genes = open(OutDirectory + "/Anco-12hr-Shift/Anco-12hr-Shift-Redundant-Genes.txt", 'a')

    outfile4 = open(OutDirectory + "/Anco-6hr-Shift/Anco-6hr-Shift.txt", 'a')
    outfile4.write('%s\n' % fileName)
    outfile4genes = open(OutDirectory + "/Anco-6hr-Shift/Anco-6hr-Shift-Redundant-Genes.txt", 'a')

    outfile5 = open(OutDirectory + "/Kafe-12hr-Shift/Kafe-12hr-Shift.txt", 'a')
    outfile5.write('%s\n' % fileName)
    outfile5genes = open(OutDirectory + "/Kafe-12hr-Shift/Kafe-12hr-Shift-Redundant-Genes.txt", 'a')

    outfile6 = open(OutDirectory + "/Kafe-6hr-Shift/Kafe-6hr-Shift.txt", 'a')
    outfile6.write('%s\n' % fileName)
    outfile6genes = open(OutDirectory + "/Kafe-6hr-Shift/Kafe-6hr-Shift-Redundant-Genes.txt", 'a')

    with open(input_file) as f:
        filedata = f.read().splitlines()
        KafeArth = []
        AncoArth = []
        KafeAnco = []
        for line in filedata[3:]:
            splitline = re.split(r'\t', line.rstrip())
            if splitline[0][:4] == "Kala" and splitline[1][:4] != "Anco" \
                and splitline[1][:4] != "Kala" and splitline[1][:3] != "Aco":
                KafeArth.append(splitline)
            elif (splitline[0][:4] == "Anco" or splitline[0][:3] == "Aco") and splitline[1][:4] != "Anco" \
                and splitline[1][:4] != "Kala" and splitline[1][:3] != "Aco":
                AncoArth.append(splitline)
            elif splitline[0][:4] == "Kala" and (splitline[1][:4] == "Anco" or splitline[1][:3] == "Aco"):
                KafeAnco.append(splitline)

        for KAr in KafeArth:
            for AnAr in AncoArth:
                for KAn in KafeAnco:
                    if KAr[0] == KAn[0] and KAr[1] == AnAr[1] and AnAr[0] == KAn[1] \
                            and float(KAr[2]) >= threshold and float(KAn[2]) >= threshold and float(AnAr[2]) >= \
                            threshold:

                        ### kafe and anco shift together  ###

                        #  test if anco and kafe shifted together 10-14 hours from arth
                        if (10 <= int(KAr[4]) <= 14) and (10 <= int(AnAr[4]) <= 14) \
                                and (int(KAn[4]) <= 3 or int(KAn[4]) >= 22):
                            outfile1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                           (KAr[0], AnAr[0], KAr[1], KAr[4],
                                            KAr[2],KAr[3], AnAr[4],AnAr[2],AnAr[3],
                                            KAn[4],KAn[2],KAn[3]))
                            outfile1genes.write('%s\t%s\t%s\n' % (KAr[0],KAn[1],KAr[1]))

                        # test if anco and kafe shifted together 4-8 hours from arth
                        elif ((4 <= int(KAr[4]) <= 8) or (16 <= int(KAr[4]) <= 20)) \
                                and ((4 <= int(AnAr[4]) <= 8) or (16 <= int(AnAr[4]) <= 20)) \
                                and (int(KAn[4]) <= 3 or int(KAn[4]) >= 22):
                            outfile2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                           (KAr[0], AnAr[0], KAr[1], KAr[4],
                                            KAr[2],KAr[3], AnAr[4],AnAr[2],AnAr[3],
                                            KAn[4],KAn[2],KAn[3]))
                            outfile2genes.write('%s\t%s\t%s\n' % (KAr[0],KAn[1],KAr[1]))

                        ### only anco shifts ###

                        # test if just anco shifted 10-14 hours from arth and kafe
                        elif (10 <= int(AnAr[4]) <= 14) and (10 <= int(KAn[4]) <= 14) \
                                and (int(KAr[4]) <= 3 or int(KAr[4]) >= 22):
                            outfile3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                           (KAr[0], AnAr[0], KAr[1], KAr[4],
                                            KAr[2],KAr[3], AnAr[4],AnAr[2],AnAr[3],
                                            KAn[4],KAn[2],KAn[3]))
                            outfile3genes.write('%s\t%s\t%s\n' % (KAr[0],KAn[1],KAr[1]))

                        # test if just anco shifted 4-8 hours from arth and kafe
                        elif ((4 <= int(AnAr[4]) <= 8) or (16 <= int(AnAr[4]) <= 20)) \
                                and ((4 <= int(KAn[4]) <= 8) or (16 <= int(KAn[4]) <= 20)) \
                                and (int(KAr[4]) <= 3 or int(KAr[4]) >= 22):
                            outfile4.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                           (KAr[0], AnAr[0], KAr[1], KAr[4],
                                            KAr[2],KAr[3], AnAr[4],AnAr[2],AnAr[3],
                                            KAn[4],KAn[2],KAn[3]))
                            outfile4genes.write('%s\t%s\t%s\n' % (KAr[0],KAn[1],KAr[1]))

                        ### only kafe shifts ###

                        # test if just kafe shifted 10-14 hours from arth and anco
                        elif (10 <= int(KAr[4]) <= 14) and (10 <= int(KAn[4]) <= 14) \
                                and (int(AnAr[4]) <= 3 or int(AnAr[4]) >= 22):
                            outfile5.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                           (KAr[0], AnAr[0], KAr[1], KAr[4],
                                            KAr[2],KAr[3], AnAr[4],AnAr[2],AnAr[3],
                                            KAn[4],KAn[2],KAn[3]))
                            outfile5genes.write('%s\t%s\t%s\n' % (KAr[0],KAn[1],KAr[1]))

                        # test if just kafe shifted 4-8 hours from arth and anco
                        elif ((4 <= int(KAr[4]) <= 8) or (16 <= int(KAr[4]) <= 20)) \
                                and ((4 <= int(KAn[4]) <= 8) or (16 <= int(KAn[4]) <= 20)) \
                                and (int(AnAr[4]) <= 3 or int(AnAr[4]) >= 22):
                            outfile6.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                                           (KAr[0], AnAr[0], KAr[1], KAr[4],
                                            KAr[2],KAr[3], AnAr[4],AnAr[2],AnAr[3],
                                            KAn[4],KAn[2],KAn[3]))
                            outfile6genes.write('%s\t%s\t%s\n' % (KAr[0],KAn[1],KAr[1]))

    outfile1.close()
    outfile2.close()
    outfile3.close()
    outfile4.close()
    outfile5.close()
    outfile6.close()
    outfile1genes.close()
    outfile2genes.close()
    outfile3genes.close()
    outfile4genes.close()
    outfile5genes.close()
    outfile6genes.close()
