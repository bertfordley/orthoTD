def Ccorr(input_file, s, OutDirectory, InDirectory):
    '''

    Author: Rob Moseley

    Usage:
        - python Ccorr.py <input_file> <s> <OutDirectory>

    Input Arguments:

    :input_file         File containing orthologs and their expression
    :s                  Spline number
    :OutDirectory       Directory to save file


    Steps:
        - computes averages for each time point sample
        - converts data to 1-hour intervals
        - converts data to z-scores
        - runs cross correlation and spearman


    Outputs:
        - Table with spearman rank and time delay for gene pairs


    '''

    import os.path
    import pandas as pd
    import numpy as np
    from numpy.linalg import norm
    from scipy.interpolate import interp1d
    from scipy.stats import zscore
    from scipy.stats import spearmanr
    import ntpath

    def cxcorr(a, b):
        """
        For calculating the circular cross correlation between two vectors
        :param a: real or complex vector
        :param b: real or complex vector
        :return: circular correlation coefficients and respective lags
        """
        a /= norm(a)  # normalization
        b /= norm(b)  # normalization
        b = np.asarray(b)
        a = np.asarray(a)
        c = []
        for k in range(len(b)):
            cor = sum(a * b.conj().transpose())
            c.append(cor)
            # circular shift
            b = np.insert(b, 0, b[-1])
            b = b[:-1]
        a = range(len(b))  # lags
        return a, c

    fileName = os.path.splitext(ntpath.basename(input_file))[0]

    # get biorep data
    data = pd.read_csv(input_file, sep="\t")
    if len(data) < 2:
        print "No expression data present in orthogroup."
        return None
    else:
        data = data.set_index("Transcript")
        # copy last three reps to front
        data = pd.concat([data.iloc[:, -1:], data], axis=1)

        # average data
        # data = pd.concat([data.iloc[:, -3:], data], axis=1)
        # data = data.transpose()
        # data = data.groupby(np.arange(len(data)) // 3).mean()
        # data = data.transpose()
        # data.columns = ['Leaf_T06p1', 'Leaf_T08', 'Leaf_T10', 'Leaf_T12',
        #                 'Leaf_T14', 'Leaf_T16', 'Leaf_T18', 'Leaf_T20',
                        # 'Leaf_T22', 'Leaf_T00', 'Leaf_T02', 'Leaf_T04', 'Leaf_T06p2']
        # data.to_csv(fileName + "_AveragedData.txt", sep='\t')

        # cubic spline interpolation
        splineMatrix = []
        numbTP = len(data.columns)
        x = np.arange(numbTP) + 1
        for idx, row in data.iterrows():
            xx = np.linspace(1, numbTP, s)
            cs = interp1d(x, row, kind='cubic')(xx)
            splineMatrix.append(cs)

        splineMatrix = pd.DataFrame(splineMatrix)
        # remove repeated time sample (6am) and leave new 7am time sample
        splineMatrix.drop(splineMatrix.columns.values[0], axis=1, inplace=True)
        splineMatrix.columns = ['7am', '8am', '9am', '10am', '11am', '12am',
                                '1pm', '2pm', '3pm', '4pm', '5pm', '6pm',
                                '7pm', '8pm', '9pm', '10pm', '11pm', '12pm',
                                '1am', '2am', '3am', '4am', '5am', '6am']
        splineMatrix.index = data.index

        # convert to z-scores
        zscoreMatrix = splineMatrix.apply(zscore, axis=1)
        zscoreMatrix.to_csv(InDirectory + "_zscores/" + fileName + "_zScores.txt", sep='\t')
        # drop genes with zero expression
        zscoreMatrix.dropna(axis=0, how="all", inplace=True)

        # run algorithm on z-scores
        numbGenes = len(zscoreMatrix.index)
        geneNames = zscoreMatrix.index
        final_network = []
        final_timeDelays = []
        pvals = []

        if zscoreMatrix.shape[0] > 1:
            for i in range(numbGenes):
                Li = []
                maxCC = []
                for j in range(numbGenes):
                    lag, ccorr = cxcorr(zscoreMatrix.iloc[j, :], zscoreMatrix.iloc[i, :])
                    # get max absolute ccorr value and respective index within max time delay
                    maxCcorr = max(ccorr)
                    maxCcorrIdx = np.argmax(ccorr)
                    maxCC.append(maxCcorr)
                    # get lag value
                    TD = lag[maxCcorrIdx]
                    Li.append(TD)

                idxMemory = []
                geneIdx = 0
                for k in range(numbGenes):
                    if k == i:
                        geneIdx = k
                    idxMemory.append(k)

                Xi = np.zeros(shape=(len(idxMemory), zscoreMatrix.shape[1]))
                # get vector of time samples for target gene (i)
                X = np.asarray(zscoreMatrix[geneIdx:geneIdx + 1])
                Xi[[idxMemory.index(x) for x in idxMemory if x == geneIdx], :] = X
                # get vectors of time samples for potential regulators (j)
                for j in range(len(idxMemory)):
                    jIdx = idxMemory[j]
                    if jIdx != geneIdx:  # ignore target gene (i)
                        Lij = Li[jIdx]
                        xj = np.asarray(zscoreMatrix[jIdx:jIdx + 1])
                        xj = np.concatenate((xj[:, Lij:], xj[:, :Lij]), axis=1)
                        # extract time samples
                        Xi[[idxMemory.index(x) for x in idxMemory if x == jIdx], :] = xj  # align time samples
                # stack sequences row wise and convert to dataframe
                Xi = pd.DataFrame(np.vstack(Xi)).transpose()
                # transpose so genes are columns and create spearman correlaton matrix and p-val matrix (goes col by col)
                if zscoreMatrix.shape[0] > 2:
                    C, pval = spearmanr(Xi)
                    Ci = C[i,:]
                    pvali = pval[i,:]
                    # add to final network
                    final_network.append(np.asarray(Ci))
                    # add pvals
                    pvals.append(np.asarray(pvali))
                else:
                    C, pval = spearmanr(Xi[0], Xi[1])
                    # add to final network
                    final_network.append(C)
                    # add pvals
                    pvals.append(pval)
                final_timeDelays.append(np.asarray(Li))
                # adjust time delays
                # adjLi = []
                # for t in Li:
                #     if t > 12:
                #         newT = 12 -t
                #         adjLi.append(newT)
                #     else:
                #         adjLi.append(t)
                # final_timeDelays.append(np.asarray(adjLi))

        if zscoreMatrix.shape[0] > 1:
            outfile = open(OutDirectory + "/" + fileName + "_Table.txt", 'w')
            outfile.write("This table shows the Max Circular Correlation (CCorr) value, p-value, and respective Time Delay\n")
            outfile.write("The p-value roughly indicates the probability of an uncorrelated system producing datasets that "
                          "have a Spearman correlation at least as extreme as the one computed from these datasets. "
                          "The p-values are not entirely reliable but are probably reasonable for datasets larger than "
                          "500 or so.\n")
            outfile.write("Gene 1\tGene 2\tSpearman Coeff\tp-val\tTime Delay\n")
            if zscoreMatrix.shape[0] > 2:
                for i in range(len(final_network)):
                    for j in range(len(final_network)):
                        outfile.write('%s\t%s\t%.5f\t%.3f\t%s\n' % (geneNames[i], geneNames[j], final_network[i][j],
                                                                    pvals[i][j], final_timeDelays[i][j]))
            else:
                ccorrs = [[1,final_network[0]],[final_network[1],1]]
                pvalues = [[0,pvals[0]],[pvals[1],0]]
                for i in range(len(final_network)):
                    for j in range(len(final_network)):
                        outfile.write('%s\t%s\t%.5f\t%.3f\t%s\n' % (geneNames[i], geneNames[j], ccorrs[i][j],
                                                                    pvalues[i][j], final_timeDelays[i][j]))

