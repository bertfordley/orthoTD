'''
adds expression data to orthogroup files
creates new folder and files
'''


import sys, re
import pandas as pd

# assigned genes
Afile = sys.argv[1]
orthoData = pd.read_csv(Afile, sep='\s+', header=None, engine='python')
orthoData  = orthoData.fillna(value="skip", axis=0)

# unassigned genes
UAfile = sys.argv[2]
UAdata = pd.read_csv(UAfile,  sep='\s+', header=None, engine='python')

# expression data
araData = pd.read_csv("/Users/rkd/Desktop/OrthoTD/TD+/orthoFinder/expression/\
Diurnal_Expression_Arabidopsis.csv")
araData = araData.set_index('Loci')

kalaData = pd.read_csv("/Users/rkd/Desktop/OrthoTD/TD+/orthoFinder/expression/\
Kaladp_expression_FPKM_primary.csv")
kalaData = kalaData.set_index('Transcript')

pineData = pd.read_csv("/Users/rkd/Desktop/OrthoTD/TD+/orthoFinder/expression/\
Diurnal_Expression_Pineapple.csv")
pineData = pineData.set_index('Transcripts')

dfs = [araData, kalaData, pineData]
expData = pd.concat(dfs)

rx = re.compile("(?:\.|^)[^\.]*")

# assigned genes
for Aidx, row in orthoData.iterrows():
    orthoGroup = row[0][:-1]
    print "Adding expression data to orthogroup: " + orthoGroup
    OUT = open("orthoExp/" + orthoGroup + ".txt", "w")
    for gene in row[1:]:
        if gene == "skip":
            pass
        else:
            regGene = rx.findall(gene)
            for expIndex, expRow in expData.iterrows():
                if regGene[0] in expIndex:
                    geneExp = expRow.tolist()
                    geneExp = [regGene[0]] + geneExp
                    OUT.write("\t".join([str(x) for x in geneExp]) + "\n")
                    expData = expData.drop(expIndex)
                    break
    OUT.close()

# unassigned genes
for UAidx, row in UAdata[1:].iterrows():
    orthoGroup = row[0][1:]
    print "Adding expression data to orthogroup: " + orthoGroup
    gene = row[1]
    geneExp = rx.findall(gene)
    for expIndex, expRow in expData.iterrows():
        if geneExp[0] in expIndex:
            OUT = open("orthoExp/" + orthoGroup + ".txt", "w")
            expression = expRow.tolist()
            expression = [geneExp[0]] + expression
            OUT.write("\t".join([str(x) for x in expression]) + "\n")
            expData = expData.drop(expIndex)
            OUT.close()
            break

