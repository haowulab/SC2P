################################################
## demo script using part of public brain data
################################################

library(SC2P)
data(brain_scRNAseq)

## assemble expressionset
colnames(Y) <- rownames(design)
phenoData <- new("AnnotatedDataFrame", data=design)
eset <- ExpressionSet(assayData=Y, phenoData=phenoData)
eset

#### infer latent status
data <- eset2Phase(eset)

## visualize the estimated latent state
zyPlot(rownames(data)[1], data, group.name="celltype")

#### DE test
de.sc2p <- twoPhaseDE(data, design="celltype", test.which=1, offset="sf")

## report top ranked genes in phase I
topGene(de.sc2p, 1)

## report top ranked genes in phase II
topGene(de.sc2p, 2)

## top ranked genes in both phases
topGene(de.sc2p, phase="both")


###### visualize some genes
## top ranked phase I DE gene
visGene(topGene(de.sc2p, 1)$Gene.name[1], data, group.name="celltype")
## top ranked phase II DE gene
visGene(topGene(de.sc2p, 2)$Gene.name[1], data, group.name="celltype")

