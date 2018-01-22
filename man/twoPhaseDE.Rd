\name{twoPhaseDE}

\alias{twoPhaseDE}

\title{
  Two-phase differential expression test for scRNA-seq data
}
\description{
  This function performs differential expression analysis for
  singel-cell RNA-seq data in two phases: phase transition (probability
  of being expressed), and magnitude tuning	(levels of expression once
  the gene is on). Currently the function only works for the comparison
  of two biological conditions. 
													
}

\usage{
twoPhaseDE(sc2p.obj, design, test.which, low.prob = 0.99, offset = c("SC2P", "sf"))
}

\arguments{
  \item{sc2p.obj}{An object of sc2pSet. Should be the result object from
  'eset2Phase' function.}

  \item{design}{A character vector of variable names in pData(sc2p.obj) to
	provide regression covariates in DE test.}

  \item{test.which}{An integer for the index of the
	to-be-tested binary variable in argument 'design'.}
  \item{low.prob}{The cutoff used in classifying genes to Phase 2. Genes with posterior probability in Phase 2 greater than \code{low.prob} is considered in Phase 2.  
  }
  \item{offset}{
	The normalization factor, or indicator for method to compute
	normalization factor. Valid entry includes:
	(1) a character, where 'SC2P' indicate to use the	normalization
	procedure provided by SC2P. 'sf' uses the total counts. (2) a
	numeric vector, where each represent the normalization factor for a
	cell. The length of the vector must equal to the number of
	cells. (3) a numeric matrix of the same dimension as the input
	expression data, where each gene in each cell has a seperate
	normalization factor.}
}

\value{
  A data frame with following columns:
  \item{p1}{The Phase 2 proportion in Group 1. }
  \item{p2}{The Phase 2 proportion in Group 2.}
  \item{Ph1.coef}{Log Odds Ratio in testing DE in the form of phase change.}
  \item{Ph1.pval}{The p-value in testing DE in the form of phase change.}
  \item{m1}{Estiamted conditional mean expression in Group 1}
  \item{m2}{Estiamted conditional mean expression in Group 2}
  \item{Ph2.coef}{Conditional log2 fold change  }
  \item{Ph2.ci.lo}{The lower end of 95\% confidence interval for the conditional log fold change. }
  \item{Ph2.ci.hi}{The upper end of 95\% confidence interval for the conditional log fold change.}
  \item{Ph2.pval}{The p-value in testing DE in the form of magnitude change within Phase 2. }
  \item{logFC}{Marginal log2 fold change}

  The data frame is sorted by gene names alphabetically. 
}

\details{
  est2Phase needs to be called before this function to estimate the
  probabilities of being in two states. 
}

\author{
  Zhijin (Jean) Wu <zwu@stat.brown.edu>
}


\seealso{
  est2Phase, topGene, visGene
}

\examples{
\dontrun{
data(brain_scRNAseq)

## create ExpressionSet
colnames(Y) <- rownames(design)
phenoData <- new("AnnotatedDataFrame", data=design)
eset <- ExpressionSet(assayData=Y, phenoData=phenoData)

## estimate phases
data <- eset2Phase(eset)
data

##  DE test
de.sc2p <- twoPhaseDE(data, design="celltype", test.which=1, offset="sf")

## list top ranked DE genes
topGene(de.sc2p, phase=1)
topGene(de.sc2p, phase=2)
topGene(de.sc2p, phase="both")

## visualize top DE gene in phase 2
visGene(topGene(de.sc2p, 2)$Gene.name[1], data, group.name="celltype")
}
}

