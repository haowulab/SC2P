\name{eset2Phase}
\alias{eset2Phase}

\title{
  Estimate the phases of expression 
}

\description{
  This function estimates the probabilities of expression for all genes in
  all cells. Given an ExpressionSet, it returns the
  posterior probability for being in phase II (expression) for each gene
  in each cell.  
}

\usage{
eset2Phase(eset, low.prob = 0.99)
}

\arguments{
  \item{eset}{An object of ExpressionSet for the single cell RNA-seq
  data. This should be constructed using the ExpressionSet function
  provided by Biobase.}

  \item{low.prob}{The probability threshold for calling phase II. When
  the posterior probability for a gene in a cell is greater than this
  threshold, it will be called as in phase II (expressed).}

}

\details{
  The estimation procedure considers the cell-specific (for example,
  sequencing depth and marginal distribution of expression for all
  genes), and gene-specific (for example, average expression 
  level for this gene cross all cells) factors, which achieves a
  data-driven thresholding and provides better results than using an
  fixed, ad hoc threshold for determining expression state. 
}

\value{
  An object of sc2pSet, which is an extension of ExpressionSet.
}


\author{
Zhijin (Jean) Wu <zwu@stat.brown.edu>
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

## visualize the results for a gene
zyPlot(rownames(data)[1], data, group.name="celltype")
}
}

\seealso{
twoPhaseDE
}

