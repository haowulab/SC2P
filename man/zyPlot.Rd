\name{zyPlot}
\alias{zyPlot}

\title{
  Visualize the posterior probability of being expressed for a gene
}

\description{
  This functions plots the posterior probability of being expressed
  versus expressions level for a particular gene cross many cells. 
}

\usage{
zyPlot(name, sc2p.obj, group.name, log2scale = TRUE, showname = TRUE, annotation = NULL)
}

\arguments{
  \item{name}{Gene name. This should either match one of the
	\code{featureNames} of the \code{sc2p.obj}, or match a name in a
	separate \code{annotation} vector.} 
  \item{sc2p.obj}{An object of sc2pSet, result from the eset2Phase
	function.} 
  \item{group.name}{Character for the name of factors indicating the
	experimental conditions. This has to be the name of a column in
	pData(sc2p.obj) data frame.
  }
  \item{log2scale}{Logical to indicate whether to plot expression level in log2
	scale. Default is TRUE.
  }
  \item{showname}{Logical to indicate whether to show the name of genes
	as the tile of the figure.
  }
  \item{annotation}{A character vector of alternative gene names, with
	length equal the number of genes in the \code{sc2p.obj}. If left as
	default \code{NULL}, the \code{featureNames} of the \code{sc2p.obj}
	is used.}
}

\details{
In general, the probability of being expressed for a gene is higher when
the expression level is high. However, the probability is also affected
by other cell- and gene-specific factors. For example, cell with higher
total depth will have greater values for all genes belong to that
cell. Genes with higher average expression levels will have higher
background, thus the threshold for expression will be higher than genes
with lower average expression.

This function provides visualization of the expression status and
expression level. 
}

\author{
 Jamie (Yi) Zhang <yi_zhang@brown.edu>
}

\seealso{
  eset2Phase
}

\examples{
\dontrun{
data(brain_scRNAseq)

## assemble expressionset
colnames(Y) <- rownames(design)
phenoData <- new("AnnotatedDataFrame", data=design)
eset <- ExpressionSet(assayData=Y, phenoData=phenoData)

## infer latent status
data <- eset2Phase(eset)

## visualize a gene
zyPlot(rownames(data)[1], data, group.name="celltype")
}
}
