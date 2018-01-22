\name{topGene}
\alias{topGene}

\title{
List top ranked DE genes 
}
\description{
Given DE test results, provide a list of top ranked DE genes according
to user specified criteria
}
\usage{
topGene(de.obj, phase = "both", number = 20, p.value = 0.05, adjust.method = "BH", annotation = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{de.obj}{A data frame for DE test result, must be returned from twoPhaseDE.
  }
  \item{phase}{Character to indicate the phase for DE gene to be
	listed. Available options include '1' for phase I; '2' for phase II;
	and "both" for both phases. 
  }
  \item{number}{Integer for the maximum number of genes to be listed. If there are less than specified \code{number} of genes with adjusted p-value lower than \code{p.value}, all of these will be included in the table.
  }
  \item{p.value}{A cutoff in adjusted p.value. Only genes with adjusted \code{p.value} lower than this cutoff are included in the table.  
  }
  \item{adjust.method}{ The method used in the function \code{p.adjust\{stats\}}. The default is \code{"BH"} for Benjamini & Hochberg (1995). 
  }
  \item{annotation}{A vector of gene names for all genes in the \code{de.obj}. If left as \code{Null}, the 
  \code{featureNames} in the \code{de.obj} will be used.
  }
}

\details{
This function provides a summary table for the genes identified as differentially expressed. Only genes with an adjusted p-value lower than specified in \code{p.value} are included, and the maximum number is set by \code{number}. 
One may choose to extract genes with differential expression in Phase 1 only, Phase 2 only. If one is interested in  differential expression in any form, setting phase as "both" will include genes with either forms of differential expression, and the genes are sorted by significance in combined FDR. 
}

\value{
  A data frame with following columns:
  \item{Gene.name}{The gene identifier.}
  \item{p1}{Probability of Phase II in group 1.}
  \item{p2}{Probability of Phase II in group 2.}
  \item{Ph1.coef}{Test statistic for Phase change}
   \item{Ph1.fdr}{FDR for differential expression in the form of phase chage.}
\item{m1}{Estimated mean phase-2 expression in group 1}
\item{m2}{Estimated mean phase-2 expression in group 2}

  \item{Ph2.coef}{Test statistic for differential expression within Phase 2}
  \item{marLogFC}{Marginal log fold change}
  \item{Ph2.fdr}{FDR for differential expression in the form of magnitude difference in Phase 2}
  \item{Phase1}{Indicator for whether there is significant Phase change DE (based on the user specified \code{p.value} cutoff.)}
  \item{Phase2}{Indicator for whether there is significant Phase-2 DE (based on the user specified \code{p.value} cutoff.)}
  \item{Comb.fdr}{A combined FDR cnosidering both}
  \item{Category}{Indicates direction of DE. For example, "--"" down regulation in both the probability in Phase-2, and the magnitude of expression.  }
  
}

\author{
  Zhijin (Jean) Wu <zwu@stat.brown.edu>
}

\seealso{
  twoPhaseDE
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

##  DE test
de.sc2p <- twoPhaseDE(data, design="celltype", test.which=1, offset="sf")

## report top ranked genes in phase I
topGene(de.sc2p, 1)

## report top ranked genes in phase II
topGene(de.sc2p, 2)

## report top ranked genes in both phases
topGene(de.sc2p, phase="both", number=5)

}
}
