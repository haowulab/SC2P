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
  \item{number}{Integer for number of genes to be listed. 
  }
  \item{p.value}{????????????
  }
  \item{adjust.method}{
	??????????????????
  }
  \item{annotation}{
	????????????????????
  }
}

\details{
  HOW TO RANK GENES WHEN 'BOTH' IS THE OPTION????
}

\value{
  A data frame with following columns:
  
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
