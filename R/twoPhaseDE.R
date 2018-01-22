### two phase DE test.
### should allow user to provide a normalization factor matrix/vector
twoPhaseDE <- function(sc2p.obj,
                       design, ## vector of covariate names
                       test.which, ## which of design to be tested?
                       low.prob=.99,
                       offset=c("SC2P", "sf")){
    Y <- assayData(sc2p.obj)$exprs
    Z <- assayData(sc2p.obj)$Z

    ## determine offset
    if(is.character(offset)) { ## offset is character
        offset <- match.arg(offset)
        if(offset == "sf") { # use size factor (derived from total counts)
            k = colSums(Y)
            k = log2(k/median(k))
            Offset = matrix(rep(k, nrow(Y)), nrow=nrow(Y), byrow=TRUE)
        } else ## use SC2P offset
        Offset <- assayData(sc2p.obj)$Offset
    } else if (is.vector(offset)) { ## offset is a vector
        if(length(offset) != ncol(Y))
            stop("offset has wrong length. Must equal to the number of cells.")
        ## turn it into a matrix, replicate for each gene
        Offset = matrix(rep(offset, nrow(Y)), nrow=nrow(Y), byrow=TRUE)
    } else if (is.matrix(offset)) { ## offset is a matrix.
        if( any(dim(offset) != dim(Y)) )
            stop("offset has wrong dimension. Must have the same dimension as the input expression.")
        Offset = offset
    } else stop("Unrecognized option for offset.")


    X <- pData(sc2p.obj)[, design, drop=F] ## or X can be permuted

    twoPhaseDE0(Y=Y, Z=Z, X=X, Offset=Offset,
                test.which=test.which, low.prob=low.prob)
}
