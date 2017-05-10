twoPhaseDE <- function(sc2p.obj,
                       design, ## vector of covariate names
                       test.which, ## which of design to be tested?
                       low.prob=.99,
                       offset=c("SC2P", "sf")){
    Y <- assayData(sc2p.obj)$exprs
    Z <- assayData(sc2p.obj)$Z
    Offset <- assayData(sc2p.obj)$Offset

    offset = match.arg(offset)
    if(offset == "SC2P") {
        Offset = assayData(sc2p.obj)$Offset
    } else {
        k = colSums(Y)
        k = log2(k/median(k))
        Offset = matrix(rep(k, nrow(Y)), nrow=nrow(Y), byrow=TRUE)
    }
    
    X <- pData(sc2p.obj)[, design, drop=F] ## or X can be permuted
    
    twoPhaseDE0(Y=Y, Z=Z, X=X, Offset=Offset,
                test.which=test.which, low.prob=low.prob)
}
