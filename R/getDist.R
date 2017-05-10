##################################################################
## given a sc2p object, calculate the distance matrix between cells
##################################################################

getLLR1 <- function(sc2p.obj, low.prob=.99){
    post.pi <- assayData(sc2p.obj)$Z
    
    n <- ncol(post.pi)
    D <- matrix(NA, n, n)
    p.g <- rowMeans(post.pi)
    keep.g <- which(p.g > 0 & p.g < 1)
    Z <- (post.pi[keep.g, ,drop=F] > low.prob)
    p.g <- p.g[keep.g]
    for (i in 1:n){
        for (j in i:n){
            z.tmp <- Z[,i] + Z[,j]
            D[i,j] = D[j,i] <- -2*sum((z.tmp==2)*log(p.g) + (z.tmp==0)*log(1-p.g))
        }
    }
    D
}
getLLR2 <- function(sc2p.obj){
    X <- assayData(sc2p.obj)$exprs
    post.pi <- assayData(sc2p.obj)$Z
    pdta <- pData(sc2p.obj)
    fdta <- fData(sc2p.obj)
    pois.lambda <- pdta$lambda; pois.p <- pdta$p0; L <- pdta$L
    mu <- fdta$mean; sigma <- fdta$sd

      G <- nrow(X)
    n <- ncol(X)
    p.g <- rowMeans(post.pi)
    pi0 <- cbind((1-p.g)^2, ## 00
                 p.g*(1-p.g),
                 p.g*(1-p.g), ## 01, 10
                 p.g^2) ## 11
    cc0 <- pi0[, 1] + pi0[ ,4]

    den.zip = den.lnp = NA*X
    for (i in 1:n){
        den.zip[,i] <- dZinf.pois(X[,i], pois.p[i], pois.lambda[i])
        den.lnp[,i] <- dLNP2(X[,i,drop=F], mu, sigma, L[i])
    }
    LLR <- matrix(NA, n, n)
    for (i in 1:n){
        for (j in i:n){
            ## under the null
            den00_0 <- den.zip[,i] * den.zip[,j]
            den01_0 <- den.zip[,i] * den.lnp[,j]
            den10_0 <- den.lnp[,i] * den.zip[,j]
            den11_0 <- den.lnp[,i] * den.lnp[,j]
            den0 <- cbind(den00_0, den01_0, den10_0, den11_0) * pi0
            L0 <- rowSums(den0)
            
            mu.ave <- log2(rowSums(X[, c(i,j)])/sum(L[c(i,j)])+1) ## 0 when both 0
            mu.ave[is.na(mu.ave)] <- 0 # when library size is 0?
            dlnp.i <- dLNP2(X[,i], mu.ave, sigma, L[i])
            dlnp.j <- dLNP2(X[,j], mu.ave, sigma, L[j])
            
            den00_1 <- den.zip[,i] * den.zip[,j]
            den01_1 <- den.zip[,i] * dlnp.j
            den10_1 <- dlnp.i * den.zip[,j]
            den11_1 <- dlnp.i * dlnp.j

            pi1  <- cbind((1-post.pi[,i])*(1-post.pi[,j]),
                          (1-post.pi[,i])*post.pi[,j],
                          post.pi[,i]*(1-post.pi[,j]),
                          post.pi[,i]*post.pi[,j])
            cc1 <- pi1[,1] + pi1[,4]
            C <- (cc1 > cc0)
            pi1 <- pi1*C + pi0*(1-C)
            den1 <- cbind(den00_1, den01_1, den10_1, den11_1) * pi1
            L1 <- rowSums(den1)
            
            llr <- log(L1) - log(L0)
            llr[llr < 0] <- 0
            ## extreme cases are excluded
            ## llr[is.infinite(llr)] <- NA
            llr[is.infinite(llr)] <- max(llr[is.finite(llr)], na.rm=T)
            LLR[i,j] = LLR[j,i] <- sum(llr, na.rm=T) 
        }
        print(paste0("Calculating ", i, "th sample"))
    }
    LLR
}

makeDist <- function(LLR){
    di <- as.matrix(expand.grid(diag(LLR), diag(LLR)))
    d <- (LLR - di[,1]/4 -di[,2]/4)/sqrt(di[,1]/2*di[,2]/2)
    d <- pmin(d, 1); d <- pmax(d, -1)
    acos(d)/pi;
}
