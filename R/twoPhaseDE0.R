##################################################
twoPhaseDE0 <- function(Y, Z, X, Offset, test.which, 
                        low.prob=.99){
    vars <- colnames(X);  vars0 <- vars[-test.which]
    group <- X[, test.which]
    if (!is.factor(group)){ stop("The variable to be tested must be a factor") }
    group <- droplevels(group)
    X[, test.which] <- group ## just putting it back for Phase II
    Ng <- length(levels(group))
    parse1 <- parse(text= paste0("glm(yyy~", paste(vars, collapse="+"),
                        ",data=X, family=binomial)"))
    contrast <- paste0(vars[test.which], levels(group)[Ng])
    ## ##############################################
    ## phase 1: change of on rate
    Z1=( Z > low.prob)^2
    ## avgZ=tapply(1:ncol(Z), group, function(ind){
    ##     rowMeans(Z[,ind]) })  
    ## avgZ=matrix(unlist(avgZ),ncol=Ng)
    n.on=rowSums(Z1)
    ind=which(n.on > 0 & n.on < ncol(Y)) 
    if (length(vars)==1){ ## single binary variable
        DE.z <- matrix(NA, nrow=nrow(Z), ncol=4)
        rownames(DE.z) <- rownames(Y)
        DE.z[ind, ]= t(apply(Z1[ind,], 1,function(yyy){
            fit <- eval(parse1)
            ss <- summary(fit)
            coef <- ss$coef[contrast, "Estimate"]
            pval <- pchisq(ss$null.deviance - ss$deviance,
                           df=ss$df.null - ss$df.residual,
                           lower.tail=FALSE)
            c(mean(yyy[group==levels(group)[1]]),
              mean(yyy[group==levels(group)[2]]),
              coef, pval)
        }))
        colnames(DE.z) <- c("p1", "p2", "Ph1.coef", "Ph1.pval")
    } else {
        DE.z <- matrix(NA, nrow=nrow(Z), ncol=2)
        rownames(DE.z) <- rownames(Y)
        parse0 <- parse(text= paste0("glm(yyy~", paste(vars0, collapse="+"),
                            ",data=X, family=binomial)"))
        DE.z[ind, ]= t(apply(Z1[ind,], 1,function(yyy){
            fit1 <- eval(parse1); fit0 <- eval(parse0)
            ss1 <- summary(fit1); ss0 <- summary(fit0)
            coef <- ss1$coef[contrast, "Estimate"]
            pval <- pchisq(ss0$deviance - ss1$deviance,
                           df=ss0$df.residual - ss1$df.residual,
                           lower.tail=FALSE)
            c(coef, pval)
        }))
        colnames(DE.z) <- c("Ph1.coef", "Ph1.pval")
    }
    ## ################################################
    ## phase 2: conditional FC
    W <- log2(Y+1) - Offset; W[!Z1] <- NA
    modelX <- eval(parse(text=paste0("model.matrix(~", paste(vars, collapse="+"),
                             ", data=X)")))
    fit <- lmFit(W, modelX)
    fit <- eBayes(fit)
    mu1 <- fit$coef[,"(Intercept)"]
    mu2 <- rowSums(fit$coef)
    coef <- fit$coef[, contrast]
    stdev <- fit$stdev.unscaled[, contrast]
    delta <- qt(.975, fit$df.residual + fit$df.prior)*stdev*sqrt(fit$s2.post)
    ci.lo <- coef - delta; ci.hi <- coef + delta
    DE.y <- cbind(mu1, mu2, coef, ci.lo, ci.hi, fit$p.value[, contrast])
    colnames(DE.y) <- c("m1", "m2", "Ph2.coef",
                        "Ph2.ci.lo", "Ph2.ci.hi", "Ph2.pval")

    ## ################################################
    ## marginal logFC change
    sel1 <- as.numeric(group)==1
    sel2 <- as.numeric(group)==2
    logFC <- apply(log2(Y+1), 1, function(y){
        mean(y[sel2], na.rm=T) - mean(y[sel1], na.rm=T)
    })
    as.data.frame(cbind(DE.z, DE.y, logFC))
}
