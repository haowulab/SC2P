visGene <- function(names, sc2p.obj, group.name=NULL, low.prob=.99,
                    show.name=TRUE, show.points=T, all.points=T,
                    annotation=NULL) {

    Y <- assayData(sc2p.obj)$exprs
    Z <- assayData(sc2p.obj)$Z
    Offset <- assayData(sc2p.obj)$Offset
    if (!is.null(annotation)){
        rownames(Y) <- annotation; rownames(Z) <- annotation
    }
    Z1 <- (Z > low.prob)
    Y.norm <- log2(Y+1) - Offset
    W = Y.norm; W[!Z1] <- NA
    group <- droplevels(pData(sc2p.obj)[, group.name])
    fit <- lmFit(W, model.matrix(~group-1))
    fit <- eBayes(fit)
    coef <- fit$coef
    delta <- qt(0.975, fit$df.residual + fit$df.prior)*fit$stdev.unscaled*sqrt(fit$s2.post)
    ci.lo <- coef - delta
    ci.hi <- coef + delta

    ind <- match(names, rownames(Y))
    if (all(is.na(ind))) { stop("The gene names are not found") }
    if (any(is.na(ind))) { warning(paste(sum(is.na(ind)),"genes are not plotted"))}
    ind <- na.omit(ind); names <- rownames(Y)[ind]

    sel1 <- as.numeric(group)==1
    sel2 <- as.numeric(group)==2
    n1=sum(sel1); n2=sum(sel2)
    at1=1;at2=2

    for (i in ind){
        thisGene <- rownames(Y)[i]
        y1 <- Y.norm[i, sel1]; y2 <- Y.norm[i, sel2]
        w1 <- W[i, sel1]; w2 <- W[i, sel2]
        m1 <- coef[i,1]; m2 <- coef[i,2]
        z1=Z[i, sel1]; z2=Z[i, sel2]
        n1.0=mean(z1<= low.prob); n2.0=mean(z2<=low.prob)
        l1 <- ci.lo[i,1]; u1 <- ci.hi[i,1]
        l2 <- ci.lo[i,2]; u2 <- ci.hi[i,2]

        if (show.name) title <- thisGene
        else title=""
        ## length of line
        ll <- 0.15

        plot(c(at1, at2), c(m1, m2), pch=16, col=rgb(0,0,0,0),
             xlim=c(0,3), ylim=c(0, max(c(y1,y2), na.rm=T)),
             xlab="",ylab="normalized log2 counts", axes=F,
             main=title)

        barplot(-n1.0*1,width=0.1,offset=1,horiz=T,add=T,
                col=rgb(.2,.2,.2,0.2),axes=F)
        barplot(n2.0*1,width=0.1,offset=2,horiz=T,add=T,axes=F,
                col=rgb(.2,.2,.2,0.2))
        text(at1,.5,round(n1.0,2))
        text(at2,.5,round(n2.0,2))
        axis(2);box()
        axis(1,at=c(at1,at2),label=levels(group)[1:2])

        if (show.points & !all.points){
            points(cbind(rnorm(n1, at1,.1), w1),
                   pch=16, col=rgb(0, 0, 1, z1))
            points(cbind(rnorm(n2,at2,.1), w2),
                   pch=16, col=rgb(0,0,1,z2))
        }
        if (show.points & all.points){
            points(cbind(rnorm(n1, at1,.1), y1),
                   pch=16, col=rgb(0, 0, 1, z1))
            points(cbind(rnorm(n2, at2,.1), y2),
                   pch=16, col=rgb(0,0,1,z2))
        }
        segments(at1-ll, m1, at1+ll, m1, col=2, lwd=2)
        segments(at1-ll/2, l1, at1+ll/2, l1, col=2, lwd=2)
        segments(at1-ll/2, u1, at1+ll/2, u1, col=2, lwd=2)
        segments(at1, l1, at1, u1, col=2)

        segments(at2-ll, m2, at2+ll, m2, col=2, lwd=2)
        segments(at2-ll/2, l2, at2+ll/2, l2, col=2, lwd=2)
        segments(at2-ll/2, u2, at2+ll/2, u2, col=2, lwd=2)
        segments(at2, l2, at2, u2, col=2)
    }
}
