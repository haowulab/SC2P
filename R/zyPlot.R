## plot a gene's estimated posterior Z=1 prob against Y
zyPlot <- function(name, sc2p.obj, group.name,
                   log2scale=TRUE, showname=TRUE, annotation=NULL){
    if (showname) title <- name
    else title <- ""
    
    Y <- assayData(sc2p.obj)$exprs
    Z <- assayData(sc2p.obj)$Z
    group <- pData(sc2p.obj)[, group.name]
    if (is.factor(group)) group <- droplevels(group)

    if (!is.null(annotation)){
        ind <- match(name, annotation)
    }else{
        ind <- match(name, rownames(Y))
    }
    if (is.na(ind)) stop("Gene not found")
    pchs <- c(1, 16)
    names(pchs) <- levels(group)
    if (log2scale){
        plot(log2(Y[ind, ]+1), Z[ind, ], col=group, pch=pchs[group],
             xlab="log2(expression count)",
             ylab="Posterior prob of foreground",
             main=title, ylim=c(0,1))
    }
    else {
        plot(Y[ind,], Z[ind, ], col=group, pch=pchs[group],
             xlab="log2(expression count)",
             ylab="posterior prob of foreground",
             main=title, ylim=c(0,1))
    }
}

    
