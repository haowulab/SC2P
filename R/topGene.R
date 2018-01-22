## given number of row or fdr cutoff, list the genes based on phase
topGene <- function(de.obj, phase="both", number=20, p.value=0.05,
                    adjust.method="BH", annotation=NULL){
    phase <- as.character(phase)
    Gene.name <- rownames(de.obj)
    if (!is.null(annotation)){ Gene.name <- annotation }
    pval1 <- de.obj$Ph1.pval; pval2 <- de.obj$Ph2.pval
    tmp <- pmin(pval1, pval2, na.rm=T) ## folows beta(1,2) distribution
    pval <- pbeta(tmp, 1, 2) ## CDF follows uniform
    fdr1 <- p.adjust(pval1, method=adjust.method)
    fdr2 <- p.adjust(pval2, method=adjust.method)
    fdr <- p.adjust(pval, method=adjust.method)
    
    coef1 <- round(de.obj[, "Ph1.coef"],2)
    coef2 <- round(de.obj[, "Ph2.coef"],2)
    Phase1 <- ifelse(fdr1 < p.value, "Y", "N")
    Phase2 <- ifelse(fdr2 < p.value, "Y", "N")
    category <- ifelse(coef1>0 & coef2>0, "++",
                       ifelse(coef1>0 & coef2<0, "+-",
                              ifelse(coef1<0 & coef2>0, "-+", "--")))

    DF <- data.frame(Gene.name=Gene.name,
                     p1=de.obj[ ,"p1"], p2=de.obj[ ,"p2"],
                     Ph1.coef=coef1, Ph1.fdr=fdr1,  
                     m1=de.obj[ ,"m1"], m2=de.obj[ ,"m2"],
                     Ph2.coef=coef2, marLogFC  = de.obj[, "logFC"], Ph2.fdr=fdr2,
                     Phase1=Phase1,Phase2=Phase2,
                     Comb.fdr=fdr,
                     Category=category)

    ## select genes first
    sel <- switch(phase,
                  "1" = which(Phase1=="Y"),
                  "2" = which(Phase2=="Y"),
                  "both" = which(Phase1=="Y" | Phase2=="Y"))
    topn <- min(number, length(sel))
    DF <- DF[sel,] ## subset only the top ones
    
    ## then order the selected genes
    ord <- switch(phase,
                  "1"  = order(DF$Ph1.fdr, DF$Ph2.fdr),
                  "2" = order(DF$Ph2.fdr, DF$Ph1.fdr),
                  "both" = order(DF$Comb.fdr))
    glist <- DF[ord,]
    rownames(glist) <- NULL
    ## if (nrow(glist)==0) stop("No gene satisfies the criteria, try a larger FDR")
    ## round the numbers
    glist <- format(glist, digits=2)
    output <- glist[1:topn, , drop=F]
    return(output)   
}
