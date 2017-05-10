####################################
## rotated image to view matrix as is
image2=function(x) image(t(x)[,nrow(x):1])

rowVar <- function(x){
    x= sweep(x,1,rowMeans(x))
    x=rowSums(x^2)/(ncol(x)-1)
}   

dZinf.pois=function(x, p0, mu){
    (x==0)*(p0+(1-p0)*dpois(x,mu))+(x>0)*(1-p0)*dpois(x,mu)}

dLNP2 <- function(x, mu, sigma, l=1){
    x.min <- pmax(0, x-0.5)
    pnorm(log2((x+0.5)/l), mu, sigma) - pnorm(log2(x.min/l), mu, sigma)
}

### MAD using only the higher half
my.mad=function (x, center = median(x), constant = 1.4826, na.rm = FALSE)
{
    if (na.rm) 
        x <- x[!is.na(x)]
    res=x-center
    constant * median(res[res>0]) 
}

std.gene=function(X,L,mu.g,sd.g){ 
    ## L is size factor;mu.g,
    ## sd.g are gene-wise mean and sd in log scale
    X=sweep(X,2,L/1e6,FUN="/")
    X=log2(1+X)
    X=sweep(X,1,mu.g)
    X=sweep(X,1,sd.g,FUN="/")
    X}


shrink.mu=function(y,s,n){
                                        #y: gene-specific mean
                                        #s: gene -specific sd
                                        #n: gene -wise sample size
    mu.g=rep(NA,length(y))
    k=which(n>1)
    if (length(k)<length(n)) {fill=TRUE} else {fill=FALSE}
    s=s[k];y=y[k];n=n[k]

    mu0=weighted.mean(y,w=n) 
    
    s2.total=sum(s^2*(n-1))+sum(n*(y-mu0)^2)
    s2.total=s2.total/sum(n)
    
    s2.1=sum(s^2*(n-1))/sum(n)
    s2.0=s2.total-s2.1
                                        #paramters
                                        #c(mu0=mu0,s2.0=s2.0,s2.1=s2.1)

### shrink mu
    mu.sub=  (y*n/s2.1+mu0/s2.0)/(n/s2.1+1/s2.0)
    mu.g[k]=mu.sub
    if (fill) mu.g[-k]=mu0
    
    mu.g
} 


#####################
cond.var=function(Y,Z,min.p=0.99){ 
                                        # Y is matrix of values; Z is matrix of indicators
    Z=Z>min.p
    n=rowSums(Z)
    m=rowSums(Y*Z)/n #conditional Mean
    v=rowSums(sweep(Y,1,m)^2*Z)/(n-1)  
    v[n==0]=0#if n=0, negative, no var
    v[n==1]=NA# if n-1, var NA 
    cbind("cond.var"=v,"df"=n-1)
}
#####################

RobustPoi0 <- function(x){
    tt=table(x)
    n0=sum(x==0)
    if (names(tt)[1]=="0") {
        xx=as.numeric(names(tt))[-1]
        tt=as.vector(tt)[-1]
    }    
    tt=log(tt)+log(gamma(1+xx))
    ##fit the regression without N0
    beta=lm(tt~xx,weight=1/exp(xx))$coef
    mu.hat=exp(beta[2])
    p0=(n0-exp(beta[1]))/(exp(beta[1]+mu.hat)+n0-exp(beta[1]))
    if (any(p0<0)) {warning("Proportion of zero inflation is negative")}
    p0 <- pmax(p0, 0) 
    
    c("p0"=p0,"mu"=mu.hat)
}


