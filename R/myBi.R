
#### modified family definition to handle logistic regression
#### with less than 100% sensitivity.

zwulink<-function(){
       linkfun <- function(mu) .Call(C_logit_link, mu/my.offset)
       linkinv <- function(eta) .Call(C_logit_linkinv, eta)*my.offset
       mu.eta <- function(eta) .Call(C_logit_mu_eta, eta)*my.offset
       valideta <- function(eta) TRUE
  environment(linkfun) <- environment(linkinv) <- environment(mu.eta)   <- environment(valideta) <- asNamespace("stats")   
     structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
        valideta = valideta, name = "mylogit"), class = "link-glm")
   }
   

mybi<-function (link = "mylogit") ##with offsets my.offset
{
   # linktemp <- substitute(link)
   # if (!is.character(linktemp))
   #     linktemp <- deparse(linktemp)
   # okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
   # if (linktemp %in% okLinks)
   #     stats <- make.link(linktemp)
   # else if (is.character(link)) {
   #     stats <- make.link(link)
   #     linktemp <- link
   # }
   # else {
   #     if (inherits(link, "link-glm")) {
   #         stats <- link
   #         if (!is.null(stats$name))
   #             linktemp <- stats$name
   #     }
   #    else {
   #         stop(gettextf("link \"%s\" not available for binomial family; available links are %s",
   #             linktemp, paste(sQuote(okLinks), collapse = ", ")),
   #             domain = NA)
   #     }
   # }
   linktemp<-"mylogit"
   stats <- zwulink()
    variance <- function(mu) mu * (1 - mu)
    validmu <- function(mu) all(is.finite(mu)) && all(mu > 0 &
        mu < 1)
    dev.resids <- function(y, mu, wt) .Call(C_binomial_dev_resids,
        y, mu, wt)
        #dev.resids <- function(y, mu, wt) .Call(C_binomial_dev_resids,
        #y, mu*my.offset, wt)
    aic <- function(y, n, mu, wt, dev) {
        m <- if (any(n > 1))
            n
        else wt
        -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m *
            y), round(m), mu, log = TRUE))
    }
    initialize <- expression({
        if (NCOL(y) == 1) {##bernoulli
            if (is.factor(y)) y <- y != levels(y)[1L]
            n <- rep.int(1, nobs)
            y[weights == 0] <- 0
            if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
            mustart <- (weights * y + 0.5)/(weights + 1)*my.offset
            m <- weights * y
            if (any(abs(m - round(m)) > 0.001)) warning("non-integer #successes in a binomial 
            glm!")
        } else if (NCOL(y) == 2) {
            if (any(abs(y - round(y)) > 0.001)) warning("non-integer counts in a binomial glm!")
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
            mustart <- (n * y + 0.5)/(n + 1)*my.offset
        } else stop("for the 'binomial' family, y must be a vector of 0 and 1's\nor a 2 column 
        matrix where col 1 is no. successes and col 2 is no. failures")
    })
    simfun <- function(object, nsim) {
        ftd <- fitted(object)
        n <- length(ftd)
        ntot <- n * nsim
        wts <- object$prior.weights
        if (any(wts%%1 != 0))
            stop("cannot simulate from non-integer prior.weights")
        if (!is.null(m <- object$model)) {
            y <- model.response(m)
            if (is.factor(y)) {
                yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd),
                  labels = levels(y))
                split(yy, rep(seq_len(nsim), each = n))
            }
            else if (is.matrix(y) && ncol(y) == 2) {
                yy <- vector("list", nsim)
                for (i in seq_len(nsim)) {
                  Y <- rbinom(n, size = wts, prob = ftd)
                  YY <- cbind(Y, wts - Y)
                  colnames(YY) <- colnames(y)
                  yy[[i]] <- YY
                }
                yy
            }
            else rbinom(ntot, size = wts, prob = ftd)/wts
        }
        else rbinom(ntot, size = wts, prob = ftd)/wts
    }
    structure(list(family = "binomial", link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
        validmu = validmu, valideta = stats$valideta, simulate = simfun),
        class = "family")
}
environment(mybi) <- environment(glm)
environment(mybi) <- environment(glm)
