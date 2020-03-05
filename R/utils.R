
## Take an M x J matrix of detection probabilities and return a matrix
## of M x J observation probs
# Compute the cell probabilities for the observation classes
# in removal sampling.
#
# Both p and the returned matrix are M x J for M sites and J sampling occasions.

removalPiFun <- function(p){
  M <- nrow(p)
  J <- ncol(p)
  pi <- matrix(NA, M, J)
  pi[,1] <- p[,1]
  for(i in seq(from = 2, length = J - 1)) {
    pi[, i] <- pi[,i-1] / p[,i-1] * (1-p[,i-1]) * p[,i]
  }
  return(pi)
}

# p is an M x 2 matrix of detection probabilities (site x observer).
# returns an M x 3 matrix of row=(1 not 2, 2 not 1, 1 and 2).
# Compute the cell probabilities for the observation classes
# in double observer sampling.

doublePiFun <- function(p){
  M <- nrow(p)
  pi <- matrix(NA, M, 3)
  pi[,1] <- p[,1] * (1 - p[,2])
  pi[,2] <- p[,2] * (1 - p[,1])
  pi[,3] <- p[,1] * p[,2]
  return(pi)
}

# p is an M x 2 matrix of detection probabilities (site x observer).
# returns an M x 2 matrix of row=(1, 2 not 1).
# Compute the cell probabilities for the observation classes
# in double observer sampling.

depDoublePiFun <- function(p){
  M <- nrow(p)
  pi <- matrix(NA, M, 2)
  pi[,1] <- p[,1]
  pi[,2] <- p[,2]*(1-p[,1])
  return(pi)
}

genFixedNLL <- function(nll, whichFixed, fixedValues)
{
    function(params) {
        params[whichFixed] <- fixedValues
        do.call(nll, list(params))
        }
}

# nll the original negative log likelihood function
# MLE the full vector of MLE values
calc.profileCI <- function(nll, whichPar, MLE, interval, level)
{
    stopifnot(length(whichPar) == 1)
    MLEnll <- nll(MLE)
    nPar <- length(MLE)
	chsq <- qchisq(level, 1)/2
    f <- function(value) {
        fixedNLL <- genFixedNLL(nll, whichPar, value)
            mleRestricted <- optim(MLE, fixedNLL)$value
        mleRestricted - MLEnll - chsq
        }
    lower <- tryCatch(uniroot(f, c(interval[1],MLE[whichPar]))$root,
        error = function(e) -Inf)
    upper <- tryCatch(upper <- uniroot(f, c(MLE[whichPar], interval[2]))$root,
        error = function(e) Inf)

    return(c(lower,upper))
}

## link functions and their gradients
logistic <- function(x) {
  1/(1 + exp(-x))
}


logistic.grad <- function(x) {
  exp(-x)/(exp(-x)+1)^2
}


log.grad <- function(x) { # duh! (but for clarity)
  1/x
}


explink <- function(x) exp(x)

exp1 <- function(x) exp(x) + 1

identLink <- function(x) x

identLinkGrad <- function(x) 1

## use logarithms to vectorize row-wise products
## this speeds things up a LOT (vs. apply(x,1,prod))
rowProds <- function(x, na.rm = FALSE)
{
  exp(rowSums(log(x), na.rm = na.rm))
}


### track linked list of parameters using a data frame
### add row to linked list
addParm <- function(list.df, parm.name, parm.length) {
    if(parm.length > 0) {
        if(nrow(list.df) == 0) {
            last.ind <- 0
        } else {
            last.ind <- list.df$end[nrow(list.df)]
        }
        parm.df <- data.frame(parameter = parm.name, start = last.ind + 1,
                              end = last.ind + parm.length,
                              stringsAsFactors = FALSE)
        list.df <- rbind(list.df, parm.df)
    }
    return(list.df)
}


parmNames <- function(list.df) {
    npar <- list.df$end[nrow(list.df)]
    names <- character(npar)
    for(i in 1:npar) {
        which.par <- which(i >= list.df$start & i <= list.df$end)
        names[i] <- list.df$parameter[which.par]
    }
    return(names)
}


# get estimated psi from rn fit

getPsi <-
function(lam)
{
  1-exp(-lam)
}

# get estimatd p from rn fit (only for a null type model so far)

getP.bar <-
function(lam, r)
{
    K = 30
    psi <- getPsi(lam)
    pN.k <- dpois(0:K,lam)
    pY.k <- 1 - (1 - r)^(0:30)
    sum(pY.k * pN.k)
}


meanstate <- function(x) {
    K <- length(x) - 1
    sum(x*(0:K))
}

truncateToBinary <- function(y) {
    if(max(y, na.rm = TRUE) > 1) {
        y <- ifelse(y > 0, 1, 0)
        warning("Some observations were > 1.  These were truncated to 1.")
    }
    return(y)
}


getSS <- function(phi) {
    ev.length <- nrow(phi)
    ev <- tryCatch(eigen(t(phi))$vectors[,1],
                   error = function(x) rep(NA, ev.length))
    ev/sum(ev)
}


SSE <- function(fit)
{
    sse <- sum(residuals(fit)^2, na.rm=TRUE)
    return(c(SSE=sse))
}


# Generate zero-inflated Poisson

rzip <- function(n, lambda, psi) {
    x <- rpois(n, lambda)
    x[runif(n) < psi] <- 0
    x
}

nllFun<- function(object) object$nllFun

mle<- function(object) object$opt$par
