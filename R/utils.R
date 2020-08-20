
## Take an M x J matrix of detection probabilities and return a matrix
## of M x J observation probs
# Compute the cell probabilities for the observation classes
# in removal sampling.
#
# Both p and the returned matrix are M x J for M sites and J sampling occasions.

removalPiFun <- function(p){
  if(is.matrix(p)) {
    M <- nrow(p)
    J <- ncol(p)
    pi <- matrix(NA, M, J)
    pi[,1] <- p[,1]
    for(i in seq(from = 2, length = J - 1)) {
      pi[, i] <- p[,i] * apply(as.matrix(p[,1:(i-1)]),1,function(x){prod(1-x)})
    }
  }
  else {
    J<- length(p)
    pi<- rep(NA, J)
    pi[1]<- p[1]
    for (i in 2:J)
      pi[i] <- p[i] * prod(1-p[1:(i-1)])
  }
  return(pi)
}

# removalPiFun <- function(p){
#   if(is.matrix(p)) {
#     M <- nrow(p)
#     J <- ncol(p)
#     pi <- matrix(NA, M, J)
#     pi[, 1] <- p[, 1]
#     for (i in seq(from = 2, length = J - 1)) {
#       pi[, i] <- pi[, i - 1]/p[, i - 1] * (1 - p[, i - 1]) * p[, i]
#     }
#   }
#   else {
#     J<- length(p)
#     pi<- rep(NA, J)
#     pi[1]<- p[1]
#     for (i in 2:J)
#       pi[i] <- p[i] * prod(1-p[1:(i-1)])
#   }
#   return(pi)
# }

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

#Complimentary log log link
cloglog <- function(x){
  1-exp(-exp(x))
}

cloglog.grad <- function(x){
  exp(-exp(x))
}

log.grad <- function(x) {
  1/x
}

exp.grad <- function(x) exp(x)

exp1 <- function(x) exp(x) + 1

identLink <- function(x) x

identLink.grad <- function(x) 1

#Complimentary log log link
cloglog <- function(x){
  1-exp(-exp(x))
}

cloglog.grad <- function(x){
  exp(-exp(x))
}

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

#---------------------------------
formatMult <- function(df.in) {
# This convenience function converts multi-year data in long format to
# eFrameRGP Object.

  years <- sort(unique(df.in[[1]]))
  nY <- length(years)
  df.obs <- list()
  nsamp <- numeric()
  maxsamp <- max(table(df.in[[1]], df.in[[2]])) # the maximum samples/yr
  for(t in 1:nY){
    df.t <- df.in[df.in[[1]] == years[t],] # subset for current year
    df.t <- df.t[,-1] # remove year column
    df.t <- dateToObs(df.t)
    nsamp <- max(df.t$obsNum)
    if(nsamp < maxsamp) {
      newrows <- df.t[1:(maxsamp - nsamp), ] # just a placeholder
      newrows[,"obsNum"] <- ((nsamp + 1) : maxsamp)
      newrows[,3 : (ncol(df.t) - 1)] <- NA
      df.t <- rbind(df.t, newrows)
    }
    df.obs <- rbind(df.obs,cbind(year = years[t],df.t))
  }

  names(df.obs)[4] <- "y"
  scol <- names(df.obs)[2]
  df_sub <- df.obs[c("year", scol, "y", "obsNum")]
  df_sub$yind <- 1:nrow(df_sub)

  ywide <- reshape(df_sub, idvar=c(scol,"year"), timevar="obsNum",
                   direction="wide")
  ywide <- reshape(ywide, idvar=scol, timevar="year", direction="wide")
  #Reshape goofs up the order
  ywide <- ywide[order(ywide[,1]),]

  y <- unname(as.matrix(ywide[grep("^y\\.", names(ywide))]))
  yind <- as.vector(t(as.matrix(ywide[grep("yind", names(ywide))])))

  #Reorder input data frame by y-index (=no factor issues)
  #Also drop site/time/y/numObs cols
  obsvars.df <- df.obs[yind, -c(1:2, 4, ncol(df.obs)), drop=FALSE]
  rownames(obsvars.df) <- NULL

  ## check for siteCovs
  obsNum <- ncol(y)
  M <- nrow(y)
  site.inds <- matrix(1:(M*obsNum), M, obsNum, byrow = TRUE)
  siteCovs <- sapply(obsvars.df, function(x) {
    obsmat <- matrix(x, M, obsNum, byrow = TRUE)
    l.u <- apply(obsmat, 1, function(y) {
      row.u <- unique(y)
      length(row.u[!is.na(row.u)])
    })
    ## if there are 0 or 1 unique vals per row, we have a sitecov
    if(all(l.u %in% 0:1)) {
      u <- apply(obsmat, 1, function(y) {
        row.u <- unique(y)
        ## only remove NAs if there are some non-NAs.
        if(!all(is.na(row.u)))
          row.u <- row.u[!is.na(row.u)]
        row.u
      })
      u
    }
  })
  siteCovs <- as.data.frame(siteCovs[!sapply(siteCovs, is.null)])
  if(nrow(siteCovs) == 0) siteCovs <- NULL

  ## only check non-sitecovs
  obsvars.df2 <- as.data.frame(obsvars.df[, !(names(obsvars.df) %in%
                                                names(siteCovs))])
  names(obsvars.df2) <- names(obsvars.df)[!(names(obsvars.df) %in%
                                              names(siteCovs))]

  yearlySiteCovs <- sapply(obsvars.df2, function(x) {
    obsmat <- matrix(x, M*nY, obsNum/nY, byrow = TRUE)
    l.u <- apply(obsmat, 1, function(y) {
      row.u <- unique(y)
      length(row.u[!is.na(row.u)])
    })
    ## if there are 0 or 1 unique vals per row, we have a sitecov
    if(all(l.u %in% 0:1)) {
      u <- apply(obsmat, 1, function(y) {
        row.u <- unique(y)
        ## only remove NAs if there are some non-NAs.
        if(!all(is.na(row.u)))
          row.u <- row.u[!is.na(row.u)]
        row.u
      })
      u
    }
  })
  yearlySiteCovs <- as.data.frame(yearlySiteCovs[!sapply(yearlySiteCovs,
                                                         is.null)])
  if(nrow(yearlySiteCovs) == 0) yearlySiteCovs <- NULL

  # Extract siteCovs and yearlySiteCovs from obsvars
  finalobsvars.df <- as.data.frame(obsvars.df[, !(names(obsvars.df) %in%
                                                    c(names(siteCovs),
                                                      names(yearlySiteCovs)))])
  names(finalobsvars.df) <- names(obsvars.df)[!(names(obsvars.df) %in%
                                                  c(names(siteCovs),
                                                    names(yearlySiteCovs)))]

  emf <- eFrameRGP(y = y, siteCovs = siteCovs, primaryCovs = yearlySiteCovs,
                           numPrimary = nY)
  return(emf)
}

#-------------------------------------
#Inverts Hessian. Returns blank matrix with a warning on a failure.
invertHessian <- function(optimOut, nparam, SE){

  blankMat <- matrix(NA, nparam, nparam)
  if(!SE) return(blankMat)

  tryCatch(solve(optimOut$hessian),
           error=function(e){
             warning("Hessian is singular. Try providing starting values or using fewer covariates.", call.=FALSE)
             return(blankMat)
           })
}


# For multinomial Open. Calculate time intervals acknowledging gaps due to NAs
# The first column indicates is time since first primary period + 1
formatDelta <- function(d, yna)
{
  M <- nrow(yna)
  T <- ncol(yna)
  d <- d - min(d, na.rm=TRUE) + 1
  dout <- matrix(NA, M, T)
  dout[,1] <- d[,1]
  dout[,2:T] <- t(apply(d, 1, diff))
  for(i in 1:M) {
    if(any(yna[i,]) & !all(yna[i,])) { # 2nd test for simulate
      last <- max(which(!yna[i,]))
      y.in <- yna[i, 1:last]
      d.in <- d[i, 1:last]
      if(any(y.in)) {
        for(j in last:2) { # first will always be time since 1
          nextReal <- which(!yna[i, 1:(j-1)])
          if(length(nextReal) > 0)
            dout[i, j] <- d[i, j] - d[i, max(nextReal)]
          else
            dout[i, j] <- d[i, j] - 1
        }
      }
    }
  }
  return(dout)
}

#---
sd.trim <- function(x, trim=0, na.rm=FALSE, ...)
{
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(NA_real_)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(!is.numeric(trim) || length(trim) != 1)
    stop("'trim' must be numeric of length one")
  n <- length(x)
  if(trim > 0 && n > 0) {
    if(is.complex(x)) stop("trimmed sd are not defined for complex data")
    if(trim >= 0.5) return(0)
    lo <- floor(n * trim) + 1
    hi <- n + 1 - lo
    x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  sd(x)
}
