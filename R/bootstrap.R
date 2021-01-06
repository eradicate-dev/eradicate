
#' bootN
#'
#' @description
#' \code{bootN} estimates abundance from a fitted model for the region covered by
#' the sampled sites. Standard errors and confidence intervals are estimated
#' using a bootstrap procedure by resampling sites.
#'
#' @param obj A fitted model object.
#' @param B The number of bootstrap replicates (default = 50).
#'    There must by covariate values for every parameter in \code{obj}.
#' @param trim The fraction of bootstrap samples to be trimmed before computation of the mean
#' and standard errors. Samples are trimmed symmetrically from the ends of the
#' bootstrap distribution.
#' @return a \code{data.frame} giving the predictions for each spatial unit
#'  as well as the overall abundance estimate for the region with associated
#'  SE and confidence intervals.
#'
#' @examples
#'  rem<- san_nic_open$removals
#'  emf <- eFrameMNO(y=rem, numPrimary=8)
#'  mod <- remMNO(~1, ~1, ~1, ~1, data=emf)
#'  Nhat<- bootN(mod)
#'
#' @export
#'
bootN <- function(obj, ...){
  # method generic
  UseMethod("bootN", obj)
}

#' @rdname bootN
#' @export
bootN.efitMNO <- function(obj, B = 50, trim = 0.1, ...) {
  if (B <= 0) {
    return(obj)
  }

  data <- obj$data
  lamformula <- obj$lamformula
  detformula <- obj$detformula
  gamformula <- obj$gamformula
  omformula <- obj$omformula
  iotaformula <- obj$iotaformula

  D <- getDesign(data, lamformula, detformula, gamformula, omformula, iotaformula)
  removed.sites <- D$removed.sites
  if(length(removed.sites)>0)
    data <- data[-removed.sites]
  M <- numSites(data)
  T <- data$numPrimary
  boot.iter <- function() {
    sites <- sort(sample(1:M, M, replace = TRUE))
    data.b <- data[sites]
    mod <- update(obj, data = data.b, se = FALSE)
    fvals<- fitted(mod)
    fvals
  }

  boot.samples <- replicate(B, boot.iter(), simplify = FALSE)
  boot.vec<- sapply(boot.samples, as.vector)
  cell.N <- apply(boot.vec,1,mean,trim=trim)
  cell.se <- apply(boot.vec,1,sd_trim)
  .site<- rep(1:M, T)
  .season <- rep(1:T, each=M)
  cellpreds <- data.frame(N = cell.N, se = cell.se, .site = .site, .season = .season)
  fitted.vals<- fitted(obj)
  Nhat<- apply(fitted.vals, 2, sum)
  boot.sum <- sapply(boot.samples, function(x) apply(x ,2, sum))
  Nboot<- apply(boot.sum, 1, mean, trim = trim)
  boot.sd<- apply(boot.sum, 1, sd_trim, trim=trim)
  boot.lcl<- apply(boot.sum, 1, quantile, 0.025)
  boot.ucl<- apply(boot.sum, 1, quantile,0.975)

  df<- data.frame(Nhat=round(Nhat),Nboot=round(Nboot), .season = 1:T, se=round(boot.sd,1),
                  lwr=round(boot.lcl),upr=round(boot.ucl))
  list(cellpreds=cellpreds, Nhat=df)
}


