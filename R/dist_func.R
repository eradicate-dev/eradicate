#' dist_func
#'
#' @name dist_func
#'
#' @description
#' \code{dist_func} estimates a distance function from distance data collected from
#' camera traps.  The distance function is then used to estimate the effective area
#' of the camera viewshed by integration. Data are assumed to be in bins specified in
#'  \code{cutpoints} (see \code{eFrameDS} for more details).  The detection function is
#'  estimated using the function \code{ds()} from package \code{Distance} using point
#'  distance sampling. Currently no adjustments are included in the fitting.
#'
#' @usage dist_func(data, theta, keyfun = c("hn","hr") ...)
#'
#' @param data \code{eFrameDS} object containing the distance data and bin information.
#' @param formula formula for the scale parameter. Default is ~1.
#' @param theta angle of the camera viewshed in degrees
#' @param keyfun detection functions to fit (currently only half-normal ("hn") or
#' hazard rate ("hr")). No adjustments are currently permitted.
#'
#' @return a \code{efitDS} model object.
#'
#' @examples
#'  ddata<- HogDeer$ddata
#'  cutpoints<- HogDeer$cutpoints
#'  emf<- eFrameDS(ddata$distance, ddata$size, ddata$camID, cutpoints, w=12.5, bin_num=TRUE)

#'  mods <- dist_func(emf, theta=HogDeer$theta)
#'
#' @export
#'
dist_func<- function (data, formula = ~1, theta, keyfun = c("hn","hr"), ...){
  if(!is(data, "eFrameDS"))
    stop("Data is not a eFrameDS")

    ddata<- emf$ddata
    cutpoints<- emf$cutpoints
    w<- emf$w
    fits<- list()
    fit.table<- NULL
    AIC_values<- NULL
    for(kf in keyfun) {
      fm<- suppressMessages(Distance::ds(ddata, formula = formula, key=kf, transect="point",
                      truncation=w, adjustment=NULL, cutpoints=cutpoints))
      summ<- summary(fm)
      fits[[kf]]<- list(pars=summ$ds$coeff)
      aic<- AIC(fm)$AIC
      Pa<- summ$ds$average.p
      Pa.se<- summ$ds$average.p.se
      esa.pred <- predict(fm, newdata=data.frame(distance=0), esw=TRUE, se.fit=TRUE)
      esa<- esa.pred$fitted * theta/360 # effective viewshed area
      esa.se <- sqrt(esa.pred$se.fit^2 * (theta/360)^2)
      chip<- Distance::gof_ds(fm)$chisquare$chi1$p
      res<- data.frame(keyfun=kf,Pa=Pa, Pa.se=Pa.se,esa=esa,esa.se=esa.se,p=chip,AIC=aic)
      fit.table<- rbind(fit.table, res)
    }

    efit <- list(fitType = "dist_func", fits=fits, fit.table=fit.table)
    class(efit) <- c('efitDS','efit','list')
    return(efit)
}






