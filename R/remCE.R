#' remCE
#'
#' @name remCE
#'
#' @description
#' \code{remCE} fits various catch-effort models to removal data from a single site
#' (or combined data from many sites).  At least 2 removal periods
#' are required to fit these models.  A port of the \code{deplet} function
#' in the \code{fishmethods} package.
#'
#' @usage remCE(catch, effort, nboot, kbh = NULL, Nstart = NULL)
#'
#' @param catch the number of removals recorded for each period.
#' @param effort the overall effort expended in each period (i.e. trapnights)
#' @param method the depletion method. Variable Effort Models:
#' \code{l}= Leslie estimator,
#' \code{hosc} = sampling coverage estimator for homogeneous model of Chao and Chang (1999),
#' \code{bh} = generalized removal model of Otis et al 1978 (i.e. Model Mbh).
#' @param nboot the number of bootstrap resamples for estimation of standard errors for model hosc,
#' hosc, hesc and hemqle methods.
#' @param kbh the number of capture parameters for model \code{bh}. NULL  for all possible parameters
#' @param Nstart initial value for N for model \code{bh}.  automatically determined if NULL.
#'

#' @return a \code{efit} model object.
#'
#' @examples
#'  remCE(catch=y, effort=eff)
#'
#' @export
#'
remCE<- function (data, method = c("l", "hosc","bh"), nboot=500, kbh=NULL, Nstart=NULL)
{
  if (!any(method %in% c("l", "hosc", "bh")))
    stop("not a valid method.")
  x <- data$counts
  x$samp <- seq(1, length(x$catch), 1)
  x$cpue <- x$catch/x$effort
  x$K <- cumsum(x$catch) - x$catch
  x$F <- cumsum(x$effort) - x$effort
  nsam <- length(x$samp)
  R <- sum(x$catch)

  if (any(method == "l")) {
    if (length(x$catch) < 3)
      stop("Leslie method requires at least 3 observations!")
    ld <- lm(cpue ~ K, data = x)
    names(ld$coefficients)[2] <- "q"
    Nld <- -coef(summary(ld))[1]/coef(summary(ld))[2]
    qld <- -coef(summary(ld))[2]
    s2ld <- summary(ld)$sigma^2
    SEnld <- sqrt((s2ld/qld^2) * ((1/nsam) + (Nld - mean(x$K))^2/sum((x$K - mean(x$K))^2)))
    CInld <- c(Nld - qt(0.975, nsam - 2) * SEnld, Nld + qt(0.975, nsam - 2) * SEnld)
    SEqld <- coef(summary(ld))[4]
    CIqld <- c(qld - qt(0.975, nsam - 2) * SEqld, qld + qt(0.975, nsam - 2) * SEqld)
    ans <- NULL
    ans$results <- matrix(NA, 2L, 6L)
    ans$results <- rbind(cbind(round(Nld, 0),
                               round(SEnld, 1), round(CInld[1], 1), round(CInld[2], 1)),
                         cbind(as.numeric(qld), as.numeric(SEqld), CIqld[1], CIqld[2]))
    dimnames(ans$results) <- list(cbind("N", "q"),
                                  c("Estimate", "SE", "95% LCI",
                                    "95% UCI"))
    ans$gof <- matrix(NA, 2L, 1L)
    ans$gof <- rbind(round(sum((x$catch - predict(ld) * x$effort)^2/
                                 (predict(ld) * x$effort)), 1), as.numeric(summary(ld)$fstatistic[3]))
    dimnames(ans$gof) <- list(c("Chi-square", "df"),"Value")
    ans$anova <- anova(ld)
    ans$summary <- summary(ld)
    ans$residuals <- as.vector(residuals(ld))
    out<- ans
  }

  if (any(method == "hosc")) {
    dd <- x
    dd$del <- ifelse(dd$catch > 0, 1, 0)
    ld <- length(dd$catch)
    for (i in 1:as.numeric(ld - 1)) {
      dd$del[i] <- ifelse(dd$catch[i] == 0 & any(dd$catch[i:as.numeric(ld)] >
                                                   0), 1, dd$del[i])
    }
    dd <- dd[dd$del == 1, ]
    kk <- dd
    nsam <- length(dd[, 1])
    dd$Dk <- cumsum(dd$catch)
    dd$ue <- dd$catch/dd$effort
    for (i in 1:as.numeric(nsam - 1)) {
      dd$Cp[i] <- 1 - (dd$ue[i + 1]/(dd$ue[1]))
    }
    dd$check <- ifelse(dd$Cp <= 0, 1, 0)
    dd$check <- ifelse(dd$check[1] == 1, 0, ifelse(dd$check[nsam] ==
                                                     1, 0, dd$check))
    tau <- ifelse(any(dd$check == TRUE), which(dd$check ==
                                                 1) - 1, nsam - 1)
    Nsc0 <- dd$Dk[tau]/dd$Cp[tau]
    q <- dd$catch[1]/(Nsc0 * dd$effort[1])
    prob <- c(dd$catch/Nsc0, 1 - R/Nsc0)
    BootNsc0 <- data.frame(rep = NULL, Nscb = NULL, qscb = NULL)
    for (j in 1:nboot) {
      kk$catch <- rmultinom(1, Nsc0, prob)[1:nsam]
      kk$Dk <- cumsum(kk$catch)
      kk$ue <- kk$catch/kk$effort
      for (i in 1:as.numeric(nsam - 1)) {
        kk$Cp[i] <- 1 - (kk$ue[i + 1]/(kk$ue[1]))
      }
      kk$check <- ifelse(kk$Cp <= 0, 1, 0)
      kk$check <- ifelse(kk$check[1] == 1, 0, ifelse(kk$check[nsam] ==
                                                       1, 0, kk$check))
      tauk <- ifelse(any(kk$check == TRUE), which(kk$check ==
                                                    1) - 1, nsam - 1)
      Nscb <- kk$Dk[tauk]/kk$Cp[tauk]
      qscb <- kk$catch[1]/(Nscb * kk$effort[1])
      BootNsc0[j, 1] <- j
      BootNsc0[j, 2] <- Nscb
      BootNsc0[j, 3] <- qscb
    }
    seN <- sd(BootNsc0[, 2])
    seqq <- sd(BootNsc0[, 3])
    C <- exp(1.96 * sqrt(log(1 + ((seN^2)/((Nsc0 - sum(dd$catch))^2)))))
    LCn <- sum(dd$catch) + (Nsc0 - sum(dd$catch))/C
    UCn <- sum(dd$catch) + (Nsc0 - sum(dd$catch)) * C
    ans <- NULL
    ans$results <- matrix(NA, 2L, 6L)
    ans$results <- rbind(cbind(round(Nsc0, 0), round(seN,
                                                     1), round(LCn, 1), round(UCn, 1)), cbind(q, seqq,
                                                                                              NA, NA))
    dimnames(ans$results) <- list(cbind("N", "lambda"),
                                  c("Estimate", "SE", "95% LCI",
                                    "95% UCI"))
    out <- ans
  }

  if(any(method=="bh")){
    if(length(x$catch) == 2){
      if(x$catch[2] >= x$catch[1]) warning ("For the two-pass estimator,
                                            catch of first pass must be greater
                                            than catch of second pass.")
    }
    nsam<- length(x$catch)
    endk<- ifelse(is.null(kbh), as.numeric(nsam-2), kbh)
    endk<- ifelse(endk > as.numeric(nsam-2),as.numeric(nsam-2),endk)
    if(nsam==2){
      warning ("Two pass estimator and variance of Otis et al. 1978 used.")
      N2<- x$catch[1]/(1-(x$catch[2]/x$catch[1]))
      q1<- 1-(x$catch[2]/x$catch[1])
      SEO<- sqrt(((N2*(1-q1)^nsam)*(1-(1-q1)^nsam))/((1-(1-q1)^nsam)^2-(nsam^2*q1^2*(1-q1)^(nsam-1))))
      ans<- NULL
      ans$two_pass<- matrix(NA,2L,2L)
      ans$two_pass<- rbind(cbind(round(N2,1),round(SEO,2)),
                           cbind(q1,NA))
      dimnames(ans$two_pass)<- list(cbind("N","q1"),c("Estimate","SE"))
    }
    if(nsam > 2){
      pcatch1<- c(rep(0,nsam))
      T<- c(rep(0,nsam))
      qq<- seq(1,nsam,1)
      T[c(seq(1,nsam,1))]<- cumsum(x$catch[c(seq(1,nsam,1))])
      TC<-sum(x$catch)
      if(is.null(Nstart)) N<- TC
      if(!is.null(Nstart)) N<- Nstart
      ans<-NULL
      resids<-NULL
      preds<-NULL
      lls<-NULL

      for(k in 1:endk){# limited to n-p
        if(k == 1) ps<- T[nsam]/(nsam*N-sum(T[1:nsam-1]))
        if(k >= 2) ps<- rep(T[nsam]/(nsam*N-sum(T[1:nsam-1])),k)
        parms<- c(N,ps)
        model<- function(y){
          N<- y[1]
          ps<- y[2:length(y)]
          pp<-NULL
          pp<- array(NA,dim=c(nsam,k))
          for(nc in 1:ncol(pp)){
            for(nr in nc:nsam){
              if(k == nc) pp[nr,nc]<- ps[nc]*(1-ps[nc])^(nr-nc)
              if(nc < k){
                if(nr == nc) pp[nr,nc]<- ps[nc]
                if(nr > nc) pp[nr,nc]<- 1-ps[nc]
              }
            }
          }
          pps<- apply(pp,1,prod,na.rm=T)
          pcatch1<- N*pps
          pT<- sum(pcatch1)
          L1<- sum(log((seq(1,T[nsam],1)+N-T[nsam])/seq(1,T[nsam],1)))
          ratio<-NULL
          for(f in 1:nsam){
            if(x$catch[f] <= 0 & pcatch1[f] > 0) ratio[f]<- 0
            if(x$catch[f] <= 0 & pcatch1[f] <= 0) ratio[f]<- 0
            if(x$catch[f] > 0 & pcatch1[f] <= 0) ratio[f]<- 0
            if(x$catch[f] > 0 & pcatch1[f] > 0) ratio[f]<- log(x$catch[f]/pcatch1[f])
          }

          H<- sum(x$catch*ratio)
          dNTP<- ifelse(as.numeric(N-pT)<=0,0,log(N-pT))
          G<- N*log(N)-T[nsam]*log(T[nsam])-(N-T[nsam])*dNTP-L1
          G+H
        }
        lower<- c(sum(x$catch),rep(1e-15,k))
        upper<- c(sum(x$catch)*20,rep(1,k)) #upper limits of N and q
        j<-0
        # try different starting value if model doesn't fit at first
        repeat{
          results<-try(optim(parms, model, gr = NULL,lower=lower,upper=upper,method=c("L-BFGS-B"),
                             control=list(maxit=100000),hessian=TRUE),TRUE)
          if(is(results,"try-error")){
            j<-j+1
            parms<-c(parms[1]*(1+j/10),parms[2:length(parms)])
            if(j==30) break
          }
          if(!is(results,"try-error")) break
        }

        if(!is(results,"try-error")){
          cov1<- solve(results$hessian)
          # in case SE is estimated inadequately
          var<- ifelse(diag(cov1)<0,NA,diag(cov1))
          SEs<- sqrt(var)
          ps<- results$par[2:length(results$par)]
          pp<- NULL
          pp<- array(NA,dim=c(nsam,length(ps)))
          for(nc in 1:ncol(pp)){
            for(nr in nc:nsam){
              if(k == nc) pp[nr,nc]<-ps[nc]*(1-ps[nc])^(nr-nc)
              if(nc < k){
                if(nr == nc) pp[nr,nc]<-ps[nc]
                if(nr > nc) pp[nr,nc]<-1-ps[nc]
              }
            }
          }
          pps<- apply(pp,1,prod,na.rm=T)
          pcatch1<- results$par[1]*pps
          resid<- x$catch - pcatch1
          chi<- ifelse(is.nan(((x$catch-pcatch1)^2)/pcatch1),0,((x$catch-pcatch1)^2)/pcatch1)
          chi<- sum(chi)
          df<- nsam - length(results$par)
          prob<- 1-pchisq(chi,df)
          aic<- 2*results$value+2*length(results$par)
        }
        if(is(results,"try-error")){
          pcatch1<- rep(NA,nsam)
          resid<- rep(NA,nsam)
          chi<- NA
          df<- NA
          prob<- NA
          SEs<- rep(NA,length(parms))
        }
        # Create CAPTURE TABLE
        plabs<- paste("p",seq(1,length(ps),1),sep="")
        if(!is(results,"try-error")) comb<- cbind(round(results$par,4),round(SEs,4))
        if(is(results,"try-error")) comb<- cbind(rep(NA,length(parms)),round(SEs,4))
        rownames(comb)<- c("N",plabs)
        comb<- as.data.frame(comb)
        names(comb)<- c("Estimate","SE")
        comb$chi<- c(round(chi,3),rep(NA,length(ps)))
        comb$df<- c(df,rep(NA,length(ps)))
        comb$prob<- c(round(prob,4),rep(NA,length(ps)))
        # Create Residuals, Observed, Predicted Table
        resids<- cbind(resids,resid)
        preds<- cbind(preds,pcatch1)
        lls<- rbind(lls,aic)
        ans[[k]]<- comb
      }#Kloop
      prc<- as.matrix(cbind(x$catch,preds))
      klabs<- paste("k=",seq(1,endk,1),sep="")
      colnames(prc)<- c("catch",klabs)
      res<- as.matrix(resids)
      colnames(res)<- c(klabs)
      kl<- paste("k",seq(1,endk,1),sep="")
      names(ans)<- c(kl)
      aicout<- as.matrix(as.vector(lls),rownames=NA)
      colnames(aicout)<- c("AIC Value")
      rownames(aicout)<- c(seq(1,endk,1))
      ans$aic<- aicout
      ans$predicted<- prc
      ans$residuals<- res
    }
    out<- ans
  }
  out
}
