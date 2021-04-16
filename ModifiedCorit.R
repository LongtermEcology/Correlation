ModifiedCorit<-function (timser1, timser2, fc,detr=FALSE, method = "InterpolationMethod", appliedFilter = "gauss", rep=100,quant=seq(0.005,0.995,0.005),pvalue=TRUE,
                         tn = seq(from = 10, to = max(c(index(timser1),index(timser2))), by = 10), dt=1/fc*0.1, int.method = "linear", k = 5, filt.output=FALSE)
{
  n <- max(c(max(index(timser1)), max(index(timser2))))
  dofs<-mean(c(length(timser1),length(timser2)))
  if (method == "InterpolationMethod") {
    filtTimser1 <- corit::InterpolationMethod(corit::detrTimser(timser1,
                                                  detr), fc, dt, n, int.method, appliedFilter, k)
    filtTimser2 <- corit::InterpolationMethod(corit::detrTimser(timser2,
                                                  detr), fc, dt, n, int.method, appliedFilter, k)
    res <- cor.test(filtTimser1, filtTimser2, method = "pearson",
                    alternative = "two.sided", na.action = TRUE)$estimate
    drop<-which(is.na(filtTimser2) | is.na(filtTimser1))
    if(length(drop)>0){
    ft2<-filtTimser2[-drop]
    ft1<-filtTimser1[-drop]
    }else{
      ft2<-filtTimser2
      ft1<-filtTimser1
    }
    lm.filt<-lm(ft2 ~ ft1)
    slope<-coef(lm.filt)[2]/mean(ft2)*100
  }
  if (method == "DirectFiltering") {
    filtTimser1 <- corit::DirectFiltering(corit::detrTimser(timser1, detr),
                                   fc, tn, appliedFilter, k)
    filtTimser2 <- corit::DirectFiltering(corit::detrTimser(timser2, detr),
                                   fc, tn, appliedFilter, k)
    res <- cor.test(filtTimser1, filtTimser2, method = "pearson",
                    alternative = "two.sided", na.action = TRUE)$estimate
    drop<-which(is.na(filtTimser2) | is.na(filtTimser1))
    if(length(drop)>0){
      ft2<-filtTimser2[-drop]
      ft1<-filtTimser1[-drop]
    }else{
      ft2<-filtTimser2
      ft1<-filtTimser1
    }
    lm.filt<-lm(ft2 ~ ft1)
    slope<-coef(lm.filt)[2]/mean(ft2)*100
  }
  if (method == "IntegrandInterpolationMethod") {
    filtTimser1 <- corit::IntegrandInterpolationMethod(corit::detrTimser(timser1,
                                                           detr), fc, tn, appliedFilter, k)
    filtTimser2 <- corit::IntegrandInterpolationMethod(corit::detrTimser(timser2,
                                                           detr), fc, tn, appliedFilter, k)
    res <- cor.test(filtTimser1, filtTimser2, method = "pearson",
                    alternative = "two.sided", na.action = TRUE)$estimate
    drop<-which(is.na(filtTimser2) | is.na(filtTimser1))
    if(length(drop)>0){
      ft2<-filtTimser2[-drop]
      ft1<-filtTimser1[-drop]
    }else{
      ft2<-filtTimser2
      ft1<-filtTimser1
    }
    lm.filt<-lm(ft2 ~ ft1)
    slope<-coef(lm.filt)[2]/mean(ft2)*100
  }
  if (filt.output == TRUE) {
    return(list(cor = res, slope = slope,ft1 = filtTimser1, ft2 = filtTimser2))
  }
  if (filt.output == FALSE) {
    if(pvalue){
    spectral.slopes<-corit::estimateTimserSlopes(timeseries1 = timser1,timeseries2 = timser2,int.step = dt)
    quants<-corit::CorQuantilesNullHyp(timser1 = timser1,timser2 = timser2,beta.noise1 = spectral.slopes$s1,beta.noise2 = spectral.slopes$s2,detr = FALSE,rep = rep,quant = quant,method = "InterpolationMethod",appliedFilter = "gauss",fc = fc,dt = dt,int.method = 'linear')
    pval<-if(res>0){1-quant[which.min(abs(quants$Quantile[[1]]-res))]}else{quant[which.min(abs(quants$Quantile[[1]]-res))]}
    names(pval)<-"pval"
    names(slope)<-"slope"
    surr.cor<-mean(quants$corPair[[1]])
    names(surr.cor)<-"surr.cor"
    return(c(res,pval,slope,surr.cor))
    }else{
      return(c(res,slope))
    }
  }
}
