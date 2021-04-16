

### Ulrike Herzschuh, Thomas Böhmer, Xianyong Cao, Raphaël Hébert, Anne Dallmeyer, Richard Telford, Stefan Kruse
###
### Reversals in temperature-precipitation interdependencies in the Northern Hemisphere extra-tropics
###
### calculation of timescale-dependent correlation, correlation significances and slope of correlation from reconstructed proxy data
###
### 
### Correspondence to: Ulrike.Herzschuh@awi.de
###
### Last edit: 2021-03-15

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

# install the corit-package from GitHub if not done yet:

if (!require("devtools")) {   
  install.packages("devtools")   
}   
devtools::install_github("EarthSystemDiagnostics/corit")

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

library(corit)

source("/ModifiedCorit.R")

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

recon.4_0  <- read.csv("/proxy_wapls_reconstruction_70taxa_4-0ka.csv", header=TRUE)

recon.8_4  <- read.csv("/proxy_wapls_reconstruction_70taxa_8-4ka.csv", header=TRUE)

recon.12_8 <- read.csv("/proxy_wapls_reconstruction_70taxa_12-8ka.csv", header=TRUE)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### 4-0 ka:

{
  
  ids.4_0 <- unique(recon.4_0$Dataset_ID)
  
  corr.4_0  <- data.frame(Dataset_ID=NA, Long=NA, Lat=NA, corr.100=NA, pval.100=NA, corr.250=NA, pval.250=NA, corr.500=NA, pval.500=NA, corr.1000=NA, pval.1000=NA)
  slope.4_0 <- data.frame(Dataset_ID=NA, Long=NA, Lat=NA, slope.100=NA, slope.250=NA, slope.500=NA, slope.1000=NA)
  
  
  for (i in 1:length(ids.4_0)) {
    
    # subset a single site from the dataset:
    reconsubset <- subset(recon.4_0, Dataset_ID == ids.4_0[i])
    uniquereconsubset <- reconsubset[which(!duplicated(reconsubset$AWI_meanAgeBP)), ]
    
    print(paste(i,"/",length(ids.4_0),sep=""))
    
    # create zoo time-series for tailored climate parameters:
    tjul.series.4_0 <- zoo(uniquereconsubset$TJul.pann_tailored, order.by=uniquereconsubset$AWI_meanAgeBP)
    pann.series.4_0 <- zoo(uniquereconsubset$Pann.tjul_tailored, order.by=uniquereconsubset$AWI_meanAgeBP)
    
    # determine the temporal resolution of the site:
    agedist.4_0 <- diff(uniquereconsubset$AWI_meanAgeBP)
    
    # -------------------------------------------------------------------------------------------------
    
    corr.4_0[i,"Dataset_ID"]  <- unique(uniquereconsubset$Dataset_ID)
    corr.4_0[i,"Long"]        <- unique(uniquereconsubset$Long)
    corr.4_0[i,"Lat"]         <- unique(uniquereconsubset$Lat)
    
    slope.4_0[i,"Dataset_ID"]  <- unique(uniquereconsubset$Dataset_ID)
    slope.4_0[i,"Long"]        <- unique(uniquereconsubset$Long)
    slope.4_0[i,"Lat"]         <- unique(uniquereconsubset$Lat) 
    
    # -------------------------------------------------------------------------------------------------
    
    if(length(uniquereconsubset$Dataset_ID) > 5 & all(!is.na(tjul.series.4_0)) & all(!is.na(pann.series.4_0))) {
      
      # -------------------------------------------------------------------------------------------------
      
      # calculate correlation, correlation significance and slope of correlation at 100 years time-scale with 10-year interpolation time step:
      
      if(mean(agedist.4_0) < 100) {
        
        modcorit100.4_0 <- ModifiedCorit(timser1=tjul.series.4_0, timser2=pann.series.4_0, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                         fc=1/100, dt=100/(0.1*100), int.method="linear", filt.output=FALSE)
        
        corr.4_0[i,"corr.100"] <- modcorit100.4_0[1]
        corr.4_0[i,"pval.100"] <- modcorit100.4_0[2]
        
        slope.4_0[i,"slope.100"] <- modcorit100.4_0[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      # calculate correlation, correlation significance and slope of correlation at 250 years time-scale with 10-year interpolation time step:
      
      if(mean(agedist.4_0) < 250) {
        
        modcorit250.4_0 <- ModifiedCorit(timser1=tjul.series.4_0, timser2=pann.series.4_0, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                         fc=1/250, dt=250/(0.1*250), int.method="linear", filt.output=FALSE)
        
        corr.4_0[i,"corr.250"] <- modcorit250.4_0[1]
        corr.4_0[i,"pval.250"] <- modcorit250.4_0[2]
        
        slope.4_0[i,"slope.250"] <- modcorit250.4_0[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      # calculate correlation, correlation significance and slope of correlation at 500 years time-scale with 10-year interpolation time step:
      
      if(mean(agedist.4_0) < 500) {
        
        modcorit500.4_0 <- ModifiedCorit(timser1=tjul.series.4_0, timser2=pann.series.4_0, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                         fc=1/500, dt=500/(0.1*500), int.method="linear", filt.output=FALSE)
        
        corr.4_0[i,"corr.500"] <- modcorit500.4_0[1]
        corr.4_0[i,"pval.500"] <- modcorit500.4_0[2]
        
        slope.4_0[i,"slope.500"] <- modcorit500.4_0[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      # calculate correlation, correlation significance and slope of correlation at 1000 years time-scale with 10-year interpolation time step:
      
      if(mean(agedist.4_0) < 1000) {
        
        modcorit1000.4_0 <- ModifiedCorit(timser1=tjul.series.4_0, timser2=pann.series.4_0, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                          fc=1/1000, dt=1000/(0.1*1000), int.method="linear", filt.output=FALSE)
        
        corr.4_0[i,"corr.1000"] <- modcorit1000.4_0[1]
        corr.4_0[i,"pval.1000"] <- modcorit1000.4_0[2]
        
        slope.4_0[i,"slope.1000"] <- modcorit1000.4_0[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
    }
    
  } #end for
  
  
  write.csv(corr.4_0,  file="/proxy_wapls_tailored_correlation_significances_timescales_4-0ka.csv", row.names=FALSE)
  write.csv(slope.4_0, file="/proxy_wapls_tailored_slope_timescales_4-0ka.csv", row.names=FALSE)
  
  # -------------------------------------------------------------------------------------------------
  
} #end 4-0ka  

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### 8-4 ka:

{
  
  ids.8_4 <- unique(recon.8_4$Dataset_ID)
  
  corr.8_4  <- data.frame(Dataset_ID=NA, Long=NA, Lat=NA, corr.100=NA, pval.100=NA, corr.250=NA, pval.250=NA, corr.500=NA, pval.500=NA, corr.1000=NA, pval.1000=NA)
  slope.8_4 <- data.frame(Dataset_ID=NA, Long=NA, Lat=NA, slope.100=NA, slope.250=NA, slope.500=NA, slope.1000=NA)
  
  
  for (i in 1:length(ids.8_4)) {
    
    reconsubset <- subset(recon.8_4, Dataset_ID == ids.8_4[i])
    uniquereconsubset <- reconsubset[which(!duplicated(reconsubset$AWI_meanAgeBP)), ]
    
    print(paste(i,"/",length(ids.8_4),sep=""))
    
    tjul.series.8_4 <- zoo(uniquereconsubset$TJul.pann_tailored, order.by=uniquereconsubset$AWI_meanAgeBP)
    pann.series.8_4 <- zoo(uniquereconsubset$Pann.tjul_tailored, order.by=uniquereconsubset$AWI_meanAgeBP)
    
    agedist.8_4 <- diff(uniquereconsubset$AWI_meanAgeBP)
    
    # -------------------------------------------------------------------------------------------------
    
    corr.8_4[i,"Dataset_ID"]  <- unique(uniquereconsubset$Dataset_ID)
    corr.8_4[i,"Long"]        <- unique(uniquereconsubset$Long)
    corr.8_4[i,"Lat"]         <- unique(uniquereconsubset$Lat)
    
    slope.8_4[i,"Dataset_ID"]  <- unique(uniquereconsubset$Dataset_ID)
    slope.8_4[i,"Long"]        <- unique(uniquereconsubset$Long)
    slope.8_4[i,"Lat"]         <- unique(uniquereconsubset$Lat) 
    
    # -------------------------------------------------------------------------------------------------
    
    if(length(uniquereconsubset$Dataset_ID) > 5 & all(!is.na(tjul.series.8_4)) & all(!is.na(pann.series.8_4))) {
      
      if(mean(agedist.8_4) < 100) {
        
        modcorit100.8_4 <- ModifiedCorit(timser1=tjul.series.8_4, timser2=pann.series.8_4, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                         fc=1/100, dt=100/(0.1*100), int.method="linear", filt.output=FALSE)
        
        corr.8_4[i,"corr.100"] <- modcorit100.8_4[1]
        corr.8_4[i,"pval.100"] <- modcorit100.8_4[2]
        
        slope.8_4[i,"slope.100"] <- modcorit100.8_4[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      if(mean(agedist.8_4) < 250) {
        
        modcorit250.8_4 <- ModifiedCorit(timser1=tjul.series.8_4, timser2=pann.series.8_4, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                         fc=1/250, dt=250/(0.1*250), int.method="linear", filt.output=FALSE)
        
        corr.8_4[i,"corr.250"] <- modcorit250.8_4[1]
        corr.8_4[i,"pval.250"] <- modcorit250.8_4[2]
        
        slope.8_4[i,"slope.250"] <- modcorit250.8_4[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      if(mean(agedist.8_4) < 500) {
        
        modcorit500.8_4 <- ModifiedCorit(timser1=tjul.series.8_4, timser2=pann.series.8_4, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                         fc=1/500, dt=500/(0.1*500), int.method="linear", filt.output=FALSE)
        
        corr.8_4[i,"corr.500"] <- modcorit500.8_4[1]
        corr.8_4[i,"pval.500"] <- modcorit500.8_4[2]
        
        slope.8_4[i,"slope.500"] <- modcorit500.8_4[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      if(mean(agedist.8_4) < 1000) {
        
        modcorit1000.8_4 <- ModifiedCorit(timser1=tjul.series.8_4, timser2=pann.series.8_4, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                          fc=1/1000, dt=1000/(0.1*1000), int.method="linear", filt.output=FALSE)
        
        corr.8_4[i,"corr.1000"] <- modcorit1000.8_4[1]
        corr.8_4[i,"pval.1000"] <- modcorit1000.8_4[2]
        
        slope.8_4[i,"slope.1000"] <- modcorit1000.8_4[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
    }
    
  } #end for
  
  
  write.csv(corr.8_4,  file="/proxy_wapls_tailored_correlation_significances_timescales_8-4ka.csv", row.names=FALSE)
  write.csv(slope.8_4, file="/proxy_wapls_tailored_slope_timescales_8-4ka.csv", row.names=FALSE)
  
  # -------------------------------------------------------------------------------------------------
  
} #end 8-4ka  

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### 12-8 ka:

{
  
  ids.12_8 <- unique(recon.12_8$Dataset_ID)
  
  corr.12_8  <- data.frame(Dataset_ID=NA, Long=NA, Lat=NA, corr.100=NA, pval.100=NA, corr.250=NA, pval.250=NA, corr.500=NA, pval.500=NA, corr.1000=NA, pval.1000=NA)
  slope.12_8 <- data.frame(Dataset_ID=NA, Long=NA, Lat=NA, slope.100=NA, slope.250=NA, slope.500=NA, slope.1000=NA)
  
  
  for (i in 1:length(ids.12_8)) {
    
    reconsubset <- subset(recon.12_8, Dataset_ID == ids.12_8[i])
    uniquereconsubset <- reconsubset[which(!duplicated(reconsubset$AWI_meanAgeBP)), ]
    
    print(paste(i,"/",length(ids.12_8),sep=""))
    
    tjul.series.12_8 <- zoo(uniquereconsubset$TJul.pann_tailored, order.by=uniquereconsubset$AWI_meanAgeBP)
    pann.series.12_8 <- zoo(uniquereconsubset$Pann.tjul_tailored, order.by=uniquereconsubset$AWI_meanAgeBP)
    
    agedist.12_8 <- diff(uniquereconsubset$AWI_meanAgeBP)
    
    # -------------------------------------------------------------------------------------------------
    
    corr.12_8[i,"Dataset_ID"]  <- unique(uniquereconsubset$Dataset_ID)
    corr.12_8[i,"Long"]        <- unique(uniquereconsubset$Long)
    corr.12_8[i,"Lat"]         <- unique(uniquereconsubset$Lat)
    
    slope.12_8[i,"Dataset_ID"]  <- unique(uniquereconsubset$Dataset_ID)
    slope.12_8[i,"Long"]        <- unique(uniquereconsubset$Long)
    slope.12_8[i,"Lat"]         <- unique(uniquereconsubset$Lat) 
    
    # -------------------------------------------------------------------------------------------------
    
    if(length(uniquereconsubset$Dataset_ID) > 5 & all(!is.na(tjul.series.12_8)) & all(!is.na(pann.series.12_8))) {
      
      if(mean(agedist.12_8) < 100) {
        
        modcorit100.12_8 <- ModifiedCorit(timser1=tjul.series.12_8, timser2=pann.series.12_8, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                          fc=1/100, dt=100/(0.1*100), int.method="linear", filt.output=FALSE)
        
        corr.12_8[i,"corr.100"] <- modcorit100.12_8[1]
        corr.12_8[i,"pval.100"] <- modcorit100.12_8[2]
        
        slope.12_8[i,"slope.100"] <- modcorit100.12_8[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      if(mean(agedist.12_8) < 250) {
        
        modcorit250.12_8 <- ModifiedCorit(timser1=tjul.series.12_8, timser2=pann.series.12_8, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                          fc=1/250, dt=250/(0.1*250), int.method="linear", filt.output=FALSE)
        
        corr.12_8[i,"corr.250"] <- modcorit250.12_8[1]
        corr.12_8[i,"pval.250"] <- modcorit250.12_8[2]
        
        slope.12_8[i,"slope.250"] <- modcorit250.12_8[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      if(mean(agedist.12_8) < 500) {
        
        modcorit500.12_8 <- ModifiedCorit(timser1=tjul.series.12_8, timser2=pann.series.12_8, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                          fc=1/500, dt=500/(0.1*500), int.method="linear", filt.output=FALSE)
        
        corr.12_8[i,"corr.500"] <- modcorit500.12_8[1]
        corr.12_8[i,"pval.500"] <- modcorit500.12_8[2]
        
        slope.12_8[i,"slope.500"] <- modcorit500.12_8[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
      if(mean(agedist.12_8) < 1000) {
        
        modcorit1000.12_8 <- ModifiedCorit(timser1=tjul.series.12_8, timser2=pann.series.12_8, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss",
                                           fc=1/1000, dt=1000/(0.1*1000), int.method="linear", filt.output=FALSE)
        
        corr.12_8[i,"corr.1000"] <- modcorit1000.12_8[1]
        corr.12_8[i,"pval.1000"] <- modcorit1000.12_8[2]
        
        slope.12_8[i,"slope.1000"] <- modcorit1000.12_8[3]
        
      }
      
      # -------------------------------------------------------------------------------------------------
      
    }
    
  } #end for
  
  
  write.csv(corr.12_8,  file="/proxy_wapls_tailored_correlation_significances_timescales_12-8ka.csv", row.names=FALSE)
  write.csv(slope.12_8, file="/proxy_wapls_tailored_slope_timescales_12-8ka.csv", row.names=FALSE)
  
  # -------------------------------------------------------------------------------------------------
  
} #end 12-8ka  

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

