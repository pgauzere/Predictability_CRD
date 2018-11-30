###############################################################################################
############# 3. Compute Community inferred climate  ######################################
#################################################### Gauzere  ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS/")

library(rgdal)
library(tidyverse)
library(raster)


#### 1.load Birds climatic Niches ####
load('results/birdsBreedClimNiches.Rdata')

#### 2.load BBS data ####
load('data/dataframe_comtrue_occurence.Rdata')

#### 3. Compute Community Inferred Climate ####

birdsBreedClimNiches$species <- as.factor(gsub(" ", "_",unique(birdsBreedClimNiches$species)))

dfComOccNiches <- left_join(dfComOcc, birdsBreedClimNiches)

BirdsClim_means_occ   <- dfComOccNiches %>% dplyr::group_by(route, year) %>% summarize( 
                               CommTmax.ave = mean(mean.Monthly.maxTemp, na.rm = T),
                               CommTmax.sd  = sd(mean.Monthly.maxTemp, na.rm = T),
                               CommPrec.ave = mean(mean.Monthly.Prec, na.rm = T),
                               CommPrec.sd  = sd(mean.Monthly.Prec, na.rm = T))

BirdsClim_means_occ   <- plyr::ddply(dfComOccNiches, .(route, year), function(x){
  # x<- dfComOccNiches[dfComOccNiches$route=='840_35_9' & dfComOccNiches$year==1972,]
  CommTmax.ave = mean(x$mean.Monthly.maxTemp, na.rm = T)
  CommTmax.var  = var(x$mean.Monthly.maxTemp, na.rm = T)
  CommPrec.ave = mean(x$mean.Monthly.Prec, na.rm = T)
  CommPrec.var  = var(x$mean.Monthly.Prec, na.rm = T)
  return(cbind(CommTmax.ave, CommTmax.var, CommPrec.ave, CommPrec.var))}, .progress='text')

save(BirdsClim_means_occ, file="results/BirdsClim_means_occurences.Rdata")


BirdsClim_means_weighted   <- plyr::ddply(na.omit(dfComOccNiches), .(route, year), function(x){
  corrective_term = sum(x$N)  / (sum(x$N) -1)
  CommTmax.ave = sum( x$mean.Monthly.maxTemp * x$N / sum(x$N))
  CommTmax.var  = corrective_term * (sum( x$mean.Monthly.maxTemp^2 * x$N / sum(x$N) ) - CommTmax.ave^2)
  CommPrec.ave = sum( x$mean.Monthly.Prec * x$N) / sum(x$N)
  CommPrec.var  = corrective_term * (sum( x$mean.Monthly.Prec^2 * x$N / sum(x$N) ) - CommPrec.ave^2)
  return(cbind(CommTmax.ave, CommTmax.var, CommPrec.ave, CommPrec.var))}, .progress='text')

BirdsClim_means_weighted <- cwi_stratified()

plot(BirdsClim_means_occ$CommTmax.var, BirdsClim_means_weighted$CommTmax.var)

save(BirdsClim_means_weighted, file="results/BirdsClim_means_weighted.Rdata")

plot(BirdsClim_means_occ$CommTmax.ave, BirdsClim_means_weighted$CommTmax.ave)
abline(0,1, col='red')


