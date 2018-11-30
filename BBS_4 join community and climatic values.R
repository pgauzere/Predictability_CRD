######################################################################
############# 4. Join Community/Climate temperature   ################
###################################### Gauzere  ######################

setwd("~/PIERRE/TRANSIENT DYNAMICS/")

library(tidyverse)
load("results/BirdsClim_means_weighted.Rdata")
BirdsClim_means_weighted$year <- as.numeric(as.character(BirdsClim_means_weighted$year))
load("results/bbsRouteClimBreed.Rdata")


load("results/bbsRouteClimBreed_Tmean.Rdata")
bbsRouteTmaxBreed <- bbsRouteClimBreed[bbsRouteClimBreed$climvar == 'tmax', -5]
bbsCommClim <-merge(BirdsClim_means_weighted[,1:4], bbsRouteTmaxBreed, by=c('route', 'year'))

colnames(bbsCommClim)[c(1,2,7,8)] <-c('site', 'time','ClimTmax.ave', 'ClimTmax.sd') 
## we only keep site for which we have more than 5 time steps
sitesToKeep <-bbsCommClim %>% group_by(site) %>% dplyr::summarize(n = length(unique(time)))   %>% filter(n>5) %>% droplevels
CommClim2 <- bbsCommClim %>% na.omit() %>% filter(site %in% sitesToKeep$site) %>% droplevels() %>% arrange(site, time)



#### smooth Ct and Ft ####
## loess independently for Ct and Ft
library(broom)
library(mgcv)
CommClim2<- CommClim2 %>% group_by(site) %>% mutate(CommTmax.ave_loess = augment(loess(CommTmax.ave ~ time))[[3]],
                                                    CommTmax.ave_loess_se=augment(loess(CommTmax.ave ~ time))[[4]],
                                                    ClimTmax.ave_loess=augment(loess(ClimTmax.ave ~ time))[[3]],
                                                    ClimTmax.ave_loess_se=augment(loess(ClimTmax.ave ~ time))[[4]])


#### compute summary statistics ####
source("scripts/PG_Summary Statistics_2.0.R")
library(plyr)
statistics_sites <-plyr::ddply(CommClim2, plyr::.(site, Longi, Lati),function(x){
  summarizeftct(ft=x$ClimTmax.ave_loess, ct=x$CommTmax.ave_loess, ft_error = x$ClimTmax.ave_loess_se, ct_error=x$CommTmax.ave_loess_se, time=x$time, name=x$site)}, .progress = "text") %>% distinct()

save(statistics_sites, file="results/statistics_sites_BBS.Rdata")

CommClimStat2 <- left_join(CommClim2, statistics_sites)
# rm(CommClim2)

test <- CommClimStat2%>% group_by(site) %>% dplyr::summarise(nTime=length(unique(time)))
hist(test$nTime, breaks=45)
min(test$nTime)
CommClimStat2 <- CommClimStat2 %>% arrange(site, desc(time))
# CommClimStat <- droplevels(CommClimStat)

save(CommClimStat2, file="results/CommClimStat2_smoothed_BBS.Rdata")



#####alternative with Tmean ####

load("results/bbsRouteClimBreed_Tmean.Rdata")
# bbsRouteTmaxBreed <- bbsRouteClimBreed_Tmean[bbsRouteClimBreed_Tmean$climvar == 'tmean', -5]
bbsCommClim <-merge(BirdsClim_means_weighted[,1:4], bbsRouteClimBreed_Tmean, by=c('route', 'year'))

colnames(bbsCommClim)[c(1,2,7,8,9)] <-c('site', 'time','ClimPrec.ave', 'ClimTmax.ave', 'ClimTmean.ave') 

## we only keep site for which we have more than 5 time steps
sitesToKeep <- bbsCommClim %>% group_by(site) %>% dplyr::summarize(n = length(unique(time)))   %>% filter(n>5) %>% droplevels
CommClim2   <- bbsCommClim %>% na.omit() %>% filter(site %in% sitesToKeep$site) %>% droplevels() %>% arrange(site, time)

#### smooth Ct and Ft ####
## loess independently for Ct and Ft
library(broom)
library(mgcv)
CommClim2<- CommClim2 %>% group_by(site) %>% mutate(CommTmax.ave_loess = augment(loess(CommTmax.ave ~ time))[[3]],
                                                    CommTmax.ave_loess_se=augment(loess(CommTmax.ave ~ time))[[4]],
                                                    ClimTmax.ave_loess=augment(loess(ClimTmax.ave ~ time))[[3]],
                                                    ClimTmax.ave_loess_se=augment(loess(ClimTmax.ave ~ time))[[4]],
                                                    ClimTmean.ave_loess=augment(loess(ClimTmean.ave ~ time))[[3]],
                                                    ClimTmean.ave_loess_se=augment(loess(ClimTmean.ave ~ time))[[4]])

plot(CommClim2$ClimTmean.ave, CommClim2$ClimTmax.ave)
plot(CommClim2$ClimTmean.ave_loess, CommClim2$ClimTmax.ave_loess)


#### compute summary statistics ####
source("scripts/PG_Summary Statistics_2.0.R")
library(plyr)
statistics_sites <-plyr::ddply(CommClim2, plyr::.(site, Longi, Lati),function(x){
  
  summarizeftct(ft=x$ClimTmean.ave_loess, 
                ct=x$CommTmax.ave_loess, 
                ft_error = x$ClimTmean.ave_loess_se, 
                ct_error=x$CommTmax.ave_loess_se, 
                time=x$time, name=x$site)}, .progress = "text") %>% distinct()

save(statistics_sites, file="results/statistics_sites_BBS_Tmean.Rdata")

CommClimStat2 <- left_join(CommClim2, statistics_sites)
# rm(CommClim2)

test <- CommClimStat2%>% group_by(site) %>% dplyr::summarise(nTime=length(unique(time)))
hist(test$nTime, breaks=45)
min(test$nTime)
CommClimStat2 <- CommClimStat2 %>% arrange(site, desc(time))
# CommClimStat <- droplevels(CommClimStat)

save(CommClimStat2, file="results/CommClimStat2_smoothed_BBS.Rdata")

