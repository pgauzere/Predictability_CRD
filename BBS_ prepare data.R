##############################################################
############# BBS prepare data  #############################
######################################### Gauzere###########

setwd("~/PIERRE/TRANSIENT DYNAMICS/")

library(rgdal)
library(tidyverse)
library(raster)

#### 1.load and check BBS Data ####
load("~/PIERRE/JYB/reviz_Oecol/data/dataframe_comtrue.Rdata")
head(dfCom)

hist(dfCom$N[dfCom$N>0], breaks=200) #8000 !!!
hist(dfCom$N[dfCom$N<100 & dfCom$N>0], breaks=200) #better limit abundance to 100 ind.

#data limited to 100 occurence per site/year
dfComCape<-dfCom
dfComCape$N[dfComCape$N>100]<-100
hist( dfComCape$N [dfComCape$N<100 & dfCom$N>0 ], breaks=100)

#abundance per year
NperYear<- dfComCape %>% na.omit() %>% group_by(year) %>% summarize(n=sum(N))
plot(NperYear)

#site per year
SitePerYear<-dfComCape %>% na.omit() %>% group_by(year) %>% summarize(n=length(unique(route)))
plot(SitePerYear)

#N per Site
NperSite<-dfComCape %>% na.omit() %>% group_by(route) %>% summarize(n=sum(N))
hist(NperSite$n, breaks=100)
NperSite[NperSite$n==max(NperSite$n),]
hist(NperSite$n[NperSite$n!=max(NperSite$n)], breaks=100)
hist(na.omit(dfComCape[dfComCape$route=='840_14_166','N']))

dfComCape <- dfComCape[dfComCape$route!='840_14_166',]
save(dfComCape, file="data/dataframe_comtrue_cape100.Rdata")

#presence/absence data
dfComAbPr<-dfComCape
dfComAbPr$N[dfComAbPr$N>0]<-1
save(dfComAbPr, file="data/dataframe_comtrue_presenceAbsence.Rdata")

#occurence data
dfComOcc<- na.omit(dfComCape[dfComCape$N > 0,])
save(dfComOcc, file="data/dataframe_comtrue_occurence.Rdata")


#N per Site/Year
NperSiteYear<- na.omit(dfComCape) %>% group_by(route, year) %>% summarize(n=sum(N))
hist(NperSiteYear$n, breaks=100)
hist(log(NperSiteYear$n))
