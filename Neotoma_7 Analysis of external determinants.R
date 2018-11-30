###################################################
############# 7. Analysis of Predictors ###########
############################## Gauzere & Blonder ##


setwd("~/PIERRE/TRANSIENT DYNAMICS")
library(tidyverse)
library(mgcv)
library(MuMIn)
library(visreg)
source("scripts/base_NAmap.R")

# load("results/statistics_sites_neotoma_SCALED_ext.drivers.Rdata")


#### load datas ####
load("results/statistics_sites_neotoma.Rdata")
load("results/predictors.Rdata")
# load("results/biotic_predictors_CWMs.Rdata")
load("results/biotic_predictors_CWMs_TRY.Rdata")

#### merge datas
temp<-left_join(statistics_sites, predictors)
stpred <- left_join(temp, sitesCWMs)
modData<- na.omit(stpred)
modData$site <- as.factor(modData$site)
modData$was_glaciated <- as.factor(modData$was_glaciated)
modData$sqrt_humanDensity <- sqrt(modData$humanDensity)
modData$sqrt_topo25 <- sqrt(modData$topo25)
modData$sqrt_dist_refugia <- sqrt(modData$dist_refugia)


#### investigate correlation between statistics ####
library("PerformanceAnalytics")
chart.Correlation(modData[, c(
  "crossings.mean.Tmax", "crossings.max.Tmax" , "deviation.base.Tmax", 
  "deviation.absolute.Tmax" , "deviation.change.Tmax"   , "deviation.change.1Ky.Tmax",
  "num.points.not.na.Tmax"  , "lm.intercept.Tmax",   
  "lm.slope.Tmax" , "lm.r2.Tmax"    , "lm.rho.Tmax")], histogram=TRUE, pch=19)

chart.Correlation(modData[, c(
  "crossings.max.Tmax" , "deviation.absolute.Tmax" , 
 "deviation.change.1Ky.Tmax",
  "num.points.not.na.Tmax")], histogram=TRUE, pch=19)

chart.Correlation(statistics_sites[, c(
  "crossings.mean.Prec", "crossings.max.Prec" , "deviation.base.Prec", 
  "deviation.absolute.Prec" , "deviation.change.Prec"   , "deviation.change.1Ky.Prec",
  "num.points.not.na.Prec"  , "lm.intercept.Prec",   
  "lm.slope.Prec" , "lm.r2.Prec"    , "lm.rho.Prec")], histogram=TRUE, pch=19)

chart.Correlation(statistics_sites[, c(
  "crossings.max.Prec" , "deviation.absolute.Prec" , 
  "deviation.change.1Ky.Prec",
  "num.points.not.na.Prec")], histogram=TRUE, pch=19)


#### investigate correlation between predictors ####
chart.Correlation(modData[, c(
                              #"Latitude", 
                              #"Longitude",
                              "base.Tmax",
                              # "change.Tmax",
                              # "thermal.heterogeneity25",
                              "sqrt_topo25",
                              "sqrt_humanDensity",
                              # "dist_refugia",
                              "mean.SeedMass",
                              # "mean.SeedLength",
                              "mean.PlantHeight"
                              # "mean.rangeTmax"
                              )], histogram=TRUE, pch=19)


chart.Correlation(modData[, c("Latitude", 
                              "Longitude",
                              "base.Tmax",
                              "base.Prec",
                              "sqrt_topo25",
                              "sqrt_humanDensity",
                              "wind",
                              "mean.SeedMass",
                              "mean.PlantHeight")], histogram=TRUE, pch=19)


#### distributions of predictors ####
baseTmax <-
  ggplot(modData, aes(x=base.Tmax))+
  geom_histogram(fill="red")+
  theme_classic()


basePrec <-
  ggplot(modData, aes(x=base.Prec))+
  geom_histogram(fill="blue")+
  theme_classic()

topo <-
  ggplot(modData, aes(x=topo25))+
  geom_histogram(fill="brown")+
  theme_classic()+
  scale_x_continuous(trans='sqrt')

refugia <-
  ggplot(modData, aes(x=dist_refugia))+
  geom_histogram(fill="purple")+
  theme_classic()

density <-
  ggplot(modData, aes(x=humanDensity))+
  geom_histogram()+
  theme_classic()+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,20,50,100,200),trans='sqrt')


wind <-
  ggplot(modData, aes(x=wind))+
  geom_histogram(fill="darkgreen")+
  theme_classic()

seedMass <-
  ggplot(modData, aes(x=mean.SeedMass))+
  geom_histogram(fill="darkgreen")+
  theme_classic()

plantHeight <-
  ggplot(modData, aes(x=mean.PlantHeight))+
  geom_histogram(fill="darkorange")+
  theme_classic()

pdf("results/distribution_predictors.pdf", width=10)
baseTmax
basePrec
topo
refugia
density
ice
seedMass
plantHeight
dev.off()



#### maps of predictors ####
library(viridis)
baseTmax <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=base.Tmax),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="baseline max. temperature", x=-65,y=30, size=5)


worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=deviation.absolute.Tmax),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="baseline max. temperature", x=-65,y=30, size=5)

basePrec <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=base.Prec),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="baseline precipitation", x=-65,y=30, size=5)

topo <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=topo25),size=4, shape=21)+
  scale_fill_viridis(trans='sqrt')+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="sqrt. Topography", x=-65,y=30, size=5)

refugia <-
  worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=dist_refugia/1000),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="distance to LGM refugia (km)", x=-65,y=30, size=5)


density <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=humanDensity),size=4, shape=21)+
  scale_fill_viridis(trans='sqrt')+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="sqrt. Human density", x=-65,y=30, size=5)

ice <-
worldmap +
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=was_glaciated),size=4, shape=21)+
  # scale_fill_viridis(trans='sqrt')+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="Ice sheet LGM", x=-65,y=30, size=5)

Wind <-
  worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=wind),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="mean seed mass", x=-65,y=30, size=5)

seedMass <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=mean.SeedMass),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="mean seed mass", x=-65,y=30, size=5)

plantHeight <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=mean.PlantHeight),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="mean plant height", x=-65,y=30, size=5)

pdf("results/maps_predictors.pdf", width=10)
baseTmax
basePrec
topo
density
refugia
seedMass
plantHeight
dev.off()




#### crossings.max Tmax ####
hist(modData$crossings.max.Tmax, breaks=c(0.5,1.5,2.5,3.5,4.5))

## Poisson Linear model
mod.N.Tmax <- glm(crossings.max.Tmax ~ 
                    scale(base.Tmax) + 
                    scale(sqrt_topo25) +
                    scale(sqrt_humanDensity) +
                    scale(dist_refugia) + 
                    scale(mean.SeedMass) + 
                    scale(mean.PlantHeight), family=poisson, data=modData)

# plot(mod.N.Tmax)
# modData$resGLM.N.Tmax<- residuals(mod.N.Tmax)
# hist(modData$resGLM.N.Tmax, breaks=100)
# 
# worldmap + 
#   geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resGLM.N.Tmax), fill=resGLM.N.Tmax), shape=21)+
#   scale_fill_gradient2()+
#   scale_size(guide=FALSE)+
#   theme_classic()
# 
# ggplot(modData, aes(Latitude, resGLM.N.Tmax))+
#   geom_point(shape=21)+
#   theme_classic()+
#   stat_smooth()


## Poisson Linear model + spatial smoothed
modSpa.N.Tmax <- gam(crossings.max.Tmax ~ s(Latitude, Longitude, bs='sos') + base.Tmax + sqrt(topo25) + humanPresence + scale(dist_refugia) + 
                    mean.SeedMass + mean.PlantHeight, family=poisson, data=modData, na.action='na.fail')

plot(modSpa.N.Tmax)
plot(residuals(modSpa.N.Tmax), fitted(modSpa.N.Tmax))
modData$resGLM.Spatial.N.Tmax<- residuals(modSpa.N.Tmax)
hist(modData$resGLM.Spatial.N.Tmax, breaks=100)

worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resGLM.Spatial.N.Tmax), fill=resGLM.Spatial.N.Tmax), shape=21)+
  scale_fill_gradient2()+
  scale_size(guide=FALSE)+
  theme_classic()

AIC(mod.N.Tmax)
AIC(modSpa.N.Tmax)

# ##  Poisson Mixed Linear model
# modMixed.N.Tmax <- gamm(crossings.max.Tmax ~ base.Tmax + topo25 + humanPresence + scale(dist_refugia) + 
#                     mean.SeedMass + mean.PlantHeight, random=list(site=~1), family=poisson, data=modData)
# modData$resGLMixed.N.Tmax<- residuals(modMixed.N.Tmax$lme)
# plot(modMixed.N.Tmax$lme)
# hist(modData$resGLMixed.N.Tmax)
# ggplot(modData, aes(y=Latitude, x=Longitude, size= abs(resGLMixed.N.Tmax), fill=resGLMixed.N.Tmax))+
#   geom_point(shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   coord_fixed(ratio=1)+
#   theme_classic()
# 
# ##  Poisson Mixed Linear model with spatial smoothing
# modMixed.Spatial.N.Tmax <- gamm(crossings.max.Tmax ~ s(Latitude, Longitude, bs='sos') + base.Tmax + topo25 + humanPresence + scale(dist_refugia) + 
#                           mean.SeedMass + mean.PlantHeight, random=list(site=~1), family=poisson, data=modData)
# modData$resGLMixed.Spatial.N.Tmax<- residuals(modMixed.Spatial.N.Tmax$lme)
# plot(modMixed.Spatial.N.Tmax$lme)
# gam.check(modMixed.Spatial.N.Tmax$gam)
# plot(modMixed.Spatial.N.Tmax$gam, page=1)
# hist(modData$resGLMixed.Spatial.N.Tmax)
# 
# ggplot(modData, aes(y=Latitude, x=Longitude, size= abs(resGLMixed.Spatial.N.Tmax), fill=resGLMixed.Spatial.N.Tmax))+
#   geom_point(shape=21)+
#   scale_fill_gradient2()+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   coord_fixed(ratio=1)+
#   theme_classic()
# 
# ggplot(modData, aes( resGLMixed.Spatial.N.Tmax, resGLMixed.N.Tmax, fill=Latitude, size=Latitude))+
#   geom_point(shape=21)+
#   theme_classic()+
#   scale_fill_viridis()+
#   # stat_smooth()+
#   geom_abline(intercept=0, slope=1)
# 
# ggplot(modData, aes(resGLM.N.Tmax, Latitude))+
#   geom_point()+
#   theme_classic()+
#   # scale_fill_viridis()+
#   stat_smooth()+
#   geom_abline(intercept=0, slope=1)


## final model 
# Poisson Linear model + spatial smoothed
modData.N.Tmax <- modData[modData$resGLM.Spatial.N.Tmax < 1.5, ]

modSpa.N.Tmax <- gam(crossings.max.Tmax ~ s(Latitude, Longitude, bs='sos') + 
                       scale(poly(base.Tmax,2)) + 
                       scale(sqrt_topo25) +
                       scale(sqrt_humanDensity) +
                       # scale(dist_refugia) + 
                       scale(mean.SeedMass) + 
                       scale(mean.PlantHeight), family='poisson', data=modData.N.Tmax, na.action='na.fail')


plot(modSpa.N.Tmax)
plot(residuals(modSpa.N.Tmax), fitted(modSpa.N.Tmax))
hist(modSpa.N.Tmax$residuals, breaks=100)
summary(modSpa.N.Tmax)

# visreg(modSpa.N.Tmax,'base.Tmax', gg=T)+
#   ylab("n partial residuals")+
#   theme_classic()

# dr_modSpa.N.Tmax <- dredge(modSpa.N.Tmax)

# modData$resGLM.Spatial.N.Tmax<- residuals(modSpa.N.Tmax)
# worldmap + 
#   geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resGLM.Spatial.N.Tmax), fill=resGLM.Spatial.N.Tmax), shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   theme_classic()+
#   ggtitle('N.Tmax model residuals')



modSpa.N.Tmax.Visreg <- gam(crossings.max.Tmax ~ s(Latitude, Longitude, bs='sos') + 
                       poly(base.Tmax,2) + 
                       scale(sqrt_topo25) +
                       scale(sqrt_humanDensity) +
                       scale(dist_refugia) + 
                       scale(mean.SeedMass) + 
                       scale(mean.PlantHeight), family='poisson', data=modData.N.Tmax, na.action='na.fail')


visreg(modSpa.N.Tmax.Visreg,'base.Tmax', gg=T)+
  ylab("n partial residuals")+
  theme_classic()


# #### crossings.max Prec ####
# 
# ## Poisson Linear model
# mod.N.Prec <- glm(crossings.max.Prec ~ 
#                     scale(base.Prec) + 
#                     scale(sqrt_topo25) +
#                     scale(sqrt_humanDensity) +
#                     scale(dist_refugia) + 
#                     scale(mean.SeedMass) + 
#                     scale(mean.PlantHeight), family=poisson, data=modData, na.action='na.fail')
# 
# # plot(mod.N.Prec)
# 
# modData$resGLM.N.Prec<- residuals(mod.N.Prec)
# hist(modData$resGLM.N.Prec, breaks=100)
# 
# worldmap + 
#   geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resGLM.N.Prec), fill=resGLM.N.Prec), shape=21)+
#   scale_fill_gradient2()+
#   scale_size(guide=FALSE)+
#   theme_classic()
# 
# ggplot(modData, aes(topo25, resGLM.N.Prec))+
#   geom_point(shape=21)+
#   theme_classic()+
#   stat_smooth()
# 
# 
# ## Poisson Linear model + spatial smoothed
# modSpa.N.Prec <- gam(crossings.max.Prec ~ s(Latitude, Longitude, bs='sos') + base.Prec + sqrt(topo25) + sqrt(humanDensity) + scale(dist_refugia) + 
#                        mean.SeedMass + mean.PlantHeight, family='poisson', data=modData, na.action='na.fail')
# 
# plot(modSpa.N.Prec)
# plot(residuals(modSpa.N.Prec), fitted(modSpa.N.Prec))
# modData$resGLM.Spatial.N.Prec<- residuals(modSpa.N.Prec)
# hist(modData$resGLM.Spatial.N.Prec, breaks=100)
# 
# worldmap + 
#   geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resGLM.Spatial.N.Prec), fill=resGLM.Spatial.N.Prec), shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   theme_classic()+
#   ggtitle('N.Prec model residuals')
# 
# AIC(mod.N.Prec)
# AIC(modSpa.N.Prec)
# 
# 
# 
# ## final model 
# # Poisson Linear model + spatial smoothed
# modData.N.Prec <- modData
# 
# modSpa.N.Prec <- gam(crossings.max.Prec ~ s(Latitude, Longitude, bs='sos') + 
#                        scale(poly(base.Prec,2)) + 
#                        scale(sqrt_topo25) +
#                        scale(sqrt_humanDensity) +
#                        scale(dist_refugia) + 
#                        scale(mean.SeedMass) + 
#                        scale(mean.PlantHeight), family='poisson', data=modData.N.Prec, na.action='na.fail')
# 
# 
# plot(modSpa.N.Prec)
# plot(residuals(modSpa.N.Prec), fitted(modSpa.N.Prec))
# hist(modSpa.N.Prec$residuals, breaks=100)
# summary(modSpa.N.Prec)
# # visreg(modSpa.N.Prec)
# # dr_modSpa.N.Prec <- dredge(modSpa.N.Prec)
# gam.vcomp(modSpa.N.Prec,rescale=TRUE,conf.lev=.95)


#### deviation.absolute Tmax####
hist(modData$deviation.absolute.Tmax)
hist(modData$deviation.absolute.Prec)


hist(sqrt(modData$deviation.absolute.Tmax))
hist(sqrt(modData$deviation.absolute.Prec))

## Linear model
mod.devAbs.Tmax <- lm(deviation.absolute.Tmax ~ 
                        scale(base.Tmax) + 
                        scale(sqrt_topo25) +
                        scale(sqrt_humanDensity) +
                        scale(dist_refugia) + 
                        scale(mean.SeedMass) + 
                        scale(mean.PlantHeight), data=modData, na.action='na.fail')
# plot(mod.devAbs.Tmax)

modData$resLM.devAbs.Tmax<- residuals(mod.devAbs.Tmax)
hist(modData$resLM.devAbs.Tmax)

# [modData$resLM.devAbs.Tmax<10, ]
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resLM.devAbs.Tmax), fill=resLM.devAbs.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()

AIC(mod.devAbs.Tmax)

## Linear model + smoothed spatial coordinates
modSpa.devAbs.Tmax <- gam(deviation.absolute.Tmax ~ s(Latitude, Longitude, bs='sos') + base.Tmax + sqrt(topo25) + sqrt(humanDensity) + scale(dist_refugia) + 
                            mean.SeedMass + mean.PlantHeight, data=modData, na.action='na.fail')

plot(y=residuals(modSpa.devAbs.Tmax), fitted(modSpa.devAbs.Tmax))

modData$resSpa.devAbs.Tmax<- residuals(modSpa.devAbs.Tmax)
hist(modData$resSpa.devAbs.Tmax, breaks=100)

[modData$resSpa.devAbs.Tmax<5, ]
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resSpa.devAbs.Tmax), fill=resSpa.devAbs.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()

modData.devAbs.Tmax <- modData[modData$resSpa.devAbs.Tmax<5, ]

summary(modSpa.devAbs.Tmax)
AIC(modSpa.devAbs.Tmax)

## final model 
# Linear model + spatial smoothed
modData.devAbs.Tmax <- modData[modData$resSpa.devAbs.Tmax<5, ]
modSpa.devAbs.Tmax <- gam(deviation.absolute.Tmax  ~ s(Latitude, Longitude, bs='sos') + 
                       scale(poly(base.Tmax, 2)) + 
                       scale(sqrt_topo25) +
                       scale(sqrt_humanDensity) +
                       # scale(dist_refugia) + 
                       scale(mean.SeedMass) + 
                       scale(mean.PlantHeight),  data=modData.devAbs.Tmax, na.action='na.fail')


plot(modSpa.devAbs.Tmax)
plot(residuals(modSpa.devAbs.Tmax), fitted(modSpa.devAbs.Tmax))
hist(modSpa.devAbs.Tmax$residuals, breaks=100)
summary(modSpa.devAbs.Tmax)


modSpa.devAbs.Tmax.Visreg <- gam(deviation.absolute.Tmax  ~ s(Latitude, Longitude, bs='sos') + 
                            poly(base.Tmax, 2) + 
                            scale(sqrt_topo25) +
                            scale(sqrt_humanDensity) +
                            # scale(dist_refugia) + 
                            scale(mean.SeedMass) + 
                            scale(mean.PlantHeight),  data=modData.devAbs.Tmax, na.action='na.fail')

visreg(modSpa.devAbs.Tmax.Visreg,'base.Tmax', gg=T)+
  ylab("deviation absolute Tmax partial residuals")+
  theme_classic()


visreg(modSpa.devAbs.Tmax.Visreg,'mean.PlantHeight', gg=T)+
  ylab("deviation absolute Tmax partial residuals")+
  theme_classic()

# dr_modSpa.devAbs.Tmax <- dredge(modSpa.devAbs.Tmax)

modData.devAbs.Tmax$resSpa.devAbs.Tmax<- residuals(modSpa.devAbs.Tmax)

worldmap + 
  geom_point(data=modData.devAbs.Tmax, aes(y=Latitude, x=Longitude, size= abs(resSpa.devAbs.Tmax), fill=resSpa.devAbs.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()+
  ggtitle('Absolute Deviation Tmax model residuals')

# #### deviation.absolute Prec####
# hist(modData$deviation.absolute.Prec)
# hist(sqrt(modData$deviation.absolute.Prec))
# 
# plot(y=modData$deviation.absolute.Prec, modData$base.Prec)
# plot(y=modData$deviation.absolute.Tmax, modData$base.Tmax)
# 
# ## Linear model
# mod.devAbs.Prec <- lm(deviation.absolute.Prec ~ 
#                         scale(base.Prec) + 
#                         scale(sqrt_topo25) +
#                         scale(sqrt_humanDensity) +
#                         scale(dist_refugia) + 
#                         scale(mean.SeedMass) + 
#                         scale(mean.PlantHeight), data=modData, na.action='na.fail')
# # plot(mod.devAbs.Prec)
# 
# modData$resLM.devAbs.Prec<- residuals(mod.devAbs.Prec)
# hist(modData$resLM.devAbs.Prec)
# hist(modData$resLM.devAbs.Prec[modData$resLM.devAbs.Prec<400])
# 
# # [modData$resLM.devAbs.Prec<400, ]
# worldmap + 
#   geom_point(data=modData[modData$resLM.devAbs.Prec<400, ], aes(y=Latitude, x=Longitude, size= abs(resLM.devAbs.Prec), fill=resLM.devAbs.Prec), shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   theme_classic()
# 
# AIC(mod.devAbs.Prec)
# 
# ## Linear model + smoothed spatial coordinates
# modSpa.devAbs.Prec <- gam(deviation.absolute.Prec ~ s(Latitude, Longitude, bs='sos') + base.Prec + sqrt(topo25) + sqrt(humanDensity) + scale(dist_refugia) + 
#                             mean.SeedMass + mean.PlantHeight, data=modData, na.action='na.fail')
# 
# plot(y=residuals(modSpa.devAbs.Prec), fitted(modSpa.devAbs.Prec))
# 
# modData$resSpa.devAbs.Prec<- residuals(modSpa.devAbs.Prec)
# hist(modData$resSpa.devAbs.Prec, breaks=100)
# hist(modData$resSpa.devAbs.Prec[modData$resSpa.devAbs.Prec<200], breaks=100)
# 
# 
# worldmap + 
#   geom_point(data=modData[modData$resSpa.devAbs.Prec<200, ], aes(y=Latitude, x=Longitude, size= abs(resSpa.devAbs.Prec), fill=resSpa.devAbs.Prec), shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   theme_classic()+
#   ggtitle('Absolute Deviation Prec model residuals')
# 
# 
# summary(modSpa.devAbs.Prec)
# AIC(modSpa.devAbs.Prec)
# 
# 
# ## final model 
# # Linear model + spatial smoothed
# modData.devAbs.Prec <- modData[modData$resSpa.devAbs.Prec<200, ]
# 
# modSpa.devAbs.Prec <- gam(deviation.absolute.Prec  ~ s(Latitude, Longitude, bs='sos') + 
#                             scale(poly(base.Prec,2)) + 
#                             scale(sqrt_topo25) +
#                             scale(sqrt_humanDensity) +
#                             scale(dist_refugia) + 
#                             scale(mean.SeedMass) + 
#                             scale(mean.PlantHeight), data=modData.devAbs.Prec, na.action='na.fail')
# 
# 
# plot(modSpa.devAbs.Prec)
# plot(residuals(modSpa.devAbs.Prec), fitted(modSpa.devAbs.Prec))
# hist(modSpa.devAbs.Prec$residuals)
# summary(modSpa.devAbs.Prec)
# 
# visreg(modSpa.devAbs.Prec,'base.Prec', gg=T)+
#   ylab("deviation absolute Prec partial residuals")+
#   theme_classic()
# 
# # dr_modSpa.devAbs.Prec <- dredge(modSpa.devAbs.Prec)
# 
# modData.devAbs.Prec$resSpa.devAbs.Prec<- residuals(modSpa.devAbs.Prec)
# 
# worldmap + 
#   geom_point(data=modData.devAbs.Prec, aes(y=Latitude, x=Longitude, size= abs(resSpa.devAbs.Prec), fill=resSpa.devAbs.Prec), shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   theme_classic()

#### deviation.change Tmax ####
hist(modData$deviation.change.1Ky.Tmax)

## Linear model
mod.devChange.Tmax <- lm(deviation.change.1Ky.Tmax ~ 
                           scale(poly(base.Tmax,2)) + 
                           scale(sqrt_topo25) +
                           scale(sqrt_humanDensity) +
                           scale(dist_refugia) + 
                           scale(mean.SeedMass) + 
                           scale(mean.PlantHeight), data=modData, na.action='na.fail')

# plot(mod.devChange.Tmax)

modData$resLM.devChange.Tmax<- residuals(mod.devChange.Tmax)
hist(modData$resLM.devChange.Tmax, breaks=100)

# [modData$resLM.devChange.Tmax>-2, ]
worldmap + 
  geom_point(data=modData[modData$resLM.devChange.Tmax>-1.5, ], aes(y=Latitude, x=Longitude, size= abs(resLM.devChange.Tmax), fill=resLM.devChange.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()

AIC(mod.devChange.Tmax)

## Linear model + smoothed spatial coordinates
modSpa.devChange.Tmax <- gam(deviation.change.1Ky.Tmax ~ s(Latitude, Longitude, bs='sos') + base.Tmax + sqrt_topo25 + sqrt_humanDensity + scale(dist_refugia) + 
                               mean.SeedMass + mean.PlantHeight, data=modData, na.action='na.fail')
plot(y=residuals(modSpa.devChange.Tmax), fitted(modSpa.devChange.Tmax))

modData$resSpa.devChange.Tmax<- residuals(modSpa.devChange.Tmax)
hist(modData$resSpa.devChange.Tmax[modData$resSpa.devChange.Tmax>-1], breaks=100)

worldmap + 
  geom_point(data=modData[modData$resSpa.devChange.Tmax>-1, ], aes(y=Latitude, x=Longitude, size= abs(resSpa.devChange.Tmax), fill=resSpa.devChange.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()

AIC(modSpa.devChange.Tmax)


modData.devChange.Tmax <- modData[modData$resSpa.devChange.Tmax>-1, ] 


mod.devChange.Tmax <- lm(deviation.change.1Ky.Tmax ~ base.Tmax + sqrt(topo25) + sqrt(humanDensity) + scale(dist_refugia) + 
                               mean.SeedMass + mean.PlantHeight, data=modData.devChange.Tmax, na.action='na.fail')

modSpa.devChange.Tmax <- gam(deviation.change.1Ky.Tmax ~ s(Latitude, Longitude, bs='sos') + base.Tmax + sqrt(topo25) + sqrt(humanDensity) + scale(dist_refugia) + 
                               mean.SeedMass + mean.PlantHeight, data=modData.devChange.Tmax, na.action='na.fail')

hist(residuals(modSpa.devChange.Tmax), breaks=100)
summary(modSpa.devChange.Tmax)
# visreg(modSpa.devChange.Tmax)


## final model 
# Linear model + spatial smoothed
modSpa.devChange.Tmax <- gam(deviation.change.1Ky.Tmax  ~ s(Latitude, Longitude, bs='sos') + 
                            scale(poly(base.Tmax,2)) + 
                            scale(sqrt_topo25) +
                            scale(sqrt_humanDensity) +
                            # scale(dist_refugia) + 
                            scale(mean.SeedMass) + 
                            scale(mean.PlantHeight),  data=modData.devChange.Tmax, na.action='na.fail')


plot(modSpa.devChange.Tmax)
plot(residuals(modSpa.devChange.Tmax), fitted(modSpa.devChange.Tmax))
hist(modSpa.devChange.Tmax$residuals, breaks=100)
summary(modSpa.devChange.Tmax)

# dr_modSpa.devChange.Tmax <- dredge(modSpa.devChange.Tmax)

modSpa.devChange.Tmax.Visreg <- gam(deviation.change.1Ky.Tmax  ~ s(Latitude, Longitude, bs='sos') + 
                               poly(base.Tmax,2) + 
                               scale(sqrt_topo25) +
                               scale(sqrt_humanDensity) +
                               # scale(dist_refugia) + 
                               scale(mean.SeedMass) + 
                               scale(mean.PlantHeight),  data=modData.devChange.Tmax, na.action='na.fail')



visreg(modSpa.devChange.Tmax.Visreg,'base.Tmax', gg=T)+
  ylab("deviation change Tmax partial residuals")+
  theme_classic()+
  geom_hline(yintercept=0)

modData.devChange.Tmax$resSpa.devChange.Tmax<- residuals(modSpa.devChange.Tmax)

worldmap + 
  geom_point(data=modData.devChange.Tmax, aes(y=Latitude, x=Longitude, size= abs(resSpa.devChange.Tmax), fill=resSpa.devChange.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()+
  ggtitle('Deviation change Tmax model residuals')

# #### deviation.change Prec ####
# hist(modData$deviation.change.1Ky.Prec)
# 
# ## Linear model
# mod.devChange.Prec <- lm(deviation.change.1Ky.Prec ~
#                            scale(base.Prec) + 
#                            scale(sqrt_topo25) +
#                            scale(sqrt_humanDensity) +
#                            scale(dist_refugia) + 
#                            scale(mean.SeedMass) + 
#                            scale(mean.PlantHeight), data=modData, na.action='na.fail')
# 
# # plot(mod.devChange.Prec)
# 
# modData$resLM.devChange.Prec<- residuals(mod.devChange.Prec)
# hist(modData$resLM.devChange.Prec, breaks=100)
# hist(modData$resLM.devChange.Prec[modData$resLM.devChange.Prec>-40 & modData$resLM.devChange.Prec < 40], breaks=100)
# 
# # [modData$resLM.devChange.Prec>-40 & modData$resLM.devChange.Prec < 40, ]
# 
# worldmap + 
#   geom_point(data=modData[modData$resLM.devChange.Prec>-40 & modData$resLM.devChange.Prec < 40, ], aes(y=Latitude, x=Longitude, size= abs(resLM.devChange.Prec), fill=resLM.devChange.Prec), shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   theme_classic()
# 
# AIC(mod.devChange.Prec)
# 
# ## Linear model + smoothed spatial coordinates
# modSpa.devChange.Prec <- gam(deviation.change.1Ky.Prec ~ s(Latitude, Longitude, bs='sos') + base.Prec + sqrt(topo25) + sqrt(humanDensity) + scale(dist_refugia) + 
#                                mean.SeedMass + mean.PlantHeight, data=modData, na.action='na.fail')
# plot(y=residuals(modSpa.devChange.Prec), fitted(modSpa.devChange.Prec))
# 
# modData$resSpa.devChange.Prec<- residuals(modSpa.devChange.Prec)
# hist(modData$resSpa.devChange.Prec, breaks=100)
# hist(modData$resSpa.devChange.Prec[modData$resSpa.devChange.Prec>-30 & modData$resSpa.devChange.Prec < 30], breaks=100)
# 
# # [modData$resSpa.devChange.Prec>-30 & modData$resSpa.devChange.Prec < 30, ]
# worldmap + 
#   geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resSpa.devChange.Prec), fill=resSpa.devChange.Prec), shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   theme_classic()
# 
# AIC(modSpa.devChange.Prec)
# 
# 
# modData.devChange.Prec <- modData[modData$resSpa.devChange.Prec>-30 & modData$resSpa.devChange.Prec < 30, ]
# 
# ## final model 
# # Poisson Linear model + spatial smoothed
# modSpa.devChange.Prec <- gam(deviation.change.1Ky.Prec  ~ s(Latitude, Longitude, bs='sos') + 
#                                poly(base.Prec,2) + 
#                                scale(sqrt_topo25) +
#                                scale(sqrt_humanDensity) +
#                                scale(dist_refugia) + 
#                                scale(mean.SeedMass) + 
#                                scale(mean.PlantHeight),  data=modData.devChange.Prec, na.action='na.fail')
# 
# 
# plot(modSpa.devChange.Prec)
# plot(residuals(modSpa.devChange.Prec), fitted(modSpa.devChange.Prec))
# hist(modSpa.devChange.Prec$residuals, breaks=100)
# summary(modSpa.devChange.Prec)
# # visreg(modSpa.devChange.Prec)
# # dr_modSpa.devChange.Prec <- dredge(modSpa.devChange.Prec)
# 
# visreg(modSpa.devChange.Prec,'base.Prec', gg=T)+
#   ylab("deviation change Prec partial residuals")+
#   theme_classic()+
#   geom_hline(yintercept=0)
# 
# modData.devChange.Prec$resSpa.devChange.Prec<- residuals(modSpa.devChange.Prec)
# 
# # [modData$resSpa.devChange.Prec>-30 & modData$resSpa.devChange.Prec < 30, ]
# worldmap + 
#   geom_point(data=modData.devChange.Prec, aes(y=Latitude, x=Longitude, size= abs(resSpa.devChange.Prec), fill=resSpa.devChange.Prec), shape=21)+
#   scale_fill_gradient2(guide=FALSE)+
#   scale_size(guide=FALSE)+
#   theme_classic()+
#   ggtitle('Deviation change Prec model residuals')

#### summary of models ####

summary(modSpa.N.Tmax)
# summary(modSpa.N.Prec)
summary(modSpa.devAbs.Tmax)
# summary(modSpa.devAbs.Prec)
summary(modSpa.devChange.Tmax)
# summary(modSpa.devChange.Prec)

visreg(modSpa.N.Tmax)
# visreg(modSpa.N.Prec)
visreg(modSpa.devAbs.Tmax)
# visreg(modSpa.devAbs.Prec)
visreg(modSpa.devChange.Tmax)
# visreg(modSpa.devChange.Prec)



summary(mod.N.Tmax)
# summary(mod.N.Prec)
summary(mod.devAbs.Tmax)
# summary(mod.devAbs.Prec)
summary(mod.devChange.Tmax)
# summary(mod.devChange.Prec)

modPred_res <-as.data.frame(rbind(
                   round(summary(modSpa.N.Tmax)$p.table, digits=4)[-1,],
                   # round(summary(modSpa.N.Prec)$p.table, digits=4)[-1,],
                   round(summary(modSpa.devAbs.Tmax)$p.table, digits=4)[-1,],
                   # round(summary(modSpa.devAbs.Prec)$p.table, digits=4)[-1,], 
                   round(summary(modSpa.devChange.Tmax)$p.table, digits=4)[-1,]))
                   # round(summary(modSpa.devChange.Prec)$p.table, digits=4)[-1,]))

modPred_res$predictor<- rep(c("Baseline.Climate","Baseline.Climate2", "Topography", "Human.Density", "Seed.mass", "Plant.height"))
modPred_res$expl.var <- rep(c("max.state.number", "absolute.deviation","deviation.change"), each=6)
colnames(modPred_res)<- c("coefficient", "SE", "z_value", "P_value","predictor","expl.var")


modPred_res$coefficient<- as.numeric(as.character(modPred_res$coefficient))
modPred_res$SE<- as.numeric(as.character(modPred_res$SE))
modPred_res$sig <- " "
modPred_res$sig[modPred_res$P_value < 0.1] <- "."
modPred_res$sig[modPred_res$P_value < 0.05] <- "*"
modPred_res$sig[modPred_res$P_value < 0.01] <- "**"


ggplot(modPred_res, aes(x=predictor, y=coefficient, col=sig))+
  geom_pointrange(aes(ymin=coefficient-(1.96*SE), ymax=coefficient+(1.96*SE)))+
  geom_hline(yintercept = 0)+
  facet_wrap(~expl.var,nrow=1, scales='free_y')+
  scale_color_manual(values=c('grey80','grey60', 'grey30', 'black'))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60, hjust=1), axis.title.x=element_blank())

write.csv(modPred_res, file="results/modPred_results_Neotoma.csv")



modPred_res_LM <- as.data.frame(rbind(
  round(summary(mod.N.Tmax)$coeff, digits=4)[-1,],
  round(summary(mod.N.Prec)$coeff, digits=4)[-1,],
  round(summary(mod.devAbs.Tmax)$coeff, digits=4)[-1,],
  round(summary(mod.devAbs.Prec)$coeff, digits=4)[-1,], 
  round(summary(mod.devChange.Tmax)$coeff, digits=4)[-1,],
  round(summary(mod.devChange.Prec)$coeff, digits=4)[-1,]))

modPred_res_LM$predictor<- rep(c("Baseline.Climate", "Topography", "Human.Density", "Glaciated.LGM", "Seed.mass", "Plant.height"))
modPred_res_LM$expl.var <- rep(c("N.Tmax", "N.Prec", "devAbs.Tmax", "devAbs.Prec", "devChange.Tmax", "devChange.Prec"), each=6)
colnames(modPred_res_LM)<- c("coefficient", "SE", "z_value", "P_value","predictor","expl.var")


modPred_res_LM$coefficient<- as.numeric(as.character(modPred_res_LM$coefficient))
modPred_res_LM$SE<- as.numeric(as.character(modPred_res_LM$SE))

ggplot(modPred_res_LM, aes(x=predictor, y=coefficient, col=predictor))+
  geom_pointrange(aes(ymin=coefficient-(1.96*SE), ymax=coefficient+(1.96*SE)))+
  geom_hline(yintercept = 0)+
  facet_wrap(~expl.var,nrow=1, scales='free_y')+
  scale_color_brewer( palette="Dark2")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none")

write.csv(modPred_res, file="results/modPred_results.csv")

#### Variance explained by non-smooth parameters ####
modSpa.devAbs.Tmax2 <- gam(deviation.absolute.Tmax  ~ 
                            scale(poly(base.Tmax, 2)) + 
                            scale(sqrt_topo25) +
                            scale(sqrt_humanDensity) +
                            scale(mean.SeedMass) + 
                            scale(mean.PlantHeight),  data=modData.devAbs.Tmax, na.action='na.fail')

summary(modSpa.devAbs.Tmax)$dev.expl 
summary(modSpa.devAbs.Tmax2)$dev.expl 




modSpa.devChange.Tmax2 <- gam(deviation.change.1Ky.Tmax  ~ 
                             scale(poly(base.Tmax, 2)) + 
                             scale(sqrt_topo25) +
                             scale(sqrt_humanDensity) +
                             scale(mean.SeedMass) + 
                             scale(mean.PlantHeight),  data=modData.devAbs.Tmax, na.action='na.fail')

summary(modSpa.devChange.Tmax)$dev.expl 
summary(modSpa.devChange.Tmax2)$dev.expl 



modSpa.N.Tmax2 <- gam(crossings.max.Tmax  ~ 
                                scale(poly(base.Tmax, 2)) + 
                                scale(sqrt_topo25) +
                                scale(sqrt_humanDensity) +
                                scale(mean.SeedMass) + 
                                scale(mean.PlantHeight),  data=modData.devAbs.Tmax, na.action='na.fail')

summary(modSpa.N.Tmax)$dev.expl 
summary(modSpa.N.Tmax2)$dev.expl 

#### hierarchical variance partitioning ####
library(hier.part)

## N tmax
hp.N.Tmax <- hier.part(y=modData.N.Tmax$crossings.max.Tmax, xcan = modData.N.Tmax[, c("Latitude", 
                                                       "Longitude",
                                                       "base.Tmax",
                                                       "sqrt_topo25",
                                                       "sqrt_humanDensity",
                                                       "mean.SeedMass",
                                                       "mean.PlantHeight")], family = 'poisson', gof="Rsqu") 

hp.N.TmaxPlot <- hp.N.Tmax$IJ
hp.N.TmaxPlot$predictor <- c("Latitude", "Longitude","Baseline.Climate", "Topography", "Human.Density", "Seed.mass", "Plant.height")
hp.N.TmaxPlot$stat <- 'N'

plot_hp.N.Tmax<-
  ggplot(hp.N.TmaxPlot, aes(y=I*100, x=predictor))+
  geom_col()+
  ggtitle('N.Tmax')+
  ylab('% explained variance by independent effect')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())+
  ylim(0,60)
# 
# 
# ## N Prec
# hp.N.Prec <- hier.part(y=modData.N.Prec$crossings.max.Prec, xcan = modData.N.Prec[, c("Latitude", 
#                                                                                       "Longitude",
#                                                                                       "base.Prec",
#                                                                                       "sqrt_topo25",
#                                                                                       "sqrt_humanDensity",
#                                                                                       "was_glaciated",
#                                                                                       "mean.SeedMass",
#                                                                                       "mean.PlantHeight")], family = 'poisson', gof="Rsqu") 
# 
# hp.N.PrecPlot <- hp.N.Prec$IJ
# hp.N.PrecPlot$predictor <- c("Latitude", "Longitude","Baseline.Climate", "Topography", "Human.Density", "Glaciated.LGM", "Seed.mass", "Plant.height")
# 
# 
# plot_hp.N.Prec<-
#   ggplot(hp.N.PrecPlot, aes(y=I*100, x=predictor))+
#   geom_col()+
#   ggtitle('N.Prec')+
#   ylab('% explained variance by independent effect')+
#   theme_classic()+
#   theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())

## absolute deviation Tmax
hp.devAbs.Tmax <- hier.part(y=modData.devAbs.Tmax$deviation.absolute.Tmax, 
                            xcan = modData.devAbs.Tmax[, c("Latitude",
                                                           "Longitude",
                                                           "base.Tmax",
                                                           "sqrt_topo25",
                                                           "sqrt_humanDensity",
                                                           "mean.SeedMass",
                                                           "mean.PlantHeight")], gof="Rsqu") 

hp.devAbs.TmaxPlot <- hp.devAbs.Tmax$IJ
hp.devAbs.TmaxPlot$predictor <- c("Latitude", "Longitude","Baseline.Climate", "Topography", "Human.Density",  "Seed.mass", "Plant.height")
hp.devAbs.TmaxPlot$stat <- 'devAbs'

plot_hp.devAbs.Tmax<-
  ggplot(hp.devAbs.TmaxPlot, aes(y=I*100, x=predictor))+
  geom_col()+
  ggtitle('devAbs.Tmax')+
  ylab('% explained variance by independent effect')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())+
  ylim(0,60)

# 
# 
# ## absolute deviation Prec
# hp.devAbs.Prec <- hier.part(y=modData.devAbs.Prec$deviation.absolute.Prec, 
#                             xcan = modData.devAbs.Prec[, c("Latitude",
#                                                            "Longitude",
#                                                            "base.Prec",
#                                                            "sqrt_topo25",
#                                                            "sqrt_humanDensity",
#                                                            "was_glaciated",
#                                                            "mean.SeedMass",
#                                                            "mean.PlantHeight")], gof="Rsqu") 
# 
# hp.devAbs.PrecPlot <- hp.devAbs.Prec$IJ
# hp.devAbs.PrecPlot$predictor <- c("Latitude", "Longitude","Baseline.Climate", "Topography", "Human.Density", "Glaciated.LGM", "Seed.mass", "Plant.height")
# 
# plot_hp.devAbs.Prec<-
#   ggplot(hp.devAbs.PrecPlot, aes(y=I*100, x=predictor))+
#   geom_col()+
#   ggtitle('devAbs.Prec')+
#   ylab('% explained variance by independent effect')+
#   theme_classic()+
#   theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())

## deviation change Tmax
hp.devChange.Tmax <- hier.part(y=modData.devChange.Tmax$deviation.change.1Ky.Tmax, 
                            xcan = modData.devChange.Tmax[, c("Latitude",
                                                           "Longitude",
                                                           "base.Tmax",
                                                           "sqrt_topo25",
                                                           "sqrt_humanDensity",
                                                           "mean.SeedMass",
                                                           "mean.PlantHeight")], gof="Rsqu") 

hp.devChange.TmaxPlot <- hp.devChange.Tmax$IJ
hp.devChange.TmaxPlot$predictor <- c("Latitude", "Longitude","Baseline.Climate", "Topography", "Human.Density", "Seed.mass", "Plant.height")
hp.devChange.TmaxPlot$stat <- 'devChange'

plot_hp.devChange.Tmax<-
  ggplot(hp.devChange.TmaxPlot, aes(y=I*100, x=predictor))+
  geom_col()+
  ggtitle('devChange.Tmax')+
  ylab('% explained variance by independent effect')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())+
  ylim(0,60)


## deviation change Prec
# hp.devChange.Prec <- hier.part(y=modData.devChange.Prec$deviation.change.1Ky.Prec, 
#                             xcan = modData.devChange.Prec[, c("Latitude",
#                                                            "Longitude",
#                                                            "base.Prec",
#                                                            "sqrt_topo25",
#                                                            "sqrt_humanDensity",
#                                                            "was_glaciated",
#                                                            "mean.SeedMass",
#                                                            "mean.PlantHeight")], gof="Rsqu") 
# 
# 
# hp.devChange.PrecPlot <- hp.devChange.Prec$IJ
# hp.devChange.PrecPlot$predictor <- c("Latitude", "Longitude","Baseline.Climate", "Topography", "Human.Density", "Glaciated.LGM", "Seed.mass", "Plant.height")
# 
# plot_hp.devChange.Prec<-
#   ggplot(hp.devChange.PrecPlot, aes(y=I*100, x=predictor))+
#   geom_col()+
#   ggtitle('devChange.Prec')+
#   ylab('% explained variance by independent effect')+
#   theme_classic()+
#   theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())


pdf('results/figures/hier.part_predictors', height=5)
plot_hp.N.Tmax
# plot_hp.N.Prec
plot_hp.devAbs.Tmax
# plot_hp.devAbs.Prec
plot_hp.devChange.Tmax
# plot_hp.devChange.Prec
dev.off()



hpPlot <- rbind(hp.devAbs.TmaxPlot, hp.devChange.TmaxPlot, hp.N.TmaxPlot)
library(RColorBrewer)

ggplot(hpPlot, aes(y=I*100, x=stat, fill=predictor))+
  geom_col(position='stack')+
  ylab('% explained variance by independent effect')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=60, hjust=1), legend.position = "bottom",
         axis.title.y=element_blank())+
  coord_flip()+
  scale_fill_brewer(palette="Dark2", type='qual')
  
  
