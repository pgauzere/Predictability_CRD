###################################################
############# 7. Analysis of Predictors ###########
############################## Gauzere & Blonder ##


setwd("~/PIERRE/TRANSIENT DYNAMICS ")
library(ggplot2)
library(plyr)
library(dplyr)
library(mgcv)
library(hier.part)
library(colorRamps)
library(rje)
# load("results/statistics_sites_neotoma_SCALED_ext.drivers.Rdata")


#### test relationships between external drivers and statistics ####
load("results/statistics_sites_neotoma.Rdata")
load("results/predictors.Rdata")
load("results/biotic_predictors_CWMs.Rdata")

temp<-join(statistics_sites, predictors)
stpred <- left_join(temp, CWMs)
modData<- na.omit(stpred)

library(hier.part)
hpMod.N.Tmax <- hier.part(y=stpred$crossings.max.Tmax, 
                          xcan=stpred[,c("Latitude", "Longitude","base.Tmax","change.Tmax", "thermal.heterogeneity25","humanPresence", "mean.seed.mass", "mean.whole.plant.height")],
                          gof="logLik")

hpMod.absDev.Tmax <- hier.part(y=stpred$deviation.absolute.Tmax, 
                               xcan=stpred[,c("Latitude", "Longitude","base.Tmax","change.Tmax", "thermal.heterogeneity25","humanPresence", "mean.seed.mass", "mean.whole.plant.height")],
                               gof="RMSPE")

hpMod.DevChange.Tmax <- hier.part(y=stpred$deviation.change.1Ky.Tmax, 
                                  xcan=stpred[,c("Latitude", "Longitude","base.Tmax","change.Tmax", "thermal.heterogeneity25","humanPresence", "mean.seed.mass", "mean.whole.plant.height")],
                                  gof="RMSPE")

TesthpMod <- rand.hp(y=modData$crossings.max.Tmax, 
                     xcan=modData[,c("Latitude", "Longitude","base.Tmax", "change.Tmax", "thermal.heterogeneity25","humanPresence")], 
                     gof="Rsqu")

hpPlot <- hpMod.N.Tmax$IJ %>% mutate(predictors = rownames(.)) %>% gather("nature","explained_variance", 1:2)
hpMod.N.Tmax$IJ$predictor <- rownames(hpPlot)
ggplot(hpMod.N.Tmax$IJ, aes(x=predictor, y=Total))+
  geom_bar(stat="identity")

pdf("results/hier.part.Externals.pdf", width=12)
hpMod.N.Tmax <- hier.part(y=modData$crossings.mean.Tmax, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
hpMod.N.Prec <- hier.part(y=modData$crossings.mean.Prec, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
hpMod.AbsDev.Tmax <- hier.part(y=modData$deviation.absolute.Tmax, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
hpMod.AbsDev.Prec <- hier.part(y=modData$deviation.absolute.Prec, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
dev.off()



world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,65), xlim=c(-130, -50))+
  theme(legend.position="bottom")


stSites <- statistics_sites_external_scaled

pairs(stSites[,c(4,5,7,9,11,13,16,18,22,26,27,30,31)])

stSites$biome <- droplevels(stSites$biome)
ggpairs(stSites, columns=c(4,5,7,9,11,13,16,18,22,26,27,30,31))

,
        diag=list(continuous="density",   discrete="boxplot"), axisLabels="show")



modData <- na.omit(stpred)
mod <- lm(log(crossings.max.Tmax) ~  Latitude + Longitude + base.Tmax +  change.Tmax +  thermal.heterogeneity25 + humanPresence + mean.seed.mass, data=modData)
modData$res <- residuals(mod)
modData$fit <- fitted(mod)
plot(mod)

summary(mod)

modGamm <- gamm(log(crossings.max.Tmax) ~  s(Latitude, Longitude) + base.Tmax +  change.Tmax +  thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=modData)
modData$resGamm <- residuals(modGamm$lme)
modData$fitGamm <- fitted(modGamm$lme)

summary(modGamm$lme)

ggplot(modData, aes(y=fitGamm, x=resGamm, col=res))+
  geom_point(size=4)+
  theme_bw()+
  scale_color_gradientn(colours=matlab.like2(50))


ggplot(modData, aes(Longitude, Latitude, col=resGamm))+
  geom_point(size=4)+
  theme_bw()+
  scale_color_gradientn(colours=matlab.like2(50))

ggplot(modData, aes(y=log1p(baseline), x=res, col=res))+
  geom_point(size=4)+
  theme_bw()+
  scale_color_gradientn(colours=matlab.like2(50))

ggplot(modData, aes(x=vegLGM, y=resGamm))+
  geom_boxplot(varwidth = F)+
  theme_bw()


modGamm <- gamm(deviation.absolute.Prec ~  s(Latitude, Longitude) + biome, random=list(site=~1), data=modData)
summary(modGamm$lme)

hpMod.N.Tmax <- hier.part(y=modData$crossings.mean.Tmax, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
TesthpMod <- rand.hp(y=modData$crossings.mean.Tmax, xcan=modData[,c("Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")

hpPlot <- hpMod$IJ %>% mutate(predictors = rownames(.)) %>% gather("nature","explained_variance", 1:2)
hpMod$IJ$predictor <- rownames(hpMod)
ggplot(hpMod$IJ, aes(x=predictor, y=Total))+
  geom_bar(stat="identity")

pdf("results/hier.part.Externals.pdf", width=12)
hpMod.N.Tmax <- hier.part(y=modData$crossings.mean.Tmax, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
hpMod.N.Prec <- hier.part(y=modData$crossings.mean.Prec, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
hpMod.AbsDev.Tmax <- hier.part(y=modData$deviation.absolute.Tmax, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
hpMod.AbsDev.Prec <- hier.part(y=modData$deviation.absolute.Prec, xcan=modData[,c("num.points.not.na.Tmax","Latitude", "Longitude","topo150","biome","vegLGM", "baseline")], gof="Rsqu")
dev.off()

ggplot(stpred, aes(x=humanDensity, y=crossings.max.Tmax))+
  geom_point()+
  # geom_point(aes(y=crossings.mean.Prec), col="blue")+
  # stat_smooth(method="lm", formula= y~ x + x^2 )+
  # stat_smooth(aes(y=crossings.mean.Prec),method="lm", col='blue')+
  theme_bw()+
  scale_y_log10()

ggplot(stSites, aes(x=biome, y=crossings.mean.Tmax))+
  geom_boxplot(varwidth=T)+
  # geom_point(aes(y=crossings.mean.Prec), col="blue")+
  # stat_smooth(method="lm", formula= y~ x + x^2 )+
  # stat_smooth(aes(y=crossings.mean.Prec),method="lm", col='blue')+
  theme_bw()

ggplot(stSites, aes(x=num.points.not.na.Tmax, y=crossings.max.Prec))+
  geom_point(col="red")+
  # geom_point(aes(y=crossings.mean.Prec), col="blue")+
  stat_smooth(method="lm", )+
  # stat_smooth(aes(y=crossings.mean.Prec),method="lm", col='blue')+
  theme_bw()
