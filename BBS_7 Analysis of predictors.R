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
load("results/statistics_sites_BBS.Rdata")
load("results/predictors_BBS.Rdata")

#### merge datas
modData <- na.omit(left_join(statistics_sites, predictors))
modData$sqrt_topo25 <- sqrt(modData$topo25)

#### investigate correlation between statistics ####
library("PerformanceAnalytics")
chart.Correlation(modData[, c(
  "crossings.mean", "crossings.max" , "deviation.base", 
  "deviation.absolute" , "deviation.change"   , "deviation.change.1Ky",
  "num.points.not.na"  , "lm.intercept",   
  "lm.slope" , "lm.r2"    , "lm.rho")], histogram=TRUE, pch=19)

chart.Correlation(modData[, c(
  "crossings.max" , "deviation.absolute" , 
 "deviation.change.1Ky",
  "num.points.not.na")], histogram=TRUE, pch=19)


#### investigate correlation between predictors ####
chart.Correlation(modData[, c("Latitude", 
                              "Longitude",
                              "base.Tmax",
                              # "change.Tmax",
                              "sqrt_topo25",
                              "HII",
                              "BodyMass", 
                              'Migration')], histogram=TRUE, pch=19)

#### distributions of predictors ####
baseTmax <-
  ggplot(modData, aes(x=base.Tmax))+
  geom_histogram(fill="red")+
  theme_classic()

topo <-
  ggplot(modData, aes(x=topo25))+
  geom_histogram(fill="brown")+
  theme_classic()+
  scale_x_continuous(trans='sqrt')

hii <-
  ggplot(modData, aes(x=HII))+
  geom_histogram()+
  theme_classic()
  
bodymass <-
  ggplot(modData, aes(x=BodyMass))+
  geom_histogram(fill="darkgreen")+
  theme_classic()

migration <-
  ggplot(modData, aes(x=Migration))+
  geom_histogram(fill="darkorange")+
  theme_classic()

pdf("results/distribution_predictors_BBS.pdf", width=10)
baseTmax
topo
hii
bodymass
migration
dev.off()

#### maps of predictors ####
library(viridis)
library(rje)
library(RColorBrewer)
library(gridExtra)
world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,50), xlim=c(-125, -65))+
  theme(legend.position="bottom")

baseTmax_map <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=base.Tmax),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="baseline max. temp.", x=-60,y=30, size=4)

topo_map <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=topo25),size=4, shape=21)+
  scale_fill_viridis(trans='sqrt')+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="sqrt. Topography", x=-70,y=30, size=4)

hii_map <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=HII),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="Human Influence Index", x=-70,y=30, size=4)

bodymass_map <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=BodyMass),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="mean body mass", x=-70,y=30, size=4)

migration_map <-
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=Migration),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")+
  annotate("text", label="% migrants", x=-70,y=30, size=4)

pdf("results/maps_predictors_BBS.pdf", width=10)
baseTmax
topo
hii
bodymass
migration
dev.off()

#### crossings.max Tmax ####
hist(modData$crossings.max, breaks=c(0.5,1.5,2.5,3.5,4.5))

## Poisson Linear model
mod.N.Tmax <- glm(crossings.max ~ 
                    scale(base.Tmax) + 
                    scale(sqrt_topo25) +
                    scale(HII) +
                    scale(BodyMass) + 
                    scale(Migration), family='poisson', data=modData)

plot(mod.N.Tmax)
modData$resGLM.N.Tmax<- residuals(mod.N.Tmax)
hist(modData$resGLM.N.Tmax, breaks=100)

worldmap +
  geom_point(data=modData, aes(y=Latitude, x=Longitude, size= resGLM.N.Tmax, fill=resGLM.N.Tmax), shape=21)+
  scale_fill_gradient2()+
  scale_size(guide=FALSE)+
  theme_classic()

ggplot(modData, aes(Latitude, resGLM.N.Tmax))+
  geom_point(shape=21)+
  theme_classic()+
  stat_smooth()

## Poisson Linear model + spatial smoothed
modSpa.N.Tmax <- gam(crossings.max ~ s(Latitude, Longitude, bs='sos') + 
                       scale(base.Tmax,2) + 
                       scale(sqrt_topo25) +
                       scale(HII) +
                       scale(BodyMass) + 
                       scale(Migration), family='poisson', data=modData, na.action='na.fail')

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
# modData.N.Tmax <- modData[modData$resGLM.Spatial.N.Tmax < 1.5, ]

modSpa.N.Tmax <- gam(crossings.max ~ s(Latitude, Longitude, bs='sos') + 
                       scale(poly(base.Tmax,2)) + 
                       scale(sqrt_topo25) +
                       scale(HII) +
                       scale(BodyMass) + 
                       scale(Migration), family='poisson', data=modData, na.action='na.fail')


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

#### deviation.absolute Tmax####
hist(modData$deviation.absolute)
hist(sqrt(modData$deviation.absolute))

## Linear model
mod.devAbs.Tmax <- lm(deviation.absolute ~   
                        scale(base.Tmax) + 
                        scale(sqrt_topo25) +
                        scale(HII) +
                        scale(BodyMass) + 
                        scale(Migration), data=modData, na.action='na.fail')
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
modSpa.devAbs.Tmax <- gam(deviation.absolute ~ s(Latitude, Longitude, bs='sos') +
                            scale(base.Tmax) + 
                            scale(sqrt_topo25) +
                            scale(HII) +
                            scale(BodyMass) + 
                            scale(Migration), data=modData, na.action='na.fail')

plot(y=residuals(modSpa.devAbs.Tmax), fitted(modSpa.devAbs.Tmax))

modData$resSpa.devAbs.Tmax<- residuals(modSpa.devAbs.Tmax)
hist(modData$resSpa.devAbs.Tmax, breaks=100)

[modData$resSpa.devAbs.Tmax<5, ]
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resSpa.devAbs.Tmax), fill=resSpa.devAbs.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()


summary(modSpa.devAbs.Tmax)
AIC(modSpa.devAbs.Tmax)

## final model 
# Linear model + spatial smoothed
modData.devAbs.Tmax <- modData[modData$resSpa.devAbs.Tmax > -2 & modData$resSpa.devAbs.Tmax<2, ]
modSpa.devAbs.Tmax <- gam(deviation.absolute  ~ s(Latitude, Longitude, bs='sos') +
                            scale(poly(base.Tmax,2)) + 
                            scale(sqrt_topo25) +
                            scale(HII) +
                            scale(BodyMass) + 
                            scale(Migration),  data=modData, na.action='na.fail')


plot(modSpa.devAbs.Tmax)
plot(residuals(modSpa.devAbs.Tmax), fitted(modSpa.devAbs.Tmax))
hist(modSpa.devAbs.Tmax$residuals, breaks=100)
summary(modSpa.devAbs.Tmax)

visreg(modSpa.devAbs.Tmax,'base.Tmax', gg=T)+
  ylab("deviation absolute Tmax partial residuals")+
  theme_classic()


# dr_modSpa.devAbs.Tmax <- dredge(modSpa.devAbs.Tmax)

modData$resSpa.devAbs.Tmax<- residuals(modSpa.devAbs.Tmax)

worldmap + 
  geom_point(data=modData.devAbs.Tmax, aes(y=Latitude, x=Longitude, size= abs(resSpa.devAbs.Tmax), fill=resSpa.devAbs.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()+
  ggtitle('Absolute Deviation Tmax model residuals')


#### deviation.change Tmax ####
hist(modData$deviation.change.1Ky)

## Linear model
mod.devChange.Tmax <- lm(deviation.change.1Ky ~ 
                           scale(base.Tmax) + 
                           scale(sqrt_topo25) +
                           scale(HII) +
                           scale(BodyMass) + 
                           scale(Migration), data=modData, na.action='na.fail')

plot(mod.devChange.Tmax)

modData$resLM.devChange.Tmax<- residuals(mod.devChange.Tmax)
hist(modData$resLM.devChange.Tmax, breaks=100)

# [modData$resLM.devChange.Tmax>-2, ]
worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resLM.devChange.Tmax), fill=resLM.devChange.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()

AIC(mod.devChange.Tmax)

## Linear model + smoothed spatial coordinates
modSpa.devChange.Tmax <- gam(deviation.change.1Ky ~ s(Latitude, Longitude, bs='sos') +
                               scale(base.Tmax) + 
                               scale(sqrt_topo25) +
                               scale(HII) +
                               scale(BodyMass) + 
                               scale(Migration),
                             data=modData, na.action='na.fail')
plot(y=residuals(modSpa.devChange.Tmax), fitted(modSpa.devChange.Tmax))

modData$resSpa.devChange.Tmax<- residuals(modSpa.devChange.Tmax)
hist(modData$resSpa.devChange.Tmax, breaks=100)

worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, size= abs(resSpa.devChange.Tmax), fill=resSpa.devChange.Tmax), shape=21)+
  scale_fill_gradient2(guide=FALSE)+
  scale_size(guide=FALSE)+
  theme_classic()

AIC(modSpa.devChange.Tmax)



## final model 
# Linear model + spatial smoothed
modSpa.devChange.Tmax <- gam(deviation.change.1Ky  ~ s(Latitude, Longitude, bs='sos') + 
                               scale(poly(base.Tmax,2)) + 
                               scale(sqrt_topo25) +
                               scale(HII) +
                               scale(BodyMass) + 
                               scale(Migration),  data=modData, na.action='na.fail')


plot(modSpa.devChange.Tmax)
plot(residuals(modSpa.devChange.Tmax), fitted(modSpa.devChange.Tmax))
hist(modSpa.devChange.Tmax$residuals, breaks=100)
summary(modSpa.devChange.Tmax)

# dr_modSpa.devChange.Tmax <- dredge(modSpa.devChange.Tmax)
visreg(modSpa.devChange.Tmax)

visreg(modSpa.devChange.Tmax,'base.Tmax', gg=T)+
  ylab("deviation change Tmax partial residuals")+
  theme_classic()+
  geom_hline(yintercept=0)

#### summary of models ####

summary(modSpa.N.Tmax)
summary(modSpa.devAbs.Tmax)
summary(modSpa.devChange.Tmax)

visreg(modSpa.N.Tmax)
visreg(modSpa.devAbs.Tmax)
visreg(modSpa.devChange.Tmax)

summary(mod.N.Tmax)
summary(mod.devAbs.Tmax)
summary(mod.devChange.Tmax)

modPred_res <-as.data.frame(rbind(
                   round(summary(modSpa.N.Tmax)$p.table, digits=4)[-1,],
                   round(summary(modSpa.devAbs.Tmax)$p.table, digits=4)[-1,], 
                   round(summary(modSpa.devChange.Tmax)$p.table, digits=4)[-1,]), row.names = 1:18)

modPred_res$predictor<- factor(rep(c("Baseline.Climate","Baseline.Climate2", "Topography", "Human.Influence", "Body.mass", "Migration"),3), 
                               levels=c("Baseline.Climate","Baseline.Climate2", "Topography", "Human.Influence", "Body.mass", "Migration"))
modPred_res$expl.var <- rep(c("max.state.number", "absolute.deviation","deviation.change"), each=6)
colnames(modPred_res)<- c("coefficient", "SE", "z_value", "P_value","predictor","expl.var")


# modPred_res$coefficient<- as.numeric(as.character(modPred_res$coefficient))
# modPred_res$SE<- as.numeric(as.character(modPred_res$SE))
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

levels(modPred_res_BBS$predictor)

write.csv(modPred_res, file="results/modPred_results_BBS.csv", row.names = F)


modPred_res_LM <- as.data.frame(rbind(
  round(summary(mod.N.Tmax)$coeff, digits=4)[-1,],
  round(summary(mod.devAbs.Tmax)$coeff, digits=4)[-1,], 
  round(summary(mod.devChange.Tmax)$coeff, digits=4)[-1,]))

modPred_res_LM$predictor<- rep(c("Baseline.Climate", "Topography", "Human.Influence", "Body.mass", "Migration"))
modPred_res_LM$expl.var <- rep(c("N", "devAbs","devChange"), each=5)
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



#### Variance explained by non-smooth parameters ####
modSpa.devAbs.Tmax2 <- gam(deviation.absolute  ~ 
                             scale(poly(base.Tmax,2)) + 
                             scale(sqrt_topo25) + 
                             scale(HII) + 
                             scale(BodyMass) + 
                             scale(Migration),  data=modData.devAbs.Tmax, na.action='na.fail')

summary(modSpa.devAbs.Tmax)$dev.expl 
summary(modSpa.devAbs.Tmax2)$dev.expl 




modSpa.devChange.Tmax2 <- gam(deviation.change.1Ky  ~ 
                                scale(poly(base.Tmax,2)) + 
                                scale(sqrt_topo25) + 
                                scale(HII) + 
                                scale(BodyMass) + 
                                scale(Migration),  data=modData.devAbs.Tmax, na.action='na.fail')

summary(modSpa.devChange.Tmax)$dev.expl 
summary(modSpa.devChange.Tmax2)$dev.expl 


modSpa.N.Tmax2 <- gam(crossings.max  ~ 
                        scale(poly(base.Tmax,2)) + 
                        scale(sqrt_topo25) + 
                        scale(HII) + 
                        scale(BodyMass) + 
                        scale(Migration),  data=modData.devAbs.Tmax, na.action='na.fail')

summary(modSpa.N.Tmax)$dev.expl 
summary(modSpa.N.Tmax2)$dev.expl 


#### hierarchical variance partitioning ####
library(hier.part)

## N tmax
hp.N.Tmax <- hier.part(y=modData$crossings.max, xcan = modData[, c("Latitude", 
                                                       "Longitude",
                                                       "base.Tmax",
                                                       "sqrt_topo25",
                                                       "HII",
                                                       "BodyMass",
                                                       "Migration")], family = 'poisson', gof="Rsqu") 

hp.N.TmaxPlot <- hp.N.Tmax$IJ
hp.N.TmaxPlot$predictor <- c("Latitude", 
                             "Longitude",
                             "base.Tmax",
                             "sqrt_topo25",
                             "HII",
                             "BodyMass",
                             "Migration")
hp.N.TmaxPlot$stat <- "N"


plot_hp.N.Tmax<-
  ggplot(hp.N.TmaxPlot, aes(y=I*100, x=predictor))+
  geom_col()+
  ggtitle('N.Tmax')+
  ylab('% explained variance by independent effect')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())

## absolute deviation Tmax
hp.devAbs.Tmax <- hier.part(y=modData$deviation.absolute, 
                            xcan = modData[, c("Latitude", 
                                               "Longitude",
                                               "base.Tmax",
                                               "sqrt_topo25",
                                               "HII",
                                               "BodyMass",
                                               "Migration")], gof="Rsqu") 

hp.devAbs.TmaxPlot <- hp.devAbs.Tmax$IJ
hp.devAbs.TmaxPlot$predictor <- c("Latitude", 
                                  "Longitude",
                                  "base.Tmax",
                                  "sqrt_topo25",
                                  "HII",
                                  "BodyMass",
                                  "Migration")

hp.devAbs.TmaxPlot$stat <- "devAbs"

plot_hp.devAbs.Tmax<-
  ggplot(hp.devAbs.TmaxPlot, aes(y=I*100, x=predictor))+
  geom_col()+
  ggtitle('devAbs.Tmax')+
  ylab('% explained variance by independent effect')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())


## deviation change Tmax
hp.devChange.Tmax <- hier.part(y=modData$deviation.change.1Ky, 
                            xcan = modData[, c("Latitude",
                                                           "Longitude",
                                                           "base.Tmax",
                                                           "sqrt_topo25",
                                                           "HII",
                                                           "BodyMass",
                                                           "Migration")], gof="Rsqu") 

hp.devChange.TmaxPlot <- hp.devChange.Tmax$IJ
hp.devChange.TmaxPlot$predictor <- c("Latitude", 
                                  "Longitude",
                                  "base.Tmax",
                                  "sqrt_topo25",
                                  "HII",
                                  "BodyMass",
                                  "Migration")
hp.devChange.TmaxPlot$stat <- "devChange"

plot_hp.devChange.Tmax<-
  ggplot(hp.devChange.TmaxPlot, aes(y=I*100, x=predictor))+
  geom_col()+
  ggtitle('devChange.Tmax')+
  ylab('% explained variance by independent effect')+
  theme_classic()+
  theme(axis.text.x=element_text(angle=60, hjust=1),legend.position = "none", axis.title.x=element_blank())

# 
# ## deviation change Prec
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
# 

pdf('results/hier.part_predictors_BBS', height=5)
plot_hp.N.Tmax
plot_hp.devAbs.Tmax
plot_hp.devChange.Tmax
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
