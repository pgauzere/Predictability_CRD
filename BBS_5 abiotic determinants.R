###############################################################################################
############# 5. Estimate external determinants of predictability in communities  #############
#################################################### Gauzere & Blonder ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS")
library(raster)
library(tidyverse)
load("results/statistics_sites_BBS.Rdata")

# statistics_sites_external_scaled <- cbind(statistics_sites_scaled, statistics_sites_external[,25:33])
# save(statistics_sites_external_scaled, file="results/statistics_sites_neotoma_SCALED_ext.drivers.Rdata")

load("results/CommClimStat2_smoothed_BBS.Rdata")

#### Baseline climate ####
# the baseline Tmax and annual precipitation at t1 for each site
abiotic <- 
CommClimStat2 %>% 
  group_by(site, Lati, Longi) %>% 
  arrange(time) %>%
  summarise(T1=time[1],
            base.Tmax=ClimTmax.ave[1], 
            change.Tmax=lm(ClimTmax.ave ~ time)$coefficients[2]) 



#### Topography & thermal heterogeneity ####
##topography
coorSites <- statistics_sites %>% select(site, Longi, Lati)
coordinates(coorSites) <- c('Longi', 'Lati')
#put on the same geo-referencial
proj4string(coorSites) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

usaAlt <- getData('alt',country="USA", mask=TRUE)[[1]]
canAlt <- getData('alt',country="CAN", mask=TRUE)
NorthAmAltRaster <- merge(usaAlt, canAlt)
# proj4string(NorthAmAltRaster) <-  CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

plot(NorthAmAltRaster)   
plot(coorSites, add=T)

abiotic$topo25 <- raster::extract(NorthAmAltRaster,coorSites, buffer=25000, fun=max, na.rm=T) -  
  raster::extract(NorthAmAltRaster,coorSites, buffer=25000, fun=min, na.rm=T)

## thermal heterogeneity
# worldclimRaster2.5             <-getData('worldclim', var='tmean', res=2.5) %>% raster::mean()
worldclimRasterMidRight        <-getData('worldclim', var='tmean', res=0.5, lon=-80, lat=45) %>% raster::mean() 
worldclimRasterMidLeft         <-getData('worldclim', var='tmean', res=0.5, lon=-100, lat=45) %>% raster::mean() 
worldclimRasterBottomRight     <-getData('worldclim', var='tmean', res=0.5, lon=-80, lat=25) %>% raster::mean() 
worldclimRasterBottomLeft      <-getData('worldclim', var='tmean', res=0.5, lon=-100, lat=25) %>% raster::mean() 
worldclimRasterRight           <-getData('worldclim', var='tmean', res=0.5, lon=-50, lat=45) %>% raster::mean() 
worldclimRastertopRight        <-getData('worldclim', var='tmean', res=0.5, lon=-80, lat=65) %>% raster::mean() 

worldclimRaster <- mosaic(worldclimRasterMidRight ,
                          worldclimRasterMidLeft ,
                          worldclimRasterBottomRight , 
                          worldclimRasterBottomLeft,
                          worldclimRasterRight,
                          worldclimRastertopRight, fun=mean)
proj4string(worldclimRaster) <-  CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

writeRaster(worldclimRaster, filename="data/worldClimRaster_Tmean_res0.5.tif", options="INTERLEAVE=BAND", overwrite=TRUE)


plot(coorSites)
plot(worldclimRaster, add=T)
plot(coorSites, add=T)
abline(h=c(25,45,60))
abline(v=c(-100,-80,-60))

abiotic$thermal.heterogeneity25 <- raster::extract(worldclimRaster,coorSites, buffer=25000, fun=sd, na.rm=T)
abiotic$thermal.heterogeneity10 <- raster::extract(worldclimRaster,coorSites, buffer=10000, fun=sd, na.rm=T)

save(abiotic, file="results/abiotic_predictors.Rdata")


## precipitation  heterogeneity
# worldclimRaster2.5             <-getData('worldclim', var='prec', res=2.5) %>% raster::mean()
worldclimRasterMidRight        <-getData('worldclim', var='prec', res=0.5, lon=-80, lat=45) %>% raster::mean() 
worldclimRasterMidLeft         <-getData('worldclim', var='prec', res=0.5, lon=-100, lat=45) %>% raster::mean() 
worldclimRasterBottomRight     <-getData('worldclim', var='prec', res=0.5, lon=-80, lat=25) %>% raster::mean() 
worldclimRasterBottomLeft      <-getData('worldclim', var='prec', res=0.5, lon=-100, lat=25) %>% raster::mean() 
worldclimRasterRight           <-getData('worldclim', var='prec', res=0.5, lon=-50, lat=45) %>% raster::mean() 
worldclimRastertopRight        <-getData('worldclim', var='prec', res=0.5, lon=-80, lat=65) %>% raster::mean() 

worldclimRaster_prec <- mosaic(worldclimRasterMidRight ,
                          worldclimRasterMidLeft ,
                          worldclimRasterBottomRight , 
                          worldclimRasterBottomLeft,
                          worldclimRasterRight,
                          worldclimRastertopRight, fun=mean)
proj4string(worldclimRaster_prec) <-  CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

writeRaster(worldclimRaster_prec, filename="data/worldClimRaster_Prec_res0.5.tif", options="INTERLEAVE=BAND", overwrite=TRUE)


plot(coorSites)
plot(worldclimRaster_prec, add=T)
plot(coorSites, add=T)
abline(h=c(25,45,60))
abline(v=c(-100,-80,-60))

abiotic$precip.heterogeneity25 <- raster::extract(worldclimRaster_prec,coorSites, buffer=25000, fun=sd, na.rm=T)
abiotic$precip.heterogeneity10 <- raster::extract(worldclimRaster_prec,coorSites, buffer=10000, fun=sd, na.rm=T)

# save(abiotic, file="results/abiotic_predictors.Rdata")




#### Human Influence Index ####
## for data download, see http://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-influence-index-geographic/data-download
hii<- raster('data/hii_n_america_grid/hii_n_amer/w001001.adf')
plot(hii, xlim=c(-130,-60), ylim=c(25,55))
plot(coorSites, add=T)
# proj4string(hii)
# proj4string(coorSites)
# #put on the same geo-referencial
# proj4string(coorSites) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

abiotic$HII <- raster::extract(hii, coorSites, buffer=25000, fun=mean, na.rm=T)

predictors<- abiotic

colnames(predictors)[2:4] <- c('Latitude', 'Longitude', 'year')

save(predictors, file="results/predictors_BBS.Rdata")

#### plots ####
#### plot distribution of abiotic predictors ####
load("results/predictors.Rdata")
library(rje)
library(RColorBrewer)
library(gridExtra)
hist_baseT <- 
  ggplot(predictors, aes(x=base.Tmax))+
  geom_histogram(fill="red")+
  theme_classic()

hist_changeT <- 
  ggplot(predictors, aes(x=change.Tmax))+
  geom_histogram(fill="red")+
  theme_classic()+  
  geom_vline(xintercept = 0, linetype=3)

hist_topo<- 
  ggplot(predictors, aes(x=topo25))+
  geom_histogram(fill="blue")+
  theme_classic()

hist_HII<-
  ggplot(predictors, aes(x=HII))+
  geom_histogram(fill="blue", bins=100)+
  scale_x_continuous(trans="log1p")+
  theme_classic()


#### map the predictors ####
library(rje)
library(RColorBrewer)
library(gridExtra)
world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,50), xlim=c(-125, -65))+
  theme(legend.position="bottom")

baseT <-
  worldmap + 
  geom_point(data=abiotic, aes(Longi, Lati, fill=base.Tmax),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")

changeT <-
  worldmap + 
  geom_point(data=abiotic, aes(Longi, Lati, fill=change.Tmax),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")

topo <-
  worldmap + 
  geom_point(data=abiotic, aes(Longi, Lati, fill=topo25),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100), trans='sqrt')+
  theme_classic()+
  theme(legend.position="bottom")

hii <-
  worldmap + 
    geom_point(data=predictors, aes(Longi, Lati, fill=HII),size=3 , shape=21)+
    scale_fill_gradientn(colours=cubeHelix(100))+
    theme_classic()+
    theme(legend.position="bottom")


pdf("results/maps_abiotic_pred.pdf", width=16)
plot( baseT)
plot( changeT)
plot( topo)
plot( humanDensity)
dev.off()


stats_predictors <- left_join(statistics_sites, predictors)
write.csv(stats_predictors, file="results/statistics_sites_neotoma_predictors.csv")

library(viridis)
ggplot(stats_predictors, aes(x=base.Tmax, y=deviation.absolute))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=0)

ggplot(stats_predictors, aes(x=topo25, y=deviation.absolute))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=0)+
  scale_x_continuous(trans='sqrt')+
  stat_smooth(method='lm', formula = y~poly(x,2))+
  coord_cartesian(ylim=c(0,8))
summary(lm(deviation.absolute ~ poly(sqrt(topo25),2), data=stats_predictors))

ggplot(stats_predictors, aes(x=HII, y=deviation.absolute))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=0)
summary(lm(deviation.absolute~ HII, data=stats_predictors))

##
ggplot(stats_predictors, aes(x=base.Tmax, y=deviation.change))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=0)+
  stat_smooth(method='lm', formula=y~poly(x,2))
summary(lm(deviation.change ~ poly(base.Tmax,2), data=stats_predictors))


ggplot(stats_predictors, aes(x=topo25, y=deviation.change))+
  geom_point()+
  theme_classic()+
  geom_hline(yintercept=0)+
  scale_x_continuous(trans='sqrt')+
  stat_smooth(method='lm')
summary(lm(deviation.change ~ topo25, data=stats_predictors))


ggplot(stats_predictors, aes(x=HII, y=deviation.change))+
  geom_point()+
  theme_classic()+
  stat_smooth(method='lm')+
  # stat_quantile()+
  geom_hline(yintercept=0)

summary(lm(deviation.change ~ HII, data=stats_predictors))

##
ggplot(stats_predictors, aes(y=base.Tmax, x=as.factor(crossings.max)))+
  geom_boxplot()+
  theme_classic()+
  coord_flip()

ggplot(stats_predictors, aes(y=topo25, x=as.factor(crossings.max)))+
  geom_boxplot()+
  theme_classic()+
  scale_y_continuous(trans='sqrt')+
  coord_flip()

ggplot(stats_predictors, aes(y=HII, x=as.factor(crossings.max)))+
  geom_boxplot()+
  theme_classic()+
  coord_flip()

