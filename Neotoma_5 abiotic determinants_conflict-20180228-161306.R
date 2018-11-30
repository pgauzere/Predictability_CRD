###############################################################################################
############# 5. Estimate external determinants of predictability in communities  #############
#################################################### Gauzere & Blonder ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS ")
library(raster)
library(ggplot2)
library(plyr)
library(dplyr)
# load("results/statistics_sites_neotoma_ext.drivers.Rdata")

load("results/statistics_sites_neotoma.Rdata")

# statistics_sites_external_scaled <- cbind(statistics_sites_scaled, statistics_sites_external[,25:33])
# save(statistics_sites_external_scaled, file="results/statistics_sites_neotoma_SCALED_ext.drivers.Rdata")

load("results/paleo_community_climato_distributions.Rdata")
load("results/sitesClim_buffer50km_neotoma.Rdata")
load("results/CommClimStat_smoothed.Rdata")

#### Baseline climate ####
# the baseline Tmax and annual precipitation at t1 for each site
abiotic <- 
CommClimStat %>% 
  filter(site %in% statistics_sites$site) %>%
  group_by(site, Latitude, Longitude, Altitude) %>% 
  arrange(desc(time)) %>%
  summarise(T1=time[1],
            base.Tmax=ClimTmax.ave[1], 
            base.Prec=ClimPrec.ave[1], 
            change.Tmax=lm(ClimTmax.ave ~ time$coefficients)[2]*(-1000),
            change.Prec=lm(ClimPrec.ave ~ time$coefficients)[2]*(-1000)) 



#### Topography & thermal heterogeneity ####
##topography
coorSites <- statistics_sites %>% select(site, Longitude, Latitude)
coordinates(coorSites) <- c('Longitude', 'Latitude')
#put on the same geo-referencial
proj4string(coorSites) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

usaAlt <- getData('alt',country="USA", mask=TRUE)[[1]]
canAlt <- getData('alt',country="CAN", mask=TRUE)
NorthAmAltRaster <- merge(usaAlt, canAlt)
proj4string(NorthAmAltRaster) <-  CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

plot(NorthAmAltRaster)   
plot(coorSites, add=T)

abiotic$topo25 <- raster::extract(NorthAmAltRaster,coorSites, buffer=25000, fun=max, na.rm=T) - # ~100km radius buffer
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

save(abiotic, file="results/abiotic_predictors.Rdata")




#### Human Settlement ####

##download HYDE 3.2 database from ftp server
plot(raster("popc_10000BC.asc", ymn=45, ymx=50, xmn=-100, xmx=105))
plot(raster("popd_10000BC.asc"))

times <- as.character( seq(from=1000, to=10000, by=1000)) # Define all needed time steps

ZIPfilenames <- as.list(paste(rep(c("lower","baseline", "upper"), time=10),"/zip/",times, "BC_pop.zip", sep="")) # create the list of filenames

names(ZIPfilenames)<- ZIPfilenames

#for all information on the HYDE 3.2 data used, see  ftp://ftp.pbl.nl/hyde/hyde3.2/readme_release_HYDE3.2.1.txt
hydeRasters <-
  llply(ZIPfilenames, function(x){    #iterate through filename list
    temp<- tempfile() #create temporary file
    download.file(url=as.character(paste("ftp://ftp.pbl.nl/hyde/hyde3.2/",x, sep="")), destfile=temp)    #download file from ftp using my filename list
    return(raster(unzip(temp)[[1]], crs=CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")))  #unzip and directly import raster to avoid stocking too much data
    unlink(temp)
  }, .progress="text")

hydeStack <- stack(hydeRasters)
plot(hydeStack)
writeRaster(hydeStack, filename="data/pop.count.Holocene.hyde3.21.rasters.tif", options="INTERLEAVE=BAND", overwrite=TRUE)




# load("results/statistics_sites_neotoma.Rdata")
coorSites <- statistics_sites %>% dplyr::select(site, Longitude, Latitude)
coordinates(coorSites) <- c('Longitude', 'Latitude')
#put on the same geo-referencial
proj4string(coorSites) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")


sitesHyde <- data.frame(coorSites$site,raster::extract(hydeStack, coorSites))
colnames(sitesHyde) <- gsub(".zip", "", colnames(sitesHyde) )
colnames(sitesHyde) <- gsub("BC_pop", "", colnames(sitesHyde) )

library(tidyr)
sitesHydeL <- sitesHyde %>% gather(key="name", "pop", 2:31)
sitesHydeL <- data.frame(sitesHydeL, ldply(strsplit(sitesHydeL$name,".",fixed=TRUE)))

sitesHydeL[,4] <- as.factor(sitesHydeL[,4])
sitesHydeL[,5] <- as.numeric(sitesHydeL[,5])
colnames(sitesHydeL)[c(1,4,5)] <- c("site", "scenario", "time")

save(sitesHydeL, file="results/sitesHydesvaluesLong.Rdata")

humanPop<-   sitesHydeL %>% group_by(site) %>% summarise(pop.ave=mean(pop,na.rm=T)) 

predictors <- abiotic

predictors$humanDensity  <- humanPop$pop.ave
predictors$humanPresence <- "human_presence"
predictors$humanPresence[predictors$humanDensity==0] <- "human_absence"

save(predictors, file="results/predictors.Rdata")

#### plots ####
#### plot distribution of statistics ####

library(rje)
library(RColorBrewer)
library(gridExtra)
hist_baseT <- 
  ggplot(predictors, aes(x=base.Tmax))+
  geom_histogram(fill="red")+
  theme_classic()

hist_baseP <- 
  ggplot(predictors, aes(x=base.Prec))+
  geom_histogram(fill="blue")+
  theme_classic()

hist_changeT <- 
  ggplot(predictors, aes(x=change.Tmax))+
  geom_histogram(fill="red")+
  theme_classic()+  
  geom_vline(xintercept = 0, linetype=3)

hist_changeP <- 
  ggplot(predictors, aes(x=change.Prec))+
  geom_histogram(fill="blue")+
  theme_classic()+
  geom_vline(xintercept = 0, linetype=3)

hist_topo<- 
  ggplot(predictors, aes(x=topo25))+
  geom_histogram(fill="blue")+
  theme_classic()


hist_heteroT<- 
  ggplot(predictors, aes(x=thermal.heterogeneity25))+
  geom_histogram(fill="blue")+
  theme_classic()

hist_heteroP<- 
  ggplot(predictors, aes(x=precip.heterogeneity25))+
  geom_histogram(fill="blue")+
  theme_classic()

hist_humanDensity<-
  ggplot(predictors, aes(x=humanDensity))+
  geom_histogram(fill="blue", bins=100)+
  scale_x_continuous(trans="log1p")+
  theme_classic()


pdf("results/distribution_predictors.pdf", width=16)
grid.arrange( hist_baseT,  hist_baseP, nrow=1)
grid.arrange( hist_changeT ,  hist_changeP, nrow=1)
grid.arrange( hist_heteroP,  hist_heteroT, nrow=1)
grid.arrange( hist_topo ,  hist_humanDensity, nrow=1)

dev.off()



#### map the statistics ####
library(rje)
library(RColorBrewer)
library(gridExtra)
world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,60), xlim=c(-115, -55))+
  theme(legend.position="bottom")

baseT <-
  worldmap + 
  geom_point(data=abiotic, aes(Longitude, Latitude, fill=base.Tmax),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")

baseP <-
  worldmap + 
  geom_point(data=abiotic, aes(Longitude, Latitude, fill=base.Prec),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")

changeT <-
  worldmap + 
  geom_point(data=abiotic, aes(Longitude, Latitude, fill=change.Tmax),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")

changeP <-
  worldmap + 
  geom_point(data=abiotic, aes(Longitude, Latitude, fill=change.Prec),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")

topo <-
  worldmap + 
  geom_point(data=abiotic, aes(Longitude, Latitude, fill=topo25),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")

heteroT <-
  worldmap + 
  geom_point(data=abiotic, aes(Longitude, Latitude, fill=thermal.heterogeneity25),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")


heteroP <-
  worldmap + 
  geom_point(data=abiotic, aes(Longitude, Latitude, fill=precip.heterogeneity25),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")

humanDensity <-
worldmap + 
    geom_point(data=predictors, aes(Longitude, Latitude, fill=humanDensity, col=humanPresence),size=3 , shape=21)+
    scale_fill_gradientn(trans="log1p",colours=cubeHelix(100))+
  scale_color_manual(values=c("red", 'black'))+
    theme_classic()+
    theme(legend.position="bottom")

  
  
pdf("results/maps_abiotic_pred.pdf", width=16)
grid.arrange( baseT,  baseP, nrow=1)
grid.arrange( changeT ,  changeP, nrow=1)
grid.arrange( heteroP ,  heteroT, nrow=1)
grid.arrange( topo ,  humanDensity, nrow=1)
dev.off()

library(viridis)
stats_predictors <- join(statistics_sites, predictors)
write.csv(stats_predictors, file="results/statistics_sites_neotoma_predictors.csv")

ggplot(stats_predictors, aes(x=base.Tmax, y=deviation.base.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  # geom_abline(intercept=0, slope=1)+
  scale_color_viridis(option="magma")


ggplot(stats_predictors, aes(x=topo25, y=deviation.absolute.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  # scale_x_continuous(trans='log1p')+
  # geom_abline(intercept=0, slope=1)+
  scale_color_viridis(option="magma")



ggplot(stats_predictors, aes(x=base.Prec, y=deviation.base.Prec))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  # geom_abline(intercept=0, slope=1)+
  scale_color_viridis(option="magma")


ggplot(stats_predictors, aes(x=base.Tmax, y=deviation.change.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  # geom_abline(intercept=0, slope=1)+
  scale_color_viridis(option="magma")

ggplot(stats_predictors, aes(x=base.Prec, y=deviation.change.Prec))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  # geom_abline(intercept=0, slope=1)+
  scale_color_viridis(option="magma")


ggplot(stats_predictors, aes(x=thermal.heterogeneity25, y=lm.rho.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  # geom_abline(intercept=0, slope=1)+
  scale_color_viridis(option="magma")

ggplot(stats_predictors, aes(x=thermal.heterogeneity25, y=precip.heterogeneity25))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  # geom_abline(intercept=0, slope=1)+
  scale_color_viridis(option="magma")



ggplot(stats_predictors, aes(x=change.Tmax, y=deviation.change.1Ky.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  geom_abline(intercept=0, slope=1)+
  scale_color_gradientn(colours=cubeHelix(100))

subtest <-test[test$num.points.not.na.Tmax>20,]

ggplot(test, aes(x=abs(change.Tmax), y=crossings.max.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  scale_color_gradientn(colours=cubeHelix(100))+
  theme_classic()
  geom_hline(yintercept=0)+
  # geom_abline(intercept=0, slope=1)+
  stat_smooth(method="lm")+

ggplot(test, aes(x=abs(change.Prec), y=crossings.mean.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  scale_color_gradientn(colours=cubeHelix(100))
  geom_hline(yintercept=0)+
  geom_abline(intercept=0, slope=1)+

ggplot(test, aes(x=base.Tmax, y=lm.rho.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  scale_color_gradientn(colours=cubeHelix(100))


ggplot(test, aes(x=change.Tmax, y=num.points.not.na.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  scale_color_gradientn(colours=cubeHelix(100))


ggplot(stats_predictors, aes(x=base.Tmax, y=change.Tmax))+
  geom_point(aes(size=num.points.not.na.Tmax,col=-num.points.not.na.Tmax))+
  theme_classic()+
  geom_hline(yintercept=0)+
  scale_color_gradientn(colours=cubeHelix(100))



library(mgcv)
mod  <- gamm(deviation.base.Tmax ~  s(Latitude, Longitude) + base.Tmax + change.Tmax + thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=stats_predictors)
mod2 <- gamm(deviation.change.1Ky.Tmax ~  s(Latitude, Longitude) + base.Tmax + change.Tmax + thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=stats_predictors)
mod3 <- gamm(lm.rho.Tmax  ~  s(Latitude, Longitude) + base.Tmax + change.Tmax + thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=stats_predictors)
mod4 <- gamm(crossings.mean.Tmax  ~  s(Latitude, Longitude) + base.Tmax + change.Tmax + thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=stats_predictors)

mod  <- gamm(deviation.base.Prec ~  s(Latitude, Longitude) + base.Prec + change.Prec + thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=stats_predictors)
mod2 <- gamm(deviation.change.1Ky.Prec ~  s(Latitude, Longitude) + base.Prec + change.Prec + thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=stats_predictors)
mod3 <- gamm(lm.rho.Prec  ~  s(Latitude, Longitude) + base.Prec + change.Prec + thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=stats_predictors)
mod4 <- gamm(crossings.mean.Prec  ~  s(Latitude, Longitude) + base.Prec + change.Prec + thermal.heterogeneity25 + log1p(humanDensity), random=list(site=~1), data=stats_predictors)



stsites <- na.omit(stats_predictors)
library(hier.part)

hpCrs <- hier.part(y=stsites$crossings.mean.Tmax, xcan=stsites[, c('Latitude', 'Longitude', 'base.Tmax', 'change.Tmax', 'thermal.heterogeneity25', 'humanDensity')], gof="Rsqu")
test_hpCrs <- rand.hp(y=stsites$crossings.mean.Tmax, xcan=stsites[, c('Latitude', 'Longitude', 'base.Tmax', 'change.Tmax', 'thermal.heterogeneity25', 'humanDensity')])


hpDeviationAbs <- hier.part(y=stsites$deviation.absolute.Tmax, xcan=stsites[, c('Latitude', 'Longitude', 'base.Tmax', 'change.Tmax', 'thermal.heterogeneity25', 'humanDensity')], gof="Rsqu")
hpDeviationChange <- hier.part(y=stsites$deviation.change.1Ky.Tmax, xcan=stsites[, c('Latitude', 'Longitude', 'base.Tmax', 'change.Tmax', 'thermal.heterogeneity25', 'humanDensity')], gof="Rsqu")
hpRho <- hier.part(y=stsites$lm.rho.Tmax, xcan=stsites[, c('Latitude', 'Longitude', 'base.Tmax', 'change.Tmax', 'thermal.heterogeneity25', 'humanDensity')], gof="Rsqu")

hpCrs <- hier.part(y=stsites$crossings.mean.Prec, xcan=stsites[, c('Latitude', 'Longitude', 'base.Prec', 'change.Prec', 'thermal.heterogeneity25', 'humanDensity')], gof="Rsqu")
hpDeviationAbs <- hier.part(y=stsites$deviation.absolute.Prec, xcan=stsites[, c('Latitude', 'Longitude', 'base.Prec', 'change.Prec', 'thermal.heterogeneity25', 'humanDensity')], gof="Rsqu")
hpDeviationChange <- hier.part(y=stsites$deviation.change.1Ky.Prec, xcan=stsites[, c('Latitude', 'Longitude', 'base.Prec', 'change.Prec', 'thermal.heterogeneity25', 'humanDensity')], gof="Rsqu")
hpRho <- hier.part(y=stsites$lm.rho.Prec, xcan=stsites[, c('Latitude', 'Longitude', 'base.Prec', 'change.Prec', 'thermal.heterogeneity25', 'humanDensity')], gof="Rsqu")

#### Ice sheet, Biomes and Vegetation at the LGM ####
load("results/statistics_sites_neotoma.Rdata")
library(maptools)
library(sp)

##read vegetation type raster at LGM from CLIMAP Project Members, 1984. The last interglacial ocean. Quat. Res. 21: 123-224. 
lgmVegRaster <-raster("data/LGM vegetation Raster.tiff")
proj4string(lgmVegRaster) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")
spplot(lgmVegRaster)

##read biome at LGM from Ray, N. and J. M. Adams. 2001. A GIS-based Vegetation Map of the World at the Last Glacial Maximum (25,000-15,000 BP). Internet Archaeology 11.
lgmVeg<- readShapeSpatial('data/LGM_vegetation/northam.shp',proj4string=CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0"))
lgmVeg@data <- join(lgmVeg@data, read.csv("data/LGM_vegetation/lgm_veg.txt", sep="\t", h=F, col.names= c("VEG_ID", "biome")), by="VEG_ID")
spplot(lgmVeg, "biome")

coorSites <- statistics_sites %>% dplyr::select(site, Longitude, Latitude)
coordinates(coorSites) <- c('Longitude', 'Latitude')
#put on the same geo-referencial
proj4string(coorSites) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

plot(lgmVegRaster, xlim=c(-140,50))
plot(lgmVeg, add=T)
plot(coorSites, add=T)  

## extract informations from shapeFile with function over from package sp
sitesBiomeType   <- data.frame(coorSites@data,coordinates(coorSites), over(coorSites,lgmVeg))
sitesBiomeType   <- join(sitesBiomeType, read.csv("data/LGM_vegetation/lgm_veg.txt", sep="\t", h=F, col.names= c("VEG_ID", "biome")), by="VEG_ID")

statistics_sites <- join(statistics_sites, sitesBiomeType) 

## extract informations from raster 
statistics_sites$vegLGM <- raster::extract(lgmVegRaster, coorSites,na.rm=T)
test <- join(data.frame("code2"=substr(statistics_sites$vegLGM, 1,1)),
read.csv("data/legend_2-LGM vegetation Raster.tiff.csv",h=F, col.names = c("code2", "veg.type")))
statistics_sites$vegLGM <- test$veg.type

##map and explore external determinants 
library(rje)
library(RColorBrewer)

world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,65), xlim=c(-120, -50))


## biome
biome <- 
worldmap + 
  geom_point(data=statistics_sites_external, aes(Longitude, Latitude, fill=biome),size=3, shape=21, col="white")+
  # scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  scale_fill_brewer(palette="Dark2")+
  theme(legend.position=c(0.85,0.2))

ggplot(statistics_sites, aes(x=biome, y=deviation.absolute.Prec))+
  # geom_jitter(aes(col=biome), shape=21)+
  geom_boxplot(aes(fill=biome),varwidth = T)+
  theme_classic()+
  scale_fill_brewer(palette="Dark2", guide=F)+
  coord_flip()

ggplot(statistics_sites, aes(x=biome, y=crossings.mean.Tmax))+
  # geom_jitter(aes(col=biome), shape=21)+
  geom_boxplot(aes(fill=biome),varwidth = T)+
  theme_classic()+
  scale_fill_brewer(palette="Dark2", guide=F)+
  coord_flip()

library(mgcv)
mod <- gamm(crossings.mean.Prec ~  s(Latitude, Longitude) + biome, random=list(site=~1), data=statistics_sites)
stsites <- na.omit(statistics_sites)
hpCrs <- hier.part(y=stsites$crossings.mean, xcan=stsites[, c("Longitude", "Latitude", "biome", "vegLGM", "thermal.gradient200")])
hpDeviationAbs <- hier.part(y=stsites$deviation.absolute, xcan=stsites[, c("Longitude", "Latitude", "biome", "vegLGM", "thermal.gradient200")])

## vegetation type
vegLGM <- 
  worldmap + 
  geom_point(data=statistics_sites_external, aes(Longitude, Latitude, fill=vegLGM),size=2, shape=21)+
  # scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  scale_fill_brewer(palette="Dark2")+
  theme(legend.position=c(0.85,0.2))

ggplot(statistics_sites, aes(x=vegLGM, y=deviation.absolute.Tmax))+
  # geom_jitter(aes(col=vegLGM), shape=21)+
  geom_boxplot(aes(fill=vegLGM))+
  theme_classic()+
  scale_fill_brewer(palette="Dark2", guide=F)+
  coord_flip()

ggplot(statistics_sites, aes(x=vegLGM, y=deviation.absolute.Prec))+
  # geom_jitter(aes(col=vegLGM), shape=21)+
  geom_boxplot(fill=NA, varwidth = T)+
  geom_violin(aes(fill=vegLGM), scale="count")+
  theme_classic()+
  scale_fill_brewer(palette="Dark2", guide=F)+
  coord_flip()




statistics_sites$ice.sheet <- "iced"
statistics_sites$ice.sheet[statistics_sites$VEG_ID != 26] <- "non-iced"

#### Topography & thermal gradient ####
##topography
country<- data.frame(ccodes())
NA_codes <- country[country$CONTINENT == "North America",]
# test <- alply(NA_codes,1, function(x){getData('alt',country=x, mask=TRUE)}, .progress = "text")

usaAlt <- getData('alt',country="USA", mask=TRUE)[[1]]
canAlt <- getData('alt',country="CAN", mask=TRUE)
NorthAmAltRaster <- (usaAlt, canAlt)
proj4string(NorthAmAltRaster) <-  CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

plot(coorSites)
plot(NorthAmAltRaster, add=T)   
plot(coorSites, add=T)

statistics_sites$topo50 <- raster::extract(NorthAmAltRaster,coorSites, buffer=50000, fun=max, na.rm=T) - # ~100km radius buffer
  raster::extract(NorthAmAltRaster,coorSites, buffer=50000, fun=min, na.rm=T)

topo<-
worldmap + 
  geom_point(data=statistics_sites, aes(Longitude, Latitude, fill=topo150),size=2, shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  # scale_fill_gradientn(colours)cubeHelix(100)+
  theme(legend.position=c(0.85,0.2))

ggplot(statistics_sites, aes(x=topo100, y=lm.r2))+
  geom_point()+
  stat_smooth(method='lm')+
  theme_classic()

## thermal gradient
worldclimRaster2.5 <- getData('worldclim', var='tmean', res=2.5) %>% raster::mean()

worldclimRasterTopLeft <-      getData('worldclim', var='tmean', res=0.5, lon=-100, lat=45) %>% raster::mean() 
worldclimRasterBottomRight <-  getData('worldclim', var='tmean', res=0.5, lon=-80, lat=25) %>% raster::mean() 
worldclimRasterBottomLeft <-   getData('worldclim', var='tmean', res=0.5, lon=-100, lat=25) %>% raster::mean() 
worldclimRastertopRightRight <-   getData('worldclim', var='tmean', res=0.5, lon=-50, lat=45) %>% raster::mean() 
worldclimRastertopTopRightRight <-   getData('worldclim', var='tmean', res=0.5, lon=-80, lat=65) %>% raster::mean() 

worldclimRaster <- mosaic(worldclimRasterTopRight , 
                          worldclimRasterTopLeft ,
                          worldclimRasterBottomRight , 
                          worldclimRasterBottomLeft,
                          worldclimRastertopTopRight,
                          worldclimRastertopTopRightRight, fun=mean)
proj4string(worldclimRaster) <-  CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

save(worldclimRaster, file="data/worldClimRaster_Tmean_res0.5.Rdata")

plot(coorSites)
plot(worldclimRaster)
     , add=T)
plot(coorSites, add=T)

statistics_sites$thermal.gradient100 <- raster::extract(worldclimRaster,coorSites, buffer=100000, fun=max, na.rm=T) -
  raster::extract(worldclimRaster,coorSites, buffer=100000, fun=min, na.rm=T)

statistics_sites$thermal.gradient200 <- raster::extract(worldclimRaster,coorSites, buffer=200000, fun=max, na.rm=T) -
  raster::extract(worldclimRaster,coorSites, buffer=200000, fun=min, na.rm=T)

worldmap + 
  geom_point(data=statistics_sites, aes(Longitude, Latitude, fill=thermal.gradient200),size=2, shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  # scale_fill_gradientn(colours)cubeHelix(100)+
  theme(legend.position=c(0.85,0.2))

ggplot(statistics_sites, aes(x=thermal.gradient200, y=crossings.mean))+
  geom_point()+
  stat_smooth(method='lm')+
  theme_classic()

#### Human Settlement ####

##download HYDE 3.2 database from ftp server
plot(raster("popc_10000BC.asc", ymn=45, ymx=50, xmn=-100, xmx=105))
plot(raster("popd_10000BC.asc"))

times <- as.character( seq(from=1000, to=10000, by=1000))
ZIPfilenames <- as.list(paste(rep(c("lower","baseline", "upper"), time=10),"/zip/",times, "BC_pop.zip", sep=""))
names(ZIPfilenames)<- ZIPfilenames
#for all information on the HYDE 3.2 data used, see  ftp://ftp.pbl.nl/hyde/hyde3.2/readme_release_HYDE3.2.1.txt
hydeRasters <- 
  llply(ZIPfilenames, function(x){
  temp<- tempfile()
  download.file(url=as.character(paste("ftp://ftp.pbl.nl/hyde/hyde3.2/",x, sep="")), destfile=temp)
  return(raster(unzip(temp)[[1]], crs=CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")))
  unlink(temp)
  }, .progress="text")


save(hydeRasters, file="data/pop.count.Holocene.hyde3.21.rasters.Rdata")


load("data/pop.count.Holocene.hyde3.21.rasters.Rdata")
hydeStack <- stack(hydeRasters)
plot(hydeStack)

# load("results/statistics_sites_neotoma.Rdata")
coorSites <- statistics_sites %>% dplyr::select(site, Longitude, Latitude)
coordinates(coorSites) <- c('Longitude', 'Latitude')
#put on the same geo-referencial
proj4string(coorSites) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")


sitesHyde <- data.frame(coorSites$site,raster::extract(hydeStack, coorSites))
colnames(sitesHyde) <- gsub(".zip", "", colnames(sitesHyde) )
colnames(sitesHyde) <- gsub("BC_pop", "", colnames(sitesHyde) )

library(tidyr)
sitesHydeL <- sitesHyde %>% gather(key="name", "pop", 2:31)
sitesHydeL <- data.frame(sitesHydeL, ldply(strsplit(sitesHydeL$name,".",fixed=TRUE)))

sitesHydeL[,4] <- as.factor(sitesHydeL[,4])
sitesHydeL[,5] <- as.numeric(sitesHydeL[,5])
colnames(sitesHydeL)[c(1,4,5)] <- c("site", "scenario", "time")

save(sitesHydeL, file="results/sitesHydesvaluesLong.Rdata")

humanPop <- sitesHydeL %>% group_by(site, scenario) %>% summarise(pop.ave=mean(pop,na.rm=T)) %>% spread(scenario, pop.ave)

statistics_sites_external<- join(statistics_sites, humanPop)
statistics_sites_external$humanDensity <- statistics_sites_external$baseline/85 #because human pop per 85km^2
save(statistics_sites_external[,-c(31,32,33)], file="results/statistics_sites_neotoma_ext.drivers.Rdata")


human<- 
  worldmap + 
  geom_point(data=statistics_sites_external, aes(Longitude, Latitude, fill=humanDensity),size=2, shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100), trans="sqrt")+
  theme_classic()+
  # scale_fill_gradientn(colours)cubeHelix(100)+
  theme(legend.position=c(0.85,0.2))
