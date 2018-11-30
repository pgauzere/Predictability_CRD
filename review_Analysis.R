#### correlation between Tmax and Tmin for modern cliamte change ####
setwd("~/PIERRE/TRANSIENT DYNAMICS")

library(prism)
library(tidyverse)
library(sp)

options(prism.path = "~/DATA_prism/")
prismFiles <- ls_prism_data(name = TRUE)
prismFiles <- separate(data=prismFiles, col=product_name,into=c('month','blk', 'year', 'res', 'type'), sep=c(3,5,9,29))[,c(1,2,4,6)]
head(prismFiles)
unique(prismFiles$type)

setwd("~/DATA_prism/")
prismFilesWanted <- prismFiles %>% filter(year >= 1966 & year <= 2011 & 
                                            month %in% c('May', 'Jun', 'Jul') & 
                                            files %in%  Sys.glob("*.bil"))

prismStack <- raster::stack(prismFilesWanted$files)

# create BBS spatial points data frame
setwd("~/PIERRE/TRANSIENT DYNAMICS")
load('data/dataframe_comtrue_occurence.Rdata')
coordBBS <- read.table('data/coords_807_rt.txt', h=T)
colnames(coordBBS)[1] <- 'route'
bbsRouteCoord <- dfComOcc %>% distinct(route) %>% left_join(coordBBS)
coordinates(bbsRouteCoord) <- c('Longi', 'Lati')
proj4string(bbsRouteCoord) <- proj4string(prismStack)
plot(prismStack[[550]])
points(bbsRouteCoord)

####### extract on BBS site point
bbsRouteClim_temporary <-raster::extract(prismStack, bbsRouteCoord, df=T) #extract clim values
# bbsRouteClim_temporary <-as.data.frame(bbsRouteClim_temporary)
bbsRouteClim_temporary <- data.frame(as.data.frame(bbsRouteCoord),bbsRouteClim_temporary)
bbsRouteClim <- gather(bbsRouteClim_temporary, key='ref', value='value', 5:ncol(bbsRouteClim_temporary))
bbsRouteClim <- separate(data=bbsRouteClim, ref,into=c('dataset','climvar','stable','res', 'date', 'ext'), sep="_")[,c(1,2,3,6,9,11)] %>%
  separate(date,into=c('year', 'month'), sep=4)

bbsRouteClimWide <- bbsRouteClim %>% spread(climvar, value)

bbsRouteClimWide$MinMaxDiff <- bbsRouteClimWide$tmax - bbsRouteClimWide$tmin

library(viridis)
ggplot(bbsRouteClimWide, aes(tmax, tmin, col=MinMaxDiff))+
  geom_point()+
  # scale_color_viridis()+
  geom_abline(slope=1, intercept=0)+
  theme_classic()

cor.test(bbsRouteClimWide$tmax, bbsRouteClimWide$tmin)

ggplot(bbsRouteClimWide, aes(y=tmean, x=as.numeric(year)))+
  # geom_point(alpha=0.01)+
  theme_classic()+
  # geom_vline(xintercept = c(-8000))+
  stat_smooth(col='black')

ggplot(bbsRouteClimWide, aes(y=tmax, x=as.numeric(year)))+
  # geom_point(alpha=0.01)+
  theme_classic()+
  # geom_vline(xintercept = c(-8000))+
  stat_smooth(col='black')

ggplot(bbsRouteClimWide, aes(y=tmin, x=as.numeric(year)))+
  # geom_point(alpha=0.01)+
  theme_classic()+
  # geom_vline(xintercept = c(-8000))+
  stat_smooth(col='black')

library(mgcv)
# mod <- gamm(tmax ~ s(as.numeric(year)), random=list(route=~1), data=bbsRouteClimWide)  
# gam.check(mod$gam)
# summary(mod$gam)
# plot(mod$gam)
# 
# mod_tmean <- gamm(tmean ~ s(as.numeric(year)), random=list(route=~1), data=bbsRouteClimWide)  
# gam.check(mod_tmean$gam)
# summary(mod$gam)
# plot(mod_tmean$gam)
# 
# mod_tmean <- gamm(tmin ~ s(as.numeric(year)), random=list(route=~1), data=bbsRouteClimWide)  
# gam.check(mod_tmean$gam)
# summary(mod$gam)
# plot(mod_tmean$gam)


bbsRouteClimBreed_Tmin <- bbsRouteClim %>% group_by(route, Lati, Longi, year, climvar) %>% summarize( breed_mean = mean(value, na.rm=T)) 
# bbsRouteClimBreed <- ddply(bbsRouteClim, .(route, year, climvar), summarize,  val = mean(value, na.rm=T)) 

bbsRouteDiff <- bbsRouteClimBreed_Tmin %>% spread(key="climvar", value="breed_mean")  %>% mutate(MaxMinDiff = tmax-tmin) %>% 
  group_by(route, Lati, Longi) %>% summarize(meanRouteDiff = mean(MaxMinDiff), sdRouteDiff=sd(MaxMinDiff)) 

ggplot(bbsRouteDiff, aes(y=Lati, x=Longi))+
  geom_point(aes(col=meanRouteDiff), size=2)+
  scale_color_viridis()+
  theme_classic()

# save(bbsRouteDiff, file="results/bbsRouteMaxMinTempDiff.Rdata")

ggplot(bbsRouteClimWide, aes(y=tmin, x=tmax))+
  geom_hex(bins=100)+
  geom_abline(intercept=0, slope=1, col="black")+
  theme_classic()+
  stat_smooth(method="lm", col="red", linetype=2)+
  scale_fill_viridis()

test<- cor.test(bbsRouteClimWide$tmin, bbsRouteClimWide$tmax)

mod<- lm(bbsRouteClimWide$tmin ~ bbsRouteClimWide$tmax)
summary(mod)

####### Explore differences 
colnames(sitesClim)[1:3] <- c( "SiteName", 'site', "Altitude")
Average_sitesClim <- sitesClim %>% group_by(SiteName, site, Latitude, Longitude) %>%
  summarize(meanRouteDiff = mean(MinMaxDiff, na.rm=T))


load("results/statistics_sites_neotoma.Rdata")
stats_diff <- left_join(statistics_sites, Average_sitesClim)

###### models
library(mgcv)
## Linear model + smoothed spatial coordinates

mod.devabs <- gam(deviation.absolute.Tmax ~ s(Latitude, Longitude, bs='sos') + meanRouteDiff, 
                  data=stats_diff, na.action='na.fail')
plot(y=residuals(mod.devabs), fitted(mod.devabs))
summary(mod.devabs)

devabs<-
visreg(mod.devabs, "meanRouteDiff", type="conditional", gg=T)+
  labs(y= "deviation absolute", x="Tmax - Tmin (degree celsius)")


mod.devchange <- gam(deviation.change.Tmax ~ s(Latitude, Longitude, bs='sos') + meanRouteDiff, 
                    data=stats_diff, na.action='na.fail')
plot(y=residuals(mod.devchange), x=fitted(mod.devchange))
summary(mod.devchange)
devchange<-
visreg(mod.devchange, "meanRouteDiff", type="conditional", gg=T)+
  labs(y= "deviation change", x="Tmax - Tmin (degree celsius)")

mod.devcross <- lm(crossings.max.Tmax ~ meanRouteDiff, 
                    data=stats_diff, na.action='na.pass')
summary(mod.devcross)

mod.devcross <- gam(crossings.max.Tmax ~ s(Latitude, Longitude, bs='sos') +
                       meanRouteDiff, 
                     data=stats_diff, family="poisson", na.action='na.fail')
summary(mod.devcross)

crossings<-
visreg(mod.devcross, "meanRouteDiff", type="conditional", gg=T)+
  labs(y= "maximum state number", x="Tmax - Tmin (degree celsius)")


png("results/figures/review_Tmax_Tmin_scatterplots.png",width =12, height=4, units="in", res=200)
plot_grid(devabs, devchange, crossings, ncol=3,  axis='bottom')
dev.off()


#### correlation between Tmax and Tmin for Paleo ####
setwd("~/PIERRE/TRANSIENT DYNAMICS")
library(tidyverse)
library(raster)
library(broom)


#### get climate paleo time-series data ####
##define directory
dir.paleoclim <- 'data/blois data/CCSM3_500_With_PaleoShorelines'
climvars <- c('tmax_year_ave','tmin_year_ave','prcp_year_ave')

##get rasters
getrasters <- function(climvar,times=seq(0,22000,by=500))
{
  r <- stack(lapply(times,function(time) {raster(sprintf("%s/%dBP/%s.tif",dir.paleoclim,time,climvar))}))
  names(r) <- paste("X",times,sep="")
  return(r)
}
climdata <- lapply(climvars, getrasters)

##separate rasters
TmaxRaster <- climdata[[1]]
proj4string(TmaxRaster) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

TminRaster <- climdata[[2]]
proj4string(TminRaster) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

PrecRaster <- climdata[[3]]
proj4string(PrecRaster) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

## get spatial info on neotoma sites
latlon <- read.csv('data/blois data/pollen old/All.site.data-withagemodel-finalv2.csv')[,c("SiteName","Handle", "Longitude", "Latitude", "Altitude")]

#Convert to spatial points data frame
coordinates(latlon) <- c('Longitude', 'Latitude')
#put on the same geo-referencial
proj4string(latlon) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

#Plot to check spatial match
plot(TmaxRaster[[1]], main="Average Tmax")
plot(latlon, add=T)

plot(TminRaster[[1]], main="Average Tmin")
plot(latlon, add=T)

plot(PrecRaster[[1]], main="Pecipitation")
plot(latlon, add=T)

## extract mean tmax and prec value for each site and time
sitesClim<- data.frame(latlon$SiteName, 
                       latlon$Handle, 
                       latlon$Altitude, 
                       coordinates(latlon),
                       raster::extract(TmaxRaster, latlon,na.rm=T),
                       raster::extract(TminRaster, latlon,na.rm=T),
                       raster::extract(PrecRaster, latlon,na.rm=T))   #extract clim values


sitestmax <- sitesClim %>% gather(key="time", "ClimTmax.ave", 6:50) %>% 
  dplyr::select(1,2,3,4,5,96,97) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))

sitestmin <- sitesClim %>% gather(key="time", "ClimTmin.ave", 51:95) %>% 
  dplyr::select(1,2,3,4,5,96,97) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))

sitesprec <- sitesClim %>% gather(key="time", "ClimPrec.ave", 96:140) %>% 
  dplyr::select(1,2,3,4,5,96,97) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))

sitesClim <- sitestmax

sitesClim$ClimTmin.ave <- sitestmin$ClimTmin.ave
sitesClim$ClimPrec.ave <- sitesprec$ClimPrec.ave

sitesClim$MinMaxDiff <- sitesClim$ClimTmax.ave - sitesClim$ClimTmin.ave

library(viridis)
ggplot(sitesClim, aes(ClimTmax.ave, ClimTmin.ave, col=time))+
  # geom_point()+
  geom_line(aes(group=latlon.SiteName))+
  scale_color_viridis()+
  geom_abline(slope=1, intercept=0)+
  theme_classic()

cor.test(sitesClim$ClimTmax.ave, sitesClim$ClimTmin.ave)


mod<-lm(sitesClim$ClimTmax.ave ~ sitesClim$ClimTmin.ave)
summary(mod)

ggplot(sitesClim, aes(y=ClimTmin.ave, x=-time))+
  geom_point(alpha=0.01)+
  geom_point(aes(y=ClimTmax.ave), col="red", alpha=0.01)+
  theme_classic()+
  # geom_vline(xintercept = c(-8000))+
  stat_smooth(col='black')+
  stat_smooth(col='red', aes(y=ClimTmax.ave))

ggplot(sitesClim, aes(y=ClimTmax.ave, x=-time))+
  geom_point(alpha=0.01)+
  theme_classic()+
  # geom_vline(xintercept = c(-8000))+
  stat_smooth(col='black')

ggplot(bbsRouteClimWide, aes(y=tmin, x=as.numeric(year)))+
  # geom_point(alpha=0.01)+
  theme_classic()+
  # geom_vline(xintercept = c(-8000))+
  stat_smooth(col='black')



####### Explore differences 
load("results/statistics_sites_neotoma.Rdata")
colnames(bbsRouteDiff)[1]<- "site"

stats_diff <- left_join(statistics_sites, sitesClim[, c(1,2)])

density <- 
  ggplot(stats_diff[stats_diff$deviation.absolute < 8,], aes(meanRouteDiff, fill=..x..))+
  geom_histogram(col="black")+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlGnBu")+
  # geom_vline(xintercept = mean(stats_diff$meanRouteDiff))+
  # geom_vline(xintercept = c(mean(stats_diff$meanRouteDiff) - sd(stats_diff$meanRouteDiff),
  # mean(stats_diff$meanRouteDiff) + sd(stats_diff$meanRouteDiff)), linetype="dashed")+
  xlab("Tmax - Tmin (degree celsius)")
library(maps)
library(ggplot2)
library(rje)

world <-map_data("world") 

worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,52), xlim=c(-125, -65))+
  theme(legend.position="bottom")

map <-
  worldmap + 
  geom_point(data=stats_diff, aes(Longi, Lati , fill=meanRouteDiff),size=3 , shape=21)+
  scale_fill_distiller(palette = "YlGnBu")+
  theme_classic()+
  theme(legend.position="none")
# annotate("text", label="Maximal temperature", x=-70,y=30, size=5, col="red")

devabs<-
  ggplot(stats_diff[stats_diff$deviation.absolute < 8,], aes(x=deviation.absolute, y=meanRouteDiff))+
  geom_pointrange(aes(ymin=meanRouteDiff - sdRouteDiff, ymax=meanRouteDiff + sdRouteDiff))+
  theme_classic()
devchange<-
  ggplot(stats_diff, aes(x=deviation.change.1Ky, y=meanRouteDiff))+
  geom_pointrange(aes(ymin=meanRouteDiff - sdRouteDiff, ymax=meanRouteDiff + sdRouteDiff))+
  theme_classic()
crossings<-
  ggplot(stats_diff, aes(x=crossings.max, y=meanRouteDiff))+
  geom_pointrange(aes(ymin=meanRouteDiff - sdRouteDiff, ymax=meanRouteDiff + sdRouteDiff))+
  theme_classic()

library(cowplot)
png("results/figures/review_Tmax_Tmin_map.png", width =12, height=4, units="in", res=200)
plot_grid(density, map, ncol=2,  axis='bottom', rel_widths =  c(1/3, 2/3))
dev.off()



###### models
library(mgcv)
## Linear model + smoothed spatial coordinates

mod.devabs <- gam(deviation.absolute ~ s(Lati, Longi, bs='sos') + meanRouteDiff, 
                  data=stats_diff, na.action='na.fail')
plot(y=residuals(mod.devabs), fitted(mod.devabs))
summary(mod.devabs)

devabs<-
  visreg(mod.devabs, "meanRouteDiff", type="conditional", gg=T)+
  labs(y= "deviation absolute", x="Tmax - Tmin (degree celsius)")


mod.devchange <- lm(deviation.change ~ meanRouteDiff, 
                    data=stats_diff, na.action='na.fail')
plot(y=residuals(mod.devchange), fitted(mod.devchange))
summary(mod.devchange)

mod.devchange <- gam(deviation.change ~ s(Lati, Longi, bs='sos') +
                       meanRouteDiff, 
                     data=stats_diff,
                     correlation= corSpher(form = ~ Lati + Longi), na.action='na.fail')
plot(y=residuals(mod.devchange), fitted(mod.devchange))
summary(mod.devchange)
devchange<-
  visreg(mod.devchange, "meanRouteDiff", type="conditional", gg=T)+
  labs(y= "deviation change", x="Tmax - Tmin (degree celsius)")



mod.devcross <- lm(crossings.max ~ meanRouteDiff, 
                   data=stats_diff, na.action='na.fail')
summary(mod.devcross)

mod.devcross <- gam(crossings.max ~ s(Lati, Longi, bs='sos') +
                      meanRouteDiff, 
                    data=stats_diff,
                    correlation= corSpher(form = ~ Lati + Longi), na.action='na.fail')
summary(mod.devcross)

crossings<-
  visreg(mod.devcross, "meanRouteDiff", type="conditional", gg=T)+
  labs(y= "maximum state number", x="Tmax - Tmin (degree celsius)")


png("results/figures/review_Tmax_Tmin_scatterplots.png",width =12, height=4, units="in", res=200)
plot_grid(devabs, devchange, crossings, ncol=3,  axis='bottom')
dev.off()


###### Birds traits not considered in the study #######

setwd("~/PIERRE/TRANSIENT DYNAMICS/")
library(tidyverse)

#### 1. read BBS species traits
traits <- read.table('data/BBS_species_traits.csv', h=T, sep=';')
#choose traits of interest
traitsOfInterest <- traits[,c("species","diet")]

ggplot(traitsOfInterest, aes(diet))+
  geom_histogram(stat="count")

load('data/dataframe_comtrue_occurence.Rdata')

#### 2. Compute Community mean Values 

dfComOccTraits <- left_join(dfComOcc, traitsOfInterest)
ggplot(dfComOccTraits, aes(diet))+
  geom_histogram(stat="count")

library(plyr)
BirdsTraits_means_weighted   <- plyr::ddply(na.omit(dfComOccTraits), .(route), function(x){
  
  seed = (sum(x$N[x$diet=='seed'])/sum(x$N))*100
  omnivore = (sum(x$N[x$diet=='omnivore'])/sum(x$N))*100
  invertebrate = (sum(x$N[x$diet=='invertebrate'])/sum(x$N))*100
  vertfishscav = (sum(x$N[x$diet=='vertfishscav'])/sum(x$N))*100
  
return(cbind(seed, omnivore, invertebrate, vertfishscav))}, .progress='text')

colnames(BirdsTraits_means_weighted)[1]<- 'site'

#### 4. test effect of different diet on statistics
load("results/statistics_sites_BBS.Rdata")

#### merge datas
modData <- na.omit(left_join(statistics_sites, BirdsTraits_means_weighted))
colnames(modData)[2:3]<- c("Longitude", "Latitude")

##### map diets
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

seed <-
  worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=seed),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")

omnivore <-
  worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=omnivore),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")

invertebrate <-
  worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=invertebrate),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")

vertfishscav <-
  worldmap + 
  geom_point(data=modData, aes(y=Latitude, x=Longitude, fill=vertfishscav),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="bottom")
library(cowplot)
png("results/figures/review_diet_maps.png",width =9, height=20, units="in", res=200)
plot_grid(seed, omnivore, invertebrate, vertfishscav, ncol=1,  axis='left')
dev.off()

#### test effect on statistics

#DEVABS
library(mgcv)
modSpa.devAbs.Tmax <- gam(deviation.absolute  ~ s(Latitude, Longitude, bs='sos') +
                            scale(seed) + 
                            scale(omnivore) +
                            scale(invertebrate) +
                            scale(vertfishscav),  data=modData, na.action='na.fail')
summary(modSpa.devAbs.Tmax)
library(visreg)
visreg(modSpa.devAbs.Tmax,"invertebrate", gg=T)+
  ylab("deviation absolute seed partial residuals")+
  theme_classic()

#DEVCHANGE
modSpa.devAbs.Tmax <- gam(deviation.absolute  ~ s(Latitude, Longitude, bs='sos') +
                            scale(seed) + 
                            scale(omnivore) +
                            scale(invertebrate) +
                            scale(vertfishscav),  data=modData, na.action='na.fail')
summary(modSpa.devAbs.Tmax)
visreg(modSpa.devAbs.Tmax,'seed', gg=T)+
  ylab("deviation absolute seed partial residuals")+
  theme_classic()

modSpa.devChange.Tmax <- gam(deviation.change.1Ky  ~ s(Latitude, Longitude, bs='sos') + 
                               scale(seed) + 
                               scale(omnivore) +
                               scale(invertebrate) +
                               scale(vertfishscav),  data=modData, na.action='na.fail')
hist(modSpa.devChange.Tmax$residuals, breaks=100)
summary(modSpa.devChange.Tmax)
visreg(modSpa.devChange.Tmax,'invertebrate', gg=T)+
  ylab("deviation change Tmax partial residuals")+
  theme_classic()+
  geom_hline(yintercept=0)

# N
modSpa.N.Tmax <- gam(crossings.max ~ s(Latitude, Longitude, bs='sos') + 
                       scale(seed) + 
                       scale(omnivore) +
                       scale(invertebrate) +
                       scale(vertfishscav),  data=modData, na.action='na.fail')
hist(modSpa.N.Tmax$residuals, breaks=100)
summary(modSpa.N.Tmax)


#### merge datas
load("results/predictors_BBS.Rdata")
modData2 <- na.omit(left_join(modData, predictors))
modData2$sqrt_topo25 <- sqrt(modData2$topo25)

##### correlation between diet at community scale 
library("PerformanceAnalytics")
chart.Correlation(modData2[, c("Latitude", 
                               "Longitude",
                               "base.Tmax",
                               "sqrt_topo25",
                               "HII",
                               "BodyMass", 
                               'Migration', "seed", 
                               "invertebrate",
                               "omnivore", 
                               "vertfishscav")], histogram=TRUE, pch=19)

chart.Correlation(modData2[, c("BodyMass", 
                               'Migration', "seed", 
                               "invertebrate",
                               "omnivore", 
                               "vertfishscav")], histogram=TRUE, pch=19)

#DEVABS
modSpa.devAbs.Tmax<- gam(deviation.absolute  ~ s(Latitude, Longitude, bs='sos') +
                            scale(poly(base.Tmax,2)) + 
                            scale(sqrt_topo25) +
                            scale(HII) +
                            scale(BodyMass) + 
                            scale(Migration)+
                            scale(invertebrate),  data=modData2, na.action='na.fail', REML=F)
hist(modSpa.devAbs.Tmax$residuals, breaks=100)
summary(modSpa.devAbs.Tmax_invertebrate)
BIC(modSpa.devAbs.Tmax_invertebrate)
AIC(modSpa.devAbs.Tmax_invertebrate)


#DEVABS
modSpa.devAbs.Tmax <- gam(deviation.absolute  ~ s(Latitude, Longitude, bs='sos') +
                            scale(poly(base.Tmax,2)) + 
                            scale(sqrt_topo25) +
                            scale(HII) +
                            scale(BodyMass) + 
                            scale(Migration),  data=modData2, na.action='na.fail', REML=F)
hist(modSpa.devAbs.Tmax$residuals, breaks=100)
summary(modSpa.devAbs.Tmax)
AIC(modSpa.devAbs.Tmax) - AIC(modSpa.devAbs.Tmax_invertebrate)
BIC(modSpa.devAbs.Tmax) - BIC(modSpa.devAbs.Tmax_invertebrate)

modSpa.devAbs.Tmax <- gam(deviation.absolute  ~ s(Latitude, Longitude, bs='sos') +
                            poly(base.Tmax,2) + 
                            scale(sqrt_topo25) +
                            scale(HII) +
                            scale(BodyMass) + 
                            scale(Migration),  data=modData2, na.action='na.fail')
hist(modSpa.devAbs.Tmax$residuals, breaks=100)
summary(modSpa.devAbs.Tmax)

visreg(modSpa.devAbs.Tmax, "BodyMass", type="conditional", gg=T)+
  labs(y= "deviation absolute", x="mean bodymass")


modSpa.devAbs.Tmax <- gam(deviation.absolute  ~ s(Latitude, Longitude, bs='sos') +
                            scale(poly(base.Tmax,2)) + 
                            scale(sqrt_topo25) +
                            scale(HII) +
                            scale(BodyMass) + 
                            scale(Migration)+
                            scale(invertebrate),  data=modData2, na.action='na.fail')
hist(modSpa.devAbs.Tmax$residuals, breaks=100)
summary(modSpa.devAbs.Tmax)


modSpa.devAbs.Tmax <- gam(deviation.absolute  ~ s(Latitude, Longitude, bs='sos') +
                            scale(base.Tmax) + 
                            scale(sqrt_topo25) +
                            scale(HII) +
                            scale(BodyMass) + 
                            scale(Migration)+
                            scale(invertebrate),  data=modData2, na.action='na.fail')
visreg(modSpa.devAbs.Tmax)

visreg(modSpa.devAbs.Tmax, "invertebrate", type="conditional", gg=T)+
  labs(y= "deviation absolute", x="% insectivorous species")

#DEVCHANGE
modSpa.devChange.Tmax <- gam(deviation.change.1Ky  ~ s(Latitude, Longitude, bs='sos') + 
                               scale(poly(base.Tmax,2)) + 
                               scale(sqrt_topo25) +
                               scale(HII) +
                               scale(BodyMass) + 
                               scale(Migration)+
                               scale(invertebrate),  data=modData2, na.action='na.fail')

hist(modSpa.devChange.Tmax$residuals, breaks=100)
summary(modSpa.devChange.Tmax)
visreg(modSpa.devChange.Tmax, "invertebrate", type="conditional", gg=T)+
  labs(y= "deviation change", x="% insectivorous species")

# N
modSpa.N.Tmax <- gam(crossings.max  ~ s(Latitude, Longitude, bs='sos') + 
                               scale(poly(base.Tmax,2)) + 
                               scale(sqrt_topo25) +
                               scale(HII) +
                               scale(BodyMass) + 
                               scale(Migration)+
                               scale(invertebrate),  data=modData2,  family='poisson', na.action='na.fail')

hist(modSpa.N.Tmax$residuals, breaks=100)
summary(modSpa.N.Tmax)
visreg(modSpa.N.Tmax, "invertebrate", type="conditional", gg=T, scale="linear")+
  labs(y= "state number", x="% insectivorous species")

#### Summary of models

summary(modSpa.devAbs.Tmax)
summary(modSpa.devChange.Tmax)
summary(modSpa.N.Tmax)

visreg(modSpa.devAbs.Tmax)
visreg(modSpa.devChange.Tmax)
visreg(modSpa.N.Tmax)


modPred_res <-as.data.frame(rbind(
  round(summary(modSpa.N.Tmax)$p.table, digits=4)[-1,],
  round(summary(modSpa.devAbs.Tmax)$p.table, digits=4)[-1,], 
  round(summary(modSpa.devChange.Tmax)$p.table, digits=4)[-1,]), row.names = 1:18)

modPred_res$predictor<- factor(rep(c("Baseline.Climate","Baseline.Climate2", "Topography", "Human.Influence", "Body.mass", "Migration", "Insectivorous"),3), 
                               levels=c("Baseline.Climate","Baseline.Climate2", "Topography", "Human.Influence", "Body.mass", "Migration", "Insectivorous"))
modPred_res$expl.var <- rep(c("max.state.number", "absolute.deviation","deviation.change"), each=7)
colnames(modPred_res)<- c("coefficient", "SE", "z_value", "P_value","predictor","expl.var")


# modPred_res$coefficient<- as.numeric(as.character(modPred_res$coefficient))
# modPred_res$SE<- as.numeric(as.character(modPred_res$SE))
modPred_res$sig <- " "
modPred_res$sig[modPred_res$P_value < 0.1] <- "."
modPred_res$sig[modPred_res$P_value < 0.05] <- "*"
modPred_res$sig[modPred_res$P_value < 0.01] <- "**"

modPred_res_BBS <- modPred_res
colnames(modPred_res_BBS)<- c("coefficient", "SE", "z_value", "P_value","predictor","expl.var")
modPred_res_BBS$predictor <- factor(modPred_res_BBS$predictor, levels=c("Baseline.Climate","Baseline.Climate2", "Topography", "Human.Influence", "Body.mass", "Migration", "Invertebrate"))
modPred_res_BBS$expl.var<- rep(c("max.state.number", "absolute.deviation","deviation.change"), each=7)

modPred_res_BBS$sig <- "ns"
modPred_res_BBS$sig[modPred_res_BBS$P_value < 0.1] <- "."
modPred_res_BBS$sig[modPred_res_BBS$P_value < 0.05] <- "* P<0.05"
modPred_res_BBS$sig[modPred_res_BBS$P_value < 0.01] <- "** P<0.01"


stat_names<- list('absolute.deviation'= expression(paste("Absolute deviation ", bar(Lambda))), 
                  'deviation.change'=expression(paste("Deviation change ", italic(d), Lambda)), 
                  'max.state.number'= expression(paste("Max. state number ", italic(n))))

stat_labeller <- function(variable,value){return(stat_names[value])}

plot_modPred_res_BBS<- 
  ggplot(modPred_res_BBS, aes(x=predictor, y=coefficient, col=sig))+
  geom_pointrange(aes(ymin=coefficient-(1.96*SE), ymax=coefficient+(1.96*SE)))+
  geom_hline(yintercept = 0)+
  facet_wrap(~expl.var,nrow=1, scales='free_y', labeller=stat_labeller)+
  # scale_color_brewer(palette = "YlGnBu")+
  # scale_color_manual(values=c(rgb(161,218,180,maxColorValue = 255), rgb(37,52,148,maxColorValue = 255),rgb(65,182,196,maxColorValue = 255)))+
  # scale_color_manual(values=c('grey80','grey60', 'grey30', 'black'))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())+
  scale_x_discrete(labels=c("Baseline Climate",expression(paste("Baseline Climate"^"2")), "Topography", "Human Influence", "Body mass", "Migration"))+
  ylab("Birds - Model coefficients")




ggplot(modPred_res_BBS, aes(x=predictor, y=coefficient, col=sig))+
  geom_pointrange(aes(ymin=coefficient-(1.96*SE), ymax=coefficient+(1.96*SE)))+
  geom_hline(yintercept = 0)+
  facet_wrap(~expl.var,nrow=1, scales='free_y', labeller=stat_labeller)+
  scale_color_manual(values=c('grey80','grey60', 'grey30', 'black'))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=60, hjust=1), axis.title.x=element_blank())

levels(modPred_res_BBS$predictor)

write.csv(modPred_res, file="results/modPred_results_BBS.csv", row.names = F)




##### sampling bias ##########

##BEN IDEA !!
occurrences_bienWithWorldclim <- readRDS("results/occurrences_bienWithWorldclim.Rdata")

library(tidyverse)
library(colorRamps)
library(rje)
library(ggplot2)
library(ggExtra)

)
species_list <- unique(occurrences_bienWithWorldclim$alltaxa)

species_convex_niches <- data.frame("species"=species_list, obs_mean=999, obs_sd=999, convex_mean=999, convex_sd=999, n_occurences=999)

for(i in 1:length(species_list)){
  sp<-species_list[i]
  print(sp)
  occurences_sp <- occurrences_bienWithWorldclim %>% filter(alltaxa==sp)
  sp_Temp <- na.omit(occurences_sp$Tmax)
  CIinfsp <- mean(sp_Temp)-1.96*sd(sp_Temp)
  CIsupsp <- mean(sp_Temp) + 1.96*sd(sp_Temp)
  sp_TempBorne <- sp_Temp[sp_Temp > CIinfsp & sp_Temp < CIsupsp]
  sp_distrib <- data.frame("realDistribution"=sp_TempBorne, "convexDistribution"=runif(length(sp_TempBorne), min=CIinfAbies, max=CIsupAbies))
  
  species_convex_niches[i, 2:6] <- c(mean(sp_TempBorne), 
                                     sd(sp_TempBorne), 
                                     mean(runif(length(sp_TempBorne), min=CIinfsp, max=CIsupsp)),
                                     sd(runif(length(sp_TempBorne), min=CIinfsp, max=CIsupsp)), 
                                     length(sp_Temp))
}


ggplot(species_convex_niches, aes(y=obs_mean, x=convex_mean))+
  geom_errorbar(aes(ymin=obs_mean-obs_sd, ymax=obs_mean+obs_sd))+
  geom_errorbarh(aes(xmin=convex_mean-convex_sd, xmax=convex_mean+convex_sd))+
  geom_point(shape=21, fill="white", aes(size=n_occurences))+
  geom_abline(interecept0, slope=1, col="red")+
  theme_classic()


cor.test(species_convex_niches$obs_mean, species_convex_niches$convex_mean)







world <-map_data("world") 

worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  theme(legend.position="bottom")+
  coord_map("albers", parameters = c(-100, -100))


occurences <- data.frame(longi=round(occurencesdf$longitude,1), lati=round(occurencesdf$latitude,1))
grid<- distinct(occurences)

Abies <- sample_n(na.omit(occurencesdf[occurencesdf$genus == "Abies",c("latitude", "longitude")]), 500)
Abies$z <- 0
contour(as.matrix(Abies))

#create convex hull
ch <- chull(Abies)

worldmap+
geom_point(data=Abies, aes(longitude, latitude))+
  theme_classic()+
  geom_path(data=Abies[c(ch, ch[1]),], aes(longitude, latitude), col="red")



inter1 <- Abies

inter1 <-read.table("~/inter1.txt")
plot(inter1[,c("lat", "long")])

#add a category (required for later rasterizing/polygonizing)
inter <- cbind(inter1, cat = rep(1L, nrow(inter1)),stringsAsFactors = FALSE)


#convert to spatial points
coordinates(inter1) = ~ longitude + latitude
plot(inter1)




#gridify your set of points
gridded(inter1) <- TRUE
#convert to raster
library(raster)
r <- raster(inter1)

library(raster)
# set up an 'empty' raster, here via an extent object derived from your data
e <- extent(inter1)
plot(e)
plot(inter1, add=T)

gridded(e) <- T

r<- rasterize(inter1)

r <- raster(e, ncol=10, nrow=2)
# or r <- raster(xmn=, xmx=,  ...

# you need to provide a function 'fun' for when there are multiple points per cell
r <- rasterize(inter1, r, Abies, fun=mean)
plot(x)

#convert raster to polygons
sp = rasterToPolygons(r, dissolve = T)



#addition transformation to distinguish well the set of polygons
polys <- slot(sp@polygons[[1]], "Polygons")

spline.poly <- function(xy, vertices, k=3, ...) {
  # Assert: xy is an n by 2 matrix with n >= k.
  
  # Wrap k vertices around each end.
  n <- dim(xy)[1]
  if (k >= 1) {
    data <- rbind(xy[(n-k+1):n,], xy, xy[1:k, ])
  } else {
    data <- xy
  }
  
  # Spline the x and y coordinates.
  data.spline <- spline(1:(n+2*k), data[,1], n=vertices, ...)
  x <- data.spline$x
  x1 <- data.spline$y
  x2 <- spline(1:(n+2*k), data[,2], n=vertices, ...)$y
  
  # Retain only the middle part.
  cbind(x1, x2)[k < x & x <= n+k, ]
}

output <- SpatialPolygons(
  Srl = lapply(1:length(polys),
               function(x){
                 p <- polys[[x]]
                 
                 #applying spline.poly function for smoothing polygon edges
                 px <- slot(polys[[x]], "coords")[,1]
                 py <- slot(polys[[x]], "coords")[,2]
                 bz <- spline.poly(slot(polys[[x]], "coords"),100, k=3)
                 bz <- rbind(bz, bz[1,])
                 slot(p, "coords") <- bz               
                 
                 # create Polygons object
                 poly <- Polygons(list(p), ID = x)
                 return(poly)
               }),
  proj4string = CRS("+init=epsg:4326")
)

#plot
plot(inter1)
plot(sp, add=T,border = "gray", lwd = 2) #polygonize result
plot(output, add=T,border = "red",  add = TRUE) #smoothed polygons




ggbase <- ggplot(occurencesdf[occurencesdf$longitude < 0 & occurencesdf$latitude < 90,], aes(y=latitude, x=longitude))+
  geom_bin2d(bins=200)+
  theme_classic()

ggbase +
  scale_fill_viridis_c(trans="log10")


grid$lati <- as.factor(grid$lati)
grid$longi <- as.factor(grid$longi)

ggplot(occure)

test <- occurences %>% group_by(lati, longi) %>% summarize(n.obs = nrow())

scale_fill_continuous()
