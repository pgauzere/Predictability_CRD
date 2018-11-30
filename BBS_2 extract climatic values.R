setwd("~/PIERRE/TRANSIENT DYNAMICS")

library(prism)
library(tidyverse)
library(sp)

options(prism.path = "data/prism")
prismFiles <- ls_prism_data(name = TRUE)
prismFiles <- separate(data=prismFiles, col=product_name,into=c('month','blk', 'year', 'res', 'type'), sep=c(3,5,9,29))[,c(1,2,4,6)]
head(prismFiles)
unique(prismFiles$type)

setwd("data/prism/")
prismFilesWanted <- prismFiles %>% filter(year >= 1966 & year <= 2011 & 
                                            month %in% c('May', 'Jun', 'Jul') & 
                                            files %in%  Sys.glob("*.bil"))

prismStack <- raster::stack(prismFilesWanted$files)

setwd("~/PIERRE/TRANSIENT DYNAMICS")

# create BBS spatial points data frame
load('data/dataframe_comtrue_occurence.Rdata')
coordBBS <- read.table('data/coords_807_rt.txt', h=T)
colnames(coordBBS)[1] <- 'route'
bbsRouteCoord <- dfComOcc %>% distinct(route) %>% left_join(coordBBS)
coordinates(bbsRouteCoord) <- c('Longi', 'Lati')
proj4string(bbsRouteCoord) <- proj4string(prismStack)
plot(prismStack[[1]])
points(bbsRouteCoord)

#### extract on BBS site point ####
bbsRouteClim_temporary <-raster::extract(prismStack, bbsRouteCoord, df=T) #extract clim values
# bbsRouteClim_temporary <-as.data.frame(bbsRouteClim_temporary)
bbsRouteClim_temporary <- data.frame(as.data.frame(bbsRouteCoord),bbsRouteClim_temporary)
bbsRouteClim <- gather(bbsRouteClim_temporary, key='ref', value='value', 5:ncol(bbsRouteClim_temporary))
bbsRouteClim <- separate(data=bbsRouteClim, ref,into=c('dataset','climvar','stable','res', 'date', 'ext'), sep="_")[,c(1,2,3,6,9,11)] %>%
  separate(date,into=c('year', 'month'), sep=4)

bbsRouteClimWide <- bbsRouteClim %>% spread(climvar, value)

ggplot(bbsRouteClimWide, aes(tmax, tmean, col=year))+
  geom_point()

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

library(mgcv)
mod <- gamm(tmax ~ s(as.numeric(year)), random=list(route=~1), data=bbsRouteClimWide)  
gam.check(mod$gam)
summary(mod$gam)
plot(mod$gam)

mod_tmean <- gamm(tmean ~ s(as.numeric(year)), random=list(route=~1), data=bbsRouteClimWide)  
gam.check(mod_tmean$gam)
summary(mod$gam)
plot(mod_tmean$gam)


bbsRouteClimBreed_Tmean <- bbsRouteClim %>% group_by(route, Lati, Longi, year, climvar) %>% summarize( breed_mean = mean(value, na.rm=T)) 
# bbsRouteClimBreed <- ddply(bbsRouteClim, .(route, year, climvar), summarize,  val = mean(value, na.rm=T)) 
bbsRouteClimBreed_Tmean <- bbsRouteClimBreed_Tmean %>% spread(key=climvar, value=breed_mean)
save(bbsRouteClimBreed_Tmean, file="results/bbsRouteClimBreed_Tmean.Rdata")

#### extract on 50km buffers around each BBS site point ####

bbsRouteClim_temporary_mean <-raster::extract(prismStack, bbsRouteCoord, buffer=50000, small=T, fun = mean, na.rm=T) #extract clim values
save(bbsRouteClim_temporary_mean, file='results/bbsRouteClim_temporary_mean.Rdata')
bbsRouteClim_temporary_sd <-raster::extract(prismStack, bbsRouteCoord, buffer=50000, small=T, fun = sd, na.rm=T) #extract clim values
save(bbsRouteClim_temporary_sd, file='results/bbsRouteClim_temporary_sd.Rdata')

bbsRouteClim_temporary_mean <-as.data.frame(bbsRouteClim_temporary_mean)
bbsRouteClim_temporary_mean <- cbind(bbsRouteCoord,bbsRouteClim_temporary_mean)


bbsRouteClim_temporary_sd <-as.data.frame(bbsRouteClim_temporary_sd)
bbsRouteClim_temporary_sd <- cbind(bbsRouteCoord,bbsRouteClim_temporary_sd)

# bbsRouteClim_temporary$route <- bbsRouteCoord$route
# bbsRouteClim_temporary_sd <- data.frame(bbsRouteCoord$route,
#                                         extract(prismStack, bbsRouteCoord, buffer=50000, fun=sd, na.rm=T, small=T)) #extract clim values

bbsRouteClim_mean <- gather(bbsRouteClim_temporary_mean, key='ref', value='mean', 4:ncol(bbsRouteClim_temporary_mean))
bbsRouteClim_mean <- separate(data=bbsRouteClim_mean, ref,into=c('dataset','climvar','stable','res', 'date', 'ext'), sep="_")[,c(1,2,3,5,8,10)] %>%
  separate(date,into=c('year', 'month'), sep=4)

bbsRouteClim_sd <- gather(bbsRouteClim_temporary_sd, key='ref', value='sd', 4:ncol(bbsRouteClim_temporary_sd))
bbsRouteClim_sd <- separate(data=bbsRouteClim_sd, ref,into=c('dataset','climvar','stable','res', 'date', 'ext'), sep="_")[,c(1,2,3,5,8,10)] %>%
  separate(date,into=c('year', 'month'), sep=4)

bbsRouteClim<- left_join(bbsRouteClim_mean, bbsRouteClim_sd)
bbsRouteClim$route   <- as.factor(bbsRouteClim$route)
bbsRouteClim$climvar <- as.factor(bbsRouteClim$climvar)
bbsRouteClim$year    <- as.numeric(bbsRouteClim$year)
bbsRouteClim$month   <- as.factor(bbsRouteClim$month)
bbsRouteClim$mean   <- as.numeric(bbsRouteClim$mean)
bbsRouteClim$sd   <- as.numeric(bbsRouteClim$sd)

bbsRouteClimBreed <- bbsRouteClim %>% group_by(route, Lati, Longi, year, climvar) %>% summarize( breed_mean = mean(mean, na.rm=T), breed_sd = sd(sd, na.rm=T)) 
# bbsRouteClimBreed <- ddply(bbsRouteClim, .(route, year, climvar), summarize,  val = mean(value, na.rm=T)) 

head(bbsRouteClimBreed)

bbsRouteClimBreed <- bbsRouteClimBreed %>% spread(key=climvar, value=val)
save(bbsRouteClimBreed, file="results/bbsRouteClimBreed.Rdata")




