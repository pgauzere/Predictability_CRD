library(raster)
library(maps)
library(plyr)
library(tidyverse)

setwd("~/PIERRE/TRANSIENT DYNAMICS /")

#### 1.Create community DF from pollen Neotoma datas ####
#### get community data ####
#neotoma data presence absence directory
dir.presabs <- 'data/blois data/output_for_BenBlonder/all data by site'
filenames.presabs <-dir(dir.presabs,pattern='*PA\\.csv',full.names=T)


filenames.presabs <-dir(dir.presabs,pattern='*PA\\.csv',full.names=T)
files.presabs <- lapply(filenames.presabs, read.csv)
names(files.presabs) <- gsub(dir.presabs,"",filenames.presabs,fixed=T)
names(files.presabs) <- gsub("/","",names(files.presabs),fixed=T)
names(files.presabs) <- gsub(".interp.data.PA.csv","",names(files.presabs),fixed=T)

files.presabs <- lapply(files.presabs, function(x) { 
  row.names(x) <- x[,"sample.ages"]
  x <- x[,2:ncol(x)]
  names(x) <- gsub(".type","",names(x),fixed=T)
  names(x) <- gsub("."," ",names(x),fixed=T)
  return(x)
} )

presAbs <- files.presabs %>% ldply(add_rownames)
colnames(presAbs)[1:2] <- c("site", "time")
presAbs$site <- as.factor(presAbs$site)
presAbs$time <- as.numeric(presAbs$time)

paleoOcc <- presAbs %>% gather("taxon", "occurence", 3:105) %>% filter(occurence == 1)


#### 2.Merge Climatic data from Paleo Clim CCSM3 to presence absence data from Neotoma ####

## 3.1 get clim data for each site each time
##define directory
dir.paleoclim <- 'data/blois data/CCSM3_500_With_PaleoShorelines'
climvars <- c('tmax_year_ave','prcp_year_ave')

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
PrecRaster <- climdata[[2]]
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
plot(PrecRaster[[1]], main="Pecipitation")
plot(latlon, add=T)

## extract mean tmax and prec value for each site and time
sitesClim<- data.frame(latlon$SiteName, 
                       latlon$Handle, 
                       latlon$Altitude, 
                       coordinates(latlon),
                       raster::extract(TmaxRaster, latlon,na.rm=T),
                       raster::extract(PrecRaster, latlon,na.rm=T))   #extract clim values


sitestmax <- sitesClim %>% gather(key="time", "ClimTmax.ave", 6:50) %>% 
  dplyr::select(1,2,3,4,5,51,52) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))

sitesprec <- sitesClim %>% gather(key="time", "ClimPrec.ave", 51:95) %>% 
  dplyr::select(1,2,3,4,5,51,52) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))

sitesClim <- sitestmax
sitesClim$ClimPrec.ave <- sitesprec$ClimPrec.ave

colnames(sitesClim)[1:3] <-  c("SiteName","site", "altitude")
save(sitesClim, file="results/sitesClim_neotoma.Rdata")

load("results/sitesClim_neotoma.Rdata")

nrow(sitesClim[which(is.na(sitesClim$ClimTmax.ave)==T),])

load("results/sitesClim_buffer50km_neotoma.Rdata")
nrow(sitesClimBuff[which(is.na(sitesClimBuff$ClimTmax.ave)==T),])

sitesClim[which(is.na(sitesClim$ClimTmax.ave)==T), c("ClimTmax.ave", "ClimPrec.ave")] <- 
  sitesClimBuff[which(is.na(sitesClim$ClimTmax.ave)==T), c("ClimTmax.ave", "ClimPrec.ave")]
nrow(sitesClim[which(is.na(sitesClim$ClimTmax.ave)==T),])


## 3.2 Get niceh bu joining climatic data and species data 
PaleoClimSpecies_site_time <- join(paleoOcc, sitesClim) %>% filter(occurence==1)
PaleoClimSpecies_site_time$species <- as.factor(PaleoClimSpecies_site_time$species)
PaleoClimSpecies_site_time <-PaleoClimSpecies_site_time%>% filter(occurence==1)
unique(PaleoClimSpecies_site_time$species)

save(PaleoClimSpecies_site_time, file="results/PaleoClimSpecies_site_time.Rdata")

load("results/PaleoClimSpecies_site_time.Rdata")

write.csv(PaleoClimSpecies_site_time, file="results/niches_paleo_all_raw.csv")

Paleo_niches_summary <- ddply(na.omit(PaleoClimSpecies_site_time), .(species), summarise, 
                        mean.Monthly.maxTemp = mean(ClimTmax.ave),
                        sd.Monthly.maxTemp = sd(ClimTmax.ave),
                        mean.Monthly.Prec = mean(ClimPrec.ave),
                        sd.Monthly.Prec = sd(ClimPrec.ave), .progress="text")

write.csv(Paleo_niches_summary, file="results/niches_paleo_summary.csv")

hist(PaleoClimSpecies_site_time$ClimTmax.ave)
hist(PaleoClimSpecies_site_time$ClimPrec.ave)






nRowSampling<- ddply(PaleoClimSpecies_site_time, .(species, time), function(x){nrow(x)})
hist(nRowSampling$V1)
mean(nRowSampling$V1)


nsamples=max(nRowSampling$V1)
PaleoClimSpecies_site_time_SAMPLED <- ddply(PaleoClimSpecies_site_time, .(species, time),  function(x) {
    # x<- PaleoClimSpecies_site_time[PaleoClimSpecies_site_time$taxon=="Abies" & PaleoClimSpecies_site_time$time==0, ]
    x[sample(1:nrow(x),nsamples,replace=T),c("ClimTmax.ave", "ClimPrec.ave","time")]
}, .progress="text")


write.csv(PaleoClimSpecies_site_time_SAMPLED, file="results/niches_paleo_all_sampled.csv")


Paleo_niches_summary_SAMPLED<-ddply(na.omit(PaleoClimSpecies_site_time_SAMPLED), .(species), summarise, 
                                    mean.Monthly.maxTemp = mean(ClimTmax.ave),
                                    sd.Monthly.maxTemp = sd(ClimTmax.ave),
                                    mean.Monthly.Prec = mean(ClimPrec.ave),
                                    sd.Monthly.Prec = sd(ClimPrec.ave), .progress="text")

plot(Paleo_niches_summary$mean.Monthly.maxTemp, Paleo_niches_summary_SAMPLED$mean.Monthly.maxTemp)

write.csv(Paleo_niches_summary_SAMPLED, file="results/niches_paleo_summary_sampled.csv")



