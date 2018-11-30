
###############################################################################################
############# 4. Compute community inferred climate and merge with paleo_climate  #############
#################################################### Gauzere & Blonder ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS")
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(raster)
library(broom)


#### get climate paleo time-series data ####

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
  select(1,2,3,4,5,51,52) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))

sitesprec <- sitesClim %>% gather(key="time", "ClimPrec.ave", 51:95) %>% 
  select(1,2,3,4,5,51,52) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))

sitesClim <- sitestmax
sitesClim$ClimPrec.ave <- sitesprec$ClimPrec.ave

colnames(sitesClim)[1:3] <-  c("SiteName","site", "altitude")
save(sitesClim, file="results/sitesClim_neotoma.Rdata")

load("results/sitesClim_neotoma.Rdata")
## extract mean and sd of tmax and prec within 50000m buffer(~ 0.5 degrees)

ClimTmax.ave <- raster::extract(TmaxRaster, latlon, buffer=50000,fun=mean,na.rm=T,small=T)
ClimPrec.ave <- raster::extract(PrecRaster, latlon, buffer=50000,fun=mean,na.rm=T,small=T)
ClimTmax.sd  <- raster::extract(TmaxRaster, latlon, buffer=50000,fun=sd,na.rm=T,small=T)
ClimPrec.sd  <- raster::extract(PrecRaster, latlon, buffer=50000,fun=sd,na.rm=T,small=T)
save(ClimTmax.ave, file="data/ClimTmax.ave.Buff.Rdata")
save(ClimPrec.ave, file="data/ClimPrec.ave.Buff.Rdata")
save(ClimTmax.sd, file="data/ClimTmax.sd.Buff.Rdata")
save(ClimPrec.sd, file="data/ClimPrec.sd.Buff.Rdata")

ClimTmax.ave <- as.data.frame(ClimTmax.ave) %>% mutate(site= latlon$Handle) %>% gather(key="time", "ClimTmax.ave",1:45) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))
ClimPrec.ave <- as.data.frame(ClimPrec.ave) %>% mutate(site= latlon$Handle) %>% gather(key="time", "ClimPrec.ave",1:45) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))
ClimTmax.sd <- as.data.frame(ClimTmax.sd) %>% mutate(site= latlon$Handle) %>% gather(key="time", "ClimTmax.sd",1:45) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))
ClimPrec.sd <- as.data.frame(ClimPrec.sd) %>% mutate(site= latlon$Handle) %>% gather(key="time", "ClimPrec.sd",1:45) %>% 
  mutate(time = as.numeric(gsub("X", "",time)))

colnames(latlon@data)[2] <- "site"
test <- merge(cbind(latlon@data,coordinates(latlon)), ClimTmax.ave, by="site")
test <- merge(test, ClimPrec.ave, by=c("site", "time"))
test <- merge(test, ClimTmax.sd, by=c("site", "time"))
test <- merge(test, ClimPrec.sd, by=c("site", "time"))
sitesClimBuff <- test
save(sitesClimBuff, file="results/sitesClim_buffer50km_neotoma.Rdata")

load("results/sitesClim_buffer50km_neotoma.Rdata")

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
save(paleoOcc, file="results/PaleoOcc.Rdata")

#### Compite Community Inferred Climate ####
load("results/PaleoOcc.Rdata")
# we decided to only use niches from comptemporary BIEN data
# niches_all_summary <- read.csv('results/niches_all_summary_sampled.csv')

niches_bien_summary <- read.csv('results/niches_summary_bien_worldclim.csv')


paleoCommClim_means   <- ddply(paleoOcc, .(site, time), summarize, 
                               CommTmax.ave = mean(niches_bien_summary$mean.Monthly.maxTemp[niches_bien_summary$alltaxa %in% taxon]),
                               CommTmax.sd  = sd(niches_bien_summary$mean.Monthly.maxTemp[niches_bien_summary$alltaxa %in% taxon]),
                               CommPrec.ave = mean(niches_bien_summary$mean.Monthly.Prec[niches_bien_summary$alltaxa %in% taxon]),
                               CommPrec.sd  = sd(niches_bien_summary$mean.Monthly.Prec[niches_bien_summary$alltaxa %in% taxon]),
                               .progress="text")    
save(paleoCommClim_means, file="results/paleo_community_climato_mean.Rdata")

load("results/paleo_community_climato_mean.Rdata")

# test <- paleoOcc %>% filter(site=="3PINES")
# x <- paleoOcc %>% filter(site=="ZUEHL")
# test <- niches_bien[niches_bien$taxon %in% x$taxon, c("Tmax","Prec", "taxon")]
# hist(niches_bien$Tmax)
# 
# paleoCommClim_distribs <- ddply(paleoOcc, .(site, time), summarize, 
#                                 CommTmax.ave = mean(niches_all$Tmax[niches_all$taxon %in% taxon]),
#                                 CommTmax.sd  = sd(niches_all$Tmax[niches_all$taxon %in% taxon]),
#                                 CommPrec.ave = mean(niches_all$Prec[niches_all$taxon %in% taxon]),
#                                 CommPrec.sd  = sd(niches_all$Prec[niches_all$taxon %in% taxon]),
#                                 .progress="text")    
# save(paleoCommClim_distribs, file="results/paleo_community_climato_distributions.Rdata")

# plot(paleoCommClim_distribs$CommTmax.ave, paleoCommClim_means$CommTmax.ave)
# abline(0,1)

#### join community and climate datasets ####
load("results/paleo_community_climato_mean.Rdata")
load("results/sitesClim_buffer50km_neotoma.Rdata")
# CommClim <- merge(paleoCommClim_distribs, sitesClimBuff, by=c("site", "time"))
CommClim <- merge(paleoCommClim_means, sitesClimBuff, by=c("site", "time"))

## we only keep site for which we have more than 5 time steps
sitesToKeep <-CommClim %>% na.omit() %>% group_by(site) %>% dplyr::summarize(n = length(unique(time)))   %>% filter(n>5) %>% droplevels
CommClim2 <- CommClim %>% na.omit() %>% filter(site %in% sitesToKeep$site) %>% droplevels() %>% arrange(site, desc(time))

rm("CommClim")
rm("paleoCommClim_mean") 
rm("sitesClim")
rm("sitesToKeep")


#### smooth Ct and Ft ####
## loess independently for Ct and Ft

# CommClim2 <- ddply(CommClim2, .(site),mutate, "CommTmax.ave_loess"=augment(loess(CommTmax.ave ~ time))[[3]], .progress = "text")
# CommClim2 <- ddply(CommClim2, .(site),mutate, "CommTmax.ave_loess_se"=augment(loess(CommTmax.ave ~ time))[[4]], .progress = "text")
# 
# CommClim2 <- ddply(CommClim2, .(site),mutate, "ClimTmax.ave_loess"=augment(loess(ClimTmax.ave ~ time))[[3]], .progress = "text")
# CommClim2 <- ddply(CommClim2, .(site),mutate, "ClimTmax.ave_loess_se"=augment(loess(ClimTmax.ave ~ time))[[4]], .progress = "text")
# 
# CommClim2 <- ddply(CommClim2, .(site),mutate, "CommPrec.ave_loess"=augment(loess(CommPrec.ave ~ time))[[3]], .progress = "text")
# CommClim2 <- ddply(CommClim2, .(site),mutate, "CommPrec.ave_loess_se"=augment(loess(CommPrec.ave ~ time))[[4]], .progress = "text")
# 
# CommClim2 <- ddply(CommClim2, .(site),mutate, "ClimPrec.ave_loess"=augment(loess(ClimPrec.ave ~ time))[[3]], .progress = "text")
# CommClim2 <- ddply(CommClim2, .(site),mutate, "ClimPrec.ave_loess_se"=augment(loess(ClimPrec.ave ~ time))[[4]], .progress = "text")


CommClim2<- CommClim2 %>% group_by(site) %>% mutate(CommTmax.ave_loess = augment(loess(CommTmax.ave ~ time))[[3]],
                                                CommTmax.ave_loess_se=augment(loess(CommTmax.ave ~ time))[[4]],
                                                ClimTmax.ave_loess=augment(loess(ClimTmax.ave ~ time))[[3]],
                                                ClimTmax.ave_loess_se=augment(loess(ClimTmax.ave ~ time))[[4]],
                                                CommPrec.ave_loess=augment(loess(CommPrec.ave ~ time))[[3]],
                                                CommPrec.ave_loess_se=augment(loess(CommPrec.ave ~ time))[[4]],
                                                ClimPrec.ave_loess=augment(loess(ClimPrec.ave ~ time))[[3]],
                                                ClimPrec.ave_loess_se=augment(loess(ClimPrec.ave ~ time))[[4]])
                                                

## or bezier curve on phase diagram
# CommClim2 <- ddply(CommClim2, .(site),mutate, 
#                    "CommTmax.ave_bezier"=bezier(t=seq(0,1,len=length(x$CommTmax.ave)), p=x[,c("ClimTmax.ave", "CommTmax.ave")], deg=(length(x$CommTmax.ave)-1))[,2], .progress = "text")
# 
# CommClim2 <- ddply(CommClim2, .(site),mutate, 
#                    "ClimTmax.ave_bezier"=bezier(t=seq(0,1,len=length(x$CommTmax.ave)), p=x[,c("ClimTmax.ave", "CommTmax.ave")], deg=(length(x$CommTmax.ave)-1))[,1], .progress = "text")
# 
# CommClim2 <- ddply(CommClim2, .(site),mutate, 
#                    "CommPrec.ave_bezier"=bezier(t=seq(0,1,len=length(x$CommPrec.ave)), p=x[,c("ClimPrec.ave", "CommPrec.ave")], deg=(length(x$CommPrec.ave)-1))[,2], .progress = "text")
# 
# CommClim2 <- ddply(CommClim2, .(site),mutate, 
#                    "ClimPrec.ave_bezier"=bezier(t=seq(0,1,len=length(x$CommPrec.ave)), p=x[,c("ClimPrec.ave", "CommPrec.ave")], deg=(length(x$CommPrec.ave)-1))[,1], .progress = "text")

#### compute summary statistics ####
source("scripts/PG_Summary Statistics_2.0.R")
library(plyr)
statistics_sites_Tmax <-ddply(CommClim2, .(site, SiteName, Altitude, Longitude, Latitude),function(x){
  summarizeftct(ft=x$ClimTmax.ave_loess, ct=x$CommTmax.ave_loess, ft_error = x$ClimTmax.ave_loess_se, ct_error=x$CommTmax.ave_loess_se, time=x$time, name=x$site)}, .progress = "text") %>% distinct()

colnames(statistics_sites_Tmax)[7:18] <- paste(colnames(statistics_sites_Tmax)[7:18], "Tmax", sep=".")

statistics_sites_Prec <-ddply(CommClim2, .(site, SiteName, Altitude, Longitude, Latitude),function(x){
  summarizeftct(ft=x$ClimPrec.ave_loess, ct=x$CommPrec.ave_loess, ft_error = x$ClimPrec.ave_loess_se, ct_error=x$CommPrec.ave_loess_se, time=x$time, name=x$site)}, .progress = "text") %>% distinct()
colnames(statistics_sites_Prec)[7:18] <- paste(colnames(statistics_sites_Prec)[7:18], "Prec", sep=".")

statistics_sites <- join(statistics_sites_Tmax, statistics_sites_Prec)

ggplot(statistics_sites, aes(deviation.base.Tmax, deviation.absolute.Tmax))+
  geom_point()

## with scales values 
# CommClim2$scaled_ClimTmax.ave <- scale(CommClim2$ClimTmax.ave)
# CommClim2$scaled_CommTmax.ave <- scale(CommClim2$CommTmax.ave)
# CommClim2$scaled_ClimPrec.ave <- scale(CommClim2$ClimPrec.ave)
# CommClim2$scaled_CommPrec.ave <- scale(CommClim2$CommPrec.ave)
# 
# source("scripts/saved scripts/Data S3 - Summary Statistics.R")
# statistics_sites_Tmax <-ddply(CommClim2, .(site, SiteName, Altitude, Longitude, Latitude),function(x){
#   summarizeftct(ft=x$scaled_ClimTmax.ave, ct=x$scaled_CommTmax.ave, name=x$site)}, .progress = "text") %>% distinct()
# colnames(statistics_sites_Tmax)[7:15] <- paste(colnames(statistics_sites_Tmax)[7:15], "Tmax", sep=".")
# 
# statistics_sites_Prec <-ddply(CommClim2, .(site, SiteName, Altitude, Longitude, Latitude),function(x){
#   summarizeftct(ft=x$scaled_ClimPrec.ave, ct=x$scaled_CommPrec.ave, name=x$site)}, .progress = "text") %>% distinct()
# colnames(statistics_sites_Prec)[7:15] <- paste(colnames(statistics_sites_Prec)[7:15], "Prec", sep=".")

##long format
statistics_long <- gather(statistics_sites, key="statistic", value="value", 7:26, factor_key = TRUE)
save(statistics_long, file="results/statistics_sites_long_neotoma.Rdata")

test <- left_join(CommClim2, statistics_sites_Tmax)
CommClimStat2 <- left_join(test, statistics_sites_Prec)


# plot(x=CommClimStat$crossings.mean.Tmax, y=CommClimStat2$crossings.mean.Tmax)
# plot(x=CommClimStat$crossings.max.Tmax, y=CommClimStat2$crossings.max.Tmax)

test <- CommClimStat2%>% group_by(site) %>% dplyr::summarise(nTime=length(unique(time)))
hist(test$nTime, breaks=45)
min(test$nTime)
CommClimStat2 <- CommClimStat2 %>% arrange(site, desc(time))
# CommClimStat <- droplevels(CommClimStat)

#add t0-first time step and tmax-last time step to commclimstat
period_study <-CommClimStat2 %>% group_by(site) %>% dplyr::summarize(t0 = min(-time), tmax = max(-time), timespan = max(-time) - min(-time))
hist(period_study$t0, breaks=50)
hist(period_study$tmax, breaks=50)
hist(period_study$timespan, breaks=50)
plot(period_study$t0, period_study$timespan)
CommClimStat2 <- left_join(CommClimStat2, period_study) %>% arrange(t0)
sites <- unique(CommClimStat2$site)

save(CommClimStat2, file="results/CommClimStat2_smoothed.Rdata")

# rm(CommClim2)
statistics_sites <- unique(CommClimStat2[, -c(2,3,4,5,6,11:23)])
save(statistics_sites, file="results/statistics_sites_neotoma.Rdata")
