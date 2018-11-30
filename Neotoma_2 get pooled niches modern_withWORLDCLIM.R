###############################################################################################
############# 2. Estimate modern pooled niches of plants  #####################################
#################################################### Gauzere & Blonder ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS /")

library(raster)
library(maps)
library(BIEN)
library(dplyr)
library(plyr)
library(tidyr)


#### 1.Create the list of taxa from pollen Neotoma datas ####
#neotoma data presence absence directory
dir.presabs <- 'data/blois data/output_for_BenBlonder/all data by site'
filenames.presabs <-dir(dir.presabs,pattern='*PA\\.csv',full.names=T)

#create a list wih elements as sites, containing presence/absence of trees
files.presabs <- filenames.presabs %>% llply(read.csv) 

#renames list elements with site names
names(files.presabs) <- gsub(dir.presabs,"",filenames.presabs,fixed=T) 
names(files.presabs) <- gsub("/","",names(files.presabs),fixed=T)
names(files.presabs) <- gsub(".interp.data.PA.csv","",names(files.presabs),fixed=T)

#transform list to DF
presAbs <- files.presabs %>% ldply

#read data on sites (lat, lon, age, etc...)
latlon <- read.csv('data/blois data/pollen old/All.site.data-withagemodel-finalv2.csv')
climvars <- c('tmax_year_ave','prcp_year_ave')

alltaxa <- 
  presAbs[-c(1,2)] %>% colnames() %>% 
  gsub(pattern=".",replacement=" ", fixed=T) %>% 
  gsub(pattern=" type",replacement="", fixed=T)

taxondf <- 
  as.data.frame(alltaxa) %>%
    separate(alltaxa, into=c("genus", "species"), fill="right", remove=F)%>%
    mutate(ref=1)

taxondf$ref[taxondf$species != "NA"] <- 2

#### 2.Extract occurences of the taxa list from BIEN ####
occurences_all <- 
  alply(taxondf, 1, function(x){
    if(x$ref==1){
      BIEN_occurrence_genus(x$genus) }
    else{
      BIEN_occurrence_species(x$species)}},
    .progress = "text")
names(occurrences_all) <- alltaxa

save(occurrences_all,file='results/occurrences_bien.Rdata')

#transform list to DF
occurencesdf <- ldply(occurences_all)

#### 3.Merge Climatic data from Wordclim to spatialized occurrences from BIEN ####

## 3.1 get Worldclim Data
worldclimTmaxRaster <- getData('worldclim', var='tmax', res=2.5) %>% raster::mean()
worldclimPrecRaster <- getData('worldclim', var='prec', res=2.5) %>% raster::mean()
save(worldclimTmaxRaster, file="data/worldclimTmaxRaster.Rdata")
save(worldclimPrecRaster, file="data/worldclimPrecRaster.Rdata")
mycrs <- worldclimTmaxRaster@crs@projargs


## 3.2 get coordinates of occurrences
load("results/occurrences_DF_bien.Rdata")

occurences_bien_coordinates <- 
  distinct(occurencesdf[, c("latitude", "longitude")]) %>%
  na.omit() %>%
  filter(longitude < 0 & latitude < 80 )

#Convert points to spatial points data frame
coordinates(occurences_bien_coordinates) <- c('longitude', 'latitude')
#put on the same geo-referencial
proj4string(occurences_bien_coordinates) <- CRS(mycrs)

#Plot to check spatial match
plot(worldclimTmaxRaster, main="Average Tmax")
plot(occurences_bien_coordinates, add=T)

bien_worldclim <- data.frame(coordinates(occurences_bien_coordinates), 
                            raster::extract(worldclimTmaxRaster, occurences_bien_coordinates), 
                            raster::extract(worldclimPrecRaster, occurences_bien_coordinates))   #extract clim values

occurrences_bienWithWorldclim <- join(unique(occurencesdf[1:9], bien_worldclim)
colnames(occurrences_bienWithWorldclim)[10:11] <- c("TmaxWorldClim", "Prec")
occurrences_bienWithWorldclim$Tmax <- occurrences_bienWithWorldclim$TmaxWorldClim/10
occurrences_bienWithWorldclim$Prec <- occurrences_bienWithWorldclim$Prec*12

occurrences_bienWithWorldclim <- occurrences_bienWithWorldclim %>% distinct(alltaxa, latitude, longitude, Prec, Tmax)

saveRDS(occurrences_bienWithWorldclim, file="results/occurrences_bienWithWorldclim.Rdata")



#### 4.Explore and summarise modern niches ####
occurrences_bienWithWorldclim <- readRDS("results/occurrences_bienWithWorldclim.Rdata")

library(colorRamps)
library(rje)
library(ggplot2)
library(ggExtra)

niches_summary <- ddply(na.omit(occurrences_bienWithWorldclim), .(alltaxa), summarise, 
                        mean.Monthly.maxTemp = mean(Tmax),
                        sd.Monthly.maxTemp = sd(Tmax),
                        mean.Monthly.Prec = mean(Prec),
                        sd.Monthly.Prec = sd(Prec), .progress="text")

write.csv(niches_summary, 'results/niches_summary_bien_worldclim.csv',row.names=F)
niches_summary <- read.csv('results/niches_summary_bien_worldclim.csv')



x <- droplevels(occurrences_bienWithWorldclim %>% filter(alltaxa=="Abies"))

# World map, using geom_path instead of geom_polygon
world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey") +
  coord_fixed(ratio=1, ylim=c(-55,80), xlim=c(-170, -30))

pdf("results/modern_niche_maps.pdf")
d_ply(occurrences_bienWithWorldclim, .(alltaxa), function(x){

  print(
    worldmap+
    geom_point(data=x, aes(longitude, latitude,fill=Tmax), shape=21)+
    scale_fill_gradientn(colours=matlab.like2(50))+
    theme_classic()+
    theme(legend.justification=c(0,0), legend.position = c(0,0))+
    labs(title=as.character(x$alltaxa[1]), x="longitude", y="latitude")
  )
  }, .progress="text")
dev.off()

pdf("results/modern_niche_climaticSpace.pdf")
d_ply(na.omit(occurrences_bienWithWorldclim), .(alltaxa), function(x){
# p <- 
  
  print(
  ggplot(data=x, aes(x=Tmax,y=Prec)) + 
  # stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
  # stat_density2d(aes(fill = ..level..), geom = "polygon")+
  geom_point(shape=21, alpha=0)+
  geom_bin2d(bins=50)+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_gradientn(colours=cubeHelix(200))+
  geom_errorbar(data=niches_summary[niches_summary$alltaxa==as.character(x$alltaxa[1]),], aes(x=mean.Monthly.maxTemp,y=mean.Monthly.Prec,ymin=mean.Monthly.Prec   -sd.Monthly.Prec,   ymax=mean.Monthly.Prec + sd.Monthly.Prec),col="red", width=0)+
  geom_errorbarh(data=niches_summary[niches_summary$alltaxa==as.character(x$alltaxa[1]),],aes(y=mean.Monthly.Prec,x=mean.Monthly.maxTemp,xmin=mean.Monthly.maxTemp-sd.Monthly.maxTemp,xmax=mean.Monthly.maxTemp + sd.Monthly.maxTemp),col="red", width=0)+
  labs(title=as.character(x$alltaxa[1]), y="Yearly precipitation", x="Mean monthly maximum temperature")+
  lims(x=c(min(occurrences_bienWithWorldclim$Tmax,na.rm = T), max(occurrences_bienWithWorldclim$Tmax,na.rm = T)),
       y=c(min(occurrences_bienWithWorldclim$Prec,na.rm = T), max(occurrences_bienWithWorldclim$Prec,na.rm = T)))  
  )             
# ggMarginal(p, type = "histogram")
  }, .progress="text")
dev.off() 

