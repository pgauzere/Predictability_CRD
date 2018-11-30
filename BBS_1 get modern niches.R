###############################################################################################
############# 1. Estimate modern pooled niches of Birds  ######################################
#################################################### Gauzere & Blonder ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS/")

library(rgdal)
library(tidyverse)
library(raster)

### 1.read distribution range of BBS species from Birdlife Datas

#create the taxalist to subset the Birdlife geodatabase in QGIS 
load("data/dataframe_comtrue_cape100.Rdata")
head(dfCom)
taxalist_BirdlLife <- gsub("_", " ", as.character(unique(dfCom$species)))
write.table(taxalist_BirdlLife, file="BBS_Birdlife_taaxalist.txt", sep=",", row.names = F)      

# In QGIS, open the birdlife geodatabse, filter as follows : # "PRESENCE"=1  AND "SEASONAL"  IN (1,2)  AND "SCINAME" IN ( "BBS_Birdlife_taaxalist.txt" ) ##copy-paste the text direclty
# Then save the filtered file as ESRI "birdlife_BBS_breeding"

# Read the Birdlife distribution ranges of species in BBS data
ranges <- readOGR(dsn="data/birdlife_BBS_breeding/", layer='birdlife_BBS_breeding')
rangesTaxalist <- as.character(unique(ranges$SCINAME))

rangesDF<- ranges@data[, c(1,2,3,5,6,7)]
rangesDF$breedDistrib<-'YES'
rangesDF %>% group_by(SCINAME) %>% summarise(length(SEASONAL))

rangesDF[rangesDF$SCINAME == 'Alectoris chukar', ]

length(unique(ranges$SCINAME))
unique(ranges$SCINAME)
unique(ranges$SEASONAL)
unique(ranges$DATE_)

# 
# nrow(unique(ranges@data[,c(2,7)]))
# nrow(unique(ranges@data[,c(2,3,7)]))
# 
# #subset breeding rangesE
# breedRanges <- ranges[ranges$SEASONAL==2,]
# breedRangesTaxalist <- as.character(breedRanges$SCINAME)
# 
# #subset full-year resident ranges
# residRanges <- ranges[ranges$SEASONAL==1,]
# residRangesTaxalist <- unique(residRanges$SCINAME)
# proj4string(residRanges)

### 2 get Worldclim raster Data
worldclimTmaxRaster <- getData('worldclim', var='tmax', res=2.5)
breedTmax <- worldclimTmaxRaster[[5:7]] %>% raster::mean() #May, June, July
# yearTmax  <- worldclimTmaxRaster        %>% raster::mean()

worldclimTminRaster <- getData('worldclim', var='tmin', res=2.5)
breedTmin <- worldclimTminRaster[[5:7]] %>% raster::mean() #May, June, July
# yearTmax  <- worldclimTmaxRaster        %>% raster::mean()

worldclimPrecRaster <- getData('worldclim', var='prec', res=2.5) 
breedPrec <-worldclimPrecRaster[[5:7]] %>% raster::mean() #May, June, July
# yearPrec <- worldclimPrecRaster        %>% raster::mean()

climStack <- stack(c(breedTmax, breedTmin, breedPrec))

# save(worldclimTmaxRaster, file="data/worldclimTmaxRaster.Rdata")
# save(worldclimPrecRaster, file="data/worldclimPrecRaster.Rdata")
# mycrs <- worldclimTmaxRaster@crs@projargs
rm(breedPrec)
rm(breedTmax)
rm(worldclimPrecRaster)
rm(worldclimTmaxRaster)
rm(rangesDF)
gc()

#### 3. Extract Climatic niche from Wordclim to distribution ranges from Birdlife ####
birdsBreedClimNiches <- data.frame('species'=rangesTaxalist, 
                                   'mean.Monthly.maxTemp'=999,
                                   'mean.Monthly.minTemp'=999,
                                   'mean.Monthly.Prec'=999 ,
                                   'sd.Monthly.maxTemp'=999 ,
                                   'sd.Monthly.minTemp'=999 ,
                                   'sd.Monthly.Prec'=999)

for(i in 1:160){
  print(paste(i, rangesTaxalist[i]))
  spRange<- aggregate(ranges[ranges$SCINAME== rangesTaxalist[i],])
  birdsBreedClimNiches[i, 2:4] <- raster::extract(climStack, spRange, fun=mean, na.rm=T)
  birdsBreedClimNiches[i, 5:7] <- raster::extract(climStack, spRange, fun=sd, na.rm=T)
  gc()
}                          


for(i in 161:395){
  print(paste(i, rangesTaxalist[i]))
  spRange<- aggregate(ranges[ranges$SCINAME== rangesTaxalist[i],])
  birdsBreedClimNiches[i, 2:4] <- raster::extract(climStack, spRange, fun=mean, na.rm=T)
  birdsBreedClimNiches[i, 5:7] <- raster::extract(climStack, spRange, fun=sd, na.rm=T)
  gc()
}                          

birdsBreedClimNiches <- birdsBreedClimNiches[-160,]


hist(birdsBreedClimNiches$mean.Monthly.maxTemp)
hist(birdsBreedClimNiches$mean.Monthly.Prec)


ggplot(birdsBreedClimNiches, aes(y=mean.Monthly.maxTemp, x=mean.Monthly.Prec))+
  geom_point()+
  geom_errorbar(aes(ymin= mean.Monthly.maxTemp - sd.Monthly.maxTemp, ymax= mean.Monthly.maxTemp + sd.Monthly.maxTemp), alpha=0.2)+
  geom_errorbarh(aes(xmin= mean.Monthly.Prec - sd.Monthly.Prec, xmax= mean.Monthly.Prec + sd.Monthly.Prec), alpha=0.2)+
  theme_classic()

birdsBreedClimNiches$mean.Monthly.maxTemp <- birdsBreedClimNiches$mean.Monthly.maxTemp/10
birdsBreedClimNiches$sd.Monthly.maxTemp <- birdsBreedClimNiches$sd.Monthly.maxTemp/10

birdsBreedClimNiches <- birdsBreedClimNiches[-160,]#delete common murre which doesnt work


hist(birdsBreedClimNiches$mean.Monthly.maxTemp)
birdsBreedClimNiches[birdsBreedClimNiches$mean.Monthly.maxTemp < 16,]
birdsBreedClimNiches[birdsBreedClimNiches$mean.Monthly.maxTemp > 30,]


save(birdsBreedClimNiches , file='results/birdsBreedClimNiches.Rdata')




