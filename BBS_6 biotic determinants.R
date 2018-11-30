###############################################################################################
############# 6. Estimate modern pooled niches of plants  #####################################
#################################################### Gauzere & Blonder ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS/")


#### 1. read BBS species traits ####
traits <- read.table('data/BBS_species_traits.csv', h=T, sep=';')
#choose traits of interest
traitsOfInterest <- traits[,c(1,3,8)]
colnames(traitsOfInterest) <- c('species', 'body_mass','migration')

#### 2.load BBS data ####
load('data/dataframe_comtrue_occurence.Rdata')

#### 3. Compute Community meam Values ####

dfComOccTraits <- left_join(dfComOcc, traitsOfInterest)

BirdsTraits_means_weighted   <- plyr::ddply(na.omit(dfComOccTraits), .(route, year), function(x){
  BodyMass = sum(x$body_mass * x$N / sum(x$N))
  Migration= (sum(x$N[x$migration=='migrant'])/sum(x$N))*100
  return(cbind(BodyMass, Migration))}, .progress='text')


colnames(BirdsTraits_means_weighted)[1]<- 'site'
BirdsTraits_means_weighted$year <- as.numeric(as.character(BirdsTraits_means_weighted$year))

load('results/predictors_BBS.Rdata')

predictors <- merge(predictors, BirdsTraits_means_weighted, by=c('site', 'year'))

save(predictors, file="results/predictors_BBS.Rdata")


# World map, using geom_path instead of geom_polygon
world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,60), xlim=c(-115, -55))+
  theme(legend.position="bottom")

worldmap + 
  geom_point(data=test, aes(Longitude, Latitude, fill=Migration),size=2, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  # scale_fill_gradientn(colours)cubeHelix(100)+
  theme(legend.position=c(0.85,0.2))

worldmap + 
  geom_point(data=test, aes(Longitude, Latitude, fill=BodyMass),size=2, shape=21)+
  scale_fill_viridis(trans="log")+
  theme_classic()+
  # scale_fill_gradientn(colours)cubeHelix(100)+
  theme(legend.position=c(0.85,0.2))

worldmap + 
  geom_point(data=test, aes(Longitude, Latitude, fill=mean.seed.mass),size=2, shape=21)+
  scale_fill_viridis(trans="log10")+
  theme_classic()+
  # scale_fill_gradientn(colours)cubeHelix(100)+
  theme(legend.position=c(0.85,0.2))



