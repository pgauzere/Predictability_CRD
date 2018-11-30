###############################################################################################
############# 6. Estimate modern pooled niches of plants  #####################################
#################################################### Gauzere & Blonder ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS /")

library(raster)
library(maps)
library(BIEN)
library(plyr)

load("data/filledXBackTrans_spmeans.Rdata")

#### 1.Create the list of taxa from pollen Neotoma datas ####
#neotoma data presence absence directory
dir.presabs <- 'data/blois data/output_for_BenBlonder/all data by site'
filenames.presabs <-dir(dir.presabs,pattern='*PA\\.csv',full.names=T)

#create a list wih elements as sites, containing presence/absence of trees
files.presabs <- filenames.presabs %>% plyr::llply(read.csv) 

#renames list elements with site names
names(files.presabs) <- gsub(dir.presabs,"",filenames.presabs,fixed=T) 
names(files.presabs) <- gsub("/","",names(files.presabs),fixed=T)
names(files.presabs) <- gsub(".interp.data.PA.csv","",names(files.presabs),fixed=T)

#transform list to DF
presAbs <- files.presabs %>% plyr::ldply(add_rownames)
colnames(presAbs)[1:2] <- c("site", "time")
presAbs$site <- as.factor(presAbs$site)
presAbs$time <- as.numeric(presAbs$time)

paleoOcc <- presAbs %>% gather("taxon", "occurence", 3:105) %>% filter(occurence == 1)


alltaxa <- 
  presAbs[-c(1,2)] %>% colnames() %>% 
  gsub(pattern=".",replacement=" ", fixed=T) %>% 
  gsub(pattern=" type",replacement="", fixed=T)

taxondf <- 
  as.data.frame(alltaxa) %>%
    separate(alltaxa, into=c("genus", "species"), fill="right", remove=F)%>%
    mutate(ref=1)

taxondf$ref[taxondf$species != "NA"] <- 2

#### 2.Extract traits of the taxa list from BIEN ####

#check available traits
BIEN_trait_list() 

#choose traits of interest
traitsOfInterest <- c("whole plant height","seed length", "seed mass","whole plant dispersal syndrome")

#extract traits of interest for each taxa
traits_bien <- 
  plyr::alply(taxondf, 1, function(x){
    if(x$ref==1){
      BIEN_trait_traitbygenus(genus=x$genus, trait=traitsOfInterest) }
    else{
      BIEN_trait_traitbyspecies(x$species, traitsOfInterest)}},
    .progress = "text")
names(traits_bien) <- alltaxa

#transform list to DF
traitsBienDF <- ldply(traits_bien)
traitsBienDF$trait_name <- as.factor(traitsBienDF$trait_name)
traitsBienDF$trait_value <- as.numeric(traitsBienDF$trait_value)

datasCheck <- ddply(traitsBienDF, .(alltaxa, trait_name),  summarize, "numberOfValues"=length(trait_value)) %>% spread(trait_name, numberOfValues)
## only seed mass and plant heights traits seems to have enough data

# traitsBienDF <- traitsBienDF[traitsBienDF$trait_name %in% traitsOfInterest,][,c(1,2,3,7,8,9)]

ggplot(traitsBienDF[traitsBienDF$genus=="Abies",], aes(x=trait_value))+
  geom_histogram()+
  facet_wrap(~trait_name)

pdf(file="results/species_traits_distribution.pdf", 16,9)
d_ply(traitsBienDF[traitsBienDF$trait_name %in% c("seed mass", "whole plant height"),], .(alltaxa), function(x){
  # x<- traitsBienDF[traitsBienDF$alltaxa=="Abies",]
  print(x$alltaxa[1])
print(
  ggplot(x, 
       aes(x=trait_value))+
  geom_histogram()+ 
  theme_classic()+
  labs(title=as.character(x$alltaxa[1]))+
    facet_wrap(~trait_name)
)
}, .progress="text")
dev.off()  


save(traitsBienDF,file='results/traitsBienDF.Rdata')

#### Compute CWM of traits of interests ####
load("results/traitsBienDF.Rdata")

test <- paleoOcc %>% filter(site=="3PINES")
x <- paleoOcc %>% filter(site=="3PINES")
# test <- niches_bien[niches_bien$taxon %in% x$taxon, c("Tmax","Prec", "taxon")]

CommTraits_distribs <- ddply(paleoOcc, .(site, time), function(x){
  ddply(traitsBienDF[traitsBienDF$alltaxa %in% x$taxon,],
        .(trait_name),summarize, "cwm"=mean(trait_value, na.rm=T), "cwv"=sd(trait_value, na.rm=T))}, .progress="text")

save(CommTraits_distribs, file="results/paleo_community_traits_distributions.Rdata")

load("results/paleo_community_traits_distributions.Rdata")


## only seed mass and plant heights traits seems to have enough data

CommTraits_distribs_cwm_wide <- CommTraits_distribs %>% spread(trait_name, cwm)
colnames(CommTraits_distribs_cwm_wide)[4:6] <- c("cwm_seed_length", "cwm_seed_mass", "cwm_whole_plant_height")

CommTraits_distribs_cwv_wide <- CommTraits_distribs %>% spread(trait_name, cwv)
colnames(CommTraits_distribs_cwv_wide)[4:6] <- c("cwv_seed_length", "cwv_seed_mass", "cwv_whole_plant_height")



# paleoCommClim_means   <- ddply(paleoOcc, .(site, time), summarize, 
#                                CommTmax.ave = mean(niches_bien_summary$mean.Monthly.maxTemp[niches_bien_summary$alltaxa %in% taxon]),
#                                CommTmax.sd  = sd(niches_bien_summary$mean.Monthly.maxTemp[niches_bien_summary$alltaxa %in% taxon]),
#                                CommPrec.ave = mean(niches_bien_summary$mean.Monthly.Prec[niches_bien_summary$alltaxa %in% taxon]),
#                                CommPrec.sd  = sd(niches_bien_summary$mean.Monthly.Prec[niches_bien_summary$alltaxa %in% taxon]),
#                                .progress="text")    
# save(paleoCommClim_means, file="results/paleo_community_climato_mean.Rdata")


library(viridis)


ggplot(CommTraits_distribs, aes(y=cwm, x=cwv))+
  geom_bin2d()+ 
  theme_classic()+
  # labs(title=as.character(x$alltaxa[1]))+
  facet_wrap(~trait_name)+
  scale_fill_viridis()+
  scale_x_continuous(trans="log")+
  scale_y_continuous(trans="log")


ggplot(CommTraits_distribs, aes(y=cwm, x=cwv))+
  geom_bin2d()+ 
  theme_classic()+
  # labs(title=as.character(x$alltaxa[1]))+
  facet_wrap(~trait_name)+
  scale_fill_viridis()+
  scale_x_continuous(trans="log")+
  scale_y_continuous(trans="log")


ggplot(CommTraits_distribs_cwm_wide, aes(y=cwm_seed_mass, x=cwm_whole_plant_height))+
  geom_point()+ 
  theme_classic()+
  # labs(title=as.character(x$alltaxa[1]))+
  # facet_wrap(~trait_name)+
  scale_fill_viridis()+
  scale_x_continuous(trans="log")+
  scale_y_continuous(trans="log")

CommTraits_distribs_cwm_wide


SitesTraits_distribs <- ddply(CommTraits_distribs, .(site), function(x){
  ddply(traitsBienDF[traitsBienDF$alltaxa %in% x$taxon,],
        .(trait_name),summarize, "cwm"=mean(trait_value, na.rm=T), "cwv"=sd(trait_value, na.rm=T))}, .progress="text")

save(SitesTraits_distribs, file="results/sites_traits_distributions.Rdata")

load("results/sites_traits_distributions.Rdata")
load("results/statistics_sites_neotoma.Rdata")

CWMs <- SitesTraits_distribs %>% select(-cwv) %>% spread(trait_name, cwm) 
colnames(CWMs)[2:4] <- c("mean.seed.length", "mean.seed.mass", "mean.whole.plant.height") 

CWVs <- SitesTraits_distribs %>% select(-cwm) %>% spread(trait_name, cwv) 
colnames(CWVs)[2:4] <- c("sd.seed.length", "sd.seed.mass", "sd.whole.plant.height") 

save(CWMs, file="results/biotic_predictors_CWMs.Rdata")

ggplot(CWMs, aes(mean.seed.mass,mean.whole.plant.height))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_viridis()

ggplot(CWVs, aes(sd.seed.mass,sd.whole.plant.height))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_viridis()


test <- join(statistics_sites, CWMs)
test <- join(test,  CWVs)


# World map, using geom_path instead of geom_polygon
world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,60), xlim=c(-115, -55))+
  theme(legend.position="bottom")

worldmap + 
  geom_point(data=test, aes(Longitude, Latitude, fill=mean.whole.plant.height),size=2, shape=21)+
  scale_fill_viridis(trans="sqrt")+
  theme_classic()+
  # scale_fill_gradientn(colours)cubeHelix(100)+
  theme(legend.position=c(0.85,0.2))

worldmap + 
  geom_point(data=test, aes(Longitude, Latitude, fill=mean.seed.length),size=2, shape=21)+
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

[test$mean.seed.mass>0.1 & test$mean.seed.mass<100]


ggplot(test[test$mean.seed.mass<150,], aes(x=mean.seed.mass,y=deviation.absolute.Tmax, fill=deviation.absolute.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log1p")+
  # scale_x_continuous(trans="log1p")+
  scale_fill_viridis()+
  # stat_smooth(method="lm")+
  theme_classic()


ggplot(test[test$mean.whole.plant.height<8,], aes(x=mean.whole.plant.height, y=deviation.base.Tmax, fill=deviation.absolute.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log1p")+
  # scale_x_continuous(trans="log1p")+
  scale_fill_viridis()+
  stat_smooth(method="lm")+
  theme_classic()

ggplot(test[test$mean.seed.mass<150,], aes(x=mean.seed.mass,y=mean.whole.plant.height, fill=deviation.absolute.Tmax))+
  geom_jitter(shape=21, size=3)+
  scale_y_continuous(trans="log")+
  # scale_x_continuous(trans="log")+
  scale_fill_viridis()+
  # stat_smooth(method="lm")+
  theme_classic()

ggplot(test[test$mean.seed.mass<150,], aes(x=mean.seed.mass,y=mean.seed.length, fill=deviation.absolute.Tmax))+
  geom_jitter(shape=21, size=3)+
  scale_y_continuous(trans="log")+
  scale_y_continuous(trans="log")+
  scale_fill_viridis()+
  # stat_smooth(method="lm")+
  theme_classic()


ggplot(test[test$mean.seed.mass<150,], aes(y=crossings.mean.Tmax,x=mean.seed.mass, fill=crossings.mean.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log")+
  # scale_x_continuous(trans="log")+
  scale_fill_viridis()+
  # stat_smooth(method="lm")+
  theme_classic()

ggplot(test[test$mean.seed.mass<150,], aes(y=crossings.mean.Tmax, x=mean.seed.length, fill=crossings.mean.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log")+
  # scale_x_continuous(trans="log")+
  scale_fill_viridis()+
  # stat_smooth(method="lm")+
  theme_classic()


ggplot(test[test$mean.whole.plant.height<8,], aes(y=crossings.mean.Tmax, x=mean.whole.plant.height, fill=crossings.mean.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log")+
  scale_x_continuous(trans="log")+
  scale_fill_viridis()+
  # stat_smooth(method="lm")+
  theme_classic()

ggplot(test[test$mean.seed.mass<150,], aes(x=mean.seed.mass,y=lm.rho.Tmax, fill=lm.rho.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log1p")+
  # scale_x_continuous(trans="log1p")+
  scale_fill_viridis()+
  stat_smooth(method="lm")+
  theme_classic()

ggplot(test, aes(x=mean.seed.length,y=lm.rho.Tmax, fill=lm.rho.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log1p")+
  # scale_x_continuous(trans="log1p")+
  scale_fill_viridis()+
  stat_smooth(method="lm")+
  theme_classic()

ggplot(test[test$mean.whole.plant.height<8,], aes(x=mean.whole.plant.height,y=lm.rho.Tmax, fill=lm.rho.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log1p")+
  # scale_x_continuous(trans="log1p")+
  scale_fill_viridis()+
  stat_smooth(method="lm")+
  theme_classic()

ggplot(test[test$mean.whole.plant.height<8,], aes(x=mean.whole.plant.height, y=crossings.max.Tmax, fill=deviation.absolute.Tmax))+
  geom_jitter(shape=21, size=3)+
  # scale_y_continuous(trans="log1p")+
  # scale_x_continuous(trans="log1p")+
  scale_fill_viridis()+
  stat_smooth(method="lm")+
  theme_classic()



