setwd("~/PIERRE/TRANSIENT DYNAMICS")
library(tidyverse)

traits <- readRDS("data/filledXBackTrans_spmeans.Rdata")

##pairs correlation of traits
library("PerformanceAnalytics")
chart.Correlation(traits[, c("X.PlantHeight.","X.SeedMass.","X.Seed.length.", "X.Dispersal.unit.length.")]
                  , histogram=TRUE)


#### 1.Create the list of taxa used in previous analysis ####
taxalist <- read.csv("results/niches_all_summary.csv")[,1]

#### 2.Extract Genus-level traits values  ####
subtraits <- traits %>% 
  filter(Genus %in% taxalist | Genus %in% c('Juniperus', 'Ostrya', 'Rumex') | AccSpeciesName %in% c('Juniperus thuja', 'Ostrya carpinus', 'Rumex oxyria')) %>%
  select(AccSpeciesName, Genus, X.PlantHeight., X.SeedMass., X.Seed.length.)%>% 
  group_by(Genus) %>%
  summarise(PlantHeight = mean(X.PlantHeight.), SeedMass = mean(X.SeedMass.), SeedLength = mean(X.Seed.length.),  
            PlantHeight.sd = sd(X.PlantHeight.), SeedMass.sd = sd(X.SeedMass.), SeedLength.sd = sd(X.Seed.length.))

ggplot(subtraits, aes(y=PlantHeight, x=SeedMass, col=SeedLength))+
  geom_pointrange(aes(ymin=PlantHeight-PlantHeight.sd, ymax=PlantHeight+PlantHeight.sd))+
  geom_errorbarh(aes(xmin=SeedMass-SeedMass.sd, xmax=SeedMass+SeedMass.sd))+
  theme_classic()


#### 3. Compute species climate tolerance as niche breath  ####
load('results/niches_all_sampled.Rdata')

nicheRange <- niches_all_sampled %>% group_by(taxon) %>% summarise(range.Tmax = sd(Tmax),
                                                           range.Prec = sd(Prec))

nicheRange$Genus <- as.character(nicheRange$taxon)

nicheRange[nicheRange$taxon == 'Juniperus.Thuja', 'Genus'] <- 'Juniperus'
nicheRange[nicheRange$taxon == 'Ostrya.Carpinus', "Genus"] <- 'Ostrya'
nicheRange[nicheRange$taxon == 'Rumex.Oxyria', "Genus"] <- 'Rumex'
nicheRange[nicheRange$taxon == 'Planera.aquatica', "Genus"] <- 'Planera'
nicheRange[nicheRange$taxon == 'Chamaedaphne.calyculata', "Genus"] <- 'Chamaedaphne'



# test2 <- niches_all %>% group_by(taxon) %>% summarise(range.Tmax = 2 * 1.96 * (sd(Tmax)  /sqrt(length(Tmax))),
#                                                      range.Prec = 2 * 1.96 * (sd(.$Prec)/sqrt(length(Tmax))))
# 
# test3 <- niches_all %>% group_by(taxon) %>% summarise(range.Tmax = max(Tmax) - min(Tmax),
#                                                      range.Prec = max(Prec) - min(Prec))



#### 4. Compute CWM of traits of interests for each neotoma community ####
load("results/PaleoOcc.Rdata")
paleoOcc$Genus <- paleoOcc$taxon
paleoOcc[paleoOcc$taxon == 'Juniperus Thuja', "Genus"] <- 'Juniperus'
paleoOcc[paleoOcc$taxon == 'Ostrya Carpinus', "Genus"] <- 'Ostrya'
paleoOcc[paleoOcc$taxon == 'Rumex Oxyria', "Genus"] <- 'Rumex'
paleoOcc[paleoOcc$taxon == 'Planera aquatica', "Genus"] <- 'Planera'
paleoOcc[paleoOcc$taxon == 'Chamaedaphne calyculata', "Genus"] <- 'Chamaedaphne'

sort(unique(paleoOcc$Genus))

cbind(sort(unique(paleoOcc$Genus)), sort(unique(nicheRange$Genus)))

test <- inner_join(paleoOcc, nicheRange) 
paleoOccTraits <- left_join(test, subtraits)

paleoOccTraits[which(is.na(paleoOccTraits$range.Tmax)==T),]

sitesCWMs <- paleoOccTraits %>% group_by(site) %>% summarise(
                                                      mean.PlantHeight =mean(PlantHeight), 
                                                      mean.SeedMass    =mean(SeedMass),
                                                      mean.SeedLength  =mean(SeedLength),
                                                      mean.rangeTmax  =mean(range.Tmax),
                                                      mean.rangePrec  =mean(range.Prec),
                                                      sd.PlantHeight =sd(PlantHeight), 
                                                      sd.SeedMass    =sd(SeedMass),
                                                      sd.SeedLength  =sd(SeedLength),
                                                      sd.rangeTmax  =sd(range.Tmax),
                                                      sd.rangePrec  =sd(range.Prec))

ggplot(sitesCWMs, aes(y=mean.PlantHeight, x=mean.SeedMass))+
  geom_errorbar(aes(ymin=mean.PlantHeight-sd.PlantHeight, ymax=mean.PlantHeight + sd.PlantHeight))+
  geom_errorbarh(aes(xmin=mean.SeedMass-sd.SeedMass, xmax=mean.SeedMass+sd.SeedMass))+
  geom_point(aes(col=mean.SeedLength))+
  theme_classic()

ggplot(sitesCWMs, aes(y=mean.rangeTmax, x=mean.rangePrec))+
  geom_errorbar(aes(ymin=mean.rangeTmax-sd.rangeTmax, ymax=mean.rangeTmax + sd.rangeTmax))+
  geom_errorbarh(aes(xmin=mean.rangePrec-sd.rangePrec, xmax=mean.rangePrec+sd.rangePrec))+
  geom_point(aes(col=mean.SeedLength))+
  theme_classic()


save(sitesCWMs, file="results/biotic_predictors_CWMs_TRY.Rdata")
