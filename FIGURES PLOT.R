setwd("~/PIERRE/TRANSIENT DYNAMICS")
library(tidyverse)
library(gridExtra)
library(viridis)
library(cowplot)
load('results/CommClimStat2_smoothed_BBS.Rdata')
CommClimStat2_BBS<- CommClimStat2

load("results/statistics_sites_BBS.Rdata")
statistics_sites_BBS <- statistics_sites
colnames(statistics_sites_BBS)[1:3] <- c('site', 'Longitude', 'Latitude') 
statistics_sites_BBS$dataset <- 'BBS'

load('results/CommClimStat2_smoothed.Rdata')
CommClimStat2_Neotoma<- CommClimStat2

load("results/statistics_sites_neotoma.Rdata")
statistics_sites_Neotoma <- statistics_sites
colnames(statistics_sites_Neotoma)[c(7,9,11,13)] <- c("crossings.max", "deviation.absolute", "deviation.change.1Ky", "num.points.not.na")
statistics_sites_Neotoma$dataset <- 'Neotoma'


statMerged<- rbind(statistics_sites_BBS[, c(1,2,3,6,8,10,12,17)], as.data.frame(statistics_sites_Neotoma[, c(1,4,5,7,9,11,13,33)]))

load("results/predictors_BBS.Rdata")
predictors_BBS<- predictors
load("results/predictors.Rdata")
load("results/biotic_predictors_CWMs_TRY.Rdata")

predMerged<-rbind(as.data.frame(predictors[,c(1,2,3,6,10)]), predictors_BBS[,c(1,3,4,5,7)])



#### merge datas
modData <- left_join(statMerged, predMerged)


##### FIGURE 2 #####
sitesPlots <- list()
# "Rockyhock Bay", 
sites_examples_neotoma<- c("Big Pond",  "Otisville", "Alexis Lake", "West Hawk Lake")
for( i in 1:length(sites_examples_neotoma)){
  x<- CommClimStat2[CommClimStat2$SiteName==sites_examples_neotoma[i],]
  sitesPlots[[i]]<-
    ggplot(x, aes(y=CommTmax.ave, x=ClimTmax.ave))+
    geom_errorbarh(aes(xmin=ClimTmax.ave-ClimTmax.sd, xmax=ClimTmax.ave+ClimTmax.sd), col="grey")+
    geom_errorbar(aes(ymin=CommTmax.ave-CommTmax.sd, ymax=CommTmax.ave+CommTmax.sd), col="grey")+
    geom_path(col="grey")+
    geom_point(shape=21, fill="white", col="grey", size=3)+
    geom_path(aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess), col='red')+
    geom_abline(slope=1, intercept)+
    theme_classic()+
    stat_smooth(aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess), method="lm", fill=NA, col="black", linetype=3, fullrange = TRUE)+
    geom_text(data=x[c(1,nrow(x)),], aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess,label=paste("-",time, sep="")),nudge_x = 0, nudge_y = -1, col="black", fontface = "bold")+
    labs(title=as.character(x$SiteName[1]),
         subtitle=paste("abs.Dev =", round(x$deviation.absolute.Tmax[1],2), " | ",
                        "dev.change =", round(x$deviation.change.Tmax[1],2), " | ",
                        "n =", round(x$crossings.max.Tmax[1],2)),
         y="Community inferred maximum temperature (°C)",
         x="Maximum temperature (°C)")+
    coord_fixed(ratio=1)+
    lims(
      x=c(min(c(x$CommTmax.ave-x$CommTmax.sd, x$ClimTmax.ave-x$ClimTmax.sd)), max(c(x$CommTmax.ave+x$CommTmax.sd, x$ClimTmax.ave+x$ClimTmax.sd))),
      y=c(min(c(x$CommTmax.ave-x$CommTmax.sd, x$ClimTmax.ave-x$ClimTmax.sd)), max(c(x$CommTmax.ave+x$CommTmax.sd, x$ClimTmax.ave+x$ClimTmax.sd))))
}

# "840_25_10",
# "840_14_71"
# "840_91_57"; "840_7_27" loop !
# 840_46_42
sites_examples_BBS<- c("840_14_182","840_72_30", "840_52_28", "840_88_22")
for( i in 1:length(sites_examples_BBS)){
  x<- CommClimStat2_BBS[CommClimStat2_BBS$site==sites_examples_BBS[i],]
  sitesPlots[[length(sites_examples_neotoma)+i]]<-
    ggplot(x, aes(y=CommTmax.ave, x=ClimTmax.ave))+
    geom_errorbarh(aes(xmin=ClimTmax.ave-ClimTmax.sd, xmax=ClimTmax.ave+ClimTmax.sd), col="grey")+
    geom_errorbar(aes(ymin=CommTmax.ave-CommTmax.sd, ymax=CommTmax.ave+CommTmax.sd), col="grey")+
    geom_path(col="grey")+
    geom_point(shape=21, fill="white", col="grey", size=3)+
    geom_path(aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess), col='red')+
    geom_abline(slope=1, intercept)+
    theme_classic()+
    stat_smooth(aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess), method="lm", fill=NA, col="black", linetype=3, fullrange = TRUE)+
    geom_text(data=x[c(1,nrow(x)),], aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess,label=paste(time, sep="")),nudge_x = 0, nudge_y = -1, col="black", fontface = "bold")+
    labs(title=as.character(x$site[1]),
         subtitle=paste("abs.Dev =", round(x$deviation.absolute[1],2), " | ",
                        "dev.change =", round(x$deviation.change[1],2), " | ",
                        "n =", round(x$crossings.max[1],2)),
         y="Community inferred maximum temperature (°C)",
         x="Maximum temperature (°C)")+
    coord_fixed(ratio=1)+
    lims(
      x=c(min(c(x$CommTmax.ave-x$CommTmax.sd, x$ClimTmax.ave-x$ClimTmax.sd)), max(c(x$CommTmax.ave+x$CommTmax.sd, x$ClimTmax.ave+x$ClimTmax.sd))),
      y=c(min(c(x$CommTmax.ave-x$CommTmax.sd, x$ClimTmax.ave-x$ClimTmax.sd)), max(c(x$CommTmax.ave+x$CommTmax.sd, x$ClimTmax.ave+x$ClimTmax.sd))))
}

# ggsave("results/figures/Figure.2.pdf", arrangeGrob(grobs = sitesPlots, ncol=length(sitesPlots)/2, nrow=2), width =14, height = 8)
# ggsave("results/figures/Figure.2.png", arrangeGrob(grobs = sitesPlots, ncol=length(sitesPlots)/2, nrow=2), width =14, height = 8)




pdf("results/figures/Figure.2.pdf", width =14, height= 8)
plot_grid(plotlist = sitesPlots,  ncol=4,  axis='bottom', labels="auto")
dev.off()

png("results/figures/Figure.2.png", width =14, height= 8, units="in", res=200)
plot_grid(plotlist = sitesPlots,  ncol=4,  axis='bottom', labels="auto")
dev.off()



##### FIGURE 3 #####


## create background maps
library(maps)
world <-map_data("world") 
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  theme(legend.position="bottom")+
  coord_map("albers", parameters = c(-100, -100),  ylim=c(22,63), xlim=c(-130, -50))

##Deviation absolute
Map_devAbs_Plants <-
worldmap + 
  geom_point(data=statistics_sites_Neotoma[statistics_sites_Neotoma$deviation.absolute<20,], size=3, alpha=0.6,aes(Longitude, Latitude , fill=deviation.absolute), shape=21)+
  # scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlGnBu")+
  labs(x = "Longitude", y="Latitude")
  # scale_fill_viridis()

# annotate("text", label="Maximal temperatur

hist_devAbs_Plants <- 
  ggplot(statistics_sites_Neotoma[statistics_sites_Neotoma$deviation.absolute<20,], aes(x=deviation.absolute, fill=..x..))+
  geom_histogram(col="black")+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlGnBu")+
  labs(x =  expression(paste("Absolute deviation ",bar(Delta))), y="Plants - Number of sites")
  # scale_fill_viridis()

Map_devAbs_Birds <-
  worldmap + 
  geom_point(data=statistics_sites_BBS[statistics_sites_BBS$deviation.absolute<10,], size=3, alpha=0.6, aes(Longitude, Latitude , fill=deviation.absolute), shape=21)+
  # scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlGnBu")+
  labs(x = "Longitude", y="Latitude")
# annotate("text", label="Maximal temperatur

hist_devAbs_Birds <- 
  ggplot(statistics_sites_BBS[statistics_sites_BBS$deviation.absolute<10,], aes(x=deviation.absolute, fill=..x..))+
  geom_histogram(col='black')+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlGnBu")+
  labs(x =  expression(paste("Absolute deviation ", bar(Lambda))), y="Birds - Number of sites")

##Deviation change
Map_devChange_Plants <-
  worldmap + 
  geom_point(data=statistics_sites_Neotoma, aes(Longitude, Latitude , fill=deviation.change.1Ky), size=3, alpha=0.6,shape=21)+
  # scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "BrBG", limits = c(-1,1)*max(abs(statistics_sites_Neotoma$deviation.change.1Ky)))+
  labs(x = "Longitude", y="Latitude")
# annotate("text", label="Maximal temperatur

hist_devChange_Plants <- 
  ggplot(statistics_sites_Neotoma[statistics_sites_Neotoma$deviation.change.1Ky>-4,], aes(x=deviation.change.1Ky, fill=..x..))+
  geom_histogram(col='black')+
  theme_classic()+
  scale_fill_distiller(palette = "BrBG", limits = c(-1,1)*max(abs(statistics_sites_Neotoma$deviation.change.1Ky)))+
  theme(legend.position="none")+
  geom_vline(xintercept=0)+
  labs(x =  expression(paste("Deviation change ", italic(d), Lambda)), y="Plants - Number of sites")


Map_devChange_Birds <-
  worldmap + 
  geom_point(data=statistics_sites_BBS, aes(Longitude, Latitude , fill=deviation.change.1Ky), size=3, alpha=0.6, shape=21)+
  # scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "BrBG", limits = c(-1,1)*max(abs(statistics_sites_BBS$deviation.change.1Ky)))+
  labs(x = "Longitude", y="Latitude")
# annotate("text", label="Maximal temperatur

hist_devChange_Birds <- 
  ggplot(statistics_sites_BBS, aes(x=deviation.change.1Ky, fill=..x..))+
  geom_histogram(col='black')+
  theme_classic()+
  # scale_fill_gradient2()+
  scale_fill_distiller(palette = "BrBG", limits = c(-1,1)*max(abs(statistics_sites_BBS$deviation.change.1Ky)))+
  theme(legend.position="none")+
  geom_vline(xintercept=0)+
  labs(x =  expression(paste("Deviation change ", italic(d), Lambda)), y="Birds - Number of sites")


# grid.arrange(Map_devChange_Birds, hist_devChange_Plants, Map_devChange_Plants, hist_devChange_Birds)


statistics_sites_Neotoma[statistics_sites_Neotoma$deviation.absolute,20,]
##state number
Map_N_Plants <-
  worldmap + 
  geom_point(data=na.omit(statistics_sites_Neotoma), aes(Longitude, Latitude , fill=as.integer(crossings.max)),size=3, alpha=0.6, shape=21)+
  # scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlOrRd")+
  labs(x = "Longitude", y="Latitude")
# annotate("text", label="Maximal temperatur

hist_N_Plants <- 
  ggplot(statistics_sites_Neotoma, aes(x=as.integer(crossings.max), fill=..x..))+
  geom_histogram(binwidth=1, col='black')+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlOrRd")+
  labs(x =  expression(paste("Maximum state number ", italic(n))), y="Plants - Number of sites")


Map_N_Birds <-
  worldmap + 
  geom_point(data=statistics_sites_BBS, aes(Longitude, Latitude , fill=as.integer(crossings.max)),size=3, alpha=0.6, shape=21)+
  # scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlOrRd")+
  labs(x = "Longitude", y="Latitude")
# annotate("text", label="Maximal temperatur

hist_N_Birds <- 
  ggplot(statistics_sites_BBS, aes(x=as.integer(crossings.max), fill=..x..))+
  geom_histogram(binwidth=1, col='black')+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_distiller(palette = "YlOrRd")+
  labs(x =  expression(paste("Maximum state number ", italic(n))), y="Birds - Number of sites")


#### ALL #####
library(cowplot)
png("results/figures/Figure.3.png",width =10, height=16.5, units="in", res=200)
plot_grid(hist_devAbs_Plants, Map_devAbs_Plants, hist_devAbs_Birds,   Map_devAbs_Birds, 
          hist_devChange_Plants, Map_devChange_Plants, hist_devChange_Birds,  Map_devChange_Birds, 
          hist_N_Plants,Map_N_Plants, hist_N_Birds, Map_N_Birds,  ncol=2,  axis='bottom',  
          rel_widths =  rep(c( 1/3, 2/3), 3), labels=c("a", "", "b", "", "c", "", "d", "","e", "", "f"), hjust=-0.1)

dev.off()


pdf("results/figures/Figure.3.pdf", width =10, height=16.5)
plot_grid(hist_devAbs_Plants, Map_devAbs_Plants, hist_devAbs_Birds,   Map_devAbs_Birds, 
          hist_devChange_Plants, Map_devChange_Plants, hist_devChange_Birds,  Map_devChange_Birds, 
          hist_N_Plants,Map_N_Plants, hist_N_Birds, Map_N_Birds,  ncol=2,  axis='bottom',  
          rel_widths =  rep(c( 1/3, 2/3), 3), labels=c("a", "", "b", "", "c", "", "d", "","e", "", "f"), hjust=-0.1)

dev.off()



##### FIGURE 4 #####
modPred_res_BBS <- read.csv("results/modPred_results_BBS.csv")
colnames(modPred_res_BBS)<- c("coefficient", "SE", "z_value", "P_value","predictor","expl.var")
modPred_res_BBS$predictor <- factor(modPred_res_BBS$predictor, levels=c("Baseline.Climate","Baseline.Climate2", "Topography", "Human.Influence", "Body.mass", "Migration"))
modPred_res_BBS$expl.var<- rep(c("max.state.number", "absolute.deviation","deviation.change"), each=6)

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
  scale_color_manual(values=c(rgb(161,218,180,maxColorValue = 255), rgb(37,52,148,maxColorValue = 255),rgb(65,182,196,maxColorValue = 255)))+
  # scale_color_manual(values=c('grey80','grey60', 'grey30', 'black'))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())+
  scale_x_discrete(labels=c("Baseline Climate",expression(paste("Baseline Climate"^"2")), "Topography", "Human Influence", "Body mass", "Migration"))+
  ylab("Birds - Model coefficients")

modPred_res_Neotoma <- read.csv("results/modPred_results_Neotoma.csv")
modPred_res_Neotoma$predictor <- factor(modPred_res_Neotoma$predictor, levels=c("Baseline.Climate","Baseline.Climate2", "Topography","Human.Density","Plant.height", "Seed.mass"))
levels(modPred_res_Neotoma$predictor)
modPred_res_Neotoma$sig <- "ns"
modPred_res_Neotoma$sig[modPred_res_Neotoma$P_value < 0.1] <- "P<0.1"
modPred_res_Neotoma$sig[modPred_res_Neotoma$P_value < 0.05] <- "* P<0.05"
modPred_res_Neotoma$sig[modPred_res_Neotoma$P_value < 0.01] <- "** P<0.01"


plot_modPred_res_Neotoma<-
  ggplot(modPred_res_Neotoma, aes(x=predictor, y=coefficient, col=sig))+
  geom_pointrange(aes(ymin=coefficient-(1.96*SE), ymax=coefficient+(1.96*SE)))+
  geom_hline(yintercept = 0)+
  facet_wrap(~expl.var,nrow=1, scales='free_y', labeller=stat_labeller)+
  # scale_color_brewer(palette = "YlGnBu")+
  scale_color_manual(values=c(rgb(161,218,180,maxColorValue = 255), rgb(37,52,148,maxColorValue = 255),rgb(65,182,196,maxColorValue = 255)))+
  # scale_color_manual(values=c('grey80','grey60', 'grey30', 'black'))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank())+
  scale_x_discrete(labels=c("Baseline Climate",expression(paste("Baseline Climate"^"2")), "Topography", "Human Density","Plant height", "Seed mass"))+
  ylab("Plants - Model coefficients")


library(cowplot)
png("results/figures/Figure.4.png", width =7, height=8, units="in", res=200)
print(plot_grid(plot_modPred_res_Neotoma, plot_modPred_res_BBS, ncol=1, labels = "auto"))
dev.off()
pdf("results/figures/Figure.4.pdf", width =7, height=8)
print(plot_grid(plot_modPred_res_Neotoma, plot_modPred_res_BBS, ncol=1, labels = "auto"))
dev.off()          


#### supplementaru Figure S1 ####
ntime_neotoma <- CommClimStat2_Neotoma %>% group_by(site, Latitude, Longitude) %>% summarise(timespan=max(-time)-min(-time))

map_neo<-
  worldmap+
  geom_point(data=ntime_neotoma , aes(x=Longitude, y=Latitude,fill=timespan), shape=21, size=3)+
  scale_fill_viridis()+
  theme(legend.position=c("right"))

nsite_neotoma <- CommClimStat2_Neotoma %>% group_by(time) %>% summarise(n.site=length(unique(site)))
time_neo<-
  ggplot(nsite_neotoma, aes(y=n.site, x=-time, fill=n.site))+
  geom_col(col='black')+
  theme_classic()+
  scale_fill_viridis(option='magma')

ntime_bbs <- CommClimStat2_BBS %>% group_by(site, Lati, Longi) %>% summarise(timespan=max(time)-min(time))

map_bbs<- 
  worldmap+
  geom_point(data=ntime_bbs , aes(x=Longi, y=Lati, fill=timespan), shape=21, size=3)+
  scale_fill_viridis()+
  theme(legend.position=c("right"))

nsite_bbs <- CommClimStat2_BBS %>% group_by(time) %>% summarise(n.site=length(unique(site)))

time_bbs<-
  ggplot(nsite_bbs, aes(y=n.site, x=time, fill=n.site))+
  geom_col(col='black')+
  theme_classic()+
  scale_fill_viridis(option='magma')


png("results/figures/supp_FigS1.png", width =12, height=6, units="in", res=200)
print(plot_grid(map_neo, map_bbs, time_neo, time_bbs, ncol=2,  labels = "auto"))
dev.off()

pdf("results/figures/supp_FigS1.pdf", width =12, height=6)
print(plot_grid(map_neo, map_bbs, time_neo, time_bbs, ncol=2, labels = "auto"))
dev.off()    

#### supplementaru Figure S2 Contemporary and paleo plant thermal niches ####
niches_bien_summary <- read.csv('results/niches_summary_bien_worldclim.csv')
niches_paleo_summary <- read.csv('results/niches_all_summary.csv')
niches_bien_summary <- niches_bien_summary[niches_bien_summary$alltaxa %in% niches_paleo_summary$taxon,]
niches_paleo_summary <- niches_paleo_summary[niches_paleo_summary$taxon %in% niches_bien_summary$alltaxa,]

cor.test(niches_bien_summary$mean.Monthly.maxTemp, niches_paleo_summary$mean.Monthly.maxTemp)


#### supplementary figure occ VS abundance ####
plot(BirdsClim_means_occ$CommTmax.ave, BirdsClim_means_weighted$CommTmax.ave)

Birds_Temp<-data.frame(BirdsClim_means_occ[1:4])
colnames(Birds_Temp)[3:4]<- paste(colnames(Birds_Temp)[3:4], "occurence", sep="_")
Birds_Temp$CommTmax.ave_abundance<- BirdsClim_means_weighted$CommTmax.ave
Birds_Temp$CommTmax.var_abundance<- BirdsClim_means_weighted$CommTmax.var

ggplot(Birds_Temp, aes(CommTmax.ave_abundance, CommTmax.ave_occurence))+
  geom_point(shape=21, alpha=0.3)+
  # geom_bin2d(bins=120)+
  coord_fixed(ratio=1)+
  geom_abline(intercept=0, slope=1, col="red")+
  theme_classic()+
  labs(x="Community Inferred Temperature based on occurences", y="Abundance weighted Community Inferred Temperature" )

ggplot(Birds_Temp, aes(CommTmax.ave_abundance))+
  geom_histogram()+
  # coord_fixed(ratio=1)+
  theme_classic()

ggplot(Birds_Temp, aes(CommTmax.ave_occurence))+
  geom_histogram()+
  # coord_fixed(ratio=1)+
  theme_classic()

cor(Birds_Temp$CommTmax.ave_abundance, Birds_Temp$CommTmax.ave_occurence)
cor.test(Birds_Temp$CommTmax.ave_abundance, Birds_Temp$CommTmax.ave_occurence)

#### supplementaru Figure SX Modern niches and Neotoma sites ####

# World map, using geom_path instead of geom_polygon
world <- map_data("world")
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey") +
  coord_fixed(ratio=1, ylim=c(-55,80), xlim=c(-170, -30))


occurrences_bienWithWorldclim <- readRDS("results/occurrences_bienWithWorldclim.Rdata")
x <- droplevels(occurrences_bienWithWorldclim %>% filter(alltaxa=="Abies"))

load("results/statistics_sites_neotoma.Rdata")

pdf("results/modern_niche_maps_with_SITES.pdf")
d_ply(occurrences_bienWithWorldclim, .(alltaxa), function(x){
  
  print(
    worldmap+
      geom_point(data=x, aes(longitude, latitude,fill=Tmax), shape=21)+
      geom_point(data=statistics_sites_Neotoma, aes(Longitude, Latitude), col='red', shape=3, size=1, alpha=0.5)+
      scale_fill_viridis("max. temperature (°C)")+
      theme_classic()+
      theme(legend.justification=c(0,0), legend.position = c(0,0))+
      labs(title=as.character(x$alltaxa[1]), x="Longitude", y="Latitude")
  )
}, .progress="text")
dev.off()


#### supplementaru Figure SX predictors ####
##BBS
#### distributions of predictors ####
baseTmax <-
  ggplot(predictors_BBS, aes(x=base.Tmax, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis()+
  theme(legend.position="none")

topo <-
  ggplot(predictors_BBS, aes(x=topo25, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis(trans="sqrt")+
  theme(legend.position="none")+
  scale_x_continuous(trans="sqrt")

hii <-
  ggplot(predictors_BBS, aes(x=HII, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis()+
  theme(legend.position="none")

bodymass <-
  ggplot(predictors_BBS, aes(x=BodyMass, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis()+
  theme(legend.position="none")

migration <-
  ggplot(predictors_BBS, aes(x=Migration, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis()+
  theme(legend.position="none")


#### maps of predictors ####

## create bacground map
library(maps)
world <-map_data("world") 
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  theme(legend.position="bottom")+
  coord_map("albers", parameters = c(-100, -100),  ylim=c(22,63), xlim=c(-130, -50))


baseTmax_map <-
  worldmap + 
  geom_point(data=predictors_BBS, aes(y=Latitude, x=Longitude, fill=base.Tmax),size=4, shape=21)+
  scale_fill_viridis()+
  theme(legend.position="none")+
  annotate("text", label="baseline max. temp.", x=-60,y=30, size=4)

topo_map <-
  worldmap + 
  geom_point(data=predictors_BBS, aes(y=Latitude, x=Longitude, fill=topo25),size=4, shape=21)+
  scale_fill_viridis(trans='sqrt')+
  theme(legend.position="none")+
  annotate("text", label="sqrt. Topography", x=-70,y=30, size=4)

hii_map <-
  worldmap + 
  geom_point(data=predictors_BBS, aes(y=Latitude, x=Longitude, fill=HII),size=4, shape=21)+
  scale_fill_viridis()+
  theme(legend.position="none")+
  annotate("text", label="Human Influence Index", x=-70,y=30, size=4)

bodymass_map <-
  worldmap + 
  geom_point(data=predictors_BBS, aes(y=Latitude, x=Longitude, fill=BodyMass),size=4, shape=21)+
  scale_fill_viridis()+
  theme(legend.position="none")+
  annotate("text", label="mean body mass", x=-70,y=30, size=4)

migration_map <-
  worldmap + 
  geom_point(data=predictors_BBS, aes(y=Latitude, x=Longitude, fill=Migration),size=4, shape=21)+
  scale_fill_viridis()+
  theme(legend.position="none")+
  annotate("text", label="% migrants", x=-70,y=30, size=4)

library(cowplot)

png("results/figures/supp_predictors_BBS.png", width =10, height=16, units="in", res=200)

print(plot_grid(baseTmax, baseTmax_map, 
                topo, topo_map,
                hii, hii_map,
                bodymass, bodymass_map,
                migration, migration_map, ncol=2,  axis='bottom',  
                rel_widths =  rep(c( 1/3, 2/3), 3), labels="auto", hjust=-0.1))

dev.off()

pdf("results/figures/supp_FigS1.pdf", width =12, height=6)
print(plot_grid(map_neo, map_bbs, time_neo, time_bbs, ncol=2, labels = "auto"))
dev.off()    

##Neotoma
#### distributions of predictors_Neotoma ####
load("results/predictors_Neotoma.Rdata")
load("results/biotic_predictors_Neotoma_CWMs_TRY.Rdata")
predictors_Neotoma<- left_join(predictors, sitesCWMs)

baseTmax <-
  ggplot(predictors_Neotoma, aes(x=base.Tmax, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis()+
  theme(legend.position="none")

topo <-
  ggplot(predictors_Neotoma, aes(x=topo25, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis(trans="sqrt")+
  theme(legend.position="none")+
  scale_x_continuous()

dens <-
  ggplot(predictors_Neotoma, aes(x=humanDensity, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis(trans="sqrt")+
  theme(legend.position="none")+
  scale_x_continuous(trans="sqrt")

seedMass <-
  ggplot(sitesCWMs, aes(x=mean.SeedMass, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis()+
  theme(legend.position="none")

plantHeight <-
  ggplot(sitesCWMs, aes(x=mean.PlantHeight, fill=..x..))+
  geom_histogram()+
  scale_fill_viridis()+
  theme(legend.position="none")


#### maps of predictors_Neotoma ####

## create bacground map
library(maps)
world <-map_data("world") 
worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  theme(legend.position="bottom")+
  coord_map("albers", parameters = c(-100, -100),  ylim=c(22,63), xlim=c(-130, -50))


baseTmax_map <-
  worldmap + 
  geom_point(data=predictors_Neotoma, aes(y=Latitude, x=Longitude, fill=base.Tmax),size=4, shape=21)+
  scale_fill_viridis()+
  theme_classic()+
  theme(legend.position="none")+
  annotate("text", label="baseline max. temp.", x=-60,y=30, size=4)

topo_map <-
  worldmap + 
  geom_point(data=predictors_Neotoma, aes(y=Latitude, x=Longitude, fill=topo25),size=4, shape=21)+
  scale_fill_viridis(trans='sqrt')+
  theme_classic()+
  theme(legend.position="none")+
  annotate("text", label="sqrt. Topography", x=-70,y=30, size=4)

dens_map <-
  worldmap + 
  geom_point(data=predictors_Neotoma, aes(y=Latitude, x=Longitude, fill=humanDensity),size=4, shape=21)+
  scale_fill_viridis(trans="sqrt")+
  theme_classic()+
  theme(legend.position="none")+
  annotate("text", label="Human Density", x=-70,y=30, size=4)

seedMass_map <-
  worldmap + 
  geom_point(data=predictors_Neotoma, aes(y=Latitude, x=Longitude, fill=mean.SeedMass),size=4, shape=21)+
  scale_fill_viridis()+
  theme(legend.position="none")+
  annotate("text", label="mean body mass", x=-70,y=30, size=4)

plantHeight_map <-
  worldmap + 
  geom_point(data=predictors_Neotoma, aes(y=Latitude, x=Longitude, fill=mean.PlantHeight),size=4, shape=21)+
  scale_fill_viridis()+

  theme(legend.position="none")+
  annotate("text", label="% migrants", x=-70,y=30, size=4)

library(cowplot)

png("results/figures/supp_predictors_Neotoma.png", width =10, height=16, units="in", res=200)

print(plot_grid(baseTmax, baseTmax_map, 
                topo, topo_map,
                dens, dens_map,
                seedMass, seedMass_map,
                plantHeight, plantHeight_map, ncol=2,  axis='bottom',  
                rel_widths =  rep(c( 1/3, 2/3), 3), labels="auto", hjust=-0.1))

dev.off()

pdf("results/figures/supp_FigS1.pdf", width =12, height=6)
print(plot_grid(map_neo, map_bbs, time_neo, time_bbs, ncol=2, labels = "auto"))
dev.off()    