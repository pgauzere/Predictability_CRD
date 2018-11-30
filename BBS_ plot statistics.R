
###############################################################################################
#################### Plot Phase diagrams, statistics maps and distributions  ##################
#################################################### Gauzere & Blonder ########################

setwd("~/PIERRE/TRANSIENT DYNAMICS")
library(tidyverse)
library(viridis)
library(gridExtra)
library(cowplot)

load("results/statistics_sites_BBS.Rdata")
load("results/CommClimStat2_smoothed_BBS.Rdata")

CommClimStat2$CommTmax.sd <- sqrt(CommClimStat2$CommTmax.var)

#### plot community response diagrams ####

library(plyr)
# x<- CommClimStat2[CommClimStat2$site=="840_25_24",]


#add t0-first time step and tmax-last time step to commclimstat
period_study <-CommClimStat2 %>% group_by(site) %>% dplyr::summarize(t0 = min(time), tmax = max(time))
hist(period_study$t0, breaks=50)
CommClimStat2 <- left_join(CommClimStat2, period_study) %>% arrange(t0)

sites <- unique(CommClimStat2$site)

# x<- CommClimStat2_BBS[CommClimStat2_BBS$site=="840_66_17",]

pdf("results/figures/phaseDiagrams_BBS_means_smoothed_stats2_lmSmoothed.pdf", width=12, height=8)
d_ply(na.omit(CommClimStat2_BBS), .(site), function(x){
  print(x$site[1])

crd<-
  ggplot(na.omit(x), aes(y=CommTmax.ave, x=ClimTmax.ave))+
    geom_errorbarh(aes(xmin=ClimTmax.ave-ClimTmax.sd, xmax=ClimTmax.ave+ClimTmax.sd), col="grey")+
    geom_errorbar(aes(ymin=CommTmax.ave-CommTmax.sd, ymax=CommTmax.ave+CommTmax.sd), col="grey")+
    geom_path(col="grey")+
    geom_point(shape=21, fill="white", col="grey", size=3)+
    geom_path(aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess), col='red')+
    geom_abline(slope=1, intercept)+
    theme_classic() + 
    stat_smooth(aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess), method="lm", fill=NA, col="black", linetype=3, fullrange = TRUE)+
    geom_text(data=x[c(1,nrow(x)),], aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess,label=paste(time, sep="")),nudge_x = 0, nudge_y = -0.5, col="black", fontface = "bold")+
    labs(
         subtitle=paste("absolute deviation =", round(x$deviation.absolute[1],2), " | ",
                        "deviation change =", round(x$deviation.change[1],2), " | ",
                        "maximum state number =", round(x$crossings.max[1],2)),
         y="Community inferred maximum temperature (°C)", 
         x="Observed maximum temperature (°C)")+
    coord_fixed(ratio=1)+
    lims(
      x=c(min(c(x$CommTmax.ave-x$CommTmax.sd, x$ClimTmax.ave-x$ClimTmax.sd)), max(c(x$CommTmax.ave+x$CommTmax.sd, x$ClimTmax.ave+x$ClimTmax.sd))),
      y=c(min(c(x$CommTmax.ave-x$CommTmax.sd, x$ClimTmax.ave-x$ClimTmax.sd)), max(c(x$CommTmax.ave+x$CommTmax.sd, x$ClimTmax.ave+x$ClimTmax.sd))))

timeComm<-
  ggplot(na.omit(x), aes(y=CommTmax.ave, x=time))+
  # geom_errorbarh(aes(xmin=ClimTmax.ave-ClimTmax.sd, xmax=ClimTmax.ave+ClimTmax.sd), col="grey")+
  geom_path(col="grey")+
  geom_point(shape=21, fill="white", col="grey", size=3)+
  geom_ribbon(aes(ymin=CommTmax.ave_loess-CommTmax.ave_loess_se, ymax=CommTmax.ave_loess+CommTmax.ave_loess_se), fill="red", alpha=0.1)+
  geom_path(aes(y=CommTmax.ave_loess, x=time), col='red')+
  geom_abline(slope=1, intercept)+
  theme_classic() +
  labs(y="Community inferred maximum temperature (°C)", 
       x="")

timeClim<-
  ggplot(na.omit(x), aes(y=ClimTmax.ave, x=time))+
    # geom_errorbarh(aes(xmin=ClimTmax.ave-ClimTmax.sd, xmax=ClimTmax.ave+ClimTmax.sd), col="grey")+
    # geom_errorbar(aes(ymin=ClimTmax.ave-ClimTmax.sd, ymax=ClimTmax.ave+ClimTmax.sd), col="grey")+
    geom_path(col="grey")+
    geom_point(shape=21, fill="white", col="grey", size=3)+
    geom_ribbon(aes(ymin=ClimTmax.ave_loess-ClimTmax.ave_loess_se, ymax=ClimTmax.ave_loess+ClimTmax.ave_loess_se), fill="red", alpha=0.1)+
    geom_path(aes(y=ClimTmax.ave_loess, x=time), col='red')+
    geom_abline(slope=1, intercept)+
    theme_classic() +
  labs(y="Observed maximum temperature (°C)", 
       x="Time (years)")  

left_col <- 
  plot_grid(timeComm, timeClim, ncol=1,  align="v")

plot<- plot_grid(left_col, crd, ncol=2,  align="h",rel_widths =c( 1, 2))
title <- ggdraw() + draw_label(paste("site", x$site[1], sep=" : "), fontface='bold')

print(plot_grid(title, plot, ncol=1, rel_heights=c(0.1, 1)))

}, .progress="text")
dev.off() 

#### stacked CRD ####


commclim.Tmax<-
  ggplot(CommClimStat2, aes(y=CommTmax.ave, x=ClimTmax.ave))+
  geom_errorbar(aes(ymin=CommTmax.ave-CommTmax.sd, ymax=CommTmax.ave+CommTmax.sd), alpha=0.01)+
  geom_path(aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess, group=site, col=t0))+
  theme_classic()+
  stat_smooth(method='gam', formula=y~s(x,k=3))+
  scale_color_viridis(guide=F)+
  # scale_color_discrete(guide=F)+
  # facet_wrap(~period2)+
  labs(y="Community inferred maximum temperature (°C)", 
       x="Maximum temperature (°C)")+
  # coord_fixed(ratio=1)+
  geom_abline(slope=1, intercept=0)


climTmax<-
  ggplot(CommClimStat2, aes(y=ClimTmax.ave, x=time))+
  # geom_point()+
  # geom_line(aes(group=site, col=t0))+
  theme_classic()+
  scale_color_viridis()+
  theme(legend.position="bottom")+
  # geom_vline(xintercept = c(-16000,-8000))+
  stat_smooth(col='black')

grid.arrange(commclim.Tmax,climTmax, ncol=1)



##pairs correlation of statistics


#### plot distribution of statistics ####

hist_devAbs.t <- 
  ggplot(statistics_sites, aes(x=deviation.absolute ))+
  geom_histogram(fill="red")+
  theme_classic()+
  geom_vline(xintercept = 0)

hist_devChange.t <- 
  ggplot(statistics_sites, aes(x=deviation.change.1Ky ))+
  geom_histogram(fill="red")+
  theme_classic()+
  geom_vline(xintercept = 0)

hist_crossMeans.t <- 
  ggplot(statistics_sites, aes(x=crossings.mean ))+
  geom_histogram(breaks=c(1,1.5,2,2.5,3,4),fill="red")+
  theme_classic()+
  geom_vline(xintercept = 1)+
  geom_vline(xintercept = c(1.5,2,2.5,3), linetype=3)

hist_crossMax.t <- 
  ggplot(statistics_sites, aes(x=crossings.max ))+
  geom_histogram(breaks=c(0.5,1.5,2.5,3.5,4.5),fill="red")+
  theme_classic()

hist_slope.t <-
  ggplot(statistics_sites, aes(x=asinh(lm.slope *10)))+
  geom_histogram(fill="red", bins=100)+
  theme_classic()+
  geom_vline(xintercept = asinh(10))+
  geom_vline(xintercept = 0, linetype=3)

hist_lmRho.t <-
  ggplot(statistics_sites, aes(x=lm.rho ))+
  geom_histogram(fill="red", boundary=1)+
  theme_classic()+
  geom_vline(xintercept = 1)+
  geom_vline(xintercept = 0, linetype=3)

hist_Npoints <- 
  ggplot(statistics_sites, aes(x=num.points.not.na ))+
  geom_histogram()+
  theme_classic()

pdf("results/distribution_summary_statistics_BBS.pdf", width=7)
print(hist_devAbs.t)
print(hist_devChange.t)
print(hist_crossMax.t)
# grid.arrange(hist_lmRho.t , hist_lmRho.p, nrow=1)
dev.off()


#### map the statistics ####
library(maps)
library(ggplot2)
library(rje)

world <-map_data("world") 

worldmap <-
  ggplot() + 
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
  coord_fixed(ratio=1, ylim=c(25,52), xlim=c(-125, -65))+
  theme(legend.position="bottom")

devAbs.t <-
  worldmap + 
  geom_point(data=statistics_sites, aes(Longi, Lati , fill=deviation.absolute),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")
  # annotate("text", label="Maximal temperature", x=-70,y=30, size=5, col="red")

devChange.t <- 
  worldmap + 
  geom_point(data=statistics_sites, aes(Longi, Lati , fill=deviation.change.1Ky),size=3 , shape=21)+
  scale_fill_gradient2()+
  theme_classic()+
  theme(legend.position="bottom")
  # annotate("text", label="Maximal temperature", x=-65,y=30, size=5, col="red")

crossMeans.t <- 
  worldmap + 
  geom_point(data=statistics_sites, aes(Longi, Lati , fill=crossings.mean),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")
  # annotate("text", label="Maximal temperature", x=-65,y=30, size=5, col="red")

crossMax.t <- 
  worldmap + 
  geom_point(data=statistics_sites, aes(Longi, Lati , fill=crossings.max),size=3 , shape=21)+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position="bottom")
  # annotate("text", label="Maximal temperature", x=-65,y=30, size=5, col="red")

lmRho.t <-
  worldmap + 
  geom_point(data=statistics_sites, aes(Longi, Lati , fill=lm.rho),size=3 , shape=21)+
  scale_fill_gradient2()+
  theme_classic()+
  theme(legend.position="bottom")
  # annotate("text", label="Maximal temperature", x=-65,y=30, size=5, col="red")
 
# npoint <- 
#   worldmap + 
#   geom_point(data=statistics_sites, aes(Longi, Lati , fill=num.points.not.na.Tmax),size=3, alpha=0.7, shape=21)+
#   scale_fill_gradientn(colours=cubeHelix(100))+
#   theme_classic()+
#   theme(legend.position="bottom")

pdf("results/maps_summary_statistics_BBS.pdf", width=16)
print( devAbs.t )
print( devChange.t)
print( crossMax.t)
dev.off()


#### interpolation ####

library(raster)   #Load the Raster Library
library(maps)
library(sp)
library(plotKML)
library(automap)

us <- raster::getData('GADM', country= "USA", level=0) #see ccodes()
can <- raster::getData('GADM', country= "CAN", level=0) #see ccodes()
e <- extent(-115, -55, 25, 60)
northAmerica <- crop(cbind(us, can), e)

northAmerica_grid <- raster(vect2rast(northAmerica,cell.size=1))

spatial_stats <- statistics_sites
coordinates(spatial_stats) <- ~ Longitude + Latitude
proj4string(spatial_stats) <- proj4string(northAmerica_grid)

kridge_stats <-autoKrige(formula = deviation.absolute.Tmax~1, input_data = spatial_stats)\, new_data = northAmerica_grid)
plot(kridge_stats)

k_stats<- raster(kridge_stats$krige_output)

# now use the mask function
rr <- mask(k_stats, northAmerica)

# plot, and overlay:
plot(rr)
plot(northAmerica,add=TRUE)

france.fort <- fortify(france, region='OBJECTID')

#### correlations between statistics ####

##pairs correlation of statistics
library("PerformanceAnalytics")
chart.Correlation(statistics_sites[, c(
 "crossings.mean", "crossings.max" , "deviation.base", 
"deviation.absolute" , "deviation.change.1Ky",
"num.points.not.na"  , "lm.intercept",   
"lm.slope" , "lm.r2"    , "lm.rho")], histogram=TRUE, pch=19)



#### general plot ####
library(viridis)

period_study <-CommClimStat2 %>% group_by(site) %>% dplyr::summarize(t0 = min(-time), tmax = max(-time))
hist(period_study$t0, breaks=50)
CommClimStat2 <- left_join(CommClimStat2, period_study)


# commclim<-
  ggplot(CommClimStat2, aes(y=CommTmax.ave, x=ClimTmax.ave))+
  geom_errorbar(aes(ymin=CommTmax.ave-CommTmax.sd, ymax=CommTmax.ave+CommTmax.sd), alpha=0.01)+
  geom_path(aes(y=CommTmax.ave_loess, x=ClimTmax.ave_loess, group=site, 
                col=deviation.absolute))+
  theme_classic()+
  stat_smooth(method='gam', formula=y~s(x, k=3), col='blue')+
  scale_color_viridis(guide=F)+
  # scale_color_gradient2()+
  # scale_color_discrete(guide=F)+
  # facet_wrap(~period2)+
  theme(legend.position="bottom")+
  labs(y="Community inferred maximum temperature (°C)", 
       x="Maximum temperature (°C)")+
  # coord_fixed(ratio=1)+
  geom_abline(slope=1, intercept=0)

clim<-
  ggplot(CommClimStat2 , aes(y=ClimTmax.ave, x=time))+
  # geom_point()+
  geom_line(aes(group=site), alpha=0.05)+
  theme_classic()+
  scale_color_viridis()+
  theme(legend.position="bottom")+
  # geom_vline(xintercept = c(-8000))+
  stat_smooth( col='black')
  facet_wrap(~period2, scales='free_x')

mod <- gamm(ClimTmax.ave ~ s(time), random=list(site=~1), data=CommClimStat2)  
gam.check(mod$gam)
summary(mod$gam)
plot(mod$gam)

  
grid.arrange(commclim.SeedMass,commclim.PlantHeight, clim)