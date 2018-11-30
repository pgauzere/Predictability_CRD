#### 

setwd("~/PIERRE/TRANSIENT DYNAMICS")

library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)
library(rje)


niches_bien_all  <- na.omit(readRDS("results/occurrences_bienWithWorldclim.Rdata")[,-c(2,3)])
unique(niches_bien_all$alltaxa)
niches_bien_all$alltaxa <- as.factor(gsub(" ", ".", niches_bien_all$alltaxa))
colnames(niches_bien_all)[1] <- "taxon"
niches_bien_all$time <- 2000
sort(unique(niches_bien_all$taxon))


niches_paleo_all <- read.csv('results/niches_paleo_all_raw.csv')[,c(4,10,11,3)]
colnames(niches_paleo_all)[c(1,2,3)] <- c("taxon", "Tmax", "Prec")
niches_paleo_all$time <- -niches_paleo_all$time
niches_paleo_all <- niches_paleo_all[niches_paleo_all$taxon %in% niches_bien_all$taxon,]
niches_paleo_all$taxon <- droplevels(niches_paleo_all$taxon)
sort(unique(niches_paleo_all$taxon))

niches_paleo_all_sampled <- na.omit(read.csv('results/niches_paleo_all_sampled.csv')[,-1])
colnames(niches_paleo_all_sampled)[c(1,2,3)] <- c("taxon", "Tmax", "Prec")
niches_paleo_all_sampled$time <- -niches_paleo_all_sampled$time
niches_paleo_all_sampled <- niches_paleo_all_sampled[niches_paleo_all_sampled$taxon %in% niches_bien_all$taxon,]
niches_paleo_all_sampled$taxon <- droplevels(niches_paleo_all_sampled$taxon)
sort(unique(niches_paleo_all_sampled$taxon))

niches_bien_all <- niches_bien_all[niches_bien_all$taxon %in% niches_paleo_all$taxon,]
niches_bien_all$taxon <- droplevels(niches_bien_all$taxon)
sort(unique(niches_bien_all$taxon))

niches_all <- na.omit(rbind(niches_bien_all[, c('taxon', 'Tmax', 'Prec', 'time')], niches_paleo_all))

niches_all_paleosampled <- na.omit(rbind(niches_bien_all[, c('taxon', 'Tmax', 'Prec', 'time')], niches_paleo_all_sampled))
save(niches_all_paleosampled, file="results/niches_all_paleosampled.Rdata")

pdf("results/niche_time_raw.pdf",width=7, height=4)
d_ply(na.omit(niches_all), .(taxon), function(x){
print(
  ggplot(x, aes(x=time, y=Tmax))+
  # stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
  geom_point(shape=21, alpha=0.2)+
  geom_hline(aes(yintercept = mean(x$Tmax)), col="red")+
  geom_hline(aes(yintercept = mean(x$Tmax)-sd(x$Tmax)), col="red", linetype=5)+
  geom_hline(aes(yintercept = mean(x$Tmax)+sd(x$Tmax)), col="red", linetype=5)+
  # geom_point(alpha=0.3, shape=21)+
  stat_smooth(col="blue")+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position = "none")+
  coord_cartesian(ylim=c(-20,35), xlim=c(-21000, 2000))+
  ggtitle(as.character(x$taxon[1]))
)},.progress="text")
dev.off()

pdf("results/niche_space_raw.pdf",width=8, height=8)
d_ply(na.omit(niches_all), .(taxon), function(x){
  
  # x <- niches_all[niches_all$taxon == "Alnus",]
  xpal <- x[x$time < 1,]
  xbien<- x[x$time > 1,]
  
  print(
    ggplot()+
      # stat_density2d(data=x, aes(x=Prec, y=Tmax), col="black")+
      # stat_density2d(data=xpal, aes(x=Prec, y=Tmax), col="red")+
      # stat_density2d(data=xbien, aes(x=Prec, y=Tmax), col="blue")+
      geom_point(data=xbien, aes(x=Prec, y=Tmax),alpha=0.05, col="blue", shape=21)+
      geom_point(data=xpal, aes(x=Prec, y=Tmax),alpha=0.05, col="red", shape=21)+
      geom_hline(data=xpal,aes(yintercept = mean(xpal$Tmax)), col="red", linetype=5)+
      geom_hline(data=xbien, aes(yintercept = mean(xbien$Tmax)), col="blue", linetype=5)+
      geom_vline(data=xpal, aes(xintercept = mean(xpal$Prec)), col="red", linetype=5)+
      geom_vline(data=xbien, aes(xintercept = mean(xbien$Prec)), col="blue", linetype=5)+
      theme_classic()+
      # theme(legend.position = "none")+
      coord_cartesian(ylim=c(min(niches_all$Tmax),max(niches_all$Tmax)), xlim=c(min(niches_all$Prec), max(niches_all$Prec)))+
      ggtitle(as.character(x$taxon[1]))
    
  )},.progress="text")
dev.off()

xpal <- niches_paleo_all[niches_paleo_all$taxon=="Alnus",]
xbien<- niches_bien_all[niches_bien_all$taxon=="Alnus",]


ggplot()+
  # stat_density2d(data=x, aes(x=Prec, y=Tmax), col="black")+
  stat_density2d(data=xpal, aes(x=Prec, y=Tmax), col="red")+
  stat_density2d(data=xbien, aes(x=Prec, y=Tmax), col="blue")+
  geom_point(data=xbien, aes(x=Prec, y=Tmax),alpha=0.05, col="blue", shape=21)+
  geom_point(data=xpal, aes(x=Prec, y=Tmax),alpha=0.05, col="red", shape=21)+
  scale_alpha_continuous()+
  scale_fill_gradientn(colours=cubeHelix(100))+
  theme_classic()+
  theme(legend.position = "none")
  # coord_cartesian(ylim=c(min(niches_all$Tmax),max(niches_all$Tmax)), xlim=c(min(niches_all$Prec),max(niches_all$Prec)))
  ggtitle(as.character(x$taxon[1]))

paleo <- 
  ggplot(data=x[x$time<1,], aes(x=time, y=Tmax))+
  # geom_point(alpha=0.3, shape=21)+
  stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
  # stat_density2d(data=x[x$time>0,],geom = "raster", aes(fill = ..density..), contour = FALSE)+
  geom_hline(aes(yintercept = mean(Tmax)), col="red")+
  geom_hline(aes(yintercept = mean(Tmax)-sd(Tmax)), col="red", linetype=5)+
  geom_hline(aes(yintercept = mean(Tmax)+sd(Tmax)), col="red", linetype=5)+
  stat_smooth(col="black")+
  scale_fill_gradientn(colours=cubeHelix(100))+
  # coord_cartesian(ylim=c(min(x$Tmax),max(x$Tmax)))+
  theme_classic()

modern <- 
  ggplot(data=x[x$time>0,], aes(x=time, y=Tmax))+
  # geom_point(alpha=0.3, shape=21)+
  stat_density2d(geom = "raster", aes(fill = ..density..), contour = FALSE)+
  # stat_density2d(data=x[x$time>0,],geom = "raster", aes(fill = ..density..), contour = FALSE)+
  geom_hline(aes(yintercept = mean(Tmax)), col="red")+
  geom_hline(aes(yintercept = mean(Tmax)-sd(Tmax)), col="red", linetype=5)+
  geom_hline(aes(yintercept = mean(Tmax)+sd(Tmax)), col="red", linetype=5)+
  stat_smooth(col="black")+
  scale_fill_gradientn(colours=cubeHelix(100))+
  coord_cartesian(ylim=c(min(x$Tmax),max(x$Tmax)))+
  theme_classic()

gridExtra::grid.arrange(paleo, modern, nrow=1)



#### compare summarised niches ####

Paleo_niches_summary_SAMPLED  <-read.csv('results/niches_paleo_summary_sampled.csv')[-1]
Paleo_niches_summary_SAMPLED$type <- "paleo_sampled"

Paleo_niches_summary          <-ddply(na.omit(niches_paleo_all), .(taxon), summarise, 
                                mean.Monthly.maxTemp = mean(Tmax),
                                sd.Monthly.maxTemp = sd(Tmax),
                                mean.Monthly.Prec = mean(Prec),
                                sd.Monthly.Prec = sd(Prec), .progress="text")
Paleo_niches_summary$type <- "paleo"

Bien_niches_summary     <-ddply(na.omit(niches_bien_all), .(taxon), summarise, 
                               mean.Monthly.maxTemp = mean(Tmax),
                               sd.Monthly.maxTemp = sd(Tmax),
                               mean.Monthly.Prec = mean(Prec),
                               sd.Monthly.Prec = sd(Prec), .progress="text")
Bien_niches_summary$type  <- "modern"

niches_summary_long <- rbind(Paleo_niches_summary, Bien_niches_summary[Bien_niches_summary$taxon %in% Paleo_niches_summary$taxon,])



ggplot(niches_summary_long, aes(x=mean.Monthly.maxTemp))+
  geom_histogram(aes(fill=type))+
  theme_classic()+
  geom_vline(xintercept= mean(niches_summary_long$mean.Monthly.maxTemp[niches_summary_long$type=="paleo"]), col='blue')+
  geom_vline(xintercept= mean(niches_summary_long$mean.Monthly.maxTemp[niches_summary_long$type=="modern"]), col='red')


plot(niches_summary_long$mean.Monthly.maxTemp[niches_summary_long$type=="paleo"], 
     niches_summary_long$mean.Monthly.maxTemp[niches_summary_long$type=="modern"])
abline(a=0, b=1)

niches_summary_wide <- data.frame(Bien_niches_summary[Bien_niches_summary$taxon %in% Paleo_niches_summary$taxon,], Paleo_niches_summary)
colnames(niches_summary_wide)[8:11] <- gsub(".1", "_paleo", colnames(niches_summary_wide)[8:11])
colnames(niches_summary_wide)[2:5] <- paste(colnames(niches_summary_wide)[2:5],"_modern", sep="")

library(smatr)

m_sma_tmx <- sma(mean.Monthly.maxTemp_modern ~ mean.Monthly.maxTemp_paleo, data=niches_summary_wide,  slope.test=1, elev.test=0)
m_sma_pre <- sma(mean.Monthly.Prec_modern    ~ mean.Monthly.Prec_paleo,    data=niches_summary_wide,  elev.test=0)
m_sma_tmx
m_sma_pre

ggplot(niches_summary_wide, aes(x=mean.Monthly.Prec_modern, y=mean.Monthly.Prec_paleo))+
  geom_errorbar( aes(ymin=mean.Monthly.Prec_paleo-sd.Monthly.Prec_paleo,   ymax=mean.Monthly.Prec_paleo + sd.Monthly.Prec_paleo), col="lightgrey")+
  geom_errorbarh(aes(xmin=mean.Monthly.Prec_modern-sd.Monthly.Prec_modern, xmax=mean.Monthly.Prec_modern + sd.Monthly.Prec_modern), col="lightgrey")+
  geom_point(shape=21, fill="white")+
  stat_smooth(method="lm")+
  geom_abline(slope=1, intercept=0, col="red")+
  theme_classic()+
  labs(title="Average monthly precipitation", x='Contemporary data (mm)',y='Paleo data (mm)')+
  coord_fixed(ratio=1)

ggplot(niches_summary_wide, aes(x=mean.Monthly.maxTemp_modern, y=mean.Monthly.maxTemp_paleo))+
  geom_errorbar( aes(ymin=mean.Monthly.maxTemp_paleo -sd.Monthly.maxTemp_paleo,      ymax=mean.Monthly.maxTemp_paleo + sd.Monthly.maxTemp_paleo), col="lightgrey")+
  geom_errorbarh(aes(xmin=mean.Monthly.maxTemp_modern-sd.Monthly.maxTemp_modern,     xmax=mean.Monthly.maxTemp_modern + sd.Monthly.maxTemp_modern), col="lightgrey")+
  geom_point(shape=21, fill="white")+
  stat_smooth(method="lm", fullrange = T)+
  geom_abline(slope=1, intercept=0, col="red")+
  theme_classic()+
  labs(x='Maximum temp. niche Contemporary data (°C)',y='Maximum temp. niche Paleo data (°C)')+
  coord_fixed(ratio=1)


#### sample, merge and summarize niches ####

species <- unique(niches_paleo_all$taxon)
nsamples=1000
niches_all_sampled <- data.frame("taxon"=rep(species, each=nsamples*2), "period"=rep(c("paleo", "modern"), each=nsamples), "Tmax"=NA, "Prec"=NA)

for(sp in species){
  print(sp)
  niches_paleo_all_SP <- niches_paleo_all_sampled[niches_paleo_all_sampled$taxon==sp,]
  niches_bien_all_SP <- niches_bien_all[niches_bien_all$taxon==sp,]
  
  niches_all_sampled[niches_all_sampled$taxon==sp, c("Tmax", "Prec")] <-
    rbind(
      niches_paleo_all_SP[sample(1:nrow(niches_paleo_all_SP),nsamples,replace=T), c("Tmax", "Prec")],
      niches_bien_all_SP[sample(1:nrow(niches_bien_all_SP),nsamples,replace=T), c("Tmax", "Prec")]
    )
  }

save(niches_all_sampled, file="results/niches_all_sampled.Rdata")
niches_all_summary_sampled <-ddply(niches_all_sampled, .(taxon), summarise, 
                                    mean.Monthly.maxTemp = mean(Tmax),
                                    sd.Monthly.maxTemp = sd(Tmax),
                                    mean.Monthly.Prec = mean(Prec),
                                    sd.Monthly.Prec = sd(Prec), .progress="text")
hist(niches_all_summary_sampled$mean.Monthly.maxTemp)



niches_all_summary <-ddply(niches_all, .(taxon), summarise, 
                                   mean.Monthly.maxTemp = mean(Tmax),
                                   sd.Monthly.maxTemp = sd(Tmax),
                                   mean.Monthly.Prec = mean(Prec),
                                   sd.Monthly.Prec = sd(Prec), .progress="text")
niches_all_summary$taxon <- as.factor(niches_all_summary$taxon)
niches_all_summary <- niches_all_summary[niches_all_summary$taxon %in% niches_all_summary_sampled$taxon,]
hist(niches_all_summary$mean.Monthly.maxTemp)


plot(niches_all_summary$mean.Monthly.maxTemp, niches_all_summary_sampled$mean.Monthly.maxTemp)
abline(0,1)

plot(niches_all_summary$mean.Monthly.maxTemp, Bien_niches_summary$mean.Monthly.maxTemp)
plot(niches_all_summary_sampled$mean.Monthly.maxTemp, Bien_niches_summary$mean.Monthly.maxTemp)

write.csv(niches_all_summary, 'results/niches_all_summary.csv',row.names=F)
write.csv(niches_all_summary_sampled, 'results/niches_all_summary_sampled.csv',row.names=F)
