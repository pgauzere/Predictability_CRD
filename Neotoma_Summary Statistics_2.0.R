allcrossings <- function(ft, ct, ct_error, ft_error, length.out=100)
{
  df <- na.omit(cbind(ft, ct, ct_error, ft_error))
  
  # plot<-
  # ggplot(as.data.frame(df), aes(y=ct, x=ft))+
  #   geom_path()+
  #   geom_errorbar(aes(ymin=ct - ct_error, ymax=ct+ct_error))
  
  if (nrow(df)<=1)
  {
    return(data.frame(seq_x=NA,counts=NA))
  }
  else
  {
    segments <- matrix(NA, nrow=length(ft)-1,ncol=8)
    colnames(segments)<- c('x0', 'y0', 'x1', 'y1', 'x0_error', 'y0_error', 'x1_error', 'y1_error')
    for (i in 1:(length(ft)-1))
    {
      segments[i,] <- c(ft[i], ct[i], ft[i+1], ct[i+1], ft_error[i], ct_error[i], ft_error[i+1], ct_error[i+1])
    }	
    
    x0 <- segments[,1]
    x1 <- segments[,3]
    
    seq_x <- seq(min(c(x0,x1),na.rm=T), max(c(x0,x1),na.rm=T), length.out=length.out)
    
    counts <- sapply(seq_x, function(x_asymptote){
      
      
      # x_asymptote <- seq_x[60]
      # plot2 <-
      #   plot+
      #   geom_vline(xintercept = x_asymptote)
      
      xsmall <- apply(segments, 1, function(row) { min(c(row[1],row[3]),na.rm=T)  })
      xbig <- apply(segments, 1, function(row) { max(c(row[1],row[3]),na.rm=T)  })
      
      isin <- (x_asymptote >= xsmall) & (x_asymptote <= xbig)
      count <- sum(as.numeric(isin), na.rm=T)
      
      if(count>1){
        segs <- segments[isin==T  ,1:4]
        errors <- segments[isin==T,5:8]
        
        testDiff <- rep("NA", count-1)
        
        for(i in seq(1:(count-1))){
          
          segA<- c(
            'xmin' = min((segs-errors)[i,c(1,3)]),
            'ymin' = min((segs-errors)[i,c(2,4)]),
            'xmax' = max((segs+errors)[i,c(1,3)]),
            'ymax' = max((segs+errors)[i,c(2,4)]), 
            'ymean'= mean(segs[i,c(2,4)]))
          
          segB<- c(
            'xmin' = min((segs-errors)[i+1,c(1,3)]),
            'ymin' = min((segs-errors)[i+1,c(2,4)]),
            'xmax' = max((segs+errors)[i+1,c(1,3)]),
            'ymax' = max((segs+errors)[i+1,c(2,4)]),
            'ymean'= mean(segs[i+1,c(2,4)]))
          
          segA_upper_segment <-  segA['ymean'] > segB['ymean']
          
          # plot2+
          #   geom_hline(yintercept=segA['ymin'], linetype=3, col="red")+
          #   geom_hline(yintercept=segA['ymax'], linetype=3, col="red")+
          #   geom_hline(yintercept=segA['ymean'], col="red")+
          #   geom_hline(yintercept=segB['ymin'], linetype=3, col="blue")+
          #   geom_hline(yintercept=segB['ymax'], linetype=3, col="blue")+
          #   geom_hline(yintercept=segB['ymean'], col="blue")
          
          if(segA_upper_segment==T){
            testDiff[i] <- segA['ymin'] > segB['ymax']
          } else{
            testDiff[i] <- segA['ymax'] < segB['ymin']
          }          
        }
        
        return(count - length(which(testDiff==FALSE)))
        # return(testDiff)
        
      } else{
        return(count)
      }})
    return(data.frame(seq_x,counts))
  }}

summarizeftct <- function(ft, ct, ft_error, ct_error, time, name)
{	
  
  num.points = length(ft)
  num.points.not.na <- length(na.omit(ft))
  
  crossings.distribution <- allcrossings(ft, ct, ft_error, ct_error)
  
  crossings.mean <- mean(crossings.distribution$counts)
  crossings.max <- max(crossings.distribution$counts)
  
  deviation <- ft-ct
  deviation.change  <- rep(NA,num.points) 
  for(i in 2:num.points){deviation.change[i] <- abs(deviation[i]) - abs(deviation[i-1])}
  
  deviation.base       <- abs(deviation[1])
  deviation.absolute   <- abs(mean(deviation,na.rm=T))
  deviation.change     <- abs(deviation[num.points]) - abs(deviation[1])
  deviation.change.1Ky <- ((abs(deviation[num.points]) - abs(deviation[1]))/ (abs(time[num.points] - time[1])))*1000
  
  
  lm.intercept <- NA
  lm.slope <- NA
  lm.se <- NA
  lm.r2 <- NA
  lm.rho <- NA
  if (num.points.not.na >= 3)
  {
    m_lm <- lm(ct~ft)
    lm.intercept <- coef(m_lm)[1]
    lm.slope <- coef(m_lm)[2]
    lm.se <- summary(m_lm)$coefficient[2,2]
    lm.pval <- summary(m_lm)$coefficient[2,4]
    lm.r2 <- summary(m_lm)$r.squared
    lm.rho <- cor(ft,ct,use="complete")
  }
	
	return(data.frame(name=as.character(name), crossings.mean=crossings.mean, crossings.max=crossings.max, deviation.base=deviation.base, deviation.absolute=deviation.absolute,
	                  deviation.change=deviation.change, deviation.change.1Ky=deviation.change.1Ky,num.points= num.points, num.points.not.na= num.points.not.na, 
	                  lm.intercept=lm.intercept, lm.slope=lm.slope, 
	                  lm.se=lm.se, 
	                  lm.pval=lm.pval,
	                  lm.r2=lm.r2, lm.rho=lm.rho))
}

