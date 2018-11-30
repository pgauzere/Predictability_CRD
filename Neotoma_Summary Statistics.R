count_crossings <- function(x_asymptote, segments) # assumes segments is in the form [x0 y0 x1 y1]
{
  x0 <- segments[,1]
  x1 <- segments[,3]
  
  xsmall <- apply(segments, 1, function(row) { min(c(row[1],row[3]),na.rm=T)  })
  xbig <- apply(segments, 1, function(row) { max(c(row[1],row[3]),na.rm=T)  })

  isin <- (x_asymptote >= xsmall) & (x_asymptote <= xbig)

  count <- sum(as.numeric(isin), na.rm=T)
  
  return(count)
}

allcrossings <- function(ft, ct, length.out=100)
{
	df <- na.omit(cbind(ft, ct))

	if (nrow(df)<=1)
	{
	  return(data.frame(seq_x=NA,counts=NA))
	}
	else
	{
	  segments <- matrix(NA, nrow=length(ft)-1,ncol=4)
	  for (i in 1:(length(ft)-1))
	  {
	    segments[i,] <- c(ft[i],ct[i],ft[i+1],ct[i+1])
	  }	
	  
	  x0 <- segments[,1]
	  x1 <- segments[,3]
	  
	  seq_x <- seq(min(c(x0,x1),na.rm=T), max(c(x0,x1),na.rm=T), length.out=length.out)
	  
	  counts <- sapply(seq_x, count_crossings, segments=segments)
	  
	  return(data.frame(seq_x,counts))
	}
}

summarizeftct <- function(ft, ct, time, name)
{	
  
  num.points = length(ft)
  num.points.not.na <- length(na.omit(ft))
  
  crossings.distribution <- allcrossings(ft, ct)
  
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
  lm.r2 <- NA
  lm.rho <- NA
  if (num.points.not.na >= 3)
  {
    m_lm <- lm(ct~ft)
    lm.intercept <- coef(m_lm)[1]
    lm.slope <- coef(m_lm)[2]
    
    lm.r2 <- summary(m_lm)$r.squared
    lm.rho <- cor(ft,ct,use="complete")
  }
	
	return(data.frame(name=as.character(name), crossings.mean=crossings.mean, crossings.max=crossings.max, deviation.base=deviation.base, deviation.absolute=deviation.absolute,
	                  deviation.change=deviation.change, deviation.change.1Ky=deviation.change.1Ky,num.points= num.points, num.points.not.na= num.points.not.na, lm.intercept=lm.intercept, lm.slope=lm.slope, lm.r2=lm.r2, lm.rho=lm.rho))
}
