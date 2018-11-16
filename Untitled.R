# 
# Simple multivariate two sample test based on a certain distance metric 
# defined in the wavelet transformed representation of the two samples
# 

library(viridis)
library(plot3D)
library(mvtnorm)
library(wavethresh)
library(scales)

# This L2 norm is computed between the wavelet powerspectrum, so it should be insensitive to translations and rotations
computeWaveletPowerSpectrumL2norm <- function(test.data.3D.A, test.data.3D.B) {
  
  # 3d wavelet decomposition
  tdwd3D2A <- threshold.wd3D(wd3D(test.data.3D.A, filter.number=4), by.level=FALSE)
  tdwd3D2B <- threshold.wd3D(wd3D(test.data.3D.B, filter.number=4), by.level=FALSE)
  
  # Get the wavelet powerspectrum
  pwrA <- get3DwaveletPowerSpectrum(tdwd3D2A)
  pwrB <- get3DwaveletPowerSpectrum(tdwd3D2B)
  
  # Compute and return the distance
  return( sqrt( sum( pwrA$power - pwrB$power)^2 )  )
}

get3DwaveletPowerForLevel <- function(tdwd3D, level=0, newQuartzWindow=FALSE) {
  reconsForlevel <- wr3D(getWvForLevel(tdwd3D, level=level))
  return( sum(reconsForlevel^2) )
}

plot3DwaveletPowerSpectrum <- function(tdwd3D, newQuartzWindow=FALSE) {
  if(newQuartzWindow) {
    quartz()
  }
  mPwrSpec <- get3DwaveletPowerSpectrum(tdwd3D)
  plot(mPwrSpec[,1], mPwrSpec[,2], xlab="Wavelet Scale", ylab="Total Power", main="Total Power Per Scale", type="n")
  grid()
  points(mPwrSpec[,1], mPwrSpec[,2], pch=19, cex=0.5)
  lines(mPwrSpec[,1], mPwrSpec[,2], col="red")
}

get3DwaveletPowerSpectrum <- function(tdwd3D) {
  mPwrSpec <- matrix(ncol=2,nrow=nlevelsWT(tdwd3D))
  for(i in 1:nlevelsWT(tdwd3D)) {
    mPwrSpec[i,] <- c(i-1, get3DwaveletPowerForLevel(tdwd3D, level=(i-1)))
  }
  return(data.frame(scale=mPwrSpec[,1], power=mPwrSpec[,2]))
}

getWvForLevel <- function(wd3Dobj, levelToKeep, verbose=TRUE) {
  copywd3Dobj <- wd3Dobj
  for(i in 1:nlevelsWT(wd3Dobj)) {
    if( (i-1) != levelToKeep ) {
      if(verbose) {
        cat(paste("Inserting in level ",(i-1),"\n"))
      }
      dimX <- (2^(i-1))
      copywd3Dobj <- putD(copywd3Dobj, v=list(a=array(rep(0, times=(dimX^3)), dim=c(dimX, dimX, dimX)), lev=(i-1), block="GGG") )
      
      copywd3Dobj <- putD(copywd3Dobj, v=list(a=array(rep(0, times=(dimX^3)), dim=c(dimX, dimX, dimX)), lev=(i-1), block="GGH") )
      
      copywd3Dobj <- putD(copywd3Dobj, v=list(a=array(rep(0, times=(dimX^3)), dim=c(dimX, dimX, dimX)), lev=(i-1), block="GHG") )
      
      copywd3Dobj <- putD(copywd3Dobj, v=list(a=array(rep(0, times=(dimX^3)), dim=c(dimX, dimX, dimX)), lev=(i-1), block="GHH") )
      
      copywd3Dobj <- putD(copywd3Dobj, v=list(a=array(rep(0, times=(dimX^3)), dim=c(dimX, dimX, dimX)), lev=(i-1), block="HGG") )
      
      copywd3Dobj <- putD(copywd3Dobj, v=list(a=array(rep(0, times=(dimX^3)), dim=c(dimX, dimX, dimX)), lev=(i-1), block="HGH") )
      
      copywd3Dobj <- putD(copywd3Dobj, v=list(a=array(rep(0, times=(dimX^3)), dim=c(dimX, dimX, dimX)), lev=(i-1), block="HHG") )
      
      if( (i-1) == 0 ) {
        copywd3Dobj <- putD(copywd3Dobj, v=list(a=array(rep(0, times=(dimX^3)), dim=c(dimX, dimX, dimX)), lev=(i-1), block="HHH") )
      }
    }
  }
  return(copywd3Dobj)
}

plot3DwaveletLevel <- function(tdwd3D, level=0, newQuartzWindow=FALSE, logPlot=FALSE) {
  if(newQuartzWindow) {
    quartz()
  }
  nRowsInData <- nrow(tdwd3D$a)
  if(logPlot==FALSE) {
    slice3D(x=1:nRowsInData, y=1:nRowsInData, z=1:nRowsInData, colvar=wr3D(getWvForLevel(tdwd3D, level=level)), zs=c(1, nRowsInData/2), xs=c(1, nRowsInData/2), ys=c(nRowsInData/2, nRowsInData), col= viridis(100, alpha=0.5, option="B", direction=1), main=paste("3D Wavelet Level", level) )
  } else {
    slice3D(x=1:nRowsInData, y=1:nRowsInData, z=1:nRowsInData, colvar=log10(rescale(wr3D(getWvForLevel(tdwd3D, level=level)), to=c(5e-1, 1)) ), zs=c(1, nRowsInData/2), xs=c(1, nRowsInData/2), ys=c(nRowsInData/2, nRowsInData), col= viridis(100, alpha=0.5, option="B", direction=1), main=paste("3D Wavelet Level", level) )
  }
}

createToyExampleAndDecomposeAtAllLevels <- function(outFolder="~/Desktop/2018-09-Crete-CRP5/Projects/TwoSampleTestWaveletCoefficientDistance/Images/", nPixSide=64) {
  
  # generate some data
  test.data.3D <- array(rnorm(nPixSide*nPixSide*nPixSide)+5*exp(-((1:nPixSide-nPixSide/2)^2)/(2*10)), dim=c(nPixSide,nPixSide,nPixSide))
  
  # plot it
  png(paste(outFolder,"ToyNoisedLine3Dslices.png",sep=""))
  slice3D(x=1:nrow(test.data.3D),y=1:nrow(test.data.3D),z=1:nrow(test.data.3D),colvar=test.data.3D, zs=c(1, nrow(test.data.3D)/2), xs=c(1, nrow(test.data.3D)/2), ys=c(nrow(test.data.3D)/2,nrow(test.data.3D)), col= viridis(100, alpha=0.5, option="B", direction=1) )
  dev.off()
  
  # 3d wavelet decomposition
  tdwd3D2 <- wd3D(test.data.3D, filter.number=4)
  
  # plot it
  for(i in 1:nlevelsWT(tdwd3D2)) {
    png(paste(outFolder,"ToyNoisedLine3Dslices-res",(i-1),".png",sep=""), width=640, height=480)
    plot3DwaveletLevel(tdwd3D2, (i-1))
    dev.off()
  }
  
  # Power Spectrum Plot
  png(paste(outFolder,"ToyNoisedLine3Dslices-PowerSpectrum.png",sep=""))
  plot3DwaveletPowerSpectrum(tdwd3D2)
  dev.off()
}

createToyDemoForPowerL2norm <- function(nPixSide=64) {
  
  # generate some data
  test.data.3D_distributionA1 <- array(rnorm(nPixSide*nPixSide*nPixSide)+5*exp(-((1:nPixSide-nPixSide/2)^2)/(2*10)), dim=c(nPixSide,nPixSide,nPixSide))
  test.data.3D_distributionA2 <- array(rnorm(nPixSide*nPixSide*nPixSide)+5*exp(-((1:nPixSide-nPixSide/2)^2)/(2*10)), dim=c(nPixSide,nPixSide,nPixSide))
  test.data.3D_distributionB  <- array(rnorm(nPixSide*nPixSide*nPixSide), dim=c(nPixSide,nPixSide,nPixSide))
  
  distA1A2 <- computeWaveletPowerSpectrumL2norm(test.data.3D_distributionA1, test.data.3D_distributionA2)
  distA1B  <- computeWaveletPowerSpectrumL2norm(test.data.3D_distributionA1, test.data.3D_distributionB)
  
  # Report the result of the test
  cat(paste("\n------------------------------------------------\n"))
  cat(paste("Comparing two draws from the same distribution \n"))
  cat(paste("Distance : ",distA1A2,"\n\n"))
  cat(paste("Comparing two draws from different distributions \n"))
  cat(paste("Distance : ",distA1B,"\n"))
  cat(paste("------------------------------------------------\n"))
}