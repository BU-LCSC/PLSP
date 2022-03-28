library(raster)
library(rgdal)
library(gdalUtils)
library(rjson)
library(scales)

numSite <- 32


########################################
params <- fromJSON(file='~/PLSP_Parameters.json')
source(params$setup$rFunctions)


########################################
## Get site name and image directory
geojsonDir <- params$setup$geojsonDir

strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0('~/chunks/',strSite)
print(ckDir)



########################################
## 
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks

cx <- (imgBase@extent@xmax+imgBase@extent@xmin)/2
cy <- (imgBase@extent@ymax+imgBase@extent@ymin)/2

# Point Shape file
utm_crs = imgBase@crs
site <- data.frame(1,cx,cy)
colnames(site) <- c('id','xcrd','ycrd')
xy   <- site[,c(2,3)]
bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=utm_crs)

pixNum <- extract(imgNum,bb)


##
ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)

blue  <- band1[pixNum%%chunk,]
green <- band2[pixNum%%chunk,]
red   <- band3[pixNum%%chunk,]
nir   <- band4[pixNum%%chunk,]
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

blue <- blue/10000; green <- green/10000; red <- red/10000; nir <- nir/10000
  
pheno_pars <- params$phenology_parameters
qa_pars    <- params$qa_parameters
  
VI     <- 2.5*(nir - red) / (nir + 2.4*red + 1)
vi     <- 2.5*(nir - red) / (nir + 2.4*red + 1)

######################################################
## Spline

# Spikes
spikes <- CheckSpike_MultiBand(blue, red, vi, dates, pheno_pars)
vi[spikes] <- NA
  
#
dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
vi_dorm <- quantile(vi[dormIms & vi>0],probs=pheno_pars$dormantQuantile,na.rm=T)   #Calc vi dormant value
  
vi[vi < vi_dorm] <- vi_dorm
  
# Gap fill
splineStart <- as.Date(as.Date(paste0(phenYrs,'-01-01')) - pheno_pars$splineBuffer) 
numDaysFit  <-  365 + (pheno_pars$splineBuffer * 2)    
splineEnd <- splineStart+(numDaysFit-1)
  
daysVec <- 1:numDaysFit
inYear <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)

#Determine gaps that require filling
gDates <- dates[!is.na(vi)]  
dDiff <- diff(gDates) > pheno_pars$gapLengthToFill     #Gaps greater than 20 days will be filled
dStart <- gDates[c(dDiff,FALSE)]
dEnd <- gDates[c(FALSE,dDiff)]

#Locate gaps in date vector
all_dates <- seq(min(splineStart), max(splineEnd), by="day")
  
fill_locations <- matrix(FALSE,length(all_dates))
for (d in 1:length(dStart)) {
  fill_locations[all_dates >= dStart[d] & all_dates < dEnd[d]] <- TRUE}
  
  fill_dates <- all_dates[fill_locations]
  
yToDo <- 1:length(phenYrs)
yrsWithGaps <- c()
for (y in yToDo) {
  pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
  if (sum(pred_dates %in% fill_dates) > 0) {yrsWithGaps <- c(yrsWithGaps,phenYrs[y])}
}
  
#If there are gaps to be filled, then we will spline all years. 
#If not, just spline product years
yrs <- phenYrs
  
numYrs <- length(yrs)
daysVec <- 1:numDaysFit
vecLength <- numDaysFit*numYrs
  
#First, we will fit splines to each year invidually
#To line up observations from each year, we will create a matrix for vi and each band (numDaysFit x numYears)
smoothMat <- matrix(NA, numDaysFit, numYrs)
maskMat <- matrix(0, numDaysFit, numYrs)
fillMat <- smoothMat
baseWeights <- maskMat
  
for (y in 1:numYrs) {
    #Use try statement, because we don't want to stop processing if only an error in one year
    try({
      
      dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(vi)   
      
      dateSub <- dates[dateRange]
      viSub <- vi[dateRange]
      
      #Get weights
      weights <- matrix(1,length(dateSub))
      
      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
      
      #Assign weights and run cubic spline
      smoothed <- Smooth_VI(viSub, dateSub, pred_dates, weights, pheno_pars, vi_dorm)
      
      
      #Mask spline in gaps, and before/after first/last image
      maskMat[fill_locations[all_dates %in% pred_dates],y] <- 1    #Mask spline in gaps
      maskMat[pred_dates < dateSub[1],y] <- 1                      #Mask spline before first image and after last image
      maskMat[pred_dates > dateSub[length(dateSub)],y] <- 1
      
      #Mask spline in the buffer years (only interested in comparing splines in target year)
      maskMat[format(pred_dates,'%Y') != yrs[y],y]  <- 1
      
      fillDs <- pred_dates %in% dateSub
      
      smoothMat[,y] <- smoothed
      baseWeights[fillDs,y] <- weights
      fillMat[fillDs,y] <- viSub
      
    },silent=TRUE)
  }
  
xs <- rep(daysVec,numYrs)
ys <- matrix(fillMat,vecLength)
ysGood <- !is.na(ys)
baseW <- matrix(baseWeights,vecLength)   #Base Weights are 1=clear observation
  
smoothMat_Masked <- smoothMat
maskMat <- as.logical(maskMat)
smoothMat_Masked[maskMat] <- NA
  
#Loop through years, compare spline to other years, weight each year based on similarity, fit spline, calculate phenology
weightArray <- calculateWeights(smoothMat_Masked, numDaysFit, numYrs, pheno_pars) 
  
prevYear <- daysVec <= pheno_pars$splineBuffer
inYear <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)
nextYear <- daysVec > (pheno_pars$splineBuffer+365)

for(y in 1:5){
  pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
  
  
  if (yrs[y] %in% yrsWithGaps) {
    
    indPrev <- y-1; indPrev[indPrev<1] <- 1
    indNext <- y+1; indNext[indNext>numYrs] <- numYrs
    
    weights <- rbind(weightArray[prevYear,,indPrev],
                     weightArray[inYear,,y],
                     weightArray[nextYear,,indNext])
    
    #Where are the gaps?
    toFill <- fill_locations[all_dates %in% pred_dates]
    
    weights[!toFill,] <- 0     #Set weight to zero for observations that aren't in a gap
    weights[,y] <- 1           #Set weights in target year to 1
    
    
    #Now that we have weights, calculate phenology
    #######################
    weights <- matrix(weights,vecLength) * baseW   #Multiple weights by base weight
    theInds <- ysGood & weights > 0
    xs_sub <- xs[theInds]; w_sub <- weights[theInds]
    smoothed_vi <- Smooth_VI(ys[theInds], xs_sub, daysVec, w_sub, pheno_pars, vi_dorm)  #Fit spline
    
  } else {
    
    #Variables needed for next steps if the above gap filling was not done
    theInds <- matrix(FALSE,length(ysGood))
    theInds[((y-1)*numDaysFit+1):(y*numDaysFit)] <- TRUE
    xs_sub <- xs[theInds]; w_sub <- baseW[theInds]
    
    smoothed_vi <- smoothMat[,y]   #if no gaps to fill, just use existing spline
  }
  
  
  # Number of clear observation
  filled_vi <- fillMat[,y]
  filled_vi[baseWeights[,y] < 1] <- NA    #If weight is less than 1, implies it is a snow-fill, and we don't want to count snow-filled as a valid observation. So set to NA.
  numObs <- sum(!is.na(filled_vi) & inYear)   #Number of observations in year
  
  #
  viSub   <- filled_vi[!is.na(filled_vi)]
  dateSub <- pred_dates[!is.na(filled_vi)]
}
  
 

######################################################
par(oma=c(2,2,2,2),mar=c(4,4,2,4),mgp=c(2.5,1,0))
plot(dates,VI,
     xlim=c(as.Date('2017-01-15'),as.Date('2021-12-15')),
     ylim=c(0,0.315),pch=19,col=rgb(0,0,0,0.2),cex=0.8,
     cex.lab=1.5,cex.axis=1.3,
     xlab='Dates',ylab='EVI2')
for(y in 1:5){
  all_d <- seq(min(splineStart[y]), max(splineEnd[y]), by="day")
  points(all_d[186:550],smoothMat[186:550,y],typ='l',lwd=4)
}

