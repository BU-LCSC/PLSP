#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 01: A script for PlanetScope image process
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 24


###############################
params <- fromJSON(file='~/PLSP_Parameters.json')
source(params$setup$rFunctions)

params$setup$neon  <- FALSE
params$setup$amflx <- FALSE


########################################
## Get site name, image directory and coordinate
geojsonDir <- params$setup$geojsonDir

siteInfo <- GetSiteInfo(numSite,geojsonDir,params)

imgDir <- siteInfo[[1]]; strSite <- siteInfo[[2]]
print(paste(strSite,';',imgDir))

cLong <- siteInfo[[3]];cLat <- siteInfo[[4]]
print(paste(cLong,';',cLat))



########################################
## List of files
dfileSR   <- list.files(path=imgDir,pattern=glob2rx('*MS_SR*.tif'),recursive=T)
dfileUDM  <- list.files(path=imgDir,pattern=glob2rx('*_DN_udm*.tif'),recursive=T)
dfileUDM2 <- list.files(path=imgDir,pattern=glob2rx('*_udm2*.tif'),recursive=T)

fileSR   <- list.files(path=imgDir,pattern=glob2rx('*MS_SR*.tif'),recursive=T,full.names=T)
fileUDM  <- list.files(path=imgDir,pattern=glob2rx('*_DN_udm*.tif'),recursive=T,full.names=T)
fileUDM2 <- list.files(path=imgDir,pattern=glob2rx('*_udm2*.tif'),recursive=T,full.names=T)


# Dates
# yy <- substr(dfileSR,8,9)
# mm <- substr(dfileSR,10,11)
# dd <- substr(dfileSR,12,13)
yy <- substr(dfileSR,58,59)
mm <- substr(dfileSR,60,61)
dd <- substr(dfileSR,62,63)
dates_all <- as.Date(paste(mm,'/',dd,'/',yy,sep=''),'%m/%d/%y')
dates <- unique(dates_all)

print(length(dates))



####
## Image process
# Create output directory for base image
outDir <- paste0(params$setup$outDir,strSite)
if (!dir.exists(outDir)) {dir.create(outDir)}  

## Create site shapefile and base image
siteWin <- GetSiteShp(fileSR,cLong,cLat)
imgBase <- GetBaseImg(fileSR,siteWin,outDir,save=T)


##
registerDoMC(params$setup$numCores)

# Create output directory for mosaiced images
outDir <- paste0(params$setup$outDir,strSite,'/mosaic')
if (!dir.exists(outDir)) {dir.create(outDir)}  

foreach(dd=1:length(dates)) %dopar% {
  
  # Find images for a specific date
  imgMulti <- which(substr(dfileSR,56,63)==paste0(substr(dates[dd],1,4),substr(dates[dd],6,7),substr(dates[dd],9,10)))
  
  # Find images that have all 4 PlanetSceop bands
  imgVaild <- c()
  for(mm in 1:length(imgMulti)){
    log <- try({
      img <- raster(fileSR[imgMulti[mm]])
      img <- crop(img,siteWin)
      },silent=T)
    if(inherits(log,'try-error')){ 
      next  
    }else{
      numBand <- nbands(raster(fileSR[imgMulti[mm]]))
      if(numBand==4){
        imgVaild <- c(imgVaild,imgMulti[mm])   
      }
    }
  }
  
  # If the number of images that have 4 bands, load them and create a mosaiced image 
  if(length(imgVaild)>0){
    
    imgB <- vector('list',length(imgVaild))
    
    for(mm in 1:length(imgVaild)){
      ii <- imgVaild[mm]
      
      img <- raster(fileSR[ii])
      numBand <- nbands(img)
      
      str <- substr(dfileSR[ii],56,77)
      
      imgP <- vector('list',numBand)
      for(i in 1:numBand){
        imgT <- raster(fileSR[ii],band=i)
        imgT <- crop(imgT,siteWin)
        
        if(length(which(substr(dfileUDM,56,77)==str))==1 & length(which(substr(dfileUDM2,56,77)==str))==0){
          log <- try({
            udmT <- raster(fileUDM[which(substr(dfileUDM,56,77)==str)])
            udmT <- crop(udmT,siteWin)
          },silent=T)
          if(inherits(log,'try-error')){ 
            next
          }else{
            imgT[udmT>0] <- NA  
          }
        }else if(length(which(substr(dfileUDM,56,77)==str))==1 & length(which(substr(dfileUDM2,56,77)==str))==1){
          log <- try({
            udmT  <- raster(fileUDM[which(substr(dfileUDM,56,77)==str)])
            udmT  <- crop(udmT,siteWin)
          },silent=T)
          if(inherits(log,'try-error')){ 
            imgT <- imgT
          }else{
            imgT[udmT>0] <- NA
          }
          log <- try({
            udm2T <- raster(fileUDM2[which(substr(dfileUDM2,56,77)==str)])
            udm2T <- crop(udm2T,siteWin)  
          },silent=T)
          if(inherits(log,'try-error')){ 
            imgT <- imgT
          }else{
            imgT[udm2T!=1] <- NA
          }
        }
        imgP[[i]] <- imgT
      }
      
      imgB[[mm]] <- brick(imgP)
    }
    
    temp1 <- vector('list',(length(imgVaild)+1))
    temp2 <- vector('list',(length(imgVaild)+1))
    temp3 <- vector('list',(length(imgVaild)+1))
    temp4 <- vector('list',(length(imgVaild)+1))
    for(i in 1:length(imgVaild)){
      temp1[[i]] <- raster(imgB[[i]],1)
      temp2[[i]] <- raster(imgB[[i]],2)
      temp3[[i]] <- raster(imgB[[i]],3)
      temp4[[i]] <- raster(imgB[[i]],4)
    }
    temp1[[(length(imgVaild)+1)]] <- imgBase
    temp2[[(length(imgVaild)+1)]] <- imgBase
    temp3[[(length(imgVaild)+1)]] <- imgBase
    temp4[[(length(imgVaild)+1)]] <- imgBase
    
    for(i in 1:length(imgVaild)){
      log <- try(compareRaster(temp1[[i]],imgBase,extent=F,rowcol=F),silent=T)
      if(inherits(log,'try-error')){
        temp1[[i]] <- projectRaster(temp1[[i]],imgBase)    
      }
      log <- try(compareRaster(temp2[[i]],imgBase,extent=F,rowcol=F),silent=T)
      if(inherits(log,'try-error')){
        temp2[[i]] <- projectRaster(temp2[[i]],imgBase)    
      }
      log <- try(compareRaster(temp3[[i]],imgBase,extent=F,rowcol=F),silent=T)
      if(inherits(log,'try-error')){
        temp3[[i]] <- projectRaster(temp3[[i]],imgBase)    
      }
      log <- try(compareRaster(temp4[[i]],imgBase,extent=F,rowcol=F),silent=T)
      if(inherits(log,'try-error')){
        temp4[[i]] <- projectRaster(temp4[[i]],imgBase)    
      }
    }
    temp1$fun <- mean; temp2$fun <- mean; temp3$fun <- mean; temp4$fun <- mean
    temp1$na.rm <- T; temp2$na.rm <- T; temp3$na.rm <- T; temp4$na.rm <- T
    rast1 <- do.call(mosaic,temp1)  
    rast2 <- do.call(mosaic,temp2)  
    rast3 <- do.call(mosaic,temp3)  
    rast4 <- do.call(mosaic,temp4)  
    
    # Brick bands
    Rast <- brick(rast1,rast2,rast3,rast4)
    
    # Save
    
    outFile    <- paste0(outDir,'/',substr(dates[dd],1,4),substr(dates[dd],6,7),substr(dates[dd],9,10),'_cliped_mosaic.tif')
    writeRaster(Rast, filename=outFile, format="GTiff", overwrite=TRUE)

    print(outFile)
  } 
} 

# Check the length of output
print(length(list.files(path=outDir)))

print(length(dates))
