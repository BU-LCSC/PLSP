#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 01_1: An ancillary script for creating water mask
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


library(raster)
library(rgdal)
library(gdalUtils)

library(rjson)


###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
ck      <- as.numeric(substr(args[3],4,5))
# numSite <- 81; ck <- 1


###############################
# Load parameters and functions
params <- fromJSON(file='~/PLSP_Parameters.json')
source(params$setup$rFunctions)


########################################
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

imgDir <- paste0(params$setup$outDir,strSite,'/mosaic')
print(imgDir)

# Load base image
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))



########################################
## Water mask
# Load water mask data
# Pekel et al., 2016. High-resolution mapping of global surface water and its long-term changes. Nature, 540(7633), pp.418-422.
waterFiles <- list.files(params$setup$water,full.names=T)

waterTile <- c()
j <- 1
for(i in 1:length(waterFiles)){
  try({
    aa <- raster(waterFiles[i])

    pp <- as(extent(imgBase), "SpatialPolygons")
    crs(pp) <- crs(imgBase)
    pp <- spTransform(pp,crs(aa))

    bb <- intersect(aa,pp)
    if(length(bb)>0){
      waterTile[j] <- i
      j <- j + 1
    }
  },silent=TRUE)
}

print(length(waterTile))


##
if(length(waterTile)==1){
  waterRaster      <- raster(waterFiles[waterTile[1]])
  waterExtent      <- projectExtent(imgBase,crs(waterRaster))
  waterRaterCroped <- crop(waterRaster,waterExtent)
  waterReprojected <- projectRaster(waterRaterCroped,imgBase,method='ngb')
  
}else{
  waterRaster <- vector('list',length(waterTile))
  for(i in 1:length(waterTile)){
    tempRaster       <- raster(waterFiles[waterTile[i]])
    tempExtent       <- projectExtent(imgBase,crs(tempRaster))
    tempRasterCroped <- crop(tempRaster,tempExtent)

    waterRaster[[i]] <- tempRasterCroped
  }
  waterRaster$fun <- max
  waterRaster$na.rm <- T
  waterRaster <- do.call(mosaic,waterRaster)
  waterReprojected <- projectRaster(waterRaster,imgBase,method='ngb')
}



########################################
if(ck==34){
  cdLeng <- 3301:nrow(waterReprojected)
}else{
  cdLeng <- ((ck-1)*100+1):(ck*100)
}
# cdLeng <- 1

coords <- matrix(NA,(length(cdLeng)*ncol(waterReprojected)),2); cc <- 1 
for(j in cdLeng){
  for(i in 1:ncol(waterReprojected)){
    coords[cc,] <- c(waterReprojected@extent@xmin+((i-1)*3+1.5),waterReprojected@extent@ymax-((j-1)*3+1.5))
    cc <- cc + 1
  }
}
coords <- as.data.frame(coords)
xy     <- coords[,c(1,2)]
bb     <- SpatialPointsDataFrame(coords=xy,data=coords,proj4string=crs(waterReprojected))

waterSmoothedVal <- extract(waterReprojected,bb,buffer=15)


##
waterVal <- c()
for(j in 1:length(waterSmoothedVal)){
  vals <- waterSmoothedVal[[j]]
  waterVal <- c(waterVal,max(vals,na.rm=T))
  
  if(j%%10000==0) print(j)
}


## Save
outDir <- paste0(params$setup$outDir,strSite,'/chunk_water_30')
if (!dir.exists(outDir)) {dir.create(outDir)}  

numCk <- sprintf('%03d',ck)
outFile    <- paste0(outDir,'/chunk_water_',numCk,'.rda')
save(waterVal, file=outFile)



########################################
files <- list.files(outDir,pattern=glob2rx('*rda'),full.names=T)

if(length(files)==34){
  waterVals <- c()
  for(i in 1:length(files)){
    load(files[i])
    waterVals <- c(waterVals,waterVal)
  }
  
  ## Save
  print(length(waterVals)==length(imgBase))
  if(length(waterVals)==length(imgBase)){
    waterValsRaster <- setValues(imgBase,waterVals)
    
    outFile    <- paste0(params$setup$outDir,strSite,'/water_mask_30_1.tif')
    writeRaster(waterValsRaster, filename=outFile, format="GTiff", overwrite=TRUE)
  }  
}

