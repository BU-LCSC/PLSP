#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 02: A script for PlanetScope image process; save mosaiced images into chunks
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Submitted to the job for each site with shell script:
# #!/bin/bash
# echo Submitting $1
# R --vanilla < ~/02_make_chunks.R $1
#
# example submission command using default parameters:
# qsub -V -pe omp 28 -l h_rt=12:00:00 run_02.sh numSite
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(raster)
library(rgdal)
library(gdalUtils)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

########################################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 51


########################################
params <- fromJSON(file='~/PLSP_Parameters.json')
source(params$setup$rFunctions)


########################################
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

imgDir <- paste0(params$setup$outDir,strSite,'/mosaic')
print(imgDir)


########################################
dfiles <- list.files(path=imgDir,pattern=glob2rx('*mosaic.tif'))
files  <- list.files(path=imgDir,pattern=glob2rx('*mosaic.tif'),full.names=T)

# Dates
yy <- substr(dfiles,3,4)
mm <- substr(dfiles,5,6)
dd <- substr(dfiles,7,8)
dates <- as.Date(paste(mm,'/',dd,'/',yy,sep=''),'%m/%d/%y')

print(length(dates))


# Divide into chunks
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))

numCk <- params$setup$numChunks
chunk <- length(imgBase)%/%numCk

# Output directory
ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
if (!dir.exists(ckDir)) {dir.create(ckDir)}

# Directory for temporal outputs
ckDirTemp <- paste0(params$setup$outDir,strSite,'/chunk/temp')
if (!dir.exists(ckDirTemp)) {dir.create(ckDirTemp)}



########################################
# For each image corresponded to each date,
# divede images into chunks, and save them as temporal files
registerDoMC(params$setup$numCores)

foreach(i=1:length(dates)) %dopar%{
  band1 <- values(raster(files[i],1))
  band2 <- values(raster(files[i],2))
  band3 <- values(raster(files[i],3))
  band4 <- values(raster(files[i],4))
  
  foreach(cc=1:numCk) %dopar%{
    ckNum <- sprintf('%03d',cc)
    dirTemp <- paste0(ckDirTemp,ckNum)
    if (!dir.exists(dirTemp)) {dir.create(dirTemp)}
    
    if(cc==numCk){
      chunks <- c((chunk*(cc-1)+1):length(imgBase))
    }else{
      chunks <- c((chunk*(cc-1)+1):(chunk*cc))
    }
    b1 <- band1[chunks]
    b2 <- band2[chunks]
    b3 <- band3[chunks]
    b4 <- band4[chunks]  
    
    save(b1,b2,b3,b4,file=paste0(dirTemp,'/',yy[i],mm[i],dd[i],'.rda'))
  }
}


########################################
# Load files for each chunk, merge then, and save 
foreach(cc=1:numCk) %dopar%{
  ckNum <- sprintf('%03d',cc)
  dirTemp <- paste0(ckDirTemp,ckNum)
  files <- list.files(dirTemp,full.names=T)
  
  if(cc==numCk){
    chunks <- c((chunk*(cc-1)+1):length(imgBase))
  }else{
    chunks <- c((chunk*(cc-1)+1):(chunk*cc))
  }
  
  if(length(files)==length(dates)){
    band1 <- matrix(NA,length(chunks),length(dates))
    band2 <- matrix(NA,length(chunks),length(dates))
    band3 <- matrix(NA,length(chunks),length(dates))
    band4 <- matrix(NA,length(chunks),length(dates))
    for(i in 1:length(dates)){
      load(files[i])
    
      band1[,i] <- b1
      band2[,i] <- b2
      band3[,i] <- b3
      band4[,i] <- b4
    }
    # Save
    save(band1,band2,band3,band4,dates,
         file=paste0(ckDir,'/chunk_',ckNum,'.rda'))
  }else{
    print('No good!')
  }
}
  
  
########################################
## Remove temporary files
system(paste0('rm -r ',ckDirTemp))



