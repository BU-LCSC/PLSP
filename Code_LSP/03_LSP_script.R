#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 03: A script for estimating phenometrics
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Submitted to the job for a single chunk for a site with shell script:
# #!/bin/bash
# echo Submitting $1
# R --vanilla < ~/03_LSP_script.R $1
#
# example submission command using default parameters:
# qsub -V -l h_rt=12:00:00 run_03.sh numSite chunk
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

numSite <- as.numeric(substr(args[3],1,3)) # site
cc      <- as.numeric(substr(args[3],4,6)) # chunk number; 1-200
# numSite <- 5; cc <- 50


########################################
## Load parameters
params <- fromJSON(file='~/PLSP_Parameters.json')
source(params$setup$rFunctions)


########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
print(ckDir)

## Load chunk image
ckNum <- sprintf('%03d',cc)
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)


## Load water mask
waterRater <- raster(paste0(params$setup$outDir,strSite,'/water_mask_30_1.tif'))

numCk <- params$setup$numChunks
chunk <- length(waterRater)%/%numCk
if(cc==numCk){
  chunks <- c((chunk*(cc-1)+1):length(waterRater))
}else{
  chunks <- c((chunk*(cc-1)+1):(chunk*cc))
}
waterMask <- values(waterRater)[chunks]


##########################################
# Estimate phenometrics
numPix <- dim(band1)[1]
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

pheno_mat <- matrix(NA,numPix,24*length(phenYrs))

for (i in 1:numPix){

  pheno_mat[i,] <- DoPhenologyPlanet(band1[i,],band2[i,],band3[i,],band4[i,],dates,phenYrs,params,waterMask[i])
  
  if(i%%10000==0) print(i)
}


# Save outputs
ckPheDir <- paste0(params$setup$outDir,strSite,'/chunk_phe')
if (!dir.exists(ckPheDir)) {dir.create(ckPheDir)}

save(pheno_mat,file=paste0(ckPheDir,'/chunk_phe_',ckNum,'.rda'))

