#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 04: A script for saving data layers into GeoTiff format
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

numSite <- as.numeric(args[3])
# numSite <- 2


########################################
params <- fromJSON(file='~/PLSP_Parameters.json')
source(params$setup$rFunctions)

productTable <- read.csv(params$setup$productTable,header=T,stringsAsFactors = F)
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr


########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckPheDir <- paste0(params$setup$outDir,strSite,'/chunk_phe_gap')
print(ckPheDir)

files <- list.files(path=ckPheDir,pattern=glob2rx('*.rda'),full.names=T)
print(length(files))



#Get all images to process
########################################
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks



# Save
pheDir <- paste0(params$setup$workDir,'Product_010/',strSite)
if (!dir.exists(pheDir)) {dir.create(pheDir)}


for(yToDo in 1:length(phenYrs)){
  
  l01 <- matrix(NA,numPix,1);l02 <- matrix(NA,numPix,1);l03 <- matrix(NA,numPix,1);l04 <- matrix(NA,numPix,1)
  l05 <- matrix(NA,numPix,1);l06 <- matrix(NA,numPix,1);l07 <- matrix(NA,numPix,1);l08 <- matrix(NA,numPix,1)
  l09 <- matrix(NA,numPix,1);l10 <- matrix(NA,numPix,1);l11 <- matrix(NA,numPix,1);l12 <- matrix(NA,numPix,1)
  l13 <- matrix(NA,numPix,1);l14 <- matrix(NA,numPix,1);l15 <- matrix(NA,numPix,1);l16 <- matrix(NA,numPix,1)
  l17 <- matrix(NA,numPix,1);l18 <- matrix(NA,numPix,1);l19 <- matrix(NA,numPix,1);l20 <- matrix(NA,numPix,1)
  l21 <- matrix(NA,numPix,1);l22 <- matrix(NA,numPix,1);l23 <- matrix(NA,numPix,1);l24 <- matrix(NA,numPix,1)
  
  for(i in 1:numChunks){
    cc <- sprintf('%03d',i)
    cfile <- paste0(ckPheDir,'/chunk_phe_',cc,'.rda') 
    log <- try(load(cfile),silent=F)
    if (inherits(log, 'try-error')) next 
    
    if(i==numChunks){chunks <- c((chunk*(i-1)+1):numPix)
    }else{chunks <- c((chunk*(i-1)+1):(chunk*i))}
    
    chunkStart <- chunks[1];  chunkEnd <- chunks[length(chunks)]
    
    l01[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+1)];l02[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+2)];l03[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+3)]
    l04[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+4)];l05[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+5)];l06[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+6)]
    l07[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+7)];l08[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+8)];l09[chunkStart:chunkEnd,] <- pheno_mat[, (24*(yToDo-1)+9)]
    l10[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+10)];l11[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+11)];l12[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+12)]
    l13[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+13)];l14[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+14)];l15[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+15)]
    l16[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+16)];l17[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+17)];l18[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+18)]
    l19[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+19)];l20[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+20)];l21[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+21)]
    l22[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+22)];l23[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+23)];l24[chunkStart:chunkEnd,] <- pheno_mat[,(24*(yToDo-1)+24)]
    
  }
  
  r01 <- setValues(imgBase,l01); r02 <- setValues(imgBase,l02); r03 <- setValues(imgBase,l03); r04 <- setValues(imgBase,l04)
  r05 <- setValues(imgBase,l05); r06 <- setValues(imgBase,l06); r07 <- setValues(imgBase,l07); r08 <- setValues(imgBase,l08)
  r09 <- setValues(imgBase,l09); r10 <- setValues(imgBase,l10); r11 <- setValues(imgBase,l11); r12 <- setValues(imgBase,l12)
  r13 <- setValues(imgBase,l13); r14 <- setValues(imgBase,l14); r15 <- setValues(imgBase,l15); r16 <- setValues(imgBase,l16)
  r17 <- setValues(imgBase,l17); r18 <- setValues(imgBase,l18); r19 <- setValues(imgBase,l19); r20 <- setValues(imgBase,l20)
  r21 <- setValues(imgBase,l21); r22 <- setValues(imgBase,l22); r23 <- setValues(imgBase,l23); r24 <- setValues(imgBase,l24)
  
  
  # Save
  pheDirYear <- paste0(pheDir,'/',phenYrs[yToDo])
  if (!dir.exists(pheDirYear)) {dir.create(pheDirYear)}
  
  writeRaster(r01,filename=paste0(pheDirYear,'/01_',phenYrs[yToDo],'_',productTable$short_name[ 1],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r02,filename=paste0(pheDirYear,'/02_',phenYrs[yToDo],'_',productTable$short_name[ 2],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r03,filename=paste0(pheDirYear,'/03_',phenYrs[yToDo],'_',productTable$short_name[ 3],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r04,filename=paste0(pheDirYear,'/04_',phenYrs[yToDo],'_',productTable$short_name[ 4],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r05,filename=paste0(pheDirYear,'/05_',phenYrs[yToDo],'_',productTable$short_name[ 5],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r06,filename=paste0(pheDirYear,'/06_',phenYrs[yToDo],'_',productTable$short_name[ 6],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r07,filename=paste0(pheDirYear,'/07_',phenYrs[yToDo],'_',productTable$short_name[ 7],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r08,filename=paste0(pheDirYear,'/08_',phenYrs[yToDo],'_',productTable$short_name[ 8],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r09,filename=paste0(pheDirYear,'/09_',phenYrs[yToDo],'_',productTable$short_name[ 9],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r10,filename=paste0(pheDirYear,'/10_',phenYrs[yToDo],'_',productTable$short_name[10],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r11,filename=paste0(pheDirYear,'/11_',phenYrs[yToDo],'_',productTable$short_name[11],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r12,filename=paste0(pheDirYear,'/12_',phenYrs[yToDo],'_',productTable$short_name[12],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r13,filename=paste0(pheDirYear,'/13_',phenYrs[yToDo],'_',productTable$short_name[13],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r14,filename=paste0(pheDirYear,'/14_',phenYrs[yToDo],'_',productTable$short_name[14],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r15,filename=paste0(pheDirYear,'/15_',phenYrs[yToDo],'_',productTable$short_name[15],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r16,filename=paste0(pheDirYear,'/16_',phenYrs[yToDo],'_',productTable$short_name[16],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r17,filename=paste0(pheDirYear,'/17_',phenYrs[yToDo],'_',productTable$short_name[17],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r18,filename=paste0(pheDirYear,'/18_',phenYrs[yToDo],'_',productTable$short_name[18],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r19,filename=paste0(pheDirYear,'/19_',phenYrs[yToDo],'_',productTable$short_name[19],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r20,filename=paste0(pheDirYear,'/20_',phenYrs[yToDo],'_',productTable$short_name[20],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r21,filename=paste0(pheDirYear,'/21_',phenYrs[yToDo],'_',productTable$short_name[21],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r22,filename=paste0(pheDirYear,'/22_',phenYrs[yToDo],'_',productTable$short_name[22],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r23,filename=paste0(pheDirYear,'/23_',phenYrs[yToDo],'_',productTable$short_name[23],'.tif'), format="GTiff", overwrite=TRUE)
  writeRaster(r24,filename=paste0(pheDirYear,'/24_',phenYrs[yToDo],'_',productTable$short_name[24],'.tif'), format="GTiff", overwrite=TRUE)
  
  print(yToDo)
}









