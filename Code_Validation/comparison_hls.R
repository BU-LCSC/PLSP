#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# An script for comparison with MSLSP data
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


library(raster)
library(rgdal)
library(gdalUtils)

library(RColorBrewer)
library(viridis)

########################################
args <- commandArgs()
print(args)

ss <- as.numeric(substr(args[3],1,3)) # site
yy <- as.numeric(substr(args[3],4,7)) # year
vv <- as.numeric(substr(args[3],8,9)) # variable
# ss <- 13; yy <- 2017; vv <- 1

vari <- c('OGI','50PCGI','OGMx','OGD','50PCGD','OGMn')


options(warn=-1)


########################################
imgs <- list.files(params$setup$workDir,pattern=glob2rx(paste0('*',yy,'_',vari[vv],'.tif')),full.names=T,recursive=T)

rastP <- raster(imgs[ss])

site <- unlist(strsplit(imgs[ss],'/'))[6]
print(site)

pp <- as(extent(rastP), "SpatialPolygons")
crs(pp) <- crs(rastP)
pts <- spsample(pp, n = 3000, type = "random")



########################################
shpHLS <- shapefile('~/sentinel2_tiles_north_america_Albers.shp')

ptsHLS <- spTransform(pts,crs(shpHLS))
tileHLS <- intersect(shpHLS,ptsHLS)

if(length(tileHLS$Name)>0){
  path <- paste0('~/product_v011/',tileHLS$Name[1])
  file <- list.files(path,pattern=glob2rx(paste0('*',yy,'*')),full.names=T)
  
  baseImg <- raster(paste0('~/AWS_results/tiles/2017/',tileHLS$Name[1],'/LSP_',tileHLS$Name[1],'_2017_50PCGD.tif'))
  
  rastH <- raster(file,varname=vari[vv])
  ptsH <- spTransform(ptsHLS,crs(baseImg))
  
  
  ########################################
  valH2 <- extract(rastH,ptsH,buffer=45)
  valP2 <- extract(rastP,pts,buffer=45)
  
  dat <- matrix(NA,length(pts),4)
  for(i in 1:length(pts)){
    dat[i,7]  <- round(median(valH2[[i]],na.rm=T))
    dat[i,8]  <- round(median(valP2[[i]],na.rm=T))
    dat[i,9]  <- round(mean(valH2[[i]],na.rm=T))
    dat[i,10] <- round(mean(valP2[[i]],na.rm=T))
  }
  
  save(dat,file=paste0('~/data/comp_hls/ext/',vari[vv],'_',yy,'_',site,'.rda'))
}




########################################
vari <- c('OGI','50PCGI','OGMx','OGD','50PCGD','OGMn')
path <- '~/data/comp_hls/ext'

pheVal <- vector('list',length(vari))
for(vv in 1:6){

  gup <- list.files(path,pattern=glob2rx(paste0(vari[vv],'*')),full.names=T)

  temp <- c()
  for(i in 1:length(gup)){
    load(gup[i])
    temp <- rbind(temp,dat[,c(7:10)])
  }
  temp <- na.omit(temp)
  temp <- temp[sample(1:dim(temp)[1],100000),]

  pheVal[[vv]] <- temp

  print(vv)
}


###
  setwd('~/data_paper/figure/')
  png(filename=paste0('1to1_hls_ext_',j,'.png'),width=12,height=7.52,unit='in',res=600)

  par(mfrow=c(2,3),oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
  spec <- rev(brewer.pal(11,'Spectral'))
  # spec <- viridis(100)
  mycolRamp = colorRampPalette(c('White',spec))
  for(i in 1:6){
    x1 <- pheVal[[i]][,j+1]
    y1 <- pheVal[[i]][,  j]

    smoothScatter(x1,y1,xlim=c(-100,500),ylim=c(-100,500),
                  colramp=mycolRamp,nbin=850,nrpoints=0,
                  transformation=function(x)x^.6,axe=F,
                  xlab='PLSP (day of year)',ylab='MSLSP (day of year)',
                  cex.lab=1.8)
    axis(1,at=seq(-100,500,100),cex.axis=1.6)
    axis(2,at=seq(-100,500,100),cex.axis=1.6)
    abline(0,1,lty=5)


    reg <- formatC(round(cor(x1,y1,use='na.or.complete'),2), digits=2,format="fg", flag="#")
    rmse <- round(sqrt(mean((x1-y1)^2,na.rm=T)))
    bias <- round(mean(x1-y1,na.rm=T),1)

    text(250,50,expression(paste(italic(r),' =',sep='')),cex=1.7,pos=4)
    text(310,50,reg,cex=1.7,pos=4)
    text(250,5,paste('RMSE = ',rmse,sep=''),cex=1.7,pos=4)
    text(250,-40,paste('Bias = ',bias,sep=''),cex=1.7,pos=4)
    text(250,-85,'n = 100000',cex=1.7,pos=4)

    text(-100,460,vari[i],cex=1.7,pos=4,font=2)
  }

  dev.off()




