library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)

library(rjson)

library(RColorBrewer)
library(viridis)

########################################
args <- commandArgs()
print(args)

ss <- as.numeric(args[3]) # site
# ss <- 103

vari <- c('50PCGI','50PCGD')


########################################
Rast <- vector('list',10); ind <- 1
for(vv in 1:2){
  for(yy in 2017:2021){
    imgs <- list.files(params$setup$outDir,pattern=glob2rx(paste0('*',yy,'_',vari[vv],'.tif')),full.names=T,recursive=T)  
    
    Rast[[ind]] <- raster(imgs[ss])    
    
    print(ind)
    ind  <- ind + 1 
  }
}

siteStr <- unlist(strsplit(imgs[ss],'/'))[6]
print(siteStr)


siteFull <- read.csv('~/site_list.csv')
afcd <- siteFull$Site_amerifulx[which(siteFull$Site==siteStr)]
vgty <- siteFull$Veg[which(siteFull$Site==siteStr)]

print(afcd)
print(vgty)




########################################
# Get site name, image directory and coordinate
rast <- Rast[[1]]
cx <- (rast@extent@xmax+rast@extent@xmin)/2
cy <- (rast@extent@ymax+rast@extent@ymin)/2
  
# Point Shape file
utm_crs = rast@crs
site <- data.frame(1,cx,cy)
colnames(site) <- c('id','xcrd','ycrd')
xy   <- site[,c(2,3)]
bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=utm_crs)

geog_crs = CRS("+proj=longlat +datum=WGS84")
bb   <- spTransform(bb,geog_crs)
  
  
# Find associate PhenoCam site
pcMeta <- list.files('~/PhenoCam_V3/data_record_1',pattern='*_meta.json',full.names=T)

siteWPC <- c()
for(i in 1:length(pcMeta)){
  meta <- fromJSON(file=pcMeta[i])
  
  # Point Shape file
  geog_crs = CRS("+proj=longlat +datum=WGS84")
  site <- data.frame(1,meta$phenocam_site$lon,meta$phenocam_site$lat)
  colnames(site) <- c('id','xcrd','ycrd')
  xy   <- site[,c(2,3)]
  bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)

  utm_crs = rast@crs
  bb   <- spTransform(bb,utm_crs)
  
  intPt <- intersect(bb,rast)  

  if(length(intPt)>0){
    siteWPC <- c(siteWPC,meta$phenocam_site$sitename)
    
    siteMeta <- c(afcd,
                  meta$phenocam_site$sitename,
                  meta$phenocam_site$long_name,
                  meta$phenocam_site$lat,
                  meta$phenocam_site$lon,
                  meta$phenocam_site$elevation,
                  meta$phenocam_site$primary_veg_type,
                  meta$phenocam_site$secondary_veg_type,
                  meta$phenocam_site$date_start,
                  meta$phenocam_site$date_end,
                  meta$phenocam_site$landcover_igbp,
                  meta$phenocam_site$site_acknowledgements)
    save(siteMeta,file=paste0('~/data/comp_phenocam/site_pc_list/pc_',afcd,'_',meta$phenocam_site$sitename,'.rda'))
  } 
  
  if(i%%100==0) print(i)
}


print(length(siteWPC))


# Get phenometric
if(length(siteWPC)>0){


  for(i in 1:length(siteWPC)){
    fJson <- list.files('~/PhenoCam_V3/data_record_1',pattern=glob2rx(paste0(siteWPC[i],'_*json')),full.names=T)
    fMetr <- list.files('~/PhenoCam_V3/data_record_5',pattern=glob2rx(paste0(siteWPC[i],'_*3day_transition_dates.csv')),full.names=T)
    
    if(length(fMetr)>0){
      
      meta <- fromJSON(file=fJson)
      datP <- c()
      for(j in 1:length(fMetr)){
        tempP <- read.csv(fMetr[j],skip=16)
        datP  <- rbind(datP,tempP)
      }
      
      ## Get PhenoCam data
      dat1 <- sort(as.Date(datP[datP$direction=='rising'&datP$gcc_value=='gcc_90',8]))
      dat2 <- sort(as.Date(datP[datP$direction=='falling'&datP$gcc_value=='gcc_90',8]))
      
      
      # Get Planet data
      # Point Shape file
      geog_crs = CRS("+proj=longlat +datum=WGS84")
      site <- data.frame(1,meta$phenocam_site$lon,meta$phenocam_site$lat)
      colnames(site) <- c('id','xcrd','ycrd')
      xy   <- site[,c(2,3)]
      bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
      
      utm_crs = rast@crs
      bb   <- spTransform(bb,utm_crs)
      
      site <- data.frame(1,bb@bbox[1],(bb@bbox[2]+30))
      colnames(site) <- c('id','xcrd','ycrd')
      xy   <- site[,c(2,3)]
      bb1  <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=utm_crs)
      
      # Get values
      dat3 <- matrix(NA,5,1)
      dat4 <- matrix(NA,5,1)
      for(j in 1:5){
        rast <- Rast[[j]]
        dat3[j] <- as.Date(round(median(unlist(extract(rast,bb1,buffer=7.5)),na.rm=T)) + as.numeric(as.Date(paste0(2015+j,'-12-31'))),origin='1970-1-1')
      }
      for(j in 1:5){
        rast <- Rast[[j+5]]
        dat4[j] <- as.Date(round(median(unlist(extract(rast,bb1,buffer=7.5)),na.rm=T)) + as.numeric(as.Date(paste0(2015+j,'-12-31'))),origin='1970-1-1')
      }
      
      save(afcd,vgty,meta,
           dat1,dat2,dat3,dat4,
           file=paste0('~/data/comp_phenocam/ext/pc_',siteStr,'_',siteWPC[i],'.rda'))
  
    }

   
    print(i)
  }


}



########################################
#
#
########################################
# Load files
files <- list.files('~/data/comp_phenocam/ext',pattern=glob2rx('*.rda'),full.names=T)


datMet1 <- c();  datVal1 <- c(); datAF1 <- c()
datMet2 <- c();  datVal2 <- c(); datAF2 <- c()
for(i in 1:length(files)){
  load(files[i])

  dat1 <- as.numeric(dat1)
  dat2 <- as.numeric(dat2)

  for(yy in 1:4){
    if(length(which(abs(dat1-dat3[yy])<100))>0){
      temp    <- c(dat1[which.min(abs(dat1-dat3[yy]))],dat3[yy])
      datVal1  <- rbind(datVal1,temp)
      datMet1  <- rbind(datMet1,vgty)
      datAF1   <- rbind(datAF1,afcd)
    }
  }

  for(yy in 1:4){
    if(length(which(abs(dat2-dat4[yy])<100))>0){
      temp    <- c(dat2[which.min(abs(dat2-dat4[yy]))],dat4[yy])
      datVal2  <- rbind(datVal2,temp)
      datMet2  <- rbind(datMet2,vgty)
      datAF1   <- rbind(datAF1,afcd)
    }
  }
}
aa1 <- as.Date(datVal1[,1],origin='1970-1-1')
aa2 <- as.Date(datVal1[,2],origin='1970-1-1')
aa11 <- as.numeric(strftime(as.Date(datVal1[,1],origin='1970-1-1'), format = "%j"))
aa22 <- as.numeric(strftime(as.Date(datVal1[,2],origin='1970-1-1'), format = "%j"))
aa11[which(substr(aa1,1,4)!=substr(aa2,1,4))] <- aa11[which(substr(aa1,1,4)!=substr(aa2,1,4))]-365

bb1 <- as.Date(datVal2[,1],origin='1970-1-1')
bb2 <- as.Date(datVal2[,2],origin='1970-1-1')
bb11 <- as.numeric(strftime(as.Date(datVal2[,1],origin='1970-1-1'), format = "%j"))
bb22 <- as.numeric(strftime(as.Date(datVal2[,2],origin='1970-1-1'), format = "%j"))
bb22[which(substr(bb1,1,4)!=substr(bb2,1,4))] <- bb22[which(substr(bb1,1,4)!=substr(bb2,1,4))]-365

mycol <- brewer.pal(7,'Set1')
mycol <- c(mycol[1:5],mycol[7])
mycol <- brewer.pal(12,'Paired')

vgt <- sort(unique(c(datMet1,datMet2)))
siteList <- sort(unique(c(datAF1,datAF2)))


vgt1 <- matrix(NA,length(datMet1),1)
vgt2 <- matrix(NA,length(datMet2),1)
for(i in 1:length(vgt1)){
  if(datMet1[i]=='CRO'){       vgt1[i]=1
  }else if(datMet1[i]=='DBF'){ vgt1[i]=2
  }else if(datMet1[i]=='MF' ){ vgt1[i]=3
  }else if(datMet1[i]=='ENF'){ vgt1[i]=4
  }else if(datMet1[i]=='GRA'){ vgt1[i]=6
  }else if(datMet1[i]=='CSH'){ vgt1[i]=5
  }else if(datMet1[i]=='OSH'){ vgt1[i]=7
  }else if(datMet1[i]=='SAV'){ vgt1[i]=8
  }else if(datMet1[i]=='WSA'){ vgt1[i]=9
  }else if(datMet1[i]=='WET'){ vgt1[i]=10
  }else if(datMet1[i]=='EBF'){ vgt1[i]=11
  }else if(datMet1[i]=='CVM'){ vgt1[i]=12}
}
for(i in 1:length(vgt2)){
  if(datMet2[i]=='CRO'){       vgt2[i]=1
  }else if(datMet2[i]=='DBF'){ vgt2[i]=2
  }else if(datMet2[i]=='MF' ){ vgt2[i]=3
  }else if(datMet2[i]=='ENF'){ vgt2[i]=4
  }else if(datMet2[i]=='GRA'){ vgt2[i]=6
  }else if(datMet2[i]=='CSH'){ vgt2[i]=5
  }else if(datMet2[i]=='OSH'){ vgt2[i]=7
  }else if(datMet2[i]=='SAV'){ vgt2[i]=8
  }else if(datMet2[i]=='WSA'){ vgt2[i]=9
  }else if(datMet2[i]=='WET'){ vgt2[i]=10
  }else if(datMet2[i]=='EBF'){ vgt2[i]=11
  }else if(datMet2[i]=='CVM'){ vgt2[i]=12}
}



###
setwd('~/figure/')
png(filename='1to1_pc_ext.png',width=12,height=6,unit='in',res=600)

par(mfrow=c(1,2),oma=c(1,1,2,2),mar=c(4,4,1,1),mgp=c(2.7,0.8,0))
plot(aa22,aa11,xlim=c(-50,400),ylim=c(-50,400),pch=21,bg=mycol[vgt1],axe=F,
     xlab='PLSP (day of year)',ylab='PhenoCam (day of year)',
     cex=1.3,cex.lab=1.7)
abline(0,1,lty=5); box()
axis(1,at=seq(-100,500,100),cex.axis=1.5)
axis(2,at=seq(-100,500,100),cex.axis=1.5)

x1 <- aa22; y1 <- aa11
reg <- formatC(round(cor(x1,y1,use='na.or.complete'),2), digits=2,format="fg", flag="#")
rmse <- round(sqrt(mean((x1-y1)^2,na.rm=T)))
bias <- round(mean(x1-y1,na.rm=T),1)

text(160, 50,expression(paste(italic(r),' =',sep='')),cex=1.5,pos=4)
text(200, 50,reg,cex=1.5,pos=4)
text(160, 20,paste('RMSE = ',rmse,sep=''),cex=1.5,pos=4)
text(160,-10,paste('Bias = ',bias,sep=''),cex=1.5,pos=4)
text(160,-40,paste('n = ',length(aa22),sep=''),cex=1.5,pos=4)
text(-50,380,'50PCGI',cex=1.7,pos=4,font=2)

legend('bottomright',c('CRO','DBF','MF','ENF','GRA','CSH','OSH','SAV','WSA','WET','EBF','CVM'),
       pch=21,bty='n',pt.bg=mycol[c(1:4,6,5,7:12)],cex=1.3,pt.cex=1.3)


plot(bb11,bb22,xlim=c(-50,400),ylim=c(-50,400),pch=21,bg=mycol[vgt2],axe=F,
     xlab='PLSP (day of year)',ylab='PhenoCam (day of year)',
     cex=1.3,cex.lab=1.7)
abline(0,1,lty=5); box()
axis(1,at=seq(-100,500,100),cex.axis=1.5)
axis(2,at=seq(-100,500,100),cex.axis=1.5)

x1 <- bb11; y1 <- bb22
reg <- formatC(round(cor(x1,y1,use='na.or.complete'),2), digits=2,format="fg", flag="#")
rmse <- round(sqrt(mean((x1-y1)^2,na.rm=T)))
bias <- round(mean(x1-y1,na.rm=T),1)

text(230, 50,expression(paste(italic(r),' =',sep='')),cex=1.5,pos=4)
text(270, 50,reg,cex=1.5,pos=4)
text(230, 20,paste('RMSE = ',rmse,sep=''),cex=1.5,pos=4)
text(230,-10,paste('Bias = ',bias,sep=''),cex=1.5,pos=4)
text(230,-40,paste('n = ',length(bb22),sep=''),cex=1.5,pos=4)
text(-50,380,'50PCGD',cex=1.7,pos=4,font=2)

dev.off()

