#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# A High Spatial Resolution Land Surface Phenology Dataset for AmeriFlux and NEON Sites
#
# 05: A script for saving data layers into netCDF format
# 
# Author: Minkyu Moon; moon.minkyu@gmail.com
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Submitted to the job for a single chunk for a site with shell script:
# #!/bin/bash
# echo Submitting $1
# R --vanilla < ~/05_netCDF.R $1
#
# example submission command using default parameters:
# qsub -V -pe omp 2 -l h_rt=12:00:00 run_05.sh numSite
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

library(raster)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(maptools)
library(rasterVis)
library(ncdf4)

########################################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
# numSite <- 5

########################################
params <- fromJSON(file='~/PLSP_Parameters.json')
source(params$setup$rFunctions)

########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)


########################################
inBase <- paste0(params$setup$workDir,'Product_GeoTiff/')
outBase <- paste0(params$setup$workDir,'Product_netCDF/')

sites <- list.dirs(inBase,recursive=F,full.names=F)
site  <- sites[ss]
print(site)

########################################
# Get product layers info
productTable <- read.csv('~/PLSP_Layers.csv',header=T,stringsAsFactors = F)
# Get a base image to pull raster info from     
baseImage <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))

########################################
# Get extent, and then define pixel centers in the x and y direction
ext = extent(baseImage)
res = res(baseImage)[1]
x = seq(ext[1]+res/2,ext[2]-res/2, res)
y = seq(ext[3]+res/2,ext[4]-res/2, res)
#Define dimensions for netCDF file
dimx = ncdim_def(name = 'x', longname = 'x coordinate', units='m', vals = as.double(x))
dimy = ncdim_def(name = 'y', longname = 'y coordinate', units='m', vals = rev(as.double(y)))


########################################
# Do a loop for 5 years; 2017-2021
for(yy in 1:5){
  # Load files
  files <- list.files(paste0(inBase,site,'/',(2016+yy)),pattern=glob2rx(paste('*.tif',sep='')),full.names=T)
  
  outFold <- paste0(outBase, site, '/')
  if(!dir.exists(outFold)) {dir.create(outFold)}
  outFile <- paste0(outFold,'PLSP_',(2016+yy),'.nc')
        
  # Loop through all the layers, and create a variable for each 
  results      <-vector("list", dim(productTable)[1]+1) 
  results[[1]] <- ncvar_def("transverse_mercator","",list(),prec="char")
  
  for(i in 1:dim(productTable)[1]){
      lyr <- productTable[i,]   # Pull the info for this layer from the productTable 
          
      if (lyr$data_type == 'Int16') {precision <- "short"}  #All are int16, so this isn't necessary
          
      # Create the variable, add to the list. Define the short_name, units, fill_value, long_name, and precision from the product table 
      results[[i+1]] <- ncvar_def(lyr$short_name, lyr$units, list(dimx,dimy), lyr$fill_value, lyr$long_name, prec=precision, compression=2)  
  }
        
  # Now create the netCDF file with the defined variables
  if (file.exists(outFile)) {file.remove(outFile)}
  ncout <- nc_create(outFile,results,force_v4=T)
        
  # Now loop through the layers again, this time actually 
  # writing the image data to the file
  for (i in 1:dim(productTable)[1]) {
      lyr <- productTable[i,]
          
      # Open the data          
      mat <- matrix(values(raster(files[i],varname=lyr$short_name)),length(x),length(y))
      # Put the image into the file
      ncvar_put(ncout,results[[i+1]], mat) 
      
      # Fill in the attributes for the layer from the product table 
      ncatt_put(ncout,lyr$short_name,"scale",lyr$scale)
      ncatt_put(ncout,lyr$short_name,"offset",lyr$offset)
      ncatt_put(ncout,lyr$short_name,"data_type",lyr$data_type)
      ncatt_put(ncout,lyr$short_name,"valid_min",lyr$valid_min)
      ncatt_put(ncout,lyr$short_name,"valid_max",lyr$valid_max)
      
      print(paste((2016+yy),';',i))
  }
  
  ## Write the projection info for the transverse_mercator variable
  # Get projection in wkt format
  wkt <- showWKT(projection(baseImage), morphToESRI = FALSE)  
  # Need to pull the central meridian from the wkt
  spt <- unlist(strsplit(gsub(']','',wkt),','))
  central_meridian    <- as.numeric(spt[which(spt == "PARAMETER[\"central_meridian\"")+1])
  
  # Fill in the info.
  ncatt_put(ncout,"transverse_mercator","long_name","CRS definition")
  ncatt_put(ncout,"transverse_mercator","grid_mapping_name","transverse_mercator")
  ncatt_put(ncout,"transverse_mercator","longitude_of_central_meridian",central_meridian)
  ncatt_put(ncout,"transverse_mercator","false_easting",5e+05)
  ncatt_put(ncout,"transverse_mercator","false_northing",0)
  ncatt_put(ncout,"transverse_mercator","latitude_of_projection_origin",0)
  ncatt_put(ncout,"transverse_mercator","scale_factor_at_central_meridian",0.9996)
  ncatt_put(ncout,"transverse_mercator","longitude_of_prime_meridian",0)
  ncatt_put(ncout,"transverse_mercator","semi_major_axis",6378137)
  ncatt_put(ncout,"transverse_mercator","inverse_flattening",298.257223563)
  ncatt_put(ncout,"transverse_mercator","GeoTransform",paste(ext[1],res,0,ext[4],0,-res))
  ncatt_put(ncout,"transverse_mercator","spatial_ref",gsub("\\", "", wkt, fixed=TRUE))
  
  ## Define global attributes
  ncatt_put(ncout,0,"title","Land Surface Phenology from PlanetScope (PLSP)")
  ncatt_put(ncout,0,"product_version","v001")
  ncatt_put(ncout,0,"summary","A High Spatial Resolution Land Surface Phenology from PlancetScope for AmeriFlux and NEON sites")
  ncatt_put(ncout,0,"software_repository","git@github.com:BU-LCSC/PLSP.git")

  ncatt_put(ncout,0,"creator_name","Land Cover & Surface Climate Group, Department of Earth & Environment, Boston University")
  ncatt_put(ncout,0,"creator_type","group")
  ncatt_put(ncout,0,"creator_email","mkmoon@bu.edu")
  ncatt_put(ncout,0,"creator_institution","Boston University")
  
  ncatt_put(ncout,0,"contributor_name", "Minkyu Moon, Andrew R. Richardson, Thomas Milliman, Mark A. Friedl")
  ncatt_put(ncout,0,"contributor_role", "Developer, Co-Investigator, Collabolator, Principal Investigator")
  ncatt_put(ncout,0,"acknowledgement","This work was supported by NASA grant #80NSSC18K0334 and by NSF award #1702627.")
  
  #Put additional attributes on coordinates
  ncatt_put(ncout,"x","axis","projection_x_coordinate")
  ncatt_put(ncout,"y","axis","projection_y_coordinate")
  
  
  ## Close the file
  nc_close(ncout)
}

