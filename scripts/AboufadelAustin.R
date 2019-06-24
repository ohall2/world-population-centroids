library(rgdal)
library(maptools)
library(raster)
library(foreach)
library(doMC)
cores <- 48
registerDoMC(cores)

##set year to run script for
year <- 1991 #

#FOLDER SETTINGS
wd <- '/data/centroid/'  # Working directory
outFolder <- '/data/centroid/out/Aboufadel_june2018'  # Output folder
dmsp_folder <-  '/data/centroid/raster/' #folder with light dmsp tif
regions_folder <- '/data/centroid/region/New' #folder with adm region shpes



#CONSTANTS
gcs <- '+proj=longlat +a=6371000 +b=6371000 +no_defs' # EPSG 4035 proj4string

#FUNCTIONS
rad2deg <- function(rad) {(rad*180)/(pi)}
deg2rad <- function(deg) {(deg*pi)/(180)}

####

# Set the working directory
setwd(wd)


####RASTER FILES

fileList <- list.files(dmsp_folder, pattern="\\.tif$", full.names=FALSE)
regions <- c('gadm28_adm0_gcs','gadm28_adm1_gcs','gadm28_adm2_gcs')

## convert year to file index
i <- 1+year-1992

# Loop for each adm level
for (x in 1:length(regions)){ #
  
  region <- readOGR(dsn=regions_folder,layer=regions[x]) #load shape
  
  light <- raster(paste(dmsp_folder,fileList[i],sep='')) #load night time lights data
  
  print(paste(Sys.time(),'---',regions[x],'filelist:',i,'of',length(fileList))) #printing to console
  
  outName <- paste(substr(fileList[i],1,9), '_', regions[x], 'AboufadelA', sep='') # creating a unique output name
  
  ##check if outfile allready exist, then skip to next.
  if(!file.exists(paste(outFolder,'/',outName,'.shp',sep=''))){  
    # run for each element in shapefile by using parallell processing, independent of order, 
    # larger shape area will take longer
    out.list <- foreach(id=1:length(region[,1]),.inorder=FALSE) %dopar% {
      #extract adm region
      area <- region[id,] 
      ext <- extent(area)
      a <- intersect(ext,extent(light)) #check if this should be processed. only if within area chosen
      if(!is.null(a)){
        ## print progress, will not show in GUI mode (e.g. Rstudio)
        ##print(paste(Sys.time(),'Processing:',id,'of',length(region[,1])))
        # mask light to area
        light.crop <- mask(crop(light,extent(area),snap='out'),area)
        #make sure we have light values
        if(sum(light.crop[,],na.rm=T)!=0){ 
          # convert to points
          light.points <- rasterToPoints(light.crop,spatial=TRUE) 
          # Convert to dataframe
          light.data <- as.data.frame(light.points)
          
          colnr <- which(colnames(light.data) != 'x' & colnames(light.data) != 'y'  )
          
          
          # Calculate weighted x,y,z coordinates
          light.data$wixi = light.data[,colnr]*(cos(deg2rad(light.data$y))*cos(deg2rad(light.data$x)))
          light.data$wiyi = light.data[,colnr]*(cos(deg2rad(light.data$y))*sin(deg2rad(light.data$x)))
          light.data$wizi = light.data[,colnr]*(sin(deg2rad(light.data$y)))
          
          
          
          #Create output data frame
          out <- data.frame(sum(light.data[,colnr]),sum(light.data$wixi),sum(light.data$wiyi),sum(light.data$wizi))
          colnames(out) <- c("wi","wixi","wiyi","wizi")
          # Calculate center x,y,z coordindates
          out$x = out$wixi/out$wi
          out$y = out$wiyi/out$wi
          out$z = out$wizi/out$wi
          out$a = sqrt(out$x^2 + out$y^2 + out$z^2)
          out$xa = out$x/out$a
          out$ya = out$y/out$a
          
          # Convert to degrees
          #out$atan2 = atan2(out$y, out$x)
          out$Latitude = rad2deg(asin(out$z/out$a))
          out$Longitude = rad2deg(atan2(out$ya,out$xa)) #rad2deg(atan((out$y/out$a)/(out$x/out$a)))
          
          out.df <- data.frame(area@data,lon.center=out$Longitude,lat.center=out$Latitude)
        } else{ #no light values return NA
          out.df <- data.frame(area@data,lon.center=NA,lat.center=NA) 
        }
      }else{
        out.df <- data.frame(area@data,lon.center=NA,lat.center=NA)
      }
      return(out.df)
      
    }
    
    # test that length of output is same as number of regions
    if(length(out.list)!=length(region)){
      print('STOPPED RUNNING DUE TO LENGTH ISSUE')
      break
      
    } 
    
    # create data frame from list and create a SpatialPointsDataFrame object and output a shapefile
    out <- do.call(rbind.data.frame, out.list)
    nolight <- out[is.na(out$lat.center),]
    out <- out[!is.na(out$lat.center),]
    
    outshp <- SpatialPointsDataFrame(out[c('lon.center','lat.center')], out, proj4string = CRS(gcs))
    
    #write output to files both shape centroids and no value flag data table.
    writeOGR(outshp, dsn=outFolder, layer=outName, driver='ESRI Shapefile') 
    write.csv(nolight,paste(outFolder,'/',substr(fileList[i],1,9), '_', regions[x], '_AboufadelA_novalue.csv', sep=''))
    
    
    
  }
}

