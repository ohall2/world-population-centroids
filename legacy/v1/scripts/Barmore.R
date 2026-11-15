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
outFolder <- '/data/centroid/out/Barmore'  # Output folder
dmsp_folder <-  '/data/centroid/raster/' #folder with light dmsp tif
regions_folder <- '/data/centroid/region' #folder with adm region shpes


#CONSTANTS
gcs <- '+proj=longlat +a=6371000 +b=6371000 +no_defs' # EPSG 4035 proj4string

### FUNCTIONS



####RASTER FILES
setwd(wd)
#FUNCTIONS
getCenter <- function(coords, weights)
{
  lon.center <- weighted.mean(coords[[1]],weights)
  lat.center <- weighted.mean(coords[[2]],weights)  
  return(data.frame(lon.center,lat.center))
}

createSpatial <- function(center, proj)
{
  return(SpatialPoints(center,proj4string=CRS(proj)))
}

CRS_AzimuthalEquidistant <- function(lon,lat)
{
  CRS <- paste('+proj=aeqd +lat_0=',
               lat,
               ' +lon_0=',
               lon,
               ' +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs',
               sep='')
  return(CRS)
}

transformPoint2center <- function(points){
  data <- cbind(data.frame(coordinates(points)),data.frame(points)[,1])
  data <- na.omit(data)
  center <- getCenter(as.data.frame(coordinates(data[,1:2])),data[,3])  
  return(center)
}




run_Barmore <- function(light,region){
  
  # run for each element in shapefile by using parallell processing, independent of order, 
  # larger shape area will take longer
  out.list <- foreach(id=1:length(region[,1]),.inorder=FALSE) %dopar% { # 1:length(region[,1])  
    #extract adm region
    area <- region[id,] 
    
    
    ##check that region is existing and completly within the light data
    ext <- extent(area)
    a <- intersect(ext,extent(light)) #check if this should be processed. only if within area chosen
    if(!is.null(a)){ 
      
      
      
      ## print progress, will not show in GUI mode (e.g. Rstudio)
      # print(paste(Sys.time(),'Processing:',id,'of',length(region[,1])))
      
      # mask light to area
      light.crop <- mask(crop(light,extent(area),snap='out'),area)
      
      if(sum(light.crop[,],na.rm=T)!=0){ #make sure we have light values
        
        # convert to points
        light.points <- rasterToPoints(light.crop,spatial=TRUE) 
        
        # Find an estimation of a possible center based only on average lat/lon 
        center <- transformPoint2center(light.points)
        
        # Project the original data to an azimuthal equidistant projection centered on that point
        CRS <- CRS_AzimuthalEquidistant(center[1],center[2])
        
        
        light.points_proj <- sp::spTransform(light.points,CRS=crs(CRS))
        
        # Find a new center based on the weighted NS/EW average of these projected values
        center <- transformPoint2center(light.points_proj)
        
        # Iterate through this process until locate a center less than 1 meter away from previous center 
        while (abs(center[[1]]) >= 1 && abs(center[[2]]) >= 1){
          # Convert the azimuthal equidistant projection to lon/lat for use in reprojection
          center_lonlat <- as.data.frame(coordinates(spTransform(createSpatial(center,CRS),CRS(gcs))))
          # Reproject original data to azimuthal equidistant projection centered on new center estimation
          CRS <- CRS_AzimuthalEquidistant(center_lonlat[1],center_lonlat[2])
          light.points_proj <- spTransform(light.points,CRS=crs(CRS))
          # Find new center based on the weighted NS/EW average of new projected values
          center <-  transformPoint2center(light.points_proj)
        }
        # Reproject the final center to lon/lat
        center_final <- spTransform(createSpatial(center,CRS),CRS(gcs))
        
        out.df <- data.frame(area@data,coordinates(center_final),id=id)
        
      }else{ #no light values return NA
        out.df <- data.frame(area@data,lon.center=NA,lat.center=NA,id=id)
      }
    }else{
      out.df <- data.frame(area@data,lon.center=NA,lat.center=NA,id=id)
    }
    return(out.df)
  }
  
  return(out.list)
  
}


### light and shape files
fileList <- list.files(dmsp_folder, pattern="\\.tif$", full.names=FALSE)
regions <- c('gadm28_adm0_gcs','gadm28_adm1_gcs','gadm28_adm2_gcs','gadm28_adm3_gcs')


## convert year to file index
i <- 1+year-1992

# Loop for each adm level
for(x in 1:3){ #1:length(regions)){
    removeTmpFiles(h=0.1)
    gc(verbose = TRUE)
    
    outName <- paste(substr(fileList[i],1,9), '_', regions[x], '_Barmore', sep='')
    
    
    
    ##check if outfile allready exist, then skip to next.
    if(!file.exists(paste(outFolder,'/',outName,'.shp',sep=''))){
      
      
      ##LOAD DATA
      region <- readOGR(dsn=regions_folder,layer=regions[x])
      light <- raster(paste(dmsp_folder,fileList[i],sep=''))
      
      ## PRINT PROGRESS
      print(paste(Sys.time(),'filelist:',i,'of',length(fileList)))
      
      out.list <- run_Barmore(light,region)
     
      
      
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
      write.csv(nolight,paste(outFolder,'/',substr(fileList[i],1,9), '_', regions[x], '_Barmore_novalue.csv', sep=''))
      
    }
}
