rm(list=ls())
library(rgdal)
library(maptools)
library(raster)
library(sf)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpmisc)


## Load country area in km^2
land_area <- read.csv('land_area_world.csv',sep=';') %>% rename(ISO=Country.Code,landArea=X2015) %>% dplyr::select(ISO,landArea,Country.Name) 


## load aboufadel and barmore
af <- st_read(paste('out/Aboufadel_june2018/DMSP_2008_gadm28_adm0_gcsAboufadelA.shp',sep='')) %>% dplyr::select(ISO)
bm <- st_read(paste('20180620_BarmoreFinal/Centroids/DMSP_2008_gadm28_adm0_gcs_Barmore_final_flagged.shp',sep='')) %>% dplyr::select(ISO)

#make sure that they have the same length and order - return value if not.
print(length(which(af$ISO != bm$ISO)))


## calculate distance between centroids
dist <- af %>% mutate(dist = as.numeric(st_distance(.,bm,by_element = T))) %>% st_set_geometry(NULL) %>% left_join(.,land_area,by='ISO') %>% mutate(dist=log10(dist),areakm=log10(landArea)) 


## Make Plot
ggplot(dist,aes(x=areakm,y=dist,label=Country.Name )) + geom_point() + 
  xlab(bquote('Adm0 Area '*log[10]*'('*km^2*')')) + 
  ylab(bquote('Distance between centroids '*log[10]*'(m)')) + 
  stat_dens2d_filter(geom = "label_repel", keep.fraction = 0.08,size=2.5) + theme_bw() + stat_dens2d_filter(geom = "point", keep.fraction = 0.08,color='red')






