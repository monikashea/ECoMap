# set your own working directory
setwd("U:/Shea/Research/GIS/Projects/TZ_Boundary/github")

library(indicspecies)
library(maptools)
library(spatstat)
library(sp)
library(rgdal)
library(ggplot2)
library(ggspatial)
library(RColorBrewer)
library(sf)
library(tidyr)
library(raster)
library(dplyr)


################### prepare data #######################
## import species data
all_trees_pred <- read.csv(file="all_trees_pred_12Aug2019.csv",header=TRUE)

# make spatial
coords <- cbind(all_trees_pred$X,all_trees_pred$Y)
ptree_diff_2018 <- SpatialPointsDataFrame(coords,all_trees_pred)
proj4string(ptree_diff_2018) <- CRS("+init=epsg:3070") # add projection NAD83/Wisconsin Transferse Mercator

## import 6 mi grid
Grid_6mi <- readOGR(dsn="./6mi_Grid",layer="6mi_Grid")
Grid_6mi <- spTransform(Grid_6mi,CRS=CRS("+init=epsg:3070"))

## import 2 km grid
Grid_2km <- readOGR(dsn="./Supplement/2km_Grid",layer="2km_Grid")
Grid_2km <- spTransform(Grid_2km,CRS=CRS("+init=epsg:3070"))

## import 4 km grid
Grid_4km <- readOGR(dsn="./Supplement/4km_Grid",layer="4km_Grid")
Grid_4km <- spTransform(Grid_4km,CRS=CRS("+init=epsg:3070"))

## import 12 mi grid
Grid_12mi <- readOGR(dsn="./Supplement/12mi_Grid",layer="12mi_Grid")
Grid_12mi <- spTransform(Grid_12mi,CRS=CRS("+init=epsg:3070"))

## import 18 mi grid
Grid_18mi <- readOGR(dsn="./Supplement/18mi_Grid",layer="18mi_Grid")
Grid_18mi <- spTransform(Grid_18mi,CRS=CRS("+init=epsg:3070"))

## bring in WI map
wisc<-readOGR(dsn="./WI_HARN_mask",layer="WI_HARN_mask")
wisc_HARN <- spTransform(wisc,CRS=CRS("+init=epsg:3070"))


### determine which grid cell each tree falls in
# 6 mi
pg_overlay_6mi <- over(ptree_diff_2018,Grid_6mi)
pg_overlay_trees_6mi <- cbind(pg_overlay_6mi,ptree_diff_2018@data)
pg_overlay_trees_6mi <- pg_overlay_trees_6mi %>% rename(Grid_6mi=Grid_ID)

# 2 km
pg_overlay_2km <- over(ptree_diff_2018,Grid_2km)
pg_overlay_trees_2km <- cbind(pg_overlay_2km,ptree_diff_2018@data)
pg_overlay_trees_2km <- pg_overlay_trees_2km %>% rename(Grid_2km=Grid_ID)

# 4 km
pg_overlay_4km <- over(ptree_diff_2018,Grid_4km)
pg_overlay_trees_4km <- cbind(pg_overlay_4km,ptree_diff_2018@data)
pg_overlay_trees_4km <- pg_overlay_trees_4km %>% rename(Grid_4km=Grid_ID)

# 12 mi
pg_overlay_12mi <- over(ptree_diff_2018,Grid_12mi)
pg_overlay_trees_12mi <- cbind(pg_overlay_12mi,ptree_diff_2018@data)
pg_overlay_trees_12mi <- pg_overlay_trees_12mi %>% rename(Grid_12mi=Grid_ID)

# 18 mi
pg_overlay_18mi <- over(ptree_diff_2018,Grid_18mi)
pg_overlay_trees_18mi <- cbind(pg_overlay_18mi,ptree_diff_2018@data)
pg_overlay_trees_18mi <- pg_overlay_trees_18mi %>% rename(Grid_18mi=Grid_ID)


### identify ecoregion for each grid cell
# bring in zones spatial
Zones <- readOGR("./212_222_boundary","212_222_merge_extended")
Zones <- spTransform(Zones,CRS=CRS("+init=epsg:3070"))

## 6mi
# convert grid to raster
r_6mi<-raster(ncol=51,nrow=54)
extent(r_6mi)<-extent(Grid_6mi)
GridRaster_6mi<-rasterize(Grid_6mi,r_6mi,'Grid_ID')
crs(GridRaster_6mi) <- "+init=epsg:3070"

# FS Normal
GridZonesFS_6mi<-raster::extract(GridRaster_6mi,Zones,df=TRUE)
GridZonesFS_6mi$ZoneFS<-ifelse (GridZonesFS_6mi$ID<2,"South","North")
GridZonesFS_6mi<-GridZonesFS_6mi %>% rename(Grid_6mi=layer) %>% select(-ID)

## 2km
# convert grid to raster
r_2km<-raster(ncol=242,nrow=256)
extent(r_2km)<-extent(Grid_2km)
GridRaster_2km<-rasterize(Grid_2km,r_2km,'Grid_ID')
crs(GridRaster_2km) <- "+init=epsg:3070"

# FS Normal
GridZonesFS_2km<-raster::extract(GridRaster_2km,Zones,df=TRUE)
GridZonesFS_2km$ZoneFS<-ifelse (GridZonesFS_2km$ID<2,"South","North")
GridZonesFS_2km<-GridZonesFS_2km %>% rename(Grid_2km=layer) %>% select(-ID)

## 4km
# convert grid to raster
r_4km<-raster(ncol=121,nrow=128)
extent(r_4km)<-extent(Grid_4km)
GridRaster_4km<-rasterize(Grid_4km,r_4km,'Grid_ID')
crs(GridRaster_4km) <- "+init=epsg:3070"

# FS Normal
GridZonesFS_4km<-raster::extract(GridRaster_4km,Zones,df=TRUE)
GridZonesFS_4km$ZoneFS<-ifelse (GridZonesFS_4km$ID<2,"South","North")
GridZonesFS_4km<-GridZonesFS_4km %>% rename(Grid_4km=layer) %>% select(-ID)

## 12mi
# convert grid to raster
r_12mi<-raster(ncol=26,nrow=27)
extent(r_12mi)<-extent(Grid_12mi)
GridRaster_12mi<-rasterize(Grid_12mi,r_12mi,'Grid_ID')
crs(GridRaster_12mi) <- "+init=epsg:3070"

# FS Normal
GridZonesFS_12mi<-raster::extract(GridRaster_12mi,Zones,df=TRUE)
GridZonesFS_12mi$ZoneFS<-ifelse (GridZonesFS_12mi$ID<2,"South","North")
GridZonesFS_12mi<-GridZonesFS_12mi %>% rename(Grid_12mi=layer) %>% select(-ID)

## 18mi
# convert grid to raster
r_18mi<-raster(ncol=17,nrow=18)
extent(r_18mi)<-extent(Grid_18mi)
GridRaster_18mi<-rasterize(Grid_18mi,r_18mi,'Grid_ID')
crs(GridRaster_18mi) <- "+init=epsg:3070"

# FS Normal
GridZonesFS_18mi<-raster::extract(GridRaster_18mi,Zones,df=TRUE)
GridZonesFS_18mi$ZoneFS<-ifelse (GridZonesFS_18mi$ID<2,"South","North")
GridZonesFS_18mi<-GridZonesFS_18mi %>% rename(Grid_18mi=layer) %>% select(-ID)

### identify cells than fall entirely inside state boundary (excluding edge cells)
Grid_WI_6mi <- raster::extract(GridRaster_6mi,wisc_HARN,df=TRUE,weights=TRUE)
Grid_WI_6mi<-Grid_WI_6mi %>% rename(Grid_6mi=layer) %>% filter(weight==max(weight))

Grid_WI_2km <- raster::extract(GridRaster_2km,wisc_HARN,df=TRUE,weights=TRUE)
Grid_WI_2km<-Grid_WI_2km %>% rename(Grid_2km=layer) %>% filter(weight==max(weight))

Grid_WI_4km <- raster::extract(GridRaster_4km,wisc_HARN,df=TRUE,weights=TRUE)
Grid_WI_4km<-Grid_WI_4km %>% rename(Grid_4km=layer) %>% filter(weight==max(weight))

Grid_WI_12mi <- raster::extract(GridRaster_12mi,wisc_HARN,df=TRUE,weights=TRUE)
Grid_WI_12mi<-Grid_WI_12mi %>% rename(Grid_12mi=layer) %>% filter(weight==max(weight))

Grid_WI_18mi <- raster::extract(GridRaster_18mi,wisc_HARN,df=TRUE,weights=TRUE)
Grid_WI_18mi<-Grid_WI_18mi %>% rename(Grid_18mi=layer) %>% filter(weight==max(weight))

### join tree data to zone data
## 6mi
Trees_6mi <- left_join(pg_overlay_trees_6mi,GridZonesFS_6mi,by="Grid_6mi")

# remove NA, UK, RR (NA, unknown species, and very rare species with <10 points)
Trees_6mi<-Trees_6mi %>% filter(!SP_new %in% c(' ','UK','RR')) %>% filter(!is.na(SP_new))
Trees_6mi$SP_new<-as.factor(Trees_6mi$SP_new)
Trees_6mi$SP_new<-droplevels(Trees_6mi$SP_new)

# convert to species counts for each cell, remove edge cells, and add ecoregion
Trees_count_6mi <- Trees_6mi %>% group_by(SP_new,Grid_6mi) %>% summarise(count=n()) %>% spread(SP_new,count) %>% 
  semi_join(Grid_WI_6mi,by="Grid_6mi") %>% left_join(GridZonesFS_6mi,by="Grid_6mi")
Trees_count_6mi[is.na(Trees_count_6mi)]<-0 # change NA to 0

# separate matrices for each ecoregion
# north
Trees_count_north_6mi<-filter(Trees_count_6mi,ZoneFS=="North")
row.names(Trees_count_north_6mi)<-Trees_count_north_6mi$Grid_6mi
Trees_count_north_6mi<-select(Trees_count_north_6mi,-c(Grid_6mi,ZoneFS))

# south
Trees_count_south_6mi<-filter(Trees_count_6mi,ZoneFS=="South")
row.names(Trees_count_south_6mi)<-Trees_count_south_6mi$Grid_6mi
Trees_count_south_6mi<-select(Trees_count_south_6mi,-c(Grid_6mi,ZoneFS))

## 2km
Trees_2km <- left_join(pg_overlay_trees_2km,GridZonesFS_2km,by="Grid_2km")

# remove NA, UK, RR (NA, unknown species, and very rare species with <10 points)
Trees_2km<-Trees_2km %>% filter(!SP_new %in% c(' ','UK','RR')) %>% filter(!is.na(SP_new))
Trees_2km$SP_new<-as.factor(Trees_2km$SP_new)
Trees_2km$SP_new<-droplevels(Trees_2km$SP_new)

# convert to species counts for each cell, remove edge cells, and add ecoregion
Trees_count_2km <- Trees_2km %>% group_by(SP_new,Grid_2km) %>% summarise(count=n()) %>% spread(SP_new,count) %>% 
  semi_join(Grid_WI_2km,by="Grid_2km") %>% left_join(GridZonesFS_2km,by="Grid_2km")
Trees_count_2km[is.na(Trees_count_2km)]<-0 # change NA to 0

# separate matrices for each ecoregion
# north
Trees_count_north_2km<-filter(Trees_count_2km,ZoneFS=="North")
row.names(Trees_count_north_2km)<-Trees_count_north_2km$Grid_2km
Trees_count_north_2km<-select(Trees_count_north_2km,-c(Grid_2km,ZoneFS))

# south
Trees_count_south_2km<-filter(Trees_count_2km,ZoneFS=="South")
row.names(Trees_count_south_2km)<-Trees_count_south_2km$Grid_2km
Trees_count_south_2km<-select(Trees_count_south_2km,-c(Grid_2km,ZoneFS))

## 4km
Trees_4km <- left_join(pg_overlay_trees_4km,GridZonesFS_4km,by="Grid_4km")

# remove NA, UK, RR (NA, unknown species, and very rare species with <10 points)
Trees_4km<-Trees_4km %>% filter(!SP_new %in% c(' ','UK','RR')) %>% filter(!is.na(SP_new))
Trees_4km$SP_new<-as.factor(Trees_4km$SP_new)
Trees_4km$SP_new<-droplevels(Trees_4km$SP_new)

# convert to species counts for each cell, remove edge cells, and add ecoregion
Trees_count_4km <- Trees_4km %>% group_by(SP_new,Grid_4km) %>% summarise(count=n()) %>% spread(SP_new,count) %>% 
  semi_join(Grid_WI_4km,by="Grid_4km") %>% left_join(GridZonesFS_4km,by="Grid_4km")
Trees_count_4km[is.na(Trees_count_4km)]<-0 # change NA to 0

# separate matrices for each ecoregion
# north
Trees_count_north_4km<-filter(Trees_count_4km,ZoneFS=="North")
row.names(Trees_count_north_4km)<-Trees_count_north_4km$Grid_4km
Trees_count_north_4km<-select(Trees_count_north_4km,-c(Grid_4km,ZoneFS))

# south
Trees_count_south_4km<-filter(Trees_count_4km,ZoneFS=="South")
row.names(Trees_count_south_4km)<-Trees_count_south_4km$Grid_4km
Trees_count_south_4km<-select(Trees_count_south_4km,-c(Grid_4km,ZoneFS))

## 12mi
Trees_12mi <- left_join(pg_overlay_trees_12mi,GridZonesFS_12mi,by="Grid_12mi")

# remove NA, UK, RR (NA, unknown species, and very rare species with <10 points)
Trees_12mi<-Trees_12mi %>% filter(!SP_new %in% c(' ','UK','RR')) %>% filter(!is.na(SP_new))
Trees_12mi$SP_new<-as.factor(Trees_12mi$SP_new)
Trees_12mi$SP_new<-droplevels(Trees_12mi$SP_new)

# convert to species counts for each cell, remove edge cells, and add ecoregion
Trees_count_12mi <- Trees_12mi %>% group_by(SP_new,Grid_12mi) %>% summarise(count=n()) %>% spread(SP_new,count) %>% 
  semi_join(Grid_WI_12mi,by="Grid_12mi") %>% left_join(GridZonesFS_12mi,by="Grid_12mi")
Trees_count_12mi[is.na(Trees_count_12mi)]<-0 # change NA to 0

# separate matrices for each ecoregion
# north
Trees_count_north_12mi<-filter(Trees_count_12mi,ZoneFS=="North")
row.names(Trees_count_north_12mi)<-Trees_count_north_12mi$Grid_12mi
Trees_count_north_12mi<-select(Trees_count_north_12mi,-c(Grid_12mi,ZoneFS))

# south
Trees_count_south_12mi<-filter(Trees_count_12mi,ZoneFS=="South")
row.names(Trees_count_south_12mi)<-Trees_count_south_12mi$Grid_12mi
Trees_count_south_12mi<-select(Trees_count_south_12mi,-c(Grid_12mi,ZoneFS))

## 18mi
Trees_18mi <- left_join(pg_overlay_trees_18mi,GridZonesFS_18mi,by="Grid_18mi")

# remove NA, UK, RR (NA, unknown species, and very rare species with <10 points)
Trees_18mi<-Trees_18mi %>% filter(!SP_new %in% c(' ','UK','RR')) %>% filter(!is.na(SP_new))
Trees_18mi$SP_new<-as.factor(Trees_18mi$SP_new)
Trees_18mi$SP_new<-droplevels(Trees_18mi$SP_new)

# convert to species counts for each cell, remove edge cells, and add ecoregion
Trees_count_18mi <- Trees_18mi %>% group_by(SP_new,Grid_18mi) %>% summarise(count=n()) %>% spread(SP_new,count) %>% 
  semi_join(Grid_WI_18mi,by="Grid_18mi") %>% left_join(GridZonesFS_18mi,by="Grid_18mi")
Trees_count_18mi[is.na(Trees_count_18mi)]<-0 # change NA to 0

# separate matrices for each ecoregion
# north
Trees_count_north_18mi<-filter(Trees_count_18mi,ZoneFS=="North")
row.names(Trees_count_north_18mi)<-Trees_count_north_18mi$Grid_18mi
Trees_count_north_18mi<-select(Trees_count_north_18mi,-c(Grid_18mi,ZoneFS))

# south
Trees_count_south_18mi<-filter(Trees_count_18mi,ZoneFS=="South")
row.names(Trees_count_south_18mi)<-Trees_count_south_18mi$Grid_18mi
Trees_count_south_18mi<-select(Trees_count_south_18mi,-c(Grid_18mi,ZoneFS))

##################### ISA sensitivity analysis ####################
### testing sensitivity of number of cells (this will be affected by cell size)
## ISA from multiple random draws
ISA_PLS<-function(north,south,nsample){
  # randomly select nsample cells from each ecoregion
  Trees_N_1<-sample_n(north,nsample,replace=FALSE)
  Trees_S_1<-sample_n(south,nsample,replace=FALSE)
  Trees_NS_1<-rbind(Trees_N_1,Trees_S_1)
  groups = c(rep(1,nsample),rep(2,nsample))
  Trees_NS_indval_1 = multipatt(Trees_NS_1,groups,control=how(nperm=9999))
  Trees_NS_indval_pvalue_1<-as.data.frame(Trees_NS_indval_1[8])
  Trees_NS_indval_pvalue_1$SP<-row.names(Trees_NS_indval_pvalue_1)
  species<-select(Trees_NS_indval_pvalue_1,SP)
  sign<-select(Trees_NS_indval_pvalue_1,sign.index)
  pvalue<-select(Trees_NS_indval_pvalue_1,sign.p.value)
  
  # formatting
  pvalue_data<-as.data.frame(pvalue)
  row.names(pvalue_data)<-unlist(species[1])
  
  # number of times each species was a significant indicator species
  pvalue_sig<-ifelse(pvalue_data<=0.05,1,0)
  pvalue_sig[is.na(pvalue_sig)]<-0
  pvalue_sig_sum<-rowSums(pvalue_sig)
  
  # number of times each species was NA indicator species (i.e. not enough of a difference between groups to test)
  pvalue_na<-ifelse(is.na(pvalue_data),1,0)
  pvalue_na_sum<-rowSums(pvalue_na)
  
  # formatting
  sign_data<-as.data.frame(sign)
  
  # northern indicator species
  sign_N<-ifelse(sign_data==1,1,0)
  sign_N_sig<-sign_N*pvalue_sig
  sign_N_sig[is.na(sign_N_sig)]<-0
  sign_N_sum<-rowSums(sign_N_sig)
  
  # southern indicator species
  sign_S<-ifelse(sign_data==2,1,0)
  sign_S_sig<-sign_S*pvalue_sig
  sign_S_sig[is.na(sign_S_sig)]<-0
  sign_S_sum<-rowSums(sign_S_sig)
  
  # both north and south
  sign_both<-ifelse(sign_data==3,1,0)
  
  sign_B_sig<-sign_both*pvalue_sig
  sign_both[is.na(sign_both)]<-0
  sign_B_sum<-rowSums(sign_both)
  
  # formatting
  results<-data.frame(pvalue_sig,sign_N_sig,sign_S_sig,sign_both,pvalue_na)
  colnames(results)<-c('pvalue_sig','sign_N_sig','sign_S_sig','sign_both','pvalue_na')
  results$SP_new <- row.names(results)
  results <- results %>% mutate(NSB=sign_N_sig+sign_S_sig+sign_both)
  results<-results
}


######## 6 mi #########
## n = 10 
niter<-50
IndSppList_10_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,10)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_10_6mi[[i]]<-IndSpp
}

IndSpp_10_6mi<-do.call("rbind",IndSppList_10_6mi)
NSpp_10_6mi<-IndSpp_10_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_10_6mi<-IndSpp_10_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_10_6mi<-IndSpp_10_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_10_6mi <- rbind(NSpp_10_6mi,SSpp_10_6mi,BSpp_10_6mi)
IndSpp_all_10_6mi <-IndSpp_all_10_6mi %>% mutate(n=10,Grid="6mi")

## n = 20 
niter<-50
IndSppList_20_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,20)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_20_6mi[[i]]<-IndSpp
}

IndSpp_20_6mi<-do.call("rbind",IndSppList_20_6mi)
NSpp_20_6mi<-IndSpp_20_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_20_6mi<-IndSpp_20_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_20_6mi<-IndSpp_20_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_20_6mi <- rbind(NSpp_20_6mi,SSpp_20_6mi,BSpp_20_6mi)
IndSpp_all_20_6mi <-IndSpp_all_20_6mi %>% mutate(n=20,Grid="6mi")

## n = 30 
niter<-50
IndSppList_30_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,30)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_30_6mi[[i]]<-IndSpp
}

IndSpp_30_6mi<-do.call("rbind",IndSppList_30_6mi)
NSpp_30_6mi<-IndSpp_30_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_30_6mi<-IndSpp_30_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_30_6mi<-IndSpp_30_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_30_6mi <- rbind(NSpp_30_6mi,SSpp_30_6mi,BSpp_30_6mi)
IndSpp_all_30_6mi <-IndSpp_all_30_6mi %>% mutate(n=30,Grid="6mi")

## n = 40 
niter<-50
IndSppList_40_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,40)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_40_6mi[[i]]<-IndSpp
}

IndSpp_40_6mi<-do.call("rbind",IndSppList_40_6mi)
NSpp_40_6mi<-IndSpp_40_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_40_6mi<-IndSpp_40_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_40_6mi<-IndSpp_40_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_40_6mi <- rbind(NSpp_40_6mi,SSpp_40_6mi,BSpp_40_6mi)
IndSpp_all_40_6mi <-IndSpp_all_40_6mi %>% mutate(n=40,Grid="6mi")

## n = 50 
niter<-50
IndSppList_50_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_50_6mi[[i]]<-IndSpp
}

IndSpp_50_6mi<-do.call("rbind",IndSppList_50_6mi)
NSpp_50_6mi<-IndSpp_50_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_50_6mi<-IndSpp_50_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_50_6mi<-IndSpp_50_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_50_6mi <- rbind(NSpp_50_6mi,SSpp_50_6mi,BSpp_50_6mi)
IndSpp_all_50_6mi <-IndSpp_all_50_6mi %>% mutate(n=50,Grid="6mi")
  
## n = 60 
niter<-50
IndSppList_60_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,60)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_60_6mi[[i]]<-IndSpp
}

IndSpp_60_6mi<-do.call("rbind",IndSppList_60_6mi)
NSpp_60_6mi<-IndSpp_60_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_60_6mi<-IndSpp_60_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_60_6mi<-IndSpp_60_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_60_6mi <- rbind(NSpp_60_6mi,SSpp_60_6mi,BSpp_60_6mi)
IndSpp_all_60_6mi <-IndSpp_all_60_6mi %>% mutate(n=60,Grid="6mi")

## n = 70 
niter<-50
IndSppList_70_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,70)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_70_6mi[[i]]<-IndSpp
}

IndSpp_70_6mi<-do.call("rbind",IndSppList_70_6mi)
NSpp_70_6mi<-IndSpp_70_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_70_6mi<-IndSpp_70_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_70_6mi<-IndSpp_70_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_70_6mi <- rbind(NSpp_70_6mi,SSpp_70_6mi,BSpp_70_6mi)
IndSpp_all_70_6mi <-IndSpp_all_70_6mi %>% mutate(n=70,Grid="6mi")

## n = 80 
niter<-50
IndSppList_80_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,80)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_80_6mi[[i]]<-IndSpp
}

IndSpp_80_6mi<-do.call("rbind",IndSppList_80_6mi)
NSpp_80_6mi<-IndSpp_80_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_80_6mi<-IndSpp_80_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_80_6mi<-IndSpp_80_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_80_6mi <- rbind(NSpp_80_6mi,SSpp_80_6mi,BSpp_80_6mi)
IndSpp_all_80_6mi <-IndSpp_all_80_6mi %>% mutate(n=80,Grid="6mi")

## n = 90 
niter<-50
IndSppList_90_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,90)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_90_6mi[[i]]<-IndSpp
}

IndSpp_90_6mi<-do.call("rbind",IndSppList_90_6mi)
NSpp_90_6mi<-IndSpp_90_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_90_6mi<-IndSpp_90_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_90_6mi<-IndSpp_90_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_90_6mi <- rbind(NSpp_90_6mi,SSpp_90_6mi,BSpp_90_6mi)
IndSpp_all_90_6mi <-IndSpp_all_90_6mi %>% mutate(n=90,Grid="6mi")

## n = 100 
niter<-50
IndSppList_100_6mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,100)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_100_6mi[[i]]<-IndSpp
}

IndSpp_100_6mi<-do.call("rbind",IndSppList_100_6mi)
NSpp_100_6mi<-IndSpp_100_6mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_100_6mi<-IndSpp_100_6mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_100_6mi<-IndSpp_100_6mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_100_6mi <- rbind(NSpp_100_6mi,SSpp_100_6mi,BSpp_100_6mi)
IndSpp_all_100_6mi <-IndSpp_all_100_6mi %>% mutate(n=100,Grid="6mi")

######## 2km #########
## n = 10 
niter<-50
IndSppList_10_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,10)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_10_2km[[i]]<-IndSpp
}

IndSpp_10_2km<-do.call("rbind",IndSppList_10_2km)
NSpp_10_2km<-IndSpp_10_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_10_2km<-IndSpp_10_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_10_2km<-IndSpp_10_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_10_2km <- rbind(NSpp_10_2km,SSpp_10_2km,BSpp_10_6mi)
IndSpp_all_10_2km <-IndSpp_all_10_2km %>% mutate(n=10,Grid="2km")

## n = 20 
niter<-50
IndSppList_20_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,20)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_20_2km[[i]]<-IndSpp
}

IndSpp_20_2km<-do.call("rbind",IndSppList_20_2km)
NSpp_20_2km<-IndSpp_20_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_20_2km<-IndSpp_20_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_20_2km<-IndSpp_20_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_20_2km <- rbind(NSpp_20_2km,SSpp_20_2km,BSpp_20_2km)
IndSpp_all_20_2km <-IndSpp_all_20_2km %>% mutate(n=20,Grid="2km")

## n = 30 
niter<-50
IndSppList_30_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,30)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_30_2km[[i]]<-IndSpp
}

IndSpp_30_2km<-do.call("rbind",IndSppList_30_2km)
NSpp_30_2km<-IndSpp_30_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_30_2km<-IndSpp_30_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_30_2km<-IndSpp_30_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_30_2km <- rbind(NSpp_30_2km,SSpp_30_2km,BSpp_30_2km)
IndSpp_all_30_2km <-IndSpp_all_30_2km %>% mutate(n=30,Grid="2km")

## n = 40 
niter<-50
IndSppList_40_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,40)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_40_2km[[i]]<-IndSpp
}

IndSpp_40_2km<-do.call("rbind",IndSppList_40_2km)
NSpp_40_2km<-IndSpp_40_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_40_2km<-IndSpp_40_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_40_2km<-IndSpp_40_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_40_2km <- rbind(NSpp_40_2km,SSpp_40_2km,BSpp_40_2km)
IndSpp_all_40_2km <-IndSpp_all_40_2km %>% mutate(n=40,Grid="2km")

## n = 50 
niter<-50
IndSppList_50_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_50_2km[[i]]<-IndSpp
}

IndSpp_50_2km<-do.call("rbind",IndSppList_50_2km)
NSpp_50_2km<-IndSpp_50_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_50_2km<-IndSpp_50_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_50_2km<-IndSpp_50_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_50_2km <- rbind(NSpp_50_2km,SSpp_50_2km,BSpp_50_2km)
IndSpp_all_50_2km <-IndSpp_all_50_2km %>% mutate(n=50,Grid="2km")

## n = 60 
niter<-50
IndSppList_60_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,60)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_60_2km[[i]]<-IndSpp
}

IndSpp_60_2km<-do.call("rbind",IndSppList_60_2km)
NSpp_60_2km<-IndSpp_60_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_60_2km<-IndSpp_60_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_60_2km<-IndSpp_60_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_60_2km <- rbind(NSpp_60_2km,SSpp_60_2km,BSpp_60_2km)
IndSpp_all_60_2km <-IndSpp_all_60_2km %>% mutate(n=60,Grid="2km")

## n = 70 
niter<-50
IndSppList_70_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,70)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_70_2km[[i]]<-IndSpp
}

IndSpp_70_2km<-do.call("rbind",IndSppList_70_2km)
NSpp_70_2km<-IndSpp_70_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_70_2km<-IndSpp_70_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_70_2km<-IndSpp_70_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_70_2km <- rbind(NSpp_70_2km,SSpp_70_2km,BSpp_70_2km)
IndSpp_all_70_2km <-IndSpp_all_70_2km %>% mutate(n=70,Grid="2km")

## n = 80 
niter<-50
IndSppList_80_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,80)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_80_2km[[i]]<-IndSpp
}

IndSpp_80_2km<-do.call("rbind",IndSppList_80_2km)
NSpp_80_2km<-IndSpp_80_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_80_2km<-IndSpp_80_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_80_2km<-IndSpp_80_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_80_2km <- rbind(NSpp_80_2km,SSpp_80_2km,BSpp_80_2km)
IndSpp_all_80_2km <-IndSpp_all_80_2km %>% mutate(n=80,Grid="2km")

## n = 90 
niter<-50
IndSppList_90_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,90)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_90_2km[[i]]<-IndSpp
}

IndSpp_90_2km<-do.call("rbind",IndSppList_90_2km)
NSpp_90_2km<-IndSpp_90_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_90_2km<-IndSpp_90_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_90_2km<-IndSpp_90_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_90_2km <- rbind(NSpp_90_2km,SSpp_90_2km,BSpp_90_2km)
IndSpp_all_90_2km <-IndSpp_all_90_2km %>% mutate(n=90,Grid="2km")

## n = 100 
niter<-50
IndSppList_100_2km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,100)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_100_2km[[i]]<-IndSpp
}

IndSpp_100_2km<-do.call("rbind",IndSppList_100_2km)
NSpp_100_2km<-IndSpp_100_2km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_100_2km<-IndSpp_100_2km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_100_2km<-IndSpp_100_2km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_100_2km <- rbind(NSpp_100_2km,SSpp_100_2km,BSpp_100_2km)
IndSpp_all_100_2km <-IndSpp_all_100_2km %>% mutate(n=100,Grid="2km")

######## 4km #########
## n = 10 
niter<-50
IndSppList_10_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,10)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_10_4km[[i]]<-IndSpp
}

IndSpp_10_4km<-do.call("rbind",IndSppList_10_4km)
NSpp_10_4km<-IndSpp_10_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_10_4km<-IndSpp_10_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_10_4km<-IndSpp_10_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_10_4km <- rbind(NSpp_10_4km,SSpp_10_4km,BSpp_10_4km)
IndSpp_all_10_4km <-IndSpp_all_10_4km %>% mutate(n=10,Grid="4km")

## n = 20 
niter<-50
IndSppList_20_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,20)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_20_4km[[i]]<-IndSpp
}

IndSpp_20_4km<-do.call("rbind",IndSppList_20_4km)
NSpp_20_4km<-IndSpp_20_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_20_4km<-IndSpp_20_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_20_4km<-IndSpp_20_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_20_4km <- rbind(NSpp_20_4km,SSpp_20_4km,BSpp_20_4km)
IndSpp_all_20_4km <-IndSpp_all_20_4km %>% mutate(n=20,Grid="4km")

## n = 30 
niter<-50
IndSppList_30_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,30)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_30_4km[[i]]<-IndSpp
}

IndSpp_30_4km<-do.call("rbind",IndSppList_30_4km)
NSpp_30_4km<-IndSpp_30_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_30_4km<-IndSpp_30_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_30_4km<-IndSpp_30_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_30_4km <- rbind(NSpp_30_4km,SSpp_30_4km,BSpp_30_4km)
IndSpp_all_30_4km <-IndSpp_all_30_4km %>% mutate(n=30,Grid="4km")

## n = 40 
niter<-50
IndSppList_40_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,40)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_40_4km[[i]]<-IndSpp
}

IndSpp_40_4km<-do.call("rbind",IndSppList_40_4km)
NSpp_40_4km<-IndSpp_40_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_40_4km<-IndSpp_40_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_40_4km<-IndSpp_40_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_40_4km <- rbind(NSpp_40_4km,SSpp_40_4km,BSpp_40_4km)
IndSpp_all_40_4km <-IndSpp_all_40_4km %>% mutate(n=40,Grid="4km")

## n = 50 
niter<-50
IndSppList_50_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_50_4km[[i]]<-IndSpp
}

IndSpp_50_4km<-do.call("rbind",IndSppList_50_4km)
NSpp_50_4km<-IndSpp_50_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_50_4km<-IndSpp_50_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_50_4km<-IndSpp_50_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_50_4km <- rbind(NSpp_50_4km,SSpp_50_4km,BSpp_50_4km)
IndSpp_all_50_4km <-IndSpp_all_50_4km %>% mutate(n=50,Grid="4km")

## n = 60 
niter<-50
IndSppList_60_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,60)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_60_4km[[i]]<-IndSpp
}

IndSpp_60_4km<-do.call("rbind",IndSppList_60_4km)
NSpp_60_4km<-IndSpp_60_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_60_4km<-IndSpp_60_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_60_4km<-IndSpp_60_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_60_4km <- rbind(NSpp_60_4km,SSpp_60_4km,BSpp_60_4km)
IndSpp_all_60_4km <-IndSpp_all_60_4km %>% mutate(n=60,Grid="4km")

## n = 70 
niter<-50
IndSppList_70_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,70)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_70_4km[[i]]<-IndSpp
}

IndSpp_70_4km<-do.call("rbind",IndSppList_70_4km)
NSpp_70_4km<-IndSpp_70_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_70_4km<-IndSpp_70_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_70_4km<-IndSpp_70_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_70_4km <- rbind(NSpp_70_4km,SSpp_70_4km,BSpp_70_4km)
IndSpp_all_70_4km <-IndSpp_all_70_4km %>% mutate(n=70,Grid="4km")

## n = 80 
niter<-50
IndSppList_80_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,80)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_80_4km[[i]]<-IndSpp
}

IndSpp_80_4km<-do.call("rbind",IndSppList_80_4km)
NSpp_80_4km<-IndSpp_80_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_80_4km<-IndSpp_80_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_80_4km<-IndSpp_80_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_80_4km <- rbind(NSpp_80_4km,SSpp_80_4km,BSpp_80_4km)
IndSpp_all_80_4km <-IndSpp_all_80_4km %>% mutate(n=80,Grid="4km")

## n = 90 
niter<-50
IndSppList_90_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,90)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_90_4km[[i]]<-IndSpp
}

IndSpp_90_4km<-do.call("rbind",IndSppList_90_4km)
NSpp_90_4km<-IndSpp_90_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_90_4km<-IndSpp_90_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_90_4km<-IndSpp_90_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_90_4km <- rbind(NSpp_90_4km,SSpp_90_4km,BSpp_90_4km)
IndSpp_all_90_4km <-IndSpp_all_90_4km %>% mutate(n=90,Grid="4km")

## n = 100 
niter<-50
IndSppList_100_4km<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,100)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_100_4km[[i]]<-IndSpp
}

IndSpp_100_4km<-do.call("rbind",IndSppList_100_4km)
NSpp_100_4km<-IndSpp_100_4km %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_100_4km<-IndSpp_100_4km %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_100_4km<-IndSpp_100_4km %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_100_4km <- rbind(NSpp_100_4km,SSpp_100_4km,BSpp_100_4km)
IndSpp_all_100_4km <-IndSpp_all_100_4km %>% mutate(n=100,Grid="4km")

######## 12mi #########
## n = 10 
niter<-50
IndSppList_10_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,10)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_10_12mi[[i]]<-IndSpp
}

IndSpp_10_12mi<-do.call("rbind",IndSppList_10_12mi)
NSpp_10_12mi<-IndSpp_10_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_10_12mi<-IndSpp_10_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_10_12mi<-IndSpp_10_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_10_12mi <- rbind(NSpp_10_12mi,SSpp_10_12mi,BSpp_10_12mi)
IndSpp_all_10_12mi <-IndSpp_all_10_12mi %>% mutate(n=10,Grid="12mi")

## n = 20 
niter<-50
IndSppList_20_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,20)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_20_12mi[[i]]<-IndSpp
}

IndSpp_20_12mi<-do.call("rbind",IndSppList_20_12mi)
NSpp_20_12mi<-IndSpp_20_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_20_12mi<-IndSpp_20_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_20_12mi<-IndSpp_20_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_20_12mi <- rbind(NSpp_20_12mi,SSpp_20_12mi,BSpp_20_12mi)
IndSpp_all_20_12mi <-IndSpp_all_20_12mi %>% mutate(n=20,Grid="12mi")

## n = 30 
niter<-50
IndSppList_30_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,30)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_30_12mi[[i]]<-IndSpp
}

IndSpp_30_12mi<-do.call("rbind",IndSppList_30_12mi)
NSpp_30_12mi<-IndSpp_30_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_30_12mi<-IndSpp_30_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_30_12mi<-IndSpp_30_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_30_12mi <- rbind(NSpp_30_12mi,SSpp_30_12mi,BSpp_30_12mi)
IndSpp_all_30_12mi <-IndSpp_all_30_12mi %>% mutate(n=30,Grid="12mi")

## n = 40 
niter<-50
IndSppList_40_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,40)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_40_12mi[[i]]<-IndSpp
}

IndSpp_40_12mi<-do.call("rbind",IndSppList_40_12mi)
NSpp_40_12mi<-IndSpp_40_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_40_12mi<-IndSpp_40_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_40_12mi<-IndSpp_40_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_40_12mi <- rbind(NSpp_40_12mi,SSpp_40_12mi,BSpp_40_12mi)
IndSpp_all_40_12mi <-IndSpp_all_40_12mi %>% mutate(n=40,Grid="12mi")

## n = 50 
niter<-50
IndSppList_50_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_50_12mi[[i]]<-IndSpp
}

IndSpp_50_12mi<-do.call("rbind",IndSppList_50_12mi)
NSpp_50_12mi<-IndSpp_50_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_50_12mi<-IndSpp_50_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_50_12mi<-IndSpp_50_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_50_12mi <- rbind(NSpp_50_12mi,SSpp_50_12mi,BSpp_50_12mi)
IndSpp_all_50_12mi <-IndSpp_all_50_12mi %>% mutate(n=50,Grid="12mi")

## n = 60 
niter<-50
IndSppList_60_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,60)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_60_12mi[[i]]<-IndSpp
}

IndSpp_60_12mi<-do.call("rbind",IndSppList_60_12mi)
NSpp_60_12mi<-IndSpp_60_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_60_12mi<-IndSpp_60_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_60_12mi<-IndSpp_60_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_60_12mi <- rbind(NSpp_60_12mi,SSpp_60_12mi,BSpp_60_12mi)
IndSpp_all_60_12mi <-IndSpp_all_60_12mi %>% mutate(n=60,Grid="12mi")

## n = 70 
niter<-50
IndSppList_70_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,70)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_70_12mi[[i]]<-IndSpp
}

IndSpp_70_12mi<-do.call("rbind",IndSppList_70_12mi)
NSpp_70_12mi<-IndSpp_70_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_70_12mi<-IndSpp_70_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_70_12mi<-IndSpp_70_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_70_12mi <- rbind(NSpp_70_12mi,SSpp_70_12mi,BSpp_70_12mi)
IndSpp_all_70_12mi <-IndSpp_all_70_12mi %>% mutate(n=70,Grid="12mi")

## n = 80 
niter<-50
IndSppList_80_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,80)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_80_12mi[[i]]<-IndSpp
}

IndSpp_80_12mi<-do.call("rbind",IndSppList_80_12mi)
NSpp_80_12mi<-IndSpp_80_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_80_12mi<-IndSpp_80_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_80_12mi<-IndSpp_80_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_80_12mi <- rbind(NSpp_80_12mi,SSpp_80_12mi,BSpp_80_12mi)
IndSpp_all_80_12mi <-IndSpp_all_80_12mi %>% mutate(n=80,Grid="12mi")

## n = 90 
niter<-50
IndSppList_90_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,90)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_90_12mi[[i]]<-IndSpp
}

IndSpp_90_12mi<-do.call("rbind",IndSppList_90_12mi)
NSpp_90_12mi<-IndSpp_90_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_90_12mi<-IndSpp_90_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_90_12mi<-IndSpp_90_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_90_12mi <- rbind(NSpp_90_12mi,SSpp_90_12mi,BSpp_90_12mi)
IndSpp_all_90_12mi <-IndSpp_all_90_12mi %>% mutate(n=90,Grid="12mi")

## n = 100 
niter<-50
IndSppList_100_12mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,100)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_100_12mi[[i]]<-IndSpp
}

IndSpp_100_12mi<-do.call("rbind",IndSppList_100_12mi)
NSpp_100_12mi<-IndSpp_100_12mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_100_12mi<-IndSpp_100_12mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_100_12mi<-IndSpp_100_12mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_100_12mi <- rbind(NSpp_100_12mi,SSpp_100_12mi,BSpp_100_12mi)
IndSpp_all_100_12mi <-IndSpp_all_100_12mi %>% mutate(n=100,Grid="12mi")

######## 18mi #########
## n = 10 
niter<-50
IndSppList_10_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,10)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_10_18mi[[i]]<-IndSpp
}

IndSpp_10_18mi<-do.call("rbind",IndSppList_10_18mi)
NSpp_10_18mi<-IndSpp_10_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_10_18mi<-IndSpp_10_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_10_18mi<-IndSpp_10_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_10_18mi <- rbind(NSpp_10_18mi,SSpp_10_18mi,BSpp_10_18mi)
IndSpp_all_10_18mi <-IndSpp_all_10_18mi %>% mutate(n=10,Grid="18mi")

## n = 20 
niter<-50
IndSppList_20_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,20)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_20_18mi[[i]]<-IndSpp
}

IndSpp_20_18mi<-do.call("rbind",IndSppList_20_18mi)
NSpp_20_18mi<-IndSpp_20_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_20_18mi<-IndSpp_20_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_20_18mi<-IndSpp_20_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_20_18mi <- rbind(NSpp_20_18mi,SSpp_20_18mi,BSpp_20_18mi)
IndSpp_all_20_18mi <-IndSpp_all_20_18mi %>% mutate(n=20,Grid="18mi")

## n = 30 
niter<-50
IndSppList_30_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,30)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_30_18mi[[i]]<-IndSpp
}

IndSpp_30_18mi<-do.call("rbind",IndSppList_30_18mi)
NSpp_30_18mi<-IndSpp_30_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_30_18mi<-IndSpp_30_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_30_18mi<-IndSpp_30_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_30_18mi <- rbind(NSpp_30_18mi,SSpp_30_18mi,BSpp_30_18mi)
IndSpp_all_30_18mi <-IndSpp_all_30_18mi %>% mutate(n=30,Grid="18mi")

## n = 40 
niter<-50
IndSppList_40_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,40)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_40_18mi[[i]]<-IndSpp
}

IndSpp_40_18mi<-do.call("rbind",IndSppList_40_18mi)
NSpp_40_18mi<-IndSpp_40_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_40_18mi<-IndSpp_40_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_40_18mi<-IndSpp_40_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_40_18mi <- rbind(NSpp_40_18mi,SSpp_40_18mi,BSpp_40_18mi)
IndSpp_all_40_18mi <-IndSpp_all_40_18mi %>% mutate(n=40,Grid="18mi")

## n = 50 
niter<-50
IndSppList_50_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_50_18mi[[i]]<-IndSpp
}

IndSpp_50_18mi<-do.call("rbind",IndSppList_50_18mi)
NSpp_50_18mi<-IndSpp_50_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_50_18mi<-IndSpp_50_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_50_18mi<-IndSpp_50_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_50_18mi <- rbind(NSpp_50_18mi,SSpp_50_18mi,BSpp_50_18mi)
IndSpp_all_50_18mi <-IndSpp_all_50_18mi %>% mutate(n=50,Grid="18mi")

## n = 60 
niter<-50
IndSppList_60_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,60)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_60_18mi[[i]]<-IndSpp
}

IndSpp_60_18mi<-do.call("rbind",IndSppList_60_18mi)
NSpp_60_18mi<-IndSpp_60_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_60_18mi<-IndSpp_60_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_60_18mi<-IndSpp_60_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_60_18mi <- rbind(NSpp_60_18mi,SSpp_60_18mi,BSpp_60_18mi)
IndSpp_all_60_18mi <-IndSpp_all_60_18mi %>% mutate(n=60,Grid="18mi")

## n = 70 
niter<-50
IndSppList_70_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,70)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_70_18mi[[i]]<-IndSpp
}

IndSpp_70_18mi<-do.call("rbind",IndSppList_70_18mi)
NSpp_70_18mi<-IndSpp_70_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_70_18mi<-IndSpp_70_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_70_18mi<-IndSpp_70_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_70_18mi <- rbind(NSpp_70_18mi,SSpp_70_18mi,BSpp_70_18mi)
IndSpp_all_70_18mi <-IndSpp_all_70_18mi %>% mutate(n=70,Grid="18mi")

## n = 80 
niter<-50
IndSppList_80_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,80)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_80_18mi[[i]]<-IndSpp
}

IndSpp_80_18mi<-do.call("rbind",IndSppList_80_18mi)
NSpp_80_18mi<-IndSpp_80_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_80_18mi<-IndSpp_80_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_80_18mi<-IndSpp_80_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_80_18mi <- rbind(NSpp_80_18mi,SSpp_80_18mi,BSpp_80_18mi)
IndSpp_all_80_18mi <-IndSpp_all_80_18mi %>% mutate(n=80,Grid="18mi")

## n = 90 
niter<-50
IndSppList_90_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,90)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_90_18mi[[i]]<-IndSpp
}

IndSpp_90_18mi<-do.call("rbind",IndSppList_90_18mi)
NSpp_90_18mi<-IndSpp_90_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_90_18mi<-IndSpp_90_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_90_18mi<-IndSpp_90_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_90_18mi <- rbind(NSpp_90_18mi,SSpp_90_18mi,BSpp_90_18mi)
IndSpp_all_90_18mi <-IndSpp_all_90_18mi %>% mutate(n=90,Grid="18mi")

## n = 100 
niter<-50
IndSppList_100_18mi<-vector("list",niter)

for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,100)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  IndSppList_100_18mi[[i]]<-IndSpp
}

IndSpp_100_18mi<-do.call("rbind",IndSppList_100_18mi)
NSpp_100_18mi<-IndSpp_100_18mi %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_100_18mi<-IndSpp_100_18mi %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
BSpp_100_18mi<-IndSpp_100_18mi %>% filter(sign_both==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="Both")
IndSpp_all_100_18mi <- rbind(NSpp_100_18mi,SSpp_100_18mi,BSpp_100_18mi)
IndSpp_all_100_18mi <-IndSpp_all_100_18mi %>% mutate(n=100,Grid="18mi")

### combine all results
IndSpp_all <- rbind(IndSpp_all_10_6mi,IndSpp_all_10_2km,IndSpp_all_10_4km,IndSpp_all_10_12mi,IndSpp_all_10_18mi,
                    IndSpp_all_20_6mi,IndSpp_all_20_2km,IndSpp_all_20_4km,IndSpp_all_20_12mi,IndSpp_all_20_18mi,
                    IndSpp_all_30_6mi,IndSpp_all_30_2km,IndSpp_all_30_4km,IndSpp_all_30_12mi,IndSpp_all_30_18mi,
                    IndSpp_all_40_6mi,IndSpp_all_40_2km,IndSpp_all_40_4km,IndSpp_all_40_12mi,IndSpp_all_40_18mi,
                    IndSpp_all_50_6mi,IndSpp_all_50_2km,IndSpp_all_50_4km,IndSpp_all_50_12mi,IndSpp_all_50_18mi,
                    IndSpp_all_60_6mi,IndSpp_all_60_2km,IndSpp_all_60_4km,IndSpp_all_60_12mi,IndSpp_all_60_18mi,
                    IndSpp_all_70_6mi,IndSpp_all_70_2km,IndSpp_all_70_4km,IndSpp_all_70_12mi,
                    IndSpp_all_80_6mi,IndSpp_all_80_2km,IndSpp_all_80_4km,IndSpp_all_80_12mi,
                    IndSpp_all_90_6mi,IndSpp_all_90_2km,IndSpp_all_90_4km,IndSpp_all_90_12mi,
                    IndSpp_all_100_6mi,IndSpp_all_100_2km,IndSpp_all_100_4km,IndSpp_all_100_12mi)


## get data ready for display
IndSpp_all_spread <- IndSpp_all %>% mutate(percent=count/50) %>% select(-count) %>% spread(SP_new,percent)
IndSpp_all_spread[is.na(IndSpp_all_spread)]<-0
IndSpp_all_long <- IndSpp_all_spread %>% gather("SP_new","percent",4:40)

# split by grid size
IndSpp_6mi_long <- IndSpp_all_long %>% filter(Grid=='6mi') %>% mutate(Zone=factor(Zone,levels=c("Both","North","South")))
IndSpp_2km_long <- IndSpp_all_long %>% filter(Grid=='2km') %>% mutate(Zone=factor(Zone,levels=c("Both","North","South")))
IndSpp_4km_long <- IndSpp_all_long %>% filter(Grid=='4km') %>% mutate(Zone=factor(Zone,levels=c("Both","North","South")))
IndSpp_12mi_long <- IndSpp_all_long %>% filter(Grid=='12mi') %>% mutate(Zone=factor(Zone,levels=c("Both","North","South")))
IndSpp_18mi_long <- IndSpp_all_long %>% filter(Grid=='18mi') %>% mutate(Zone=factor(Zone,levels=c("Both","North","South")))


### plots
# 6mi
FigS8c <- ggplot(data=IndSpp_6mi_long,aes(x=n,y=percent))+
  geom_bar(stat="identity",aes(fill=Zone))+
  facet_wrap(~SP_new,ncol = 8)+
  scale_fill_manual(values=c("gray","#2C7BB6","#D7191C"))+
  scale_x_continuous(breaks=c(10,40,70,100))+
  scale_y_continuous(breaks=c(0.0,0.5,1.0))+
  xlab("Number of samples")+
  ylab("Proportion of iterations")+
  theme_bw()+
  theme(strip.background = element_blank())+
  ggtitle("c. 9.7-kilometer grid")+
  theme(plot.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        axis.text.x = element_text(angle = 270,vjust = 0.5,hjust=0)
        )
  
# 2km
FigS8a <- ggplot(data=IndSpp_2km_long,aes(x=n,y=percent))+
  geom_bar(stat="identity",aes(fill=Zone))+
  facet_wrap(~SP_new,ncol = 8)+
  scale_fill_manual(values=c("gray","#2C7BB6","#D7191C"))+
  scale_x_continuous(breaks=c(10,40,70,100))+
  scale_y_continuous(breaks=c(0.0,0.5,1.0))+
  xlab("Number of samples")+
  ylab("Proportion of iterations")+
  theme_bw()+
  theme(strip.background = element_blank())+
  ggtitle("a. 2.0-kilometer grid")+
  theme(plot.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        axis.text.x = element_text(angle = 270,vjust = 0.5,hjust=0))

# 4km
FigS8b <- ggplot(data=IndSpp_4km_long,aes(x=n,y=percent))+
  geom_bar(stat="identity",aes(fill=Zone))+
  facet_wrap(~SP_new,ncol = 8)+
  scale_fill_manual(values=c("gray","#2C7BB6","#D7191C"))+
  scale_x_continuous(breaks=c(10,40,70,100))+
  scale_y_continuous(breaks=c(0.0,0.5,1.0))+
  xlab("Number of samples")+
  ylab("Proportion of iterations")+
  theme_bw()+
  theme(strip.background = element_blank())+
  ggtitle("b. 4.0-kilometer grid")+
  theme(plot.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        axis.text.x = element_text(angle = 270,vjust = 0.5,hjust=0))

# 12mi
FigS8d <- ggplot(data=IndSpp_12mi_long,aes(x=n,y=percent))+
  geom_bar(stat="identity",aes(fill=Zone))+
  facet_wrap(~SP_new,ncol = 8)+
  scale_fill_manual(values=c("gray","#2C7BB6","#D7191C"))+
  scale_x_continuous(breaks=c(10,40,70,100))+
  scale_y_continuous(breaks=c(0.0,0.5,1.0))+
  xlab("Number of samples")+
  ylab("Proportion of iterations")+
  theme_bw()+
  theme(strip.background = element_blank())+
  ggtitle("d. 19.3-kilometer grid")+
  theme(plot.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        axis.text.x = element_text(angle = 270,vjust = 0.5,hjust=0))

# 18mi
FigS8e <- ggplot(data=IndSpp_18mi_long,aes(x=n,y=percent))+
  geom_bar(stat="identity",aes(fill=Zone))+
  facet_wrap(~SP_new,ncol = 8)+
  scale_fill_manual(values=c("gray","#2C7BB6","#D7191C"))+
  scale_x_continuous(breaks=c(10,40,70,100))+
  scale_y_continuous(breaks=c(0.0,0.5,1.0))+
  xlab("Number of samples")+
  ylab("Proportion of iterations")+
  theme_bw()+
  theme(strip.background = element_blank())+
  ggtitle("e. 29.0-kilometer grid")+
  theme(plot.title = element_text(size = 12),
        plot.title.position = "plot",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 11),
        strip.text = element_text(size = 11),
        axis.text.x = element_text(angle = 270,vjust = 0.5,hjust=0))

#### save plots
dpi=300
ggsave("./Output/FigS8a.png", plot=FigS8a, width = 6.5, height = 4, dpi = dpi, units="in")
ggsave("./Output/FigS8b.png", plot=FigS8b, width = 7, height = 4.5, dpi = dpi, units="in")
ggsave("./Output/FigS8c.png", plot=FigS8c, width = 7, height = 4.5, dpi = dpi, units="in")
ggsave("./Output/FigS8d.png", plot=FigS8d, width = 7, height = 4.5, dpi = dpi, units="in")
ggsave("./Output/FigS8e.png", plot=FigS8e, width = 7, height = 4.5, dpi = dpi, units="in")




##################### ISA and mapping ####################
## ISA from multiple random draws
ISA_PLS<-function(north,south,nsample){
  # randomly select nsample cells from each ecoregion
  Trees_N_1<-sample_n(north,nsample,replace=FALSE)
  Trees_S_1<-sample_n(south,nsample,replace=FALSE)
  Trees_NS_1<-rbind(Trees_N_1,Trees_S_1)
  groups = c(rep(1,nsample),rep(2,nsample))
  Trees_NS_indval_1 = multipatt(Trees_NS_1,groups,control=how(nperm=9999))
  Trees_NS_indval_pvalue_1<-as.data.frame(Trees_NS_indval_1[8])
  Trees_NS_indval_pvalue_1$SP<-row.names(Trees_NS_indval_pvalue_1)
  species<-select(Trees_NS_indval_pvalue_1,SP)
  sign<-select(Trees_NS_indval_pvalue_1,sign.index)
  pvalue<-select(Trees_NS_indval_pvalue_1,sign.p.value)
  
  # formatting
  pvalue_data<-as.data.frame(pvalue)
  row.names(pvalue_data)<-unlist(species[1])
  
  # number of times each species was a significant indicator species
  pvalue_sig<-ifelse(pvalue_data<=0.05,1,0)
  pvalue_sig[is.na(pvalue_sig)]<-0
  pvalue_sig_sum<-rowSums(pvalue_sig)
  
  # number of times each species was NA indicator species (i.e. not enough of a difference between groups to test)
  pvalue_na<-ifelse(is.na(pvalue_data),1,0)
  pvalue_na_sum<-rowSums(pvalue_na)
  
  # formatting
  sign_data<-as.data.frame(sign)
  
  # northern indicator species
  sign_N<-ifelse(sign_data==1,1,0)
  sign_N_sig<-sign_N*pvalue_sig
  sign_N_sig[is.na(sign_N_sig)]<-0
  sign_N_sum<-rowSums(sign_N_sig)
  
  # southern indicator species
  sign_S<-ifelse(sign_data==2,1,0)
  sign_S_sig<-sign_S*pvalue_sig
  sign_S_sig[is.na(sign_S_sig)]<-0
  sign_S_sum<-rowSums(sign_S_sig)
  
  # both north and south
  sign_both<-ifelse(sign_data==3,1,0)
  
  #sign_B_sig<-sign_both*pvalue_sig
  sign_both[is.na(sign_both)]<-0
  sign_B_sum<-rowSums(sign_both)
  
  # formatting
  results<-data.frame(pvalue_sig,sign_N_sig,sign_S_sig,sign_both,pvalue_na)
  colnames(results)<-c('pvalue_sig','sign_N_sig','sign_S_sig','sign_both','pvalue_na')
  results$SP_new <- row.names(results)
  results <- results %>% mutate(NSB=sign_N_sig+sign_S_sig+sign_both)
  results<-results
  }


#### mapping ####
### 6 mi
niter<-50
IndSppList_6mi<-vector("list",niter)
NS_ratio_6mi <- vector("list",niter)
NS_logratio_6mi <- vector("list",niter)
NS_ratiotrans_6mi <- vector("list",niter)
TZ4_loess_6mi <- vector("list",niter)
TZ4_cutoff_6mi <- vector("list",niter)
TZ_6mi <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_6mi,Trees_count_south_6mi,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_6mi[[i]]<-IndSpp
  
  #### intensity
  ## processing
  TreesSub<-filter(Trees,SP_new %in% c(IndSpp$SP_new))
  coords<-SpatialPoints(TreesSub[,c("X","Y")])
  TreesSp<-SpatialPointsDataFrame(coords,TreesSub)
  
  # separate datasets for each region
  North<-TreesSp[TreesSp$SP_new %in% c(NSpp$SP_new),]
  South<-TreesSp[TreesSp$SP_new %in% c(SSpp$SP_new),]
  North$N = 1
  South$N = 0
  
  # set window for anlaysis
  wisc_win<-as.owin(wisc_HARN)
  
  # convert spatial points to point pattern dataset
  North.ppp<-ppp(North$X, North$Y, window=wisc_win, marks=North$N)
  South.ppp<-ppp(South$X, South$Y, window=wisc_win, marks=South$N)
  
  # count number of points in each dataset
  nNorth<-npoints(North.ppp)
  nSouth<-npoints(South.ppp)
  
  # set bandwidths
  bw.tz4<-4000
  
  ### calculate intensity
  North.dens4 <- density(North.ppp,bw.tz4,eps=1600)
  South.dens4 <- density(South.ppp,bw.tz4,eps=1600)
  
  TZratio0<-as(North.dens4,"SpatialGridDataFrame")
  names(TZratio0)<-"North"
  TZratio0$South<-as(South.dens4,"SpatialGridDataFrame")$v
  TZratio4<-as(TZratio0,"SpatialPixelsDataFrame")
  TZratio4$TZratio<-TZratio4$North/TZratio4$South
  TZratio4$logratio<-log(TZratio4$TZratio)
  TZratio4$ratio_transform<-ifelse(TZratio4$logratio<0,1/TZratio4$TZratio,TZratio4$TZratio)
  
  # store ratio and logratio from each run
  TZratioGrid <- eval.im(North.dens4/South.dens4)
  TZlogratioGrid <- eval.im(log(North.dens4/South.dens4))
  TZratiotransGrid <- eval.im(ifelse(TZlogratioGrid<0,1/TZratioGrid,TZratioGrid))
  NS_ratio_6mi[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_6mi[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_6mi[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.1, truncation point = 20000 #############################
  ############## span = 50 bins ###############
  TZratio4_nc1<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins<-hist(TZratio4_nc1$ratio_transform,breaks=seq(0,max(TZratio4_nc1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins<-data.frame(TZratio4_ratio_halfbins$mids)
  TZ4_ratio_halfbins$counts<-TZratio4_ratio_halfbins$counts
  colnames(TZ4_ratio_halfbins)<-c("bin","counts")
  TZ4_ratio_halfbins<-slice(TZ4_ratio_halfbins,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit <- loess(TZ4_ratio_halfbins$counts~TZ4_ratio_halfbins$bin,span = (50/nrow(TZ4_ratio_halfbins)))
  plot(TZ4_ratio_halfbins$bin[1:100], TZ4_ratio_halfbins$counts[1:100])
  lines(TZ4_ratio_halfbins$bin[1:100],TZ4_halfbin_fit$fitted[1:100])
  
  # local slope
  TZ4_loess_slope <- data.frame(diff(TZ4_halfbin_fit$fitted)/diff(TZ4_ratio_halfbins$bin)) %>% bind_cols(TZ4_ratio_halfbins[2:nrow(TZ4_ratio_halfbins),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_6mi[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_6mi[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_6mi[[i]]<-TZ_4
}

m_6mi<-do.call(bind,TZ_6mi)

# take mean of all ratio maps
### take mean
NS_ratio_table_6mi<-matrix(unlist(NS_ratio_6mi),ncol=niter)
NS_logratio_table_6mi<-matrix(unlist(NS_logratio_6mi),ncol=niter)
NS_ratiotrans_table_6mi<-matrix(unlist(NS_ratiotrans_6mi),ncol=niter)

NS_ratiotrans_mean_6mi<-rowMeans(NS_ratiotrans_table_6mi)

TZratio0_6mi<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_6mi)<-"North"
TZratio0_6mi$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_6mi$mean_ratiotrans<-NS_ratiotrans_mean_6mi

# percent southern
N_pct_mean_6mi <- NS_ratio_mean_6mi/(NS_ratio_mean_6mi+1)
S_pct_mean_6mi <- 1-N_pct_mean_6mi

TZratio0_6mi$N_pct <- N_pct_mean_6mi
TZratio0_6mi$S_pct <- S_pct_mean_6mi

## make spatial pixel data frame
TZratio_mean_6mi<-as(TZratio0_6mi,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_6mi,a=0.2) # set formatting for the 50 lines

TZratio_mean_6mi@data$cut_mean <- cut(TZratio_mean_6mi@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("    1","    2","    4","    7","   12","   20","   35","   60","   90","  150","  200","> 200")
SupFig2A <-spplot(TZratio_mean_6mi,"cut_mean", sp.layout=list(iterline,wisc_HARN),
              col.regions=rev(brewer.pal(11,"RdYlBu")),
              colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
              par.settings = list(axis.line = list(col =  'transparent')),
              auto.key = list(title = "Relative ratio"),
              main="9.6 kilometer grid")

SupFig2A

## Average ratio line 
avg_cutoff_6mi<-mean(unlist(TZ4_cutoff_6mi))
sd_cutoff_6mi<-sd(unlist(TZ4_cutoff_6mi))

# map of average cutoff on average ratio map
TZratio4_mean_image_6mi<-as.image.SpatialGridDataFrame(TZratio_mean_6mi["mean_ratiotrans"])
TZcontour_mean_4_6mi<-contourLines(TZratio4_mean_image_6mi,levels=avg_cutoff_6mi)
TZ_mean_4_6mi<-ContourLines2SLDF(TZcontour_mean_4_6mi)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_6mi@data$cut_mean <- cut(TZratio_mean_6mi@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S2A <-spplot(TZratio_mean_6mi,"cut_mean", sp.layout=list(TZ_mean_4_6mi,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"))

mean_plot_S2A

### 2km
niter<-50
IndSppList_2km<-vector("list",niter)
NS_ratio_2km <- vector("list",niter)
NS_logratio_2km <- vector("list",niter)
NS_ratiotrans_2km <- vector("list",niter)
TZ4_loess_2km <- vector("list",niter)
TZ4_cutoff_2km <- vector("list",niter)
TZ_2km <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_2km,Trees_count_south_2km,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_2km[[i]]<-IndSpp
  
  #### intensity
  ## processing
  TreesSub<-filter(Trees,SP_new %in% c(IndSpp$SP_new))
  coords<-SpatialPoints(TreesSub[,c("X","Y")])
  TreesSp<-SpatialPointsDataFrame(coords,TreesSub)
  
  # separate datasets for each region
  North<-TreesSp[TreesSp$SP_new %in% c(NSpp$SP_new),]
  South<-TreesSp[TreesSp$SP_new %in% c(SSpp$SP_new),]
  North$N = 1
  South$N = 0
  
  # set window for anlaysis
  wisc_win<-as.owin(wisc_HARN)
  
  # convert spatial points to point pattern dataset
  North.ppp<-ppp(North$X, North$Y, window=wisc_win, marks=North$N)
  South.ppp<-ppp(South$X, South$Y, window=wisc_win, marks=South$N)
  
  # count number of points in each dataset
  nNorth<-npoints(North.ppp)
  nSouth<-npoints(South.ppp)
  
  # set bandwidths
  bw.tz4<-4000
  
  ### calculate intensity
  North.dens4 <- density(North.ppp,bw.tz4,eps=1600)
  South.dens4 <- density(South.ppp,bw.tz4,eps=1600)
  
  TZratio0<-as(North.dens4,"SpatialGridDataFrame")
  names(TZratio0)<-"North"
  TZratio0$South<-as(South.dens4,"SpatialGridDataFrame")$v
  TZratio4<-as(TZratio0,"SpatialPixelsDataFrame")
  TZratio4$TZratio<-TZratio4$North/TZratio4$South
  TZratio4$logratio<-log(TZratio4$TZratio)
  TZratio4$ratio_transform<-ifelse(TZratio4$logratio<0,1/TZratio4$TZratio,TZratio4$TZratio)
  
  # store ratio and logratio from each run
  TZratioGrid <- eval.im(North.dens4/South.dens4)
  TZlogratioGrid <- eval.im(log(North.dens4/South.dens4))
  TZratiotransGrid <- eval.im(ifelse(TZlogratioGrid<0,1/TZratioGrid,TZratioGrid))
  NS_ratio_2km[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_2km[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_2km[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.1, truncation point = 20000 #############################
  ############## span = 50 bins ###############
  TZratio4_nc1<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins<-hist(TZratio4_nc1$ratio_transform,breaks=seq(0,max(TZratio4_nc1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins<-data.frame(TZratio4_ratio_halfbins$mids)
  TZ4_ratio_halfbins$counts<-TZratio4_ratio_halfbins$counts
  colnames(TZ4_ratio_halfbins)<-c("bin","counts")
  TZ4_ratio_halfbins<-slice(TZ4_ratio_halfbins,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit <- loess(TZ4_ratio_halfbins$counts~TZ4_ratio_halfbins$bin,span = (50/nrow(TZ4_ratio_halfbins)))
  plot(TZ4_ratio_halfbins$bin[1:100], TZ4_ratio_halfbins$counts[1:100])
  lines(TZ4_ratio_halfbins$bin[1:100],TZ4_halfbin_fit$fitted[1:100])
  
  # local slope
  TZ4_loess_slope <- data.frame(diff(TZ4_halfbin_fit$fitted)/diff(TZ4_ratio_halfbins$bin)) %>% bind_cols(TZ4_ratio_halfbins[2:nrow(TZ4_ratio_halfbins),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_2km[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_2km[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_2km[[i]]<-TZ_4
}

m_2km<-do.call(bind,TZ_2km)

# take mean of all ratio maps
### take mean
NS_ratio_table_2km<-matrix(unlist(NS_ratio_2km),ncol=niter)
NS_logratio_table_2km<-matrix(unlist(NS_logratio_2km),ncol=niter)
NS_ratiotrans_table_2km<-matrix(unlist(NS_ratiotrans_2km),ncol=niter)

NS_ratiotrans_mean_2km<-rowMeans(NS_ratiotrans_table_2km)

TZratio0_2km<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_2km)<-"North"
TZratio0_2km$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_2km$mean_ratiotrans<-NS_ratiotrans_mean_2km

# percent southern
N_pct_mean_2km <- NS_ratio_mean_2km/(NS_ratio_mean_2km+1)
S_pct_mean_2km <- 1-N_pct_mean_2km

TZratio0_2km$N_pct <- N_pct_mean_2km
TZratio0_2km$S_pct <- S_pct_mean_2km

## make spatial pixel data frame
TZratio_mean_2km<-as(TZratio0_2km,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_2km,a=0.2) # set formatting for the 50 lines

TZratio_mean_2km@data$cut_mean <- cut(TZratio_mean_2km@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("    1","    2","    4","    7","   12","   20","   35","   60","   90","  150","  200","> 200")
SupFig2B <-spplot(TZratio_mean_2km,"cut_mean", sp.layout=list(iterline,wisc_HARN),
                  col.regions=rev(brewer.pal(11,"RdYlBu")),
                  colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                  par.settings = list(axis.line = list(col =  'transparent')),
                  auto.key = list(title = "Relative ratio"),
                  main="2 kilometer grid")

SupFig2B

## Average ratio line 
avg_cutoff_2km<-mean(unlist(TZ4_cutoff_2km))
sd_cutoff_2km<-sd(unlist(TZ4_cutoff_2km))

# map of average cutoff on average ratio map
TZratio4_mean_image_2km<-as.image.SpatialGridDataFrame(TZratio_mean_2km["mean_ratiotrans"])
TZcontour_mean_4_2km<-contourLines(TZratio4_mean_image_2km,levels=avg_cutoff_2km)
TZ_mean_4_2km<-ContourLines2SLDF(TZcontour_mean_4_2km)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_2km@data$cut_mean <- cut(TZratio_mean_2km@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S2B <-spplot(TZratio_mean_2km,"cut_mean", sp.layout=list(TZ_mean_4_2km,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"))

mean_plot_S2B


### 4km
niter<-50
IndSppList_4km<-vector("list",niter)
NS_ratio_4km <- vector("list",niter)
NS_logratio_4km <- vector("list",niter)
NS_ratiotrans_4km <- vector("list",niter)
TZ4_loess_4km <- vector("list",niter)
TZ4_cutoff_4km <- vector("list",niter)
TZ_4km <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_4km,Trees_count_south_4km,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_4km[[i]]<-IndSpp
  
  #### intensity
  ## processing
  TreesSub<-filter(Trees,SP_new %in% c(IndSpp$SP_new))
  coords<-SpatialPoints(TreesSub[,c("X","Y")])
  TreesSp<-SpatialPointsDataFrame(coords,TreesSub)
  
  # separate datasets for each region
  North<-TreesSp[TreesSp$SP_new %in% c(NSpp$SP_new),]
  South<-TreesSp[TreesSp$SP_new %in% c(SSpp$SP_new),]
  North$N = 1
  South$N = 0
  
  # set window for anlaysis
  wisc_win<-as.owin(wisc_HARN)
  
  # convert spatial points to point pattern dataset
  North.ppp<-ppp(North$X, North$Y, window=wisc_win, marks=North$N)
  South.ppp<-ppp(South$X, South$Y, window=wisc_win, marks=South$N)
  
  # count number of points in each dataset
  nNorth<-npoints(North.ppp)
  nSouth<-npoints(South.ppp)
  
  # set bandwidths
  bw.tz4<-4000
  
  ### calculate intensity
  North.dens4 <- density(North.ppp,bw.tz4,eps=1600)
  South.dens4 <- density(South.ppp,bw.tz4,eps=1600)
  
  TZratio0<-as(North.dens4,"SpatialGridDataFrame")
  names(TZratio0)<-"North"
  TZratio0$South<-as(South.dens4,"SpatialGridDataFrame")$v
  TZratio4<-as(TZratio0,"SpatialPixelsDataFrame")
  TZratio4$TZratio<-TZratio4$North/TZratio4$South
  TZratio4$logratio<-log(TZratio4$TZratio)
  TZratio4$ratio_transform<-ifelse(TZratio4$logratio<0,1/TZratio4$TZratio,TZratio4$TZratio)
  
  # store ratio and logratio from each run
  TZratioGrid <- eval.im(North.dens4/South.dens4)
  TZlogratioGrid <- eval.im(log(North.dens4/South.dens4))
  TZratiotransGrid <- eval.im(ifelse(TZlogratioGrid<0,1/TZratioGrid,TZratioGrid))
  NS_ratio_4km[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_4km[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_4km[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.1, truncation point = 20000 #############################
  ############## span = 50 bins ###############
  TZratio4_nc1<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins<-hist(TZratio4_nc1$ratio_transform,breaks=seq(0,max(TZratio4_nc1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins<-data.frame(TZratio4_ratio_halfbins$mids)
  TZ4_ratio_halfbins$counts<-TZratio4_ratio_halfbins$counts
  colnames(TZ4_ratio_halfbins)<-c("bin","counts")
  TZ4_ratio_halfbins<-slice(TZ4_ratio_halfbins,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit <- loess(TZ4_ratio_halfbins$counts~TZ4_ratio_halfbins$bin,span = (50/nrow(TZ4_ratio_halfbins)))
  plot(TZ4_ratio_halfbins$bin[1:100], TZ4_ratio_halfbins$counts[1:100])
  lines(TZ4_ratio_halfbins$bin[1:100],TZ4_halfbin_fit$fitted[1:100])
  
  # local slope
  TZ4_loess_slope <- data.frame(diff(TZ4_halfbin_fit$fitted)/diff(TZ4_ratio_halfbins$bin)) %>% bind_cols(TZ4_ratio_halfbins[2:nrow(TZ4_ratio_halfbins),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_4km[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_4km[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_4km[[i]]<-TZ_4
}

m_4km<-do.call(bind,TZ_4km)

# take mean of all ratio maps
### take mean
NS_ratio_table_4km<-matrix(unlist(NS_ratio_4km),ncol=niter)
NS_logratio_table_4km<-matrix(unlist(NS_logratio_4km),ncol=niter)
NS_ratiotrans_table_4km<-matrix(unlist(NS_ratiotrans_4km),ncol=niter)

NS_ratiotrans_mean_4km<-rowMeans(NS_ratiotrans_table_4km)

TZratio0_4km<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_4km)<-"North"
TZratio0_4km$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_4km$mean_ratiotrans<-NS_ratiotrans_mean_4km

# percent southern
N_pct_mean_4km <- NS_ratio_mean_4km/(NS_ratio_mean_4km+1)
S_pct_mean_4km <- 1-N_pct_mean_4km

TZratio0_4km$N_pct <- N_pct_mean_4km
TZratio0_4km$S_pct <- S_pct_mean_4km

## make spatial pixel data frame
TZratio_mean_4km<-as(TZratio0_4km,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_4km,a=0.2) # set formatting for the 50 lines

TZratio_mean_4km@data$cut_mean <- cut(TZratio_mean_4km@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("    1","    2","    4","    7","   12","   20","   35","   60","   90","  150","  200","> 200")
SupFig2C <-spplot(TZratio_mean_4km,"cut_mean", sp.layout=list(iterline,wisc_HARN),
                  col.regions=rev(brewer.pal(11,"RdYlBu")),
                  colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                  par.settings = list(axis.line = list(col =  'transparent')),
                  auto.key = list(title = "Relative ratio"),
                  main="4 kilometer grid")

SupFig2C

## Average ratio line 
avg_cutoff_4km<-mean(unlist(TZ4_cutoff_4km))
sd_cutoff_4km<-sd(unlist(TZ4_cutoff_4km))

# map of average cutoff on average ratio map
TZratio4_mean_image_4km<-as.image.SpatialGridDataFrame(TZratio_mean_4km["mean_ratiotrans"])
TZcontour_mean_4_4km<-contourLines(TZratio4_mean_image_4km,levels=avg_cutoff_4km)
TZ_mean_4_4km<-ContourLines2SLDF(TZcontour_mean_4_4km)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_4km@data$cut_mean <- cut(TZratio_mean_4km@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S2C <-spplot(TZratio_mean_4km,"cut_mean", sp.layout=list(TZ_mean_4_4km,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"))

mean_plot_S2C

### 12mi
niter<-50
IndSppList_12mi<-vector("list",niter)
NS_ratio_12mi <- vector("list",niter)
NS_logratio_12mi <- vector("list",niter)
NS_ratiotrans_12mi <- vector("list",niter)
TZ4_loess_12mi <- vector("list",niter)
TZ4_cutoff_12mi <- vector("list",niter)
TZ_12mi <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_12mi,Trees_count_south_12mi,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_12mi[[i]]<-IndSpp
  
  #### intensity
  ## processing
  TreesSub<-filter(Trees,SP_new %in% c(IndSpp$SP_new))
  coords<-SpatialPoints(TreesSub[,c("X","Y")])
  TreesSp<-SpatialPointsDataFrame(coords,TreesSub)
  
  # separate datasets for each region
  North<-TreesSp[TreesSp$SP_new %in% c(NSpp$SP_new),]
  South<-TreesSp[TreesSp$SP_new %in% c(SSpp$SP_new),]
  North$N = 1
  South$N = 0
  
  # set window for anlaysis
  wisc_win<-as.owin(wisc_HARN)
  
  # convert spatial points to point pattern dataset
  North.ppp<-ppp(North$X, North$Y, window=wisc_win, marks=North$N)
  South.ppp<-ppp(South$X, South$Y, window=wisc_win, marks=South$N)
  
  # count number of points in each dataset
  nNorth<-npoints(North.ppp)
  nSouth<-npoints(South.ppp)
  
  # set bandwidths
  bw.tz4<-4000
  
  ### calculate intensity
  North.dens4 <- density(North.ppp,bw.tz4,eps=1600)
  South.dens4 <- density(South.ppp,bw.tz4,eps=1600)
  
  TZratio0<-as(North.dens4,"SpatialGridDataFrame")
  names(TZratio0)<-"North"
  TZratio0$South<-as(South.dens4,"SpatialGridDataFrame")$v
  TZratio4<-as(TZratio0,"SpatialPixelsDataFrame")
  TZratio4$TZratio<-TZratio4$North/TZratio4$South
  TZratio4$logratio<-log(TZratio4$TZratio)
  TZratio4$ratio_transform<-ifelse(TZratio4$logratio<0,1/TZratio4$TZratio,TZratio4$TZratio)
  
  # store ratio and logratio from each run
  TZratioGrid <- eval.im(North.dens4/South.dens4)
  TZlogratioGrid <- eval.im(log(North.dens4/South.dens4))
  TZratiotransGrid <- eval.im(ifelse(TZlogratioGrid<0,1/TZratioGrid,TZratioGrid))
  NS_ratio_12mi[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_12mi[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_12mi[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.1, truncation point = 20000 #############################
  ############## span = 50 bins ###############
  TZratio4_nc1<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins<-hist(TZratio4_nc1$ratio_transform,breaks=seq(0,max(TZratio4_nc1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins<-data.frame(TZratio4_ratio_halfbins$mids)
  TZ4_ratio_halfbins$counts<-TZratio4_ratio_halfbins$counts
  colnames(TZ4_ratio_halfbins)<-c("bin","counts")
  TZ4_ratio_halfbins<-slice(TZ4_ratio_halfbins,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit <- loess(TZ4_ratio_halfbins$counts~TZ4_ratio_halfbins$bin,span = (50/nrow(TZ4_ratio_halfbins)))
  plot(TZ4_ratio_halfbins$bin[1:100], TZ4_ratio_halfbins$counts[1:100])
  lines(TZ4_ratio_halfbins$bin[1:100],TZ4_halfbin_fit$fitted[1:100])
  
  # local slope
  TZ4_loess_slope <- data.frame(diff(TZ4_halfbin_fit$fitted)/diff(TZ4_ratio_halfbins$bin)) %>% bind_cols(TZ4_ratio_halfbins[2:nrow(TZ4_ratio_halfbins),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_12mi[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_12mi[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_12mi[[i]]<-TZ_4
}

m_12mi<-do.call(bind,TZ_12mi)

# take mean of all ratio maps
### take mean
NS_ratio_table_12mi<-matrix(unlist(NS_ratio_12mi),ncol=niter)
NS_logratio_table_12mi<-matrix(unlist(NS_logratio_12mi),ncol=niter)
NS_ratiotrans_table_12mi<-matrix(unlist(NS_ratiotrans_12mi),ncol=niter)

NS_ratiotrans_mean_12mi<-rowMeans(NS_ratiotrans_table_12mi)

TZratio0_12mi<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_12mi)<-"North"
TZratio0_12mi$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_12mi$mean_ratiotrans<-NS_ratiotrans_mean_12mi

# percent southern
N_pct_mean_12mi <- NS_ratio_mean_12mi/(NS_ratio_mean_12mi+1)
S_pct_mean_12mi <- 1-N_pct_mean_12mi

TZratio0_12mi$N_pct <- N_pct_mean_12mi
TZratio0_12mi$S_pct <- S_pct_mean_12mi

## make spatial pixel data frame
TZratio_mean_12mi<-as(TZratio0_12mi,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_12mi,a=0.2) # set formatting for the 50 lines

TZratio_mean_12mi@data$cut_mean <- cut(TZratio_mean_12mi@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("    1","    2","    4","    7","   12","   20","   35","   60","   90","  150","  200","> 200")
SupFig2D <-spplot(TZratio_mean_12mi,"cut_mean", sp.layout=list(iterline,wisc_HARN),
                  col.regions=rev(brewer.pal(11,"RdYlBu")),
                  colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                  par.settings = list(axis.line = list(col =  'transparent')),
                  auto.key = list(title = "Relative ratio"),
                  main="19.3 kilometer grid")

SupFig2D

## Average ratio line 
avg_cutoff_12mi<-mean(unlist(TZ4_cutoff_12mi))
sd_cutoff_12mi<-sd(unlist(TZ4_cutoff_12mi))

# map of average cutoff on average ratio map
TZratio4_mean_image_12mi<-as.image.SpatialGridDataFrame(TZratio_mean_12mi["mean_ratiotrans"])
TZcontour_mean_4_12mi<-contourLines(TZratio4_mean_image_12mi,levels=avg_cutoff_12mi)
TZ_mean_4_12mi<-ContourLines2SLDF(TZcontour_mean_4_12mi)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_12mi@data$cut_mean <- cut(TZratio_mean_12mi@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S2D <-spplot(TZratio_mean_12mi,"cut_mean", sp.layout=list(TZ_mean_4_12mi,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"))

mean_plot_S2D


### 18mi
niter<-50
IndSppList_18mi<-vector("list",niter)
NS_ratio_18mi <- vector("list",niter)
NS_logratio_18mi <- vector("list",niter)
NS_ratiotrans_18mi <- vector("list",niter)
TZ4_loess_18mi <- vector("list",niter)
TZ4_cutoff_18mi <- vector("list",niter)
TZ_18mi <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north_18mi,Trees_count_south_18mi,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_18mi[[i]]<-IndSpp
  
  #### intensity
  ## processing
  TreesSub<-filter(Trees,SP_new %in% c(IndSpp$SP_new))
  coords<-SpatialPoints(TreesSub[,c("X","Y")])
  TreesSp<-SpatialPointsDataFrame(coords,TreesSub)
  
  # separate datasets for each region
  North<-TreesSp[TreesSp$SP_new %in% c(NSpp$SP_new),]
  South<-TreesSp[TreesSp$SP_new %in% c(SSpp$SP_new),]
  North$N = 1
  South$N = 0
  
  # set window for anlaysis
  wisc_win<-as.owin(wisc_HARN)
  
  # convert spatial points to point pattern dataset
  North.ppp<-ppp(North$X, North$Y, window=wisc_win, marks=North$N)
  South.ppp<-ppp(South$X, South$Y, window=wisc_win, marks=South$N)
  
  # count number of points in each dataset
  nNorth<-npoints(North.ppp)
  nSouth<-npoints(South.ppp)
  
  # set bandwidths
  bw.tz4<-4000
  
  ### calculate intensity
  North.dens4 <- density(North.ppp,bw.tz4,eps=1600)
  South.dens4 <- density(South.ppp,bw.tz4,eps=1600)
  
  TZratio0<-as(North.dens4,"SpatialGridDataFrame")
  names(TZratio0)<-"North"
  TZratio0$South<-as(South.dens4,"SpatialGridDataFrame")$v
  TZratio4<-as(TZratio0,"SpatialPixelsDataFrame")
  TZratio4$TZratio<-TZratio4$North/TZratio4$South
  TZratio4$logratio<-log(TZratio4$TZratio)
  TZratio4$ratio_transform<-ifelse(TZratio4$logratio<0,1/TZratio4$TZratio,TZratio4$TZratio)
  
  # store ratio and logratio from each run
  TZratioGrid <- eval.im(North.dens4/South.dens4)
  TZlogratioGrid <- eval.im(log(North.dens4/South.dens4))
  TZratiotransGrid <- eval.im(ifelse(TZlogratioGrid<0,1/TZratioGrid,TZratioGrid))
  NS_ratio_18mi[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_18mi[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_18mi[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.1, truncation point = 20000 #############################
  ############## span = 50 bins ###############
  TZratio4_nc1<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins<-hist(TZratio4_nc1$ratio_transform,breaks=seq(0,max(TZratio4_nc1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins<-data.frame(TZratio4_ratio_halfbins$mids)
  TZ4_ratio_halfbins$counts<-TZratio4_ratio_halfbins$counts
  colnames(TZ4_ratio_halfbins)<-c("bin","counts")
  TZ4_ratio_halfbins<-slice(TZ4_ratio_halfbins,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit <- loess(TZ4_ratio_halfbins$counts~TZ4_ratio_halfbins$bin,span = (50/nrow(TZ4_ratio_halfbins)))
  plot(TZ4_ratio_halfbins$bin[1:100], TZ4_ratio_halfbins$counts[1:100])
  lines(TZ4_ratio_halfbins$bin[1:100],TZ4_halfbin_fit$fitted[1:100])
  
  # local slope
  TZ4_loess_slope <- data.frame(diff(TZ4_halfbin_fit$fitted)/diff(TZ4_ratio_halfbins$bin)) %>% bind_cols(TZ4_ratio_halfbins[2:nrow(TZ4_ratio_halfbins),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_18mi[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_18mi[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_18mi[[i]]<-TZ_4
}

m_18mi<-do.call(bind,TZ_18mi)

# take mean of all ratio maps
### take mean
NS_ratio_table_18mi<-matrix(unlist(NS_ratio_18mi),ncol=niter)
NS_logratio_table_18mi<-matrix(unlist(NS_logratio_18mi),ncol=niter)
NS_ratiotrans_table_18mi<-matrix(unlist(NS_ratiotrans_18mi),ncol=niter)

NS_ratiotrans_mean_18mi<-rowMeans(NS_ratiotrans_table_18mi)

TZratio0_18mi<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_18mi)<-"North"
TZratio0_18mi$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_18mi$mean_ratiotrans<-NS_ratiotrans_mean_18mi

## make spatial pixel data frame
TZratio_mean_18mi<-as(TZratio0_18mi,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_18mi,a=0.2) # set formatting for the 50 lines

TZratio_mean_18mi@data$cut_mean <- cut(TZratio_mean_18mi@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("    1","    2","    4","    7","   12","   20","   35","   60","   90","  150","  200","> 200")
SupFig2E <-spplot(TZratio_mean_18mi,"cut_mean", sp.layout=list(iterline,wisc_HARN),
                  col.regions=rev(brewer.pal(11,"RdYlBu")),
                  colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                  par.settings = list(axis.line = list(col =  'transparent')),
                  auto.key = list(title = "Relative ratio"),
                  main="29 kilometer grid")

SupFig2E

## Average ratio line 
avg_cutoff_18mi<-mean(unlist(TZ4_cutoff_18mi))
sd_cutoff_18mi<-sd(unlist(TZ4_cutoff_18mi))

# map of average cutoff on average ratio map
TZratio4_mean_image_18mi<-as.image.SpatialGridDataFrame(TZratio_mean_18mi["mean_ratiotrans"])
TZcontour_mean_4_18mi<-contourLines(TZratio4_mean_image_18mi,levels=avg_cutoff_18mi)
TZ_mean_4_18mi<-ContourLines2SLDF(TZcontour_mean_4_18mi)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_18mi@data$cut_mean <- cut(TZratio_mean_18mi@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S2E <-spplot(TZratio_mean_18mi,"cut_mean", sp.layout=list(TZ_mean_4_18mi,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"))

mean_plot_S2E


plot(wisc_HARN)
plot(Grid_6mi,add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

plot(wisc_HARN)
plot(Grid_2km,add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

plot(wisc_HARN)
plot(Grid_4km,add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

plot(wisc_HARN)
plot(Grid_12mi,add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

plot(wisc_HARN)
plot(Grid_18mi,add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')


##### plot mean lines together #####
plot(wisc_HARN)
plot(TZ_mean_4_2km, col="red",add=TRUE)
plot(TZ_mean_4_4km,col="orange",add=TRUE)
plot(TZ_mean_4_6mi,col="green",add=TRUE)
plot(TZ_mean_4_12mi,col="blue",add=TRUE)
plot(TZ_mean_4_18mi,col="purple",add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

#### save image #####
### convert to sf
wisc_HARN_sf <- as(wisc_HARN,"sf")

## combine lines
# add field with variable
TZ_mean_4_2km$gridsize <- "2km"
TZ_mean_4_4km$gridsize <- "4km"
TZ_mean_4_6mi$gridsize <- "9.6km"
TZ_mean_4_12mi$gridsize <- "19.3km"
TZ_mean_4_18mi$gridsize <- "29.0km"

# combine FS lines
TZ_mean_4_grid_all <- rbind(TZ_mean_4_2km,TZ_mean_4_4km,TZ_mean_4_6mi,TZ_mean_4_12mi,TZ_mean_4_18mi)
crs(TZ_mean_4_grid_all) <- "+init=epsg:3070"
TZ_mean_4_grid_all$gridsize <- factor(TZ_mean_4_grid_all$gridsize,levels = c("2km","4km","9.6km","19.3km","29.0km"))
TZ_mean_4_grid_all_sf <- as(TZ_mean_4_grid_all,"sf")

### plot
FigS3 <- ggplot(data=wisc_HARN_sf)+
  geom_sf(fill="white")+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1, "cm"), #how tall the arrow should be
                         width= unit(0.5, "cm"), 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_orienteering)+
  geom_sf(data=TZ_mean_4_grid_all_sf,aes(col=gridsize,fill=gridsize))+
  scale_color_manual(breaks=c("2km","4km","9.6km","19.3km","29.0km"),
                     values=alpha(c("2km"="red","4km"="orange","9.6km"="green","19.3km"="blue","29.0km"="purple"),0.8),
                     labels=c("2.0 km","4.0 km","9.6 km","19.3 km","29.0 km"),guide="none")+
  scale_fill_manual(values=c("2km"="red","4km"="orange","9.6km"="green","19.3km"="blue","29.0km"="purple"),
                    labels=c("2.0 km","4.0 km","9.6 km","19.3 km","29.0 km"))+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=10))

#### save plot
dpi=600
ggsave("./Output/FigS3.jpg", plot=FigS3, width = 4000/dpi, height = 4000/dpi, dpi = dpi)


