# set your own working directory to where the files are stored on your computer
setwd("U:/Shea/Research/GIS/Projects/TZ_Boundary/github")

library(indicspecies)
library(maptools)
library(spatstat)
library(sp)
library(rgdal)
library(RColorBrewer)
library(ggplot2)
library(ggspatial)
library(sf)
library(tidyr)
library(raster)
library(dplyr)


##################### move boundary to determine impact on results #############################

################### prepare data #######################
## create polygons for zones
# bring in WI map
wisc<-readOGR(dsn="./WI_HARN_mask",layer="WI_HARN_mask")
wisc_HARN <- spTransform(wisc,CRS=CRS("+init=epsg:3070"))

## import 6 mi grid
Grid_6mi <- readOGR(dsn="./6mi_Grid",layer="6mi_Grid")
Grid_6mi <- spTransform(Grid_6mi,CRS=CRS("+init=epsg:3070"))

# make overall polygon
wisc_coords<-wisc_HARN@bbox
wisc_poly<-as(raster::extent(wisc_coords),"SpatialPolygons")
proj4string(wisc_poly)<-"+init=epsg:3070"

wisc_coords_df<-as.data.frame(wisc_coords)
x_expand<-wisc_coords_df[1,] %>% mutate(max2=max+(max-min),min2=min-(max-min)) %>% mutate(max=max2,min=min2) %>% select(min,max)

# separate polygon for N
wisc_coords_N<- as.data.frame(wisc_coords) %>% 
  slice(2) %>% 
  mutate(avg=(min+max)/2) %>%
  select(-min) %>%
  rename(min=avg) %>% 
  mutate(max=(max+((max-min)*2))) %>% 
  bind_rows(x_expand) %>% 
  arrange(desc(max))
row.names(wisc_coords_N)<-c("x","y")
wisc_coords_N<-as.matrix(wisc_coords_N)
wisc_poly_N<-as(raster::extent(wisc_coords_N),"SpatialPolygons")
proj4string(wisc_poly_N)<-"+init=epsg:3070"
zone<-as.data.frame("north")
colnames(zone)<-"zone"
wisc_polyDF_N<-SpatialPolygonsDataFrame(wisc_poly_N,zone) # add data on zone

# separate polygon for S
wisc_coords_S<- as.data.frame(wisc_coords) %>% 
  slice(2) %>% 
  mutate(avg=(min+max)/2) %>% 
  select(-max) %>% 
  rename(max=avg) %>% 
  mutate(min=(min-((max-min)*2))) %>% 
  bind_rows(x_expand) %>% 
  arrange(desc(max))
row.names(wisc_coords_S)<-c("x","y")
wisc_coords_S<-as.matrix(wisc_coords_S)
wisc_poly_S<-as(raster::extent(wisc_coords_S),"SpatialPolygons")
proj4string(wisc_poly_S)<-"+init=epsg:3070"
zone<-as.data.frame("south",colnames("zone"))
colnames(zone)<-"zone"
wisc_polyDF_S<-SpatialPolygonsDataFrame(wisc_poly_S,zone) # add data on zone

# join polygons
wisc_list<-list(wisc_polyDF_N,wisc_polyDF_S)
makeUniform<-function(SPDF){
  pref<-substitute(SPDF)  #just putting the file name in front.
  newSPDF<-spChFIDs(SPDF,as.character(paste(pref,rownames(as(SPDF,"data.frame")),sep="_")))
  return(newSPDF)
} ### not sure why this is necessary but it makes the "rotate zones" section work
newIDs<-lapply(wisc_list,function(x) makeUniform(x))
wisc_split<-do.call(rbind,newIDs)

# shift zones N/S
wisc_split_50kS<-raster::shift(wisc_split,0,-50000)
wisc_split_50kN<-raster::shift(wisc_split,0,50000)

# rotate zones
wisc_split_45deg<-elide(wisc_split,rotate=45,center=apply(wisc_coords,1,mean))
wisc_split_135deg<-elide(wisc_split,rotate=135,center=apply(wisc_coords,1,mean))

proj4string(wisc_split_45deg) <- CRS("+init=epsg:3070")
proj4string(wisc_split_135deg) <- CRS("+init=epsg:3070")

## bring in Forest Service ecoregion boudnary
Zones <- readOGR("./212_222_boundary","212_222_merge_extended")
Zones <- spTransform(Zones,CRS=CRS("+init=epsg:3070"))

# shift zones south and north
Zones_50kS<-raster::shift(Zones,0,-50000)
Zones_50kN<-raster::shift(Zones,0,50000)

# add FS zones table (NOTE: this was greated in ArcGIS to determine if cells split by ecoregion boundary had the majority of their area in N or S)
Grid_zone<-read.csv(file="6mi_Grid_ecoregions_area.csv",header=TRUE)

# select grid ecoregions with largest area (this is only for the original boundary location)
Grid_ecoregion <- Grid_zone %>% group_by(Grid_ID) %>% summarise(max=max(F_AREA)) %>% left_join(Grid_zone,by="Grid_ID") %>% 
  filter(max==F_AREA) %>% select(Grid_ID,Ecoregion)

# add to spatial grid
Grid_6mi <- merge(Grid_6mi,Grid_ecoregion,by="Grid_ID")

# convert grid to raster
r<-raster(ncol=51,nrow=54)
extent(r)<-extent(Grid_6mi)
GridRaster<-rasterize(Grid_6mi,r,'Grid_ID')
crs(GridRaster) <- "+init=epsg:3070"

## join grid id to zones
# Normal
GridZones<-raster::extract(GridRaster,wisc_split,df=TRUE)
GridZones$Zone<-ifelse (GridZones$ID<2,"North","South")
GridZones<-GridZones %>% rename(Grid_ID=layer) %>% select(-ID)

# Shift N
GridZones_50kN<-raster::extract(GridRaster,wisc_split_50kN,df=TRUE)
GridZones_50kN$Zone_50kN<-ifelse (GridZones_50kN$ID<2,"North","South")
GridZones_50kN<-GridZones_50kN %>% rename(Grid_ID=layer) %>% select(-ID)

# Shift S
GridZones_50kS<-raster::extract(GridRaster,wisc_split_50kS,df=TRUE)
GridZones_50kS$Zone_50kS<-ifelse (GridZones_50kS$ID<2,"North","South")
GridZones_50kS<-GridZones_50kS %>% rename(Grid_ID=layer) %>% select(-ID)

# Rotate 45
GridZones_45deg<-raster::extract(GridRaster,wisc_split_45deg,df=TRUE)
GridZones_45deg$Zone_45deg<-ifelse (GridZones_45deg$ID<2,"North","South")
GridZones_45deg<-GridZones_45deg %>% rename(Grid_ID=layer) %>% select(-ID)

# Rotate 135
GridZones_135deg<-raster::extract(GridRaster,wisc_split_135deg,df=TRUE)
GridZones_135deg$Zone_135deg<-ifelse (GridZones_135deg$ID<2,"North","South")
GridZones_135deg<-GridZones_135deg %>% rename(Grid_ID=layer) %>% select(-ID)

# FS Normal
GridZonesFS<-raster::extract(GridRaster,Zones,df=TRUE)
GridZonesFS$ZoneFS<-ifelse (GridZonesFS$ID<2,"South","North")
GridZonesFS<-GridZonesFS %>% rename(Grid_ID=layer) %>% select(-ID)

# FS Shift S
GridZonesFS_50kS<-raster::extract(GridRaster,Zones_50kS,df=TRUE)
GridZonesFS_50kS$ZoneFS_50kS<-ifelse (GridZonesFS_50kS$ID<2,"South","North")
GridZonesFS_50kS<-GridZonesFS_50kS %>% rename(Grid_ID=layer) %>% select(-ID)

# FS Shift N
GridZonesFS_50kN<-raster::extract(GridRaster,Zones_50kN,df=TRUE)
GridZonesFS_50kN$ZoneFS_50kN<-ifelse (GridZonesFS_50kN$ID<2,"South","North")
GridZonesFS_50kN<-GridZonesFS_50kN %>% rename(Grid_ID=layer) %>% select(-ID)

GridList<-list(GridZones,GridZones_50kN,GridZones_50kS,GridZones_45deg,GridZones_135deg,GridZonesFS,GridZonesFS_50kN,GridZonesFS_50kS)
GridZones_all<-Reduce(function(x,y) merge(x,y,by="Grid_ID",all=TRUE), GridList)


### identify cells than fall entirely inside state boundary (excluding edge cells)
Grid_WI_6mi <- raster::extract(GridRaster,wisc_HARN,df=TRUE,weights=TRUE)
Grid_WI_6mi<-Grid_WI_6mi %>% rename(Grid_ID=layer) %>% filter(weight==max(weight))


## import species data
all_trees_pred <- read.csv(file="all_trees_pred_12Aug2019.csv",header=TRUE)

# make spatial
coords <- cbind(all_trees_pred$X,all_trees_pred$Y)
ptree_diff_2018 <- SpatialPointsDataFrame(coords,all_trees_pred)
proj4string(ptree_diff_2018) <- CRS("+init=epsg:3070") # add projection NAD83/Wisconsin Transferse Mercator

## determine which trees belong to each grid cell
pg_overlay <- over(ptree_diff_2018,Grid_6mi)
Trees <- cbind(pg_overlay,ptree_diff_2018@data)

# remove NA, UK, RR (NA, unknown species, and very rare species with <10 points)
Trees<-Trees %>% filter(!SP_new %in% c(' ','UK','RR')) %>% filter(!is.na(SP_new))
Trees$SP_new<-as.factor(Trees$SP_new)
Trees$SP_new<-droplevels(Trees$SP_new)

# convert to species counts for each cell, remove edge cells, and add ecoregion
Trees_count_ecoregion <- Trees %>% group_by(SP_new,Grid_ID) %>% summarise(count=n()) %>% spread(SP_new,count) %>% 
  semi_join(Grid_WI_6mi,by="Grid_ID") %>% left_join(GridZones_all,by="Grid_ID")
Trees_count_ecoregion[is.na(Trees_count_ecoregion)]<-0 # change NA to 0


# separate matrices for each ecoregion
# north normal
Trees_count_N_hor<-filter(Trees_count_ecoregion,Zone=="North")
row.names(Trees_count_N_hor)<-Trees_count_N_hor$Grid_ID
Trees_count_N_hor<-select(Trees_count_N_hor,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# south normal
Trees_count_S_hor<-filter(Trees_count_ecoregion,Zone=="South")
row.names(Trees_count_S_hor)<-Trees_count_S_hor$Grid_ID
Trees_count_S_hor<-select(Trees_count_S_hor,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# north 50kN
Trees_count_N_50kN<-filter(Trees_count_ecoregion,Zone_50kN=="North")
row.names(Trees_count_N_50kN)<-Trees_count_N_50kN$Grid_ID
Trees_count_N_50kN<-select(Trees_count_N_50kN,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# south 50kN
Trees_count_S_50kN<-filter(Trees_count_ecoregion,Zone_50kN=="South")
row.names(Trees_count_S_50kN)<-Trees_count_S_50kN$Grid_ID
Trees_count_S_50kN<-select(Trees_count_S_50kN,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# north 50kS
Trees_count_N_50kS<-filter(Trees_count_ecoregion,Zone_50kS=="North")
row.names(Trees_count_N_50kS)<-Trees_count_N_50kS$Grid_ID
Trees_count_N_50kS<-select(Trees_count_N_50kS,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# south 50kS
Trees_count_S_50kS<-filter(Trees_count_ecoregion,Zone_50kS=="South")
row.names(Trees_count_S_50kS)<-Trees_count_S_50kS$Grid_ID
Trees_count_S_50kS<-select(Trees_count_S_50kS,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# north 45deg
Trees_count_N_45deg<-filter(Trees_count_ecoregion,Zone_45deg=="North")
row.names(Trees_count_N_45deg)<-Trees_count_N_45deg$Grid_ID
Trees_count_N_45deg<-select(Trees_count_N_45deg,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# south 45deg
Trees_count_S_45deg<-filter(Trees_count_ecoregion,Zone_45deg=="South")
row.names(Trees_count_S_45deg)<-Trees_count_S_45deg$Grid_ID
Trees_count_S_45deg<-select(Trees_count_S_45deg,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# north 135deg
Trees_count_N_135deg<-filter(Trees_count_ecoregion,Zone_135deg=="North")
row.names(Trees_count_N_135deg)<-Trees_count_N_135deg$Grid_ID
Trees_count_N_135deg<-select(Trees_count_N_135deg,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# south 135deg
Trees_count_S_135deg<-filter(Trees_count_ecoregion,Zone_135deg=="South")
row.names(Trees_count_S_135deg)<-Trees_count_S_135deg$Grid_ID
Trees_count_S_135deg<-select(Trees_count_S_135deg,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# north normal FS
Trees_count_FS_N<-filter(Trees_count_ecoregion,ZoneFS=="North")
row.names(Trees_count_FS_N)<-Trees_count_FS_N$Grid_ID
Trees_count_FS_N<-select(Trees_count_FS_N,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# south normal FS
Trees_count_FS_S<-filter(Trees_count_ecoregion,ZoneFS=="South")
row.names(Trees_count_FS_S)<-Trees_count_FS_S$Grid_ID
Trees_count_FS_S<-select(Trees_count_FS_S,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# north 50kN FS
Trees_count_FS_N_50kN<-filter(Trees_count_ecoregion,ZoneFS_50kN=="North")
row.names(Trees_count_FS_N_50kN)<-Trees_count_FS_N_50kN$Grid_ID
Trees_count_FS_N_50kN<-select(Trees_count_FS_N_50kN,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# south 50kN FS
Trees_count_FS_S_50kN<-filter(Trees_count_ecoregion,ZoneFS_50kN=="South")
row.names(Trees_count_FS_S_50kN)<-Trees_count_FS_S_50kN$Grid_ID
Trees_count_FS_S_50kN<-select(Trees_count_FS_S_50kN,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# north 50kS FS
Trees_count_FS_N_50kS<-filter(Trees_count_ecoregion,ZoneFS_50kS=="North")
row.names(Trees_count_FS_N_50kS)<-Trees_count_FS_N_50kS$Grid_ID
Trees_count_FS_N_50kS<-select(Trees_count_FS_N_50kS,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))

# south 50kS FS
Trees_count_FS_S_50kS<-filter(Trees_count_ecoregion,ZoneFS_50kS=="South")
row.names(Trees_count_FS_S_50kS)<-Trees_count_FS_S_50kS$Grid_ID
Trees_count_FS_S_50kS<-select(Trees_count_FS_S_50kS,-c(Grid_ID,Zone,Zone_50kS,Zone_50kN,Zone_45deg,Zone_135deg,ZoneFS,ZoneFS_50kN,ZoneFS_50kS))


##################### ISA and mapping #####################
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
######## horizontal line ##########
niter<-20
IndSppList_hor<-vector("list",niter)
NS_ratio_hor <- vector("list",niter)
NS_logratio_hor <- vector("list",niter)
NS_ratiotrans_hor <- vector("list",niter)
TZ4_loess_hor <- vector("list",niter)
TZ4_cutoff_hor <- vector("list",niter)
TZ_hor <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_N_hor,Trees_count_S_hor,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_hor[[i]]<-IndSpp
  
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
  NS_ratio_hor[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_hor[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_hor[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_hor[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_hor[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_hor[[i]]<-TZ_4
}

m_hor<-do.call(bind,TZ_hor)

# take mean of all ratio maps
NS_ratio_table_hor<-matrix(unlist(NS_ratio_hor),ncol=niter)
NS_logratio_table_hor<-matrix(unlist(NS_logratio_hor),ncol=niter)
NS_ratiotrans_table_hor<-matrix(unlist(NS_ratiotrans_hor),ncol=niter)

NS_ratio_mean_hor<-rowMeans(NS_ratio_table_hor)
NS_ratio_sd_hor<-apply(NS_ratio_table_hor,1,sd)

NS_logratio_mean_hor<-rowMeans(NS_logratio_table_hor)
NS_logratio_sd_hor<-apply(NS_logratio_table_hor,1,sd)

NS_ratiotrans_mean_hor<-rowMeans(NS_ratiotrans_table_hor)
NS_ratiotrans_sd_hor<-apply(NS_ratiotrans_table_hor,1,sd)

TZratio0_hor<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_hor)<-"North"
TZratio0_hor$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_hor$mean_ratio<-NS_ratio_mean_hor
TZratio0_hor$sd_ratio<-NS_ratio_sd_hor

TZratio0_hor$mean_logratio<-NS_logratio_mean_hor
TZratio0_hor$sd_logratio<-NS_logratio_sd_hor

TZratio0_hor$mean_ratiotrans<-NS_ratiotrans_mean_hor
TZratio0_hor$sd_ratiotrans<-NS_ratiotrans_sd_hor

# make spatial pixel data frame
TZratio_mean_hor<-as(TZratio0_hor,"SpatialPixelsDataFrame")

# map ratio up to 100
brks<-c(1,seq(10,100,length=10),5e17)
iterline<-list("sp.lines",m_hor,a=0.2)
spplot(TZratio_mean_hor,"mean_ratiotrans", sp.layout=iterline,main="50 samples ISA mean ratio, horizontal line",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))

########## 50k North line ###########
niter<-20
IndSppList_50kN<-vector("list",niter)
NS_ratio_50kN <- vector("list",niter)
NS_logratio_50kN <- vector("list",niter)
NS_ratiotrans_50kN <- vector("list",niter)
TZ4_loess_50kN <- vector("list",niter)
TZ4_cutoff_50kN <- vector("list",niter)
TZ_50kN <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_N_50kN,Trees_count_S_50kN,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_50kN[[i]]<-IndSpp
  
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
  NS_ratio_50kN[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_50kN[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_50kN[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ######################### binwidth = 0.1, truncation point = 20000 #############################
  ############## span = 50 bins ###############
  TZratio4_nc1<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins<-hist(TZratio4_nc1$ratio_transform,breaks=seq(0,max(TZratio4_nc1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins<-data.frame(TZratio4_ratio_halfbins$mids)
  TZ4_ratio_halfbins$counts<-TZratio4_ratio_halfbins$counts
  colnames(TZ4_ratio_halfbins)<-c("bin","counts")
  TZ4_ratio_halfbins<-slice(TZ4_ratio_halfbins,-c(1:10))
  
  # use loess to fit a line
  #TZ4_halfbin_fit <- loess(TZ4_ratio_halfbins$counts~TZ4_ratio_halfbins$bin)
  #plot(TZ4_ratio_halfbins$bin[1:100], TZ4_ratio_halfbins$counts[1:100])
  #lines(TZ4_ratio_halfbins$bin[1:100],TZ4_halfbin_fit$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit <- loess(TZ4_ratio_halfbins$counts~TZ4_ratio_halfbins$bin,span = (50/nrow(TZ4_ratio_halfbins)))
  plot(TZ4_ratio_halfbins$bin[1:100], TZ4_ratio_halfbins$counts[1:100])
  lines(TZ4_ratio_halfbins$bin[1:100],TZ4_halfbin_fit$fitted[1:100])
  
  # local slope
  TZ4_loess_slope <- data.frame(diff(TZ4_halfbin_fit$fitted)/diff(TZ4_ratio_halfbins$bin)) %>% bind_cols(TZ4_ratio_halfbins[2:nrow(TZ4_ratio_halfbins),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_50kN[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_50kN[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_50kN[[i]]<-TZ_4
}

m_50kN<-do.call(bind,TZ_50kN)

# take mean of all ratio maps
NS_ratio_50kN_table<-matrix(unlist(NS_ratio_50kN),ncol=niter)
NS_logratio_50kN_table<-matrix(unlist(NS_logratio_50kN),ncol=niter)
NS_ratiotrans_50kN_table<-matrix(unlist(NS_ratiotrans_50kN),ncol=niter)

NS_ratio_50kN_mean<-rowMeans(NS_ratio_50kN_table)
NS_ratio_50kN_sd<-apply(NS_ratio_50kN_table,1,sd)

NS_logratio_50kN_mean<-rowMeans(NS_logratio_50kN_table)
NS_logratio_50kN_sd<-apply(NS_logratio_50kN_table,1,sd)

NS_ratiotrans_50kN_mean<-rowMeans(NS_ratiotrans_50kN_table)
NS_ratiotrans_50kN_sd<-apply(NS_ratiotrans_50kN_table,1,sd)

TZratio0_50kN<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_50kN)<-"North"
TZratio0_50kN$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_50kN$mean_ratio<-NS_ratio_50kN_mean
TZratio0_50kN$sd_ratio<-NS_ratio_50kN_sd

TZratio0_50kN$mean_logratio<-NS_logratio_50kN_mean
TZratio0_50kN$sd_logratio<-NS_logratio_50kN_sd

TZratio0_50kN$mean_ratiotrans<-NS_ratiotrans_50kN_mean
TZratio0_50kN$sd_ratiotrans<-NS_ratiotrans_50kN_sd

# make spatial pixel data frame
TZratio_mean_50kN<-as(TZratio0_50kN,"SpatialPixelsDataFrame")

# map ratio up to 100
brks<-c(1,seq(10,100,length=10),5e17)
iterline_50kN<-list("sp.lines",m_50kN,a=0.2)
spplot(TZratio_mean_50kN,"mean_ratiotrans", sp.layout=iterline_50kN,main="50 samples ISA mean ratio, horizontal line shifted 50 km north",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))

## average TZ line
avg_cutoff_50kN<-mean(unlist(TZ4_cutoff_50kN))
sd_cutoff_50kN<-sd(unlist(TZ4_cutoff_50kN))

# map of average cutoff on average ratio map
TZratio4_mean_image_50kN<-as.image.SpatialGridDataFrame(TZratio_mean_50kN["mean_ratiotrans"])
TZcontour_mean_4_50kN<-contourLines(TZratio4_mean_image_50kN,levels=avg_cutoff_50kN)
TZ_mean_4_50kN<-ContourLines2SLDF(TZcontour_mean_4_50kN)

brks<-c(1,seq(10,100,length=10),5e17)
spplot(TZratio_mean_50kN,"mean_ratiotrans", sp.layout=TZ_mean_4_50kN,main="50 samples ISA, mean cutoff at mean ratio, horizontal line shifted 50 km north",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))



####### 50k South line #########
niter<-20
IndSppList_50kS<-vector("list",niter)
NS_ratio_50kS <- vector("list",niter)
NS_logratio_50kS <- vector("list",niter)
NS_ratiotrans_50kS <- vector("list",niter)
TZ4_loess_50kS <- vector("list",niter)
TZ4_cutoff_50kS <- vector("list",niter)
TZ_50kS <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_N_50kS,Trees_count_S_50kS,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_50kS[[i]]<-IndSpp
  
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
  NS_ratio_50kS[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_50kS[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_50kS[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_50kS[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_50kS[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_50kS[[i]]<-TZ_4
}

m_50kS<-do.call(bind,TZ_50kS)

# take mean of all ratio maps
NS_ratio_50kS_table<-matrix(unlist(NS_ratio_50kS),ncol=niter)
NS_logratio_50kS_table<-matrix(unlist(NS_logratio_50kS),ncol=niter)
NS_ratiotrans_50kS_table<-matrix(unlist(NS_ratiotrans_50kS),ncol=niter)

NS_ratio_50kS_mean<-rowMeans(NS_ratio_50kS_table)
NS_ratio_50kS_sd<-apply(NS_ratio_50kS_table,1,sd)

NS_logratio_50kS_mean<-rowMeans(NS_logratio_50kS_table)
NS_logratio_50kS_sd<-apply(NS_logratio_50kS_table,1,sd)

NS_ratiotrans_50kS_mean<-rowMeans(NS_ratiotrans_50kS_table)
NS_ratiotrans_50kS_sd<-apply(NS_ratiotrans_50kS_table,1,sd)

TZratio0_50kS<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_50kS)<-"North"
TZratio0_50kS$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_50kS$mean_ratio<-NS_ratio_50kS_mean
TZratio0_50kS$sd_ratio<-NS_ratio_50kS_sd

TZratio0_50kS$mean_logratio<-NS_logratio_50kS_mean
TZratio0_50kS$sd_logratio<-NS_logratio_50kS_sd

TZratio0_50kS$mean_ratiotrans<-NS_ratiotrans_50kS_mean
TZratio0_50kS$sd_ratiotrans<-NS_ratiotrans_50kS_sd

# make spatial pixel data frame
TZratio_mean_50kS<-as(TZratio0_50kS,"SpatialPixelsDataFrame")

# map ratio up to 100
brks<-c(1,seq(10,100,length=10),5e17)
iterline_50kS<-list("sp.lines",m_50kS,a=0.2)
spplot(TZratio_mean_50kS,"mean_ratiotrans", sp.layout=iterline_50kS,main="50 samples ISA mean ratio, horizontal line shifted 50 km south",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))

## average TZ line
avg_cutoff_50kS<-mean(unlist(TZ4_cutoff_50kS))
sd_cutoff_50kS<-sd(unlist(TZ4_cutoff_50kS))

# map of average cutoff on average ratio map
TZratio4_mean_image_50kS<-as.image.SpatialGridDataFrame(TZratio_mean_50kS["mean_ratiotrans"])
TZcontour_mean_4_50kS<-contourLines(TZratio4_mean_image_50kS,levels=avg_cutoff_50kS)
TZ_mean_4_50kS<-ContourLines2SLDF(TZcontour_mean_4_50kS)

brks<-c(1,seq(10,100,length=10),5e17)
spplot(TZratio_mean_50kS,"mean_ratiotrans", sp.layout=TZ_mean_4_50kS,main="50 samples ISA, mean cutoff at mean ratio, FW ecoregion boundary shifted 50 km south",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))


######### 45 deg line #########
niter<-20
IndSppList_45deg<-vector("list",niter)
NS_ratio_45deg <- vector("list",niter)
NS_logratio_45deg <- vector("list",niter)
NS_ratiotrans_45deg <- vector("list",niter)
TZ4_loess_45deg <- vector("list",niter)
TZ4_cutoff_45deg <- vector("list",niter)
TZ_45deg <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_N_45deg,Trees_count_S_45deg,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_45deg[[i]]<-IndSpp
  
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
  NS_ratio_45deg[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_45deg[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_45deg[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_45deg[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_45deg[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_45deg[[i]]<-TZ_4
}

m_45deg<-do.call(bind,TZ_45deg)

# take mean of all ratio maps
NS_ratio_45deg_table<-matrix(unlist(NS_ratio_45deg),ncol=niter)
NS_logratio_45deg_table<-matrix(unlist(NS_logratio_45deg),ncol=niter)
NS_ratiotrans_45deg_table<-matrix(unlist(NS_ratiotrans_45deg),ncol=niter)

NS_ratio_45deg_mean<-rowMeans(NS_ratio_45deg_table)
NS_ratio_45deg_sd<-apply(NS_ratio_45deg_table,1,sd)

NS_logratio_45deg_mean<-rowMeans(NS_logratio_45deg_table)
NS_logratio_45deg_sd<-apply(NS_logratio_45deg_table,1,sd)

NS_ratiotrans_45deg_mean<-rowMeans(NS_ratiotrans_45deg_table)
NS_ratiotrans_45deg_sd<-apply(NS_ratiotrans_45deg_table,1,sd)

TZratio0_45deg<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_45deg)<-"North"
TZratio0_45deg$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_45deg$mean_ratio<-NS_ratio_45deg_mean
TZratio0_45deg$sd_ratio<-NS_ratio_45deg_sd

TZratio0_45deg$mean_logratio<-NS_logratio_45deg_mean
TZratio0_45deg$sd_logratio<-NS_logratio_45deg_sd

TZratio0_45deg$mean_ratiotrans<-NS_ratiotrans_45deg_mean
TZratio0_45deg$sd_ratiotrans<-NS_ratiotrans_45deg_sd

# make spatial pixel data frame
TZratio_mean_45deg<-as(TZratio0_45deg,"SpatialPixelsDataFrame")

# map ratio up to 100
brks<-c(1,seq(10,100,length=10),5e17)
iterline_45deg<-list("sp.lines",m_45deg,a=0.2)
spplot(TZratio_mean_45deg,"mean_ratiotrans", sp.layout=iterline_45deg,main="50 samples ISA mean ratio, horizontal line rotated 45 degrees",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))

## average TZ line
avg_cutoff_45deg<-mean(unlist(TZ4_cutoff_45deg))
sd_cutoff_45deg<-sd(unlist(TZ4_cutoff_45deg))

# map of average cutoff on average ratio map
TZratio4_mean_image_45deg<-as.image.SpatialGridDataFrame(TZratio_mean_45deg["mean_ratiotrans"])
TZcontour_mean_4_45deg<-contourLines(TZratio4_mean_image_45deg,levels=avg_cutoff_45deg)
TZ_mean_4_45deg<-ContourLines2SLDF(TZcontour_mean_4_45deg)

brks<-c(1,seq(10,100,length=10),5e17)
spplot(TZratio_mean_45deg,"mean_ratiotrans", sp.layout=TZ_mean_4_45deg,main="50 samples ISA, mean cutoff at mean ratio, horizontal line rotated 45 degrees",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))


############ 135 deg line ############
niter<-20
IndSppList_135deg<-vector("list",niter)
NS_ratio_135deg <- vector("list",niter)
NS_logratio_135deg <- vector("list",niter)
NS_ratiotrans_135deg <- vector("list",niter)
TZ4_loess_135deg <- vector("list",niter)
TZ4_cutoff_135deg <- vector("list",niter)
TZ_135deg <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_N_135deg,Trees_count_S_135deg,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_135deg[[i]]<-IndSpp
  
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
  NS_ratio_135deg[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_135deg[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_135deg[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_135deg[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_135deg[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_135deg[[i]]<-TZ_4
}

m_135deg<-do.call(bind,TZ_135deg)

# take mean of all ratio maps
NS_ratio_135deg_table<-matrix(unlist(NS_ratio_135deg),ncol=niter)
NS_logratio_135deg_table<-matrix(unlist(NS_logratio_135deg),ncol=niter)
NS_ratiotrans_135deg_table<-matrix(unlist(NS_ratiotrans_135deg),ncol=niter)

NS_ratio_135deg_mean<-rowMeans(NS_ratio_135deg_table)
NS_ratio_135deg_sd<-apply(NS_ratio_135deg_table,1,sd)

NS_logratio_135deg_mean<-rowMeans(NS_logratio_135deg_table)
NS_logratio_135deg_sd<-apply(NS_logratio_135deg_table,1,sd)

NS_ratiotrans_135deg_mean<-rowMeans(NS_ratiotrans_135deg_table)
NS_ratiotrans_135deg_sd<-apply(NS_ratiotrans_135deg_table,1,sd)

TZratio0_135deg<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_135deg)<-"North"
TZratio0_135deg$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_135deg$mean_ratio<-NS_ratio_135deg_mean
TZratio0_135deg$sd_ratio<-NS_ratio_135deg_sd

TZratio0_135deg$mean_logratio<-NS_logratio_135deg_mean
TZratio0_135deg$sd_logratio<-NS_logratio_135deg_sd

TZratio0_135deg$mean_ratiotrans<-NS_ratiotrans_135deg_mean
TZratio0_135deg$sd_ratiotrans<-NS_ratiotrans_135deg_sd

# make spatial pixel data frame
TZratio_mean_135deg<-as(TZratio0_135deg,"SpatialPixelsDataFrame")

# map ratio up to 100
brks<-c(1,seq(10,100,length=10),5e17)
iterline_135deg<-list("sp.lines",m_135deg,a=0.2)
spplot(TZratio_mean_135deg,"mean_ratiotrans", sp.layout=iterline_135deg,main="50 samples ISA mean ratio, horizontal line rotated 135 degrees",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))

## average TZ line
avg_cutoff_135deg<-mean(unlist(TZ4_cutoff_135deg))
sd_cutoff_135deg<-sd(unlist(TZ4_cutoff_135deg))

# map of average cutoff on average ratio map
TZratio4_mean_image_135deg<-as.image.SpatialGridDataFrame(TZratio_mean_135deg["mean_ratiotrans"])
TZcontour_mean_4_135deg<-contourLines(TZratio4_mean_image_135deg,levels=avg_cutoff_135deg)
TZ_mean_4_135deg<-ContourLines2SLDF(TZcontour_mean_4_135deg)

brks<-c(1,seq(10,100,length=10),5e17)
spplot(TZratio_mean_135deg,"mean_ratiotrans", sp.layout=TZ_mean_4_135deg,main="50 samples ISA, mean cutoff at mean ratio, horizontal line rotated 135 degrees",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))


########### FS line ##############
niter<-20
IndSppList_FS<-vector("list",niter)
NS_ratio_FS <- vector("list",niter)
NS_logratio_FS <- vector("list",niter)
NS_ratiotrans_FS <- vector("list",niter)
TZ4_loess_FS <- vector("list",niter)
TZ4_cutoff_FS <- vector("list",niter)
TZ_FS <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_FS_N,Trees_count_FS_S,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_FS[[i]]<-IndSpp
  
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
  NS_ratio_FS[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_FS[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_FS[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_FS[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_FS[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_FS[[i]]<-TZ_4
}

m_FS<-do.call(bind,TZ_FS)

# take mean of all ratio maps
NS_ratio_FS_table<-matrix(unlist(NS_ratio_FS),ncol=niter)
NS_logratio_FS_table<-matrix(unlist(NS_logratio_FS),ncol=niter)
NS_ratiotrans_FS_table<-matrix(unlist(NS_ratiotrans_FS),ncol=niter)

NS_ratio_FS_mean<-rowMeans(NS_ratio_FS_table)
NS_ratio_FS_sd<-apply(NS_ratio_FS_table,1,sd)

NS_logratio_FS_mean<-rowMeans(NS_logratio_FS_table)
NS_logratio_FS_sd<-apply(NS_logratio_FS_table,1,sd)

NS_ratiotrans_FS_mean<-rowMeans(NS_ratiotrans_FS_table)
NS_ratiotrans_FS_sd<-apply(NS_ratiotrans_FS_table,1,sd)

TZratio0_FS<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_FS)<-"North"
TZratio0_FS$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_FS$mean_ratio<-NS_ratio_FS_mean
TZratio0_FS$sd_ratio<-NS_ratio_FS_sd

TZratio0_FS$mean_logratio<-NS_logratio_FS_mean
TZratio0_FS$sd_logratio<-NS_logratio_FS_sd

TZratio0_FS$mean_ratiotrans<-NS_ratiotrans_FS_mean
TZratio0_FS$sd_ratiotrans<-NS_ratiotrans_FS_sd

# make spatial pixel data frame
TZratio_mean_FS<-as(TZratio0_FS,"SpatialPixelsDataFrame")

# map ratio up to 100
brks<-c(1,seq(10,100,length=10),5e17)
iterline_FS<-list("sp.lines",m_FS,a=0.2)
spplot(TZratio_mean_FS,"mean_ratiotrans", sp.layout=iterline_FS,main="50 samples ISA mean ratio, FS ecoregion boundary",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))

## average TZ line
avg_cutoff_FS<-mean(unlist(TZ4_cutoff_FS))
sd_cutoff_FS<-sd(unlist(TZ4_cutoff_FS))

# map of average cutoff on average ratio map
TZratio4_mean_image_FS<-as.image.SpatialGridDataFrame(TZratio_mean_FS["mean_ratiotrans"])
TZcontour_mean_4_FS<-contourLines(TZratio4_mean_image_FS,levels=avg_cutoff_FS)
TZ_mean_4_FS<-ContourLines2SLDF(TZcontour_mean_4_FS)

brks<-c(1,seq(10,100,length=10),5e17)
spplot(TZratio_mean_FS,"mean_ratiotrans", sp.layout=TZ_mean_4_FS,main="50 samples ISA, mean cutoff at mean ratio, FS ecoregion boundary",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))


########### FS 50 km N line ##############
niter<-20
IndSppList_FS_50kN<-vector("list",niter)
NS_ratio_FS_50kN <- vector("list",niter)
NS_logratio_FS_50kN <- vector("list",niter)
NS_ratiotrans_FS_50kN <- vector("list",niter)
TZ4_loess_FS_50kN <- vector("list",niter)
TZ4_cutoff_FS_50kN <- vector("list",niter)
TZ_FS_50kN <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_FS_N_50kN,Trees_count_FS_S_50kN,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_FS_50kN[[i]]<-IndSpp
  
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
  NS_ratio_FS_50kN[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_FS_50kN[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_FS_50kN[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_FS_50kN[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_FS_50kN[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_FS_50kN[[i]]<-TZ_4
}

m_FS_50kN<-do.call(bind,TZ_FS_50kN)

# take mean of all ratio maps
NS_ratio_FS_50kN_table<-matrix(unlist(NS_ratio_FS_50kN),ncol=niter)
NS_logratio_FS_50kN_table<-matrix(unlist(NS_logratio_FS_50kN),ncol=niter)
NS_ratiotrans_FS_50kN_table<-matrix(unlist(NS_ratiotrans_FS_50kN),ncol=niter)

NS_ratio_FS_50kN_mean<-rowMeans(NS_ratio_FS_50kN_table)
NS_ratio_FS_50kN_sd<-apply(NS_ratio_FS_50kN_table,1,sd)

NS_logratio_FS_50kN_mean<-rowMeans(NS_logratio_FS_50kN_table)
NS_logratio_FS_50kN_sd<-apply(NS_logratio_FS_50kN_table,1,sd)

NS_ratiotrans_FS_50kN_mean<-rowMeans(NS_ratiotrans_FS_50kN_table)
NS_ratiotrans_FS_50kN_sd<-apply(NS_ratiotrans_FS_50kN_table,1,sd)

TZratio0_FS_50kN<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_FS_50kN)<-"North"
TZratio0_FS_50kN$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_FS_50kN$mean_ratio<-NS_ratio_FS_50kN_mean
TZratio0_FS_50kN$sd_ratio<-NS_ratio_FS_50kN_sd

TZratio0_FS_50kN$mean_logratio<-NS_logratio_FS_50kN_mean
TZratio0_FS_50kN$sd_logratio<-NS_logratio_FS_50kN_sd

TZratio0_FS_50kN$mean_ratiotrans<-NS_ratiotrans_FS_50kN_mean
TZratio0_FS_50kN$sd_ratiotrans<-NS_ratiotrans_FS_50kN_sd

# make spatial pixel data frame
TZratio_mean_FS_50kN<-as(TZratio0_FS_50kN,"SpatialPixelsDataFrame")

# map ratio up to 100
brks<-c(1,seq(10,100,length=10),5e17)
iterline_FS_50kN<-list("sp.lines",m_FS_50kN,a=0.2)
spplot(TZratio_mean_FS_50kN,"mean_ratiotrans", sp.layout=iterline_FS_50kN,main="50 samples ISA mean ratio, FS ecoregion boundary shifted 50 km north",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))

## average TZ line
avg_cutoff_FS_50kN<-mean(unlist(TZ4_cutoff_FS_50kN))
sd_cutoff_FS_50kN<-sd(unlist(TZ4_cutoff_FS_50kN))

# map of average cutoff on average ratio map
TZratio4_mean_image_FS_50kN<-as.image.SpatialGridDataFrame(TZratio_mean_FS_50kN["mean_ratiotrans"])
TZcontour_mean_4_FS_50kN<-contourLines(TZratio4_mean_image_FS_50kN,levels=avg_cutoff_FS_50kN)
TZ_mean_4_FS_50kN<-ContourLines2SLDF(TZcontour_mean_4_FS_50kN)

brks<-c(1,seq(10,100,length=10),5e17)
spplot(TZratio_mean_FS_50kN,"mean_ratiotrans", sp.layout=TZ_mean_4_FS_50kN,main="50 samples ISA, mean cutoff at mean ratio, FS ecoregion boundary shifted 50 km north",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))


########## FS 50k South line ##########
niter<-20
IndSppList_FS_50kS<-vector("list",niter)
NS_ratio_FS_50kS <- vector("list",niter)
NS_logratio_FS_50kS <- vector("list",niter)
NS_ratiotrans_FS_50kS <- vector("list",niter)
TZ4_loess_FS_50kS <- vector("list",niter)
TZ4_cutoff_FS_50kS <- vector("list",niter)
TZ_FS_50kS <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_FS_N_50kS,Trees_count_FS_S_50kS,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_FS_50kS[[i]]<-IndSpp
  
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
  NS_ratio_FS_50kS[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_FS_50kS[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_FS_50kS[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_FS_50kS[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_FS_50kS[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_FS_50kS[[i]]<-TZ_4
}

m_FS_50kS<-do.call(bind,TZ_FS_50kS)

# take mean of all ratio maps
NS_ratio_FS_50kS_table<-matrix(unlist(NS_ratio_FS_50kS),ncol=niter)
NS_logratio_FS_50kS_table<-matrix(unlist(NS_logratio_FS_50kS),ncol=niter)
NS_ratiotrans_FS_50kS_table<-matrix(unlist(NS_ratiotrans_FS_50kS),ncol=niter)

NS_ratio_FS_50kS_mean<-rowMeans(NS_ratio_FS_50kS_table)
NS_ratio_FS_50kS_sd<-apply(NS_ratio_FS_50kS_table,1,sd)

NS_logratio_FS_50kS_mean<-rowMeans(NS_logratio_FS_50kS_table)
NS_logratio_FS_50kS_sd<-apply(NS_logratio_FS_50kS_table,1,sd)

NS_ratiotrans_FS_50kS_mean<-rowMeans(NS_ratiotrans_FS_50kS_table)
NS_ratiotrans_FS_50kS_sd<-apply(NS_ratiotrans_FS_50kS_table,1,sd)

TZratio0_FS_50kS<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_FS_50kS)<-"North"
TZratio0_FS_50kS$South<-as(South.dens4,"SpatialGridDataFrame")$v

TZratio0_FS_50kS$mean_ratio<-NS_ratio_FS_50kS_mean
TZratio0_FS_50kS$sd_ratio<-NS_ratio_FS_50kS_sd

TZratio0_FS_50kS$mean_logratio<-NS_logratio_FS_50kS_mean
TZratio0_FS_50kS$sd_logratio<-NS_logratio_FS_50kS_sd

TZratio0_FS_50kS$mean_ratiotrans<-NS_ratiotrans_FS_50kS_mean
TZratio0_FS_50kS$sd_ratiotrans<-NS_ratiotrans_FS_50kS_sd

# make spatial pixel data frame
TZratio_mean_FS_50kS<-as(TZratio0_FS_50kS,"SpatialPixelsDataFrame")

# map ratio up to 100
brks<-c(1,seq(10,100,length=10),5e17)
iterline_FS_50kS<-list("sp.lines",m_FS_50kS,a=0.2)
spplot(TZratio_mean_FS_50kS,"mean_ratiotrans", sp.layout=iterline_FS_50kS,main="50 samples ISA mean ratio, FS ecoregion boundary shifted 50 km south",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))

## average TZ line
avg_cutoff_FS_50kS<-mean(unlist(TZ4_cutoff_FS_50kS))
sd_cutoff_FS_50kS<-sd(unlist(TZ4_cutoff_FS_50kS))

# map of average cutoff on average ratio map
TZratio4_mean_image_FS_50kS<-as.image.SpatialGridDataFrame(TZratio_mean_FS_50kS["mean_ratiotrans"])
TZcontour_mean_4_FS_50kS<-contourLines(TZratio4_mean_image_FS_50kS,levels=avg_cutoff_FS_50kS)
TZ_mean_4_FS_50kS<-ContourLines2SLDF(TZcontour_mean_4_FS_50kS)

brks<-c(1,seq(10,100,length=10),5e17)
spplot(TZratio_mean_FS_50kS,"mean_ratiotrans", sp.layout=TZ_mean_4_FS_50kS,main="50 samples ISA, mean cutoff at mean ratio, FS ecoregion boundary shifted 50 km south",at=brks,col.regions=rev(brewer.pal(11,"RdYlBu")))


#### plots of boundaries ####
plot(wisc_HARN,border="grey50",lwd=2,main="FS line position")
plot(m_FS,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(Zones,border="blue",add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2,main="50 km north FS line position")
plot(m_FS_50kN,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(Zones_50kN,border="blue",add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2,main="50km S FS line position")
plot(m_FS_50kS,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(Zones_50kS,border="blue",add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2,main="Horizontal line position")
plot(m_hor,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(wisc_split,border="blue",add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2,main="50 km north horizontal line position")
plot(m_50kN,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(wisc_split_50kN,border="blue",add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2,main="50 km south horizontal line position")
plot(m_50kS,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(wisc_split_50kS,border="blue",add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2,main="45 deg rotation line position")
plot(m_45deg,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(wisc_split_45deg,border="blue",add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2,main="135 deg rotation line position")
plot(m_135deg,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(wisc_split_135deg,border="blue",add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2,main="Overlay TZ line from different straight line starting positions")
plot(m_hor,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(m_50kN,col=rgb(0,0,1,alpha=0.2),add=TRUE)
plot(m_50kS,col=rgb(1,0,0,alpha=0.2),add=TRUE)
plot(m_45deg,col=rgb(0,1,1,alpha=0.2),add=TRUE)
plot(m_135deg,col=rgb(1,0,1,alpha=0.2),add=TRUE)

plot(wisc_HARN,border="grey50",lwd=2)
plot(m_FS,col=rgb(0,0,0,alpha=0.3),add=TRUE)
plot(m_FS_50kN,col=rgb(0,0,1,alpha=0.2),add=TRUE)
plot(m_FS_50kS,col=rgb(1,0,0,alpha=0.2),add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

plot(wisc_HARN,border="grey50",lwd=2,main="Overlay TZ line from different FS boundary starting positions")
plot(TZ_mean_4_FS,col=rgb(0,0,0,alpha=0.8),add=TRUE)
plot(TZ_mean_4_FS_50kN,col=rgb(0,0,1,alpha=0.7),add=TRUE)
plot(TZ_mean_4_FS_50kS,col=rgb(1,0,0,alpha=0.7),add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')


#### save figures
theme_set(theme_bw())

### convert to sf
wisc_HARN_sf <- as(wisc_HARN,"sf")

## combine lines
# add field with variable
m_FS$init <- "FS"
m_FS_50kN$init <- "FS_50kN"
m_FS_50kS$init <- "FS_50kS"
m_hor$init <- "hor"
m_45deg$init <- "45deg"
m_135deg$init <- "135deg"

# combine FS lines
m_FS_all <- rbind(m_FS, m_FS_50kN, m_FS_50kS)
crs(m_FS_all) <- "+init=epsg:3070"
m_FS_all_sf <- as(m_FS_all,"sf")

# combine straight lines
m_sl_all <- rbind(m_hor,m_45deg,m_135deg)
crs(m_sl_all) <- "+init=epsg:3070"
m_sl_all$init <- factor(m_sl_all$init, levels=c("hor","45deg","135deg"))
m_sl_all_sf <- as(m_sl_all,"sf")

## plots
# FS lines
Fig4_1 <- ggplot(data=wisc_HARN_sf)+
  geom_sf(fill="white")+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1, "cm"), #how tall the arrow should be
                         width= unit(0.5, "cm"), 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_orienteering)+
  geom_sf(data=m_FS_all_sf,aes(col=init,fill=init))+
  scale_color_manual(values=alpha(c("FS"="black",
                              "FS_50kN" = "blue",
                              "FS_50kS" = "red"),0.2),
                     labels=c("Original boundary","Shifted 50 km north","Shifted 50 km south"),guide="none")+
  scale_fill_manual(values=c("FS"="black",
                                    "FS_50kN" = "blue",
                                    "FS_50kS" = "red"),
                    labels=c("Original boundary","Shifted 50 km north","Shifted 50 km south"))+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")

# straight lines
Fig4_2 <- ggplot(data=wisc_HARN_sf)+
  geom_sf(fill="white")+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1, "cm"), #how tall the arrow should be
                         width= unit(0.5, "cm"), 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_orienteering)+
  geom_sf(data=m_sl_all_sf,aes(col=init,fill=init))+
  scale_color_manual(breaks=c("hor","45deg","135deg"),values=alpha(c("hor"="black",
                                    "45deg" = "cyan",
                                    "135deg" = "magenta"),0.2),
                     labels=c("Original boundary","Shifted 50 km north","Shifted 50 km south"),guide="none")+
  scale_fill_manual(values=c("hor"="black",
                             "45deg" = "cyan",
                             "135deg" = "magenta"),
                    labels=c("Horizontal line","NW to SE diagonal line","SW to NE diagonal line"))+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")

### combine initial boundary lines
# add field with variable
wisc_split$init <- "Horizontal"
wisc_split_135deg$init <- "Diagonal, 45 degrees"
wisc_split_45deg$init <- "Diagonal, 135 degrees"
Zones$init <- "FS Ecoregion Boundary"
Zones_50kN$init <- "FS Ecoregion Boundary, 50 km N"
Zones_50kS$init <- "FS Ecoregion Boundary, 50 km S"

# convert from line to polygon
wisc_split_line <- as(wisc_split,"SpatialLinesDataFrame")
wisc_split_135deg_line <- as(wisc_split_135deg,"SpatialLinesDataFrame")
wisc_split_45deg_line <- as(wisc_split_45deg,"SpatialLinesDataFrame")
Zones_line <- as(Zones,"SpatialLinesDataFrame")
Zones_50kN_line <- as(Zones_50kN,"SpatialLinesDataFrame")
Zones_50kS_line <- as(Zones_50kS,"SpatialLinesDataFrame")

## split and extract line crossing state
# overlay
wisc_split_line_short <- raster::intersect(wisc_split_line,wisc_poly)
wisc_split_135deg_line_short <- raster::intersect(wisc_split_135deg_line,wisc_poly)
wisc_split_45deg_line_short <- raster::intersect(wisc_split_45deg_line,wisc_poly)

Zones_line_short <- raster::intersect(Zones_line,wisc_HARN)
Zones_50kN_line_short <- raster::intersect(Zones_50kN_line,wisc_HARN)
Zones_50kS_line_short <- raster::intersect(Zones_50kS_line,wisc_HARN)

plot(wisc_HARN)  
plot(wisc_split_line_short,add=TRUE)
plot(wisc_split_135deg_line_short,col="magenta",add=TRUE)
plot(wisc_split_45deg_line_short,col="cyan",add=TRUE)

plot(wisc_HARN)
plot(Zones_line_short,add=TRUE)
plot(Zones_50kN_line_short,col="blue",add=TRUE)
plot(Zones_50kS_line_short,col="red",add=TRUE)
  
# FS boundary
Zones_sl <- rbind(Zones_line_short,Zones_50kN_line_short,Zones_50kS_line_short)
crs(Zones_sl) <- "+init=epsg:3070"
Zones_sl_sf <- as(Zones_sl,"sf")


Fig4_3 <- ggplot(data=wisc_HARN_sf)+
  geom_sf(fill="white")+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1, "cm"), #how tall the arrow should be
                         width= unit(0.5, "cm"), 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_orienteering)+
  geom_sf(data=Zones_sl_sf,aes(col=init,fill=init))+
  scale_color_manual(breaks=c("FS Ecoregion Boundary","FS Ecoregion Boundary, 50 km N","FS Ecoregion Boundary, 50 km S"),values=alpha(c("FS Ecoregion Boundary"="black",
                                                                                                          "FS Ecoregion Boundary, 50 km N" = "blue",
                                                                                                          "FS Ecoregion Boundary, 50 km S" = "red")),
                     labels=c("FS Ecoregion Boundary","FS Ecoregion Boundary, 50 km N","FS Ecoregion Boundary, 50 km S"),guide="none")+
  scale_fill_manual(values=c("FS Ecoregion Boundary"="black",
                             "FS Ecoregion Boundary, 50 km N" = "blue",
                             "FS Ecoregion Boundary, 50 km S" = "red"),
                    labels=c("Original boundary","Shifted 50 km north","Shifted 50 km south"))+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")


Fig4_3

# straight lines
wisc_split_sl <- rbind(wisc_split_line_short,wisc_split_45deg_line_short,wisc_split_135deg_line_short)
crs(wisc_split_sl) <- "+init=epsg:3070"
wisc_split_sl_sf <- as(wisc_split_sl,"sf")


Fig4_4 <- ggplot(data=wisc_HARN_sf)+
  geom_sf(fill="white")+
  annotation_scale(location = "bl", width_hint = 0.25) +
  geom_sf(data=wisc_split_sl_sf,aes(col=init,fill=init))+
  scale_color_manual(breaks=c("Horizontal","Diagonal, 45 degrees","Diagonal, 135 degrees"),
                     values=alpha(c("Horizontal"="black","Diagonal, 45 degrees" = "magenta","Diagonal, 135 degrees" = "cyan")),
                     labels=c("Horizontal","Diagonal, 45 degrees","Diagonal, 135 degrees"),guide="none")+
  scale_fill_manual(breaks=c("Horizontal","Diagonal, 45 degrees","Diagonal, 135 degrees"),
                              values=c("Horizontal"="black",
                             "Diagonal, 45 degrees" = "magenta",
                             "Diagonal, 135 degrees" = "cyan"),
                    labels=c("Horizontal line","SW to NE diagonal line","NW to SE diagonal line"))+
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())+
  theme(legend.title = element_blank(),
        legend.text=element_text(size=10),
        legend.position = "bottom",legend.direction = "vertical")


Fig4_4



#### save plots
dpi=600
ggsave("./Output/Fig4_1.jpg", plot=Fig4_1, width = 4000/dpi, height = 4000/dpi, dpi = dpi)
ggsave("./Output/Fig4_2.jpg", plot=Fig4_2, width = 4000/dpi, height = 4000/dpi, dpi = dpi)
ggsave("./Output/Fig4_3.jpg", plot=Fig4_3, width = 4000/dpi, height = 4000/dpi, dpi = dpi)
ggsave("./Output/Fig4_4.jpg", plot=Fig4_4, width = 4000/dpi, height = 4000/dpi, dpi = dpi)



