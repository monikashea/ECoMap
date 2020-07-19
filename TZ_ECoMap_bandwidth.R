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

## bring in WI map to set boundaries for density mapping
wisc<-readOGR(dsn="./WI_HARN_mask",layer="WI_HARN_mask")
wisc_HARN <- spTransform(wisc,CRS=CRS("+init=epsg:3070"))

## import 6 mi grid for ISA
Grid_6mi <- readOGR(dsn="./6mi_Grid",layer="6mi_Grid")
Grid_6mi <- spTransform(Grid_6mi,CRS=CRS("+init=epsg:3070"))

## determine which trees belong to each grid cell
pg_overlay <- over(ptree_diff_2018,Grid_6mi)
pg_overlay_trees <- cbind(pg_overlay,ptree_diff_2018@data)
pg_overlay_trees <- pg_overlay_trees %>% rename(Grid_6mi=Grid_ID)

## data on area of each section of grid cell x ecoregion boundary intersection (done in ArcGIS, could also be done in R, just didn't have the skill at the time)
Grid_zone<-read.csv(file="6mi_Grid_ecoregions_area.csv",header=TRUE)
Grid_zone<- Grid_zone %>% rename(Grid_6mi=Grid_ID) %>% mutate(ZoneFS=ifelse(Ecoregion==212,"North","South"))

# select grid ecoregions with largest area
Grid_ecoregion <- Grid_zone %>% group_by(Grid_6mi) %>% summarise(max=max(F_AREA)) %>% left_join(Grid_zone,by="Grid_6mi") %>% 
  filter(max==F_AREA) %>% select(Grid_6mi,ZoneFS)

### identify grid cells on edges of state boundary (or in lakes)
# bring in zones spatial
Zones <- readOGR("./212_222_boundary","212_222_merge_extended")
Zones <- spTransform(Zones,CRS=CRS("+init=epsg:3070"))

# convert grid to raster
r_6mi<-raster(ncol=51,nrow=54)
extent(r_6mi)<-extent(Grid_6mi)
GridRaster_6mi<-rasterize(Grid_6mi,r_6mi,'Grid_ID')
crs(GridRaster_6mi) <- "+init=epsg:3070"

# filter cells than fall entirely inside state boundary (excluding edge cells)
Grid_WI_6mi <- raster::extract(GridRaster_6mi,wisc_HARN,df=TRUE,weights=TRUE)
Grid_WI_6mi<-Grid_WI_6mi %>% rename(Grid_6mi=layer) %>% filter(weight==max(weight))


### join tree data to zone data
Trees <- left_join(pg_overlay_trees,Grid_ecoregion,by="Grid_6mi")

# remove NA, UK, RR (NA, unknown species, and very rare species with <10 points)
Trees<-Trees %>% filter(!SP_new %in% c(' ','UK','RR')) %>% filter(!is.na(SP_new))
Trees$SP_new<-as.factor(Trees$SP_new)
Trees$SP_new<-droplevels(Trees$SP_new)

# convert to species counts for each cell, remove edge cells, and add ecoregion
Trees_count <- Trees %>% group_by(SP_new,Grid_6mi) %>% summarise(count=n()) %>% spread(SP_new,count) %>% 
  semi_join(Grid_WI_6mi,by="Grid_6mi") %>% left_join(Grid_ecoregion,by="Grid_6mi")
Trees_count[is.na(Trees_count)]<-0 # change NA to 0

## separate matrices for each ecoregion
# north
Trees_count_north<-filter(Trees_count,ZoneFS=="North")
row.names(Trees_count_north)<-Trees_count_north$Grid_6mi
Trees_count_north<-select(Trees_count_north,-c(Grid_6mi,ZoneFS))

# south
Trees_count_south<-filter(Trees_count,ZoneFS=="South")
row.names(Trees_count_south)<-Trees_count_south$Grid_6mi
Trees_count_south<-select(Trees_count_south,-c(Grid_6mi,ZoneFS))


##################### ISA and mapping #####################
## prepare data for ISA
north <- Trees_count_north %>% mutate(Grid=row.names(.))
south <- Trees_count_south %>% mutate(Grid=row.names(.))

# ISA function
ISA_PLS<-function(north,south,nsample){
  # randomly select nsample cells from each ecoregion
  Trees_N_1<-sample_n(north,nsample,replace=FALSE)
  Trees_S_1<-sample_n(south,nsample,replace=FALSE)
  Trees_NS_1<-rbind(Trees_N_1,Trees_S_1)
  groups = c(rep(1,nsample),rep(2,nsample))
  
  # ISA
  Trees_NS_indval_1 = multipatt(Trees_NS_1,groups,control=how(nperm=9999))
  
  ## summarize results
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

############## bandwidth = 2 km ################
#### mapping ####
### iterations
niter<-50
IndSppList_bw2<-vector("list",niter)
NS_ratio_bw2 <- vector("list",niter)
NS_logratio_bw2 <- vector("list",niter)
NS_ratiotrans_bw2 <- vector("list",niter)
TZ4_loess_bw2 <- vector("list",niter)
TZ4_cutoff_bw2 <- vector("list",niter)
TZ_bw2 <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north,Trees_count_south,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_bw2[[i]]<-IndSpp
  
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
  bw.tz4<-2000
  
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
  NS_ratio_bw2[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_bw2[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_bw2[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_bw2[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_bw2[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_bw2[[i]]<-TZ_4
}

m_bw2<-do.call(bind,TZ_bw2)

## Relative ratio of northern and southern trees in Wisconsin
# take mean of all ratio maps
NS_ratiotrans_table_bw2<-matrix(unlist(NS_ratiotrans_bw2),ncol=niter)
NS_ratiotrans_mean_bw2<-rowMeans(NS_ratiotrans_table_bw2)

# make new SGDF with mean ratiotrans
TZratio0_bw2<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_bw2)<-"North"
TZratio0_bw2$South<-as(South.dens4,"SpatialGridDataFrame")$v
TZratio0_bw2$mean_ratiotrans<-NS_ratiotrans_mean_bw2

# convert to spatial pixel data frame
TZratio_mean_bw2<-as(TZratio0_bw2,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e21) # set brakes so scale is readable
iterline<-list("sp.lines",m_bw2,a=0.2) # set formatting for the 50 lines

TZratio_mean_bw2@data$cut_mean <- cut(TZratio_mean_bw2@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
FigS4A <-spplot(TZratio_mean_bw2,"cut_mean", sp.layout=list(iterline,wisc_HARN),
       col.regions=rev(brewer.pal(11,"RdYlBu")),
       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
       par.settings = list(axis.line = list(col =  'transparent')),
       auto.key = list(title = "Relative ratio"),
       main="Bandwidth = 2 km")

FigS4A


## Average ratio line
avg_cutoff_bw2<-mean(unlist(TZ4_cutoff_bw2))
sd_cutoff_bw2<-sd(unlist(TZ4_cutoff_bw2))

# map of average cutoff on average ratio map
TZratio4_mean_image_bw2<-as.image.SpatialGridDataFrame(TZratio_mean_bw2["mean_ratiotrans"])
TZcontour_mean_4_bw2<-contourLines(TZratio4_mean_image_bw2,levels=avg_cutoff_bw2)
TZ_mean_4_bw2<-ContourLines2SLDF(TZcontour_mean_4_bw2)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e19) # set brakes so scale is readable
TZratio_mean_bw2@data$cut_mean <- cut(TZratio_mean_bw2@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S4A <-spplot(TZratio_mean_bw2,"cut_mean", sp.layout=list(TZ_mean_4_bw2,wisc_HARN),
              col.regions=rev(brewer.pal(11,"RdYlBu")),
              colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
              par.settings = list(axis.line = list(col =  'transparent')),
              auto.key = list(title = "Relative ratio"))

mean_plot_S4A


############## bandwidth = 3 km ################
#### mapping ####
### iterations
niter<-50
IndSppList_bw3<-vector("list",niter)
NS_ratio_bw3 <- vector("list",niter)
NS_logratio_bw3 <- vector("list",niter)
NS_ratiotrans_bw3 <- vector("list",niter)
TZ4_loess_bw3 <- vector("list",niter)
TZ4_cutoff_bw3 <- vector("list",niter)
TZ_bw3 <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north,Trees_count_south,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_bw3[[i]]<-IndSpp
  
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
  bw.tz4<-3000
  
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
  NS_ratio_bw3[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_bw3[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_bw3[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_bw3[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_bw3[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_bw3[[i]]<-TZ_4
}

m_bw3<-do.call(bind,TZ_bw3)

## Relative ratio of northern and southern trees in Wisconsin
# take mean of all ratio maps
NS_ratiotrans_table_bw3<-matrix(unlist(NS_ratiotrans_bw3),ncol=niter)
NS_ratiotrans_mean_bw3<-rowMeans(NS_ratiotrans_table_bw3)

# make new SGDF with mean ratiotrans
TZratio0_bw3<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_bw3)<-"North"
TZratio0_bw3$South<-as(South.dens4,"SpatialGridDataFrame")$v
TZratio0_bw3$mean_ratiotrans<-NS_ratiotrans_mean_bw3

# convert to spatial pixel data frame
TZratio_mean_bw3<-as(TZratio0_bw3,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_bw3,a=0.2) # set formatting for the 50 lines

TZratio_mean_bw3@data$cut_mean <- cut(TZratio_mean_bw3@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
FigS4B <-spplot(TZratio_mean_bw3,"cut_mean", sp.layout=list(iterline,wisc_HARN),
              col.regions=rev(brewer.pal(11,"RdYlBu")),
              colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
              par.settings = list(axis.line = list(col =  'transparent')),
              auto.key = list(title = "Relative ratio"),
              main="Bandwidth = 3 km")

FigS4B



## Average ratio line 
avg_cutoff_bw3<-mean(unlist(TZ4_cutoff_bw3))
sd_cutoff_bw3<-sd(unlist(TZ4_cutoff_bw3))

# map of average cutoff on average ratio map
TZratio4_mean_image_bw3<-as.image.SpatialGridDataFrame(TZratio_mean_bw3["mean_ratiotrans"])
TZcontour_mean_4_bw3<-contourLines(TZratio4_mean_image_bw3,levels=avg_cutoff_bw3)
TZ_mean_4_bw3<-ContourLines2SLDF(TZcontour_mean_4_bw3)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_bw3@data$cut_mean <- cut(TZratio_mean_bw3@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S4B <-spplot(TZratio_mean_bw3,"cut_mean", sp.layout=list(TZ_mean_4_bw3,wisc_HARN),
                   col.regions=rev(brewer.pal(11,"RdYlBu")),
                   colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                   par.settings = list(axis.line = list(col =  'transparent')),
                   auto.key = list(title = "Relative ratio"))

mean_plot_S4B


############## bandwidth = 4 km ################
#### mapping ####
### iterations
niter<-50
IndSppList_bw4<-vector("list",niter)
NS_ratio_bw4 <- vector("list",niter)
NS_logratio_bw4 <- vector("list",niter)
NS_ratiotrans_bw4 <- vector("list",niter)
TZ4_loess_bw4 <- vector("list",niter)
TZ4_cutoff_bw4 <- vector("list",niter)
TZ_bw4 <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north,Trees_count_south,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_bw4[[i]]<-IndSpp
  
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
  NS_ratio_bw4[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_bw4[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_bw4[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_bw4[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_bw4[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_bw4[[i]]<-TZ_4
}

m_bw4<-do.call(bind,TZ_bw4)

## Relative ratio of northern and southern trees in Wisconsin
# take mean of all ratio maps
NS_ratiotrans_table_bw4<-matrix(unlist(NS_ratiotrans_bw4),ncol=niter)
NS_ratiotrans_mean_bw4<-rowMeans(NS_ratiotrans_table_bw4)

# make new SGDF with mean ratiotrans
TZratio0_bw4<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_bw4)<-"North"
TZratio0_bw4$South<-as(South.dens4,"SpatialGridDataFrame")$v
TZratio0_bw4$mean_ratiotrans<-NS_ratiotrans_mean_bw4

# convert to spatial pixel data frame
TZratio_mean_bw4<-as(TZratio0_bw4,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_bw4,a=0.2) # set formatting for the 50 lines

TZratio_mean_bw4@data$cut_mean <- cut(TZratio_mean_bw4@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
FigS4C <-spplot(TZratio_mean_bw4,"cut_mean", sp.layout=list(iterline,wisc_HARN),
              col.regions=rev(brewer.pal(11,"RdYlBu")),
              colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
              par.settings = list(axis.line = list(col =  'transparent')),
              auto.key = list(title = "Relative ratio"),
              main="Bandwidth = 4 km")

FigS4C


## Average ratio line 
avg_cutoff_bw4<-mean(unlist(TZ4_cutoff_bw4))
sd_cutoff_bw4<-sd(unlist(TZ4_cutoff_bw4))

# map of average cutoff on average ratio map
TZratio4_mean_image_bw4<-as.image.SpatialGridDataFrame(TZratio_mean_bw4["mean_ratiotrans"])
TZcontour_mean_4_bw4<-contourLines(TZratio4_mean_image_bw4,levels=avg_cutoff_bw4)
TZ_mean_4_bw4<-ContourLines2SLDF(TZcontour_mean_4_bw4)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_bw4@data$cut_mean <- cut(TZratio_mean_bw4@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S4C <-spplot(TZratio_mean_bw4,"cut_mean", sp.layout=list(TZ_mean_4_bw4,wisc_HARN),
                   col.regions=rev(brewer.pal(11,"RdYlBu")),
                   colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                   par.settings = list(axis.line = list(col =  'transparent')),
                   auto.key = list(title = "Relative ratio"))

mean_plot_S4C


############## bandwidth = 5 km ################
#### mapping ####
### iterations
niter<-50
IndSppList_bw5<-vector("list",niter)
NS_ratio_bw5 <- vector("list",niter)
NS_logratio_bw5 <- vector("list",niter)
NS_ratiotrans_bw5 <- vector("list",niter)
TZ4_loess_bw5 <- vector("list",niter)
TZ4_cutoff_bw5 <- vector("list",niter)
TZ_bw5 <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north,Trees_count_south,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_bw5[[i]]<-IndSpp
  
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
  bw.tz4<-5000
  
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
  NS_ratio_bw5[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_bw5[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_bw5[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_bw5[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_bw5[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_bw5[[i]]<-TZ_4
}

m_bw5<-do.call(bind,TZ_bw5)

## Relative ratio of northern and southern trees in Wisconsin
# take mean of all ratio maps
NS_ratiotrans_table_bw5<-matrix(unlist(NS_ratiotrans_bw5),ncol=niter)
NS_ratiotrans_mean_bw5<-rowMeans(NS_ratiotrans_table_bw5)

# make new SGDF with mean ratiotrans
TZratio0_bw5<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_bw5)<-"North"
TZratio0_bw5$South<-as(South.dens4,"SpatialGridDataFrame")$v
TZratio0_bw5$mean_ratiotrans<-NS_ratiotrans_mean_bw5

# convert to spatial pixel data frame
TZratio_mean_bw5<-as(TZratio0_bw5,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e21) # set brakes so scale is readable
iterline<-list("sp.lines",m_bw5,a=0.2) # set formatting for the 50 lines

TZratio_mean_bw5@data$cut_mean <- cut(TZratio_mean_bw5@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
FigS4D <-spplot(TZratio_mean_bw5,"cut_mean", sp.layout=list(iterline,wisc_HARN),
              col.regions=rev(brewer.pal(11,"RdYlBu")),
              colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
              par.settings = list(axis.line = list(col =  'transparent')),
              auto.key = list(title = "Relative ratio"),
              main="Bandwidth = 5 km")

FigS4D


## Average ratio line 
avg_cutoff_bw5<-mean(unlist(TZ4_cutoff_bw5))
sd_cutoff_bw5<-sd(unlist(TZ4_cutoff_bw5))

# map of average cutoff on average ratio map
TZratio4_mean_image_bw5<-as.image.SpatialGridDataFrame(TZratio_mean_bw5["mean_ratiotrans"])
TZcontour_mean_4_bw5<-contourLines(TZratio4_mean_image_bw5,levels=avg_cutoff_bw5)
TZ_mean_4_bw5<-ContourLines2SLDF(TZcontour_mean_4_bw5)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_bw5@data$cut_mean <- cut(TZratio_mean_bw5@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S4D <-spplot(TZratio_mean_bw5,"cut_mean", sp.layout=list(TZ_mean_4_bw5,wisc_HARN),
                   col.regions=rev(brewer.pal(11,"RdYlBu")),
                   colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                   par.settings = list(axis.line = list(col =  'transparent')),
                   auto.key = list(title = "Relative ratio"))

mean_plot_S4D


############## bandwidth = 6 km ################
#### mapping ####
### iterations
niter<-50
IndSppList_bw6<-vector("list",niter)
NS_ratio_bw6 <- vector("list",niter)
NS_logratio_bw6 <- vector("list",niter)
NS_ratiotrans_bw6 <- vector("list",niter)
TZ4_loess_bw6 <- vector("list",niter)
TZ4_cutoff_bw6 <- vector("list",niter)
TZ_bw6 <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north,Trees_count_south,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList_bw6[[i]]<-IndSpp
  
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
  bw.tz4<-6000
  
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
  NS_ratio_bw6[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio_bw6[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans_bw6[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
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
  TZ4_loess_bw6[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_bw6[[i]]<-TZ4_halfbin_cutoff
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  TZ_bw6[[i]]<-TZ_4
}

m_bw6<-do.call(bind,TZ_bw6)

## Relative ratio of northern and southern trees in Wisconsin
# take mean of all ratio maps
NS_ratiotrans_table_bw6<-matrix(unlist(NS_ratiotrans_bw6),ncol=niter)
NS_ratiotrans_mean_bw6<-rowMeans(NS_ratiotrans_table_bw6)

# make new SGDF with mean ratiotrans
TZratio0_bw6<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0_bw6)<-"North"
TZratio0_bw6$South<-as(South.dens4,"SpatialGridDataFrame")$v
TZratio0_bw6$mean_ratiotrans<-NS_ratiotrans_mean_bw6

# convert to spatial pixel data frame
TZratio_mean_bw6<-as(TZratio0_bw6,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_bw6,a=0.2) # set formatting for the 50 lines

TZratio_mean_bw6@data$cut_mean <- cut(TZratio_mean_bw6@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
FigS4E <-spplot(TZratio_mean_bw6,"cut_mean", sp.layout=list(iterline,wisc_HARN),
              col.regions=rev(brewer.pal(11,"RdYlBu")),
              colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
              par.settings = list(axis.line = list(col =  'transparent')),
              auto.key = list(title = "Relative ratio"),
              main="Bandwidth = 6 km")

FigS4E


## Average ratio line 
avg_cutoff_bw6<-mean(unlist(TZ4_cutoff_bw6))
sd_cutoff_bw6<-sd(unlist(TZ4_cutoff_bw6))

# map of average cutoff on average ratio map
TZratio4_mean_image_bw6<-as.image.SpatialGridDataFrame(TZratio_mean_bw6["mean_ratiotrans"])
TZcontour_mean_4_bw6<-contourLines(TZratio4_mean_image_bw6,levels=avg_cutoff_bw6)
TZ_mean_4_bw6<-ContourLines2SLDF(TZcontour_mean_4_bw6)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean_bw6@data$cut_mean <- cut(TZratio_mean_bw6@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot_S4E <-spplot(TZratio_mean_bw6,"cut_mean", sp.layout=list(TZ_mean_4_bw6,wisc_HARN),
                   col.regions=rev(brewer.pal(11,"RdYlBu")),
                   colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                   par.settings = list(axis.line = list(col =  'transparent')),
                   auto.key = list(title = "Relative ratio"))

mean_plot_S4E


##### plot mean lines together #####
plot(wisc_HARN)
plot(TZ_mean_4_bw2, col="red",add=TRUE)
plot(TZ_mean_4_bw3,col="orange",add=TRUE)
plot(TZ_mean_4_bw4,col="green",add=TRUE)
plot(TZ_mean_4_bw5,col="blue",add=TRUE)
plot(TZ_mean_4_bw6,col="purple",add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

#### save image #####
### convert to sf
wisc_HARN_sf <- as(wisc_HARN,"sf")

## combine lines
# add field with variable
TZ_mean_4_bw2$bw <- "2km"
TZ_mean_4_bw3$bw <- "3km"
TZ_mean_4_bw4$bw <- "4km"
TZ_mean_4_bw5$bw <- "5km"
TZ_mean_4_bw6$bw <- "6km"

# combine FS lines
TZ_mean_4_bw_all <- rbind(TZ_mean_4_bw2,TZ_mean_4_bw3,TZ_mean_4_bw4,TZ_mean_4_bw5,TZ_mean_4_bw6)
crs(TZ_mean_4_bw_all) <- "+init=epsg:3070"
TZ_mean_4_bw_all$bw <- factor(TZ_mean_4_bw_all$bw,levels = c("2km","3km","4km","5km","6km"))
TZ_mean_4_bw_all_sf <- as(TZ_mean_4_bw_all,"sf")

### plot
FigS4 <- ggplot(data=wisc_HARN_sf)+
  geom_sf(fill="white")+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1, "cm"), #how tall the arrow should be
                         width= unit(0.5, "cm"), 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_orienteering)+
  geom_sf(data=TZ_mean_4_bw_all_sf,aes(col=bw,fill=bw))+
  scale_color_manual(breaks=c("2km","3km","4km","5km","6km"),
                     values=alpha(c("2km"="red","3km"="orange","4km"="green","5km"="blue","6km"="purple"),0.8),
                     labels=c("2 km","3 km","4 km","5 km","6 km"),guide="none")+
  scale_fill_manual(values=c("2km"="red","3km"="orange","4km"="green","5km"="blue","6km"="purple"),
                    labels=c("2 km","3 km","4 km","5 km","6 km"))+
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
ggsave("./Output/FigS4.jpg", plot=FigS4, width = 4000/dpi, height = 4000/dpi, dpi = dpi)






