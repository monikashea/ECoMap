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
library(TTR)


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

############################# testing different truncation points and bin sizes ##############################

#### mapping ####
### iterations
niter<-50
IndSppList<-vector("list",niter)
NS_ratio <- vector("list",niter)
NS_logratio <- vector("list",niter)
NS_ratiotrans <- vector("list",niter)
TZ4_loess_nc5 <- vector("list",niter)
TZ4_cutoff_nc5 <- vector("list",niter)
TZ_nc5 <- vector("list",niter)
TZ4_loess_co2 <- vector("list",niter)
TZ4_cutoff_co2 <- vector("list",niter)
TZ_co2 <- vector("list",niter)
TZ4_loess_co1 <- vector("list",niter)
TZ4_cutoff_co1 <- vector("list",niter)
TZ_co1 <- vector("list",niter)
TZ4_loess_co05 <- vector("list",niter)
TZ4_cutoff_co05 <- vector("list",niter)
TZ_co05 <- vector("list",niter)
TZ4_loess_nc1 <- vector("list",niter)
TZ4_cutoff_nc1 <- vector("list",niter)
TZ_nc1 <- vector("list",niter)
TZ4_loess_co2bs1 <- vector("list",niter)
TZ4_cutoff_co2bs1 <- vector("list",niter)
TZ_co2bs1 <- vector("list",niter)
TZ4_loess_co1bs1 <- vector("list",niter)
TZ4_cutoff_co1bs1 <- vector("list",niter)
TZ_co1bs1 <- vector("list",niter)
TZ4_loess_co05bs1 <- vector("list",niter)
TZ4_cutoff_co05bs1 <- vector("list",niter)
TZ_co05bs1 <- vector("list",niter)
TZ4_loess_nc10 <- vector("list",niter)
TZ4_cutoff_nc10 <- vector("list",niter)
TZ_nc10 <- vector("list",niter)
TZ4_loess_co2bs10 <- vector("list",niter)
TZ4_cutoff_co2bs10 <- vector("list",niter)
TZ_co2bs10 <- vector("list",niter)
TZ4_loess_co1bs10 <- vector("list",niter)
TZ4_cutoff_co1bs10 <- vector("list",niter)
TZ_co1bs10 <- vector("list",niter)
TZ4_loess_co05bs10 <- vector("list",niter)
TZ4_cutoff_co05bs10 <- vector("list",niter)
TZ_co05bs10 <- vector("list",niter)


set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north,Trees_count_south,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList[[i]]<-IndSpp
  
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
  NS_ratio[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.5, truncation point = 20000 #############################
  TZratio4_nc5<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc5<-hist(TZratio4_nc5$ratio_transform,breaks=seq(0,max(TZratio4_nc5$ratio_transform)+0.5,by=0.5),plot=FALSE)
  TZ4_ratio_halfbins_nc5<-data.frame(TZratio4_ratio_halfbins_nc5$mids)
  TZ4_ratio_halfbins_nc5$counts<-TZratio4_ratio_halfbins_nc5$counts
  colnames(TZ4_ratio_halfbins_nc5)<-c("bin","counts")
  TZ4_ratio_halfbins_nc5<-slice(TZ4_ratio_halfbins_nc5,-c(1,2))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_nc5 <- loess(TZ4_ratio_halfbins_nc5$counts~TZ4_ratio_halfbins_nc5$bin)
  #plot(TZ4_ratio_halfbins_nc5$bin[1:100], TZ4_ratio_halfbins_nc5$counts[1:100])
  #lines(TZ4_ratio_halfbins_nc5$bin[1:100],TZ4_halfbin_fit_nc5$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc5 <- loess(TZ4_ratio_halfbins_nc5$counts~TZ4_ratio_halfbins_nc5$bin,span = (20/nrow(TZ4_ratio_halfbins_nc5)))
  plot(TZ4_ratio_halfbins_nc5$bin[1:100], TZ4_ratio_halfbins_nc5$counts[1:100])
  lines(TZ4_ratio_halfbins_nc5$bin[1:100],TZ4_halfbin_fit_nc5$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_nc5 <- data.frame(diff(TZ4_halfbin_fit_nc5$fitted)/diff(TZ4_ratio_halfbins_nc5$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc5[2:nrow(TZ4_ratio_halfbins_nc5),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc5[[i]] <- TZ4_loess_slope_nc5
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc5<-TZ4_loess_slope_nc5 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc5[[i]]<-TZ4_halfbin_cutoff_nc5
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc5<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc5)
  TZ_4_nc5<-ContourLines2SLDF(TZcontour_4_nc5)
  
  TZ_nc5[[i]]<-TZ_4_nc5
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.5, truncation point = 200 #############################
  TZratio4_co2<-TZratio4@data %>% filter(ratio_transform<=200)
  
  # get data ready
  TZratio4_ratio_halfbins_co2<-hist(TZratio4_co2$ratio_transform,breaks=seq(0,max(TZratio4_co2$ratio_transform)+0.5,by=0.5),plot=FALSE)
  TZ4_ratio_halfbins_co2<-data.frame(TZratio4_ratio_halfbins_co2$mids)
  TZ4_ratio_halfbins_co2$counts<-TZratio4_ratio_halfbins_co2$counts
  colnames(TZ4_ratio_halfbins_co2)<-c("bin","counts")
  TZ4_ratio_halfbins_co2<-slice(TZ4_ratio_halfbins_co2,-c(1,2))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co2 <- loess(TZ4_ratio_halfbins_co2$counts~TZ4_ratio_halfbins_co2$bin)
  #plot(TZ4_ratio_halfbins_co2$bin[1:100], TZ4_ratio_halfbins_co2$counts[1:100])
  #lines(TZ4_ratio_halfbins_co2$bin[1:100],TZ4_halfbin_fit_co2$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_co2 <- loess(TZ4_ratio_halfbins_co2$counts~TZ4_ratio_halfbins_co2$bin,span = (20/nrow(TZ4_ratio_halfbins_co2)))
  plot(TZ4_ratio_halfbins_co2$bin[1:100], TZ4_ratio_halfbins_co2$counts[1:100])
  lines(TZ4_ratio_halfbins_co2$bin[1:100],TZ4_halfbin_fit_co2$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co2 <- data.frame(diff(TZ4_halfbin_fit_co2$fitted)/diff(TZ4_ratio_halfbins_co2$bin)) %>% bind_cols(TZ4_ratio_halfbins_co2[2:nrow(TZ4_ratio_halfbins_co2),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co2[[i]] <- TZ4_loess_slope_co2
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co2<-TZ4_loess_slope_co2 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co2[[i]]<-TZ4_halfbin_cutoff_co2
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co2<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co2)
  TZ_4_co2<-ContourLines2SLDF(TZcontour_4_co2)
  
  TZ_co2[[i]]<-TZ_4_co2
  
  ######################### binwidth = 0.5, truncation point = 100 #############################
  TZratio4_co1<-TZratio4@data %>% filter(ratio_transform<=100)
  
  # get data ready
  TZratio4_ratio_halfbins_co1<-hist(TZratio4_co1$ratio_transform,breaks=seq(0,max(TZratio4_co1$ratio_transform)+0.5,by=0.5),plot=FALSE)
  TZ4_ratio_halfbins_co1<-data.frame(TZratio4_ratio_halfbins_co1$mids)
  TZ4_ratio_halfbins_co1$counts<-TZratio4_ratio_halfbins_co1$counts
  colnames(TZ4_ratio_halfbins_co1)<-c("bin","counts")
  TZ4_ratio_halfbins_co1<-slice(TZ4_ratio_halfbins_co1,-c(1,2))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co1 <- loess(TZ4_ratio_halfbins_co1$counts~TZ4_ratio_halfbins_co1$bin)
  #plot(TZ4_ratio_halfbins_co1$bin[1:100], TZ4_ratio_halfbins_co1$counts[1:100])
  #lines(TZ4_ratio_halfbins_co1$bin[1:100],TZ4_halfbin_fit_co1$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_co1 <- loess(TZ4_ratio_halfbins_co1$counts~TZ4_ratio_halfbins_co1$bin,span = (20/nrow(TZ4_ratio_halfbins_co1)))
  plot(TZ4_ratio_halfbins_co1$bin[1:100], TZ4_ratio_halfbins_co1$counts[1:100])
  lines(TZ4_ratio_halfbins_co1$bin[1:100],TZ4_halfbin_fit_co1$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co1 <- data.frame(diff(TZ4_halfbin_fit_co1$fitted)/diff(TZ4_ratio_halfbins_co1$bin)) %>% bind_cols(TZ4_ratio_halfbins_co1[2:nrow(TZ4_ratio_halfbins_co1),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co1[[i]] <- TZ4_loess_slope_co1
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co1<-TZ4_loess_slope_co1 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co1[[i]]<-TZ4_halfbin_cutoff_co1
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co1<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co1)
  TZ_4_co1<-ContourLines2SLDF(TZcontour_4_co1)
  
  TZ_co1[[i]]<-TZ_4_co1
  
  ######################### binwidth = 0.5, truncation point = 50 #############################
  TZratio4_co05<-TZratio4@data %>% filter(ratio_transform<=50)
  
  # get data ready
  TZratio4_ratio_halfbins_co05<-hist(TZratio4_co05$ratio_transform,breaks=seq(0,max(TZratio4_co05$ratio_transform)+0.5,by=0.5),plot=FALSE)
  TZ4_ratio_halfbins_co05<-data.frame(TZratio4_ratio_halfbins_co05$mids)
  TZ4_ratio_halfbins_co05$counts<-TZratio4_ratio_halfbins_co05$counts
  colnames(TZ4_ratio_halfbins_co05)<-c("bin","counts")
  TZ4_ratio_halfbins_co05<-slice(TZ4_ratio_halfbins_co05,-c(1,2))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co05 <- loess(TZ4_ratio_halfbins_co05$counts~TZ4_ratio_halfbins_co05$bin)
  #plot(TZ4_ratio_halfbins_co05$bin[1:100], TZ4_ratio_halfbins_co05$counts[1:100])
  #lines(TZ4_ratio_halfbins_co05$bin[1:100],TZ4_halfbin_fit_co05$fitted[1:100])
 
  # change smoothing parameter
  TZ4_halfbin_fit_co05 <- loess(TZ4_ratio_halfbins_co05$counts~TZ4_ratio_halfbins_co05$bin,span = (20/nrow(TZ4_ratio_halfbins_co05)))
  plot(TZ4_ratio_halfbins_co05$bin[1:100], TZ4_ratio_halfbins_co05$counts[1:100])
  lines(TZ4_ratio_halfbins_co05$bin[1:100],TZ4_halfbin_fit_co05$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co05 <- data.frame(diff(TZ4_halfbin_fit_co05$fitted)/diff(TZ4_ratio_halfbins_co05$bin)) %>% bind_cols(TZ4_ratio_halfbins_co05[2:nrow(TZ4_ratio_halfbins_co05),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co05[[i]] <- TZ4_loess_slope_co05
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co05<-TZ4_loess_slope_co05 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co05[[i]]<-TZ4_halfbin_cutoff_co05
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co05<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co05)
  TZ_4_co05<-ContourLines2SLDF(TZcontour_4_co05)
  
  TZ_co05[[i]]<-TZ_4_co05
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.1, truncation point = 20000 #############################
  TZratio4_nc1<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1<-hist(TZratio4_nc1$ratio_transform,breaks=seq(0,max(TZratio4_nc1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1<-data.frame(TZratio4_ratio_halfbins_nc1$mids)
  TZ4_ratio_halfbins_nc1$counts<-TZratio4_ratio_halfbins_nc1$counts
  colnames(TZ4_ratio_halfbins_nc1)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1<-slice(TZ4_ratio_halfbins_nc1,-c(1:10))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_nc1 <- loess(TZ4_ratio_halfbins_nc1$counts~TZ4_ratio_halfbins_nc1$bin)
  #plot(TZ4_ratio_halfbins_nc1$bin[1:100], TZ4_ratio_halfbins_nc1$counts[1:100])
  #lines(TZ4_ratio_halfbins_nc1$bin[1:100],TZ4_halfbin_fit_nc1$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1 <- loess(TZ4_ratio_halfbins_nc1$counts~TZ4_ratio_halfbins_nc1$bin,span = (100/nrow(TZ4_ratio_halfbins_nc1)))
  plot(TZ4_ratio_halfbins_nc1$bin[1:100], TZ4_ratio_halfbins_nc1$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1$bin[1:100],TZ4_halfbin_fit_nc1$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_nc1 <- data.frame(diff(TZ4_halfbin_fit_nc1$fitted)/diff(TZ4_ratio_halfbins_nc1$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1[2:nrow(TZ4_ratio_halfbins_nc1),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1[[i]] <- TZ4_loess_slope_nc1
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1<-TZ4_loess_slope_nc1 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1[[i]]<-TZ4_halfbin_cutoff_nc1
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1)
  TZ_4_nc1<-ContourLines2SLDF(TZcontour_4_nc1)
  
  TZ_nc1[[i]]<-TZ_4_nc1
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.1, truncation point = 200 #############################
  TZratio4_co2bs1<-TZratio4@data %>% filter(ratio_transform<=200)
  
  # get data ready
  TZratio4_ratio_halfbins_co2bs1<-hist(TZratio4_co2bs1$ratio_transform,breaks=seq(0,max(TZratio4_co2bs1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_co2bs1<-data.frame(TZratio4_ratio_halfbins_co2bs1$mids)
  TZ4_ratio_halfbins_co2bs1$counts<-TZratio4_ratio_halfbins_co2bs1$counts
  colnames(TZ4_ratio_halfbins_co2bs1)<-c("bin","counts")
  TZ4_ratio_halfbins_co2bs1<-slice(TZ4_ratio_halfbins_co2bs1,-c(1:10))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co2bs1 <- loess(TZ4_ratio_halfbins_co2bs1$counts~TZ4_ratio_halfbins_co2bs1$bin)
  #plot(TZ4_ratio_halfbins_co2bs1$bin[1:100], TZ4_ratio_halfbins_co2bs1$counts[1:100])
  #lines(TZ4_ratio_halfbins_co2bs1$bin[1:100],TZ4_halfbin_fit_co2bs1$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_co2bs1 <- loess(TZ4_ratio_halfbins_co2bs1$counts~TZ4_ratio_halfbins_co2bs1$bin,span = (100/nrow(TZ4_ratio_halfbins_co2bs1)))
  plot(TZ4_ratio_halfbins_co2bs1$bin[1:100], TZ4_ratio_halfbins_co2bs1$counts[1:100])
  lines(TZ4_ratio_halfbins_co2bs1$bin[1:100],TZ4_halfbin_fit_co2bs1$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co2bs1 <- data.frame(diff(TZ4_halfbin_fit_co2bs1$fitted)/diff(TZ4_ratio_halfbins_co2bs1$bin)) %>% bind_cols(TZ4_ratio_halfbins_co2bs1[2:nrow(TZ4_ratio_halfbins_co2bs1),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co2bs1[[i]] <- TZ4_loess_slope_co2bs1
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co2bs1<-TZ4_loess_slope_co2bs1 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co2bs1[[i]]<-TZ4_halfbin_cutoff_co2bs1
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co2bs1<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co2bs1)
  TZ_4_co2bs1<-ContourLines2SLDF(TZcontour_4_co2bs1)
  
  TZ_co2bs1[[i]]<-TZ_4_co2bs1
  
  ######################### binwidth = 0.1, truncation point = 100 #############################
  TZratio4_co1bs1<-TZratio4@data %>% filter(ratio_transform<=100)
  
  # get data ready
  TZratio4_ratio_halfbins_co1bs1<-hist(TZratio4_co1bs1$ratio_transform,breaks=seq(0,max(TZratio4_co1bs1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_co1bs1<-data.frame(TZratio4_ratio_halfbins_co1bs1$mids)
  TZ4_ratio_halfbins_co1bs1$counts<-TZratio4_ratio_halfbins_co1bs1$counts
  colnames(TZ4_ratio_halfbins_co1bs1)<-c("bin","counts")
  TZ4_ratio_halfbins_co1bs1<-slice(TZ4_ratio_halfbins_co1bs1,-c(1:10))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co1bs1 <- loess(TZ4_ratio_halfbins_co1bs1$counts~TZ4_ratio_halfbins_co1bs1$bin)
  #plot(TZ4_ratio_halfbins_co1bs1$bin[1:100], TZ4_ratio_halfbins_co1bs1$counts[1:100])
  #lines(TZ4_ratio_halfbins_co1bs1$bin[1:100],TZ4_halfbin_fit_co1bs1$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_co1bs1 <- loess(TZ4_ratio_halfbins_co1bs1$counts~TZ4_ratio_halfbins_co1bs1$bin,span = (100/nrow(TZ4_ratio_halfbins_co1bs1)))
  plot(TZ4_ratio_halfbins_co1bs1$bin[1:100], TZ4_ratio_halfbins_co1bs1$counts[1:100])
  lines(TZ4_ratio_halfbins_co1bs1$bin[1:100],TZ4_halfbin_fit_co1bs1$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co1bs1 <- data.frame(diff(TZ4_halfbin_fit_co1bs1$fitted)/diff(TZ4_ratio_halfbins_co1bs1$bin)) %>% bind_cols(TZ4_ratio_halfbins_co1bs1[2:nrow(TZ4_ratio_halfbins_co1bs1),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co1bs1[[i]] <- TZ4_loess_slope_co1bs1
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co1bs1<-TZ4_loess_slope_co1bs1 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co1bs1[[i]]<-TZ4_halfbin_cutoff_co1bs1
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co1bs1<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co1bs1)
  TZ_4_co1bs1<-ContourLines2SLDF(TZcontour_4_co1bs1)
  
  TZ_co1bs1[[i]]<-TZ_4_co1bs1
  
  ######################### binwidth = 0.1, truncation point = 50 #############################
  TZratio4_co05bs1<-TZratio4@data %>% filter(ratio_transform<=50)
  
  # get data ready
  TZratio4_ratio_halfbins_co05bs1<-hist(TZratio4_co05bs1$ratio_transform,breaks=seq(0,max(TZratio4_co05bs1$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_co05bs1<-data.frame(TZratio4_ratio_halfbins_co05bs1$mids)
  TZ4_ratio_halfbins_co05bs1$counts<-TZratio4_ratio_halfbins_co05bs1$counts
  colnames(TZ4_ratio_halfbins_co05bs1)<-c("bin","counts")
  TZ4_ratio_halfbins_co05bs1<-slice(TZ4_ratio_halfbins_co05bs1,-c(1:10))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co05bs1 <- loess(TZ4_ratio_halfbins_co05bs1$counts~TZ4_ratio_halfbins_co05bs1$bin)
  #plot(TZ4_ratio_halfbins_co05bs1$bin[1:100], TZ4_ratio_halfbins_co05bs1$counts[1:100])
  #lines(TZ4_ratio_halfbins_co05bs1$bin[1:100],TZ4_halfbin_fit_co05bs1$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_co05bs1 <- loess(TZ4_ratio_halfbins_co05bs1$counts~TZ4_ratio_halfbins_co05bs1$bin,span = (100/nrow(TZ4_ratio_halfbins_co05bs1)))
  plot(TZ4_ratio_halfbins_co05bs1$bin[1:100], TZ4_ratio_halfbins_co05bs1$counts[1:100])
  lines(TZ4_ratio_halfbins_co05bs1$bin[1:100],TZ4_halfbin_fit_co05bs1$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co05bs1 <- data.frame(diff(TZ4_halfbin_fit_co05bs1$fitted)/diff(TZ4_ratio_halfbins_co05bs1$bin)) %>% bind_cols(TZ4_ratio_halfbins_co05bs1[2:nrow(TZ4_ratio_halfbins_co05bs1),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co05bs1[[i]] <- TZ4_loess_slope_co05bs1
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co05bs1<-TZ4_loess_slope_co05bs1 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co05bs1[[i]]<-TZ4_halfbin_cutoff_co05bs1
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co05bs1<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co05bs1)
  TZ_4_co05bs1<-ContourLines2SLDF(TZcontour_4_co05bs1)
  
  TZ_co05bs1[[i]]<-TZ_4_co05bs1
  
  ## getting raw data from histogram bins
  ######################### binwidth = 1, truncation point = 20000 #############################
  TZratio4_nc10<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc10<-hist(TZratio4_nc10$ratio_transform,breaks=seq(0,max(TZratio4_nc10$ratio_transform)+1,by=1),plot=FALSE)
  TZ4_ratio_halfbins_nc10<-data.frame(TZratio4_ratio_halfbins_nc10$mids)
  TZ4_ratio_halfbins_nc10$counts<-TZratio4_ratio_halfbins_nc10$counts
  colnames(TZ4_ratio_halfbins_nc10)<-c("bin","counts")
  TZ4_ratio_halfbins_nc10<-slice(TZ4_ratio_halfbins_nc10,-c(1))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_nc10 <- loess(TZ4_ratio_halfbins_nc10$counts~TZ4_ratio_halfbins_nc10$bin)
  #plot(TZ4_ratio_halfbins_nc10$bin[1:100], TZ4_ratio_halfbins_nc10$counts[1:100])
  #lines(TZ4_ratio_halfbins_nc10$bin[1:100],TZ4_halfbin_fit_nc10$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc10 <- loess(TZ4_ratio_halfbins_nc10$counts~TZ4_ratio_halfbins_nc10$bin,span = (10/nrow(TZ4_ratio_halfbins_nc10)))
  plot(TZ4_ratio_halfbins_nc10$bin[1:100], TZ4_ratio_halfbins_nc10$counts[1:100])
  lines(TZ4_ratio_halfbins_nc10$bin[1:100],TZ4_halfbin_fit_nc10$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_nc10 <- data.frame(diff(TZ4_halfbin_fit_nc10$fitted)/diff(TZ4_ratio_halfbins_nc10$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc10[2:nrow(TZ4_ratio_halfbins_nc10),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc10[[i]] <- TZ4_loess_slope_nc10
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc10<-TZ4_loess_slope_nc10 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc10[[i]]<-TZ4_halfbin_cutoff_nc10
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc10<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc10)
  TZ_4_nc10<-ContourLines2SLDF(TZcontour_4_nc10)
  
  TZ_nc10[[i]]<-TZ_4_nc10
  
  ## getting raw data from histogram bins
  ######################### binwidth = 1, truncation point = 200 #############################
  TZratio4_co2bs10<-TZratio4@data %>% filter(ratio_transform<=200)
  
  # get data ready
  TZratio4_ratio_halfbins_co2bs10<-hist(TZratio4_co2bs10$ratio_transform,breaks=seq(0,max(TZratio4_co2bs10$ratio_transform)+1,by=1),plot=FALSE)
  TZ4_ratio_halfbins_co2bs10<-data.frame(TZratio4_ratio_halfbins_co2bs10$mids)
  TZ4_ratio_halfbins_co2bs10$counts<-TZratio4_ratio_halfbins_co2bs10$counts
  colnames(TZ4_ratio_halfbins_co2bs10)<-c("bin","counts")
  TZ4_ratio_halfbins_co2bs10<-slice(TZ4_ratio_halfbins_co2bs10,-c(1))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co2bs10 <- loess(TZ4_ratio_halfbins_co2bs10$counts~TZ4_ratio_halfbins_co2bs10$bin)
  #plot(TZ4_ratio_halfbins_co2bs10$bin[1:100], TZ4_ratio_halfbins_co2bs10$counts[1:100])
  #lines(TZ4_ratio_halfbins_co2bs10$bin[1:100],TZ4_halfbin_fit_co2bs10$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_co2bs10 <- loess(TZ4_ratio_halfbins_co2bs10$counts~TZ4_ratio_halfbins_co2bs10$bin,span = (10/nrow(TZ4_ratio_halfbins_co2bs10)))
  plot(TZ4_ratio_halfbins_co2bs10$bin[1:100], TZ4_ratio_halfbins_co2bs10$counts[1:100])
  lines(TZ4_ratio_halfbins_co2bs10$bin[1:100],TZ4_halfbin_fit_co2bs10$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co2bs10 <- data.frame(diff(TZ4_halfbin_fit_co2bs10$fitted)/diff(TZ4_ratio_halfbins_co2bs10$bin)) %>% bind_cols(TZ4_ratio_halfbins_co2bs10[2:nrow(TZ4_ratio_halfbins_co2bs10),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co2bs10[[i]] <- TZ4_loess_slope_co2bs10
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co2bs10<-TZ4_loess_slope_co2bs10 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co2bs10[[i]]<-TZ4_halfbin_cutoff_co2bs10
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co2bs10<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co2bs10)
  TZ_4_co2bs10<-ContourLines2SLDF(TZcontour_4_co2bs10)
  
  TZ_co2bs10[[i]]<-TZ_4_co2bs10
  
  ######################### binwidth = 1, truncation point = 100 #############################
  TZratio4_co1bs10<-TZratio4@data %>% filter(ratio_transform<=100)
  
  # get data ready
  TZratio4_ratio_halfbins_co1bs10<-hist(TZratio4_co1bs10$ratio_transform,breaks=seq(0,max(TZratio4_co1bs10$ratio_transform)+1,by=1),plot=FALSE)
  TZ4_ratio_halfbins_co1bs10<-data.frame(TZratio4_ratio_halfbins_co1bs10$mids)
  TZ4_ratio_halfbins_co1bs10$counts<-TZratio4_ratio_halfbins_co1bs10$counts
  colnames(TZ4_ratio_halfbins_co1bs10)<-c("bin","counts")
  TZ4_ratio_halfbins_co1bs10<-slice(TZ4_ratio_halfbins_co1bs10,-c(1))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co1bs10 <- loess(TZ4_ratio_halfbins_co1bs10$counts~TZ4_ratio_halfbins_co1bs10$bin)
  #plot(TZ4_ratio_halfbins_co1bs10$bin[1:100], TZ4_ratio_halfbins_co1bs10$counts[1:100])
  #lines(TZ4_ratio_halfbins_co1bs10$bin[1:100],TZ4_halfbin_fit_co1bs10$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_co1bs10 <- loess(TZ4_ratio_halfbins_co1bs10$counts~TZ4_ratio_halfbins_co1bs10$bin,span = (10/nrow(TZ4_ratio_halfbins_co1bs10)))
  plot(TZ4_ratio_halfbins_co1bs10$bin[1:100], TZ4_ratio_halfbins_co1bs10$counts[1:100])
  lines(TZ4_ratio_halfbins_co1bs10$bin[1:100],TZ4_halfbin_fit_co1bs10$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co1bs10 <- data.frame(diff(TZ4_halfbin_fit_co1bs10$fitted)/diff(TZ4_ratio_halfbins_co1bs10$bin)) %>% bind_cols(TZ4_ratio_halfbins_co1bs10[2:nrow(TZ4_ratio_halfbins_co1bs10),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co1bs10[[i]] <- TZ4_loess_slope_co1bs10
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co1bs10<-TZ4_loess_slope_co1bs10 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co1bs10[[i]]<-TZ4_halfbin_cutoff_co1bs10
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co1bs10<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co1bs10)
  TZ_4_co1bs10<-ContourLines2SLDF(TZcontour_4_co1bs10)
  
  TZ_co1bs10[[i]]<-TZ_4_co1bs10
  
  ######################### binwidth = 1, truncation point = 50 #############################
  TZratio4_co05bs10<-TZratio4@data %>% filter(ratio_transform<=50)
  
  # get data ready
  TZratio4_ratio_halfbins_co05bs10<-hist(TZratio4_co05bs10$ratio_transform,breaks=seq(0,max(TZratio4_co05bs10$ratio_transform)+1,by=1),plot=FALSE)
  TZ4_ratio_halfbins_co05bs10<-data.frame(TZratio4_ratio_halfbins_co05bs10$mids)
  TZ4_ratio_halfbins_co05bs10$counts<-TZratio4_ratio_halfbins_co05bs10$counts
  colnames(TZ4_ratio_halfbins_co05bs10)<-c("bin","counts")
  TZ4_ratio_halfbins_co05bs10<-slice(TZ4_ratio_halfbins_co05bs10,-c(1))
  
  # use loess to fit a line
  #TZ4_halfbin_fit_co05bs10 <- loess(TZ4_ratio_halfbins_co05bs10$counts~TZ4_ratio_halfbins_co05bs10$bin)
  #plot(TZ4_ratio_halfbins_co05bs10$bin[1:100], TZ4_ratio_halfbins_co05bs10$counts[1:100])
  #lines(TZ4_ratio_halfbins_co05bs10$bin[1:100],TZ4_halfbin_fit_co05bs10$fitted[1:100])
  
  # change smoothing parameter
  TZ4_halfbin_fit_co05bs10 <- loess(TZ4_ratio_halfbins_co05bs10$counts~TZ4_ratio_halfbins_co05bs10$bin,span = (10/nrow(TZ4_ratio_halfbins_co05bs10)))
  plot(TZ4_ratio_halfbins_co05bs10$bin[1:100], TZ4_ratio_halfbins_co05bs10$counts[1:100])
  lines(TZ4_ratio_halfbins_co05bs10$bin[1:100],TZ4_halfbin_fit_co05bs10$fitted[1:100])
  
  # local slope
  TZ4_loess_slope_co05bs10 <- data.frame(diff(TZ4_halfbin_fit_co05bs10$fitted)/diff(TZ4_ratio_halfbins_co05bs10$bin)) %>% bind_cols(TZ4_ratio_halfbins_co05bs10[2:nrow(TZ4_ratio_halfbins_co05bs10),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_co05bs10[[i]] <- TZ4_loess_slope_co05bs10
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_co05bs10<-TZ4_loess_slope_co05bs10 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_co05bs10[[i]]<-TZ4_halfbin_cutoff_co05bs10
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_co05bs10<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_co05bs10)
  TZ_4_co05bs10<-ContourLines2SLDF(TZcontour_4_co05bs10)
  
  TZ_co05bs10[[i]]<-TZ_4_co05bs10
}

m_nc5<-do.call(bind,TZ_nc5)
m_co2<-do.call(bind,TZ_co2)
m_co1<-do.call(bind,TZ_co1)
m_co05<-do.call(bind,TZ_co05)
m_nc1<-do.call(bind,TZ_nc1)
m_co2bs1<-do.call(bind,TZ_co2bs1)
m_co1bs1<-do.call(bind,TZ_co1bs1)
m_co05bs1<-do.call(bind,TZ_co05bs1)
m_nc10<-do.call(bind,TZ_nc10)
m_co2bs10<-do.call(bind,TZ_co2bs10)
m_co1bs10<-do.call(bind,TZ_co1bs10)
m_co05bs10<-do.call(bind,TZ_co05bs10)

TZ4_loess_nc5_df <- do.call(bind,TZ4_loess_nc5)
TZ4_loess_co2_df <- do.call(bind,TZ4_loess_co2)
TZ4_loess_co1_df <- do.call(bind,TZ4_loess_co1)
TZ4_loess_co05_df <- do.call(bind,TZ4_loess_co05)
TZ4_loess_nc1_df <- do.call(bind,TZ4_loess_nc1)
TZ4_loess_co2bs1_df <- do.call(bind,TZ4_loess_co2bs1)
TZ4_loess_co1bs1_df <- do.call(bind,TZ4_loess_co1bs1)
TZ4_loess_co05bs1_df <- do.call(bind,TZ4_loess_co05bs1)
TZ4_loess_nc10_df <- do.call(bind,TZ4_loess_nc10)
TZ4_loess_co2bs10_df <- do.call(bind,TZ4_loess_co2bs10)
TZ4_loess_co1bs10_df <- do.call(bind,TZ4_loess_co1bs10)
TZ4_loess_co05bs10_df <- do.call(bind,TZ4_loess_co05bs10)

low_nc1<-TZ4_loess_nc1_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))
low_nc5<-TZ4_loess_nc5_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))

TZ4_loess_nc5_df$ma <- runMean(TZ4_loess_nc5_df$counts,10)
TZ4_loess_co2_df$ma <- runMean(TZ4_loess_co2_df$counts,10)
TZ4_loess_co1_df$ma <- runMean(TZ4_loess_co1_df$counts,10)
TZ4_loess_co05_df$ma <- runMean(TZ4_loess_co05_df$counts,10)
TZ4_loess_nc1_df$ma <- runMean(TZ4_loess_nc1_df$counts,10)
TZ4_loess_co2bs1_df$ma <- runMean(TZ4_loess_co2bs1_df$counts,10)
TZ4_loess_co1bs1_df$ma <- runMean(TZ4_loess_co1bs1_df$counts,10)
TZ4_loess_co05bs1_df$ma <- runMean(TZ4_loess_co05bs1_df$counts,10)
TZ4_loess_nc10_df$ma <- runMean(TZ4_loess_nc10_df$counts,10)
TZ4_loess_co2bs10_df$ma <- runMean(TZ4_loess_co2bs10_df$counts,10)
TZ4_loess_co1bs10_df$ma <- runMean(TZ4_loess_co1bs10_df$counts,10)
TZ4_loess_co05bs10_df$ma <- runMean(TZ4_loess_co05bs10_df$counts,10)

low_ma_nc5<-TZ4_loess_nc5_df %>% filter(ma<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))

## Relative ratio of northern and southern trees in Wisconsin
# take mean of all ratio maps
NS_ratiotrans_table<-matrix(unlist(NS_ratiotrans),ncol=niter)
NS_ratiotrans_mean<-rowMeans(NS_ratiotrans_table)

# make new SGDF with mean ratiotrans
TZratio0<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0)<-"North"
TZratio0$South<-as(South.dens4,"SpatialGridDataFrame")$v
TZratio0$mean_ratiotrans<-NS_ratiotrans_mean

# convert to spatial pixel data frame
TZratio_mean<-as(TZratio0,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_nc5,a=0.2) # set formatting for the 50 lines

TZratio_mean@data$cut_mean <- cut(TZratio_mean@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
Fig_nc5 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(iterline,wisc_HARN),
                col.regions=rev(brewer.pal(11,"RdYlBu")),
                colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                par.settings = list(axis.line = list(col =  'transparent')),
                auto.key = list(title = "Relative ratio"),
                main="Cutoff = 20000, bin = 0.5")

Fig_nc5

iterline<-list("sp.lines",m_co2,a=0.2) # set formatting for the 50 lines

Fig_co2 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(iterline,wisc_HARN),
                 col.regions=rev(brewer.pal(11,"RdYlBu")),
                 colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                 par.settings = list(axis.line = list(col =  'transparent')),
                 auto.key = list(title = "Relative ratio"),
                 main="Cutoff = 200, bin = 0.5")

Fig_co2

## ratio summaries
TZ_cutoff_nc5_ratio <- data.frame(unlist(TZ4_cutoff_nc5)) %>% mutate(cutoff=20000, bin=0.5) %>% rename(cp=1)
TZ_cutoff_co2_ratio <- data.frame(unlist(TZ4_cutoff_co2)) %>% mutate(cutoff=200, bin=0.5) %>% rename(cp=1)
TZ_cutoff_co1_ratio <- data.frame(unlist(TZ4_cutoff_co1)) %>% mutate(cutoff=100, bin=0.5) %>% rename(cp=1)
TZ_cutoff_co05_ratio <- data.frame(unlist(TZ4_cutoff_co05)) %>% mutate(cutoff=50, bin=0.5) %>% rename(cp=1)
TZ_cutoff_nc1_ratio <- data.frame(unlist(TZ4_cutoff_nc1)) %>% mutate(cutoff=20000, bin=0.1) %>% rename(cp=1)
TZ_cutoff_co2bs1_ratio <- data.frame(unlist(TZ4_cutoff_co2bs1)) %>% mutate(cutoff=200, bin=0.1) %>% rename(cp=1)
TZ_cutoff_co1bs1_ratio <- data.frame(unlist(TZ4_cutoff_co1bs1)) %>% mutate(cutoff=100, bin=0.1) %>% rename(cp=1)
TZ_cutoff_co05bs1_ratio <- data.frame(unlist(TZ4_cutoff_co05bs1)) %>% mutate(cutoff=50, bin=0.1) %>% rename(cp=1)
TZ_cutoff_nc10_ratio <- data.frame(unlist(TZ4_cutoff_nc10)) %>% mutate(cutoff=20000, bin=1) %>% rename(cp=1)
TZ_cutoff_co2bs10_ratio <- data.frame(unlist(TZ4_cutoff_co2bs10)) %>% mutate(cutoff=200, bin=1) %>% rename(cp=1)
TZ_cutoff_co1bs10_ratio <- data.frame(unlist(TZ4_cutoff_co1bs10)) %>% mutate(cutoff=100, bin=1) %>% rename(cp=1)
TZ_cutoff_co05bs10_ratio <- data.frame(unlist(TZ4_cutoff_co05bs10)) %>% mutate(cutoff=50, bin=1) %>% rename(cp=1)

TZ_cutoff_ratio_all <- bind_rows(TZ_cutoff_nc5_ratio,TZ_cutoff_co2_ratio,TZ_cutoff_co1_ratio,TZ_cutoff_co05_ratio,
                                 TZ_cutoff_nc1_ratio,TZ_cutoff_co2bs1_ratio,TZ_cutoff_co1bs1_ratio,TZ_cutoff_co05bs1_ratio,
                                 TZ_cutoff_nc10_ratio,TZ_cutoff_co2bs10_ratio,TZ_cutoff_co1bs10_ratio,TZ_cutoff_co05bs10_ratio) %>% mutate(cutoff=as.factor(cutoff))

TZ_cutoff_ratio_all$bin_label <- ifelse(TZ_cutoff_ratio_all$bin==0.1,"Bin = 0.1",(ifelse(TZ_cutoff_ratio_all$bin==0.5,"Bin = 0.5","Bin = 1.0")))

### boxplot of CP results
ggplot(data=TZ_cutoff_ratio_all,aes(x=cutoff,y=cp))+
  geom_boxplot()+
  facet_wrap(~bin_label)+
  theme_bw()+
  ylab("Critical point")+
  xlab("Truncation ratio value")+
  theme(strip.background = element_blank())






## Average ratio line
avg_cutoff_nc5<-mean(unlist(TZ4_cutoff_nc5))
sd_cutoff_nc5<-sd(unlist(TZ4_cutoff_nc5))

avg_cutoff_co2<-mean(unlist(TZ4_cutoff_co2))
sd_cutoff_co2<-sd(unlist(TZ4_cutoff_co2))

avg_cutoff_co1<-mean(unlist(TZ4_cutoff_co1))
sd_cutoff_co1<-sd(unlist(TZ4_cutoff_co1))

avg_cutoff_co05<-mean(unlist(TZ4_cutoff_co05))
sd_cutoff_co05<-sd(unlist(TZ4_cutoff_co05))

avg_cutoff_nc1<-mean(unlist(TZ4_cutoff_nc1))
sd_cutoff_nc1<-sd(unlist(TZ4_cutoff_nc1))

avg_cutoff_co2bs1<-mean(unlist(TZ4_cutoff_co2bs1))
sd_cutoff_co2bs1<-sd(unlist(TZ4_cutoff_co2bs1))

avg_cutoff_nc10<-mean(unlist(TZ4_cutoff_nc10))
sd_cutoff_nc10<-sd(unlist(TZ4_cutoff_nc10))

avg_cutoff_co2bs10<-mean(unlist(TZ4_cutoff_co2bs10))
sd_cutoff_co2bs10<-sd(unlist(TZ4_cutoff_co2bs10))

# map of average cutoff on average ratio map
TZratio4_mean_image_nc5<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc5<-contourLines(TZratio4_mean_image_nc5,levels=avg_cutoff_nc5)
TZ_mean_4_nc5<-ContourLines2SLDF(TZcontour_mean_4_nc5)

TZratio4_mean_image_co2<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_co2<-contourLines(TZratio4_mean_image_co2,levels=avg_cutoff_co2)
TZ_mean_4_co2<-ContourLines2SLDF(TZcontour_mean_4_co2)

TZratio4_mean_image_co1<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_co1<-contourLines(TZratio4_mean_image_co1,levels=avg_cutoff_co1)
TZ_mean_4_co1<-ContourLines2SLDF(TZcontour_mean_4_co1)

TZratio4_mean_image_co05<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_co05<-contourLines(TZratio4_mean_image_co05,levels=avg_cutoff_co05)
TZ_mean_4_co05<-ContourLines2SLDF(TZcontour_mean_4_co05)

TZratio4_mean_image_nc1<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1<-contourLines(TZratio4_mean_image_nc1,levels=avg_cutoff_nc1)
TZ_mean_4_nc1<-ContourLines2SLDF(TZcontour_mean_4_nc1)

TZratio4_mean_image_co2bs1<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_co2bs1<-contourLines(TZratio4_mean_image_co2bs1,levels=avg_cutoff_co2bs1)
TZ_mean_4_co2bs1<-ContourLines2SLDF(TZcontour_mean_4_co2bs1)

TZratio4_mean_image_nc10<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc10<-contourLines(TZratio4_mean_image_nc10,levels=avg_cutoff_nc10)
TZ_mean_4_nc10<-ContourLines2SLDF(TZcontour_mean_4_nc10)

TZratio4_mean_image_co2bs10<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_co2bs10<-contourLines(TZratio4_mean_image_co2bs10,levels=avg_cutoff_co2bs10)
TZ_mean_4_co2bs10<-ContourLines2SLDF(TZcontour_mean_4_co2bs10)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean@data$cut_mean <- cut(TZratio_mean@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")

mean_plot_nc5 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc5,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"),
                       main = "bin size = 0.5, trunction point = 20,000")

mean_plot_co2 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_co2,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"),
                       main = "bin size = 0.5, trunction point = 200")

mean_plot_co1 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_co1,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"),
                       main = "bin size = 0.5, trunction point = 100")

mean_plot_co05 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_co05,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"),
                       main = "bin size = 0.5, trunction point = 50")

mean_plot_nc1 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"),
                       main = "bin size = 0.1, trunction point = 20,000")

mean_plot_co2bs1 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_co2bs1,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 200")

mean_plot_nc10 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc10,wisc_HARN),
                       col.regions=rev(brewer.pal(11,"RdYlBu")),
                       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                       par.settings = list(axis.line = list(col =  'transparent')),
                       auto.key = list(title = "Relative ratio"),
                       main = "bin size = 1, trunction point = 20,000")

mean_plot_co2bs10 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_co2bs10,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 1, trunction point = 200")


mean_plot_nc5
mean_plot_co2
mean_plot_co1
mean_plot_co05
mean_plot_nc1
mean_plot_co2bs1
mean_plot_nc10
mean_plot_co2bs10


########################### going with bin size = 0.1 and truncation point = 20,000 ###############################

########################### testing different spans ##########################

#### mapping ####
### iterations
niter<-50
IndSppList<-vector("list",niter)
NS_ratio <- vector("list",niter)
NS_logratio <- vector("list",niter)
NS_ratiotrans <- vector("list",niter)
TZ4_loess_nc1_20 <- vector("list",niter)
TZ4_cutoff_nc1_20 <- vector("list",niter)
TZ_nc1_20 <- vector("list",niter)
TZ4_loess_resid_nc1_20_list <- vector("list",niter)
TZ4_loess_nc1_30 <- vector("list",niter)
TZ4_cutoff_nc1_30 <- vector("list",niter)
TZ_nc1_30 <- vector("list",niter)
TZ4_loess_resid_nc1_30_list <- vector("list",niter)
TZ4_loess_nc1_40 <- vector("list",niter)
TZ4_cutoff_nc1_40 <- vector("list",niter)
TZ_nc1_40 <- vector("list",niter)
TZ4_loess_resid_nc1_40_list <- vector("list",niter)
TZ4_loess_nc1_50 <- vector("list",niter)
TZ4_cutoff_nc1_50 <- vector("list",niter)
TZ_nc1_50 <- vector("list",niter)
TZ4_loess_resid_nc1_50_list <- vector("list",niter)
TZ4_loess_nc1_60 <- vector("list",niter)
TZ4_cutoff_nc1_60 <- vector("list",niter)
TZ_nc1_60 <- vector("list",niter)
TZ4_loess_resid_nc1_60_list <- vector("list",niter)
TZ4_loess_nc1_70 <- vector("list",niter)
TZ4_cutoff_nc1_70 <- vector("list",niter)
TZ_nc1_70 <- vector("list",niter)
TZ4_loess_resid_nc1_70_list <- vector("list",niter)
TZ4_loess_nc1_80 <- vector("list",niter)
TZ4_cutoff_nc1_80 <- vector("list",niter)
TZ_nc1_80 <- vector("list",niter)
TZ4_loess_resid_nc1_80_list <- vector("list",niter)
TZ4_loess_nc1_90 <- vector("list",niter)
TZ4_cutoff_nc1_90 <- vector("list",niter)
TZ_nc1_90 <- vector("list",niter)
TZ4_loess_resid_nc1_90_list <- vector("list",niter)

set.seed(6542)
for (i in 1:niter){
  #### ISA
  results<-ISA_PLS(Trees_count_north,Trees_count_south,50)
  
  # subset of trees that are identifies as northern, southern, or both
  IndSpp<-results[which(results[,7]>0),]
  NSpp<-IndSpp[which(IndSpp[,2]>0),]
  SSpp<-IndSpp[which(IndSpp[,3]>0),]
  IndSppList[[i]]<-IndSpp
  
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
  NS_ratio[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ## getting raw data from histogram bins
  ######################### binwidth = 0.1, truncation point = 20000 #############################
  ############## span = 20 bins ###############
  TZratio4_nc1_20<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1_20<-hist(TZratio4_nc1_20$ratio_transform,breaks=seq(0,max(TZratio4_nc1_20$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1_20<-data.frame(TZratio4_ratio_halfbins_nc1_20$mids)
  TZ4_ratio_halfbins_nc1_20$counts<-TZratio4_ratio_halfbins_nc1_20$counts
  colnames(TZ4_ratio_halfbins_nc1_20)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1_20<-slice(TZ4_ratio_halfbins_nc1_20,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1_20 <- loess(TZ4_ratio_halfbins_nc1_20$counts~TZ4_ratio_halfbins_nc1_20$bin,span = (20/nrow(TZ4_ratio_halfbins_nc1_20)))
  plot(TZ4_ratio_halfbins_nc1_20$bin[1:100], TZ4_ratio_halfbins_nc1_20$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1_20$bin[1:100],TZ4_halfbin_fit_nc1_20$fitted[1:100])
  
  TZ4_loess_resid_nc1_20 <- data.frame(TZ4_ratio_halfbins_nc1_20$counts,TZ4_halfbin_fit_nc1_20$residual) %>% rename(count=1,residual=2)
  TZ4_loess_resid_nc1_20_list[[i]] <- TZ4_loess_resid_nc1_20
  
  # local slope
  TZ4_loess_slope_nc1_20 <- data.frame(diff(TZ4_halfbin_fit_nc1_20$fitted)/diff(TZ4_ratio_halfbins_nc1_20$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1_20[2:nrow(TZ4_ratio_halfbins_nc1_20),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1_20[[i]] <- TZ4_loess_slope_nc1_20
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1_20<-TZ4_loess_slope_nc1_20 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1_20[[i]]<-TZ4_halfbin_cutoff_nc1_20
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1_20<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1_20)
  TZ_4_nc1_20<-ContourLines2SLDF(TZcontour_4_nc1_20)
  
  TZ_nc1_20[[i]]<-TZ_4_nc1_20
  
  ############## span = 30 bins ###############
  TZratio4_nc1_30<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1_30<-hist(TZratio4_nc1_30$ratio_transform,breaks=seq(0,max(TZratio4_nc1_30$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1_30<-data.frame(TZratio4_ratio_halfbins_nc1_30$mids)
  TZ4_ratio_halfbins_nc1_30$counts<-TZratio4_ratio_halfbins_nc1_30$counts
  colnames(TZ4_ratio_halfbins_nc1_30)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1_30<-slice(TZ4_ratio_halfbins_nc1_30,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1_30 <- loess(TZ4_ratio_halfbins_nc1_30$counts~TZ4_ratio_halfbins_nc1_30$bin,span = (30/nrow(TZ4_ratio_halfbins_nc1_30)))
  plot(TZ4_ratio_halfbins_nc1_30$bin[1:100], TZ4_ratio_halfbins_nc1_30$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1_30$bin[1:100],TZ4_halfbin_fit_nc1_30$fitted[1:100])
  
  TZ4_loess_resid_nc1_30 <- data.frame(TZ4_ratio_halfbins_nc1_30$counts,TZ4_halfbin_fit_nc1_30$residual) %>% rename(count=1,residual=2)
  TZ4_loess_resid_nc1_30_list[[i]] <- TZ4_loess_resid_nc1_30
  
  # local slope
  TZ4_loess_slope_nc1_30 <- data.frame(diff(TZ4_halfbin_fit_nc1_30$fitted)/diff(TZ4_ratio_halfbins_nc1_30$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1_30[2:nrow(TZ4_ratio_halfbins_nc1_30),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1_30[[i]] <- TZ4_loess_slope_nc1_30
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1_30<-TZ4_loess_slope_nc1_30 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1_30[[i]]<-TZ4_halfbin_cutoff_nc1_30
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1_30<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1_30)
  TZ_4_nc1_30<-ContourLines2SLDF(TZcontour_4_nc1_30)
  
  TZ_nc1_30[[i]]<-TZ_4_nc1_30
  
  ############## span = 40 bins ###############
  TZratio4_nc1_40<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1_40<-hist(TZratio4_nc1_40$ratio_transform,breaks=seq(0,max(TZratio4_nc1_40$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1_40<-data.frame(TZratio4_ratio_halfbins_nc1_40$mids)
  TZ4_ratio_halfbins_nc1_40$counts<-TZratio4_ratio_halfbins_nc1_40$counts
  colnames(TZ4_ratio_halfbins_nc1_40)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1_40<-slice(TZ4_ratio_halfbins_nc1_40,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1_40 <- loess(TZ4_ratio_halfbins_nc1_40$counts~TZ4_ratio_halfbins_nc1_40$bin,span = (40/nrow(TZ4_ratio_halfbins_nc1_40)))
  plot(TZ4_ratio_halfbins_nc1_40$bin[1:100], TZ4_ratio_halfbins_nc1_40$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1_40$bin[1:100],TZ4_halfbin_fit_nc1_40$fitted[1:100])
  
  TZ4_loess_resid_nc1_40 <- data.frame(TZ4_ratio_halfbins_nc1_40$counts,TZ4_halfbin_fit_nc1_40$residual) %>% rename(count=1,residual=2)
  TZ4_loess_resid_nc1_40_list[[i]] <- TZ4_loess_resid_nc1_40
  
  # local slope
  TZ4_loess_slope_nc1_40 <- data.frame(diff(TZ4_halfbin_fit_nc1_40$fitted)/diff(TZ4_ratio_halfbins_nc1_40$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1_40[2:nrow(TZ4_ratio_halfbins_nc1_40),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1_40[[i]] <- TZ4_loess_slope_nc1_40
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1_40<-TZ4_loess_slope_nc1_40 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1_40[[i]]<-TZ4_halfbin_cutoff_nc1_40
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1_40<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1_40)
  TZ_4_nc1_40<-ContourLines2SLDF(TZcontour_4_nc1_40)
  
  TZ_nc1_40[[i]]<-TZ_4_nc1_40
  
  ############## span = 50 bins ###############
  TZratio4_nc1_50<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1_50<-hist(TZratio4_nc1_50$ratio_transform,breaks=seq(0,max(TZratio4_nc1_50$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1_50<-data.frame(TZratio4_ratio_halfbins_nc1_50$mids)
  TZ4_ratio_halfbins_nc1_50$counts<-TZratio4_ratio_halfbins_nc1_50$counts
  colnames(TZ4_ratio_halfbins_nc1_50)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1_50<-slice(TZ4_ratio_halfbins_nc1_50,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1_50 <- loess(TZ4_ratio_halfbins_nc1_50$counts~TZ4_ratio_halfbins_nc1_50$bin,span = (50/nrow(TZ4_ratio_halfbins_nc1_50)))
  plot(TZ4_ratio_halfbins_nc1_50$bin[1:100], TZ4_ratio_halfbins_nc1_50$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1_50$bin[1:100],TZ4_halfbin_fit_nc1_50$fitted[1:100])
  
  TZ4_loess_resid_nc1_50 <- data.frame(TZ4_ratio_halfbins_nc1_50$counts,TZ4_halfbin_fit_nc1_50$residual) %>% rename(count=1,residual=2)
  TZ4_loess_resid_nc1_50_list[[i]] <- TZ4_loess_resid_nc1_50
  
  # local slope
  TZ4_loess_slope_nc1_50 <- data.frame(diff(TZ4_halfbin_fit_nc1_50$fitted)/diff(TZ4_ratio_halfbins_nc1_50$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1_50[2:nrow(TZ4_ratio_halfbins_nc1_50),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1_50[[i]] <- TZ4_loess_slope_nc1_50
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1_50<-TZ4_loess_slope_nc1_50 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1_50[[i]]<-TZ4_halfbin_cutoff_nc1_50
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1_50<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1_50)
  TZ_4_nc1_50<-ContourLines2SLDF(TZcontour_4_nc1_50)
  
  TZ_nc1_50[[i]]<-TZ_4_nc1_50
  
  ############## span = 60 bins ###############
  TZratio4_nc1_60<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1_60<-hist(TZratio4_nc1_60$ratio_transform,breaks=seq(0,max(TZratio4_nc1_60$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1_60<-data.frame(TZratio4_ratio_halfbins_nc1_60$mids)
  TZ4_ratio_halfbins_nc1_60$counts<-TZratio4_ratio_halfbins_nc1_60$counts
  colnames(TZ4_ratio_halfbins_nc1_60)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1_60<-slice(TZ4_ratio_halfbins_nc1_60,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1_60 <- loess(TZ4_ratio_halfbins_nc1_60$counts~TZ4_ratio_halfbins_nc1_60$bin,span = (60/nrow(TZ4_ratio_halfbins_nc1_60)))
  plot(TZ4_ratio_halfbins_nc1_60$bin[1:100], TZ4_ratio_halfbins_nc1_60$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1_60$bin[1:100],TZ4_halfbin_fit_nc1_60$fitted[1:100])
  
  TZ4_loess_resid_nc1_60 <- data.frame(TZ4_ratio_halfbins_nc1_60$counts,TZ4_halfbin_fit_nc1_60$residual) %>% rename(count=1,residual=2)
  TZ4_loess_resid_nc1_60_list[[i]] <- TZ4_loess_resid_nc1_60
  
  # local slope
  TZ4_loess_slope_nc1_60 <- data.frame(diff(TZ4_halfbin_fit_nc1_60$fitted)/diff(TZ4_ratio_halfbins_nc1_60$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1_60[2:nrow(TZ4_ratio_halfbins_nc1_60),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1_60[[i]] <- TZ4_loess_slope_nc1_60
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1_60<-TZ4_loess_slope_nc1_60 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1_60[[i]]<-TZ4_halfbin_cutoff_nc1_60
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1_60<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1_60)
  TZ_4_nc1_60<-ContourLines2SLDF(TZcontour_4_nc1_60)
  
  TZ_nc1_60[[i]]<-TZ_4_nc1_60
  
  ############## span = 70 bins ###############
  TZratio4_nc1_70<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1_70<-hist(TZratio4_nc1_70$ratio_transform,breaks=seq(0,max(TZratio4_nc1_70$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1_70<-data.frame(TZratio4_ratio_halfbins_nc1_70$mids)
  TZ4_ratio_halfbins_nc1_70$counts<-TZratio4_ratio_halfbins_nc1_70$counts
  colnames(TZ4_ratio_halfbins_nc1_70)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1_70<-slice(TZ4_ratio_halfbins_nc1_70,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1_70 <- loess(TZ4_ratio_halfbins_nc1_70$counts~TZ4_ratio_halfbins_nc1_70$bin,span = (70/nrow(TZ4_ratio_halfbins_nc1_70)))
  plot(TZ4_ratio_halfbins_nc1_70$bin[1:100], TZ4_ratio_halfbins_nc1_70$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1_70$bin[1:100],TZ4_halfbin_fit_nc1_70$fitted[1:100])
  
  TZ4_loess_resid_nc1_70 <- data.frame(TZ4_ratio_halfbins_nc1_70$counts,TZ4_halfbin_fit_nc1_70$residual) %>% rename(count=1,residual=2)
  TZ4_loess_resid_nc1_70_list[[i]] <- TZ4_loess_resid_nc1_70
  
  # local slope
  TZ4_loess_slope_nc1_70 <- data.frame(diff(TZ4_halfbin_fit_nc1_70$fitted)/diff(TZ4_ratio_halfbins_nc1_70$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1_70[2:nrow(TZ4_ratio_halfbins_nc1_70),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1_70[[i]] <- TZ4_loess_slope_nc1_70
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1_70<-TZ4_loess_slope_nc1_70 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1_70[[i]]<-TZ4_halfbin_cutoff_nc1_70
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1_70<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1_70)
  TZ_4_nc1_70<-ContourLines2SLDF(TZcontour_4_nc1_70)
  
  TZ_nc1_70[[i]]<-TZ_4_nc1_70
  
  ############## span = 80 bins ###############
  TZratio4_nc1_80<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1_80<-hist(TZratio4_nc1_80$ratio_transform,breaks=seq(0,max(TZratio4_nc1_80$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1_80<-data.frame(TZratio4_ratio_halfbins_nc1_80$mids)
  TZ4_ratio_halfbins_nc1_80$counts<-TZratio4_ratio_halfbins_nc1_80$counts
  colnames(TZ4_ratio_halfbins_nc1_80)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1_80<-slice(TZ4_ratio_halfbins_nc1_80,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1_80 <- loess(TZ4_ratio_halfbins_nc1_80$counts~TZ4_ratio_halfbins_nc1_80$bin,span = (80/nrow(TZ4_ratio_halfbins_nc1_80)))
  plot(TZ4_ratio_halfbins_nc1_80$bin[1:100], TZ4_ratio_halfbins_nc1_80$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1_80$bin[1:100],TZ4_halfbin_fit_nc1_80$fitted[1:100])
  
  TZ4_loess_resid_nc1_80 <- data.frame(TZ4_ratio_halfbins_nc1_80$counts,TZ4_halfbin_fit_nc1_80$residual) %>% rename(count=1,residual=2)
  TZ4_loess_resid_nc1_80_list[[i]] <- TZ4_loess_resid_nc1_80
  
  # local slope
  TZ4_loess_slope_nc1_80 <- data.frame(diff(TZ4_halfbin_fit_nc1_80$fitted)/diff(TZ4_ratio_halfbins_nc1_80$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1_80[2:nrow(TZ4_ratio_halfbins_nc1_80),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1_80[[i]] <- TZ4_loess_slope_nc1_80
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1_80<-TZ4_loess_slope_nc1_80 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1_80[[i]]<-TZ4_halfbin_cutoff_nc1_80
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1_80<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1_80)
  TZ_4_nc1_80<-ContourLines2SLDF(TZcontour_4_nc1_80)
  
  TZ_nc1_80[[i]]<-TZ_4_nc1_80
  
  ############## span = 90 bins ###############
  TZratio4_nc1_90<-TZratio4@data %>% filter(ratio_transform<=20000)
  
  # get data ready
  TZratio4_ratio_halfbins_nc1_90<-hist(TZratio4_nc1_90$ratio_transform,breaks=seq(0,max(TZratio4_nc1_90$ratio_transform)+0.1,by=0.1),plot=FALSE)
  TZ4_ratio_halfbins_nc1_90<-data.frame(TZratio4_ratio_halfbins_nc1_90$mids)
  TZ4_ratio_halfbins_nc1_90$counts<-TZratio4_ratio_halfbins_nc1_90$counts
  colnames(TZ4_ratio_halfbins_nc1_90)<-c("bin","counts")
  TZ4_ratio_halfbins_nc1_90<-slice(TZ4_ratio_halfbins_nc1_90,-c(1:10))
  
  # change smoothing parameter
  TZ4_halfbin_fit_nc1_90 <- loess(TZ4_ratio_halfbins_nc1_90$counts~TZ4_ratio_halfbins_nc1_90$bin,span = (90/nrow(TZ4_ratio_halfbins_nc1_90)))
  plot(TZ4_ratio_halfbins_nc1_90$bin[1:100], TZ4_ratio_halfbins_nc1_90$counts[1:100])
  lines(TZ4_ratio_halfbins_nc1_90$bin[1:100],TZ4_halfbin_fit_nc1_90$fitted[1:100])
  
  TZ4_loess_resid_nc1_90 <- data.frame(TZ4_ratio_halfbins_nc1_90$counts,TZ4_halfbin_fit_nc1_90$residual) %>% rename(count=1,residual=2)
  TZ4_loess_resid_nc1_90_list[[i]] <- TZ4_loess_resid_nc1_90
  
  # local slope
  TZ4_loess_slope_nc1_90 <- data.frame(diff(TZ4_halfbin_fit_nc1_90$fitted)/diff(TZ4_ratio_halfbins_nc1_90$bin)) %>% bind_cols(TZ4_ratio_halfbins_nc1_90[2:nrow(TZ4_ratio_halfbins_nc1_90),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess_nc1_90[[i]] <- TZ4_loess_slope_nc1_90
  
  # NOTE: where pred_diff_bin_diff > -1 is the "inflection point" chosen ratio for TZ
  TZ4_halfbin_cutoff_nc1_90<-TZ4_loess_slope_nc1_90 %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff_nc1_90[[i]]<-TZ4_halfbin_cutoff_nc1_90
  
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4_nc1_90<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff_nc1_90)
  TZ_4_nc1_90<-ContourLines2SLDF(TZcontour_4_nc1_90)
  
  TZ_nc1_90[[i]]<-TZ_4_nc1_90
}


m_nc1_20<-do.call(bind,TZ_nc1_20)
m_nc1_30<-do.call(bind,TZ_nc1_30)
m_nc1_40<-do.call(bind,TZ_nc1_40)
m_nc1_50<-do.call(bind,TZ_nc1_50)
m_nc1_60<-do.call(bind,TZ_nc1_60)
m_nc1_70<-do.call(bind,TZ_nc1_70)
m_nc1_80<-do.call(bind,TZ_nc1_80)
m_nc1_90<-do.call(bind,TZ_nc1_90)


TZ4_loess_nc1_20_df <- do.call(bind,TZ4_loess_nc1_20)
TZ4_loess_nc1_30_df <- do.call(bind,TZ4_loess_nc1_30)
TZ4_loess_nc1_40_df <- do.call(bind,TZ4_loess_nc1_40)
TZ4_loess_nc1_50_df <- do.call(bind,TZ4_loess_nc1_50)
TZ4_loess_nc1_60_df <- do.call(bind,TZ4_loess_nc1_60)
TZ4_loess_nc1_70_df <- do.call(bind,TZ4_loess_nc1_70)
TZ4_loess_nc1_80_df <- do.call(bind,TZ4_loess_nc1_80)
TZ4_loess_nc1_90_df <- do.call(bind,TZ4_loess_nc1_90)

low_nc1_20<-TZ4_loess_nc1_20_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))
low_nc1_30<-TZ4_loess_nc1_30_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))
low_nc1_40<-TZ4_loess_nc1_40_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))
low_nc1_50<-TZ4_loess_nc1_50_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))
low_nc1_60<-TZ4_loess_nc1_60_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))
low_nc1_70<-TZ4_loess_nc1_70_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))
low_nc1_80<-TZ4_loess_nc1_80_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))
low_nc1_90<-TZ4_loess_nc1_90_df %>% filter(counts<25) %>% group_by(iteration) %>% summarise(bin_min=min(bin))


## Relative ratio of northern and southern trees in Wisconsin
# take mean of all ratio maps
NS_ratiotrans_table<-matrix(unlist(NS_ratiotrans),ncol=niter)
NS_ratiotrans_mean<-rowMeans(NS_ratiotrans_table)

# make new SGDF with mean ratiotrans
TZratio0<-as(North.dens4,"SpatialGridDataFrame")
names(TZratio0)<-"North"
TZratio0$South<-as(South.dens4,"SpatialGridDataFrame")$v
TZratio0$mean_ratiotrans<-NS_ratiotrans_mean

# convert to spatial pixel data frame
TZratio_mean<-as(TZratio0,"SpatialPixelsDataFrame")

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
iterline<-list("sp.lines",m_nc1_20,a=0.2) # set formatting for the 50 lines

TZratio_mean@data$cut_mean <- cut(TZratio_mean@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
Fig_nc1_20 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(iterline,wisc_HARN),
                    col.regions=rev(brewer.pal(11,"RdYlBu")),
                    colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                    par.settings = list(axis.line = list(col =  'transparent')),
                    auto.key = list(title = "Relative ratio"),
                    main="Cutoff = 20000, bin = 0.1, span = 20")

Fig_nc1_20

iterline<-list("sp.lines",m_nc1_50,a=0.2) # set formatting for the 50 lines

Fig_nc1_50 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(iterline,wisc_HARN),
                    col.regions=rev(brewer.pal(11,"RdYlBu")),
                    colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                    par.settings = list(axis.line = list(col =  'transparent')),
                    auto.key = list(title = "Relative ratio"),
                    main="Cutoff = 20000, bin = 0.1, span = 50")

Fig_nc1_50

## ratio summaries
TZ_cutoff_nc1_20_ratio <- data.frame(unlist(TZ4_cutoff_nc1_20)) %>% mutate(span=20) %>% rename(cp=1)
TZ_cutoff_nc1_30_ratio <- data.frame(unlist(TZ4_cutoff_nc1_30)) %>% mutate(span=30) %>% rename(cp=1)
TZ_cutoff_nc1_40_ratio <- data.frame(unlist(TZ4_cutoff_nc1_40)) %>% mutate(span=40) %>% rename(cp=1)
TZ_cutoff_nc1_50_ratio <- data.frame(unlist(TZ4_cutoff_nc1_50)) %>% mutate(span=50) %>% rename(cp=1)
TZ_cutoff_nc1_60_ratio <- data.frame(unlist(TZ4_cutoff_nc1_60)) %>% mutate(span=60) %>% rename(cp=1)
TZ_cutoff_nc1_70_ratio <- data.frame(unlist(TZ4_cutoff_nc1_70)) %>% mutate(span=70) %>% rename(cp=1)
TZ_cutoff_nc1_80_ratio <- data.frame(unlist(TZ4_cutoff_nc1_80)) %>% mutate(span=80) %>% rename(cp=1)
TZ_cutoff_nc1_90_ratio <- data.frame(unlist(TZ4_cutoff_nc1_90)) %>% mutate(span=90) %>% rename(cp=1)


TZ_cutoff_ratio_all <- bind_rows(TZ_cutoff_nc1_20_ratio,TZ_cutoff_nc1_30_ratio,TZ_cutoff_nc1_40_ratio,TZ_cutoff_nc1_50_ratio,
                                 TZ_cutoff_nc1_60_ratio,TZ_cutoff_nc1_70_ratio,TZ_cutoff_nc1_80_ratio,TZ_cutoff_nc1_90_ratio,) %>% mutate(span=as.factor(span))

ggplot(data=TZ_cutoff_ratio_all,aes(x=span,y=cp))+
  geom_boxplot()


#################### boxplot ##############################

## Average ratio line
avg_cutoff_nc1_20<-mean(unlist(TZ4_cutoff_nc1_20))
sd_cutoff_nc1_20<-sd(unlist(TZ4_cutoff_nc1_20))

avg_cutoff_nc1_30<-mean(unlist(TZ4_cutoff_nc1_30))
sd_cutoff_nc1_30<-sd(unlist(TZ4_cutoff_nc1_30))

avg_cutoff_nc1_40<-mean(unlist(TZ4_cutoff_nc1_40))
sd_cutoff_nc1_40<-sd(unlist(TZ4_cutoff_nc1_40))

avg_cutoff_nc1_50<-mean(unlist(TZ4_cutoff_nc1_50))
sd_cutoff_nc1_50<-sd(unlist(TZ4_cutoff_nc1_50))

avg_cutoff_nc1_60<-mean(unlist(TZ4_cutoff_nc1_60))
sd_cutoff_nc1_60<-sd(unlist(TZ4_cutoff_nc1_60))

avg_cutoff_nc1_70<-mean(unlist(TZ4_cutoff_nc1_70))
sd_cutoff_nc1_70<-sd(unlist(TZ4_cutoff_nc1_70))

avg_cutoff_nc1_80<-mean(unlist(TZ4_cutoff_nc1_80))
sd_cutoff_nc1_80<-sd(unlist(TZ4_cutoff_nc1_80))

avg_cutoff_nc1_90<-mean(unlist(TZ4_cutoff_nc1_90))
sd_cutoff_nc1_90<-sd(unlist(TZ4_cutoff_nc1_90))

# map of average cutoff on average ratio map
TZratio4_mean_image_nc1_20<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1_20<-contourLines(TZratio4_mean_image_nc1_20,levels=avg_cutoff_nc1_20)
TZ_mean_4_nc1_20<-ContourLines2SLDF(TZcontour_mean_4_nc1_20)

TZratio4_mean_image_nc1_30<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1_30<-contourLines(TZratio4_mean_image_nc1_30,levels=avg_cutoff_nc1_30)
TZ_mean_4_nc1_30<-ContourLines2SLDF(TZcontour_mean_4_nc1_30)

TZratio4_mean_image_nc1_40<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1_40<-contourLines(TZratio4_mean_image_nc1_40,levels=avg_cutoff_nc1_40)
TZ_mean_4_nc1_40<-ContourLines2SLDF(TZcontour_mean_4_nc1_40)

TZratio4_mean_image_nc1_50<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1_50<-contourLines(TZratio4_mean_image_nc1_50,levels=avg_cutoff_nc1_50)
TZ_mean_4_nc1_50<-ContourLines2SLDF(TZcontour_mean_4_nc1_50)

TZratio4_mean_image_nc1_60<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1_60<-contourLines(TZratio4_mean_image_nc1_60,levels=avg_cutoff_nc1_60)
TZ_mean_4_nc1_60<-ContourLines2SLDF(TZcontour_mean_4_nc1_60)

TZratio4_mean_image_nc1_70<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1_70<-contourLines(TZratio4_mean_image_nc1_70,levels=avg_cutoff_nc1_70)
TZ_mean_4_nc1_70<-ContourLines2SLDF(TZcontour_mean_4_nc1_70)

TZratio4_mean_image_nc1_80<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1_80<-contourLines(TZratio4_mean_image_nc1_80,levels=avg_cutoff_nc1_80)
TZ_mean_4_nc1_80<-ContourLines2SLDF(TZcontour_mean_4_nc1_80)

TZratio4_mean_image_nc1_90<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4_nc1_90<-contourLines(TZratio4_mean_image_nc1_90,levels=avg_cutoff_nc1_90)
TZ_mean_4_nc1_90<-ContourLines2SLDF(TZcontour_mean_4_nc1_90)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean@data$cut_mean <- cut(TZratio_mean@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")

mean_plot_nc1_20 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1_20,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 20,000, span = 20")

mean_plot_nc1_30 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1_30,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 20,000, span = 30")

mean_plot_nc1_40 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1_40,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 20,000, span = 40")

mean_plot_nc1_50 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1_50,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 20,000, span = 50")

mean_plot_nc1_60 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1_60,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 20,000, span = 60")

mean_plot_nc1_70 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1_70,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 20,000, span = 70")

mean_plot_nc1_80 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1_80,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 20,000, span = 80")

mean_plot_nc1_90 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4_nc1_90,wisc_HARN),
                          col.regions=rev(brewer.pal(11,"RdYlBu")),
                          colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                          par.settings = list(axis.line = list(col =  'transparent')),
                          auto.key = list(title = "Relative ratio"),
                          main = "bin size = 0.1, trunction point = 20,000, span = 90")

mean_plot_nc1_20
mean_plot_nc1_30
mean_plot_nc1_40
mean_plot_nc1_50
mean_plot_nc1_60
mean_plot_nc1_70
mean_plot_nc1_80
mean_plot_nc1_90


### plot residuals
plot(TZ4_loess_resid_nc1_20_list[[4]][,1],TZ4_loess_resid_nc1_20_list[[4]][,2])
plot(TZ4_loess_resid_nc1_30_list[[4]][,1],TZ4_loess_resid_nc1_30_list[[4]][,2])
plot(TZ4_loess_resid_nc1_40_list[[4]][,1],TZ4_loess_resid_nc1_40_list[[4]][,2])
plot(TZ4_loess_resid_nc1_50_list[[4]][,1],TZ4_loess_resid_nc1_50_list[[4]][,2])
plot(TZ4_loess_resid_nc1_60_list[[4]][,1],TZ4_loess_resid_nc1_60_list[[4]][,2])
plot(TZ4_loess_resid_nc1_70_list[[4]][,1],TZ4_loess_resid_nc1_70_list[[4]][,2])
plot(TZ4_loess_resid_nc1_80_list[[4]][,1],TZ4_loess_resid_nc1_80_list[[4]][,2])
plot(TZ4_loess_resid_nc1_90_list[[4]][,1],TZ4_loess_resid_nc1_90_list[[4]][,2])

### plot lines together
plot(wisc_HARN)
plot(TZ_mean_4_nc1_20, col="gray80",add=TRUE)
plot(TZ_mean_4_nc1_30, col= "orange",add=TRUE)
plot(TZ_mean_4_nc1_40,col="blue",add=TRUE)
plot(TZ_mean_4_nc1_50,col="red",add=TRUE)
plot(TZ_mean_4_nc1_60,col="green",add=TRUE)
plot(TZ_mean_4_nc1_70,col="purple",add=TRUE)
#plot(TZ_mean_4_nc1_80,col="gray",add=TRUE)
#plot(TZ_mean_4_nc1_90,col="purple",add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

plot(wisc_HARN)
plot(TZ_mean_4_nc1_20, col="gray85",add=TRUE)
plot(TZ_mean_4_nc1_30, col= "gray70",add=TRUE)
plot(TZ_mean_4_nc1_40,col="gray55",add=TRUE)
plot(TZ_mean_4_nc1_60,col="gray45",add=TRUE)
plot(TZ_mean_4_nc1_70,col="gray25",add=TRUE)
#plot(TZ_mean_4_nc1_80,col="gray",add=TRUE)
#plot(TZ_mean_4_nc1_90,col="purple",add=TRUE)

plot(TZ_mean_4_nc1_50,col="blue",lwd=2,add=TRUE)
scalebar(d= 100000, type='bar',divs=4,below = 'm')

#### save image #####
### convert to sf
wisc_HARN_sf <- as(wisc_HARN,"sf")

## combine lines
# add field with variable
TZ_mean_4_nc1_20$smooth <- "20"
TZ_mean_4_nc1_30$smooth <- "30"
TZ_mean_4_nc1_40$smooth <- "40"
TZ_mean_4_nc1_50$smooth <- "50"
TZ_mean_4_nc1_60$smooth <- "60"
TZ_mean_4_nc1_70$smooth <- "70"

# combine FS lines
TZ_mean_4_nc1_all <- rbind(TZ_mean_4_nc1_20,TZ_mean_4_nc1_30,TZ_mean_4_nc1_40,TZ_mean_4_nc1_50,TZ_mean_4_nc1_60,TZ_mean_4_nc1_70)
crs(TZ_mean_4_nc1_all) <- "+init=epsg:3070"
TZ_mean_4_nc1_all$smooth <- factor(TZ_mean_4_nc1_all$smooth,levels = c("20","30","40","50","60","70"))
TZ_mean_4_nc1_all_sf <- as(TZ_mean_4_nc1_all,"sf")

### plot
FigS6 <- ggplot(data=wisc_HARN_sf)+
  geom_sf(fill="white")+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true",
                         height = unit(1, "cm"), #how tall the arrow should be
                         width= unit(0.5, "cm"), 
                         pad_x = unit(0.25, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_orienteering)+
  geom_sf(data=TZ_mean_4_nc1_all_sf,aes(col=smooth,fill=smooth))+
  scale_color_manual(breaks=c("20","30","40","50","60","70"),
                     values=alpha(c("20"="gray85","30"="gray70","40"="gray55","50"="blue","60"="gray45","70"="gray25"),0.8),
                     labels=c("0.0010","0.0015","0.0020","0.0025","0.0030","0.0035"),guide="none")+
  scale_fill_manual(values=c("20"="gray85","30"="gray70","40"="gray55","50"="blue","60"="gray45","70"="gray25"),
                    labels=c("0.0010","0.0015","0.0020","0.0025","0.0030","0.0035"))+
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
ggsave("./Output/FigS6.jpg", plot=FigS6, width = 4000/dpi, height = 4000/dpi, dpi = dpi)

