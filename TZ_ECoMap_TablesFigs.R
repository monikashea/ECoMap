# set your own working directory
setwd("U:/Shea/Research/GIS/Projects/TZ_Boundary/github")

library(indicspecies)
library(maptools)
library(spatstat)
library(sp)
library(rgdal)
library(ggplot2)
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


#### mapping ####
## number of iterations; can change, but doing at least 10 is suggested
niter<-50

## create lists
IndSppList<-vector("list",niter)
NS_ratio <- vector("list",niter)
NS_logratio <- vector("list",niter)
NS_ratiotrans <- vector("list",niter)
TZ4_loess <- vector("list",niter)
TZ4_cutoff <- vector("list",niter)
TZ <- vector("list",niter)

## loop conducts multiple iterations
## this step can take a while to run

set.seed(6542)
for (i in 1:niter){
  #### ISA
  # 50 is number of samples; can adjust
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
  
  # set bandwidth; can adjust
  bw.tz4<-4000
  
  ### calculate intensity
  # eps indicates pixel size. 1600m was chosen because this matches how the data were collected, in ~1600 grid cells
  # adjust eps as appropriate for data
  North.dens4 <- density(North.ppp,bw.tz4,eps=1600)
  South.dens4 <- density(South.ppp,bw.tz4,eps=1600)
  
  ### take ratio
  # processing
  TZratio0<-as(North.dens4,"SpatialGridDataFrame")
  names(TZratio0)<-"North"
  TZratio0$South<-as(South.dens4,"SpatialGridDataFrame")$v
  TZratio4<-as(TZratio0,"SpatialPixelsDataFrame")
  
  # ratio
  TZratio4$TZratio<-TZratio4$North/TZratio4$South
  TZratio4$logratio<-log(TZratio4$TZratio)
  
  # relative ratio
  TZratio4$ratio_transform<-ifelse(TZratio4$logratio<0,1/TZratio4$TZratio,TZratio4$TZratio)
  
  # store ratio, logratio, and transformed ratio from each run
  TZratioGrid <- eval.im(North.dens4/South.dens4)
  TZlogratioGrid <- eval.im(log(North.dens4/South.dens4))
  TZratiotransGrid <- eval.im(ifelse(TZlogratioGrid<0,1/TZratioGrid,TZratioGrid))
  NS_ratio[[i]] <-as(TZratioGrid,"SpatialGridDataFrame")$v
  NS_logratio[[i]] <-as(TZlogratioGrid,"SpatialGridDataFrame")$v
  NS_ratiotrans[[i]] <-as(TZratiotransGrid,"SpatialGridDataFrame")$v
  
  ### create histogram for loess curve to determine critical point
  ## getting raw data from histogram bins
  ## binwidth = 0.1, truncation point = 20000, span = 50 bins; these parameters can be adjusted
  TZratio4_nc1<-TZratio4@data %>% filter(ratio_transform<=20000) # ratio_transform<=20000 sets truncation point, can be adjusted
  
  # get data ready
  TZratio4_ratio_halfbins<-hist(TZratio4_nc1$ratio_transform,breaks=seq(0,max(TZratio4_nc1$ratio_transform)+0.1,by=0.1),plot=FALSE) # by=0.1 sets bin size, can be adjusted
  TZ4_ratio_halfbins<-data.frame(TZratio4_ratio_halfbins$mids)
  TZ4_ratio_halfbins$counts<-TZratio4_ratio_halfbins$counts
  colnames(TZ4_ratio_halfbins)<-c("bin","counts")
  TZ4_ratio_halfbins<-slice(TZ4_ratio_halfbins,-c(1:10))
  
  # use loess to fit a line
  TZ4_halfbin_fit <- loess(TZ4_ratio_halfbins$counts~TZ4_ratio_halfbins$bin,span = (50/nrow(TZ4_ratio_halfbins))) # 50 is the number of bins being smoothed, formula sets the smoothing parameter; can be adjusted
  plot(TZ4_ratio_halfbins$bin[1:100], TZ4_ratio_halfbins$counts[1:100])
  lines(TZ4_ratio_halfbins$bin[1:100],TZ4_halfbin_fit$fitted[1:100])
  
  # local slope
  TZ4_loess_slope <- data.frame(diff(TZ4_halfbin_fit$fitted)/diff(TZ4_ratio_halfbins$bin)) %>% bind_cols(TZ4_ratio_halfbins[2:nrow(TZ4_ratio_halfbins),]) %>% rename(slope=1) %>% mutate(iteration=i)
  TZ4_loess[[i]] <- TZ4_loess_slope
  
  # NOTE: where pred_diff_bin_diff > -1 is the "leveling off point" critical point ratio for TZ margins
  TZ4_halfbin_cutoff<-TZ4_loess_slope %>% filter(slope>-1) %>% summarise(ratio=min(bin))
  TZ4_cutoff[[i]]<-TZ4_halfbin_cutoff
  
  # draw contour line on relative ratio map
  TZratio4_image<-as.image.SpatialGridDataFrame(TZratio4["ratio_transform"])
  TZcontour_4<-contourLines(TZratio4_image,levels=TZ4_halfbin_cutoff)
  TZ_4<-ContourLines2SLDF(TZcontour_4)
  
  # save contour line
  TZ[[i]]<-TZ_4
}

# ignore warnings

# compile contour lines from each permutation
m<-do.call(bind,TZ)


########### Table 1 #############
## Percent of fifty runs that each species is a significant indicator species for the north or south
# summarize indicator species
IndSpp_TZ4<-do.call("rbind",IndSppList)
NSpp_TZ4<-IndSpp_TZ4 %>% filter(sign_N_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="North")
SSpp_TZ4<-IndSpp_TZ4 %>% filter(sign_S_sig==1) %>% group_by(SP_new) %>% summarise(count=n()) %>% mutate(Zone="South")
IndSpp_all_TZ4 <- rbind(NSpp_TZ4,SSpp_TZ4)
IndSpp_NSB_TZ4 <- IndSpp_TZ4 %>% group_by(SP_new) %>% summarise(N=sum(sign_N_sig),S=sum(sign_S_sig),B=sum(sign_both))

# export
#write.csv(NSpp_TZ4,file="./Output/NSpp_TZ4_FS.csv",row.names = FALSE)
#write.csv(SSpp_TZ4,file="./Output/SSpp_TZ4_FS.csv",row.names = FALSE)
write.csv(IndSpp_NSB_TZ4,file="./Output/IndSpp_NSB_TZ4_FS.csv",row.names = FALSE)

## Table 1 was created in Excel based on the files exported above

########### Figure 1a #############
## USDA Forest Service Ecological Provinces in Wisconsin, overlaying a map of Wisconsin county boundaries. 
# Mapped in ArcGIS

########### Figure 1b ##############
## density maps for north and south
plot(North.dens4,zlim=c(5.837e-20,3.175e-06))
plot(South.dens4,zlim=c(5.837e-20,3.175e-06))

jpeg(filename = "./Output/Fig1b_Ndens_2.jpg",width=6,height=6,units="in",res=600)
plot(North.dens4,zlim=c(5.837e-20,3.175e-06))
dev.off()

jpeg(filename = "./Output/Fig1b_Sdens_2.jpg",width=6,height=6,units="in",res=600)
plot(South.dens4,zlim=c(5.837e-20,3.175e-06))
dev.off()

########### Figure 1c #############
## relative ratio plot for one iteration
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")

TZratio4_fig <- TZratio4
TZratio4_fig@data$cut <- cut(TZratio4_fig@data$ratio_transform,breaks=brks)
wisc_map <- list("sp.polygons",wisc_HARN,col="black",alpha=0.5,first=F)

Fig1c <- spplot(TZratio4_fig,"cut", sp.layout=wisc_map,
       col.regions=rev(brewer.pal(11,"RdYlBu")),
       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
       par.settings = list(axis.line = list(col =  'transparent')),
       auto.key = list(title = "Relative ratio"))

Fig1c

jpeg(filename = "./Output/Fig1c_relratio.jpg",width=6,height=6,units="in",res=600)
print(Fig1c)
dev.off()


########### Figure 1d #############
## Example histogram with critical point
TZratio4_nc1_short <- TZratio4_nc1 %>% filter(ratio_transform<=100)
  
Fig1d_data <- data.frame(TZ4_halfbin_fit$fitted) %>% rename(fitted=1) %>% bind_cols(TZ4_ratio_halfbins) %>% filter(bin<=100)

Fig1d <- ggplot(TZratio4_nc1_short,aes(x=ratio_transform))+
  geom_histogram(binwidth = 0.1,color="grey",fill="white")+
  geom_line(data=Fig1d_data,aes(x=bin,y=fitted),linetype="dashed")+
  geom_vline(xintercept=TZ4_halfbin_cutoff$ratio, linetype="dotted")+
  scale_x_continuous(breaks=c(1,seq(25,75,by=25),100))+
  theme_classic()+
  xlab("Relative ratio")+
  ylab("Number of pixels")+
  theme(text=element_text(size=10, color = "black"),
        axis.title = element_text(size=12))

Fig1d

ggsave("Fig1d_hist.jpeg", plot = Fig1d, device = "jpeg", path = "./Output",
       width = 3, height = 2.5, units = "in", dpi=600)


########### Figure 1e #############
## Example rel ratio map with contour line
m_example <-list("sp.lines",TZ_4)

Fig1e <- spplot(TZratio4_fig,"cut", sp.layout=list(wisc_map,m_example),
                col.regions=rev(brewer.pal(11,"RdYlBu")),
                colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
                par.settings = list(axis.line = list(col =  'transparent')),
                auto.key = list(title = "Relative ratio"))

Fig1e

jpeg(filename = "./Output/Fig1e_map.jpg",width=6,height=6,units="in",res=600)
print(Fig1e)
dev.off()

########### Figure 2 #############
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
iterline<-list("sp.lines",m,a=0.2) # set formatting for the 50 lines
scale1 <- list("SpatialPolygonsRescale",layout.scale.bar(),offset = c(300000,230000), scale = 100000, fill=c("transparent","black"),labels = c("0","100 km"))
s1_text0 <- list("sp.text", c(300000,230000+15000), "0", cex = 1)
s1_text1 <- list("sp.text", c(300000+100000,230000+15000), "100 km", cex = 1)

TZratio_mean@data$cut_mean <- cut(TZratio_mean@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
Fig2 <-spplot(TZratio_mean,"cut_mean", sp.layout=list(iterline,wisc_HARN,scale1,s1_text0,s1_text1),
       col.regions=rev(brewer.pal(11,"RdYlBu")),
       colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
       par.settings = list(axis.line = list(col =  'transparent')),
       auto.key = list(title = "Relative ratio"))

Fig2

jpeg(filename = "./Output/Fig2_meanRR.jpg",width=6,height=6,units="in",res=600)
print(Fig2)
dev.off()

## NOTE: added north arrow manually (in paint)

########### Figure 3 #############
## Figure created in ArcGIS using average ratio line
## Average ratio line created below, exported to ArcGIS where islands and holes were removed.
avg_cutoff<-mean(unlist(TZ4_cutoff))
sd_cutoff<-sd(unlist(TZ4_cutoff))

# map of average cutoff on average ratio map
TZratio4_mean_image<-as.image.SpatialGridDataFrame(TZratio_mean["mean_ratiotrans"])
TZcontour_mean_4<-contourLines(TZratio4_mean_image,levels=avg_cutoff)
TZ_mean_4<-ContourLines2SLDF(TZcontour_mean_4)

# map relative ratio at range of levels
brks<-c(1,2,4,7,12,20,35,60,90,150,200,5e17) # set brakes so scale is readable
TZratio_mean@data$cut_mean <- cut(TZratio_mean@data$mean_ratiotrans,breaks=brks) # convert to factor to make a readable scale

brks_label <-c("1","2","4","7","12","20","35","60","90","150","200","> 200")
mean_plot <-spplot(TZratio_mean,"cut_mean", sp.layout=list(TZ_mean_4,wisc_HARN),
              col.regions=rev(brewer.pal(11,"RdYlBu")),
              colorkey=list(height = 0.5, labels = list(at = seq(0.5, length(brks) -0.5), labels = brks_label)), # colorkey sets label names
              par.settings = list(axis.line = list(col =  'transparent')),
              auto.key = list(title = "Relative ratio"))

mean_plot

mean_FS_TZ4<-writeOGR(TZ_mean_4,'./Output/shapefiles','mean_FS_TZ4_20191227',driver="ESRI Shapefile",overwrite_layer = TRUE)
lines_FS_TZ4<-writeOGR(m,'./Output/shapefiles','lines_FS_TZ4_20191227',driver="ESRI Shapefile",overwrite_layer = TRUE)

## Edit and finish in ArcGIS

############# Figure 4 ################
##### Figure 4 is in TZ_ECoMap_movezones #####


