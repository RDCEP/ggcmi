##################################################################
#                                                                #
#   GGCMI_PHASE2_MINIMUM_CROPLAND_MASK.R                         #
#   R-Script to generate harmonized growing areas to be used     #
#   at least for GGCMI phase 2 simulations.                      #
#   Data is based on                                             #
#   * GAEZ 3.0 suitability indeces as prepared by Tom Pugh       #
#   * MIRCA2000 total cropland extent                            #
#   excluding all areas that are unsuitable according to         #
#   GAEZ but only if there is no cropland according              #
#   to MIRCA2000. A minimum threshold for MIRCA2000 is used.     #
#                                                                #
#   written by Christoph Mueller, PIK                            #
#   2015/08/20                                                   #
#                                                                #
#   modified by Jannes Breier                                    #
#   2015/10/21                                                   #
#                                                                # 
##################################################################


### settings ####

# deletes all variables
rm(list=ls(all=TRUE))

# load required libaries
require(fields)
require(ncdf4)
require(maps)

# data location
path.mask <- "/p/projects/macmit/data/GGCMI/phase2_masks/"
path.mirca <- "/p/projects/macmit/data/GGCMI/AgMIP.output/processed.150505/masks/weight/"

testing <- F

### functions ####

read.cropland <- function(fname,va="sum"){
  nf <- nc_open(fname)    																	
  data <- ncvar_get(nf,varid=va) # reading sum of all cropland??? or just of the major 4?
  nc_close(nf)
  data
}

read.total.cropland <- function(fname,va="sum"){
  nf <- nc_open(fname)  
  vars <- names(nf$var)
  data <- ncvar_get(nf,varid=vars[1])
  data[is.finite(data)] <- 0 #reset
  for(var in vars[-which(vars=="Managed_grass")]) data <- data + ncvar_get(nf,varid=var)
  nc_close(nf)
  data
}

read.suitability.mask <- function(fname,va="suitability"){
  nf <- nc_open(fname)    																	
  data <- ncvar_get(nf,varid=va) # reading boolean suitability mask (0,1)
  nc_close(nf)
  data  
}

read.landmask <- function(fname,va="mask"){
  nf <- nc_open(fname)      																
  data <- ncvar_get(nf,varid=va) # reading boolean suitability mask (0,1)
  nc_close(nf)
  data  
}

### main ####

# all.nc4 does not seem to have all data
cropland <- read.cropland(paste(path.mirca,"all.nc4",sep=""))
cropland.m <- read.cropland(paste(path.mirca,"all.nc4",sep=""),"sum_mai")
cropland.w <- read.cropland(paste(path.mirca,"all.nc4",sep=""),"sum_whe")
cropland.r <- read.cropland(paste(path.mirca,"all.nc4",sep=""),"sum_ric")
cropland.s <- read.cropland(paste(path.mirca,"all.nc4",sep=""),"sum_soy")

cropland.rf <- read.total.cropland(paste(path.mirca,"landuse.rf.nc4",sep=""))
cropland.ir <- read.total.cropland(paste(path.mirca,"landuse.ir.nc4",sep=""))
total.cropland <- cropland.rf + cropland.ir


# read 90% severe constraint masks as prepared by Tom
non.suitable <- read.suitability.mask(paste(path.mask,"gaez_soil_suitability_mask_90percent.nc4",sep=""))[,360:1]

boolean.mask <- non.suitable
boolean.mask[non.suitable==1] <- 0
boolean.mask[non.suitable==0] <- 1
boolean.mask[total.cropland>0] <- 1 #>0.001 seemed to be ineffective, so >0 was set

landmask <- read.landmask(paste(path.mirca,"../aggr/binary.mask.nc4",sep="")) 



### write NC file ####

require(ncdf4)

# missing value
mv <- 1.e20 

#write dim-vector
dim_lon <- ncdim_def("lon","degrees_east",seq(-179.75,179.75,len=360/0.5))
dim_lat <- ncdim_def("lat","degrees_north",seq(-89.75,89.75,len=180/0.5))


#define variables including the predefined dimensions
ncm <- ncvar_def("minimum cropland mask"," boolean",list(dim_lon,dim_lat),mv,
                longname="minimum cropland to be simulated in GGCMI phase 2 simulations",
                compression=9)

# create files
nf <- nc_create(paste0(path.mask,"boolean_cropmask_ggcmi_phase2.nc4"), ncm)

#
# space for further insertions of titles/comments for (global) variables
#

# adding variables to files
ncatt_put(nf,varid=0,"title","minimum cropland to be simulated in GGCMI phase 2 simulations")
ncatt_put(nf,varid=0,"author","Christoph Mueller, PIK cmueller@pik-potsdam.de")
ncatt_put(nf,varid=0,"comment 1","use this mask to mask out unsuitable areas from simulations for AgMIP GGCMI phase2")
ncatt_put(nf,varid=0,"comment 2","unsuitability is based on GAEZ; Grid-cells are marked as unsitable if they have at least 90% of their area classified as at least severe constraints (categories 3-6). Variables considered are: Rooting conditions, oxygen availability, excess salts, toxicities, workability.")
ncatt_put(nf,varid=0,"comment 3","unsuitability is neglected, if pixel is used for cropland (any crop) according to the MIRCA2000 data set")

# preparing data for NC files
mapm <- boolean.mask[,360:1] # boolean.mask - invert latitudes 
mapm[is.na(mapm)] <- mv
mapm[!is.finite(landmask[,360:1])] <- mv

# writing data to NC files
ncvar_put(nf,ncm,mapm)

# closing NC files
nc_close(nf)

if(testing){
  ### plot for testing and comparision ###
  
  cropland[!is.finite(landmask)] <- NA
  total.cropland[!is.finite(landmask)] <- NA
  cropland.rf[!is.finite(landmask)] <- NA
  cropland.ir[!is.finite(landmask)] <- NA
  boolean.mask[!is.finite(landmask)] <- NA
  non.suitable[!is.finite(landmask)] <- NA
  cropland.rf[cropland.rf==0] <- NA
  cropland.ir[cropland.ir==0] <- NA
  
  
  ## comparision between the sum of different data sources ("all.nc4" vs. "landuse.rf.nc4")
  
  png(paste0(path.mask,"testing.png"),width=2*300,height=3*300,res=300,pointsize=6)
  par(oma=c(0,0,0,0))
  split.screen(c(2,1))
  
  screen(1)
  par(mar=c(0,0,0,0))
  image(x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),
        cropland[,360:1]) 
  map(add=T,boundary=T,lty=0,col="grey60")
  
  
  screen(2)
  par(mar=c(0,0,0,0))
  image(x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),
        total.cropland[,360:1])
  map(add=T,boundary=T,lty=2,col="grey60")
  close.screen(all=T)
  dev.off()
  
  ## comparision between the consecutive steps until the final boolean mask
  
  png(paste0(path.mask,"testing2.png"),width=5*300,height=3*300,res=300,pointsize=6)
  par(oma=c(0,0,0,0))
  
  split.screen(c(2,2))
  
  screen(1)
  par(mar=c(0,0,0,0))
  image(x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),
        cropland.rf[,360:1]) 
  map(add=T,boundary=T,lty=0,lwd=0.001,resolution=1,col="grey60")
  
  
  screen(2)
  par(mar=c(0,0,0,0))
  image(x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),
        cropland.ir[,360:1])
  map(add=T,boundary=T,lwd=0.01,col="grey60")
  
  
  screen(3)
  par(mar=c(0,0,0,0))
  image(x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),
        non.suitable[,360:1])
  map(add=T,boundary=T,lwd=0.01,col="grey60")
  
  
  screen(4)
  par(mar=c(0,0,0,0))
  image(x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),
        boolean.mask[,360:1])
  map(add=T,boundary=T,lwd=0.01,col="grey60")
  
  close.screen(all=T)  
  dev.off()
}
