#require(plotrix) # not needed as I'm using a tweaked copy of the taylor.diagram() function
rm(list=ls(all=TRUE))

# turn on warnings as they occur
options(warn=1)

require(ncdf4)
require(weights)

# script to plot Taylor diagrams for GGCMI phase 1. 
# written by Christoph Mueller, PIK
# cmueller@pik-potsdam.de
# I'm using a slightly modified version of taylor.diagram() of the plotrix package
# to avoid comparison with NA values and to allow for weighted correlation analysis.
# For that I've copied the source code of taylor.diagram to GGCMI_taylor_function.R and modified it accordingly.

# However, as I need to compute correlation coefficients of the shifted time series before plotting,
# I don't use the modified taylor.diagram funktion any more (except for the basic plot)
# but compute the coordinates from the correlation coefficient and the standard deviation
# colors adjusted for color-blindness and to avoid green and blue which is used for the contour lines

path <- "/p/projects/macmit/data/GGCMI/AgMIP.output/"
path.ref <- "/p/projects/macmit/data/GGCMI/reference/"
path.ref.detrend <- paste0(path.ref,"ray-dt/")
processed.ts <- "processed"

# read FPU names
fpu.names <- read.csv(paste0(path,processed.ts,"/masks/aggr/fpu.meta.csv"),header=F)[,2]
gadm0.names <- array("missing",253)
gadm0.index <- read.csv(paste0(path,processed.ts,"/masks/aggr/gadm0.meta.csv"),header=F)[,1]
name <- read.csv(paste0(path,processed.ts,"/masks/aggr/gadm0.meta.csv"),header=F)[,2]


source("/p/projects/macmit/users/cmueller/GGCMI/GGCMI_taylor_function.R")

get_nc4_data_slice <- function(fname,crd=cr0,dtd=dt0,mpd=mp0,scend=scen0){
  nf <- nc_open(fname)
  cr <- which(strsplit(ncatt_get(nf,varid="cr")$long_name,split=", ")[[1]]==crd)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  #fpu <- which(strsplit(ncatt_get(nf,varid="fpu")$long_name,split=", ")[[1]]==cr0)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  scen <- which(strsplit(ncatt_get(nf,varid="scen")$long_name,split=", ")[[1]]==scend)
  if(length(scen)==0){
    nc_close(nf)
    cat(scend,"not available.\n")
    return(NA)
  }
  # order is cr, mp, dt, scen, time, fpu
  # test if scen dimension is > 1, otherwise object has one dimension less
  if(length(strsplit(ncatt_get(nf,varid="scen")$long_name,split=", ")[[1]])>1) {
    data.bc <- ncvar_get(nf,varid="yield_detrend")[cr,mp,dt,scen,,]
  } else {
    data.bc <- ncvar_get(nf,varid="yield_detrend")[cr,mp,dt,,]
  }
  nc_close(nf)
  data.bc
}

get_nc4_ensemble_slice <- function(fname,crd=cr0,dtd=dt0,mpd=mp0,nmd=nm0,wtd=wt0){
  nf <- nc_open(fname)
  cr <- which(strsplit(ncatt_get(nf,varid="cr")$long_name,split=", ")[[1]]==crd)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  # always uses the best scen
  #scen <- which(strsplit(ncatt_get(nf,varid="top_scens")$long_name,split=", ")[[1]]==scend)
  wt <- which(strsplit(ncatt_get(nf,varid="wt")$long_name,split=", ")[[1]]==wtd)
  # order is cr, mp, dt, scen, time, fpu
  # test if scen dimension is > 1, otherwise object has one dimension less
  data.bc <- ncvar_get(nf,varid="yield_detrend")[wt,nmd,cr,mp,dt,,]
  nc_close(nf)
  data.bc
}

get_nc4_ref_slice <- function(fname,crop,dtd=dt0,mpd=mp0,getarea=F){
  nf <- nc_open(fname)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  # order is mp, dt, time, fpu
  if(getarea){
    if(crop=="mai") { data <- ncvar_get(nf,varid="area_mai")
    } else if(crop=="whe"){ data <- ncvar_get(nf,varid="area_whe")
    } else if(crop=="ric"){ data <- ncvar_get(nf,varid="area_ric")
    } else if(crop=="soy"){ data <- ncvar_get(nf,varid="area_soy")}
  } else{
    if(crop=="mai") { data <- ncvar_get(nf,varid="yield_mai")[mp,dt,,]
    } else if(crop=="whe"){ data <- ncvar_get(nf,varid="yield_whe")[mp,dt,,]
    } else if(crop=="ric"){ data <- ncvar_get(nf,varid="yield_ric")[mp,dt,,]
    } else if(crop=="soy"){ data <- ncvar_get(nf,varid="yield_soy")[mp,dt,,]}
  }
  nc_close(nf)
  data
}

adjust.ts <- function(data,ref){ # ref is extened one year in both directions
  if(length(which(!is.na(data*ref[c(1:length(data))+1])))<10) return(ref[c(1:length(data))+1]) # arbitrarily set a minimum of 10 years with data
  p0 <- wtd.cors(data,ref[c(1:length(data))+1])
  p1 <- wtd.cors(data,ref[c(1:length(data))])
  p2 <- wtd.cors(data,ref[c(1:length(data))+2])
  m.alt <- max(p1,p2,na.rm=T)
  w.alt <- which.max(c(p1,p0,p2))-1
  if(p0 < (m.alt+shift.threshold)){
    # w.alt is shifting data by 0 or 2, default (unshifted) would be a shift of 1
    return (ref[c(1:length(data))+w.alt])
  }
  ref[c(1:length(data))+1] # if alternatives are not better than standard ts, return standard ts
}

wtd.cor.mod <- function(aa,bb,cc){
  if(do.means){
    ncols <- dim(aa)[1]
    aa <- matrix(rep(colMeans(aa,na.rm=T),ncol=ncols),byrow=F)
    aa[!is.finite(aa)] <- NA
    ncols <- dim(bb)[1]
    bb <- matrix(rep(colMeans(bb,na.rm=T),ncol=ncols),byrow=F)
    bb[!is.finite(bb)] <- NA
    ncols <- dim(cc)[1]
    cc <- matrix(rep(colMeans(cc,na.rm=T),ncol=ncols),byrow=F)
    cc[!is.finite(cc)] <- NA
  }
  wtd.cor(as.vector(aa),as.vector(bb),as.vector(cc))
}


# currently only works for agmerra
clim <- c("agmerra","wfdei.gpcc","watch")
cfirst <- c(1980,1979,1958)
clast <- c(2010,2009,2001)
# using less RAy/Iizumi data for MA detrending
rayf <- 1961
rayl <- 2008
iizumif <- 1982
iizumil <- 2005 #2006
faof <- 1961
faol <- 2012
# specify which crop mask to use for aggregation for comparison with FPU.Iizumi, FPU.Ray, FAO
agg.fao <- "dynamic_ray_mask"
agg.fao <- "fixed_mirca_mask"
agg.iizumi <- "fixed_iizumi_mask"
agg.ray <- "dynamic_ray_mask"
aggs <- c("fixed_mirca_mask","dynamic_ray_mask","fixed_iizumi_mask","fixed_spam_mask")



# loop through some other models 
ggcms <- c("pdssat","epic-boku","epic-iiasa","gepic",
           "papsim","pegasus","lpj-guess","lpjml",
           "cgms-wofost","clm-crop","epic-tamu","orchidee-crop",
           "pepic","prysbi2")
ensembles <- c("rmse","tscorr")
# loop through crops
crops <- c("mai","whe","ric","soy")
cropsi <- c("maize_major","wheat","rice_major","soybean")
cropsl <- c("Maize","Wheat","Rice","Soy")
cropsl2 <- c("Maize","Wheat","Rice","Soybean")

pchs <- c(0:6,8,11,13,15:18)

#### topproducers. 

topproducer.soy <- c("global","USA","Brazil","Argentina","China","India","Paraguay","Canada","Uruguay","Ukraine","Bolivia")
gadm0.soy <- c(0,240,32,11,48,105,175,41,242,237,27)
topproducer.ric <- c("global","China","India","Indonesia","Bangladesh","Viet Nam","Thailand","Myanmar","Philippines","Brazil","Japan")
gadm0.ric <- c(0,48,105,106,19,247,226,154,177,32,114)
topproducer.mai <- c("global","USA","China","Brazil","Argentina","Mexico","India","Ukraine","Indonesia","France","South_Africa")
gadm0.mai <- c(0,240,48,32,11,145,105,237,106,79,209)
topproducer.whe <- c("global","China","India","USA","Russia","France","Canada","Australia","Pakistan","Germany","Turkey")
gadm0.whe <- c(0,48,105,240,186,79,41,14,170,86,232)

nf <- nc_open(paste(path,processed.ts,"/masks/aggr/gadm0.mask.nc4",sep=""))
gadm.mask <- ncvar_get(nf)
nc_close(nf)

#### define settings ####

do.gadm0 <- T
do.fpu <- F
do.raw <- F
do.norm <- T
do.png <- T
ignore.y <- 1
pcor <- F
do.rescale <- F
do.pixel <- F
do.means <- T

cr0 <- "none"#"mean-scale"#"variance-scale"#"none"
dt0 <- "ma" #"none"#"quad"
mp0 <- "true" #"false" #"true"#
scen0 <- "default"
nm0 <- 1 #number of top model for the ensemble
wt0 <- "unweighted" #wheighting for ensemble (1:unweighted, 2:weighted)

if(do.means) mp0 <- "true" # if only taking means, make sure to preserve them in the first place...

sig.threshold <- 0.1
shift.threshold <- 0.2

FRESHMATTER <- 100 / c(88, 88, 87, 91) 


prefix <- paste("cr_",cr0,".dt_",dt0,".mp_",mp0,".weighted_by_FAO_prod.shifted_ts",sep="")
if(do.means) prefix <- paste0(prefix,"_national_means_only")


#### do.gadm0 ####
#### national level comparison
if(do.gadm0){
  # use biascorr/gadm0/faostat
  topps <- c("global","top10")
  for(cl in 1){
    fao.season <- if(agg.fao=="dynamic_ray_mask") c(max(rayf,cfirst[cl],faof):min(rayl,clast[cl],faol)) - faof+1 else c(max(cfirst[cl],faof):min(clast[cl],faol)) - faof+1
    clim.season.f <- if(agg.fao=="dynamic_ray_mask") c(max(rayf,faof,cfirst[cl]) :min(rayl,faol,clast[cl])) - cfirst[cl]+1 else  c(max(faof,cfirst[cl]) :min(faol,clast[cl])) - cfirst[cl]+1  
    # remove first and last years
    fao.season <- c(fao.season[1+ignore.y]:fao.season[length(fao.season)-ignore.y])
    clim.season.f <- c((clim.season.f[1+ignore.y]):(clim.season.f[length(clim.season.f)-ignore.y]))    
      
    for(cc in crops){
      fname <- paste0(path.ref,"/faostat/faostat.1961-2012.gadm0.nc4")
      # get valid gadm0 IDs
      nf <- nc_open(fname)
      gadm0.id <- ncvar_get(nf,varid="gadm0")
      nc_close(nf)
      fname <- if(agg.ray=="dynamic_ray_mask")  paste0(path.ref,"/ray/ray.1961-2008.gadm0.ray.nc4") else  paste0(path.ref,"/ray/ray.1961-2008.gadm0.fixed.nc4")
      # get valid gadm0 IDs
      nf <- nc_open(fname)
      gadm0.id.r <- ncvar_get(nf,varid="gadm0")
      nc_close(nf)
      fname <- if(agg.ray=="dynamic_ray_mask")  paste0(path.ref,"/iizumi/iizumi.1982-2006.gadm0.iizumi.nc4") else  paste0(path.ref,"/iizumi/iizumi.1982-2006.gadm0.fixed.nc4")
      # get valid gadm0 IDs
      nf <- nc_open(fname)
      gadm0.id.i <- ncvar_get(nf,varid="gadm0")
      nc_close(nf)
      # f#*k. different country selection in sim files
      fname<-paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg.fao,"/","pdssat","_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                  clast[cl],".biascorr.nc4",sep="")
      nf <- nc_open(fname)
      gadm0.id2 <- ncvar_get(nf,varid="gadm0")
      nc_close(nf)
      sset <- which(gadm0.index %in% gadm0.id & gadm0.index %in% gadm0.id2)

      fname <-  paste0(path.ref,"/faostat/faostat.1961-2012.gadm0.nc4")
      r.fao <- get_nc4_ref_slice(fname,cc)[,which(gadm0.id %in% gadm0.id2)]
      # FAO production to use as weight
      r.fao2 <- get_nc4_ref_slice(fname,cc,mpd="true")[,which(gadm0.id %in% gadm0.id2)]
      a.fao <- get_nc4_ref_slice(fname,cc,getarea=T)[,which(gadm0.id %in% gadm0.id2)]
      p.fao <- r.fao2 * a.fao
      rm(a.fao,r.fao2)
      
      fname <-  paste0(path.ref,"/ray/ray.1961-2008.gadm0.ray.nc4")
      r.ray <- get_nc4_ref_slice(fname,cc)[,which(gadm0.id.r %in% gadm0.id2)]
      fname <-  paste0(path.ref,"/iizumi/iizumi.1982-2006.gadm0.iizumi.nc4")
      r.iizumi <- get_nc4_ref_slice(fname,cc)[,which(gadm0.id.i %in% gadm0.id2)]
      # read area per gadm0
      nf <- nc_open(paste(path,processed.ts,"/masks/weight/aggs/",tolower(cropsl[which(crops==cc)]),".gadm0.nc4",sep=""))
      gadm0.id.area <- ncvar_get(nf,varid="gadm0_index")
      area_per_gadm0_rf <- ncvar_get(nf,varid="rainfed_gadm0")[which(gadm0.id.area %in% gadm0.id2)]
      area_per_gadm0_ir <- ncvar_get(nf,varid="irrigated_gadm0")[which(gadm0.id.area %in% gadm0.id2)]
      nc_close(nf)
      # use FAO production as weight
      ww.f <- p.fao[clim.season.f,]
      
      gadm0.names <- name[sset]
      gadm0.id <- 1:length(gadm0.names)
      fyear.plot <- max(cfirst[cl],faof)+ignore.y
      lyear.plot <- min(clast[cl],faol)-ignore.y

      while(length(dev.list())>0) dev.off() # make sure no figures are open to allow for 2 simultaneous figures
      for(topp in topps){
        if(do.png){
          png(paste(path,"taylor/",cc,".taylor.",topp,".gadm0.best_mask.",clim[cl],".",prefix,".",agg.fao,".",Sys.Date(),".png",sep=""),height=5*300,width=5*300,res=300,pointsize=9,type="cairo")
        } else{
          postscript(paste(path,"taylor/",cc,".taylor.",topp,".gadm0.best_mask.",clim[cl],".",prefix,".",agg.fao,".",Sys.Date(),".eps",sep=""),height=5*100,width=5*100,pointsize=22,family="NimbusSan")
        }
        par(xpd=NA)
        if(topp=="global"){
          gadm0s <- gadm0.id
        }else{
          gadm0s <- which(gadm0.names %in% if(cc=="whe") topproducer.whe[2:11] else if(cc=="mai") topproducer.mai[2:11] else if(cc=="ric") topproducer.ric[2:11] else topproducer.soy[2:11])
        }
        # plot fao in taylor diagram
        pos <- taylor.diagram2(as.vector(r.fao[c(fyear.plot:lyear.plot)-faof+1,gadm0s]),
                               as.vector(r.fao[c(fyear.plot:lyear.plot)-faof+1,gadm0s]),
                               as.vector(ww.f[,gadm0s]),ref.sd=T,
                               show.gamma=F,normalize=do.norm,col="grey60",pos.cor=F,#pos.cor=if(topp=="global") T else F,
                               main="",#paste(topp,clim[cl],cropsl[which(crops==cc)],"\n",prefix,"\n",agg.fao),
                               type="n",sig.threshold=sig.threshold)
        
        # loop through models and add points to taylor diagram
        for(gg in ggcms){
          cat(gg,": ",sep="")
          sdf <- shf <- sff <- 0
          for(agg in aggs){
            fao.season <- c(max(rayf,cfirst[cl],faof):min(rayl,clast[cl],faol)) - faof+1
            clim.season.f <- c(max(rayf,faof,cfirst[cl]) :min(rayl,faol,clast[cl])) - cfirst[cl]+1  
            # remove first and last years
            fao.season <- c(fao.season[1+ignore.y]:(fao.season[length(fao.season)-ignore.y]))
            clim.season.f <- c((clim.season.f[1+ignore.y]):(clim.season.f[length(clim.season.f)-ignore.y]))    
            # use FAO production as weight
            ww.f <- p.fao[clim.season.f,]
            
            fn <- paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.default.f <- if(is.finite(file.info(fn)$size)) get_nc4_data_slice(fname=fn,scend="default") * FRESHMATTER[which(crops==cc)] else 0
            fn <- paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.fullharm.f <- if(is.finite(file.info(fn)$size)) get_nc4_data_slice(fname=fn,scend="fullharm") * FRESHMATTER[which(crops==cc)] else 0
            fn <- paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.harmnon.f <- if(is.finite(file.info(fn)$size)) get_nc4_data_slice(fname=fn,scend="harmnon") * FRESHMATTER[which(crops==cc)] else 0
            # shift time series by +/-1 year to test if the time series should be shifted per country
            # need to shift FAO data as that allows for keeping all non-NA data points
            r.fao.d <- r.fao.f <- r.fao.h <- r.fao # resetting to original
            for(gad in gadm0s){
              if(length(data.default.f)>1) r.fao.d[fao.season,gad] <- adjust.ts(data.default.f[clim.season.f,gad],r.fao[c((fao.season[1]-1):(fao.season[length(fao.season)]+1)),gad])
              if(length(data.fullharm.f)>1) r.fao.f[fao.season,gad] <- adjust.ts(data.fullharm.f[clim.season.f,gad],r.fao[c((fao.season[1]-1):(fao.season[length(fao.season)]+1)),gad])
              if(length(data.harmnon.f)>1) r.fao.h[fao.season,gad] <- adjust.ts(data.harmnon.f[clim.season.f,gad],r.fao[c((fao.season[1]-1):(fao.season[length(fao.season)]+1)),gad])
            }
            if(agg!=aggs[1]){
              if(length(data.default.f)>1){
                buf <- wtd.cor.mod(r.fao.d[fao.season,gadm0s],
                               data.default.f[clim.season.f,gadm0s],
                               ww.f[,gadm0s])
                if(buf[1]>cd[1]){
                  cd <- buf
                  dd.f <- data.default.f[clim.season.f,gadm0s]   
                  sdf <- sqrt(wtd.var(as.vector(data.default.f[clim.season.f,gadm0s]),as.vector(ww.f[,gadm0s]),na.rm=T))
                }
              }
              if(length(data.fullharm.f)>1){
                buf <- wtd.cor.mod(r.fao.f[fao.season,gadm0s],
                                data.fullharm.f[clim.season.f,gadm0s],
                                ww.f[,gadm0s])
                if(buf[1]>cf[1]){
                  cf <- buf
                  df.f <- data.fullharm.f[clim.season.f,gadm0s] 
                  sff <- sqrt(wtd.var(as.vector(data.fullharm.f[clim.season.f,gadm0s]),as.vector(ww.f[,gadm0s]),na.rm=T))
                }
              }
              if(length(data.harmnon.f)>1){
                buf <- wtd.cor.mod(r.fao.h[fao.season,gadm0s],
                               data.harmnon.f[clim.season.f,gadm0s],
                               ww.f[,gadm0s])
                if(buf[1]>ch[1]){
                  ch <- buf
                  dh.f <- data.harmnon.f[clim.season.f,gadm0s]
                  shf <- sqrt(wtd.var(as.vector(data.harmnon.f[clim.season.f,gadm0s]),as.vector(ww.f[,gadm0s]),na.rm=T))
                }
              }
            } else {
              if(length(data.default.f)>1){
                dd.f <- data.default.f[clim.season.f,gadm0s]
                cd <- wtd.cor.mod(r.fao.d[fao.season,gadm0s],dd.f,ww.f[,gadm0s])
                sdf <- sqrt(wtd.var(as.vector(dd.f),as.vector(ww.f[,gadm0s]),na.rm=T))
              } else{
                dd.f <- NA
                cd <- -1
              }
              if(length(data.fullharm.f)>1){
                df.f <- data.fullharm.f[clim.season.f,gadm0s]
                cf <- wtd.cor.mod(r.fao.f[fao.season,gadm0s],df.f,ww.f[,gadm0s])
                sff <- sqrt(wtd.var(as.vector(df.f),as.vector(ww.f[,gadm0s]),na.rm=T))
              } else{
                df.f <- NA
                cf <- -1
              }
              if(length(data.harmnon.f)>1){
                dh.f <- data.harmnon.f[clim.season.f,gadm0s]
                ch <- wtd.cor.mod(r.fao.h[fao.season,gadm0s],dh.f,ww.f[,gadm0s])
                shf <- sqrt(wtd.var(as.vector(dh.f),as.vector(ww.f[,gadm0s]),na.rm=T))
              } else{
                dh.f <- NA
                ch <- -1
              }
            }
          }# end aggs
          #normalize sds
          if(do.norm){
            sd.r <- sqrt(wtd.var(as.vector(r.fao[fao.season,gadm0s]),as.vector(ww.f[,gadm0s]),na.rm=T))
            sdf <- sdf/sd.r
            sff <- sff/sd.r
            shf <- shf/sd.r
            sd.r <- 1
          }
          if(length(which(!is.na(dd.f)))>1) {
            if(cd[4]>sig.threshold) {
              col2 <- col2rgb("#4575b4")
              points(sdf * cd[1], sdf * sin(acos(cd[1])), pch = pchs[which(ggcms==gg)], 
                     col = rgb(col2[1,],col2[2,],col2[3,],30,maxColorValue=255), 
                     cex = 1)    
              cat("bad p value!",cf[4],"\n")
              
            } else{
              points(sdf * cd[1], sdf * sin(acos(cd[1])), pch = pchs[which(ggcms==gg)], col = "#4575b4", 
                     cex = 1) 
              cat("plotting",sdf*cd[1],sdf,cd[1],sdf*sin(acos(cd[1])),sin(acos(cd[1])),cd[4],"\n")
            }
          } else cat("skipping",gg,cc,topp,clim[cl],"default with fao\n")
          if(length(which(!is.na(df.f)))>1) {
            if(cf[4]>sig.threshold) {
              col2 <- col2rgb("#d73027")
              points(sff * cf[1], sff * sin(acos(cf[1])), pch = pchs[which(ggcms==gg)], 
                     col = rgb(col2[1,],col2[2,],col2[3,],30,maxColorValue=255), 
                     cex = 1)    
              cat("bad p value!",cf[4],"\n")
              
            } else{
              points(sff * cf[1], sff * sin(acos(cf[1])), pch = pchs[which(ggcms==gg)], col = "#d73027", 
                     cex = 1)                 
              cat("plotting",sff*cf[1],sff,cf[1],sff*sin(acos(cf[1])),sin(acos(cf[1])),cf[4],"\n")
            }
            
          } else cat("skipping",gg,cc,topp,clim[cl],"fullharm with fao\n")
          if(length(which(!is.na(dh.f)))>1) {
            if(ch[4]>sig.threshold) {
              col2 <- col2rgb("#fee090")
              points(shf * cd[1], shf * sin(acos(ch[1])), pch = pchs[which(ggcms==gg)], 
                     col = rgb(col2[1,],col2[2,],col2[3,],30,maxColorValue=255), 
                     cex = 1)    
              cat("bad p value!",ch[4],"\n")
              
            } else{
              points(shf * ch[1], shf * sin(acos(ch[1])), pch = pchs[which(ggcms==gg)], col = "#fee090", 
                     cex = 1)                 
              cat("plotting",shf*ch[1],shf,ch[1],shf*sin(acos(ch[1])),sin(acos(ch[1])),ch[4],"\n")
            }
            
          } else cat("skipping",gg,cc,topp,clim[cl],"harmnon with fao\n")          
        } # for ggcms
        

        legend(-2.9,3.9 ,legend=c(ggcms),ncol=3,
               pch=pchs,col=c(rep(1,length(ggcms))),bty="n",cex=0.8,xjust=0)
        legend(2.9,3.9 ,legend=c("default vs. FAOstat","fullharm vs. FAOstat","harm-suffN vs. FAOstat"),
               pch=c(16,16,16),col=c("#4575b4","#d73027","#fee090"),
               bty="n",cex=0.8,xjust=1)
        
        dev.off()
      }
    }
  }
  
}
