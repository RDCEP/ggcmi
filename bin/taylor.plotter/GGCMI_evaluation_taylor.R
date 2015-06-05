# issues
# why is there no mask/weight/aggs/*.gadm0.nc4 ?
# why is model ensemble for soy from 1979 not from 1980? also inconsistencies for other climate data products??
# 


require(ncdf4)
#require(plotrix) # not needed as I'm using a tweaked copy of the taylor.diagram() function
rm(list=ls(all=TRUE))

# script to plot Taylor diagrams for GGCMI phase 1. 
# written by Christoph Müller, PIK
# cmueller@pik-potsdam.de
# I'm using a slightly modified version of taylor.diagram() of the plotrix package
# to avoid comparison with NA values and to allow for weighted correlation analysis.
# For that I've copied the source code of taylor.diagram to GGCMI_taylor_function.R and modified it accordingly.

# intermediate version
# TODO
# * adjust functions for reading national level data for new nc file structure
# * adjust functions for reading detrended reference data for new nc file structure 
# * check other issues
# * do pixel-based assessment for individual countries

#path <- "/iplex/01/2011/isimip-lpj/yields/GGCMI2_outputs/WFDEI_GPCC_default/ncdf/"
path <- "D:/data/GGCMI/"
path.ref <- "D:/Dropbox/GGCMI/"
path.ref.detrend <- "D:/data/GGCMI/ray-dt/"
processed.ts <- "processed.150505"

# read FPU names
fpu.names <- read.csv("D:/data/GGCMI/fpu.meta.csv",header=F)[,2]
gadm0.names <- array("missing",253)
gadm0.index <- read.csv("D:/data/GGCMI/gadm0.meta.csv",header=F)[,1]
name <- read.csv("D:/data/GGCMI/gadm0.meta.csv",header=F)[,2]


source("D:/github/ggcmi/bin/taylor.plotter/GGCMI_taylor_function.R")

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

get_nc4_ref_slice <- function(fname,crop,dtd=dt0,mpd=mp0){
  nf <- nc_open(fname)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  # order is mp, dt, time, fpu
  if(crop=="mai") { data <- ncvar_get(nf,varid="yield_mai")[mp,dt,,]
  } else if(crop=="whe"){ data <- ncvar_get(nf,varid="yield_whe")[mp,dt,,]
  } else if(crop=="ric"){ data <- ncvar_get(nf,varid="yield_ric")[mp,dt,,]
  } else if(crop=="soy"){ data <- ncvar_get(nf,varid="yield_soy")[mp,dt,,]}
  nc_close(nf)
  data
}

# different function needed for detreneded ref data, as there are more dimensions. TODO
# get_nc4_refdt_slice <- function(fname,crd,dtd=dt0,mpd=mp0,nmd=1,){
#   #forgot what this was about. TODO
# }

# different function needed for national biascorr data, which now have additional dimensions. TODO

# loop through weather products TODO
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
agg.ray <- "fixed_mirca_mask"
agg.iizumi <- "fixed_mirca_mask"


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

#### topproducers. FPU topproducer regions are just randomly selected ones from the topproducer countries

topproducer.soy <- c("global","USA","Brazil","Argentina","China","India","Paraguay","Canada","Uruguay","Ukraine","Bolivia")
gadm0.soy <- c(0,240,32,11,48,105,175,41,242,237,27)
fpu.soy <- c(0,219,206,190,193,286,178,182,243,69,234)
topproducer.ric <- c("global","China","India","Indonesia","Bangladesh","Viet Nam","Thailand","Myanmar","Philippines","Brazil","Japan")
gadm0.ric <- c(0,48,105,106,19,247,226,154,177,32,114)
fpu.ric <- c(0,193,286,290,31,194,131,226,181,206,294)
topproducer.mai <- c("global","USA","China","Brazil","Argentina","Mexico","India","Ukraine","Indonesia","France","Canada")
gadm0.mai <- c(0,240,48,32,11,145,105,237,106,79,41)
fpu.mai <- c(0,134,281,206,190,191,286,69,290,207,182)
topproducer.whe <- c("global","China","India","USA","Russia","France","Canada","Australia","Pakistan","Germany","Turkey")
gadm0.whe <- c(0,48,105,240,186,79,41,14,170,86,232)
fpu.whe <- c(0,281,103,133,24,207,50,135,293,187,65)

nf <- nc_open(paste(path,processed.ts,"/masks/aggr/gadm0.mask.nc4",sep=""))
gadm.mask <- ncvar_get(nf)
nc_close(nf)

#### define settings ####

do.gadm0 <- T
do.fpu <- T
do.raw <- T
do.norm <- T
do.png <- T
ignore.y <- 1
pcor <- F
do.rescale <- F
do.pixel <- T

cr0 <- "none"#"mean-scale"#"variance-scale"#"none"
dt0 <- "quad" #"none"#"quad"
mp0 <- "false"
scen0 <- "default"
nm0 <- 1 #number of top model for the ensemble
wt0 <- "unweighted" #wheighting for ensemble (1:unweighted, 2:weighted)


prefix <- paste("cr_",cr0,".dt_",dt0,".mp_",mp0,".weighted_by_ir+rf_mask",sep="")
#prefix <- "cr_mean-scale.dt_quad.mp_true"

#### do.fpu ####
# use biascorr/fpu/ray
if(do.fpu){
  # use biascorr/fpu/ray
  topps <- c("global","top10")
  for(cl in 1){
    ray.season <- c(max(cfirst[cl],rayf):min(clast[cl],rayl)) - rayf+1
    iizumi.season <- c(max(cfirst[cl],iizumif):min(clast[cl],iizumil)) - iizumif+1
    clim.season.r <- c(max(rayf,cfirst[cl]) :min(rayl,clast[cl])) - cfirst[cl]+1  
    clim.season.i <- c(max(iizumif,cfirst[cl]):min(iizumil,clast[cl])) - cfirst[cl]+1 
    # remove first and last years
    ray.season <- c(ray.season[1+ignore.y]:ray.season[length(ray.season)-ignore.y])
    iizumi.season <- c((iizumi.season[1+ignore.y]):(iizumi.season[length(iizumi.season)-ignore.y]))
    clim.season.r <- c((clim.season.r[1+ignore.y]):(clim.season.r[length(clim.season.r)-ignore.y]))
    clim.season.i <- c((clim.season.i[1+ignore.y]):(clim.season.i[length(clim.season.i)-ignore.y]))
    
    for(cc in crops){
    #for(cc in "mai"){
      fname <- if(agg.ray=="dynamic_ray_mask") "D:/data/GGCMI/reference/ray/ray.1961-2008.fpu.dynamic.nc4" else "D:/data/GGCMI/reference/ray/ray.1961-2008.fpu.fixed.nc4"
      r.ray <- get_nc4_ref_slice(fname,cc)
      fname <- if(agg.iizumi =="dynamic_ray_mask") "D:/data/GGCMI/reference/iizumi/iizumi.1982-2006.fpu.dynamic.nc4" else "D:/data/GGCMI/reference/iizumi/iizumi.1982-2006.fpu.fixed.nc4"
      r.iizumi <- get_nc4_ref_slice(fname,cc)
      fyear.plot <- max(cfirst[cl],rayf,iizumif)+ignore.y
      lyear.plot <- min(clast[cl],rayl,iizumil)-ignore.y
      # for testing
      if(F){
        png(paste("D:/data/GGCMI/reference/comparison_Italy_",fyear.plot,"-",lyear.plot,"_",prefix,".png",sep=""))
        #randomcountries <- which(!is.na(r.ray[1,]))[9:10]
        # check Italy, which is FPU==full country
        randomcountries <- 105
        plot(c(fyear.plot:lyear.plot),r.ray[c(fyear.plot:lyear.plot)-rayf+1,randomcountries[1]],
             ylim=range(r.ray[c(fyear.plot:lyear.plot)-rayf+1,randomcountries],
                        r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,randomcountries],
                        r.fao[c(fyear.plot:lyear.plot)-rayf+1,90],na.rm=T),
             type="n",ylab="yield [t/ha]",main=paste("Italy",fyear.plot,"-",lyear.plot,prefix,sep=" "))
        for(rci in randomcountries){
          lines(c(fyear.plot:lyear.plot),r.ray[c(fyear.plot:lyear.plot)-rayf+1,rci],col=rci)
          lines(c(fyear.plot:lyear.plot),r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,rci],col=rci,lty=2)
          lines(c(fyear.plot:lyear.plot),r.fao[c(fyear.plot:lyear.plot)-faof+1,90],col=rci,lty=3)
        }
        legend("topleft",legend=c("Ray","Iizumi","FAO"),lty=c(1,2,3))
        dev.off()
        png(paste("D:/data/GGCMI/reference/comparison_mean_per_FPU_",fyear.plot,"-",lyear.plot,"_",prefix,".png",sep=""),
            width=5*1200,height=3*1200,res=1200,pointsize=4)
        #randomcountries <- which(!is.na(r.ray[1,]))[9:10]
        # check Italy, which is FPU==full country
        diff.ref <- colMeans(r.ray[c(fyear.plot:lyear.plot)-rayf+1,])-colMeans(r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,])
        prange <- range(r.ray[c(fyear.plot:lyear.plot)-rayf+1,],diff.ref,
                        r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,],na.rm=T)
        plot(1:dim(r.ray)[2],colMeans(r.ray[c(fyear.plot:lyear.plot)-rayf+1,]),
             ylim=prange,
             type="p",ylab="yield [t/ha]",main=paste("FPU means",fyear.plot,"-",lyear.plot,prefix,sep=" "),
             col=2,pch=3,axes=F)
        axis(2)
        #axis(1,crt=90,labels=fpu.names,at=c(1:309),cex=.4)
        text(x=seq(1,by=1,length.out=dim(r.ray)[2]),y=min(diff.ref,na.rm=T),labels=fpu.names,
             adj=c(1,.5),cex=0.3,srt=90,xpd=NA)
        
        points(1:dim(r.ray)[2],colMeans(r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,]),
               col=3,pch=3)
        points(1:dim(r.ray)[2],diff.ref,
               col=4,pch=4)
        for(i in 1:dim(r.ray)[2]) lines(c(i,i),prange,lwd=0.1)
        for(i in c(as.integer(prange[1]):as.integer(prange[2]))) lines(c(1,dim(r.ray)[2]),c(i,i),lwd=0.1,lty=2)
        legend("topleft",legend=c("Ray","Iizumi","difference"),col=c(2,3,4),pch=c(3,3,4))
        dev.off()
        
      }
      # read area per FPU
      nf <- nc_open(paste("D:/data/GGCMI/",processed.ts,"/masks/weight/aggs/",tolower(cropsl[which(crops==cc)]),".fpu.nc4",sep=""))
      area_per_fpu_rf <- ncvar_get(nf,varid="rainfed_fpu")
      area_per_fpu_ir <- ncvar_get(nf,varid="irrigated_fpu")
      nc_close(nf)
      ww.r <- matrix(rep(cbind(area_per_fpu_rf+area_per_fpu_ir),length(clim.season.r)),nrow=length(clim.season.r),byrow=T)
      ww.i <- matrix(rep(cbind(area_per_fpu_rf+area_per_fpu_ir),length(clim.season.i)),nrow=length(clim.season.i),byrow=T)
      while(length(dev.list())>0) dev.off() # make sure not figures are open to allow for 2 simultaneous figures
      for(topp in topps){
        if(do.png){
          png(paste(path,"taylor/",cc,".taylor.",topp,".fpu.",clim[cl],".",prefix,".","r_vs_",agg.ray,".i_vs_",agg.iizumi,".",Sys.Date(),".png",sep=""),height=5*300,width=5*300,res=300,pointsize=9)
        } else{
          postscript(paste(path,"taylor/",cc,".taylor.",topp,".fpu.",clim[cl],".",prefix,".","r_vs_",agg.ray,".i_vs_",agg.iizumi,".",Sys.Date(),".eps",sep=""),height=5*100,width=5*100,pointsize=22,family="NimbusSan")
        }
        par(xpd=NA)
        if(topp=="global"){
          fpus <- c(1:309)
        }else{
          # this is for country level
          #fpus <- if(cc=="whe") gadm.whe[2:11] else if(cc=="mai") gadm.mai[2:11] else if(cc=="ric") gadm.ric[2:11] else if(cc=="soy") gadm.soy[2:11]
          # this is a proxy for FPU level, just using the 10 highest productivity ones
          #fpus <- order(colMeans(r.ray),decreasing=T)[1:10]
          fpus <- if(cc=="whe") fpu.whe[2:11] else if(cc=="mai") fpu.mai[2:11] else if(cc=="ric") fpu.ric[2:11] else if(cc=="soy") fpu.soy[2:11]
        }
        # plot ray vs. iizumi in taylor diagram
        pos <- taylor.diagram2(as.vector(r.ray[c(fyear.plot:lyear.plot)-rayf+1,fpus]),
                               as.vector(r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,fpus]),
                               as.vector(ww.r[1:(lyear.plot-fyear.plot+1),fpus]),
                               ref.sd=T,
                               show.gamma=F,normalize=do.norm,col="grey60",pos.cor=pcor,#pos.cor=if(topp=="global") T else F,
                               main=if(do.png) paste(topp,clim[cl],cropsl[which(crops==cc)],"\n",prefix) else "")
        
        # loop through models and add points to taylor diagram
        for(gg in ggcms){
          cat(gg,": ",sep="")
          #if(!((cc=="ric" & gg=="papsim")|(cc=="ric" & gg=="pegasus"))){ #no rice for papsim and pegasus
            fn <- paste(path,processed.ts,"/biascorr/fpu/ray/",agg.ray,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.default.r <- if(is.finite(file.info(fn)$size>0))  get_nc4_data_slice(fname=fn,scend="default") else 0
            fn <- paste(path,processed.ts,"/biascorr/fpu/ray/",agg.ray,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.fullharm.r <- if(is.finite(file.info(fn)$size>0))  get_nc4_data_slice(fname=fn,scend="fullharm") else 0
            fn <- paste(path,processed.ts,"/biascorr/fpu/ray/",agg.ray,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.harmnon.r <- if(is.finite(file.info(fn)$size>0))  get_nc4_data_slice(fname=fn,scend="harmnon") else 0
            fn <- paste(path,processed.ts,"/biascorr/fpu/iizumi/",agg.iizumi,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.default.i <- if(is.finite(file.info(fn)$size>0))  get_nc4_data_slice(fname=fn,scend="default") else 0
            fn <- paste(path,processed.ts,"/biascorr/fpu/iizumi/",agg.iizumi,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.fullharm.i <- if(is.finite(file.info(fn)$size>0))  get_nc4_data_slice(fname=fn,scend="fullharm") else 0
            fn <- paste(path,processed.ts,"/biascorr/fpu/iizumi/",agg.iizumi,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                        if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
            data.harmnon.i <- if(is.finite(file.info(fn)$size>0))  get_nc4_data_slice(fname=fn,scend="harmnon") else 0
            dd.r <- if(length(data.default.r)>1) data.default.r[clim.season.r,fpus] else NA
            df.r <- if(length(data.fullharm.r)>1) data.fullharm.r[clim.season.r,fpus] else NA
            dh.r <- if(length(data.harmnon.r)>1) data.harmnon.r[clim.season.r,fpus] else NA
            dd.i <- if(length(data.default.i)>1) data.default.i[clim.season.i,fpus] else NA
            df.i <- if(length(data.fullharm.i)>1) data.fullharm.i[clim.season.i,fpus] else NA
            dh.i <- if(length(data.harmnon.i)>1) data.harmnon.i[clim.season.i,fpus] else NA
            
            # plot time series of models for top10
            if(topp==topps[2]){
              ref.vec <- as.vector(r.iizumi[iizumi.season,fpus])
              sim.vec.dd <- as.vector(dd.i)
              sim.vec.df <- as.vector(df.i)
              sim.vec.dh <- as.vector(dh.i)
              ww.vec <- as.vector(ww.i[,fpus])
              if(length(dim(dd.i))>0){
                png(paste(path,"timeseries/",gg,"/",gg,".",cc,".timeseries.",topp,".jointly.iizumi.",clim[cl],".",prefix,".","r_vs_",agg.ray,"i_vs_",agg.iizumi,".",Sys.Date(),".png",sep=""),height=3*300,width=5*300,res=300,pointsize=9)
                par(cex=0.6)
                prange <- range(ref.vec,sim.vec.dd,sim.vec.dh,sim.vec.df,na.rm=T)
                plot(ref.vec,ylim=prange,type="o",cex=0.6,xlab="",ylab="")
                text(x=seq(1,by=length(iizumi.season),length.out=length(fpus)),y=prange[2]+(prange[2]-prange[1])*.1,labels=fpu.names[fpus],
                     adj=0,cex=0.6,srt=20,xpd=NA)
                text(x=seq(1,by=length(iizumi.season),length.out=length(fpus)),y=prange[2]+(prange[2]-prange[1])*.075,labels=paste("weight: ",round(ww.i[1,fpus])),
                     adj=0,cex=0.6,srt=20,xpd=NA)
                text(x=seq(1,by=length(iizumi.season),
                           length.out=length(fpus)),
                     y=prange[2]+(prange[2]-prange[1])*.05,
                     labels=paste("cor: ",
                                  c(round(cor(r.iizumi[iizumi.season,fpus[1]],dd.i[,1]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[2]],dd.i[,2]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[3]],dd.i[,3]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[4]],dd.i[,4]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[5]],dd.i[,5]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[6]],dd.i[,6]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[7]],dd.i[,7]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[8]],dd.i[,8]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[9]],dd.i[,9]),3),
                                    round(cor(r.iizumi[iizumi.season,fpus[10]],dd.i[,10]),3))),
                     adj=0,cex=0.6,srt=20,xpd=NA)
                lines(sim.vec.dd,type="o",cex=0.6,lwd=0.5,col="green")
                lines(sim.vec.df,type="o",cex=0.6,lwd=0.5,col="orange")
                lines(sim.vec.dh,type="o",cex=0.6,lwd=0.5,col="yellow2")
                legend("topleft",legend=c("Iizumi",paste("default r: ",
                                                         round(cor(ref.vec,sim.vec.dd,"pairwise"),3)," weighted: ",round(wtd.cor(ref.vec,sim.vec.dd,ww.vec)[1],3)),
                                          paste("fullharm r: ",
                                                if(length(which(!is.na(sim.vec.df)))>0) round(cor(ref.vec,sim.vec.df,"pairwise"),3),
                                                " weighted: ",if(length(which(!is.na(sim.vec.df)))>0) round(wtd.cor(ref.vec,sim.vec.df,ww.vec)[1],3)),
                                          paste("harmnon r: ",if(length(which(!is.na(sim.vec.dh)))>0) round(cor(ref.vec,sim.vec.dh,"pairwise"),3),
                                                " weighted: ",if(length(which(!is.na(sim.vec.dh)))>0) round(wtd.cor(ref.vec,sim.vec.dh,ww.vec)[1],3))),
                       col=c(1,"green","orange","yellow2"),lty=1,bty="n")
                dev.off()                
              }
              if(length(dim(dd.r))>0){
                png(paste(path,"timeseries/",gg,"/",gg,".",cc,".timeseries.",topp,".jointly.ray.",clim[cl],".",prefix,".","r_vs_",agg.ray,".i_vs_",agg.iizumi,".",Sys.Date(),".png",sep=""),height=3*300,width=5*300,res=300,pointsize=9)
                par(cex=0.6)
                ref.vec <- as.vector(r.ray[ray.season,fpus])
                sim.vec.dd <- as.vector(dd.r)
                sim.vec.df <- as.vector(df.r)
                sim.vec.dh <- as.vector(dh.r)
                prange <- range(ref.vec,sim.vec.dd,sim.vec.dh,sim.vec.df,na.rm=T)
                plot(ref.vec,ylim=prange,type="o",cex=0.6,xlab="",ylab="")
                text(x=seq(1,by=length(ray.season),length.out=length(fpus)),y=prange[2]+(prange[2]-prange[1])*.1,labels=fpu.names[fpus],
                     adj=0,cex=0.6,srt=20,xpd=NA)
                text(x=seq(1,by=length(ray.season),length.out=length(fpus)),y=prange[2]+(prange[2]-prange[1])*.075,labels=paste("weight: ",round(ww.i[1,fpus])),
                     adj=0,cex=0.6,srt=20,xpd=NA)
                text(x=seq(1,by=length(ray.season),
                           length.out=length(fpus)),
                     y=prange[2]+(prange[2]-prange[1])*.05,
                     labels=paste("cor: ",
                                  c(round(cor(r.ray[ray.season,fpus[1]],dd.r[,1]),3),
                                    round(cor(r.ray[ray.season,fpus[2]],dd.r[,2]),3),
                                    round(cor(r.ray[ray.season,fpus[3]],dd.r[,3]),3),
                                    round(cor(r.ray[ray.season,fpus[4]],dd.r[,4]),3),
                                    round(cor(r.ray[ray.season,fpus[5]],dd.r[,5]),3),
                                    round(cor(r.ray[ray.season,fpus[6]],dd.r[,6]),3),
                                    round(cor(r.ray[ray.season,fpus[7]],dd.r[,7]),3),
                                    round(cor(r.ray[ray.season,fpus[8]],dd.r[,8]),3),
                                    round(cor(r.ray[ray.season,fpus[9]],dd.r[,9]),3),
                                    round(cor(r.ray[ray.season,fpus[10]],dd.r[,10]),3))),
                     adj=0,cex=0.6,srt=20,xpd=NA)
                lines(sim.vec.dd,type="o",cex=0.6,lwd=0.5,col="green")
                lines(sim.vec.df,type="o",cex=0.6,lwd=0.5,col="orange")
                lines(sim.vec.dh,type="o",cex=0.6,lwd=0.5,col="yellow2")
                legend("topleft",legend=c("ray",paste("default r: ",
                                                      round(cor(ref.vec,sim.vec.dd,"pairwise"),3)," weighted: ",round(wtd.cor(ref.vec,sim.vec.dd,ww.vec)[1],3)),
                                          paste("fullharm r: ",
                                                if(length(which(!is.na(sim.vec.df)))>0) round(cor(ref.vec,sim.vec.df,"pairwise"),3),
                                                " weighted: ",if(length(which(!is.na(sim.vec.df)))>0) round(wtd.cor(ref.vec,sim.vec.df,ww.vec)[1],3)),
                                          paste("harmnon r: ",if(length(which(!is.na(sim.vec.dh)))>0) round(cor(ref.vec,sim.vec.dh,"pairwise"),3),
                                                " weighted: ",if(length(which(!is.na(sim.vec.dh)))>0) round(wtd.cor(ref.vec,sim.vec.dh,ww.vec)[1],3))),
                       col=c(1,"green","orange","yellow2"),lty=1,bty="n")
                dev.off()                
              }
              png(paste(path,"timeseries/",gg,"/",gg,".",cc,".timeseries.",topp,".fpu.",clim[cl],".",prefix,".","r_vs_",agg.ray,".i_vs_",agg.iizumi,".",Sys.Date(),".png",sep=""),height=10*300,width=5*300,res=300,pointsize=9)
              split.screen(c(5,2))
              for(i in 1:length(fpus)){
                screen(i)
                par(cex=0.7,mar=c(2,2,2.5,0))
                if(length(which(!is.na(r.ray[ray.season,fpus[i]])))>0 | length(which(!is.na(r.iizumi[iizumi.season,fpus[i]])))){
                  drange <- range(r.ray[ray.season,fpus[i]],r.iizumi[iizumi.season,fpus[i]],
                                  if(length(dim(dd.r))>0)dd.r[,i],if(length(dim(dd.i))>0)dd.i[,i],
                                  if(length(dim(df.r))>0)df.r[,i],if(length(dim(df.i))>0)df.i[,i],
                                  if(length(dim(dh.r))>0)dh.r[,i],if(length(dim(dh.i))>0)dh.i[,i],na.rm=T)
                  plot(c(rayf:rayl)[ray.season],r.ray[ray.season,fpus[i]],type="o",ylim=drange,pch=16,xlab="",ylab="",cex=0.6,
                       main=fpu.names[fpus[i]])  
                  if(length(dim(dd.r))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.r],dd.r[,i],col=rainbow(10)[2],type="o",pch=16,cex=0.6)
                  if(length(dim(df.r))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.r],df.r[,i],col=rainbow(10)[4],type="o",pch=16,cex=0.6)
                  if(length(dim(dh.r))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.r],dh.r[,i],col=rainbow(10)[6],type="o",pch=16,cex=0.6)
                  lines(c(iizumif:iizumil)[iizumi.season],r.iizumi[iizumi.season,fpus[i]],type="o",ylim=drange,lty=2,pch=1,cex=0.6)  
                  if(length(dim(dd.i))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.i],dd.i[,i],col=rainbow(10)[2],type="o",lty=2,pch=1,cex=0.6)
                  if(length(dim(df.i))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.i],df.i[,i],col=rainbow(10)[4],type="o",lty=2,pch=1,cex=0.6)
                  if(length(dim(dh.i))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.i],dh.i[,i],col=rainbow(10)[6],type="o",lty=2,pch=1,cex=0.6)
                  par(cex=0.6)
                  legend("topleft",bty="n",legend=c("ray","iizumi","default","fullharm","harmnon"),
                         col=c(1,1,rainbow(10)[c(2,4,6)]),lty=c(1,2,1,1,1),pch=c(16,1,1,1))
                  
                }
              }
              close.screen(all=T)
              
              dev.off()
            }
            
            if(length(which(!is.na(dd.r)))>1) {
              taylor.diagram2(as.vector(r.ray[ray.season,fpus]),as.vector(dd.r),as.vector(ww.r[,fpus]),add=T,col="blue",pch=which(ggcms==gg),normalize=do.norm)
            } else cat("skipping",gg,cc,topp,clim[cl],"default with Ray\n")
            if(length(which(!is.na(dd.i)))>1) {
              taylor.diagram2(as.vector(r.iizumi[iizumi.season,fpus]),as.vector(dd.i),as.vector(ww.i[,fpus]),add=T,col="green",pch=which(ggcms==gg),normalize=do.norm)
            } else cat("skipping",gg,cc,topp,clim[cl],"default with Iizumi\n")
            if(length(which(!is.na(df.r)))>1) {
              taylor.diagram2(as.vector(r.ray[ray.season,fpus]),as.vector(df.r),as.vector(ww.r[,fpus]),add=T,col="red",pch=which(ggcms==gg),normalize=do.norm)
            } else cat("skipping",gg,cc,topp,clim[cl],"fullharm with Ray\n")
            if(length(which(!is.na(df.i)))>1) {
              taylor.diagram2(as.vector(r.iizumi[iizumi.season,fpus]),as.vector(df.i),as.vector(ww.i[,fpus]),add=T,col="orange",pch=which(ggcms==gg),normalize=do.norm)
            } else cat("skipping",gg,cc,topp,clim[cl],"fullharm with Iizumi\n")
            if(length(which(!is.na(dh.r)))>1) {
              taylor.diagram2(as.vector(r.ray[ray.season,fpus]),as.vector(dh.r),as.vector(ww.r[,fpus]),add=T,col="magenta",pch=which(ggcms==gg),normalize=do.norm)
            } else cat("skipping",gg,cc,topp,clim[cl],"harmnon with Ray\n")
            if(length(which(!is.na(dh.i)))>1) {
              taylor.diagram2(as.vector(r.iizumi[iizumi.season,fpus]),as.vector(dh.i),as.vector(ww.i[,fpus]),add=T,col="yellow2",pch=which(ggcms==gg),normalize=do.norm)
            } else cat("skipping",gg,cc,topp,clim[cl],"harmnon with Iizumi\n")
        #}# if ric and pegasus or papsim
          
        }# for ggcm
        for(ee in ensembles){
          #get_nc4_ensemble_slice <- function(fname,crd=cr0,dtd=dt0,mpd=mp0,scend=scen0,nmd=nm0,wtd=wt0){
          fn <- paste(path,processed.ts,"/modelensemble/fpu/ray/",agg.ray,"/",ee,"_",clim[cl],"_hist_",cc,"_annual_",if("soy"==cc) cfirst[cl]-1 else cfirst[cl],"_",
                      clast[cl],".ensemble.nc4",sep="")
          data.ensemble.r <- if(is.finite(file.info(fn)$size>0))  get_nc4_ensemble_slice(fname=fn,nmd=nm0,wtd=wt0) else 0
          fn <- paste(path,processed.ts,"/modelensemble/fpu/iizumi/",agg.iizumi,"/",ee,"_",clim[cl],"_hist_",cc,"_annual_",if("soy"==cc) cfirst[cl]-1 else cfirst[cl],"_",
                      clast[cl],".ensemble.nc4",sep="")
          data.ensemble.i <- if(is.finite(file.info(fn)$size>0))  get_nc4_ensemble_slice(fname=fn,nmd=nm0,wtd=wt0) else 0
          de.r <- if(length(data.ensemble.r)>1) data.ensemble.r[clim.season.r,fpus] else NA
          de.i <- if(length(data.ensemble.i)>1) data.ensemble.i[clim.season.i,fpus] else NA
          
          if(length(which(!is.na(de.r)))>1) {
            taylor.diagram2(as.vector(r.ray[ray.season,fpus]),as.vector(de.r),as.vector(ww.r[,fpus]),add=T,col="blue",pch=if(ensembles[1]==ee) 17 else 15,normalize=do.norm)
          } else cat("skipping ensemble for Ray\n")
          if(length(which(!is.na(de.i)))>1) {
            taylor.diagram2(as.vector(r.iizumi[iizumi.season,fpus]),as.vector(de.i),as.vector(ww.i[,fpus]),add=T,col="red",pch=if(ensembles[1]==ee) 17 else 15,normalize=do.norm)
          } else cat("skipping ensemble for Iizumi\n")
          
        }
        #       legend(pos* if(topp==topps[1]).9 else .8,pos* if(topp==topps[1]) 1.3 else 1.4 ,legend=c(ggcms,"Iizumi vs. Ray","default vs. Ray","fullharm vs. Ray","default vs. Iizumi", "fullharm vs. Iizumi","harmnon vs. Ray","harmnon vs. Iizumi"),
        #              pch=c(1:length(ggcms),16,16,16,16,16,16,16),col=c(rep(1,length(ggcms)),"grey60","red","blue","green","orange","magenta","yellow2"),bty="n",cex=0.7)
        legend(-2.9,3.9 ,legend=c(ggcms),
               pch=c(1:length(ggcms)),col=c(rep(1,length(ggcms))),bty="n",cex=0.7,xjust=0)
        legend(2.9,3.9 ,legend=c("default vs. Ray","fullharm vs. Ray","harmnon vs. Ray",
                                 "default vs. Iizumi", "fullharm vs. Iizumi","harmnon vs. Iizumi",
                                 "Iizumi vs. Ray","","RMSE ensemble vs.Ray","tscorr ensemble vs.Ray","RMSE ensemble vs. Iizumi","tscorr ensemble vs. Iizumi"),
               pch=c(16,16,16,16,16,16,16,1,17,15,17,15),col=c("blue","red","magenta","green","orange","yellow2","grey60","white","blue","blue","red","red"),
               bty="n",cex=0.7,xjust=1)
        
        dev.off()
      }
    }
  }
  
}


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
    #for(cc in "mai"){
      fname <- "D:/data/GGCMI/reference/faostat/faostat.1961-2012.gadm0.nc4"
      # get valid gadm0 IDs
      nf <- nc_open(fname)
      gadm0.id <- ncvar_get(nf,varid="gadm0")
      nc_close(nf)
      fname <- if(agg.ray=="dynamic_ray_mask") "D:/data/GGCMI/reference/ray/ray.1961-2008.gadm0.ray.nc4" else "D:/data/GGCMI/reference/ray/ray.1961-2008.gadm0.fixed.nc4"
      # get valid gadm0 IDs
      nf <- nc_open(fname)
      gadm0.id.r <- ncvar_get(nf,varid="gadm0")
      nc_close(nf)
      fname <- if(agg.ray=="dynamic_ray_mask") "D:/data/GGCMI/reference/iizumi/iizumi.1982-2006.gadm0.iizumi.nc4" else "D:/data/GGCMI/reference/iizumi/iizumi.1982-2006.gadm0.fixed.nc4"
      # get valid gadm0 IDs
      nf <- nc_open(fname)
      gadm0.id.i <- ncvar_get(nf,varid="gadm0")
      nc_close(nf)
      # fu*k. different country selection in sim files
      fname<-paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg.fao,"/","pdssat","_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                  clast[cl],".biascorr.nc4",sep="")
      nf <- nc_open(fname)
      gadm0.id2 <- ncvar_get(nf,varid="gadm0")
      nc_close(nf)
      sset <- which(gadm0.index %in% gadm0.id & gadm0.index %in% gadm0.id2)

      fname <- "D:/data/GGCMI/reference/faostat/faostat.1961-2012.gadm0.nc4"
      r.fao <- get_nc4_ref_slice(fname,cc)[,which(gadm0.id %in% gadm0.id2)]
      fname <- "D:/data/GGCMI/reference/ray/ray.1961-2008.gadm0.ray.nc4"
      r.ray <- get_nc4_ref_slice(fname,cc)[,which(gadm0.id.r %in% gadm0.id2)]
      fname <- "D:/data/GGCMI/reference/iizumi/iizumi.1982-2006.gadm0.iizumi.nc4"
      r.iizumi <- get_nc4_ref_slice(fname,cc)[,which(gadm0.id.i %in% gadm0.id2)]

      gadm0.names <- name[sset]
      gadm0.id <- 1:length(gadm0.names)
      fyear.plot <- max(cfirst[cl],faof)+ignore.y
      lyear.plot <- min(clast[cl],faol)-ignore.y
#       # read area per country
#       nf <- nc_open(paste("D:/data/GGCMI/",processed.ts,"/masks/weight/aggs/",tolower(cropsl[which(crops==cc)]),".fpu.nc4",sep=""))
#       area_per_fpu_rf <- ncvar_get(nf,varid="rainfed_fpu")
#       area_per_fpu_ir <- ncvar_get(nf,varid="irrigated_fpu")
#       nc_close(nf)
#       ww.r <- matrix(rep(cbind(area_per_fpu_rf+area_per_fpu_ir),length(clim.season.r)),nrow=length(clim.season.r),byrow=T)
#       ww.i <- matrix(rep(cbind(area_per_fpu_rf+area_per_fpu_ir),length(clim.season.i)),nrow=length(clim.season.i),byrow=T)
      
      if(F){ # some more checking
        fyear.plot1 <- max(cfirst[cl],faof,rayf,iizumif) + ignore.y
        lyear.plot1 <- min(clast[cl],faol,rayl,iizumil) - ignore.y
        randomcountries <- c(1:208)
        for(rci in randomcountries){
          prange <- range(r.ray[c(fyear.plot1:lyear.plot1)-rayf+1,rci],
                          r.iizumi[c(fyear.plot1:lyear.plot1)-iizumif+1,rci],
                          r.fao[c(fyear.plot1:lyear.plot1)-rayf+1,rci],na.rm=T)
          if(!is.infinite(prange[1])){
            png(paste("D:/data/GGCMI/reference/comparison_",gadm0.names[rci],"_gadm0_",fyear.plot,"-",lyear.plot,"_",prefix,Sys.Date(),".png",sep=""))
            #randomcountries <- which(!is.na(r.ray[1,]))[9:10]
            # check Italy, which is FPU==full country
            #rc2 <- which(gadm0.id.r==randomcountries)
            plot(c(fyear.plot1:lyear.plot1),r.ray[c(fyear.plot1:lyear.plot1)-rayf+1,rci],
                 ylim=prange,
                 type="n",ylab="yield [t/ha]",main=paste(gadm0.names[rci],fyear.plot1,"-",lyear.plot1,prefix,sep=" "))
            
            lines(c(fyear.plot1:lyear.plot1),r.ray[c(fyear.plot1:lyear.plot1)-rayf+1,rci],col=1)
            lines(c(fyear.plot1:lyear.plot1),r.iizumi[c(fyear.plot1:lyear.plot1)-iizumif+1,rci],col=2,lty=2)
            lines(c(fyear.plot1:lyear.plot1),r.fao[c(fyear.plot1:lyear.plot1)-faof+1,rci],col=3,lty=3)
            legend("topleft",legend=c("Ray","Iizumi","FAO"),lty=c(1,2,3),col=c(1,2,3))
            dev.off()
            
          }
        }
        png(paste("D:/data/GGCMI/reference/comparison_mean_per_GADM0_",fyear.plot,"-",lyear.plot,"_",prefix,Sys.Date(),".png",sep=""),
            width=5*1200,height=3*1200,res=1200,pointsize=4)
        #randomcountries <- which(!is.na(r.ray[1,]))[9:10]
        # check Italy, which is FPU==full country
        diff.ref <- colMeans(r.ray[c(fyear.plot:lyear.plot)-rayf+1,])-colMeans(r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,])
        prange <- range(r.ray[c(fyear.plot:lyear.plot)-rayf+1,],diff.ref,
                        r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,],na.rm=T)
        plot(1:dim(r.ray)[2],colMeans(r.ray[c(fyear.plot:lyear.plot)-rayf+1,]),
             ylim=prange,
             type="p",ylab="yield [t/ha]",main=paste("FPU means",fyear.plot,"-",lyear.plot,prefix,sep=" "),
             col=2,pch=3,axes=F)
        axis(2)
        #axis(1,crt=90,labels=fpu.names,at=c(1:309),cex=.4)
        text(x=seq(1,by=1,length.out=dim(r.ray)[2]),y=min(diff.ref,na.rm=T),labels=gadm0.names,
             adj=c(1,.5),cex=0.3,srt=90,xpd=NA)
        
        points(1:dim(r.ray)[2],colMeans(r.iizumi[c(fyear.plot:lyear.plot)-iizumif+1,]),
               col=3,pch=3)
        points(1:dim(r.ray)[2],diff.ref,
               col=4,pch=4)
        for(i in 1:dim(r.ray)[2]) lines(c(i,i),prange,lwd=0.1)
        for(i in c(as.integer(prange[1]):as.integer(prange[2]))) lines(c(1,dim(r.ray)[2]),c(i,i),lwd=0.1,lty=2)
        legend("topleft",legend=c("Ray","Iizumi","difference"),col=c(2,3,4),pch=c(3,3,4))
        dev.off()
      }
      # read area per gadm0
      #     nf <- nc_open("D:/data/GGCMI/processed/masks/weight/aggs/maize.agg.nc4")
      #     area_per_gadm0_rf <- ncvar_get(nf,varid="rainfed_gadm0")
      #     area_per_gadm0_ir <- ncvar_get(nf,varid="irrigated_gadm0")
      #     nc_close(nf)
      while(length(dev.list())>0) dev.off() # make sure no figures are open to allow for 2 simultaneous figures
      for(topp in topps){
        if(do.png){
          png(paste(path,"taylor/",cc,".taylor.",topp,".gadm0.",clim[cl],".",prefix,".",agg.fao,".",Sys.Date(),".png",sep=""),height=5*300,width=5*300,res=300,pointsize=9)
        } else{
          postscript(paste(path,"taylor/",cc,".taylor.",topp,".gadm0.",clim[cl],".",prefix,".",agg.fao,".",Sys.Date(),".eps",sep=""),height=5*100,width=5*100,pointsize=22,family="NimbusSan")
        }
        par(xpd=NA)
        if(topp=="global"){
          gadm0s <- gadm0.id
        }else{
          # this is for country level
          #gadm0s <- if(cc=="whe") gadm.whe[2:11] else if(cc=="mai") gadm.mai[2:11] else if(cc=="ric") gadm.ric[2:11] else if(cc=="soy") gadm.soy[2:11]
          # this is a proxy for gadm0 level, just using the 10 highest productivity ones
          #gadm0s <- order(colMeans(r.fao),decreasing=T)[1:10]
          gadm0s <- which(gadm0.names %in% if(cc=="whe") topproducer.whe[2:11] else if(cc=="mai") topproducer.mai[2:11] else if(cc=="ric") topproducer.ric[2:11] else topproducer.soy[2:11])
        }
        # plot fao vs. iizumi in taylor diagram
        pos <- taylor.diagram2(as.vector(r.fao[c(fyear.plot:lyear.plot)-faof+1,gadm0s]),
                               as.vector(r.fao[c(fyear.plot:lyear.plot)-faof+1,gadm0s]),ref.sd=T,
                               show.gamma=F,normalize=do.norm,col="grey60",pos.cor=F,#pos.cor=if(topp=="global") T else F,
                               main=paste(topp,clim[cl],cropsl[which(crops==cc)],"\n",prefix,"\n",agg.fao),type="n")
        
        # loop through models and add points to taylor diagram
        for(gg in ggcms){
          cat(gg,": ",sep="")
          fn <- paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg.fao,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                      if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
          data.default.f <- if(is.finite(file.info(fn)$size)) get_nc4_data_slice(fname=fn,scend="default") else 0
          fn <- paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg.fao,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                      if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
          data.fullharm.f <- if(is.finite(file.info(fn)$size)) get_nc4_data_slice(fname=fn,scend="fullharm") else 0
          fn <- paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg.fao,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                      if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
          data.harmnon.f <- if(is.finite(file.info(fn)$size)) get_nc4_data_slice(fname=fn,scend="harmnon") else 0
          dd.f <- if(length(data.default.f)>1) data.default.f[clim.season.f,gadm0s] else NA
          df.f <- if(length(data.fullharm.f)>1) data.fullharm.f[clim.season.f,gadm0s] else NA
          dh.f <- if(length(data.harmnon.f)>1) data.harmnon.f[clim.season.f,gadm0s] else NA
          
          # plot time series of models for top10
          if(topp==topps[2]){
            ref.vec <- as.vector(r.fao[fao.season,gadm0s])
            sim.vec.dd <- as.vector(dd.f)
            sim.vec.df <- as.vector(df.f)
            sim.vec.dh <- as.vector(dh.f)
            png(paste(path,"timeseries/",gg,"/",gg,".",cc,".timeseries.",topp,".gadm0.jointly.",clim[cl],".",prefix,".",agg.fao,".",Sys.Date(),".png",sep=""),height=3*300,width=5*300,res=300,pointsize=9)
            par(cex=0.6)
            prange <- range(ref.vec,sim.vec.dd,sim.vec.dh,sim.vec.df,na.rm=T)
            plot(ref.vec,ylim=prange,type="o",cex=0.6,xlab="",ylab="")
            text(x=seq(1,by=length(fao.season),length.out=length(gadm0s)),y=prange[2]*1.1,labels=gadm0.names[gadm0s],
                 adj=0,cex=0.6,srt=20,xpd=NA)
            lines(sim.vec.dd,type="o",cex=0.6,lwd=0.5,col="green")
            lines(sim.vec.df,type="o",cex=0.6,lwd=0.5,col="orange")
            lines(sim.vec.dh,type="o",cex=0.6,lwd=0.5,col="yellow2")
            legend("topleft",legend=c("Iizumi",paste("default r: ",if(length(which(!is.na(sim.vec.dd)))>0) round(cor(ref.vec,sim.vec.dd,"pairwise"),3)),
                                      paste("fullharm r: ",if(length(which(!is.na(sim.vec.df)))>0) round(cor(ref.vec,sim.vec.df,"pairwise"),3)),
                                      paste("harmnon r: ",if(length(which(!is.na(sim.vec.dh)))>0) round(cor(ref.vec,sim.vec.dh,"pairwise"),3))),
                   col=c(1,"green","orange","yellow2"),lty=1,bty="n")
            dev.off()
            png(paste(path,"timeseries/",gg,"/",gg,".",cc,".timeseries.",topp,".gadm0.",clim[cl],".",prefix,".",agg.fao,".",Sys.Date(),".png",sep=""),height=10*300,width=5*300,res=300,pointsize=9)
            split.screen(c(5,2))
            for(i in 1:length(gadm0s)){
              screen(i)
              par(cex=0.7,mar=c(2,2,2.5,0))
              if(length(which(!is.na(r.fao[fao.season,gadm0s[i]])))>0 ){
                drange <- range(r.fao[fao.season,gadm0s[i]],
                                if(length(dim(dd.f))>0)dd.f[,i],
                                if(length(dim(df.f))>0)df.f[,i],
                                if(length(dim(dh.f))>0)dh.f[,i],na.rm=T)
                plot(c(faof:faol)[fao.season],r.fao[fao.season,gadm0s[i]],type="o",ylim=drange,pch=16,xlab="",ylab="",cex=0.6,
                     main=gadm0.names[gadm0s[i]])  
                if(length(dim(dd.f))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.f],dd.f[,i],col=rainbow(10)[2],type="o",pch=16,cex=0.6)
                if(length(dim(df.f))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.f],df.f[,i],col=rainbow(10)[4],type="o",pch=16,cex=0.6)
                if(length(dim(dh.f))>0)lines(c(cfirst[cl]:clast[cl])[clim.season.f],dh.f[,i],col=rainbow(10)[6],type="o",pch=16,cex=0.6)
                par(cex=0.6)
                legend("topleft",bty="n",legend=c("fao","default","fullharm","harmnon"),
                       col=c(1,rainbow(10)[c(2,4,6)]),lty=c(1,1,1,1))
                
              }
            }
            close.screen(all=T)
            
            dev.off()
          }
          
          if(length(which(!is.na(dd.f)))>1) {
            taylor.diagram2(as.vector(r.fao[fao.season,gadm0s]),as.vector(dd.f),add=T,col="blue",pch=which(ggcms==gg),normalize=do.norm)
          } else cat("skipping",gg,cc,topp,clim[cl],"default with fao\n")
          if(length(which(!is.na(df.f)))>1) {
            taylor.diagram2(as.vector(r.fao[fao.season,gadm0s]),as.vector(df.f),add=T,col="red",pch=which(ggcms==gg),normalize=do.norm)
          } else cat("skipping",gg,cc,topp,clim[cl],"fullharm with fao\n")
          if(length(which(!is.na(dh.f)))>1) {
            taylor.diagram2(as.vector(r.fao[fao.season,gadm0s]),as.vector(dh.f),add=T,col="magenta",pch=which(ggcms==gg),normalize=do.norm)
          } else cat("skipping",gg,cc,topp,clim[cl],"harmnon with fao\n")          
        } # for ggcms
        for(ee in ensembles){
          #get_nc4_ensemble_slice <- function(fname,crd=cr0,dtd=dt0,mpd=mp0,scend=scen0,nmd=nm0,wtd=wt0){
          fn <- paste(path,processed.ts,"/modelensemble/gadm0/ray/",agg.ray,"/",ee,"_",clim[cl],"_hist_",cc,"_annual_",if("soy"==cc) cfirst[cl]-1 else cfirst[cl],"_",
                      clast[cl],".ensemble.nc4",sep="")
          data.ensemble.r <- if(is.finite(file.info(fn)$size>0))  get_nc4_ensemble_slice(fname=fn,nmd=nm0,wtd=wt0) else 0
          fn <- paste(path,processed.ts,"/modelensemble/gadm0/iizumi/",agg.ray,"/",ee,"_",clim[cl],"_hist_",cc,"_annual_",if("soy"==cc) cfirst[cl]-1 else cfirst[cl],"_",
                      clast[cl],".ensemble.nc4",sep="")
          data.ensemble.i <- if(is.finite(file.info(fn)$size>0))  get_nc4_ensemble_slice(fname=fn,nmd=nm0,wtd=wt0) else 0
          de.r <- if(length(data.ensemble.r)>1) data.ensemble.r[clim.season.r,gadm0s] else NA
          de.i <- if(length(data.ensemble.i)>1) data.ensemble.i[clim.season.i,gadm0s] else NA
          
          if(length(which(!is.na(de.r)))>1) {
            taylor.diagram2(as.vector(r.ray[ray.season,gadm0s]),as.vector(de.r),add=T,col="blue",pch=if(ensembles[1]==ee) 17 else 15,normalize=do.norm)
          } else cat("skipping ensemble for Ray\n")
          if(length(which(!is.na(de.i)))>1) {
            taylor.diagram2(as.vector(r.iizumi[iizumi.season,gadm0s]),as.vector(de.i),add=T,col="red",pch=if(ensembles[1]==ee) 17 else 15,normalize=do.norm)
          } else cat("skipping ensemble for Ray\n")
          
        }
        
#         legend(pos* if(topp==topps[1]).9 else .8,pos* if(topp==topps[1]) 1.3 else 1.4 ,legend=c(ggcms,"Iizumi vs. fao","default vs. fao","fullharm vs. fao","default vs. Iizumi", "fullharm vs. Iizumi","harmnon vs. fao","harmnon vs. Iizumi"),
#                pch=c(1:length(ggcms),16,16,16,16,16,16,16),col=c(rep(1,length(ggcms)),"grey60","red","blue","green","orange","magenta","yellow2"),bty="n",cex=0.7)
#         legend(-2.8,3.5 ,legend=c(ggcms),
#                pch=c(1:length(ggcms)),col=c(rep(1,length(ggcms))),bty="n",cex=0.7,xjust=0)
#         legend(2.8,3.5 ,legend=c("default vs. fao","fullharm vs. fao","harmnon vs. fao"),
#                pch=c(16,16,16),col=c("blue","red","magenta"),
#                bty="n",cex=0.7,xjust=1)
        legend(-2.9,3.9 ,legend=c(ggcms),
               pch=c(1:length(ggcms)),col=c(rep(1,length(ggcms))),bty="n",cex=0.7,xjust=0)
        legend(2.9,3.9 ,legend=c("default vs. FAOstat","fullharm vs. FAOstat","harmnon vs. FAOstat",
                                 "","RMSE ensemble vs.Ray","tscorr ensemble vs.Ray","RMSE ensemble vs. Iizumi","tscorr ensemble vs. Iizumi"),
               pch=c(16,16,16,1,17,15,17,15),col=c("blue","red","magenta","white","blue","blue","red","red"),
               bty="n",cex=0.7,xjust=1)
        
        dev.off()
      }
    }
  }
  
}

#### do.pixel ####
#### pixel-wise comparison
if(do.pixel){ #do pixel-wise comparison for selected countries
  prefix <- paste("cr_",cr0,".dt_",dt0,".",sep="")
  

  for(cl in 1){ # AgMERRA only
    if(do.raw & cl>1) cat("WARNING! raw works for AgMERRA only so far\n!")
    ray.season <- c(max(cfirst[cl],rayf):min(clast[cl],rayl)) - rayf+1
    iizumi.season <- c(max(cfirst[cl],iizumif):min(clast[cl],iizumil)) - iizumif+1
    clim.season.r <- c(max(rayf,cfirst[cl]) :min(rayl,clast[cl])) - max(1961,cfirst[cl])+1  # rescaled data is from 1961 onwards only, even if files are names 1958 (WATCH)
    clim.season.i <- c(max(iizumif,cfirst[cl]):min(iizumil,clast[cl])) - max(1961,cfirst[cl])+1 
    #  cat(ray.season,"\n",iizumi.season,"\n",clim.season.r,"\n",clim.season.i,"\n")
    # cut away first and last 2 years for MA
    # remove first and last years
    ray.season <- c(ray.season[1+ignore.y]:ray.season[length(ray.season)-ignore.y])
    iizumi.season <- c(iizumi.season[1+ignore.y]:iizumi.season[length(iizumi.season)-ignore.y])
    clim.season.i <- c((clim.season.i[1+ignore.y]):(clim.season.i[length(clim.season.i)-ignore.y]))    
    clim.season.r <- c((clim.season.r[1+ignore.y]):(clim.season.r[length(clim.season.r)-ignore.y]))    
    
    
    for(cc in crops){
      if(cc=="whe") {topprod <-  topproducer.whe;gadm <- gadm0.whe
      }else if(cc=="mai") {topprod <-  topproducer.mai;gadm <- gadm0.mai
      }else if(cc=="ric") {topprod <-  topproducer.ric;gadm <- gadm0.ric
      }else if(cc=="soy") {topprod <-  topproducer.soy;gadm <- gadm0.soy
      }
      #for(j in 1:length(topprod)){
      for(j in 2:11){
        topp <- topprod[j]
        gad <- gadm[j]
        png(paste(path,"taylor/",cc,".taylor.",if(!do.rescale)"detrended.",topp,".",clim[cl],".",prefix,".r_vs_",agg.ray,".i_vs_",agg.iizumi,".",Sys.Date(),".png",sep=""),height=5*300,width=5*300,res=300,pointsize=9)
        
        #ref data
        
        #nf  <- nc_open(paste(path.ref,cc,if(cc=="mai") "_weight","_ray_1961-2008.nc4",sep=""))
        nf  <- nc_open(paste(path.ref,cc,"_weight_ray_1961-2008.nc4",sep=""))
        r.ray <- ncvar_get(nf,varid=paste("yield_",cc,sep=""))
        #if(cc=="mai"){
        
        area <- ncvar_get(nf,varid=paste("area_",cc,sep=""))
        r.ray[area<2000] <- NA
        r.iizumi[area<2000] <- NA
        
        #}
        #ref <- ref[,,ray.season]
        #ref[!is.finite(ref)] <- 0
        nc_close(nf)
        # detrended Ray data
        nf  <- nc_open(paste(path.ref.detrend,cc,"_weight_ray_1961-2008.nc4",sep=""))
        r.ray.dt <- ncvar_get(nf,varid=paste("yield_",cc,sep=""))[1,3,,,] #mp=T, dt=quad(3) or ma(4)
        r.ray.dt[area<2000] <- NA
        nc_close(nf)
        
        mask <- r.ray[,,1]
        mask[] <- NA
        if(gad==gadm[1]) mask[] <- 1 else mask[gadm.mask==gad] <- 1
        
        
        nf  <- nc_open(paste("D:/data/GGCMI/iizumi/iizumi.2013JAN29.",cropsi[which(crops==cc)],".1982-2006.30min.nc4",sep=""))
        r.iizumi <- ncvar_get(nf,varid="yield50")
        # cut away 2006
        #ref2 <- ref2[,,iizumi.season]
        #ref2[!is.finite(ref2)] <- 0
        nc_close(nf)
        nf  <- nc_open(paste("D:/data/GGCMI/iizumi/detrended/iizumi.2013JAN29.",cropsi[which(crops==cc)],".1982-2006.30min.nc4",sep=""))
        r.iizumi.dt <- ncvar_get(nf,varid="yield50")[1,3,,,]#mp=T, dt=quad(3) or ma(4)
        nc_close(nf)
        
        
        
        
        par(xpd=NA)
        # plot taylor diagram
        #nonzero2 <- apply(ref,MARGIN=c(1,2),sum)
        #nonzero3 <- apply(ref2,MARGIN=c(1,2),sum)
        sset <- as.logical(r.ray)
        sset[] <- F
        #sset[which(is.finite(nonzero2) & is.finite(nonzero3) & is.finite(mask))] <- T
        sset[which(is.finite(mask))] <- T
        if(do.rescale){
          pos <- taylor.diagram2(subset(r.ray[,,c(max(iizumif,rayf):min(iizumil,rayl))-rayf+1],sset),subset(r.iizumi[,,c(max(iizumif,rayf):min(iizumil,rayl))-iizumif+1],sset),ref.sd=T,
                                 show.gamma=F,normalize=do.norm,col="grey60",pos.cor=if(topp=="global") T else F,
                                 main=paste(topp,clim[cl],cropsl[which(crops==cc)]))
        } else
        {
          pos <- taylor.diagram2(subset(r.ray.dt[,,c(max(iizumif,rayf):min(iizumil,rayl))-rayf+1],sset),subset(r.iizumi.dt[,,c(max(iizumif,rayf):min(iizumil,rayl))-iizumif+1],sset),ref.sd=T,
                                 show.gamma=F,normalize=do.norm,col="grey60",pos.cor=if(topp=="global") T else F,
                                 main=paste(topp,clim[cl],cropsl[which(crops==cc)]))
        }
        #taylor.diagram(subset(ref,sset),subset(data.rescale,sset),add=T,col="blue")
        
        # now, after comparing Ray and Iizumi, cut ref data to relevant period
        r.iizumi <- r.iizumi[,,iizumi.season]
        r.ray <- r.ray[,,ray.season]
        r.iizumi.dt <- r.iizumi.dt[,,iizumi.season]
        r.ray.dt <- r.ray.dt[,,ray.season]
        
        if(do.raw){
          # read weights for combining raw irrig and rainfed data
          nf <- nc_open(paste(path,processed.ts,"/masks/weight/landuse.ir.nc4",sep=""))
          ir <- ncvar_get(nf,varid=cropsl[which(crops==cc)])
          nc_close(nf)
          nf <- nc_open(paste(path,processed.ts,"/masks/weight/landuse.rf.nc4",sep=""))
          rf <- ncvar_get(nf,varid=cropsl[which(crops==cc)])
          nc_close(nf)
          rf[rf==0] <- NA
          ir[ir==0] <- NA        
        }
        
        for(gg in ggcms[1:8]){ # I currently only have 7 GGCMs on my local drive...
          if(!((cc=="ric" & gg=="papsim")|(cc=="ric" & gg=="pegasus"))){ #no rice for papsim and pegasus
            #           if(do.rescale){
            #           nf <- nc_open(paste(path,"processed.140821/rescaled/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
            #                               if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".rescaled.nc4",sep=""))
            #           data.rescale <- ncvar_get(nf)
            #           # sum of firr and noirr, default
            #           if(length(dim(data.rescale))==5){ # default is not always there. If there's only one dimension, it will be fullharm/harmnon
            #             dr <- data.rescale[3,1,,,clim.season.r] # Ray has 1980-2008
            #             dr2 <- data.rescale[3,1,,,clim.season.i] # Iizumi has 1982-2006
            #             #data.rescale[!is.finite(data.rescale)] <- 0
            #             #nonzero <- apply(dr,MARGIN=c(1,2),sum)
            #             sset2 <- sset <- as.logical(ref)
            #             sset[] <- F
            #             #sset[which(is.finite(nonzero) & is.finite(nonzero2) & is.finite(mask))] <- T
            #             sset[is.finite(mask)] <- T
            #             if(length(subset(ref,sset))>0 & length(subset(dr,sset))>0) {
            #               taylor.diagram2(subset(ref,sset),subset(dr,sset),add=T,col="red",pch=which(ggcms==gg),normalize=do.norm)
            #             } else cat("skipping",gg,cc,topp,clim[cl],length(subset(ref,sset)),length(subset(dr,sset)),"default with Ray\n")
            #             sset2[] <- F
            #             #sset2[which(is.finite(nonzero) & is.finite(nonzero3) & is.finite(mask))] <- T
            #             sset2[is.finite(mask)] <- T
            #             if(length(subset(ref2,sset2))>0 & length(subset(dr2,sset2))>0) {
            #               taylor.diagram2(subset(ref2,sset2),subset(dr2,sset2),add=T,col="green",pch=which(ggcms==gg),normalize=do.norm)
            #             } else cat("skipping",gg,cc,topp,clim[cl],length(subset(ref2,sset2)),length(subset(dr2,sset2)),"deffault with Iizumi\n")
            #           }
            #           # sum of firr and noirr, fullharm/harmnon
            #           if(length(dim(data.rescale))==5){
            #             dr <- data.rescale[3,2,,,clim.season.r] # Ray has 1980-2008
            #             dr2 <- data.rescale[3,2,,,clim.season.i] # Iizumi has 1982-2006
            #           } else {
            #             dr <- data.rescale[3,,,clim.season.r]
            #             dr2 <- data.rescale[3,,,clim.season.i] # Iizumi has 1982-2006, dropping 2006 for bad quality            
            #           }
            #           #data.rescale[!is.finite(data.rescale)] <- 0
            #           if(length(subset(ref,sset))>0 & length(subset(dr,sset))>0) {
            #             taylor.diagram2(subset(ref,sset),subset(dr,sset),add=T,col="blue",pch=which(ggcms==gg),normalize=do.norm)
            #           } else cat("skipping",gg,cc,topp,clim[cl],length(subset(ref,sset)),length(subset(dr,sset)),"fullharm with Ray\n")
            #           if(length(subset(ref2,sset2))>0 & length(subset(dr2,sset2))>0) {
            #             taylor.diagram2(subset(ref2,sset2),subset(dr2,sset2),add=T,col="orange",pch=which(ggcms==gg),normalize=do.norm)
            #           } else cat("skipping",gg,cc,topp,clim[cl],length(subset(ref2,sset2)),length(subset(dr2,sset2)),"fullharm with Iizumi\n")
            #           nc_close(nf)
            #           } #do.rescale
            # raw data
            if(do.raw & gg!="epic-iiasa"){
              fn <- paste(path,"raw/",gg,"_agmerra_hist_",if(gg=="lpjml" | gg=="lpj-guess")"harmnon" else "fullharm","_noirr_yield_",cc,"_annual_1980_2010.nc4",sep="")
              if(is.finite(file.info(fn)$size)){
                nf <- nc_open(fn)
                yrf <- ncvar_get(nf)
                nc_close(nf)
              } else yrf <- NA
              fn <- paste(path,"raw/",gg,"_agmerra_hist_",if(gg=="lpjml" | gg=="lpj-guess")"harmnon" else "fullharm","_firr_yield_",cc,"_annual_1980_2010.nc4",sep="")
              if(is.finite(file.info(fn)$size)){
                nf <- nc_open(fn)
                yir <- ncvar_get(nf)
                nc_close(nf)
              } else yrf <- NA
              fn <-paste(path,"raw/",gg,"_agmerra_hist_default_noirr_yield_",cc,"_annual_1980_2010.nc4",sep="")
              if(is.finite(file.info(fn)$size)){
                nf <- nc_open(fn)
                yrf2 <- ncvar_get(nf)
                nc_close(nf)
              } else yrf <- NA
              fn <- paste(path,"raw/",gg,"_agmerra_hist_default_firr_yield_",cc,"_annual_1980_2010.nc4",sep="")
              if(is.finite(file.info(fn)$size)){
                nf <- nc_open(fn)
                yir2 <- ncvar_get(nf)
                nc_close(nf)
              } else yrf <- NA
              for(i in 1:dim(yrf)[3]) {
                yrf[,,i] <- yrf[,,i]*rf #harmonized
                yir[,,i] <- yir[,,i]*ir #harmonized
                yrf2[,,i] <- yrf2[,,i]*rf #default
                yir2[,,i] <- yir2[,,i]*ir #default
              }
              raw.sum <- yrf[,,clim.season.r]+yir[,,clim.season.r]
              raw.sum2 <- yrf2[,,clim.season.r]+yir2[,,clim.season.r]
              #rm(yrf,yir,yrf2,yir2)
              #nonzero <- apply(raw.sum,MARGIN=c(1,2),sum)
              sset2 <- sset <- as.logical(r.ray)
              sset[] <- F
              #sset[which(is.finite(nonzero) & is.finite(nonzero2) & is.finite(nonzero3))] <- T
              sset[is.finite(mask)] <- T
              sset2[] <- F
              sset2[is.finite(mask)] <- T
              #if(length(subset(r.ray.dt,sset))>0 & length(subset(raw.sum,sset))>0) 
              if(length(which(is.finite(subset(r.ray.dt,sset)) & is.finite(subset(raw.sum,sset))))>0)
                taylor.diagram2(subset(r.ray.dt,sset),subset(raw.sum,sset),add=T,col=if(gg=="lpjml" | gg=="lpj-guess") "magenta" else"red",pch=which(ggcms==gg),normalize=do.norm)
              #if(length(subset(r.ray.dt,sset))>0 & length(subset(raw.sum2,sset))>0)
              if(length(which(is.finite(subset(r.ray.dt,sset)) & is.finite(subset(raw.sum2,sset))))>0)
                taylor.diagram2(subset(r.ray.dt,sset),subset(raw.sum2,sset),add=T,col="blue",pch=which(ggcms==gg),normalize=do.norm)
              # cut away first and last 2 years for Iizumi
              raw.sum <- yrf[,,clim.season.i]+yir[,,clim.season.i]
              raw.sum2 <- yrf2[,,clim.season.i]+yir2[,,clim.season.i]
              rm(yrf,yir,yrf2,yir2)
              #if(length(subset(r.iizumi.dt,sset2))>0 & length(subset(raw.sum,sset2))>0) 
              if(length(which(is.finite(subset(r.iizumi.dt,sset2)) & is.finite(subset(raw.sum,sset2))))>0)
                taylor.diagram2(subset(r.iizumi.dt,sset2),subset(raw.sum,sset2),add=T,col=if(gg=="lpjml" | gg=="lpj-guess") "yellow2" else "orange",pch=which(ggcms==gg),normalize=do.norm)
              if(length(which(is.finite(subset(r.iizumi.dt,sset2)) & is.finite(subset(raw.sum2,sset2))))>0)
                taylor.diagram2(subset(r.iizumi.dt,sset2),subset(raw.sum2,sset2),add=T,col="green",pch=which(ggcms==gg),normalize=do.norm)          
            } # do.raw
          }
          cat("done",gg,"\n")
        }#ggcm
#         legend(pos* if(1==j)0.9 else .8,pos* if(1==j) 1.15 else 1.3,legend=c(ggcms,"Iizumi vs. Ray","default vs. Ray","fullharm vs. Ray","harmnon vs. Ray","default vs. Iizumi", "fullharm vs. Iizumi","harmnon vs. Iizumi", "Ray vs. Iizumi"),
#                pch=c(1:length(ggcms),16,16,16,16,16),col=c(rep(1,length(ggcms)),"grey60","red","blue","magenta","green","orange","yellow2","grey60"),bty="n",cex=0.7)
        legend(-2.8,3.5 ,legend=c(ggcms),
               pch=c(1:length(ggcms)),col=c(rep(1,length(ggcms))),bty="n",cex=0.7,xjust=0)
        legend(2.8,3.5 ,legend=c("default vs. Ray","fullharm vs. Ray","harmnon vs. Ray","default vs. Iizumi", 
                                 "fullharm vs. Iizumi","harmnon vs. Iizumi","Iizumi vs. Ray"),
               pch=c(16,16,16,16,16,16,16),col=c("blue","red","magenta","green","orange","yellow2","grey60"),
               bty="n",cex=0.7,xjust=1)
        
        dev.off()
        
      } # topprod    
    }# crop
  }# clim
}
