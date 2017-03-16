rm(list=ls(all=TRUE))


require(ncdf4)
options(warn=1)


# script to plot heatmap diagrams for GGCMI phase 1. 
# written by Christoph Müller, PIK
# cmueller@pik-potsdam.de
# I'm using a slightly modified version of heatmap.2() of the gplots package
# to allow for a second y-axis labeling

path <- "/p/projects/macmit/data/GGCMI/AgMIP.output/"
path.out <- "/p/projects/macmit/users/cmueller/GGCMI/heatmaps/"
path.ref <- "/p/projects/macmit/data/GGCMI/reference/"
path.ref.detrend <- "/p/projects/macmit/data/GGCMI/reference/ray-dt/"
path.script <- "/p/projects/macmit/users/cmueller/GGCMI/"
processed.ts <- "processed"

source(paste0(path.script,"heatmap.2b.R"))


# only works for agmerra so far
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
sig.threshold <- 0.1
shift.threshold <- 0.2

FRESHMATTER <- 100 / c(88, 88, 87, 91) 


aggs <- c("fixed_mirca_mask","dynamic_ray_mask","fixed_iizumi_mask","fixed_spam_mask")
dts <- c("none","lin","quad","ma","ffdtr")
scens <- c("default","fullharm","harmnon")
scens2 <- c("default","fullharm","harm-suffN")

# loop through some other models 
ggcms <- c("pdssat","epic-boku","epic-iiasa","gepic",
           "papsim","pegasus","lpj-guess","lpjml",
           "cgms-wofost","clm-crop","epic-tamu","orchidee-crop",
           "pepic","prysbi2")
ensembles <- c("mean.bias","tscorr")
# loop through crops
crops <- c("mai","whe","ric","soy")
cropsi <- c("maize_major","wheat","rice_major","soybean")
cropsl <- c("Maize","Wheat","Rice","Soy")
cropsl2 <- c("Maize","Wheat","Rice","Soybean")

#### functions ####

get_nc4_multimetrics <- function(fname,crd=cr0,dtd=dt0,mpd=mp0,scend=scen0,vid="mean.bias"){
  nf <- nc_open(fname)
  cr <- which(strsplit(ncatt_get(nf,varid="cr")$long_name,split=", ")[[1]]==crd)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  scen <- which(strsplit(ncatt_get(nf,varid="scen")$long_name,split=", ")[[1]]==scend)
  if(length(scen)==0){
    nc_close(nf)
    cat(scend,"not available.\n")
    return(NA)
  }
  if(length(strsplit(ncatt_get(nf,varid="scen")$long_name,split=", ")[[1]])>1) {
    # structure unclear for 3-element dimensions, but will use first element of both (full time range and no cr)
    data <- ncvar_get(nf,varid=vid)[1,1,mp,dt,scen,]
  }else{
    data <- ncvar_get(nf,varid=vid)[1,1,mp,dt,]
  }
  nc_close(nf)
  data
}

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
    return(array(NA,dim=c(31,208)))
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


terrain.colors.reverse <- function (n, alpha = 1) 
{
  if ((n <- as.integer(n[1L])) > 0) {
    k <- n%/%2
    h <- c(4/12, 2/12, 0/12)
    s <- c(1, 1, 0)
    v <- c(0.65, 0.9, 0.95)
    out <- c(hsv(h = seq.int(h[1L], h[2L], length.out = k), s = seq.int(s[1L], 
                                                                        s[2L], length.out = k), v = seq.int(v[1L], v[2L], 
                                                                                                            length.out = k), alpha = alpha), hsv(h = seq.int(h[2L], 
                                                                                                                                                             h[3L], length.out = n - k + 1)[-1L], s = seq.int(s[2L], 
                                                                                                                                                                                                              s[3L], length.out = n - k + 1)[-1L], v = seq.int(v[2L], 
                                                                                                                                                                                                                                                               v[3L], length.out = n - k + 1)[-1L], alpha = alpha))
    out[length(out):1]
  }
  else character()
}

# color scheme
color.red <- colorRampPalette(c("red4","red","orange","lightgoldenrod","yellow"))
color.green <- colorRampPalette(c("blue4","turquoise","green3","greenyellow"))
color.red2blue <- colorRampPalette(c("red4","red","orange","lightgoldenrod","lightblue","turquoise","blue","blue4"))


red.yellow.grey.green.blue <- function(n,alpha=1)
{
  if(n%%2 != 1) n <- n+1  
  nn <- n%/%2
  mid <- n%/%2 + 1
  col.veg <- color.red2blue(n)
  col.veg[mid] <- "grey90"
  col.veg[(mid+1):n] <- color.green(nn)[nn:1]
  col.veg[1:nn] <- color.red(nn)
  col.veg
}

# rebuild moving average from GGCMI processing pipeline detrender
winave <- function(d,nstart,nend){
  nd <- length(d)
  dave <- d
  numramp <- (nend-nstart)/2
  n <- nstart - 2
  for(i in 1:nd){ # python indexing starts at 0, R indexing at 1
    if((i-1)<= numramp){
      n <- n + 2
    } else if ((i-1)>=nd-numramp){
      n <- n-2
    }
    idx0 <- max(i - as.integer(n / 2), 0)
    idx1 <- min(i + as.integer(n / 2) + 1, nd)
    dave[i] <- mean(d[idx0 : idx1],na.rm=T)
    #dave <- dave[c(as.integer(n/2):as.integer(nd-n/2))+1]
  }
  dave[c(as.integer(n/2):as.integer(nd-n/2))+1] # again, move index by 1 to account for differences in pyhton and R
}

detrend_ma <- function(y){
  dy <- y
  dy[] <- NaN
  line <- dy
  # python indexing x[2 : -2] basically removes the first 2 and the last 2 elements
  # equivalent to x[3:(length(x)-2)] in R
  ly <- length(y) # helper for more efficient indexing
  line[3:(ly-2)] <- winave(y,5,7)
  dy[3:(ly-2)] <- y[3:(ly-2)] - line[3:(ly-2)]
  dy[!is.finite(dy)] <- NA
  dy
}

# detrend_cor function using the detrend_ma function
detrended_cor <- function(x1,x2){
  # check for computation of correlation coefficient | p-value
  x12 <- x1*x2
  if(length(x12[which(is.finite(x12))])<3 | sd(x1,na.rm=T)==0 | sd(x2,na.rm=T)==0)
    return(NA)
  x1_dt <- detrend_ma(x1)
  x2_dt <- detrend_ma(x2)
  x12 <- x1_dt*x2_dt
  if(length(x12[which(is.finite(x12))])<3 | sd(x1_dt,na.rm=T)==0 | sd(x2_dt,na.rm=T)==0)
    return(NA)
  co <- cor.test(x1_dt,x2_dt, use ="na.or.complete")
  if(co$p.value>sig.threshold)
    return(NA)
  as.numeric(co$estimate)
}

mean.bias <- function(x1,x2){
  x12 <- x1*x2
  if(length(x12[which(is.finite(x12))])<1)
    return(NA)
  mean.bias <- mean((x2-x1),na.rm=TRUE)
  mean.bias
}

#### define settings ####

cr0 <- "none"#"mean-scale"#"variance-scale"#"none"
dt0 <- "ma" #"none"#"quad"
mp0 <- "true"
scen0 <- "default"
nm0 <- 1 #number of top model for the ensemble
wt0 <- "unweighted" #wheighting for ensemble (1:unweighted, 2:weighted)

# ignore one year at beginning and end of time series
ignore.y <- 1

# read GADM0 names
gadm0.names <- array("missing",253)
gadm0.index <- read.csv(paste(path,"processed.150505","/masks/aggr/gadm0.meta.csv",sep=""),header=F)[,1]
name <- read.csv(paste(path,"processed.150505","/masks/aggr/gadm0.meta.csv",sep=""),header=F)[,2]

# f#*k. different country selection in sim files
fname<-paste(path,processed.ts,"/multimetrics/gadm0/faostat/fixed_mirca_mask/","pdssat","_",clim[1],"_hist_mai_annual_",cfirst[1],"_",
             clast[1],".multimetrics.nc4",sep="")
nf <- nc_open(fname)
gadm0.id2 <- ncvar_get(nf,varid="gadm0")
nc_close(nf)
sset <- which(gadm0.index %in% gadm0.id2)
gadm0.names <- name[sset]
gadm0.id <- 1:length(gadm0.names)


topproducer.soy <- c("United States","Brazil","Argentina","China","India","Paraguay","Canada","Uruguay","Ukraine","Bolivia")
gadm0.soy <- c(240,32,11,48,105,175,41,242,237,27)
topproducer.ric <- c("China","India","Indonesia","Bangladesh","Vietnam","Thailand","Myanmar","Philippines","Brazil","Japan")
gadm0.ric <- c(48,105,106,19,247,226,154,177,32,114)
topproducer.mai <- c("United States","China","Brazil","Argentina","Mexico","India","Ukraine","Indonesia","France","South Africa")
gadm0.mai <- c(240,48,32,11,145,105,237,106,79,209)
topproducer.whe <- c("China","India","United States","Russia","France","Canada","Australia","Pakistan","Germany","Turkey")
gadm0.whe <- c(48,105,240,186,79,41,14,170,86,232)



colo <- red.yellow.grey.green.blue

prefix <- paste("cr_",cr0,".dt_",dt0,".mp_",mp0,".multimetrics.mean_bias",sep="")

for(cc in crops){
  hm0 <- array(NA,dim=c(length(gadm0.names),length(ggcms)*3,length(aggs)))
  hm02 <- array(NA,dim=c(length(gadm0.names),length(ggcms)*3,length(aggs)))
  for(cl in 1){        
    for(gg in ggcms){
      cat("doing",gg,"\n")
      for(agg in aggs){
        
        fname <- paste0(path.ref,"/faostat/faostat.1961-2012.gadm0.nc4")
        # get valid gadm0 IDs
        nf <- nc_open(fname)
        gadm0.id <- ncvar_get(nf,varid="gadm0")
        nc_close(nf)
        # f#*k. different country selection in sim files
        fname<-paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg,"/","pdssat","_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                     clast[cl],".biascorr.nc4",sep="")
        nf <- nc_open(fname)
        gadm0.id2 <- ncvar_get(nf,varid="gadm0")
        nc_close(nf)
        sset <- which(gadm0.index %in% gadm0.id & gadm0.index %in% gadm0.id2)
        
        fname <-  paste0(path.ref,"/faostat/faostat.1961-2012.gadm0.nc4")
        r.fao <- get_nc4_ref_slice(fname,cc)[,which(gadm0.id %in% gadm0.id2)]
        
        # for a fair comparison, always use the same (shorter) time frame until 2008
        fao.season <- c(max(rayf,cfirst[cl],faof):min(rayl,clast[cl],faol)) - faof+1
        clim.season.f <- c(max(rayf,faof,cfirst[cl]) :min(rayl,faol,clast[cl])) - cfirst[cl]+1
        # remove first and last years
        fao.season <- c(fao.season[1+ignore.y]:fao.season[length(fao.season)-ignore.y])
        clim.season.f <- c((clim.season.f[1+ignore.y]):(clim.season.f[length(clim.season.f)-ignore.y]))    

        fn <- paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                    if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
        if(file.exists(fn)){
          data.default.f <- get_nc4_data_slice(fname=fn,scend="default") * FRESHMATTER[which(crops==cc)]
          data.fullharm.f <- get_nc4_data_slice(fname=fn,scend="fullharm") * FRESHMATTER[which(crops==cc)]
          data.harmnon.f <- get_nc4_data_slice(fname=fn,scend="harmnon") * FRESHMATTER[which(crops==cc)]
        } else{
          data.default.f <- data.fullharm.f <- data.harmnon.f <- array(NA,dim=c(31,208))
        }
        
        if(!all(is.na(c(data.default.f,data.fullharm.f,data.harmnon.f)))){
          i1 <- (which(ggcms==gg)-1)*3
          i2 <- which(aggs==agg)
          for(gadm in 1:dim(r.fao)[2]){
            hm0[gadm,i1+1,i2] <- mean.bias(r.fao[fao.season,gadm],data.default.f[clim.season.f,gadm])
            hm0[gadm,i1+2,i2] <- mean.bias(r.fao[fao.season,gadm],data.fullharm.f[clim.season.f,gadm])
            hm0[gadm,i1+3,i2] <- mean.bias(r.fao[fao.season,gadm],data.harmnon.f[clim.season.f,gadm])
          }
        }        
      }
    }#ggcms
  }#cl
  # select best aggregation mask (min mean.bias)
  hm <- apply(hm0,c(1,2),min)
  rownames(hm) <- gadm0.names
  colnames(hm) <- as.vector(outer(scens,ggcms, paste, sep=" "))
  # I want GGCM names come first, but inverting the outer function arguments does not do the right thing...
  cna <- as.vector(outer(scens,ggcms, paste, sep=" "))
  for(gg in 1:length(ggcms))
    for(ha in 1:length(scens))
      cna[(gg-1)*length(scens)+ha] <- paste(ggcms[gg],scens2[ha])
  colnames(hm) <- cna
  rm(cna)

  # remove all GGCM-scenario combinations that are non-existent
  hm <- hm[,colSums(is.na(hm))<nrow(hm)]
  # remove all countries that are not covered by any models
  hm <- hm[rowSums(!is.na(hm))>0,]
  best <- apply(hm,MARGIN=1,function(x) x[which.min(abs(x))])
  # add empty column to separate best from rest, this is painted by sepcol so it can be anything
  hm2 <- cbind(best,best,hm)
  colnames(hm2)[2] <- ""
  best.lab <- as.vector(apply(hm,MARGIN=1,function(x) paste(names(which.min(abs(x))),format(min(abs(x),na.rm=T),nsmall=2,digits=2))))

  # all countries ####
  png(paste(path.out,"heatmap_mean.bias_all_countries_best_agg_",cc,"_",prefix,".png",sep=""),
      width=6*300,height=10*300,res=300,pointsize=5,type="cairo")
  par(xpd=F,cex=1.5)
  ra <- max(abs(hm2),na.rm=T)
  heatmap.2b(hm2,Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(best.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm)[2],rowsep=0:dim(hm)[1],sepcol="white",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(17,17),col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-ra,ra,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="mean bias")
  dev.off()

  #top10 only ####
  # there should be a smarter way of ordering things according to the topproducer list...
  tp <- get(paste("topproducer",cc,sep="."))
  selec <- NULL
  for(tt in tp) selec <- c(selec,which(rownames(hm) %in% tt))
  best.lab <- as.vector(apply(hm[selec,],MARGIN=1,function(x) paste(names(which.min(abs(x))),format(min(abs(x),na.rm=T),nsmall=2,digits=2))))
  png(paste(path.out,"heatmap_mean.bias_top10_best_agg_",cc,"_",prefix,".png",sep=""),
      width=8*300,height=5*300,res=300,pointsize=9,type="cairo")
  par(xpd=F)
  ra <- max(abs(hm2[selec,]),na.rm=T)
  heatmap.2b(hm2[selec,],Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(best.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm2[selec,])[2],rowsep=0:dim(hm2[selec,])[1],sepcol="grey60",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(1.2,4),lhei=c(1.0,.3,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12),cex=0.9),
             margins=c(16,10),left.margin=0,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-ra,ra,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="mean bias")
  dev.off()
  save(hm0,selec,r.fao,tp,tt,gadm0.names,scens,scens2,ggcms,
       file=paste0(path.out,"heatmap_best_agg_mean_bias_",cc,".Rdata"))
  
}#crops
              