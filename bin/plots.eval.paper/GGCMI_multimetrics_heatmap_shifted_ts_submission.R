rm(list=ls(all=TRUE))


require(ncdf4)
#require(zoo)

#require(gplots)

options(warn=1)


# script to plot heatmap diagrams for GGCMI phase 1. 
# written by Christoph Mueller, PIK
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
sig.threshold <- 0.1
shift.threshold <- 0.2

aggs <- c("fixed_mirca_mask","dynamic_ray_mask","fixed_iizumi_mask","fixed_spam_mask")
dts <- c("none","lin","quad","ma","ffdtr")
scens <- c("default","fullharm","harmnon")
scens2 <- c("default","fullharm","harm-suffN")

# loop through some other models 
ggcms <- c("pdssat","epic-boku","epic-iiasa","gepic",
           "papsim","pegasus","lpj-guess","lpjml",
           "cgms-wofost","clm-crop","epic-tamu","orchidee-crop",
           "pepic","prysbi2")

# loop through crops
crops <- c("mai","whe","ric","soy")
cropsi <- c("maize_major","wheat","rice_major","soybean")
cropsl <- c("Maize","Wheat","Rice","Soy")
cropsl2 <- c("Maize","Wheat","Rice","Soybean")

#### functions ####

get_nc4_multimetrics <- function(fname,crd=cr0,dtd=dt0,mpd=mp0,scend=scen0,vid="tscorr"){
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
color.green <- colorRampPalette(c("lightblue","turquoise","green","green4"))
color.red2blue <- colorRampPalette(c("red4","red","orange","lightgoldenrod","lightblue","turquoise","blue","blue4"))


red.yellow.grey.blue.green <- function(n,alpha=1)
{
  if(n%%2 != 1) n <- n+1  
  nn <- n%/%2
  mid <- n%/%2 + 1
  col.veg <- color.red2blue(n)
  col.veg[mid] <- "grey90"
  col.veg[(mid+1):n] <- color.green(nn)
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

#### define settings ####

cr0 <- "none"#"mean-scale"#"variance-scale"#"none"
dt0 <- "ma" #"none"#"quad"
mp0 <- "false"
scen0 <- "default"
nm0 <- 1 #number of top model for the ensemble
wt0 <- "unweighted" #wheighting for ensemble (1:unweighted, 2:weighted)

# ignore one year at beginning and end of time series
ignore.y <- 1

# options to shift time series
shift.ts <- c(0,-1,1)


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

colo <- red.yellow.grey.blue.green

prefix <- paste("cr_",cr0,".dt_",dt0,".mp_",mp0,".multimetrics.processing_dt2",sep="")

for(cc in crops){
  hm0 <- array(NA,dim=c(length(gadm0.names),length(ggcms)*3,length(aggs)))
  hm02 <- array(NA,dim=c(length(gadm0.names),length(ggcms)*3,length(aggs),length(shift.ts)))
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
        r.fao <- get_nc4_ref_slice(fname,cc,dtd="none",mpd="true")[,which(gadm0.id %in% gadm0.id2)]
        
        fao.season <- c(max(rayf,cfirst[cl],faof):min(rayl,clast[cl],faol)) - faof+1
        clim.season.f <- c(max(rayf,faof,cfirst[cl]) :min(rayl,faol,clast[cl])) - cfirst[cl]+1
        # remove first and last years
        fao.season <- c(fao.season[1+ignore.y]:fao.season[length(fao.season)-ignore.y])
        clim.season.f <- c((clim.season.f[1+ignore.y]):(clim.season.f[length(clim.season.f)-ignore.y]))    

        fn <- paste(path,processed.ts,"/biascorr/gadm0/faostat/",agg,"/",gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                    if(cl==2 & (gg=="pegasus" | gg=="lpjml")) clast[cl]+1 else clast[cl],".biascorr.nc4",sep="")
        if(file.exists(fn)){
          data.default.f <- get_nc4_data_slice(fname=fn,dtd="none",mpd="true",scend="default")
          data.fullharm.f <- get_nc4_data_slice(fname=fn,dtd="none",mpd="true",scend="fullharm")
          data.harmnon.f <- get_nc4_data_slice(fname=fn,dtd="none",mpd="true",scend="harmnon")
        } else{
          data.default.f <- data.fullharm.f <- data.harmnon.f <- array(NA,dim=c(31,208))
        }

        fn <- paste(path,processed.ts,
                    "/multimetrics/gadm0/faostat/",agg,"/",
                    gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                    clast[cl],".multimetrics.nc4",sep="")
        if(!is.na(file.info(fn)$size)){
          i1 <- (which(ggcms==gg)-1)*3
          i2 <- which(aggs==agg)
          for(sc in scens){
            hm0[,i1+which(scens==sc),i2] <- get_nc4_multimetrics(fn,scend=sc)
          }
          for(gadm in 1:dim(r.fao)[2]){
            hm02[gadm,i1+1,i2,1] <- detrended_cor(r.fao[fao.season,gadm],data.default.f[clim.season.f,gadm])
            hm02[gadm,i1+2,i2,1] <- detrended_cor(r.fao[fao.season,gadm],data.fullharm.f[clim.season.f,gadm])
            hm02[gadm,i1+3,i2,1] <- detrended_cor(r.fao[fao.season,gadm],data.harmnon.f[clim.season.f,gadm])
            hm02[gadm,i1+1,i2,2] <- detrended_cor(c(NA,r.fao[fao.season,gadm]),c(data.default.f[clim.season.f,gadm],NA))
            hm02[gadm,i1+2,i2,2] <- detrended_cor(c(NA,r.fao[fao.season,gadm]),c(data.fullharm.f[clim.season.f,gadm],NA))
            hm02[gadm,i1+3,i2,2] <- detrended_cor(c(NA,r.fao[fao.season,gadm]),c(data.harmnon.f[clim.season.f,gadm],NA))
            hm02[gadm,i1+1,i2,3] <- detrended_cor(c(r.fao[fao.season,gadm],NA),c(NA,data.default.f[clim.season.f,gadm]))
            hm02[gadm,i1+2,i2,3] <- detrended_cor(c(r.fao[fao.season,gadm],NA),c(NA,data.fullharm.f[clim.season.f,gadm]))
            hm02[gadm,i1+3,i2,3] <- detrended_cor(c(r.fao[fao.season,gadm],NA),c(NA,data.harmnon.f[clim.season.f,gadm]))
          }
        }        
      }
    }#ggcms
  }#cl
  # select best correlation per shift and then best agg
  hm.cor <- apply(hm02,c(1,2,3),max,na.rm=T)
  hm.cor0 <- hm02[,,,1]
  hm.cor1 <- hm02[,,,2]
  hm.cor2 <- hm02[,,,3]
  hm.cor[!is.finite(hm.cor)] <- NA
  hm.cor.mask <- apply(hm02,c(1,2,3),function(x) which.max(x)[1])
  hm.cor.agg <- apply(hm.cor,c(1,2),max,na.rm=T)
  hm.cor.agg[!is.finite(hm.cor.agg)] <- NA
  hm.cor.agg0 <- apply(hm.cor0,c(1,2),max,na.rm=T)
  hm.cor.agg0[!is.finite(hm.cor.agg0)] <- NA
  hm.cor.agg1 <- apply(hm.cor1,c(1,2),max,na.rm=T)
  hm.cor.agg1[!is.finite(hm.cor.agg1)] <- NA
  hm.cor.agg2 <- apply(hm.cor2,c(1,2),max,na.rm=T)
  hm.cor.agg2[!is.finite(hm.cor.agg2)] <- NA
  hm.cor.agg.mask <- apply(hm.cor,c(1,2),function(x) which.max(x)[1])
  # if shifted is not at least "shift.threshold" better than non-shifted, stick to non-shifted
  sub1 <- hm.cor.agg
  sub1[is.na(sub1)] <- 0
  sub2 <- hm.cor.agg0
  sub2[is.na(sub2)] <- 0
  sub2 <- sub2 + shift.threshold
  hm.cor.agg[sub1<sub2] <- hm.cor.agg0[sub1<sub2]
  rm(sub1,sub2)
  # redo heatmap of which needed shifting
  hm.cor.agg.mask[] <- NA
  hm.cor.agg.mask[hm.cor.agg>(hm.cor.agg1-0.001) & hm.cor.agg<(hm.cor.agg1+0.001)] <- -1
  hm.cor.agg.mask[hm.cor.agg>(hm.cor.agg2-0.001) & hm.cor.agg<(hm.cor.agg2+0.001)] <- 1
  hm.cor.agg.mask[hm.cor.agg>(hm.cor.agg0-0.001) & hm.cor.agg<(hm.cor.agg0+0.001)] <- 0
  
  # I want GGCM names come first, but inverting the outer function arguments does not do the right thing...
  cna <- as.vector(outer(scens,ggcms, paste, sep=" "))
  for(gg in 1:length(ggcms))
    for(ha in 1:length(scens))
      cna[(gg-1)*length(scens)+ha] <- paste(ggcms[gg],scens2[ha])
  colnames(hm.cor.agg) <- cna
  rownames(hm.cor.agg) <- gadm0.names
  colnames(hm.cor.agg0) <- cna
  rownames(hm.cor.agg0) <- gadm0.names
  colnames(hm.cor.agg1) <- cna
  rownames(hm.cor.agg1) <- gadm0.names
  colnames(hm.cor.agg2) <- cna
  rownames(hm.cor.agg2) <- gadm0.names
  colnames(hm.cor.agg.mask) <- cna
  rownames(hm.cor.agg.mask) <- gadm0.names
  rm(cna)

  hm.cor.agg0 <- hm.cor.agg0[,colSums(is.na(hm.cor.agg))<nrow(hm.cor.agg)] # use multimetric-file based selection as the others have NA for non-significant as well
  hm.cor.agg1 <- hm.cor.agg1[,colSums(is.na(hm.cor.agg))<nrow(hm.cor.agg)] # use multimetric-file based selection as the others have NA for non-significant as well
  hm.cor.agg2 <- hm.cor.agg2[,colSums(is.na(hm.cor.agg))<nrow(hm.cor.agg)] # use multimetric-file based selection as the others have NA for non-significant as well
  hm.cor.agg.mask <- hm.cor.agg.mask[,colSums(is.na(hm.cor.agg))<nrow(hm.cor.agg)] # use multimetric-file based selection as the others have NA for non-significant as well
  hm.cor.agg <- hm.cor.agg[,colSums(is.na(hm.cor.agg))<nrow(hm.cor.agg)] # use multimetric-file based selection as the others have NA for non-significant as well

  # remove all countries with no data
  hm.cor.agg0 <- hm.cor.agg0[rowSums(!is.na(hm.cor.agg))>0,]
  hm.cor.agg1 <- hm.cor.agg1[rowSums(!is.na(hm.cor.agg))>0,]
  hm.cor.agg2 <- hm.cor.agg2[rowSums(!is.na(hm.cor.agg))>0,]
  hm.cor.agg.mask <- hm.cor.agg.mask[rowSums(!is.na(hm.cor.agg))>0,]
  hm.cor.agg <- hm.cor.agg[rowSums(!is.na(hm.cor.agg))>0,]

  best.cor.agg <- apply(hm.cor.agg,MARGIN=1,max,na.rm=T)
  hm.cor.agg.b <- cbind(best.cor.agg,best.cor.agg,hm.cor.agg)
  colnames(hm.cor.agg.b)[2] <- ""
  colnames(hm.cor.agg.b)[1] <- "best"
  best.cor.agg0 <- apply(hm.cor.agg0,MARGIN=1,max,na.rm=T)
  hm.cor.agg.b0 <- cbind(best.cor.agg0,best.cor.agg0,hm.cor.agg0)
  colnames(hm.cor.agg.b0)[2] <- ""
  colnames(hm.cor.agg.b0)[1] <- "best"
  best.cor.agg1 <- apply(hm.cor.agg1,MARGIN=1,max,na.rm=T)
  hm.cor.agg.b1 <- cbind(best.cor.agg1,best.cor.agg1,hm.cor.agg1)
  colnames(hm.cor.agg.b1)[2] <- ""
  colnames(hm.cor.agg.b1)[1] <- "best"
  best.cor.agg2 <- apply(hm.cor.agg2,MARGIN=1,max,na.rm=T)
  hm.cor.agg.b2 <- cbind(best.cor.agg2,best.cor.agg2,hm.cor.agg2)
  colnames(hm.cor.agg.b2)[2] <- ""
  colnames(hm.cor.agg.b2)[1] <- "best"
  
  best.cor.agg.lab <- as.vector(apply(hm.cor.agg,MARGIN=1,function(x) paste(names(which.max(x)),format(max(x,na.rm=T),nsmall=2,digits=2))))
  best.cor.agg.lab0 <- as.vector(apply(hm.cor.agg0,MARGIN=1,function(x) paste(names(which.max(x)),format(max(x,na.rm=T),nsmall=2,digits=2))))
  best.cor.agg.lab1 <- as.vector(apply(hm.cor.agg1,MARGIN=1,function(x) paste(names(which.max(x)),format(max(x,na.rm=T),nsmall=2,digits=2))))
  best.cor.agg.lab2 <- as.vector(apply(hm.cor.agg2,MARGIN=1,function(x) paste(names(which.max(x)),format(max(x,na.rm=T),nsmall=2,digits=2))))

  # all countries ####
  png(paste(path.out,"heatmap_all_countries_which_shift_",cc,"_",prefix,".png",sep=""),
      width=6*300,height=10*300,res=300,pointsize=5,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.mask,Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.cor.agg.lab[length(best.cor.agg.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.mask)[2],rowsep=0:dim(hm.cor.agg.mask)[1],sepcol="white",
             # set key to below heatmat, see
             #http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(16,22),left.margin=8,col=c(2,3,4),cexRow=1.5,cexCol=1.5,
             breaks=seq(-1.5,1.5,by=1),key.line=2.5,symkey=F,adjCol=c(NA,0.5),
             key.xlab="which offset provides best r: -1, 0, +1 years")
  dev.off()

  png(paste(path.out,"heatmap_all_countries_best_agg_shiftedTS_",cc,"_",prefix,".png",sep=""),
      width=6*300,height=10*300,res=300,pointsize=5,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.b,Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.cor.agg.lab[length(best.cor.agg.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.b)[2],rowsep=0:dim(hm.cor.agg.b)[1],sepcol="white",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(16,22),left.margin=8,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-1,1,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="timeseries correlation factor r")
  dev.off()
  
  png(paste(path.out,"heatmap_all_countries_best_agg_shift0_",cc,"_",prefix,".png",sep=""),
      width=6*300,height=10*300,res=300,pointsize=5,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.b0,Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.cor.agg.lab0[length(best.cor.agg.lab0):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.b0)[2],rowsep=0:dim(hm.cor.agg.b0)[1],sepcol="white",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(16,22),left.margin=8,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-1,1,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="timeseries correlation factor r")
  dev.off()

  png(paste(path.out,"heatmap_all_countries_best_agg_shift-1_",cc,"_",prefix,".png",sep=""),
      width=6*300,height=10*300,res=300,pointsize=5,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.b1,Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.cor.agg.lab1[length(best.cor.agg.lab1):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.b1)[2],rowsep=0:dim(hm.cor.agg.b1)[1],sepcol="white",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(16,22),left.margin=8,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-1,1,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="timeseries correlation factor r")
  dev.off()
  
  png(paste(path.out,"heatmap_all_countries_best_agg_shift+1_",cc,"_",prefix,".png",sep=""),
      width=6*300,height=10*300,res=300,pointsize=5,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.b2,Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.cor.agg.lab2[length(best.cor.agg.lab2):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.b2)[2],rowsep=0:dim(hm.cor.agg.b2)[1],sepcol="white",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(16,22),left.margin=8,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-1,1,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="timeseries correlation factor r")
  dev.off()
  
  #top10 only ####
  #selec <- which(rownames(hm) %in% get(paste("topproducer",cc,sep=".")))
  # there should be a smarter way of ordering things according to the topproducer list...
  tp <- get(paste("topproducer",cc,sep="."))
  selec <- NULL
  for(tt in tp) selec <- c(selec,which(rownames(hm.cor.agg.b) %in% tt))

  best.lab <- as.vector(apply(hm.cor.agg.b[selec,-c(1,2)],MARGIN=1,function(x) paste(names(which.max(x)),format(max(x,na.rm=T),nsmall=2,digits=2))))
  png(paste(path.out,"heatmap_top10_best_agg_shiftedTS_",cc,"_",prefix,".png",sep=""),
      width=8*300,height=5*300,res=300,pointsize=9,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.b[selec,],Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(best.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.b[selec,])[2],rowsep=0:dim(hm.cor.agg.b[selec,])[1],sepcol="grey60",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(1.2,4),lhei=c(1.0,.3,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12),cex=0.9),
             margins=c(16,10),left.margin=0,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-1,1,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="timeseries correlation factor r")
  dev.off()

  best.lab <- as.vector(apply(hm.cor.agg.b0[selec,-c(1,2)],MARGIN=1,function(x) paste(names(which.max(x)),format(max(x,na.rm=T),nsmall=2,digits=2))))
  png(paste(path.out,"heatmap_top10_best_agg_shift0_",cc,"_",prefix,".png",sep=""),
      width=8*300,height=5*300,res=300,pointsize=9,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.b0[selec,],Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(best.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.b0[selec,])[2],rowsep=0:dim(hm.cor.agg.b0[selec,])[1],sepcol="grey60",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(1.2,4),lhei=c(1.0,.3,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12),cex=0.9),
             margins=c(16,10),left.margin=0,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-1,1,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="timeseries correlation factor r")
  dev.off()
  
  best.lab <- as.vector(apply(hm.cor.agg.b1[selec,-c(1,2)],MARGIN=1,function(x) paste(names(which.max(x)),format(max(x,na.rm=T),nsmall=2,digits=2))))
  png(paste(path.out,"heatmap_top10_best_agg_shift-1_",cc,"_",prefix,".png",sep=""),
      width=8*300,height=5*300,res=300,pointsize=9,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.b1[selec,],Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(best.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.b1[selec,])[2],rowsep=0:dim(hm.cor.agg.b1[selec,])[1],sepcol="grey60",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(1.2,4),lhei=c(1.0,.3,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12),cex=0.9),
             margins=c(16,10),left.margin=0,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-1,1,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="timeseries correlation factor r")
  dev.off()
  
  best.lab <- as.vector(apply(hm.cor.agg.b2[selec,-c(1,2)],MARGIN=1,function(x) paste(names(which.max(x)),format(max(x,na.rm=T),nsmall=2,digits=2))))
  png(paste(path.out,"heatmap_top10_best_agg_shift+1_",cc,"_",prefix,".png",sep=""),
      width=8*300,height=5*300,res=300,pointsize=9,type="cairo")
  par(xpd=F,cex=1.5)
  heatmap.2b(hm.cor.agg.b2[selec,],Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(best.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=0:dim(hm.cor.agg.b2[selec,])[2],rowsep=0:dim(hm.cor.agg.b2[selec,])[1],sepcol="grey60",
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(1.2,4),lhei=c(1.0,.3,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12),cex=0.9),
             margins=c(16,10),left.margin=0,col=colo,cexRow=1.5,cexCol=1.5,
             breaks=seq(-1,1,length.out=102),key.line=2.5,symkey=T,adjCol=c(NA,0.5),
             key.xlab="timeseries correlation factor r")
  dev.off()
  
}#crops
              