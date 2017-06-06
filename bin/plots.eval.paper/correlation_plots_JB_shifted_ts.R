# (C) 2014-2017 Potsdam Institute for Climate Impact Research (PIK),
# written by Christoph Mueller & Jannes Breier, PIK
# cmueller@pik-potsdam.de
# Licensed under GNU AGPL Version 3 <LICENSE.txt in ggcmi directory>

##########################################################
# Script for plotting maps of correlation between        #
# reference data and model output for GGCMI phase 1      #
#                                                        #
##########################################################

rm(list=ls(all=TRUE))

##########################################################
# libraries ####
##########################################################
require(fields)
require(maps)
require(ncdf4)
require(zoo)

##########################################################
# global settings ####
##########################################################
NODATA <- 1e20
ncell <- 67420


path.sim <- "/p/projects/macmit/data/GGCMI/AgMIP.output/"
#path.sim <- "/iplex/01/2014/macmit/data/GGCMI/AgMIP.output/"
#path.sim <- "/iplex/01/2014/macmit/data/GGCMI/AgMIP.output/LPJmL/AgMERRA/maize/"

path.ref <- "/p/projects/macmit/data/GGCMI/reference/ray-dt/"
path.ref.detrend <- "/p/projects/macmit/data/GGCMI/reference/ray-dt/"
path.ray.2015 <- "/p/projects/macmit/users/cmueller/GGCMI/ray2015/"
picture.path <- "/p/projects/macmit/users/cmueller/GGCMI_cormap_test/"
path.iizumi <- "/p/projects/macmit/data/GGCMI/reference/iizumi/detrended/"


# scenario: model|climate|harmonisation|crop 

ggcms <- c("pDSSAT","EPIC-Boku","EPIC-IIASA","GEPIC",
           "pAPSIM","PEGASUS","LPJ-GUESS","LPJmL",
           "CGMS-WOFOST","CLM-Crop","EPIC-TAMU","ORCHIDEE-crop",
           "PEPIC","PRYSBI2")

clim <- c("AgMERRA","WFDEI.GPCC","WATCH")
harms <- c("default", "harmnon", "fullharm")

cropsl <- c("maize","wheat","rice","soy")
cropsl2 <- c("maize","wheat","rice","soybean")
cropsr <- c("A_Maize","C_Wheat","B_Rice","D_Soybean")
cropss <- c("mai","whe","ric","soy")

# time range
years <- seq(from=1961, to=2010, by=1)
start.year <- 1980
end.year <- 2008
r.years <- seq(from=1961, to=2008, by=1)
sim.years <- seq(from=1980, to=2010, by=1)
s.r <- which(r.years==start.year)
e.r <- which(r.years==end.year)
s.sim <- which(sim.years==start.year)
e.sim <- which(sim.years==end.year)

threshold.r2 <- 0.3
threshold.r <- sqrt(threshold.r2)
shift.threshold <- 0.2

##########################################################
# functions ####
##########################################################

read.map.from.nc <- function(fn,vn){  
  vec <- array(0,dim=c(ncell,2))
  ncs <- nc_open(fn)
  data <- ncvar_get(ncs,varid=vn)
  nc_close(ncs)
  data
}

read.ray.from.nc <- function(fn,vn){
  data <- read.map.from.nc(fn,vn)
  bins <- seq(0,by=0.15,length.out=7)
  data[data==1] <- 0
  for(i in 2:7){
    data[data==i] <- bins[i-1]+(bins[i]-bins[i-1])/2
  }
  data
}

get_nc4_ref_slice <- function(fname,crop,dtd="none",mpd="true"){ 
  nf <- nc_open(fname)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  # order mp|dt|time|fpu
  if(crop=="mai") { data <- ncvar_get(nf,varid="yield_mai")[mp,dt,,,]
  } else if(crop=="whe"){ data <- ncvar_get(nf,varid="yield_whe")[mp,dt,,,]
  } else if(crop=="ric"){ data <- ncvar_get(nf,varid="yield_ric")[mp,dt,,,]
  } else if(crop=="soy"){ data <- ncvar_get(nf,varid="yield_soy")[mp,dt,,,]}
  nc_close(nf)
  data
}

# detrend <- function(x1,L){
#   #
#   t_x1 <- rollapply(c(rep(NA,(L-1)),x1),L,mean,na.rm=T)
#   x1-t_x1
# }
# 
# detrended_cor <- function(x1,x2){
#   # check for computation of correlation coefficient | p-value
#   x12 <- x1*x2
#   if(length(x12[which(is.finite(x12))])<3 | sd(x1,na.rm=T)==0 | sd(x2,na.rm=T)==0)
#     return(c(NA,NA))
#   x1_dt <- detrend(x1,3)[-1]
#   x2_dt <- detrend(x2,3)[-1] 
#   x12 <- x1_dt*x2_dt
#   if(length(x12[which(is.finite(x12))])<4)
#     return(c(NA,NA))
#   co <- cor.test(x1_dt,x2_dt, use ="na.or.complete")
#   cop1 <- cor.test(c(x1_dt,NA),c(NA,x2_dt), use ="na.or.complete")
#   com1 <- cor.test(c(NA,x1_dt),c(x2_dt,NA), use ="na.or.complete")
#   if(!is.na(cop1$estimate) & cop1$estimate > (co$estimate + shift.threshold) & (is.na(com1$estimate) | cop1$estimate > com1$estimate)){
#     return(c(as.numeric(cop1$estimate),cop1$p.value))
#   }else if(!is.na(com1$estimate) & com1$estimate > (co$estimate + shift.threshold)){
#     return(c(as.numeric(com1$estimate),com1$p.value))
#   }
#   c(as.numeric(co$estimate),co$p.value)
# }

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

detrended_cor <- function(x1,x2){
  # check for computation of correlation coefficient | p-value
  x12 <- x1*x2
  if(length(x12[which(is.finite(x12))])<3 | sd(x1,na.rm=T)==0 | sd(x2,na.rm=T)==0)
    return(c(NA,NA))
  x1_dt <- detrend_ma(x1)
  x2_dt <- detrend_ma(x2)
  x12 <- x1_dt*x2_dt
  if(length(x12[which(is.finite(x12))])<4)
    return(c(NA,NA))
  co <- cor.test(x1_dt,x2_dt, use ="na.or.complete")
  cop1 <- cor.test(c(x1_dt,NA),c(NA,x2_dt), use ="na.or.complete")
  com1 <- cor.test(c(NA,x1_dt),c(x2_dt,NA), use ="na.or.complete")
  if(!is.na(cop1$estimate) & cop1$estimate > (co$estimate + shift.threshold) & (is.na(com1$estimate) | cop1$estimate > com1$estimate)){
    return(c(as.numeric(cop1$estimate),cop1$p.value))
  }else if(!is.na(com1$estimate) & com1$estimate > (co$estimate + shift.threshold)){
    return(c(as.numeric(com1$estimate),com1$p.value))
  }
  c(as.numeric(co$estimate),co$p.value)
}

rsme <- function(x1,x2){
  x12 <- x1*x2
  if(length(x12[which(is.finite(x12))])<1)
    return(NA)
  rsme <- sqrt(mean((x1-x2)^2,na.rm=TRUE))
}

mae <- function(x1,x2,w=NA){
  if(length(x1) != length(x2)) 
    stop("Invalid argument: 'sim' & 'obs' doesn't have the same length !")
  if(is.na(w)){
    mae <- mean(abs(x1 - x2), na.rm=T)  
  } else { # weighted mae
    if(length(x1) != length(w))
      stop("Invalid argument: sim % obs don't have the same length!")
    mae <- weighted.mean(abs(x1 - x2),w,na.rm=T)
  }
  return(mae)
}
# from the package simPopulatoin, which is no longer available via CRAN
quantileWt <- function(x, weights = NULL, 
                       probs = seq(0, 1, 0.25), na.rm = TRUE) {
  # initializations
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  x <- unname(x)  # unlike 'quantile', this never returns a named vector
  if(is.null(weights)) {
    return(quantile(x, probs, na.rm=na.rm, names=FALSE, type=1))
  } else if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
  else if(length(weights) != length(x)) {
    stop("'weights' must have the same length as 'x'")
  } else if(!all(is.finite(weights))) stop("missing or infinite weights")
  if(!is.numeric(probs) || all(is.na(probs)) || 
       isTRUE(any(probs < 0 | probs > 1))) {
    stop("'probs' must be a numeric vector with values in [0,1]")
  }
  if(length(x) == 0) return(rep.int(NA, length(probs)))
  if(!isTRUE(na.rm) && any(is.na(x))) {
    stop("missing values and NaN's not allowed if 'na.rm' is not TRUE")
  }
  # sort values and weights
  ord <- order(x, na.last=NA)
  x <- x[ord]
  weights <- weights[ord]
  # some preparations
  rw <- cumsum(weights)/sum(weights)
  # obtain quantiles
  select <- sapply(probs, function(p) min(which(rw >= p)))
  q <- x[select]
  return(q)
}

# from the package simPopulation, which is no longer available via CRAN
# modified to return the list element "groups"
spBwplotStats <- function(x, weights = NULL, coef = 1.5, 
                          zeros = TRUE, do.out = TRUE) {
  # initializations
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  if(!is.numeric(coef) || length(coef) != 1 || coef < 0) {
    stop("'coef' must be a single non-negative number")
  }
  # get quantiles
  if(isTRUE(zeros)) {
    zero <- ifelse(is.na(x), FALSE, x == 0)
    x <- x[!zero]
    if(is.null(weights)) nzero <- sum(zero)
    else {
      # if 'zeros' is not TRUE, these checks are done in 'quantileWt'
      # but here we need them since we use subscripting
      if(!is.numeric(weights)) stop("'weights' must be a numeric vector")
      else if(length(weights) != length(zero)) {
        stop("'weights' must have the same length as 'x'")
      }
      # avoid infinite elements in weights
      weights[!is.finite(weights)] <- 0
      nzero <- sum(weights[zero])
      weights <- weights[!zero]
    }
  } else nzero <- NULL
  # normalize weights???
  ok <- !is.na(x)
  n <- if(is.null(weights)) sum(ok) else sum(weights[ok])
  if(n == 0) stats <- rep.int(NA, 5)
  else stats <- quantileWt(x, weights)
  iqr <- diff(stats[c(2, 4)])  # inter quartile range
  if(coef == 0) do.out <- FALSE
  else {
    if(is.na(iqr)) out <- is.infinite(x) 
    else {
      lower <- stats[2] - coef * iqr
      upper <- stats[4] + coef * iqr
      out <- ifelse(ok, x < lower | x > upper, FALSE)
    }
    if(any(out)) stats[c(1, 5)] <- range(x[!out], na.rm=TRUE)
  }
  #stats needs to be a matrix with 5 rows and 1 column to work with bxp()
  res <- list(stats=matrix(stats,5,1), n=n, nzero=nzero, 
              out=if(isTRUE(do.out)) x[out] else numeric(),
              groups=if(isTRUE(do.out) & any(out)) rep(1,length(out)) else numeric())
  class(res) <- "spBwplotStats"
  res
}


# image.legend <- function(data,zlim=(c(min(data,na.rm=T),max(data,na.rm=T)+1)),
#                          col=tim.colors(zlim[2]-zlim[1]+1),
#                          breaks=seq(min(data,na.rm=T),max(data,na.rm=T)+1,,length(col)+1),
#                          at=c(1:length(col))+min(data,na.rm=T)-0.5,
#                          labels=c(min(data,na.rm=T):max(data,na.rm=T)))
# {
#   if(length(col)!=length(breaks)-1){
#     cat("missmatching color/break scheme! col:",length(col),"break:", length(breaks),"\n")
#     return()
#   }
#   # for a correct representation of all colors, zlim needs to be max(data)+1
#   image.plot(data,zlim=zlim,legend.only=T,
#              breaks=breaks, #one more break than colors
#              col=col,
#              # default: one label per color
#              axis.args=list(at=at,labels=labels) 
#   )
# }
#########################################################
# colors #####
#########################################################

noco<- 60
hnoco <- noco/2
# color scheme
col.diff <- rainbow(noco)
r1 <- rainbow(hnoco,start=6/6,end=1/6)
col.diff[1:hnoco]<-r1
r1 <- rainbow(hnoco,start=2/6,end=4/6)
col.diff[(hnoco+1):(2*hnoco)]<-r1
r1 <- colorRampPalette(c("red4","red","orange","lightgoldenrod","yellow","green","green4"))
col.r2 <- r1(20) # 20 color steps
col.r2[21] <- "gray95" #21st color grey for non-significance
col.mods <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#ffed6f','gray95')
col.box <- c('#1b9e77','#d95f02','#7570b3')

##########################################################
# main #######
##########################################################

#   for (cl in clim){
for (cl in clim[1]){
  #for(cc in 1){ 
  for(cc in 1:length(cropsl)){
    # plot correlation of Ray vs. Iizumi ####
    nf <- nc_open(paste0(path.iizumi,"iizumi.2013JAN29.",cropsl2[cc],".1982-2006.30min.nc4"))
    iizumiref <- ncvar_get(nf,varid="yield50")[1,4,,,]
    nc_close(nf)
    nf <- nc_open(paste0(path.ref,cropss[cc],"_weight_ray_1961-2008.nc4"))
    rayref <- ncvar_get(nf,varid=paste0("yield_",cropss[cc]))[1,4,,,22:46]
    nc_close(nf)
    mapi <- array(NA,dim=c(720,360,2))
    for(i in 1:720) {
      for(j in 1:360){
        mapi[i,j,] <- detrended_cor(rayref[i,j,], iizumiref[i,j,])
      }
    }
    # coefficient of determination
    mapr2 <- mapi[,,1] * mapi[,,1]
    # 0.9 confidence interval - non signficance in grey (col.r2)
    mapr2[mapi[,,2]>0.1] <- 1.05
    save(mapr2,file=paste0(picture.path,cropsl[cc],"_",tolower(cl),
                                   "_hist_yield_",cropss[cc],"_ray_vs_iizumi_R2__shifted_ts_processed_dt.Rdata"))
    # png(paste0(picture.path,cropsl[cc],"_",tolower(cl),
    #            "_hist_yield_",cropss[cc],"_ray_vs_iizumi_R2__shifted_ts_processed_dt.png"),
    #     height=5*300,width=8*300,res=300,pointsize=6,type="cairo") 
    # par("mar"=c(2,2,2,5))
    # image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
    #       mapr2[,360:1],main="",col=col.r2,zlim=c(0,1.05))
    # grid(10,10,lwd=.5)
    # map(boundary=T,add=T)
    # image.plot(mapr2,legend.only=T,col=col.r2[1:20],zlim=c(0,1),oldstyle=TRUE)
    # #close.screen(all=T)
    # dev.off()
    
    # plot Ray R2 values ####
    rayr2 <- read.map.from.nc(paste0(path.ray.2015,"NonCategoricalFigure2",cropsr[cc],"05.nc"),"Data")
    # png(paste0(picture.path,cropsl[cc],"_",tolower(cl),
    #            "_hist_yield_",cropss[cc],"_ray_R2_shifted_ts_processed_dt.png"),
    #     height=5*300,width=8*300,res=300,pointsize=6,type="cairo") 
    # par("mar"=c(2,2,2,5))
    # image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
    #       rayr2,main="",col=col.r2,zlim=c(0,1.05))
    # grid(10,10,lwd=.5)
    # map(boundary=T,add=T)
    # image.plot(rayr2,legend.only=T,col=col.r2[1:20],zlim=c(0,1),oldstyle=TRUE)
    # #close.screen(all=T)
    # dev.off()
    rm(rayr2,mapi,mapr2,iizumiref,rayref)
    gc()
    
    tn <- paste0(picture.path,cropsl[cc],"_statistics_shifted_ts.csv")
    write("Crop,GGCM,harmonization,rsme,median R2,median sig R2,mean R2,mean sig R2, area weighted mean R2,area weighted mean sig R2,production weighted mean R2,production weighted mean sig R2,mae,# pixels with significant correlation,# pixels with R2>=0.3,area with significant correlation,share of area with significant correlation, area with R2>=0.3, share of area with R2>=0.3,production with significant correlation,share of production with significant correlation, production with R2>=0.3,share of production with R2>=0.3, # pixels with best correlation,area with best correlation,share of area with best correlation, production with best correlation, share of production with best correlation, # pixels with best rsme,area with best rsme,  share of area with best rsme, production with best rsme, share of production with best rsme",
          file=tn,sep="")
    best.cor <- best.cor2 <- best.mod <- best.mod2 <- num30 <- totaln <- totalv <- mean <- map.rsme <- best.rsme <- best.mod.rsme <- map.mae <- best.mae <- best.mod.mae <- array(NA,dim=c(720,360))
    
    for(ha in harms){ # loop through vector for default, harmnon and fullharm
    #for(ha in harms[1]){
      # maps for best correlation and best model per harmonization
      best.cor[] <- best.cor2[] <- best.mod[] <- best.mod2[] <- 
        num30[] <- totaln[] <- totalv[] <- mean[] <- map.rsme[] <- best.rsme[] <- 
        best.mod.rsme[] <- map.mae[] <- best.mae[] <- best.mod.mae[] <- NA
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_shifted_ts.Rdata")
      sn2 <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      list.stats <- list() # empty list for all GGCM-specific statistics
      list.summarydata <- list() # empty list for boxplots 
      list.data <- list()
      #if(file.exist(sn)){
      #  load(sn)
      #} else {
        for (gg in ggcms){
          #for (gg in ggcms[c(1,8)]){
          # initializing
          map.rsme[] <- NA
          #mean[] <- NA
          #totalv[] <- NA
          #totaln[] <- NA
          #num30[] <- NA
          map.mae[] <- NA
          
          ## model output
          # condition irrigated  
          fn <- paste0(path.sim,gg,"/",cl,"/",cropsl[cc],"/", tolower(gg),"_", tolower(cl),
                       "_hist_",ha, "_firr_yield_",cropss[cc], "_annual_1980_2010.nc4")
          if(!file.exists(fn)){
            # fill NA in statistics table, skip
            list.stats[[which(ggcms==gg)]] <- rep(NA,19)
            list.summarydata[[which(ggcms==gg)]] <- list()
            list.data[[which(ggcms==gg)]] <- NA
            if(gg==ggcms[length(ggcms)]){
              rsquared <- rsquared.sig <- best.cor*best.cor
              rsquared[best.mod<0] <- 0 # include non-significant correlations as zero
              rsquared.sig[best.mod<0] <- NA # include non-significant correlations as zero
              rsquared[ray.mean>0 & !is.finite(rsquared)] <- 0
              #best.mod[i,j] <- -9
              issig <- best.mod
              issig[best.mod<0] <- 0
              issig[best.mod>0] <- 1
              det.r_ens <- best.cor*best.cor
              det.r_ens[best.mod<0] <- 1.05
              list.stats[[length(ggcms)+1]] <- c(NA,median(rsquared,na.rm=T),median(rsquared.sig,na.rm=T),
                                                 mean(rsquared,na.rm=T),mean(rsquared.sig,na.rm=T),
                                                 weighted.mean(rsquared,area.irrf,na.rm=T),weighted.mean(rsquared.sig,area.irrf,na.rm=T),
                                                 weighted.mean(rsquared,production,na.rm=T),weighted.mean(rsquared.sig,production,na.rm=T),
                                                 NA,#mae(r.ray.sub,sim.yield.sub),
                                                 length(which(issig==1)),length(which(rsquared>=0.3)),
                                                 sum(area.irrf*issig,na.rm=T),sum(area.irrf*issig,na.rm=T)/sum(area.irrf,na.rm=T),
                                                 sum(area.irrf[rsquared>=0.3],na.rm=T),sum(area.irrf[rsquared>=0.3],na.rm=T)/sum(area.irrf,na.rm=T),
                                                 sum(production*issig,na.rm=T),sum(production*issig,na.rm=T)/sum(production,na.rm=T),
                                                 sum(production[rsquared>=0.3],na.rm=T),sum(production[rsquared>=0.3],na.rm=T)/sum(production,na.rm=T)
              )
              
              list.data[[length(ggcms)+1]] <- list(det.r_ens,NA,NA,NA,rsquared,rsquared.sig,
                                                   best.rsme,NA,NA,best.cor,best.mae)
              rm(det.r_ens)
              
              png(paste0(picture.path,"dummy_shifted_ts_processed_dt.png"),type="cairo")
              b1 <- boxplot(as.vector(best.rsme),plot=F)
              #rsqua <- det.r_sim
              #rsqua[cor.r_sim[,,2]>0.1] <- 0
              b2 <- boxplot(as.vector(rsquared.sig),plot=F)
              dev.off()
              b3 <- spBwplotStats(as.vector(best.mae))
              b4 <- spBwplotStats(as.vector(best.mae),as.vector(area.irrf))
              b5 <- spBwplotStats(as.vector(best.mae),as.vector(production))
              b6 <- spBwplotStats(as.vector(rsquared),as.vector(area.irrf))
              b7 <- spBwplotStats(as.vector(rsquared),as.vector(production))
              b8 <- spBwplotStats(as.vector(rsquared.sig),as.vector(area.irrf))
              b9 <- spBwplotStats(as.vector(rsquared.sig),as.vector(production))
              list.summarydata[[length(ggcms)+1]] <- list(rsme.stats=list(stats=b1$stats,n=b1$n),
                                                          r2.stats=list(stats=b2$stats,n=b2$n,out=b2$out,group=b2$group),
                                                          mae.stats=list(stats=b3$stats,n=b3$n,out=b3$out,group=b3$group),
                                                          mae_area.stats=list(stats=b4$stats,n=b4$n,out=b4$out,group=b4$group),
                                                          mae_production.stats=list(stats=b5$stats,n=b5$n,out=b5$out,group=b5$group),
                                                          r2_area.stats=list(stats=b6$stats,n=b6$n,out=b6$out,group=b6$group),
                                                          r2_production.stats=list(stats=b7$stats,n=b7$n,out=b7$out,group=b7$group),
                                                          r2_sig_area.stats=list(stats=b8$stats,n=b8$n,out=b8$out,group=b8$group),
                                                          r2_sig_production.stats=list(stats=b9$stats,n=b9$n,out=b9$out,group=b9$group))
              
            }
            next
          }
          vn <- paste0("yield_",cropss[cc])
          sim.yield.ir <- read.map.from.nc(fn,vn)
          
          # condition rainfed
          if(gg != "PRYSBI2") fn <- paste0(path.sim,gg,"/",cl,"/",cropsl[cc],"/", tolower(gg),"_", tolower(cl),
                                           "_hist_",ha, "_noirr_yield_",cropss[cc], "_annual_1980_2010.nc4")
          vn <- paste0("yield_",cropss[cc])
          sim.yield.rf <- read.map.from.nc(fn,vn)
          
          # weighted average for consistency between model output and reference data
          fn <- paste0(path.sim,"processed/masks/weight/",cropsl[cc], ".nc4")
          area.rf <- read.map.from.nc(fn,"rainfed")
          area.ir <- read.map.from.nc(fn,"irrigated")
          sim.yield.irif <- array(NA,dim=c(720,360,31))
          
          for (i in 1:dim(sim.yield.ir)[3]){
            sim.yield.irif[,,i] <- (sim.yield.ir[,,i] * area.ir + sim.yield.rf[,,i] * area.rf)/(area.ir + area.rf)
          }
          
          ## reference data
          fname <- paste0(path.ref.detrend,cropss[cc],"_weight_ray_1961-2008.nc4")
          crop <- cropss[cc]
          r.ray <- get_nc4_ref_slice(fname,crop)
          
          ## detrending and correlation computation (using function detrended_cor)
          r.ray.sub <- r.ray[,,s.r:e.r]
          sim.yield.sub <- sim.yield.irif[,,s.sim:e.sim]
          cor.r_sim <- array(NA,dim=c(720,360,2))
          
          ## area
          area.irrf <- area.ir + area.rf
          ## production
          ray.mean <- apply(r.ray.sub,c(1,2),mean,na.rm=T)
          ray.mean[!is.finite(ray.mean)] <- NA
          production <- area.irrf * ray.mean
          production[!is.finite(production)] <- NA
          
          for (i in 1:dim(cor.r_sim)[1]){
            for(j in 1:dim(cor.r_sim)[2]){
              cor.r_sim[i,j,] <- detrended_cor(r.ray.sub[i,j,], sim.yield.sub[i,j,])
              map.rsme[i,j] <- rsme(r.ray.sub[i,j,],sim.yield.sub[i,j,])
              map.mae[i,j] <- mae(r.ray.sub[i,j,],sim.yield.sub[i,j,])
              if(!all(is.na(r.ray.sub[i,j,]))){
                totaln[i,j] <- if(is.na(totaln[i,j])) 1 else totaln[i,j] + 1
              }
              if(is.na(best.rsme[i,j]) | (!is.na(map.rsme[i,j]) & best.rsme[i,j]>map.rsme[i,j])){
                best.rsme[i,j] <- map.rsme[i,j]
                best.mod.rsme[i,j] <- which(ggcms==gg)
              }
              if(is.na(best.mae[i,j]) | (!is.na(map.mae[i,j]) & best.mae[i,j]>map.mae[i,j])){
                best.mae[i,j] <- map.mae[i,j]
                best.mod.mae[i,j] <- which(ggcms==gg)
              }
              if(!is.na(cor.r_sim[i,j,2])){
                if(cor.r_sim[i,j,2]<=0.1){
                  totalv[i,j] <- if(is.na(totalv[i,j])) 1 else totalv[i,j] + 1
                  if(is.na(best.cor[i,j]) | best.cor[i,j]<cor.r_sim[i,j,1]){
                    best.cor2[i,j] <- best.cor[i,j]
                    best.mod2[i,j] <- best.mod[i,j]
                    best.cor[i,j] <- cor.r_sim[i,j,1]
                    best.mod[i,j] <- which(ggcms==gg)
                  }  
                  if(cor.r_sim[i,j,1]>threshold.r) {
                    num30[i,j] <- if(is.na(num30[i,j])) 1 else num30[i,j] + 1
                  }
                  mean[i,j] <- if(is.na(mean[i,j])) cor.r_sim[i,j,1] else mean[i,j] + cor.r_sim[i,j,1]
                } else if(is.na(best.mod[i,j])){
                  best.mod[i,j] <- -9
                }
                
              }
            }  
          }
          
          # coefficient of determination
          det.r_sim <- cor.r_sim[,,1] * cor.r_sim[,,1]
          # 0.9 confidence interval - non signficance in grey (col.r2)
          det.r_sim[cor.r_sim[,,2]>0.1] <- 1.05
          
          #  Crop,GGCM,
          # mean R2,area weighted mean R2,production weighted mean R2,
          # # pixels with significant correlation,# pixels with R2>=0.3,
          # area with significant correlation,area with R2>=0.3,
          # production with significant correlation,production with R2>=0.3,
          # # pixels with best correlation,area with best correlation,production with best correlation
          
          rsquared <- rsquared.sig <- cor.r_sim[,,1] * cor.r_sim[,,1]
          rsquared[cor.r_sim[,,2]>0.1] <- 0 # include non-significant correlations as zero
          rsquared.sig[cor.r_sim[,,2]>0.1] <- NA # include non-significant correlations as zero
          rsquared[ray.mean>0 & !is.finite(rsquared)] <- 0
          issig <- cor.r_sim[,,2]
          issig[cor.r_sim[,,2]>0.1] <- 0
          issig[cor.r_sim[,,2]>0] <- 1
          
        
          list.stats[[which(ggcms==gg)]] <- c(rsme(r.ray.sub,sim.yield.sub),median(rsquared,na.rm=T),median(rsquared.sig,na.rm=T),
                                              mean(rsquared,na.rm=T),mean(rsquared.sig,na.rm=T),
                                              weighted.mean(rsquared,area.irrf,na.rm=T),weighted.mean(rsquared.sig,area.irrf,na.rm=T),
                                              weighted.mean(rsquared,production,na.rm=T),weighted.mean(rsquared.sig,production,na.rm=T),
                                              mae(r.ray.sub,sim.yield.sub),
                                              length(which(issig==1)),length(which(rsquared>=0.3)),
                                              sum(area.irrf*issig,na.rm=T),sum(area.irrf*issig,na.rm=T)/sum(area.irrf,na.rm=T),
                                              sum(area.irrf[rsquared>=0.3],na.rm=T),sum(area.irrf[rsquared>=0.3],na.rm=T)/sum(area.irrf,na.rm=T),
                                              sum(production*issig,na.rm=T),sum(production*issig,na.rm=T)/sum(production,na.rm=T),
                                              sum(production[rsquared>=0.3],na.rm=T),sum(production[rsquared>=0.3],na.rm=T)/sum(production,na.rm=T)
          )
          
          list.data[[which(ggcms==gg)]] <- list(det.r_sim,sim.yield.sub,totaln,totalv,rsquared,rsquared.sig,
                                                map.rsme,num30,mean,issig,map.mae)
          png(paste0(picture.path,"dummy_shifted_ts_processed_dt.png"),type="cairo")
          b1 <- boxplot(as.vector(map.rsme),plot=F)
          #rsqua <- det.r_sim
          #rsqua[cor.r_sim[,,2]>0.1] <- 0
          b2 <- boxplot(as.vector(rsquared.sig),plot=F)
          dev.off()
          b3 <- spBwplotStats(as.vector(map.mae))
          b4 <- spBwplotStats(as.vector(map.mae),as.vector(area.irrf))
          b5 <- spBwplotStats(as.vector(map.mae),as.vector(production))
          b6 <- spBwplotStats(as.vector(rsquared),as.vector(area.irrf))
          b7 <- spBwplotStats(as.vector(rsquared),as.vector(production))
          b8 <- spBwplotStats(as.vector(rsquared.sig),as.vector(area.irrf))
          b9 <- spBwplotStats(as.vector(rsquared.sig),as.vector(production))
          list.summarydata[[which(ggcms==gg)]] <- list(rsme.stats=list(stats=b1$stats,n=b1$n),
                                                       r2.stats=list(stats=b2$stats,n=b2$n,out=b2$out,group=b2$group),
                                                       mae.stats=list(stats=b3$stats,n=b3$n,out=b3$out,group=b3$group),
                                                       mae_area.stats=list(stats=b4$stats,n=b4$n,out=b4$out,group=b4$group),
                                                       mae_production.stats=list(stats=b5$stats,n=b5$n,out=b5$out,group=b5$group),
                                                       r2_area.stats=list(stats=b6$stats,n=b6$n,out=b6$out,group=b6$group),
                                                       r2_production.stats=list(stats=b7$stats,n=b7$n,out=b7$out,group=b7$group),
                                                       r2_sig_area.stats=list(stats=b8$stats,n=b8$n,out=b8$out,group=b8$group),
                                                       r2_sig_production.stats=list(stats=b9$stats,n=b9$n,out=b9$out,group=b9$group))
          if(gg==ggcms[length(ggcms)]){ # do for best model ensemble
            rsquared <- rsquared.sig <- best.cor*best.cor
            rsquared[best.mod<0] <- 0 # include non-significant correlations as zero
            rsquared.sig[best.mod<0] <- NA # include non-significant correlations as zero
            rsquared[ray.mean>0 & !is.finite(rsquared)] <- 0
            #best.mod[i,j] <- -9
            issig <- best.mod
            issig[best.mod<0] <- 0
            issig[best.mod>0] <- 1
            det.r_ens <- best.cor*best.cor
            det.r_ens[best.mod<0] <- 1.05
            list.stats[[length(ggcms)+1]] <- c(NA,median(rsquared,na.rm=T),median(rsquared.sig,na.rm=T),
                                                mean(rsquared,na.rm=T),mean(rsquared.sig,na.rm=T),
                                                weighted.mean(rsquared,area.irrf,na.rm=T),weighted.mean(rsquared.sig,area.irrf,na.rm=T),
                                                weighted.mean(rsquared,production,na.rm=T),weighted.mean(rsquared.sig,production,na.rm=T),
                                                NA,#mae(r.ray.sub,sim.yield.sub),
                                                length(which(issig==1)),length(which(rsquared>=0.3)),
                                                sum(area.irrf*issig,na.rm=T),sum(area.irrf*issig,na.rm=T)/sum(area.irrf,na.rm=T),
                                                sum(area.irrf[rsquared>=0.3],na.rm=T),sum(area.irrf[rsquared>=0.3],na.rm=T)/sum(area.irrf,na.rm=T),
                                                sum(production*issig,na.rm=T),sum(production*issig,na.rm=T)/sum(production,na.rm=T),
                                                sum(production[rsquared>=0.3],na.rm=T),sum(production[rsquared>=0.3],na.rm=T)/sum(production,na.rm=T)
            )
            
            list.data[[length(ggcms)+1]] <- list(det.r_ens,NA,NA,NA,rsquared,rsquared.sig,
                                                  best.rsme,NA,NA,best.cor,best.mae)
            rm(det.r_ens)
            
            png(paste0(picture.path,"dummy_shifted_ts_processed_dt.png"),type="cairo")
            b1 <- boxplot(as.vector(best.rsme),plot=F)
            #rsqua <- det.r_sim
            #rsqua[cor.r_sim[,,2]>0.1] <- 0
            b2 <- boxplot(as.vector(rsquared.sig),plot=F)
            dev.off()
            b3 <- spBwplotStats(as.vector(best.mae))
            b4 <- spBwplotStats(as.vector(best.mae),as.vector(area.irrf))
            b5 <- spBwplotStats(as.vector(best.mae),as.vector(production))
            b6 <- spBwplotStats(as.vector(rsquared),as.vector(area.irrf))
            b7 <- spBwplotStats(as.vector(rsquared),as.vector(production))
            b8 <- spBwplotStats(as.vector(rsquared.sig),as.vector(area.irrf))
            b9 <- spBwplotStats(as.vector(rsquared.sig),as.vector(production))
            list.summarydata[[length(ggcms)+1]] <- list(rsme.stats=list(stats=b1$stats,n=b1$n),
                                                         r2.stats=list(stats=b2$stats,n=b2$n,out=b2$out,group=b2$group),
                                                         mae.stats=list(stats=b3$stats,n=b3$n,out=b3$out,group=b3$group),
                                                         mae_area.stats=list(stats=b4$stats,n=b4$n,out=b4$out,group=b4$group),
                                                         mae_production.stats=list(stats=b5$stats,n=b5$n,out=b5$out,group=b5$group),
                                                         r2_area.stats=list(stats=b6$stats,n=b6$n,out=b6$out,group=b6$group),
                                                         r2_production.stats=list(stats=b7$stats,n=b7$n,out=b7$out,group=b7$group),
                                                         r2_sig_area.stats=list(stats=b8$stats,n=b8$n,out=b8$out,group=b8$group),
                                                         r2_sig_production.stats=list(stats=b9$stats,n=b9$n,out=b9$out,group=b9$group))
            
          }
          rm(b1,b2,b3,b4,b5,b6,b7,b8,b9)
          gc()
          ## visualisation
          # coefficient of determination
          cat(gg,cropsl[cc],"R2:",range(det.r_sim,na.rm=T),range(cor.r_sim[,,1] * cor.r_sim[,,1],na.rm=T),"\n")
          save(det.r_sim,file=paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
                                             "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_R2_shifted_ts_processed_dt.Rdata"))
          png(paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
                     "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_R2_shifted_ts_processed_dt.png"),
              height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
          par("mar"=c(2,2,2,5))
          image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
                det.r_sim[,360:1],main="",col=col.r2,zlim=c(0,1.05))
          grid(10,10,lwd=.5)
          map(boundary=T,add=T)
          cat("det.r_sim",range(det.r_sim,na.rm=T),"\n")
          image.plot(det.r_sim,legend.only=T,col=col.r2[1:20],zlim=c(0,1),oldstyle=TRUE)
          #close.screen(all=T)
          dev.off()
          
          # rsme
          mapi <- map.rsme
          mapi[mapi>10] <- 10
          save(mapi,file=paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
                                        "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_rsme_shifted_ts_processed_dt.Rdata"))
          png(paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
                     "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_rsme_shifted_ts_processed_dt.png"),
              height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
          par("mar"=c(2,2,2,5))
          image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
                mapi[,360:1],main="",col=col.r2[20:1])
          grid(10,10,lwd=.5)
          map(boundary=T,add=T)
          cat("mapi",range(mapi,na.rm=T),"\n")
          image.plot(mapi,legend.only=T,col=col.r2[20:1],oldstyle=TRUE,
                     axis.args=list(at=c(0:10),labels=c(0:9,">=10")))
          #close.screen(all=T)
          dev.off()

          # mae
          mapi <- map.mae
          mapi[mapi>10] <- 10
          save(mapi,file=paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
                                        "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_mae_shifted_ts_processed_dt.Rdata"))
          png(paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
                     "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_mae_shifted_ts_processed_dt.png"),
              height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
          par("mar"=c(2,2,2,5))
          image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
                mapi[,360:1],main="",col=col.r2[20:1])
          grid(10,10,lwd=.5)
          map(boundary=T,add=T)
          cat("mapi",range(mapi,na.rm=T),"\n")
          image.plot(mapi,legend.only=T,col=col.r2[20:1],oldstyle=TRUE,
                     axis.args=list(at=c(0:10),labels=c(0:9,">=10")))
          #close.screen(all=T)
          dev.off()
          rm(mapi)
          gc()
          
          # yields
          #           sim.mean <- apply(sim.yield.sub,c(1,2),mean,na.rm=T)
          #           ra <- range(sim.mean,ray.mean,na.rm=T)
          #           diff <- sim.mean-ray.mean
          #           rad <- max(abs(diff),na.rm=T)
          #           rad <- c(-rad,rad)
          #           png(paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
          #                      "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_yielddiff_shifted_ts_processed_dt.png"),
          #               height=10*300,width=8*300,res=300,pointsize=6,type="cairo") 
          #           split.screen(c(3,1))
          #           screen(1)
          #           par("mar"=c(2,2,2,5))
          #           image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
          #                 ray.mean[,360:1],col=col.r2[20:1],main="ref",zlim=ra)
          #           grid(10,10,lwd=.5)
          #           map(boundary=T,add=T)
          #           image.plot(ray.mean,legend.only=T,col=col.r2[20:1],oldstyle=TRUE,zlim=ra)
          #           screen(2)
          #           par("mar"=c(2,2,2,5))
          #           image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
          #                 sim.mean[,360:1],col=col.r2[20:1],main="sim",zlim=ra)
          #           grid(10,10,lwd=.5)
          #           map(boundary=T,add=T)
          #           image.plot(ray.mean,legend.only=T,col=col.r2[20:1],oldstyle=TRUE,zlim=ra)
          #           screen(3)
          #           par("mar"=c(2,2,2,5))
          #           image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
          #                 diff[,360:1],col=col.diff,main="sim",zlim=rad)
          #           grid(10,10,lwd=.5)
          #           map(boundary=T,add=T)
          #           cat("diff",range(diff,na.rm=T),"\n")
          #           image.plot(diff,legend.only=T,col=col.diff,zlim=rad,oldstyle=TRUE)
          #           close.screen(all=T)
          #           dev.off()
          #           
          # p-value
          #       png(paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
          #                  "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.p_shifted_ts_processed_dt.png"),
          #           height=5*300,width=5*300,res=300,pointsize=6) 
          #       par("mar"=c(2,2,2,5))
          #       image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
          #             cor.r_sim[,360:1,2] ,main="", col=col.r2,zlim=c(0,1.05),oldstyle=TRUE)
          #       grid(10,10,lwd=.5)
          #       map(boundary=T,add=T)
          #       image.plot(cor.r_sim,legend.only=T,col=col.r2[1:20],zlim=c(0,1),oldstyle=TRUE)
          #       close.screen(all=T)
          #       dev.off()
          
          #           # correlation coefficient
          #           png(paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
          #                      "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          #               height=5*300,width=8*300,res=300,pointsize=6,type="cairo") 
          #           par("mar"=c(2,2,2,5))
          #           image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
          #                 cor.r_sim[,360:1,1] ,main="",col=col.diff,zlim=c(-1,1),oldstyle=TRUE)
          #           grid(10,10,lwd=.5)
          #           map(boundary=T,add=T)
          #           image.plot(cor.r_sim,legend.only=T,col=col.diff,zlim=c(-1,1),oldstyle=TRUE)
          #           #close.screen(all=T)
          #           dev.off()
          
        }
      save(det.r_sim,cor.r_sim,best.cor,best.cor2,best.mod,best.mod2,num30,map.rsme,
           totalv,totaln,r.ray.sub,sim.yield.sub,mean,list.stats,list.summarydata,map.mae,
           file=sn)        
      save(list.summarydata,file=sn2)        
      #}#else file.exist
      # plotting best correlation and best model
      
      # fill table with statistics
      for(gg in ggcms){
        # add "best of" metrics
        # # pixels with best correlation,area with best correlation,production with best correlation
        #    write.table(annual.az,paste(results.folder,"/LP",suffix,".txt",sep=""),append=T,row.names=F,col.names=F,quote=F,sep="\t")
        
        write.table(paste(cropsl[cc],gg,ha,paste(list.stats[[which(ggcms==gg)]],collapse=","),
                          length(which(best.mod==which(ggcms==gg))),
                          sum(area.irrf[best.mod==which(ggcms==gg)],na.rm=T),sum(area.irrf[best.mod==which(ggcms==gg)],na.rm=T)/sum(area.irrf,na.rm=T),
                          sum(production[best.mod==which(ggcms==gg)],na.rm=T),sum(production[best.mod==which(ggcms==gg)],na.rm=T)/sum(production,na.rm=T),
                          length(which(best.mod.rsme==which(ggcms==gg))),
                          sum(area.irrf[best.mod.rsme==which(ggcms==gg)],na.rm=T),sum(area.irrf[best.mod.rsme==which(ggcms==gg)],na.rm=T)/sum(area.irrf,na.rm=T),
                          sum(production[best.mod.rsme==which(ggcms==gg)],na.rm=T),sum(production[best.mod.rsme==which(ggcms==gg)],na.rm=T)/sum(production,na.rm=T),
                          collapse=",",sep=","),
                    file=tn,append=T,row.names=F,col.names=F,quote=F,sep=",")        
      }
      #write ensemble line to table
      write.table(paste(cropsl[cc],"ensemble best",ha,paste(list.stats[[length(ggcms)+1]],collapse=","),
                        length(which(best.mod>0)),
                        sum(area.irrf[best.mod>0],na.rm=T),sum(area.irrf[best.mod>0],na.rm=T)/sum(area.irrf,na.rm=T),
                        sum(production[best.mod>0],na.rm=T),sum(production[best.mod>0],na.rm=T)/sum(production,na.rm=T),
                        length(which(best.mod.rsme>0)),
                        sum(area.irrf[best.mod.rsme>0],na.rm=T),sum(area.irrf[best.mod.rsme>0],na.rm=T)/sum(area.irrf,na.rm=T),
                        sum(production[best.mod.rsme>0],na.rm=T),sum(production[best.mod.rsme>0],na.rm=T)/sum(production,na.rm=T),
                        collapse=",",sep=","),
                  file=tn,append=T,row.names=F,col.names=F,quote=F,sep=",")        
      

      if(F){
        png(paste0(picture.path,cropsl[cc],"_boxplot_rsme_",tolower(cl),"_",ha,"_1980_2010_shifted_ts_processed_dt.png"),
            height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
        par(mar=c(8,3,1,1))
        ra <- 0
        for(gg in 1:(length(ggcms)+1)){
          if(length(list.summarydata[[gg]])>0) ra <- range(ra,list.summarydata[[gg]]$rsme.stats$stats,na.rm=T)
        }
        #plot, no content
        plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
        axis(2)
        mtext("RMSE (tDM/ha)",2,line=2)
        box()
        #axis(1,labels=F)
        mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
        for(gg in 1:(length(ggcms)+1)){
          li <- list.summarydata[[gg]]
          #cat(ggcms[gg],str(li$rsme.stats),"\n")
          if(length(li)>0) bxp(li$rsme.stats,at=gg+0.5,add=T)
        }
        dev.off()
        
        png(paste0(picture.path,cropsl[cc],"_boxplot_mae_",tolower(cl),"_",ha,"_1980_2010_shifted_ts_processed_dt.png"),
            height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
        par(mar=c(8,3,1,1))
        ra <- 0
        for(gg in 1:(length(ggcms)+1)){
          if(length(list.summarydata[[gg]])>0) ra <- range(ra,list.summarydata[[gg]]$mae.stats$stats,na.rm=T)
        }
        #plot, no content
        plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
        axis(2)
        mtext("MAE (tDM/ha)",2,line=2)
        box()
        #axis(1,labels=F)
        mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
        for(gg in 1:(length(ggcms)+1)){
          li <- list.summarydata[[gg]]
          #cat(ggcms[gg],str(li$mae.stats),"\n")
          if(length(li)>0) bxp(li$mae.stats,at=gg+0.5,add=T)
        }
        dev.off()
        
        png(paste0(picture.path,cropsl[cc],"_boxplot_mae_area_weighted_",tolower(cl),"_",ha,"_1980_2010_shifted_ts_processed_dt.png"),
            height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
        par(mar=c(8,3,1,1))
        ra <- 0
        for(gg in 1:(length(ggcms)+1)){
          if(length(list.summarydata[[gg]])>0) ra <- range(ra,list.summarydata[[gg]]$mae_area.stats$stats,na.rm=T)
        }
        #plot, no content
        plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
        axis(2)
        mtext("MAE (tDM/ha)",2,line=2)
        box()
        #axis(1,labels=F)
        mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
        for(gg in 1:(length(ggcms)+1)){
          li <- list.summarydata[[gg]]
          #cat(ggcms[gg],str(li$mae_area.stats),"\n")
          if(length(li)>0) bxp(li$mae_area.stats,at=gg+0.5,add=T)
        }
        dev.off()
        
        png(paste0(picture.path,cropsl[cc],"_boxplot_mae_production_weighted_",tolower(cl),"_",ha,"_1980_2010_shifted_ts_processed_dt.png"),
            height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
        par(mar=c(8,3,1,1))
        ra <- 0
        for(gg in 1:(length(ggcms)+1)){
          if(length(list.summarydata[[gg]])>0) ra <- range(ra,list.summarydata[[gg]]$mae_production.stats$stats,na.rm=T)
        }
        #plot, no content
        plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
        axis(2)
        mtext("MAE (tDM/ha)",2,line=2)
        box()
        #axis(1,labels=F)
        mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
        for(gg in 1:(length(ggcms)+1)){
          li <- list.summarydata[[gg]]
          #cat(ggcms[gg],str(li$mae_production.stats),"\n")
          if(length(li)>0) bxp(li$mae_production.stats,at=gg+0.5,add=T)
        }
        dev.off()
        
        
        png(paste0(picture.path,cropsl[cc],"_boxplot_r2_",tolower(cl),"_",ha,"_1980_2010_shifted_ts_processed_dt.png"),
            height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
        par(mar=c(8,3,1,1))
        ra <- c(0,1)
        #plot, no content
        plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
        axis(2)
        mtext("R2 ",2,line=2)
        box()
        #axis(1,labels=F)
        mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
        for(gg in 1:(length(ggcms)+1)){
          li <- list.summarydata[[gg]]
          #cat(ggcms[gg],str(li$r2.stats),"\n")
          if(length(li)>0) bxp(li$r2.stats,at=gg+0.5,add=T)
        }
        dev.off()
        
      }
      
      mapi <- best.cor * best.cor
      mapi[best.mod==-9] <- 1.05
      save(mapi,file=paste0(picture.path,cropsl[cc],"_best_corrlation_", tolower(cl),
                                    "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata"))
      png(paste0(picture.path,cropsl[cc],"_best_corrlation_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            mapi[,360:1] ,main="",col=col.r2,zlim=c(0,1.05),oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(mapi,legend.only=T,col=col.r2[1:20],zlim=c(0,1),
                 oldstyle=TRUE,smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()
      rm(mapi)

      save(best.rsme,file=paste0(picture.path,cropsl[cc],"_best_rsme_", tolower(cl),
                                         "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_shifted_ts_processed_dt.Rdata"))
      png(paste0(picture.path,cropsl[cc],"_best_rsme_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            best.rsme[,360:1] ,main="",col=col.r2[20:1],oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(best.rsme,legend.only=T,col=col.r2[20:1],
                 oldstyle=TRUE,smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()

      mapi6 <- (best.cor * best.cor) - (best.cor2 * best.cor2)
      mapi6[best.mod==-9 | best.mod2==-9] <- 1.05

      save(mapi6,file=paste0(picture.path,cropsl[cc],"_dist_best2ndbest_corrlation_", tolower(cl),
                                     "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata"))      
      png(paste0(picture.path,cropsl[cc],"_dist_best2ndbest_corrlation_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            mapi6[,360:1] ,main="",col=col.r2,zlim=c(0,1.05),oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(mapi6,legend.only=T,col=col.r2[1:20],zlim=c(0,1),
                 oldstyle=TRUE,smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()
      rm(mapi6)
      gc()
      
      mapi2 <- best.mod
      mapi2[best.mod==-9] <- NA
      save(mapi2,file=paste0(picture.path,cropsl[cc],"_best_model_", tolower(cl),
                                     "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata"))
      png(paste0(picture.path,cropsl[cc],"_best_model_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            mapi2[,360:1] ,main="",col=col.mods,zlim=c(1,14),oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(mapi2,legend.only=T,col=col.mods,zlim=c(1,14),
                 oldstyle=TRUE,axis.args=list(at=1:14,labels=ggcms),
                 smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()
      rm(mapi2)
      gc()
      
      mapi3 <- mean/totalv*mean/totalv #R2 not R
      mapi3[!is.finite(mapi3)] <- NA
      save(mapi3,file=paste0(picture.path,cropsl[cc],"_mean_correlation_", tolower(cl),
                                     "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata"))
      png(paste0(picture.path,cropsl[cc],"_mean_correlation_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            mapi3[,360:1] ,main="",col=col.r2,zlim=c(0,1.05),oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(mapi3,legend.only=T,col=col.r2[1:20],zlim=c(0,1),
                 oldstyle=TRUE,smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()
      rm(mapi3)
      gc()
      
      mapi4 <- totalv/totaln
      mapi4[!is.finite(mapi4)] <- NA
      save(mapi4,file=paste0(picture.path,cropsl[cc],"_valid_models_", tolower(cl),
                                     "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata"))
      png(paste0(picture.path,cropsl[cc],"_valid_models_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            mapi4[,360:1] ,main="",col=col.r2,zlim=c(0,1.05),oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(mapi4,legend.only=T,col=col.r2[1:20],zlim=c(0,1),
                 oldstyle=TRUE,smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()
      rm(mapi4)
      gc()
      
      mapi5 <- num30/totalv
      mapi5[!is.finite(mapi5)] <- NA
      png(paste0(picture.path,cropsl[cc],"_valid_models_r30_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            mapi5[,360:1] ,main="",col=col.r2,zlim=c(0,1.05),oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(mapi5,legend.only=T,col=col.r2[1:20],zlim=c(0,1),
                 oldstyle=TRUE,smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()
      rm(mapi5)
      gc()

      png(paste0(picture.path,cropsl[cc],"_number_models_r30_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            num30[,360:1] ,main="",col=col.mods,zlim=c(1,14),oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(num30,legend.only=T,col=col.mods,zlim=c(1,14),
                 oldstyle=TRUE,smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()
      png(paste0(picture.path,cropsl[cc],"_number_models_sig_", tolower(cl),
                 "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.png"),
          height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
      par("mar"=c(2,2,2,12))
      image(x=seq(-179.75,179.75,len=720),y=seq(-89.75,89.75,len=360),ylim=c(-55,65),
            totalv[,360:1] ,main="",col=col.mods,zlim=c(1,14),oldstyle=TRUE)
      grid(10,10,lwd=.5)
      map(boundary=T,add=T)
      image.plot(totalv,legend.only=T,col=col.mods,zlim=c(1,14),
                 oldstyle=TRUE,smallplot=c(0.87,0.9,0.1,0.9))
      #close.screen(all=T)
      dev.off()
      
    }#harm
    # plot boxplots for all harmonization settings ####
    fn <- paste0(path.sim,"processed/masks/weight/",cropsl[cc], ".nc4")
    area.rf <- read.map.from.nc(fn,"rainfed")
    area.ir <- read.map.from.nc(fn,"irrigated")
    fname <- paste0(path.ref.detrend,cropss[cc],"_weight_ray_1961-2008.nc4")
    crop <- cropss[cc]
    r.ray <- get_nc4_ref_slice(fname,crop)
    
    ## detrending and correlation computation (using function detrended_cor)
    r.ray.sub <- r.ray[,,s.r:e.r]
    ## area
    area.irrf <- area.ir + area.rf
    ## production
    ray.mean <- apply(r.ray.sub,c(1,2),mean,na.rm=T)
    ray.mean[!is.finite(ray.mean)] <- NA
    production <- area.irrf * ray.mean
    production[!is.finite(production)] <- NA
    sum.area <- sum(area.irrf,na.rm=T)
    sum.production <- sum(production,na.rm=T)
    
    
    png(paste0(picture.path,cropsl[cc],"_boxplot_rsme_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- 0

    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        if(length(list.summarydata[[gg]])>0) ra <- range(ra,list.summarydata[[gg]]$rsme.stats$stats,list.summarydata[[gg]]$rsme.stats$out,na.rm=T)
      }
    }
    #plot, no content
    plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("RMSE (tDM/ha)",2,line=2)
    box()
    #axis(1,labels=F)
    mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$rsme.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=li$rsme.stats$n,srt=90,xpd=NA)
        }
      }
    }
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()

    png(paste0(picture.path,cropsl[cc],"_boxplot_mae_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- 0

    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        if(length(list.summarydata[[gg]])>0) ra <- range(ra,list.summarydata[[gg]]$mae.stats$stats,list.summarydata[[gg]]$mae.stats$out,na.rm=T)
      }
    }
    #plot, no content
    plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("MAE (tDM/ha)",2,line=2)
    box()
    #axis(1,labels=F)
    mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$mae.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=li$mae.stats$n,srt=90,xpd=NA)
        }
      }
    }
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()

    png(paste0(picture.path,cropsl[cc],"_boxplot_mae_areaw_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- 0

    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        if(length(list.summarydata[[gg]])>0) ra <- range(ra,list.summarydata[[gg]]$mae_area.stats$stats,list.summarydata[[gg]]$mae_area.stats$out,na.rm=T)
      }
    }
    #plot, no content
    plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("MAE (tDM/ha)",2,line=2)
    box()
    #axis(1,labels=F)
    mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$mae_area.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=format(li$mae_area.stats$n/sum.area*100,digits=3,scientific=F,nsmall=1),srt=90,xpd=NA)
        }
      }
    }
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()

    png(paste0(picture.path,cropsl[cc],"_boxplot_mae_prodw_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- 0

    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        if(length(list.summarydata[[gg]])>0) ra <- range(ra,list.summarydata[[gg]]$mae_production.stats$stats,list.summarydata[[gg]]$mae_production.stats$out,na.rm=T)
      }
    }
    #plot, no content
    plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("MAE (tDM/ha)",2,line=2)
    box()
    #axis(1,labels=F)
    mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$mae_production.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=format(li$mae_production.stats$n/sum.production,digits=3,scientific=F,nsmall=1),srt=90,xpd=NA)
        }
      }
    }
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()
    
    # invert lats to match format of other data
    rayr2 <- read.map.from.nc(paste0(path.ray.2015,"NonCategoricalFigure2",cropsr[cc],"05.nc"),"Data")[,360:1]
    png(paste0(picture.path,cropsl[cc],"_boxplot_r2_sig_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- c(0,1)

    #plot, no content
    plot(1:(length(ggcms)+5),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("R2",2,line=2)
    box()
    #axis(1,labels=F)
    #text(c(1:length(ggcms))+0.5, par("usr")[3] - 0.1, labels = ggcms, srt = 45, pos = 1, xpd = TRUE,adj=c(1,0))
    mtext(c(ggcms,"ensemble best","ensemble X Ray","Ray2015","Ray X ensemble"),1,las=2,adj=1,at=c(1:(length(ggcms)+4))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$r2.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=li$r2.stats$n,srt=90,xpd=NA)
        }
      }
      # read in ensemble map again to compute R2 distribution for area that is covered by Ray data
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_shifted_ts.Rdata")
      # avoid keeping too much junk, I'm loading this to a new environemnt, extract the map and delete the rest
      tmp.env <- new.env()
      load(sn,envir=tmp.env)
      best.cor <- get("best.cor",pos=tmp.env)
      rm(tmp.env)
      # mask all areas where Ray has no significant correlations
      best.x.ray <- best.cor*best.cor
      best.x.ray[!is.finite(rayr2)] <- NA
      if(ha==harms[1]){
        ray.x.best <- rayr2
        ray.x.best[!is.finite(best.cor)] <- NA
        b6a <- spBwplotStats(as.vector(ray.x.best),as.vector(area.irrf))
        b6p <- spBwplotStats(as.vector(ray.x.best),as.vector(production))
      }
      # computed weighted stats for next plots
      assign(paste0("b5a",which(harms==ha)),spBwplotStats(as.vector(best.x.ray),as.vector(area.irrf)))
      assign(paste0("b5p",which(harms==ha)),spBwplotStats(as.vector(best.x.ray),as.vector(production)))
      boxplot(as.vector(best.x.ray),at=length(ggcms)+2+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
      text(length(ggcms)+2+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=length(which(best.x.ray>0)),srt=90,xpd=NA)
    }
    boxplot(as.vector(rayr2),at=length(ggcms)+3+0.5,add=T,boxwex=0.2,border="black")
    text(length(ggcms)+3+0.5,ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=length(which(rayr2>0)),srt=90,xpd=NA)
    boxplot(as.vector(ray.x.best),at=length(ggcms)+4+0.5,add=T,boxwex=0.2,border="black")
    text(length(ggcms)+4+0.5,ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=length(which(ray.x.best>0)),srt=90,xpd=NA)
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()


    b4 <- spBwplotStats(as.vector(rayr2),as.vector(area.irrf))


    png(paste0(picture.path,cropsl[cc],"_boxplot_r2_areaw_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- c(0,1)

    #plot, no content
    plot(1:(length(ggcms)+5),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("R2",2,line=2)
    box()
    #axis(1,labels=F)
    #text(c(1:length(ggcms))+0.5, par("usr")[3] - 0.1, labels = ggcms, srt = 45, pos = 1, xpd = TRUE,adj=c(1,0))
    mtext(c(ggcms,"ensemble best","ensemble X Ray","Ray2015","Ray X ensemble"),1,las=2,adj=1,at=c(1:(length(ggcms)+4))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$r2_area.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=format(li$r2_area.stats$n/sum.area,digits=3,scientific=F,nsmall=1),srt=90,xpd=NA)
        }
      }
      b5 <- get(paste0("b5a",which(harms==ha)))
      bxp(b5,at=length(ggcms)+2+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
      text(length(ggcms)+2+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=format(b5$n/sum.area,digits=3,scientific=F,nsmall=1),srt=90,xpd=NA)
    }
    bxp(b4,at=length(ggcms)+3+0.5,add=T,boxwex=0.2,border="black")
    text(length(ggcms)+3+0.5,ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=format(b4$n/sum.area,digits=3,scientific=F,nsmall=1),srt=90,xpd=NA)
    bxp(b6a,at=length(ggcms)+4+0.5,add=T,boxwex=0.2,border="black")
    text(length(ggcms)+4+0.5,ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=format(b6a$n/sum.area,digits=3,scientific=F,nsmall=1),srt=90,xpd=NA)
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()

    b4 <- spBwplotStats(as.vector(rayr2),as.vector(production))
    png(paste0(picture.path,cropsl[cc],"_boxplot_r2_prodw_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=10,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- c(0,1)

    #plot, no content
    plot(1:(length(ggcms)+5),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("R2",2,line=2)
    box()
    #axis(1,labels=F)
    #text(c(1:length(ggcms))+0.5, par("usr")[3] - 0.1, labels = ggcms, srt = 45, pos = 1, xpd = TRUE,adj=c(1,0))
    mtext(c(ggcms,"ensemble best","ensemble X Ray","Ray2015","Ray X ensemble"),1,las=2,adj=1,at=c(1:(length(ggcms)+4))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$r2_production.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.3,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,
               adj=c(1,0.5),srt=90,xpd=NA,cex=0.7,
               labels=format(li$r2_production.stats$n/sum.production,digits=3,scientific=F,nsmall=1))
        }
      }
      b5 <- get(paste0("b5p",which(harms==ha)))
      bxp(b5,at=length(ggcms)+2+0.1+0.2*which(harms==ha),add=T,boxwex=0.3,border=col.box[which(harms==ha)])
      text(length(ggcms)+2+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,
           adj=c(1,0.5),srt=90,xpd=NA,cex=0.7,
           labels=format(b5$n/sum.production,digits=3,scientific=F,nsmall=1))
    }
    bxp(b4,at=length(ggcms)+3+0.5,add=T,boxwex=0.3,border="black")
    text(length(ggcms)+3+0.5,ra[2]+(ra[2]-ra[1])*0.155,
         adj=c(1,0.5),srt=90,xpd=NA,cex=0.7,
         labels=format(b4$n/sum.production,digits=3,scientific=F,nsmall=1))
    bxp(b6p,at=length(ggcms)+4+0.5,add=T,boxwex=0.3,border="black")
    text(length(ggcms)+4+0.5,ra[2]+(ra[2]-ra[1])*0.155,
         adj=c(1,0.5),srt=90,xpd=NA,cex=0.7,
         labels=format(b6p$n/sum.production,digits=3,scientific=F,nsmall=1))
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()


    png(paste0(picture.path,cropsl[cc],"_boxplot_r2_sig_areaw_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- c(0,1)

    #plot, no content
    plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("R2",2,line=2)
    box()
    #axis(1,labels=F)
    #text(c(1:length(ggcms))+0.5, par("usr")[3] - 0.1, labels = ggcms, srt = 45, pos = 1, xpd = TRUE,adj=c(1,0))
    mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$r2_sig_area.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=format(li$r2_sig_area.stats$n/sum.area,digits=3,scientific=F,nsmall=1),srt=90,xpd=NA)
        }
      }
    }
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()
    png(paste0(picture.path,cropsl[cc],"_boxplot_r2_sig_prodw_allharms_",tolower(cl),"_1980_2010_shifted_ts_processed_dt.png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo")
    par(mar=c(8,3,5,1),lwd=2)
    ra <- c(0,1)

    #plot, no content
    plot(1:(length(ggcms)+2),ylim=ra,type="n",axes=F,xlab="",ylab="")
    axis(2)
    mtext("R2",2,line=2)
    box()
    #axis(1,labels=F)
    #text(c(1:length(ggcms))+0.5, par("usr")[3] - 0.1, labels = ggcms, srt = 45, pos = 1, xpd = TRUE,adj=c(1,0))
    mtext(c(ggcms,"ensemble best"),1,las=2,adj=1,at=c(1:(length(ggcms)+1))+0.5,line=1)
    for(ha in harms){
      sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
      load(sn)
      for(gg in 1:(length(ggcms)+1)){
        li <- list.summarydata[[gg]]
        if(length(li)>0) {
          bxp(li$r2_sig_production.stats,at=gg+0.1+0.2*which(harms==ha),add=T,boxwex=0.2,border=col.box[which(harms==ha)])
          text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,adj=c(1,0.5),labels=format(li$r2_sig_production.stats$n/sum.production,digits=3,scientific=F,nsmall=1),srt=90,xpd=NA)
        }
      }
    }
    legend("bottomleft",horiz=T,legend=harms,bty="n",fill="white",border=col.box)
    dev.off()
    print(ls())
  }#crops
}#clim


