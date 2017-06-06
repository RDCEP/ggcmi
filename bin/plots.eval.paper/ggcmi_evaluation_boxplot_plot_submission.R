# (C) 2014-2017 Potsdam Institute for Climate Impact Research (PIK),
# written by Christoph Mueller, PIK
# cmueller@pik-potsdam.de
# Licensed under GNU AGPL Version 3 <LICENSE.txt in ggcmi directory>

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
harms2 <- c("default", "harm-suffN", "fullharm")

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

cl <- clim[1]

read.map.from.nc <- function(fn,vn){  
  vec <- array(0,dim=c(ncell,2))
  ncs <- nc_open(fn)
  data <- ncvar_get(ncs,varid=vn)
  nc_close(ncs)
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

# from the package simPopulation, which is no longer available via CRAN
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


for(cc in 1:length(cropsl)){
  
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

  #### use this for paper ####
  # invert lats to match format of other data
  rayr2 <- read.map.from.nc(paste0(path.ray.2015,"NonCategoricalFigure2",cropsr[cc],"05.nc"),"Data")[,360:1]
  b4 <- spBwplotStats(as.vector(rayr2),as.vector(production))  
  png(paste0(picture.path,cropsl[cc],"_boxplot_r2_prodw_allharms_",tolower(cl),"_1980_2010.png"),
      height=5*300,width=8*300,res=300,pointsize=10,type="cairo")
  par(mar=c(8,3,5,1),lwd=2)
  ra <- c(0,1)
  
  #plot, no content
  plot(1:(length(ggcms)+5),ylim=ra,type="n",axes=F,xlab="",ylab="")
  axis(2)
  mtext(expression(paste("R"^2)),2,line=2)
  box()
  mtext(c(ggcms,"ensemble best","ensemble X Ray","Ray2015","Ray X ensemble"),1,las=2,adj=1,at=c(1:(length(ggcms)+4))+0.5,line=1)
  for(ha in harms){
    sn <- paste0(picture.path,cropsl[cc],"_",tolower(cl),"_hist_",ha,"_yield_boxdata_shifted_ts.Rdata")
    load(sn)
    for(gg in 1:(length(ggcms)+1)){
      li <- list.summarydata[[gg]]
      if(length(li)>0) {
        bxp(li$r2_production.stats,at=gg+0.1+0.2*which(harms==ha),add=T,
            boxlwd=1,whisklwd=1,staplelwd=1,medlwd=3,
            outcex=0.2,outpch=20,boxwex=0.3,border=col.box[which(harms==ha)])
        text(gg+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,
             adj=c(1,0.5),srt=90,xpd=NA,cex=0.7,
             labels=format(li$r2_production.stats$n/sum.production,digits=3,scientific=F,nsmall=1))
      }
    }
    b5 <- get(paste0("b5p",which(harms==ha)))
    bxp(b5,at=length(ggcms)+2+0.1+0.2*which(harms==ha),add=T,
        boxlwd=1,whisklwd=1,staplelwd=1,medlwd=3,
        outcex=0.2,outpch=20,boxwex=0.3,border=col.box[which(harms==ha)])
    text(length(ggcms)+2+0.1+0.2*which(harms==ha),ra[2]+(ra[2]-ra[1])*0.155,
         adj=c(1,0.5),srt=90,xpd=NA,cex=0.7,
         labels=format(b5$n/sum.production,digits=3,scientific=F,nsmall=1))
  }
  bxp(b4,at=length(ggcms)+3+0.5,add=T,
      boxlwd=1,whisklwd=1,staplelwd=1,medlwd=3,
      outcex=0.2,outpch=20,boxwex=0.3,border="black")
  text(length(ggcms)+3+0.5,ra[2]+(ra[2]-ra[1])*0.155,
       adj=c(1,0.5),srt=90,xpd=NA,cex=0.7,
       labels=format(b4$n/sum.production,digits=3,scientific=F,nsmall=1))
  bxp(b6p,at=length(ggcms)+4+0.5,add=T,
      boxlwd=1,whisklwd=1,staplelwd=1,medlwd=3,
      outcex=0.2,outpch=20,boxwex=0.3,border="black")
  text(length(ggcms)+4+0.5,ra[2]+(ra[2]-ra[1])*0.155,
       adj=c(1,0.5),srt=90,xpd=NA,cex=0.7,
       labels=format(b6p$n/sum.production,digits=3,scientific=F,nsmall=1))
  legend("bottomleft",horiz=T,legend=harms2,bty="n",fill="white",border=col.box)
  dev.off()
  
}#crops
