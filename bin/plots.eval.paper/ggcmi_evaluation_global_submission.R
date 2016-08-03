rm(list=ls(all=TRUE))
options(warn=1)

##########################################################
# libraries ####
##########################################################
require(fields)
require(maps)
require(ncdf4)

##########################################################
# global settings ####
##########################################################
NODATA <- 1e20
ncell <- 67420
cr0 <- "none"
dt0 <- "ma"
mp0 <- "false"
scen0 <- "default"

shift.threshold <- 0.2

ggcms <- c("pDSSAT","EPIC-Boku","EPIC-IIASA","GEPIC",
           "pAPSIM","PEGASUS","LPJ-GUESS","LPJmL",
           "CGMS-WOFOST","CLM-Crop","EPIC-TAMU","ORCHIDEE-crop",
           "PEPIC","PRYSBI2")

clim <- c("AgMERRA","WFDEI.GPCC","WATCH")
harms <- c("default", "harmnon", "fullharm") 
harms2 <- c("default", "harm-suffN", "fullharm") 

cropsl <- c("maize","wheat","rice","soy")
cropss <- c("mai","whe","ric","soy")
FRESHMATTER <- 100 / c(88, 88, 87, 91) 

aggs <- c("fixed_mirca_mask","dynamic_ray_mask","fixed_iizumi_mask","fixed_spam_mask")

# time range
years <- c(1961:2010)
start.year <- 1982 #1980 corrected to remove the first 2 NA years
end.year <- 2006 #2010 corrected to remove the last 2 NA years and also to make sure all masks end the same (Ray is only 2008)
r.years <- c(1961:2012)
sim.years <- c(1980:2010)
s.r <- which(r.years==start.year)
e.r <- which(r.years==end.year)
s.sim <- which(sim.years==start.year)
e.sim <- which(sim.years==end.year)

path.sim <- "/p/projects/macmit/data/GGCMI/AgMIP.output/processed/biascorr/global/faostat/"
path.ens <- "/p/projects/macmit/data/GGCMI/AgMIP.output/processed/modelensemble/global/faostat/"
path.ref <- "/p/projects/macmit/data/GGCMI/AgMIP.input/other.inputs/reference/faostat/"
path.pic <- "/p/projects/macmit/users/cmueller/GGCMI/global_agg/"
#########################################################
# functions ####
#########################################################

get_nc4_data_slice <- function(fname,crd=cr0,dtd=dt0,mpd=mp0,scend=scen0){
  nf <- nc_open(fname)
  cr <- which(strsplit(ncatt_get(nf,varid="cr")$long_name,split=", ")[[1]]==crd)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  scen <- which(strsplit(ncatt_get(nf,varid="scen")$long_name,split=", ")[[1]]==scend)
  if(length(scen)==0){
    nc_close(nf)
    return(NA)
  }
  # order is cr, mp, dt, scen, time, fpu
  # test if scen dimension is > 1, otherwise object has one dimension less
  if(length(strsplit(ncatt_get(nf,varid="scen")$long_name,split=", ")[[1]])>1) {
    data.bc <- ncvar_get(nf,varid="yield_detrend")[cr,mp,dt,scen,]
  } else {
    data.bc <- ncvar_get(nf,varid="yield_detrend")[cr,mp,dt,]
  }
  nc_close(nf)
  data.bc
}

get_nc4_ensemble_slice <- function(fname,crd=cr0,dtd=dt0,mpd=mp0,scend=scen0,nm=1,wtd="unweighted"){
  nf <- nc_open(fname)
  cr <- which(strsplit(ncatt_get(nf,varid="cr")$long_name,split=", ")[[1]]==crd)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  wt <- which(strsplit(ncatt_get(nf,varid="wt")$long_name,split=", ")[[1]]==wtd)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  scen <- which(strsplit(ncatt_get(nf,varid="top_scens")$long_name,split=", ")[[1]]==scend)
  # order is wt, nm, cr, mp, dt, time
  # test if scen dimension is > 1, otherwise object has one dimension less
  data.bc <- ncvar_get(nf,varid="yield_detrend")[wt,nm,cr,mp,dt,]
  nc_close(nf)
  data.bc
}

get_nc4_ref_slice <- function(fname,crop,dtd=dt0,mpd=mp0){
  nf <- nc_open(fname)
  dt <- which(strsplit(ncatt_get(nf,varid="dt")$long_name,split=", ")[[1]]==dtd)
  mp <- which(strsplit(ncatt_get(nf,varid="mp")$long_name,split=", ")[[1]]==mpd)
  data.bc <- ncvar_get(nf,varid=paste0("yield_",crop))[mp,dt,]
  nc_close(nf)
  data.bc
}

mean.bias <- function(x1,x2){
  x12 <- x1*x2
  if(length(x12[which(is.finite(x12))])<1)
    return(NA)
  mean.bias <- mean((x2-x1),na.rm=TRUE)
  mean.bias
}

##################################################################
# colors ####
##################################################################
col.mods2 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#ffed6f','gray50',2,3)
col.mods <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#ffed6f','gray50')



##################################################################
# main ####
##################################################################
cl <- clim[1]
for(cro in cropss){
  for(ha in harms){
    data <- array(NA,dim=c(length(s.sim:e.sim),length(ggcms),length(aggs)))
    for(agg in aggs){
      ggcm.names <- NULL
      ref.data <- get_nc4_ref_slice(paste0(path.ref,"faostat.1961-2012.global.",agg,".nc4"),cro)[s.r:e.r]
      for(ggcm in ggcms){
        ggcm.names <- c(ggcm.names,ggcm)
        fn <- paste0(path.sim,agg,"/",tolower(ggcm),"_",tolower(cl),"_hist_",cro,"_annual_1980_2010.biascorr.nc4")
        if(!file.exists(fn)) next
        buf <- get_nc4_data_slice(fn,scend=ha)[s.sim:e.sim]
        if(ha==harms[1] & all(is.na(buf))){ # for the default plot, add fullharm for those that don't have a default
          buf <- get_nc4_data_slice(fn,scend=harms[3])[s.sim:e.sim]
          if(!all(is.na(buf))) ggcm.names[which(ggcm.names==ggcm)] <- paste0(ggcm," (",harms[3],")")
        }
        data[,which(ggcms==ggcm),which(aggs==agg)] <- buf *FRESHMATTER[which(cropss==cro)]
      }      
    }
    png(paste0(path.pic,"global_aggregation_best_mask_",cro,"_",ha,"_",cl,"_shifted_ts_no_ensemble.png"),
        width=8*300,height=5*300,res=300,type="cairo",pointsize=10)
    ra <- range(ref.data,data,na.rm=T)
    # add extra space at bottom for legend
    ra[1] <- ra[1] - (ra[2]-ra[1])*0.2
    plot(c(start.year:end.year),ref.data,ylim=ra,type="l",
         lwd=2.5,lty=3,ylab="t/ha",xlab="year (AD)",
         xaxp=c(start.year,end.year,end.year-start.year))
    cors <- array(NA,length(ggcms))
    rs <- array(NA,length(ggcms))
    for(gg in 1:length(ggcms)){
      cor4b <- cor4m1 <- cor4p1 <- NULL
      r4b <- r4m1 <- r4p1 <- NULL
      if(all(is.na(data[,gg,]))) next
      for(agg in 1:length(aggs)){
        # compute correlations with and without shifted time series
        buf <- cor.test(ref.data,data[,gg,agg],use="pairwise")
        cor4b <- c(cor4b,buf$estimate)
        r4b <- c(r4b,buf$p.value)
        buf <- cor.test(ref.data[-1],data[-dim(data)[1],gg,agg],use="pairwise")
        cor4m1 <- c(cor4m1,buf$estimate)
        r4m1 <- c(r4m1,buf$p.value)
        buf <- cor.test(ref.data[-length(ref.data)],data[-1,gg,agg],use="pairwise")
        cor4p1 <- c(cor4p1,buf$estimate)
        r4p1 <- c(r4p1,buf$p.value)
      } 
      cat(ggcms[gg],cor4b,cor4m1,cor4p1,"\n")
      # plot line shifted by 1 if cor is improved by at least 0.2
      if(max(cor4m1,na.rm=T) > (max(cor4b,na.rm=T)+shift.threshold) & max(cor4m1,na.rm=T)>max(cor4p1,na.rm=T)){
        cat("using SF for",ggcm.names[gg],ha,cro,max(cor4b,na.rm=T),"<",max(cor4m1,na.rm=T),"\n")
        cor4 <- cor4m1
        r4 <- r4m1
        lines(c(start.year:end.year),c(NA,data[-dim(data)[1],gg,which.max(cor4)]),col=col.mods[gg],lwd=1.5,lty=if(gg>length(ggcms)) 3 else 1)
        # adjust names for legend, add SF for "shifted forward"
        if(gg <= length(ggcms) & ggcm.names[gg] != ggcms[gg]){
          ggcm.names[gg] <- paste0(ggcms[gg]," (sf) (",harms[3],")")
        } else {
          ggcm.names[gg] <- paste0(ggcms[gg]," (sf)")
        }
      } else if (max(cor4p1,na.rm=T) > (max(cor4b,na.rm=T)+shift.threshold) & max(cor4p1,na.rm=T)>max(cor4m1,na.rm=T)){
        cat("using SB for",ggcm.names[gg],ha,cro,max(cor4b,na.rm=T),"<",max(cor4p1,na.rm=T),"\n")
        cor4 <- cor4p1
        r4 <- r4p1
        lines(c(start.year:end.year),c(data[-1,gg,which.max(cor4)],NA),col=col.mods[gg],lwd=1.5,lty=if(gg>length(ggcms)) 3 else 1)
        # adjust names for legend, add SB for "shifted backward"
        if(gg <= length(ggcms) & ggcm.names[gg] != ggcms[gg]){
          ggcm.names[gg] <- paste0(ggcms[gg]," (sb) (",harms[3],")")
        } else {
          ggcm.names[gg] <- paste0(ggcms[gg]," (sb)")
        }
      } else {
        cor4 <- cor4b
        r4 <- r4b
        lines(c(start.year:end.year),data[,gg,which.max(cor4)],col=col.mods[gg],lwd=1.5,lty=if(gg>length(ggcms)) 3 else 1)
      }
      cors[gg] <- max(cor4,na.rm=T) #cor(ref.data,data[,gg],use="pairwise")
      rs[gg] <- r4[which.max(cor4)]
    }
    # plot reference again on top of all
    lines(c(start.year:end.year),ref.data,lwd=2.5,lty=3)
    legend("bottomright",ncol=4,lwd=c(rep(2,length(ggcm.names)),2.5),lty=c(rep(1,length(ggcms)),3),
           bty="n",col=c(col.mods,1),cex=0.6,
           legend=paste(c(ggcm.names,"FAOstat"),c(format(cors,digits=1,nsmall=3),""),
                        c(ifelse(rs<0.001,"***",ifelse(rs<0.05,"**",ifelse(rs<0.1,"*","n.s."))),"")))
    dev.off()
  }
}


# mean.bias vs. r scatterplots
for(cro in cropss){
  cors <- array(NA,dim=c(length(ggcms),length(harms)))
  rs <- array(NA,dim=c(length(ggcms),length(harms)))
  mean.biases <- array(NA,dim=c(length(ggcms),length(harms)))
  for(ha in harms){
    data <- array(NA,dim=c(length(s.sim:e.sim),length(ggcms),length(aggs)))
    for(agg in aggs){
      ggcm.names <- NULL
      ref.data <- get_nc4_ref_slice(paste0(path.ref,"faostat.1961-2012.global.",agg,".nc4"),cro,mpd="true")[s.r:e.r]
      for(ggcm in ggcms){
        ggcm.names <- c(ggcm.names,ggcm)
        fn <- paste0(path.sim,agg,"/",tolower(ggcm),"_",tolower(cl),"_hist_",cro,"_annual_1980_2010.biascorr.nc4")
        if(!file.exists(fn)) next
        buf <- get_nc4_data_slice(fn,scend=ha,mp="true")[s.sim:e.sim]
        if(ha==harms[1] & all(is.na(buf))){ # for the default plot, add fullharm for those that don't have a default
          buf <- get_nc4_data_slice(fn,scend=harms[3],mp="true")[s.sim:e.sim]
          if(!all(is.na(buf))) ggcm.names[which(ggcm.names==ggcm)] <- paste0(ggcm," (",harms[3],")")
        }
        data[,which(ggcms==ggcm),which(aggs==agg)] <- buf *FRESHMATTER[which(cropss==cro)]
      }      
    }
    for(gg in 1:length(ggcms)){
      cor4b <- cor4m1 <- cor4p1 <- NULL
      r4b <- r4m1 <- r4p1 <- NULL
      mean.bias4b <- NULL
      if(all(is.na(data[,gg,]))) next
      for(agg in 1:length(aggs)){
        buf <- cor.test(ref.data,data[,gg,agg],use="pairwise")
        cor4b <- c(cor4b,buf$estimate)
        r4b <- c(r4b,buf$p.value)
        mean.bias4b <- c(mean.bias4b,mean.bias(ref.data,data[,gg,agg]))
        buf <- cor.test(ref.data[-1],data[-dim(data)[1],gg,agg],use="pairwise")
        cor4m1 <- c(cor4m1,buf$estimate)
        r4m1 <- c(r4m1,buf$p.value)
        buf <- cor.test(ref.data[-length(ref.data)],data[-1,gg,agg],use="pairwise")
        cor4p1 <- c(cor4p1,buf$estimate)
        r4p1 <- c(r4p1,buf$p.value)
      } 
      cat(ggcms[gg],cor4b,cor4m1,cor4p1,"\n")
      # plot line shifted by 1 if cor is improved by at least 0.2
      if(max(cor4m1,na.rm=T) > (max(cor4b,na.rm=T)+shift.threshold) & max(cor4m1,na.rm=T)>max(cor4p1,na.rm=T)){
        cat("using SF for",ggcm.names[gg],ha,cro,max(cor4b,na.rm=T),"<",max(cor4m1,na.rm=T),"\n")
        cor4 <- cor4m1
        r4 <- r4m1
        mean.bias4 <- mean.bias4b
        if(gg <= length(ggcms) & ggcm.names[gg] != ggcms[gg]){
          ggcm.names[gg] <- paste0(ggcms[gg]," (sf) (",harms[3],")")
        } else {
          ggcm.names[gg] <- paste0(ggcms[gg]," (sf)")
        }
      } else if (max(cor4p1,na.rm=T) > (max(cor4b,na.rm=T)+shift.threshold) & max(cor4p1,na.rm=T)>max(cor4m1,na.rm=T)){
        cat("using SB for",ggcm.names[gg],ha,cro,max(cor4b,na.rm=T),"<",max(cor4p1,na.rm=T),"\n")
        cor4 <- cor4p1
        r4 <- r4p1
        mean.bias4 <- mean.bias4b
        if(gg <= length(ggcms) & ggcm.names[gg] != ggcms[gg]){
          ggcm.names[gg] <- paste0(ggcms[gg]," (sb) (",harms[3],")")
        } else {
          ggcm.names[gg] <- paste0(ggcms[gg]," (sb)")
        }
      } else {
        cor4 <- cor4b
        r4 <- r4b
        mean.bias4 <- mean.bias4b
      }
      cors[gg,which(harms==ha)] <- max(cor4,na.rm=T) #cor(ref.data,data[,gg],use="pairwise")
      rs[gg,which(harms==ha)] <- r4[which.max(cor4)]
      mean.biases[gg,which(harms==ha)] <- mean.bias4[which.max(cor4)]
    }
  }
  # one plot for all harms, GGCMS, best masks, best r
  png(paste0(path.pic,"global_aggregation_best_mask_",cro,"_",cl,"_shifted_ts_no_ensemble_r_vs_mean.bias.png"),
      width=8*300,height=5*300,res=300,type="cairo",pointsize=10)
  ra <- range(mean.biases,na.rm=T)
  plot(1:1,ylim=ra,xlim=c(0,1.1),type="n",ylab="mean bias [t/ha]",xlab="correlation (r)")
  for(ha in 1:length(harms)){
    for(ggcm in 1:length(ggcms)){
      points(cors[ggcm,ha],mean.biases[ggcm,ha],
             pch=22+ha,col=col.mods[ggcm],bg=col.mods[ggcm])
    }
  }
  model <- lm(as.vector(mean.biases) ~ as.vector(cors))
  abline(model,lty=2,col="grey30")
  r2 <- paste(format(summary(model)$r.squared,nsmall=2,digits=1))
  pv <- summary(model)$coefficients[2,4]
  sig <- ifelse(pv<0.001,"***",ifelse(pv<0.05,"**",ifelse(pv<0.1,"*","n.s.")))
  
  legend("topright",ncol=1,pch=c(rep(15,length(ggcms)),NA,23:25),
         #bty="n",
         col=c(col.mods[1:length(ggcms)],NA,rep("grey20",3)),cex=0.8,
         legend=c(ggcms,"",harms2))
  legend("bottomright",bty="n",legend=bquote(R^2 ~ ":" ~ .(r2) ~ .(sig)),lty=2,col="grey30")
  dev.off()
}


