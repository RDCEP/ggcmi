rm(list=ls(all=TRUE))


require(ncdf4)
require(gplots)



# script to plot heatmap diagrams for GGCMI phase 1. 
# written by Christoph Müller, PIK
# cmueller@pik-potsdam.de
# I'm using a slightly modified version of heatmap.2() of the gplots package
# to allow for a second y-axis labeling

# intermediate version
# TODO
# include ensemble

#path <- "/iplex/01/2011/isimip-lpj/yields/GGCMI2_outputs/WFDEI_GPCC_default/ncdf/"
path <- "/iplex/01/2014/macmit/data/GGCMI/AgMIP.output/"
path.out <- "/iplex/01/2014/macmit/users/cmueller/GGCMI/"
path.ref <- "/iplex/01/2014/macmit/data/GGCMI/reference/"
path.ref.detrend <- "/iplex/01/2014/macmit/data/GGCMI/reference/ray-dt/"
#path <- "M:/data/GGCMI/AgMIP.output/"
path <- "D:/data/GGCMI/"
path.out <- "D:/publications/ggcmi_evaluation_paper/heatmaps/"
processed.ts <- "processed"


# loop through weather products TODO
clim <- c("agmerra","wfdei.gpcc","watch")
cfirst <- c(1980,1979,1958)
clast <- c(2010,2009,2001)

aggs <- c("fixed_mirca_mask","dynamic_ray_mask","fixed_iizumi_mask","fixed_spam_mask")
dts <- c("none","lin","quad","ma","ffdtr")
scens <- c("default","fullharm","harmnon")

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


#### define settings ####

cr0 <- "none"#"mean-scale"#"variance-scale"#"none"
dt0 <- "quad" #"none"#"quad"
mp0 <- "false"
scen0 <- "default"
nm0 <- 1 #number of top model for the ensemble
wt0 <- "unweighted" #wheighting for ensemble (1:unweighted, 2:weighted)

# read GADM0 names
gadm0.names <- array("missing",253)
gadm0.index <- read.csv(paste(path,"processed.150505","/masks/aggr/gadm0.meta.csv",sep=""),header=F)[,1]
name <- read.csv(paste(path,"processed.150505","/masks/aggr/gadm0.meta.csv",sep=""),header=F)[,2]

# fu*k. different country selection in sim files
fname<-paste(path,processed.ts,"/multimetrics/gadm0/faostat/",agg.fao,"/","pdssat","_",clim[1],"_hist_mai_annual_",cfirst[1],"_",
             clast[1],".multimetrics.nc4",sep="")
nf <- nc_open(fname)
gadm0.id2 <- ncvar_get(nf,varid="gadm0")
nc_close(nf)
sset <- which(gadm0.index %in% gadm0.id2)
gadm0.names <- name[sset]
gadm0.id <- 1:length(gadm0.names)



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

topproducer.soy <- c("United States","Brazil","Argentina","China","India","Paraguay","Canada","Uruguay","Ukraine","Bolivia")
gadm0.soy <- c(240,32,11,48,105,175,41,242,237,27)
topproducer.ric <- c("China","India","Indonesia","Bangladesh","Vietnam","Thailand","Myanmar","Philippines","Brazil","Japan")
gadm0.ric <- c(48,105,106,19,247,226,154,177,32,114)
topproducer.mai <- c("United States","China","Brazil","Argentina","Mexico","India","Ukraine","Indonesia","France","South_Africa")
gadm0.mai <- c(240,48,32,11,145,105,237,106,79,209)
topproducer.whe <- c("China","India","United States","Russia","France","Canada","Australia","Pakistan","Germany","Turkey")
gadm0.whe <- c(48,105,240,186,79,41,14,170,86,232)



colo <- terrain.colors(256)[265:1]

prefix <- paste("cr_",cr0,".dt_",dt0,".mp_",mp0,".multimetrics",sep="")

for(cc in crops){
  hm <- array(NA,dim=c(length(gadm0.names),length(ggcms)*3))
  for(cl in 1){        
    for(gg in ggcms){
      fn <- paste(path,processed.ts,
                  "/multimetrics/gadm0/faostat/fixed_mirca_mask/",
                  gg,"_",clim[cl],"_hist_",cc,"_annual_",cfirst[cl],"_",
                  clast[cl],".multimetrics.nc4",sep="")
      if(!is.na(file.info(fn)$size)){
        for(sc in scens){
          hm[,(which(ggcms==gg)-1)*3+which(scens==sc)] <- get_nc4_multimetrics(fn,scend=sc)
        }
      }
    }#ggcms
  }#cl
  rownames(hm) <- gadm0.names
  colnames(hm) <- as.vector(outer(scens,ggcms, paste, sep="."))

  # remove all GGCM-scenario combinations that are non-existent
  hm <- hm[,colSums(is.na(hm))<nrow(hm)]
  # remove all countries that are not covered by all models
  hm <- hm[rowSums(is.na(hm))==0,]
  best <- apply(hm,MARGIN=1,max)
  hm2 <- cbind(best,hm)
  colnames(hm2)[1] <- "best"
  best.lab <- as.vector(apply(hm,MARGIN=1,function(x) names(which.max(x))))
  
  # all countries ####
  png(paste(path.out,"heatmap_all_countries_",cc,"_",prefix,".png",sep=""),width=6*300,height=10*300,res=300,pointsize=5)
  par(xpd=NA,cex=1.5)
  heatmap.2b(hm2,Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(best.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=1:dim(hm)[2],rowsep=1:dim(hm)[1],sepcol="white",
             # set key to below heatmat, see
             #http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(15,12),col=colo,cexRow=1.5,cexCol=1.5)
  dev.off()
  
  #top10 only ####
  selec <- which(rownames(hm) %in% get(paste("topproducer",cc,sep=".")))
  best.lab <- as.vector(apply(hm[selec,],MARGIN=1,function(x) names(which.max(x))))
  png(paste(path.out,"heatmap_top10_",cc,"_",prefix,".png",sep=""),width=6*300,height=6*300,res=300,pointsize=5)
  par(xpd=NA,cex=1.5)
  heatmap.2b(hm2[selec,],Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(best.lab):1],
             dendrogram="none", trace="none", key=T,
             colsep=1:dim(hm2[selec,])[2],rowsep=1:dim(hm2[selec,])[1],sepcol="white",
             # set key to below heatmat, see
             #http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(15,12),col=colo,cexRow=1.5,cexCol=1.5)
  dev.off()
  
  # top 10 fasttrack only #### 
  
  selec <- which(rownames(hm) %in% get(paste("topproducer",cc,sep=".")))
  ggsel <- which(colnames(hm) %in% c("default.pdssat","default.gepic","default.epic-boku",
                                      "default.lpj-guess","default.lpjml","default.pegasus"))
  best <- apply(hm[,ggsel],MARGIN=1,max)
  best.lab <- as.vector(apply(hm[,ggsel],MARGIN=1,function(x) names(which.max(x))))[selec]
  hm2 <- cbind(best,hm[,ggsel])
  colnames(hm2)[1] <- "best"
  png(paste(path.out,"heatmap_top10_",cc,"_fasttrack_ggcm_",prefix,".png",sep=""),width=6*300,height=6*300,res=300,pointsize=5)
  par(xpd=NA,cex=1.5)
  heatmap.2b(hm2[selec,],Rowv=NULL,Colv=NULL,#keysize=0.5,
             second.label=best.lab[length(selec):1],
             dendrogram="none", trace="none", key=T,
             colsep=1:dim(hm2[selec,])[2],rowsep=1:dim(hm2[selec,])[1],sepcol="white",
             # set key to below heatmat, see
             #http://stackoverflow.com/questions/15351575/moving-color-key-in-r-heatmap-2-function-of-gplots-package
             lmat = rbind(c(0,4),c(0,3),c(2,1)),
             lwid=c(.5,4),lhei=c(.3,.1,4),
             sepwidth=c(0.01,0.01),key.par=list(mar=c(2,2,2,12)),
             margins=c(15,12),col=colo,cexRow=1.5,cexCol=1.5)
  dev.off()
  
  
}#crops
              