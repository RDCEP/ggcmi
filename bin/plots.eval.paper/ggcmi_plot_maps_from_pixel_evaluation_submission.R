# (C) 2014-2017 Potsdam Institute for Climate Impact Research (PIK),
# written by Christoph Mueller, PIK
# cmueller@pik-potsdam.de
# Licensed under GNU AGPL Version 3 <LICENSE.txt in ggcmi directory>

##########################################################
# script to plot maps from 
# correlation_plots_JB_shifted_ts.R
# for publication
# combining different maps from *.Rdata files to one 
# figure file
##########################################################
rm(list=ls(all=T))


require(fields)
require(maps)
require(ncdf4)

picture.path <- "/p/projects/macmit/users/cmueller/GGCMI_cormap_test/"
path.ray.2015 <- "/p/projects/macmit/users/cmueller/GGCMI/ray2015/"


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
NODATA <- -9999

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
col.r2 <- r1(21) # 20 color steps
col.r2[22] <- "gray50" #21st color grey for non-significance
col.mods <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#ffed6f','gray70')
col.box <- c('#1b9e77','#d95f02','#7570b3')

##########################################################
# functions ####
##########################################################

read.map.from.nc <- function(fn,vn){  
  ncs <- nc_open(fn)
  data <- ncvar_get(ncs,varid=vn)
  nc_close(ncs)
  data
}

plot.mymap1 <- function(screen.ind,sx,ex,lx,sy,ey,ly,mp,ra,lra,cx,cy,co,lco,tx,la,lo){
  screen(screen.ind)
  par(mar=c(0.1,0.1,0.1,0.1),cex=1.2)
  #cat(sx,ex,lx,sy,ey,ly,"\n\n")
  image(x=seq(sx,ex,len=lx),y=seq(sy,ey,len=ly),
        mp,zlim=ra,xlim=cx, ylim=cy, col=co,xlab="",ylab="",axes=F)
  text(-160,-30,tx,adj=c(0,0),pos=4)
  map(boundary=T,add=T,lwd=0.2,col="grey20")
  #plot(world_borders,xlim=cx,ylim=cy,add=TRUE,fg="transparent",lwd=.75,col="grey70")
  box()
}

four.map1 <- function(map1,map2,map3,map4,
                      range1,range2,range3,range4,
                      l.range1,l.range2,l.range3,l.range4,
                      col1,col2,col3,col4,
                      l.col1,l.col2,l.col3,l.col4,
                      lab1="",lab2="",lab3="",lab4="",
                      text1="A)",text2="B)",text3="C)",text4="D)",
                      l1="",l2="",l3="",l4="",
                      startx=-179.75,endx=179.75,lenx=720,starty=-89.75,endy=89.75,leny=360,
                      cxlim=c(-180,180),cylim=c(-55,65)) {
  
  # first divide up the figure region
  split.screen( rbind(c(0,1,0.08,1), c(0,1,0,0.08)))->ind1
  split.screen(rbind(c(0,0.5,0.5,1),c(0.5,1,0.5,1),c(0,0.5,0,0.5),c(0.5,1,0,0.5)),screen=ind1[1]) -> ind2
  # now divide screen into the figure region and legend colorbar on the
  # right to put a legend.
  #   split.screen( rbind(c(0,1,0.15,1), c(0,1,0,0.15)),screen=1)->ind1
  #   split.screen( rbind(c(0,1,0.15,1), c(0,1,0,0.15)),screen=2)->ind2
  #   split.screen( rbind(c(0,1,0.15,1), c(0,1,0,0.15)),screen=3)->ind3
  #   split.screen( rbind(c(0,1,0.15,1), c(0,1,0,0.15)),screen=4)->ind4
  if(max(map1,na.rm=T)>NODATA){
    plot.mymap1(ind2[1],startx,endx,lenx,starty,endy,leny,map1,range1,l.range1,cxlim,cylim,col1,l.col1,text1,lab1,l1)
  }
  if(max(map2,na.rm=T)>NODATA){
    plot.mymap1(ind2[2],startx,endx,lenx,starty,endy,leny,map2,range2,l.range2,cxlim,cylim,col2,l.col2,text2,lab2,l2)
  }
  if(max(map3,na.rm=T)>NODATA){
    plot.mymap1(ind2[3],startx,endx,lenx,starty,endy,leny,map3,range3,l.range3,cxlim,cylim,col3,l.col3,text3,lab3,l3)
  }
  if(max(map4,na.rm=T)>NODATA){
    plot.mymap1(ind2[4],startx,endx,lenx,starty,endy,leny,map4,range4,l.range4,cxlim,cylim,col4,l.col4,text4,lab4,l4)
  }
  screen(ind1[2])
  #plot.mymap1 <- function(screen.ind,sx,ex,lx,sy,ey,ly,mp,ra,lra,cx,cy,co,lco,tx,la,lo){
  
  par(mar=c(1.2,0.1,0.1,0.1),cex=1.5)
  image.plot(zlim=l.range1, col=l.col1,horizontal=T,legend.only=T,legend.width=6,
             legend.args=list( text=lab1, side=4, line=0.3, las=1, cex=1.5 ),
             axis.args=list(lwd=0.5,mgp=c(3, .1, 0),tcl=-0.2),
             smallplot=c(0.05,0.7,0.66,0.95))
  
  close.screen( all=TRUE)
}

##########################################################
# main                                             #######
##########################################################
cl <- clim[1]

for(cc in 1:length(cropsl)){
  # plot 4 maps with one figure on Ray vs. Iizumi correlation
  load(paste0(picture.path,cropsl[cc],"_",tolower(cl),
              "_hist_yield_",cropss[cc],"_ray_vs_iizumi_R2__shifted_ts_processed_dt.Rdata"))
  assign(paste0("rvi",cc),mapr2)
  rm(mapr2)
  ha <- "default"
  load(paste0(picture.path,cropsl[cc],"_best_corrlation_", tolower(cl),
              "_hist_",ha, "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata"))
  assign(paste0("bestr2",cc),mapi)
  rm(mapi)
}

# plot 4 maps with one figure on Ray vs. Iizumi correlation
png(paste0(picture.path,"allcrops_4panel_ray_vs_iizumi_R2.png"),
    height=5*300,width=8*300,res=300,pointsize=6,type="cairo") 
par("oma"=c(0,0,0,0))

four.map1(rvi1[,360:1],rvi2[,360:1],rvi3[,360:1],rvi4[,360:1],
          c(0,1.05),c(0,1.05),c(0,1.05),c(0,1.05),
          c(0,1),c(0,1),c(0,1),c(0,1),
          col.r2,col.r2,col.r2,col.r2,
          col.r2[1:21],col.r2[1:21],col.r2[1:21],col.r2[1:21],
          lab1=expression(paste("coefficient of determination (R"^2,")")))
dev.off()

# 4-panel map with R2 per harmonization setting + the Ray et al. original. ####
for(cc in 1:length(cropsl)){
  # plot Ray R2 values ####
  rayr2 <- read.map.from.nc(paste0(path.ray.2015,"NonCategoricalFigure2",cropsr[cc],"05.nc"),"Data")
  for(ha in 1:length(harms)){
    load(paste0(picture.path,cropsl[cc],"_best_corrlation_", tolower(cl),
              "_hist_",harms[ha], "_yield_",cropss[cc],"_annual_1980_2010_TEST.R_shifted_ts_processed_dt.Rdata"))
    assign(paste0("bestr2",ha),mapi)
    rm(mapi)
  }
  png(paste0(picture.path,cropss[cc],"_4panel_best_per_harm_vs_ray_R2.png"),
      height=5*300,width=8*300,res=300,pointsize=6,type="cairo") 
  par("oma"=c(0,0,0,0))
  
  four.map1(bestr21[,360:1],bestr22[,360:1],bestr23[,360:1],rayr2,
            c(0,1.05),c(0,1.05),c(0,1.05),c(0,1.05),
            c(0,1),c(0,1),c(0,1),c(0,1),
            col.r2,col.r2,col.r2,col.r2,
            col.r2[1:21],col.r2[1:21],col.r2[1:21],col.r2[1:21],
            text1="A) GGCMs default",text2="B) GGCMs fullharm",
            text3="C) GGCMs\n     harm-suffN",text4="D) Ray et al. (2015)",
            lab1=expression(paste("coefficient of determination (R"^2,")")))
  dev.off()
  rm(bestr21,bestr22,bestr23,rayr2)
}

# for individual GGCMs
for(cc in 1:length(cropsl)){
  # plot Ray R2 values ####
  rayr2 <- read.map.from.nc(paste0(path.ray.2015,"NonCategoricalFigure2",cropsr[cc],"05.nc"),"Data")
  for(gg in ggcms){
    # empty arrays in case there is no data for one or several.
    bestr21 <- bestr22 <- bestr23 <- array(NA,dim=c(720,360))
    for(ha in 1:length(harms)){
      fn <- paste0(picture.path,cropsl[cc],"_", tolower(gg),"_", tolower(cl),
                   "_hist_",harms[ha], "_yield_",cropss[cc],"_annual_1980_2010_R2_shifted_ts_processed_dt.Rdata")
      if(!file.exists(fn)) next
      load(fn)
      assign(paste0("bestr2",ha),det.r_sim)
      rm(det.r_sim)
    }
    png(paste0(picture.path,cropss[cc],"_4panel_best_per_harm_vs_ray_R2_",gg,".png"),
        height=5*300,width=8*300,res=300,pointsize=6,type="cairo") 
    par("oma"=c(0,0,0,0))
    
    four.map1(bestr21[,360:1],bestr22[,360:1],bestr23[,360:1],rayr2,
              c(0,1.05),c(0,1.05),c(0,1.05),c(0,1.05),
              c(0,1),c(0,1),c(0,1),c(0,1),
              col.r2,col.r2,col.r2,col.r2,
              col.r2[1:21],col.r2[1:21],col.r2[1:21],col.r2[1:21],
              text1=paste0("A) ",gg,"\n     default"),text2=paste0("B) ",gg,"\n     fullharm"),
              text3=paste0("C) ",gg,"\n     harm-suffN"),text4="D) Ray et al. (2015)",
              lab1=expression(paste("coefficient of determination (R"^2,")")))
    dev.off()
  }
  rm(bestr21,bestr22,bestr23,rayr2)
}
