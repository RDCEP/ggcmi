require(lubridate)
require(ncdf4)
rm(list=ls(all=TRUE))

midway <- F
cluster <- T

pix <- matrix(c(-112.25,33.25,-109.75,27.25,-99.25,19.25,-51.25,-23.25,32.75,24.25,
         33.75,14.25,75.25,15.25,88.75,25.75,5.75,51.75,-58.25,-37.75,
         75.75,30.75,75.75,22.75,89.25,43.25,-96.75,39.25,0.25,51.75,
         3.25,49.75,1.75,47.75,9.75,54.75,118.25,32.25,114.25,37.75,126.25,45.25,
         117.25,-33.75,146.25,-34.25,50.75,35.75,73.25,31.25,72.75,50.25,
         38.75,45.25,33.25,49.25,27.25,38.75,-112.75,49.75),nrow=30,byrow=T)

index <- cbind(as.integer((pix[,1]+180)/0.5+1.01),as.integer((pix[,2]+90)/0.5+1.01))


syear <- c(1950,1951,1961,1971,1981,1991,2001,2005,2011,2021,2031,2041,2051,2061,2071,2081,2091)
eyear <- c("","-1960","-1970","-1980","-1990","-2000","-2004","-2010","-2020","-2030","-2040","-2050","-2060","-2070","-2080","-2090","-2099")
years <- c(1950:2099)
ly <- 37#length(which(leap_year(years)))
# cleaned files on midway have 4 dims
extract.data <- function(fn,ind){ # filename, mask
  nc <- nc_open(fn)
  nn <- names(nc$var)
  data <- ncvar_get(nc,varid=nn[1]) # get first (and only) variable
  nc_close(nc)
  # don't know a more elegant way, so do it cumbersomely
  dat <- array(0,dim=c(dim(ind)[1],dim(data)[3:4]))
  for(i in 1:dim(ind)[1]){
    dat[i,,] <- data[ind[i,1],360-ind[i,2],,] # invert latitudes
  }
  dat
}
extract.T.data <- function(fn,ind,rcp){ # filename, mask
  
  dat <- array(0,dim=c(dim(ind)[1],365*length(years)+ly))
  pos <- 1
  for(y in 1:length(syear)){
    fn2 <- paste(fn,if(syear[y]<2003) "historical" else rcp,"_",syear[y],eyear[y],".nc4",sep="")
    nc <- nc_open(fn2)
    nn <- names(nc$var)
    cat("adding to", pos,"\n")
    for(i in 1:dim(ind)[1]){
      # get one point at the time (otherwise too memory demanding)
      data <- ncvar_get(nc,varid=nn[1],start=c(ind[i,1],360-ind[i,2],1),count=c(1,1,-1)) # length of 3rd dim of data is variable
      dat[i,c(pos:(pos+length(data)-1))] <- data
    }
    nc_close(nc)
    pos <- pos+length(data) # where to glue the new data...
  }
  dat
}

ggcms <- c("EPIC","GEPIC","IMAGE_LEITAP","LPJ-GUESS","LPJmL","pDSSAT","PEGASUS")
gcms <- "HadGEM2-ES"
rcps <- c("rcp2p6","rcp4p5","rcp6p0","rcp8p5")[c(1,4)] # currently only rcp2p6 and rcp8p5 are available in the cleaned version
irrig <- c("noirr","firr")

everything <- list()

if(midway){
  # this works on midway.rcc.uchicago.edu
  path <- "/project/ggcmi/isi1/isi1.clean/"
  for(ggcm in ggcms){
    for(gcm in gcms){
      for(rcp in rcps){
        fn <- paste(path,ggcm,"/",gcm,"/wheat/",rcp,"/noco2/",tolower(ggcm),"_",tolower(gcm),"_ssp2_noco2_yield_whe_annual_1980_2099.nc4",sep="")
        if(!is.na(file.info(fn)$size)){
          data <- extract.data(fn,index)
          for(ir in 1:length(irrig)){
            everything[[paste(ggcm,gcm,rcp,irrig[ir],"noco2",sep="_")]] <- data[,,ir]
          }        
        }
      }
    }
  }
  
  save(everything,ggcms,gcms,rcps,irrig,file="/home/chmueller/GGCMI/extracted.ggcmi.fasttrack.wheat.30pix.Rdata")
  
}

if(cluster){
  # this works on the PIK cluster
  path <- "/iplex/01/2011/isimip/inputdata_bced/HadGEM2-ES/"
  
  everything <- list()
  
  for(gcm in gcms){
    for(rcp in rcps){
      fn <- paste(path,"tas_bced_1960_1999_hadgem2-es_",sep="")
      # extract.T.data is pretty slow with all the loops and point-wise extraction from so many files.
      # But I ran into memory constraints when doing more wasteful data extraction from NC files, so I chose this...
      data <- extract.T.data(fn,index,rcp)
      everything[[paste(gcm,rcp,sep="_")]] <- data
    }
  }
  save(list = ls(all = TRUE), file = "extract.RData")
  save(everything,gcms,file="/iplex/01/2014/macmit/users/cmueller/GGCMI/extracted.daily.Tas.30pix.Rdata")
  
}
# #### testing at PIK
# fn <- "D:/data/GGCMI/epic_hadgem2-es_ssp2_noco2_yield_whe_annual_1980_2099.nc4"
# image.plot(x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),data[,,1,1])
# points(pix2,col=2)
# pix2 <- pix
# pix2[,2] <- -pix[,2]
