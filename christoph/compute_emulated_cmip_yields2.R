require(ncdf4)
require(fields)
require(reshape)
require(abind)
require(vioplot)

path.clim5 <- "/p/projects/macmit/data/CMIP/CMIP5/pcmdi_crugrid/"
path.rclim5 <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/CMIP/CMIP5/"
path.clim6 <- "/p/projects/macmit/data/CMIP/CMIP6/crugrid/"
path.rclim6 <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/CMIP/CMIP6/"
path.gs <- "/p/projects/macmit/data/GGCMI/AgMIP.input/other.inputs/"
path.lu <- "/p/projects/macmit/data/GGCMI/AgMIP.output/processed/masks/weight/"
path.prod <- "/p/projects/lpjml/reference_data/yield/monfreda_2008/Rdata/05deg/"
path.fert <- "/p/projects/macmit/data/GGCMI/AgMIP.input/other.inputs/AGMIP_NUTRIENTS.HARM.version1.1/"
path.coeff <- "/p/projects/macmit/data/GGCMI/AgMIP.output/Jim_Emulator/EMULATOR_PARAMS/"
path.ryield5 <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/CMIP/CMIP5/yields/"
path.ryield6 <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/CMIP/CMIP6/yields/"
path.riwd5 <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/CMIP/CMIP5/iwd/"
path.rprod5 <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/CMIP/CMIP5/production/"
path.riwd6 <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/CMIP/CMIP6/iwd/"
path.rprod6 <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/CMIP/CMIP6/production/"
path.figs <- "/p/projects/macmit/users/cmueller/GGCMI/phase2/emulate_CMIP/figures/"

crops.gs <- c("Maize","Rice","Soybeans","wwh","swh")
crops.gs2 <- c("mai","ric","soy","wwh","swh")
crops.nice <- c("Maize","Rice","Soybeans","Winter Wheat","Spring Wheat","all combined")
crops.param <- c("maize","rice","soy","winter_wheat","spring_wheat")
crops.param2 <- c("maize","rice","soy","winter_wheat","spring_wheat","all")
crops.fert <- c("maize","rice","soybean","wheat","wheat")
ndaymonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
ggcms <- c("CARAIB","EPIC-TAMU","GEPIC","JULES","LPJ-GUESS","LPJmL","PEPIC","PROMET","pDSSAT")
gcms5 <- c("ACCESS1-0","ACCESS1-3","BCC-CSM1-1","BCC-CSM1-1-M","BNU-ESM","CanESM2","CCSM4","CESM1-BGC","CESM1-CAM5","CESM1-CAM5-1-FV2","CESM1-WACCM",
           "CMCC-CESM","CMCC-CM","CMCC-CMS","CNRM-CM5","CSIRO-Mk3-6-0","EC-EARTH","FGOALS-g2","FIO-ESM","GFDL-CM3","GFDL-ESM2G","GFDL-ESM2M",
           "GISS-E2-H","GISS-E2-H","GISS-E2-H","GISS-E2-H-CC","GISS-E2-R","GISS-E2-R","GISS-E2-R","GISS-E2-R-CC","HadGEM2-AO","HadGEM2-CC","HadGEM2-ES","INMCM4",
           "IPSL-CM5A-LR","IPSL-CM5A-MR","IPSL-CM5B-LR","MIROC-ESM","MIROC-ESM-CHEM","MIROC5","MPI-ESM-LR","MPI-ESM-MR","MRI-CGCM3","MRI-ESM1",
           "NorESM1-M","NorESM1-ME")
gcms5 <- gcms5[-which(gcms5=="CESM1-CAM5-1-FV2")] # misses precip data for Dec 2056
gcms6 <- c("ACCESS-CM2","ACCESS-ESM1-5","BCC-CSM2-MR","CAMS-CSM1-0","CanESM5","CanESM5","CanESM5-CanOE","CESM2","CESM2-WACCM","CIESM","CNRM-CM6-1","CNRM-ESM2-1","EC-Earth3",
           "EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","FIO-ESM-2-0","GFDL-CM4","GFDL-ESM4","GISS-E2-1-G","HadGEM3-GC31-LL","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MCM-UA-1-0",
           "MIROC6","MIROC-ES2L","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NESM3","NorESM2-LM","NorESM2-MM","UKESM1-0-LL")
run5 <- rep(1,length(gcms5))
run5[which(gcms5=="CESM1-WACCM")] <- 2
run5[which(gcms5=="HadGEM2-ES")] <- 2
run5[which(gcms5=="EC-EARTH")] <- 8
run6 <- rep(1,length(gcms6))
param5 <- rep(1,length(gcms5))
param5[which(gcms5=="GISS-E2-H")[2]] <- 2
param5[which(gcms5=="GISS-E2-H")[3]] <- 3
param5[which(gcms5=="GISS-E2-R")[2]] <- 2
param5[which(gcms5=="GISS-E2-R")[3]] <- 3
param6 <- rep(1,length(gcms6))
param6[which(gcms6=="CanESM5")[2]] <- 2
param6[which(gcms6=="CanESM5-CanOE")] <- 2
param6[which(gcms6=="GISS-E2-1-G")] <- 3
file6 <- rep(1,length(gcms6))
file6[which(gcms6=="CNRM-CM6-1")] <- 2
file6[which(gcms6=="CNRM-ESM2-1")] <- 2
file6[which(gcms6=="UKESM1-0-LL")] <- 2
file6[which(gcms6=="HadGEM3-GC31-LL")] <- 3
file6[which(gcms6=="MCM-UA-1-0")] <- 2
file6[which(gcms6=="MIROC-ES2L")] <- 2
sy5 <- rep(1950,length(gcms5))
sy5[which(gcms5=="CESM1-WACCM")] <- 1955 # CESM1-WACCM only starts in 1955 -- despite the file name stating 1950...

rcps <- c("rcp26","rcp45","rcp60","rcp85")
rcps.nr <- c("2.6","4.5","6.0","8.5")
rcps.nice <- c("RCP 2.6","RCP 4.5","RCP 6.0","RCP 8.5")

ssps <- c("ssp119","ssp126","ssp245","ssp370","ssp585")[2:5]
rssps <- c("rcp19","rcp26","rcp45","rcp70","rcp85")[2:5]
ssps.nr <- c("1.9","2.6","4.5","7.0","8.5")[2:5]
ssps.nice <- c("SSP1 RCP 1.9","SSP1 RCP 2.6","SSP2 RCP 4.5","SSP3 RCP 7.0","SSP5 RCP 8.5")[2:5]

#c("whe","soy","ric","mai","mgr","mil","pea","nut","cas","rap","sgb","sug","sun")
#ENERG <- c(3340,3350,2800,3560,1,3400,3410,4140,1090,4940,700,300,3080) # kcal/kg fresh matter, dummy 1 for managed grass
#FRESHMATTER <- 100 / c(88, 91, 87, 88,100,88,90,94,35,92,24,27,93) #dummy 1 for managed grass
# maize, rice, soy, wheat, wheat
ENERG <- c(3560,2800,3350,3340,3340)# kcal/kg fresh matter
ENERG <- ENERG*1e3 #kcal/tFM
FRESHMATTER <- 100/c(88,87,91,88,88)
ENERG_DM <- ENERG*FRESHMATTER

# colors ####
col.var <- c(rgb(133,200,182,maxColorValue=255),rgb(19,134,178,maxColorValue=255),rgb(70,85,107,maxColorValue=255))
col.rcp <- c(rgb(239,209,0,maxColorValue=255),rgb(203,88,27,maxColorValue=255),rgb(139,50,57,maxColorValue=255),rgb(65,12,83,maxColorValue=255))
col.rcp2 <- c(rgb(239,209,0,45,maxColorValue=255),rgb(203,88,27,45,maxColorValue=255),rgb(139,50,57,45,maxColorValue=255),rgb(65,12,83,45,maxColorValue=255))
col.ggcms <- c(rgb(141,211,199,maxColorValue=255),rgb(255,255,179,maxColorValue=255),rgb(190,186,218,maxColorValue=255),rgb(251,128,114,maxColorValue=255),rgb(128,177,211,maxColorValue=255),rgb(253,180,98,maxColorValue=255),rgb(179,222,105,maxColorValue=255),rgb(252,205,229,maxColorValue=255),rgb(217,217,217,maxColorValue=255),rgb(188,128,189,maxColorValue=255))
col.ggcms2 <- c(rgb(141,211,199,25,maxColorValue=255),rgb(255,255,179,25,maxColorValue=255),rgb(190,186,218,25,maxColorValue=255),rgb(251,128,114,25,maxColorValue=255),rgb(128,177,211,25,maxColorValue=255),rgb(253,180,98,25,maxColorValue=255),rgb(179,222,105,25,maxColorValue=255),rgb(252,205,229,25,maxColorValue=255),rgb(217,217,217,25,maxColorValue=255),rgb(188,128,189,25,maxColorValue=255))



cr <- 2
ggcm <- 2
gcm <- 23
rcp <- 4
adap <- 0

# functions ####

# function to read NC map
readmap.nc <- function(filename,var="",lo="lon",la="lat",starttime=1){
  nc <- nc_open(filename)
  if(var=="") var <- names(nc$var)[1]
  lon <- ncvar_get(nc,lo)
  if(min(lon)>=0){
    cat("WARNING! Longitude does not contain negative values, shifting >180 by 360\n")
    lon[lon>180] <- lon[lon>180]-360
  }
  lat <- ncvar_get(nc,la)
  if(lat[1] > lat[length(lat)]){
    cat("WARNING, inverting latitudes\n")
    
  }
  if(starttime==1)
    buf <- ncvar_get(nc,var)
  else
    buf <-ncvar_get(nc,var,start=c(1,1,starttime))
  nc_close(nc)
  if(length(dim(buf))==2)
    buf <- buf[order(lon),order(lat)]
  else if(length(dim(buf))==3)
    buf <- buf[order(lon),order(lat),]
  else if(length(dim(buf))>3)
    cat("WARNING, cannot adjust lon/lat setting for 4-dim array\n")
  buf
}

writemap.nc <- function(fn,ldata,var="yield",var2=var,crop="",cropf="",
                        units="tDM/ha",mv=1e20,start.year=1980,end.year=2099,
                        title="emulator based CMIP5 yield projections",comment2=""){
  
  # NC file dimensions
  dim_lon <- ncdim_def("lon","degrees_east",seq(-179.75,179.75,len=360/0.5))
  #change order of latitudes
  dim_lat <- ncdim_def("lat","degrees_north",seq(89.75,-89.75,len=180/0.5))
  dim_time <- ncdim_def("time",paste("years since ",start.year-1,"-01-01 00:00:00",sep=""),
                        c(start.year:end.year)-start.year+1,
                        calendar = "standard")
  
  lncv <- list()
  for(i in 1:length(var)){
    # define variable
    data <- ldata[[i]]
    if(length(dim(data))==2){
      ncv <- ncvar_def(paste(var[i],"_",crop,sep=""),units,list(dim_lon,dim_lat),mv,
                       longname=paste(cropf,var2[i]),
                       compression=9)
    } else if(length(dim(data))==3){
      ncv <- ncvar_def(paste(var[i],"_",crop,sep=""),units,list(dim_lon,dim_lat,dim_time),mv,
                       longname=paste(cropf,var2[i]),
                       compression=9)
    } else {
      return(cat("wrong number of dimensions in data structure",dim(data),"\n"))
    }
    lncv[[i]] <- ncv
  }
  # create files
  #ncf <- nc_create(fn,list(ncv))
  ncf <- nc_create(fn,lncv)
  
  # commenting
  if(""!=title) ncatt_put(ncf,varid=0,"title",title)
  ncatt_put(ncf,varid=0,"experiment","AgMIP GGCMI phase 2")
  if(""!=comment2) ncatt_put(ncf,varid=0,"method",comment2)
  ncatt_put(ncf,varid=0,"author","cmueller@pik-potsdam.de")
  for(i in 1:length(ldata)){
    data <- ldata[[i]]
    data[!is.finite(data)] <- mv
    # checking data for NC files
    if(dim(data)[1]!=720)
      return(cat("wrong length of longitude dimension (1) in data structure",dim(data)[1],"\n"))
    if(dim(data)[2]!=360)
      return(cat("wrong length of latitude dimension (2) in data structure",dim(data)[2],"\n"))
    # needs 3D matrix: lon, lat, time
    if(length(dim(data))==2){
      ncvar_put(ncf,lncv[[i]],data,start=c(1,1),count=c(-1,-1))
    } else if(length(dim(data))==3){
      # writing data to NC files looping through time
      for(j in 1:dim(data)[3]){
        ncvar_put(ncf,lncv[[i]],data[,,j],start=c(1,1,j),count=c(-1,-1,1))
      }
    } 
  }
  # closing NC files
  nc_close(ncf)
}
# taken from LPJmL code
cvrtdaymonth <- function(doy){
  if(doy<0)
    return(NA)
  sum <- 0
  for(i in 1:12)
    if(doy<=sum+ndaymonth[i])
    {
      month <- i
      dayofmonth <- doy-sum
      break;
    }
  else
    sum <- sum + ndaymonth[i];
  c(month,dayofmonth)
}
compute.weights <- function(pdp,hdp){
  ved <- rep(0,12)
  #pdp <- li[[1]]
  #hdp <- li[[2]]
  pds <- cvrtdaymonth(pdp)
  hds <- cvrtdaymonth(hdp)
  if(is.na(pds[1]) | is.na(hds[1])){
    return(ved)
  } else if(pds[1]>hds[1]){ # planting in fall, harvest next year
    ved[pds[1]:12] <- ndaymonth[pds[1]:12]
    ved[1:hds[1]] <- ndaymonth[1:hds[1]]
  } else {
    ved[pds[1]:hds[1]] <- ndaymonth[pds[1]:hds[1]]
  }
  ved[pds[1]] <- ndaymonth[hds[1]]-pds[2]+1
  ved[hds[1]] <- hds[2]
  ved/sum(ved)
}

# emulator

emulator <- function(cf,C,T,W,N,irrig=F){
  if(dim(cf)[3] == 34 | (irrig & dim(cf)[3]==19)){
    #cat("34 coefficients found\n")
    if(irrig){
      Y = (cf[,,1] + cf[,,2]*C + cf[,,3]*T + cf[,,4]*N + cf[,,5]*C^2 + cf[,,6]*C*T + cf[,,7]*C*N + cf[,,8]*T^2 + 
             cf[,,9]*T*N + cf[,,10]*N^2 + cf[,,11]*C^3 + cf[,,12]*C^2*T + cf[,,13]*C^2*N + cf[,,14]*C*T^2 + 
             cf[,,15]*C*T*N + cf[,,16]*C*N^2 + cf[,,17]*T^3 + cf[,,18]*T^2*N + cf[,,19]*T*N^2)
    } else {
      Y <- (cf[,,1] + cf[,,2]*C + cf[,,3]*T + cf[,,4]*W + cf[,,5]*N + cf[,,6]*C^2 + cf[,,7]*C*T + cf[,,8]*C*W + 
              cf[,,9]*C*N + cf[,,10]*T^2 + cf[,,11]*T*W + cf[,,12]*T*N + cf[,,13]*W^2 + cf[,,14]*W*N + cf[,,15]*N^2 + 
              cf[,,16]*C^3 + cf[,,17]*C^2*T + cf[,,18]*C^2*W + cf[,,19]*C^2*N + cf[,,20]*C*T^2 + cf[,,21]*C*T*W + 
              cf[,,22]*C*T*N + cf[,,23]*C*W^2 + cf[,,24]*C*W*N + cf[,,25]*C*N^2 + cf[,,26]*T^3 + cf[,,27]*T^2*W + 
              cf[,,28]*T^2*N + cf[,,29]*T*W^2 + cf[,,30]*T*W*N + cf[,,31]*T*N^2 + cf[,,32]*W^3 + 
              cf[,,33]*W^2*N + cf[,,34]*W*N^2)
    }
  } else if(dim(cf)[3]==20 | (irrig & dim(cf)[3]==10)){
    #cat("20 coefficients found\n")
    if(irrig){
      Y = (cf[,,1] + cf[,,2]*C + cf[,,3]*T + cf[,,4]*C^2 + cf[,,5]*C*T + cf[,,6]*T^2 + cf[,,7]*C^3 + 
             cf[,,8]*C^2*T + cf[,,9]*C*T^2 + cf[,,10]*T^3)
    } else{
      Y <- (cf[,,1] + cf[,,2]*C + cf[,,3]*T + cf[,,4]*W + cf[,,5]*C^2 + cf[,,6]*C*T + cf[,,7]*C*W + cf[,,8]*T^2 +
              cf[,,9]*T*W + cf[,,10]*W^2 + cf[,,11]*C^3 + cf[,,12]*C^2*T + cf[,,13]*C^2*W + cf[,,14]*C*T^2 + 
              cf[,,15]*C*T*W + cf[,,16]*C*W^2 + cf[,,17]*T^3 + cf[,,18]*T^2*W + cf[,,19]*T*W^2 + cf[,,20]*W^3)
    }
  }
  Y
}
emulator.iwd <- function(cf,C,T,N){
  if(dim(cf)[3]==19){
    # Y = (cf[,,1] + cf[,,2]*C + cf[,,3]*T + cf[,,4]*N + cf[,,5]*C^2 + cf[,,6]*T^2 + cf[,,7]*N^2 + cf[,,8]*C*N + 
    #        cf[,,9]*T*N + cf[,,10]*C*T + cf[,,11]*T^3 + cf[,,12]*C^3 + cf[,,13]*T*N + cf[,,14]*C*N^2 + 
    #        cf[,,15]*N^2*T + cf[,,16]*T^2*N + cf[,,17]*T^2*C + cf[,,18]*C^2*T + cf[,,19]*C^2*N)
    Y = (cf[,,1] + cf[,,2]*C + cf[,,3]*T + cf[,,4]*N + cf[,,5]*C^2 + cf[,,6]*C*T + cf[,,7]*C*N + cf[,,8]*T^2 + 
           cf[,,9]*T*N + cf[,,10]*N^2 + cf[,,11]*C^3 + cf[,,12]*C^2*T + cf[,,13]*C^2*N + cf[,,14]*C*T^2 + 
           cf[,,15]*C*T*N + cf[,,16]*C*N^2 + cf[,,17]*T^3 + cf[,,18]*T^2*N + cf[,,19]*T*N^2)
  } else if(dim(cf)[3]==10){
    #cat("20 coefficients found\n")
    # Y = (cf[,,1] + cf[,,2]*C + cf[,,3]*T + cf[,,4]*C^2 + cf[,,5]*T^2 + 
    #        cf[,,6]*C*T + cf[,,7]*T^3 + cf[,,8]*C^3 +  
    #        cf[,,9]*T^2*C + cf[,,10]*C^2*T)
    Y = (cf[,,1] + cf[,,2]*C + cf[,,3]*T + cf[,,4]*C^2 + cf[,,5]*C*T + cf[,,6]*T^2 + cf[,,7]*C^3 + 
           cf[,,8]*C^2*T + cf[,,9]*C*T^2 + cf[,,10]*T^3)
  }
  Y
}

# helper functions
varall <- function(x){
  vari <- var(as.vector(x),na.rm=T)
  vari
}

replace.with.first.miss <- function(x,thres){
  if(length(x)%%2 != 0){
    cat("ERROR, vector has odd number of elements!\n")
    return(NA)
  }
  br <- as.integer(length(x)/2)
  # first half of vector is values, 2nd half is by which to split
  val <- x[1:br]
  cri <- x[c(1:br)+br]
  repl <- which(cri>=thres)
  if(length(repl)>0){
    val[cri>thres] <- val[repl[1]]
  }
  val
}

# vioplot funciont
vioplot_gs <- function(fn,A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,B6,colA1,colA2,colB1,colB2){
  A1 <- sample(A1,1e5)
  A2 <- sample(A2,1e5)
  A3 <- sample(A3,1e5)
  A4 <- sample(A4,1e5)
  A5 <- sample(A5,1e5)
  A6 <- sample(A6,1e5)
  B1 <- sample(B1,1e5)
  B2 <- sample(B2,1e5)
  B3 <- sample(B3,1e5)
  B4 <- sample(B4,1e5)
  B5 <- sample(B5,1e5)
  B6 <- sample(B6,1e5)
  if(!is.null(fn)) png(fn,width=5*300,height=5*300,res=300,pointsize=10)
  
  vioplot(list(all=A1,maize=1,rice=1,soy=1,swheat=1,wwheat=1),
          horizontal=T,at=c(1:6)-0.1,
          border=c(1,NA,NA,NA,NA,NA),col=c(colA1,NA,NA,NA,NA,NA),rectCol=c(colA2,NA,NA,NA,NA,NA),
          lineCol=c(colA2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
  vioplot(list(all=1,maize=A2,rice=1,soy=1,swheat=1,wwheat=1),
          add=T,horizontal=T,at=c(1:6)-0.1,
          border=c(NA,1,NA,NA,NA,NA),col=c(NA,colA1,NA,NA,NA,NA),rectCol=c(NA,colA2,NA,NA,NA,NA),
          lineCol=c(NA,colA2,NA,NA,NA,NA),colMed=c(NA,1,NA,NA,NA,NA),colMed2=c(NA,1,NA,NA,NA,NA))
  vioplot(list(all=1,maize=1,rice=A3,soy=1,swheat=1,wwheat=1),
          add=T,horizontal=T,at=c(1:6)-0.1,
          border=c(NA,NA,1,NA,NA,NA),col=c(NA,NA,colA1,NA,NA,NA),rectCol=c(NA,NA,colA2,NA,NA,NA),
          lineCol=c(NA,NA,colA2,NA,NA,NA),colMed=c(NA,NA,1,NA,NA,NA),colMed2=c(NA,NA,1,NA,NA,NA))
  vioplot(list(all=1,maize=1,rice=1,soy=A4,swheat=1,wwheat=1),
          add=T,horizontal=T,at=c(1:6)-0.1,
          border=c(NA,NA,NA,1,NA,NA),col=c(NA,NA,NA,colA1,NA,NA),rectCol=c(NA,NA,NA,colA2,NA,NA),
          lineCol=c(NA,NA,NA,colA2,NA,NA),colMed=c(NA,NA,NA,1,NA,NA),colMed2=c(NA,NA,NA,1,NA,NA))
  vioplot(list(all=1,maize=1,rice=1,soy=1,swheat=A5,wwheat=1),
          add=T,horizontal=T,at=c(1:6)-0.1,
          border=c(NA,NA,NA,NA,1,NA),col=c(NA,NA,NA,NA,colA1,NA),rectCol=c(NA,NA,NA,NA,colA2,NA),
          lineCol=c(NA,NA,NA,NA,colA2,NA),colMed=c(NA,NA,NA,NA,1,NA),colMed2=c(NA,NA,NA,NA,1,NA))
  vioplot(list(all=1,maize=1,rice=1,soy=1,swheat=1,wwheat=A6),
          add=T,horizontal=T,at=c(1:6)-0.1,
          border=c(NA,NA,NA,NA,NA,1),col=c(NA,NA,NA,NA,NA,colA1),rectCol=c(NA,NA,NA,NA,NA,colA2),
          lineCol=c(NA,NA,NA,NA,NA,colA2),colMed=c(NA,NA,NA,NA,NA,1),colMed2=c(NA,NA,NA,NA,NA,1))
  
  vioplot(list(all=B1,maize=1,rice=1,soy=1,swheat=1,wwheat=1),
          add=T,horizontal=T,at=c(1:6)+0.1,
          border=c(1,NA,NA,NA,NA,NA),col=c(colB1,NA,NA,NA,NA,NA),rectCol=c(colB2,NA,NA,NA,NA,NA),
          lineCol=c(colB2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
  vioplot(list(all=1,maize=B2,rice=1,soy=1,swheat=1,wwheat=1),
          add=T,horizontal=T,at=c(1:6)+0.1,
          border=c(NA,1,NA,NA,NA,NA),col=c(NA,colB1,NA,NA,NA,NA),rectCol=c(NA,colB2,NA,NA,NA,NA),
          lineCol=c(NA,colB2,NA,NA,NA,NA),colMed=c(NA,1,NA,NA,NA,NA),colMed2=c(NA,1,NA,NA,NA,NA))
  vioplot(list(all=1,maize=1,rice=B3,soy=1,swheat=1,wwheat=1),
          add=T,horizontal=T,at=c(1:6)+0.1,
          border=c(NA,NA,1,NA,NA,NA),col=c(NA,NA,colB1,NA,NA,NA),rectCol=c(NA,NA,colB2,NA,NA,NA),
          lineCol=c(NA,NA,colB2,NA,NA,NA),colMed=c(NA,NA,1,NA,NA,NA),colMed2=c(NA,NA,1,NA,NA,NA))
  vioplot(list(all=1,maize=1,rice=1,soy=B4,swheat=1,wwheat=1),
          add=T,horizontal=T,at=c(1:6)+0.1,
          border=c(NA,NA,NA,1,NA,NA),col=c(NA,NA,NA,colB1,NA,NA),rectCol=c(NA,NA,NA,colB2,NA,NA),
          lineCol=c(NA,NA,NA,colB2,NA,NA),colMed=c(NA,NA,NA,1,NA,NA),colMed2=c(NA,NA,NA,1,NA,NA))
  vioplot(list(all=1,maize=1,rice=1,soy=1,swheat=B5,wwheat=1),
          add=T,horizontal=T,at=c(1:6)+0.1,
          border=c(NA,NA,NA,NA,1,NA),col=c(NA,NA,NA,NA,colB1,NA),rectCol=c(NA,NA,NA,NA,colB2,NA),
          lineCol=c(NA,NA,NA,NA,colB2,NA),colMed=c(NA,NA,NA,NA,1,NA),colMed2=c(NA,NA,NA,NA,1,NA))
  vioplot(list(all=1,maize=1,rice=1,soy=1,swheat=1,wwheat=B6),
          add=T,horizontal=T,at=c(1:6)+0.1,
          border=c(NA,NA,NA,NA,NA,1),col=c(NA,NA,NA,NA,NA,colB1),rectCol=c(NA,NA,NA,NA,NA,colB2),
          lineCol=c(NA,NA,NA,NA,NA,colB2),colMed=c(NA,NA,NA,NA,NA,1),colMed2=c(NA,NA,NA,NA,NA,1))
  
  
  dev.off()
}
vioplot_gs2 <- function(fn,A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,B6,colA1,colA2,colB1,colB2,reduce=F){
  if(reduce){
    A1 <- sample(A1,1e5)
    A2 <- sample(A2,1e5)
    A3 <- sample(A3,1e5)
    A4 <- sample(A4,1e5)
    A5 <- sample(A5,1e5)
    A6 <- sample(A6,1e5)
    B1 <- sample(B1,1e5)
    B2 <- sample(B2,1e5)
    B3 <- sample(B3,1e5)
    B4 <- sample(B4,1e5)
    B5 <- sample(B5,1e5)
    B6 <- sample(B6,1e5)
  }
  if(!is.null(fn)) png(fn,width=5*300,height=5*300,res=300,pointsize=10)
  
  vioplot(list(all=A1,maize=A2,rice=A3,soy=A4,swheat=A5,wwheat=A6),
          horizontal=T,at=c(1:6)-0.1,
          #border=c(1,NA,NA,NA,NA,NA),col=c(colA1,NA,NA,NA,NA,NA),rectCol=c(colA2,NA,NA,NA,NA,NA),
          #lineCol=c(colA2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
          col=colA1)
  
  vioplot(list(all=B1,maize=B2,rice=B3,soy=B4,swheat=B5,wwheat=B6),
          add=T,horizontal=T,at=c(1:6)+0.1,
          #border=c(1,NA,NA,NA,NA,NA),col=c(colB1,NA,NA,NA,NA,NA),rectCol=c(colB2,NA,NA,NA,NA,NA),
          #lineCol=c(colB2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
          col=colB1,rectCol=colB2)
  
  if(!is.null(fn)) dev.off()
}

boxplot_gs <- function(fn=NULL,A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,B6,colA1,colB1,reduce=T,diff=F){
  if(reduce){
    A1 <- sample(A1,1e5)
    A2 <- sample(A2,1e5)
    A3 <- sample(A3,1e5)
    A4 <- sample(A4,1e5)
    A5 <- sample(A5,1e5)
    A6 <- sample(A6,1e5)
    B1 <- sample(B1,1e5)
    B2 <- sample(B2,1e5)
    B3 <- sample(B3,1e5)
    B4 <- sample(B4,1e5)
    B5 <- sample(B5,1e5)
    B6 <- sample(B6,1e5)
  }
  if(!is.null(fn))png(fn,width=5*300,height=5*300,res=300,pointsize=10)
  if(!diff){
    boxplot(list(A1,A2,A3,A4,A5,A6),xlim=c(0.5,7),
            horizontal=T,at=c(1:6)-0.1,axes=F,outline=F,
            #border=c(1,NA,NA,NA,NA,NA),col=c(colA1,NA,NA,NA,NA,NA),rectCol=c(colA2,NA,NA,NA,NA,NA),
            #lineCol=c(colA2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
            col=colA1)
    
    boxplot(list(B1,B2,B3,B4,B5,B6),
            add=T,horizontal=T,at=c(1:6)+0.1,axes=F,outline=F,
            #border=c(1,NA,NA,NA,NA,NA),col=c(colB1,NA,NA,NA,NA,NA),rectCol=c(colB2,NA,NA,NA,NA,NA),
            #lineCol=c(colB2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
            col=colB1)
    legend("topleft",legend=c("CMIP5","CMIP6"),fill=c(colA1,colB1),ncol=2,bty="n")
    axis(1)
    axis(2,at=c(1:6),labels=c("all","maize","rice","soybean","w.wheat","s.wheat"),las=2)
    box()
  } else {
    boxplot(list(B1-A1,B2-A2,B3-A3,B4-A4,B5-A5,B6-A6),xlim=c(0.5,7),
            horizontal=T,at=c(1:6)-0.1,axes=F,outline=F,
            #border=c(1,NA,NA,NA,NA,NA),col=c(colA1,NA,NA,NA,NA,NA),rectCol=c(colA2,NA,NA,NA,NA,NA),
            #lineCol=c(colA2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
            col=colB1)
    axis(1)
    axis(2,at=c(1:6),labels=c("all","maize","rice","soybean","w.wheat","s.wheat"),las=2)
    box()
  }
  if(!is.null(fn))dev.off()
}
boxplot_gs2 <- function(fn=NULL,A2,A3,A4,A5,A6,B2,B3,B4,B5,B6,colA1,colB1,reduce=T,diff=F,pdf=F,xl=""){
  if(reduce){
    A1 <- sample(A1,1e5)
    A2 <- sample(A2,1e5)
    A3 <- sample(A3,1e5)
    A4 <- sample(A4,1e5)
    A5 <- sample(A5,1e5)
    A6 <- sample(A6,1e5)
    B1 <- sample(B1,1e5)
    B2 <- sample(B2,1e5)
    B3 <- sample(B3,1e5)
    B4 <- sample(B4,1e5)
    B5 <- sample(B5,1e5)
    B6 <- sample(B6,1e5)
  }
  A1 <- c(A2,A3,A4,A5,A6)
  B1 <- c(B2,B3,B4,B5,B6)
  if(!is.null(fn)){
    if(pdf) {
      pdf(paste0(fn,".pdf"))
    } else{
      png(paste0(fn,".png"),width=5*300,height=5*300,res=300,pointsize=10)
    }
  }
  if(!diff){
    boxplot(list(A1,A2,A3,A4,A5,A6),xlim=c(0.5,7),
            horizontal=T,at=c(1:6)-0.1,axes=F,outline=F,
            #border=c(1,NA,NA,NA,NA,NA),col=c(colA1,NA,NA,NA,NA,NA),rectCol=c(colA2,NA,NA,NA,NA,NA),
            #lineCol=c(colA2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
            col=colA1,xlab=xl)
    
    boxplot(list(B1,B2,B3,B4,B5,B6),
            add=T,horizontal=T,at=c(1:6)+0.1,axes=F,outline=F,
            #border=c(1,NA,NA,NA,NA,NA),col=c(colB1,NA,NA,NA,NA,NA),rectCol=c(colB2,NA,NA,NA,NA,NA),
            #lineCol=c(colB2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
            col=colB1)
    legend("topleft",legend=c("CMIP5","CMIP6"),fill=c(colA1,colB1),ncol=2,bty="n")
    axis(1)
    axis(2,at=c(1:6),labels=c("all","maize","rice","soybean","w.wheat","s.wheat"),las=2)
    box()
  } else {
    boxplot(list(B1-A1,B2-A2,B3-A3,B4-A4,B5-A5,B6-A6),xlim=c(0.5,7),
            horizontal=T,at=c(1:6)-0.1,axes=F,outline=F,
            #border=c(1,NA,NA,NA,NA,NA),col=c(colA1,NA,NA,NA,NA,NA),rectCol=c(colA2,NA,NA,NA,NA,NA),
            #lineCol=c(colA2,NA,NA,NA,NA,NA),colMed=c(1,NA,NA,NA,NA,NA),colMed2=c(1,NA,NA,NA,NA,NA))
            col=colB1,xlab=xl)
    axis(1)
    axis(2,at=c(1:6),labels=c("all","maize","rice","soybean","w.wheat","s.wheat"),las=2)
    box()
  }
  if(!is.null(fn))dev.off()
}
# parallel setup ####
parallel <- TRUE

#parallel <- F

if(parallel==T) {
  library(Rmpi)  # R implementation of MPI interface
  library(doMPI) # interface for foreach construct to run in MPI parallel mode
  clu <- startMPIcluster(verbose=F) # start cluster (link R instances together)
  num.cluster <- clusterSize(clu)
  if (num.cluster > 1) {
    # we are using more than 1 CPU, so really run in parallel mode
    registerDoMPI(clu) # tells foreach to use MPI parallel mode
    print(paste("Running in parallel mode on",num.cluster,"worker nodes."))
  } else {
    registerDoSEQ() # tells foreach to use sequential mode
    print("Running in sequential mode.")
  }
} else {
  library(foreach)
  registerDoSEQ() # tells foreach to use sequential mode
  print("Running in sequential mode.")
} 

# select CMIP
do.cmip5 <- F
do.cmip6 <- T

version <- "2.0"

# read growing season data and generate monthly weighted deltas ####

# this does not work in parallel mode for some unknown reason
if(F){
  raw_deltas <- T
  capped_deltas <- T
  for(cr in 1:length(crops.gs)){
    fn <- paste0(path.gs,"AGMIP_GROWING_SEASON.HARM.version",if(crops.gs[cr] %in% c("swh","wwh")) "2.0" else "1.25","/",crops.gs[cr],"_rf_growing_season_dates_v",
                 if(crops.gs[cr] %in% c("swh","wwh")) "2" else "1.25",".nc4")
    pd <- readmap.nc(fn,"planting day")#[,360:1] # inverted latitudes compared to emulator coefficients
    hd <- readmap.nc(fn,"harvest day")#[,360:1]
    
    index1 <- which(!is.na(pd) & pd>0 & hd>0,arr.ind = TRUE)
    weights <- array(NA,dim=c(720,360,12))
    
    for(i in 1:dim(index1)[1]){
      weights[index1[i,1],index1[i,2],] <- compute.weights(pd[index1[i,1],index1[i,2]],hd[index1[i,1],index1[i,2]])
    }
    assign(paste0("weights_",crops.gs[cr]),weights)
  }
  
  # read climate data and compute growing season T and P
  
  #dT <- pd/max(pd,na.rm=T)*3
  #dW <- 1.2
  #N <- 200
  #C <- 510
  
  #path.clim <- "/p/projects/macmit/data/CMIP/CMIP5/pcmdi_merged/"
  #fn <- paste0(path.clim,gcms[gcm],"/pr_mon_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",rcps[rcp],"_195001-208412.nc")
  if(do.cmip5){
    #tabl <- array("NA",dim=c(length(gcms),length(rcps)))
    #for(rcp in 1:length(rcps)){
    for(rcp in 1){
      #for(gcm in 1:length(gcms)){
      parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
        #for(gcm in 1:length(gcms)){
        #for(gcm in 10){
        fn <- paste0(path.clim5,gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"/pr_mon_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",rcps[rcp],"_195001-209912.crugrid.nc")
        if(file.exists(fn)){
          #tabl[gcm,rcp] <- "X"
          #cat("found",fn,"\n")
          pr <- readmap.nc(fn,"pr",starttime=(1980-sy5[gcm])*12+1) #start in 1980
          fn <- paste0(path.clim5,gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"/tas_mon_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",rcps[rcp],"_195001-209912.crugrid.nc")
          tas <- readmap.nc(fn,"tas",starttime=(1980-sy5[gcm])*12+1)
          
          # compute growing season means per crop ####
          for(cr in 1:length(crops.gs)){
            #for(cr in 1){
            tas.gs <- pr.gs <- array(NA,dim=c(dim(pr)[1:2],dim(pr)[3]/12))
            weights <- get(paste0("weights_",crops.gs[cr]))
            # with weights summing up to 1, this computes average tas and average pr
            # as long as growing seasons don't change over time (they don't do so here)
            # it's fine to compute fractional changes in pr from averages rather than sums
            # as the denominator is the same in both pr estimates used for computing the fraction
            for(i in 1:dim(pr.gs)[3]){
              tas.gs[,,i] <- apply(tas[,,c(1:12)+(i-1)*12]*weights,c(1,2),sum)
              pr.gs[,,i] <- apply(pr[,,c(1:12)+(i-1)*12]*weights,c(1,2),sum)
            }
            
            tas.base <- apply(tas.gs[,,1:31],c(1,2),mean) #1980-2010 mean
            pr.base <- apply(pr.gs[,,1:31],c(1,2),mean)
            fn <- paste0(path.rclim5,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_climate.Rdata")
            save(tas.gs,pr.gs,tas.base,pr.base,file=fn)
            dT <- dP <- array(NA,dim=c(dim(tas.gs)[1:2],dim(tas.gs)[3]-30-16))
            dT2 <- dP2 <- array(NA,dim=c(dim(tas.gs)[1:2],dim(tas.gs)[3]-31))
            for(i in 32:(dim(tas.gs)[3]-15)){ #start in 2011
              dT[,,(i-31)] <- apply(tas.gs[,,c(1:31)-16+i],c(1,2),mean,na.rm=T)-tas.base
              dP[,,(i-31)] <- apply(pr.gs[,,c(1:31)-16+i],c(1,2),mean,na.rm=T)/pr.base
            }
            for(i in 32:(dim(tas.gs)[3])){ #start in 2011
              dT2[,,(i-31)] <- tas.gs[,,i]-tas.base
              dP2[,,(i-31)] <- pr.gs[,,i]/pr.base
            }
            if(raw_deltas){
              opath <- paste0(path.clim5,"../growingseason_uncapped_deltas/",rcps[rcp],"/")
              dir.create(opath,showWarnings=F,recursive=T)
              fn2 <- paste0(opath,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
              fn3 <- paste0(opath,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_moving_31year_mean_2011-2084.nc4")
              fn2y <- paste0(opath,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_annual_2011-2099_vs_baseline_mean.nc4")
              fn3y <- paste0(opath,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_annual_2011-2099_vs_baseline_mean.nc4")
              writemap.nc(fn2,list(dP[,360:1,]),var="delta precip",var2="change in 31-year precip mean vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="fraction (-)",title="CMIP5 growing season changes in precip",
                          comment2="fractional changes in precip are computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                          start.year=2011,end.year=2084)
              writemap.nc(fn2y,list(dP2[,360:1,]),var="delta precip",var2="change in annual precip vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="fraction (-)",title="CMIP5 growing season changes in precip",
                          comment2="fractional changes in precip are computed based on differences between the 1980-2010 baseline growing season mean and annual growing season means",
                          start.year=2011,end.year=2099)
              writemap.nc(fn3,list(dT[,360:1,]),var="delta tas",var2="change in 31-year tas mean vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="K",title="CMIP5 growing season changes in daily mean temperature (tas)",
                          comment2="absolute changes in tas are computed based on differences between the 1980-2010 baseline growing season mean and a moving 31-year average",
                          start.year=2011,end.year=2084)
              writemap.nc(fn3y,list(dT2[,360:1,]),var="delta tas",var2="change in annual tas vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="K",title="CMIP5 growing season changes in daily mean temperature (tas)",
                          comment2="absolute changes in tas are computed based on differences between the 1980-2010 baseline and annual growing season means",
                          start.year=2011,end.year=2099)
            }
            if(capped_deltas){
              opath <- paste0(path.clim5,"../growingseason_deltas/",rcps[rcp],"/")
              dir.create(opath,showWarnings=F,recursive=T)
              fn2 <- paste0(opath,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_pr_moving_31year_mean_2011-2084.nc4")
              fn3 <- paste0(opath,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_tas_moving_31year_mean_2011-2084.nc4")
              fn2y <- paste0(opath,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_pr_annual_2011-2099_vs_baseline_mean.nc4")
              fn3y <- paste0(opath,rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_tas_annual_2011-2099_vs_baseline_mean.nc4")
              dT[dT<(-1)] <- -1
              dT[dT>6] <- 6
              dT2[dT2<(-1)] <- -1
              dT2[dT2>6] <- 6
              dP[dP<0.5] <- 0.5
              dP[dP>1.3] <- 1.3
              dP2[dP2<0.5] <- 0.5
              dP2[dP2>1.3] <- 1.3
              
              #writemap.nc <- function(filenmae,data,var="yield",var2=var,crop="",cropf="",units="tDM/ha",mv=1e20,start.year=1980,end.year=2099)
              writemap.nc(fn2,list(dP[,360:1,]),var="delta precip",var2="change in 31-year precip mean vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="fraction (-)",title="CMIP5 growing season changes in precip",
                          comment2="fractional changes in precip are computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                          start.year=2011,end.year=2084)
              writemap.nc(fn2y,list(dP2[,360:1,]),var="delta precip",var2="change in annual precip vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="fraction (-)",title="CMIP5 growing season changes in precip",
                          comment2="fractional changes in precip are computed based on differences between the 1980-2010 baseline growing season mean and annual growing season means",
                          start.year=2011,end.year=2099)
              writemap.nc(fn3,list(dT[,360:1,]),var="delta tas",var2="change in 31-year tas mean vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="K",title="CMIP5 growing season changes in daily mean temperature (tas)",
                          comment2="absolute changes in tas are computed based on differences between the 1980-2010 baseline growing season mean and a moving 31-year average",
                          start.year=2011,end.year=2084)
              writemap.nc(fn3y,list(dT2[,360:1,]),var="delta tas",var2="change in annual tas vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="K",title="CMIP5 growing season changes in daily mean temperature (tas)",
                          comment2="absolute changes in tas are computed based on differences between the 1980-2010 baseline and annual growing season means",
                          start.year=2011,end.year=2099)
            }
          }
        }
        else{
          #cat("missing",fn,"\n")
        }
      }
    }
    # colnames(tabl) <- rcps
    # gcms_p <- gcms
    # gcms_p[which(gcms=="GISS-E2-H")[2]] <- "GISS-E2-H_p2"
    # gcms_p[which(gcms=="GISS-E2-H")[3]] <- "GISS-E2-H_p3"
    # gcms_p[which(gcms=="GISS-E2-R")[2]] <- "GISS-E2-R_p2"
    # gcms_p[which(gcms=="GISS-E2-R")[3]] <- "GISS-E2-R_p3"
    # rownames(tabl) <- gcms_p
    
    #gcm <- 1
  }
  if(do.cmip6)
  {
    for(ssp in 1:length(ssps)){
    #for(ssp in 3){
      #parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
      #for(gcm in c(17,18,21,22)){
      for(gcm in c(10,29)){
        opath <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/")
        fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
        if(file.exists(fn2)){
          cat(fn2,"already exists, skipping\n")
          next
        }
        fn <- paste0(path.clim6,gcms6[gcm],"/",ssps[ssp],"/pr_mon_",gcms6[gcm],"_",ssps[ssp],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_crugrid_1950_2099.nc")
        if(file.exists(fn)){
          #cat("found",fn,"\n")
          #tabl[gcm,rcp] <- "X"
          cat("found",fn,"\n")
          pr <- readmap.nc(fn,"pr",starttime=(1980-1950)*12+1) #start in 1980
          fn <- paste0(path.clim6,gcms6[gcm],"/",ssps[ssp],"/tas_mon_",gcms6[gcm],"_",ssps[ssp],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_crugrid_1950_2099.nc")
          tas <- readmap.nc(fn,"tas",starttime=(1980-1950)*12+1)
          
          # compute growing season means per crop ####
          for(cr in 1:length(crops.gs)){
            #for(cr in 1){
            tas.gs <- pr.gs <- array(NA,dim=c(dim(pr)[1:2],dim(pr)[3]/12))
            weights <- get(paste0("weights_",crops.gs[cr]))
            # with weights summing up to 1, this computes average tas and average pr
            # as long as growing seasons don't change over time (they don't do so here)
            # it's fine to compute fractional changes in pr from averages rather than sums
            # as the denominator is the same in both pr estimates used for computing the fraction
            for(i in 1:dim(pr.gs)[3]){
              tas.gs[,,i] <- apply(tas[,,c(1:12)+(i-1)*12]*weights,c(1,2),sum)
              pr.gs[,,i] <- apply(pr[,,c(1:12)+(i-1)*12]*weights,c(1,2),sum)
            }
            
            tas.base <- apply(tas.gs[,,1:31],c(1,2),mean) #1980-2010 mean
            pr.base <- apply(pr.gs[,,1:31],c(1,2),mean)
            fn <- paste0(path.rclim6,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_climate.Rdata")
            save(tas.gs,pr.gs,tas.base,pr.base,file=fn)
            dT <- dP <- array(NA,dim=c(dim(tas.gs)[1:2],dim(tas.gs)[3]-30-16))
            dT2 <- dP2 <- array(NA,dim=c(dim(tas.gs)[1:2],dim(tas.gs)[3]-31))
            for(i in 32:(dim(tas.gs)[3]-15)){ #start in 2011
              dT[,,(i-31)] <- apply(tas.gs[,,c(1:31)-16+i],c(1,2),mean,na.rm=T)-tas.base
              dP[,,(i-31)] <- apply(pr.gs[,,c(1:31)-16+i],c(1,2),mean,na.rm=T)/pr.base
            }
            for(i in 32:(dim(tas.gs)[3])){ #start in 2011
              dT2[,,(i-31)] <- tas.gs[,,i]-tas.base
              dP2[,,(i-31)] <- pr.gs[,,i]/pr.base
            }
            if(raw_deltas){
              opath <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/")
              dir.create(opath,showWarnings=F,recursive=T)
              fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
              fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_moving_31year_mean_2011-2084.nc4")
              fn2y <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_annual_2011-2099_vs_baseline_mean.nc4")
              fn3y <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_annual_2011-2099_vs_baseline_mean.nc4")
              writemap.nc(fn2,list(dP[,360:1,]),var="delta precip",var2="change in 31-year precip mean vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="fraction (-)",title="CMIP6 growing season changes in precip",
                          comment2="fractional changes in precip are computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                          start.year=2011,end.year=2084)
              writemap.nc(fn2y,list(dP2[,360:1,]),var="delta precip",var2="change in annual precip vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="fraction (-)",title="CMIP6 growing season changes in precip",
                          comment2="fractional changes in precip are computed based on differences between the 1980-2010 baseline growing season mean and annual growing season means",
                          start.year=2011,end.year=2099)
              writemap.nc(fn3,list(dT[,360:1,]),var="delta tas",var2="change in 31-year tas mean vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="K",title="CMIP6 growing season changes in daily mean temperature (tas)",
                          comment2="absolute changes in tas are computed based on differences between the 1980-2010 baseline growing season mean and a moving 31-year average",
                          start.year=2011,end.year=2084)
              writemap.nc(fn3y,list(dT2[,360:1,]),var="delta tas",var2="change in annual tas vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="K",title="CMIP6 growing season changes in daily mean temperature (tas)",
                          comment2="absolute changes in tas are computed based on differences between the 1980-2010 baseline and annual growing season means",
                          start.year=2011,end.year=2099)
            }
            if(capped_deltas){
              opath <- paste0(path.clim6,"../growingseason_deltas/",ssps[ssp],"/")
              dir.create(opath,showWarnings=F,recursive=T)
              fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_pr_moving_31year_mean_2011-2084.nc4")
              fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_tas_moving_31year_mean_2011-2084.nc4")
              fn2y <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_pr_annual_2011-2099_vs_baseline_mean.nc4")
              fn3y <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_tas_annual_2011-2099_vs_baseline_mean.nc4")
              dT[dT<(-1)] <- -1
              dT[dT>6] <- 6
              dT2[dT2<(-1)] <- -1
              dT2[dT2>6] <- 6
              dP[dP<0.5] <- 0.5
              dP[dP>1.3] <- 1.3
              dP2[dP2<0.5] <- 0.5
              dP2[dP2>1.3] <- 1.3
              
              #writemap.nc <- function(filenmae,data,var="yield",var2=var,crop="",cropf="",units="tDM/ha",mv=1e20,start.year=1980,end.year=2099)
              writemap.nc(fn2,list(dP[,360:1,]),var="delta precip",var2="change in 31-year precip mean vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="fraction (-)",title="CMIP6 growing season changes in precip",
                          comment2="fractional changes in precip are computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                          start.year=2011,end.year=2084)
              writemap.nc(fn2y,list(dP2[,360:1,]),var="delta precip",var2="change in annual precip vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="fraction (-)",title="CMIP6 growing season changes in precip",
                          comment2="fractional changes in precip are computed based on differences between the 1980-2010 baseline growing season mean and annual growing season means",
                          start.year=2011,end.year=2099)
              writemap.nc(fn3,list(dT[,360:1,]),var="delta tas",var2="change in 31-year tas mean vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="K",title="CMIP6 growing season changes in daily mean temperature (tas)",
                          comment2="absolute changes in tas are computed based on differences between the 1980-2010 baseline growing season mean and a moving 31-year average",
                          start.year=2011,end.year=2084)
              writemap.nc(fn3y,list(dT2[,360:1,]),var="delta tas",var2="change in annual tas vs. 1980-2010 baseline",crop=crops.gs2[cr],cropf=crops.nice[cr],
                          units="K",title="CMIP6 growing season changes in daily mean temperature (tas)",
                          comment2="absolute changes in tas are computed based on differences between the 1980-2010 baseline and annual growing season means",
                          start.year=2011,end.year=2099)
            }
          }
        }
        else{
          cat("missing",fn,"\n")
        }
      }
    }
    # colnames(tabl) <- rcps
    # gcms_p <- gcms
    # gcms_p[which(gcms=="GISS-E2-H")[2]] <- "GISS-E2-H_p2"
    # gcms_p[which(gcms=="GISS-E2-H")[3]] <- "GISS-E2-H_p3"
    # gcms_p[which(gcms=="GISS-E2-R")[2]] <- "GISS-E2-R_p2"
    # gcms_p[which(gcms=="GISS-E2-R")[3]] <- "GISS-E2-R_p3"
    # rownames(tabl) <- gcms_p
    
    #gcm <- 1
    
  }
} #if(T)

# compute emulated iwd fields
if(F){
#   do.standard <- T
#   do.lowinput <- F
#   do.highinput <- F
#   #fn <- paste0(path.clim,gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"/pr_mon_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",rcps[rcp],"_195001-208412.crugrid.nc")
#   
#   
#   # read emulator parameters and compute yield changes #### 
#   
#   for(cr in 1:length(crops.param)){
#     #for(cr in 3){
#     # read fertilizer data 
#     nfert <- readmap.nc(paste0(path.fert,"agmip_",crops.fert[cr],"_apprate_fill_NPK_0.5.nc4"),"Napprate",lo="longitude",la="latitude")
#     for(rcp in 1:length(rcps)){
#       #for(rcp in 1){
#       buf <- read.table(paste0("/p/projects/lpjml/input/scenarios/",rcps[rcp],"_CO2_1765-2200.dat"))
#       co2s <- buf$V2[which(buf$V1 %in% c(1981:2084))]
#       #for(gcm in 1:length(gcms)){
#       #for(gcm in which(gcms=="CESM1-WACCM")){
#       parloop <- foreach(gcm=c(1:length(gcms)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
#         #fn <- paste0(path.rclim,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_growingseason_climate.Rdata")
#         fnc2 <- paste0(path.clim,"../growingseason_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_growingseason_pr_moving_31year_mean_2011-2084.nc4")
#         fnc3 <- paste0(path.clim,"../growingseason_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_growingseason_tas_moving_31year_mean_2011-2084.nc4")
#         if(file.exists(fnc2) & file.exists(fnc3)){
#           dP <- readmap.nc(fnc2)
#           dT <- readmap.nc(fnc3)
#           #for(ggcm in 1:length(ggcms)){
#           #parloop <- foreach(ggcm=c(1:length(ggcms)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
#           for(ggcm in 6){
#             if(do.standard){
#               fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A0.nc4")
#               if(file.exists(fn)){
#                 cat("computing iwd for",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
#                 coeff <- readmap.nc(fn,"K_IWD")
#                 iwd_constt <- iwd_constp <- iwd_constc <- iwd <- array(NA,dim=dim(dT))
#                 #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
#                 iwd_base <- emulator.iwd(coeff,360,0,nfert)
#                 for(i in 1:dim(dT)[3]){
#                   iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],nfert)
#                   iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,nfert)
#                   iwd_constp[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],nfert)
#                   iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],nfert)
#                 }
#                 iwd[iwd<0 & !is.na(iwd)] <- 0
#                 iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
#                 iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
#                 iwd_constp[iwd_constp<0 & !is.na(iwd_constp)] <- 0
#                 iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
#                 fn <- paste0(path.riwd,"A0/",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_movingwindow.Rdata")
#                 save(iwd,iwd_constt,iwd_constp,
#                      iwd_constc,iwd_base,
#                      file=fn)
#                 opath <- paste0(path.riwd,"ncdf/A0/",rcps[rcp],"/",ggcms[ggcm],"/")
#                 dir.create(opath,showWarnings=F,recursive=T)
#                 fn2 <- paste0(opath,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_movingwindow_2011_2084.nc4")
#                 fn3 <- paste0(opath,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_baseline_1980_2010_average.nc4")
#                 writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("irrigation water demand"),crop=crops.gs2[cr],cropf=crops.nice[cr],
#                             units="mm",title="emulated crop irrigation water demand for 1980-2010 baseline",
#                             comment2="emulated crop irrigation water demand, computed based on the 1980-2010 baseline 31-year average",
#                             start.year=1995,end.year=1995)
#                 writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("irrigation water demand"),crop=crops.gs2[cr],cropf=crops.nice[cr],
#                             units="mm",title="emulated crop irrigation water demand for CMIP5 growing season changes",
#                             comment2="emulated crop irrigation water demand, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
#                             start.year=2011,end.year=2084)
#               }
#               fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A1.nc4")
#               if(file.exists(fn)){
#                 cat("computing A1 iwd for",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
#                 coeff <- readmap.nc(fn,"K_IWD")
#                 iwd_constt <- iwd_constp <- iwd_constc <- iwd <- array(NA,dim=dim(dT))
#                 #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
#                 iwd_base <- emulator.iwd(coeff,360,0,nfert)
#                 for(i in 1:dim(dT)[3]){
#                   iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],nfert)
#                   iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,nfert)
#                   iwd_constp[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],nfert)
#                   iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],nfert)
#                 }
#                 iwd[iwd<0 & !is.na(iwd)] <- 0
#                 iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
#                 iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
#                 iwd_constp[iwd_constp<0 & !is.na(iwd_constp)] <- 0
#                 iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
#                 fn <- paste0(path.riwd,"A1/",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_movingwindow.Rdata")
#                 save(iwd,iwd_constt,iwd_constp,
#                      iwd_constc,iwd_base,
#                      file=fn)
#                 opath <- paste0(path.riwd,"ncdf/A1/",rcps[rcp],"/",ggcms[ggcm],"/")
#                 dir.create(opath,showWarnings=F,recursive=T)
#                 fn2 <- paste0(opath,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_movingwindow_2011_2084.nc4")
#                 fn3 <- paste0(opath,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_baseline_1980_2010_average.nc4")
#                 writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("irrigation water demand"),crop=crops.gs2[cr],cropf=crops.nice[cr],
#                             units="mm",title="emulated crop irrigation water demand for CMIP5 growing season changes",
#                             comment2="emulated crop irrigation water demand, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
#                             start.year=2011,end.year=2084)
#                 writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("irrigation water demand"),crop=crops.gs2[cr],cropf=crops.nice[cr],
#                             units="mm",title="emulated crop irrigation water demand for 1980-2010 baseline",
#                             comment2="emulated crop irrigation water demand, computed based on the 1980-2010 baseline 31-year average",
#                             start.year=1995,end.year=1995)
#               } # else do nothing
#             }
#             if(FALSE & do.lowinput){
#               fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A0.nc4")
#               if(file.exists(fn)){
#                 cat("computing yields for",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
#                 coeff_rf <- readmap.nc(fn,"K_rf")
#                 coeff_ir <- readmap.nc(fn,"K_ir")
#                 yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
#                   yields_rf <- yields_ir <- array(NA,dim=dim(dT))
#                 #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
#                 yields_rf_base <- emulator(coeff_rf,360,0,1,10)
#                 yields_ir_base <- emulator(coeff_ir,360,0,1,10,irrig=T)
#                 #for(i in 31:dim(tas.gs)[3]){
#                 #for(i in 31:(dim(tas.gs)[3]-15)){
#                 for(i in 1:dim(dT)[3]){
#                   yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],10)
#                   yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],10,irrig=T)
#                   yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],10)
#                   yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],10,irrig=T)
#                   yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,10)
#                   yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,10,irrig=T)
#                   yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],10)
#                   yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],10,irrig=T)
#                 }
#                 yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
#                 yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
#                 yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
#                 yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
#                 yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
#                 yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
#                 yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
#                 yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
#                 yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
#                 yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
#                 fn <- paste0(path.ryield,"A0_N10/",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_movingwindow.Rdata")
#                 save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
#                      yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
#                      file=fn)
#                 opath <- paste0(path.ryield,"ncdf/A0_N10/",rcps[rcp],"/",ggcms[ggcm],"/")
#                 dir.create(opath,showWarnings=F,recursive=T)
#                 fn2 <- paste0(opath,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_movingwindow_2011_2084.nc4")
#                 fn3 <- paste0(opath,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_baseline_1980_2010_average.nc4")
#                 writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
#                             units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
#                             comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
#                             start.year=2011,end.year=2084)
#                 writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
#                             units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
#                             comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
#                             start.year=1995,end.year=1995)
#                 
#               }
#             }
#             if(FALSE & do.highinput){
#               fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A0.nc4")
#               if(file.exists(fn)){
#                 cat("computing yields for",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
#                 coeff_rf <- readmap.nc(fn,"K_rf")
#                 coeff_ir <- readmap.nc(fn,"K_ir")
#                 yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
#                   yields_rf <- yields_ir <- array(NA,dim=dim(dT))
#                 #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
#                 yields_rf_base <- emulator(coeff_rf,360,0,1,200)
#                 yields_ir_base <- emulator(coeff_ir,360,0,1,200,irrig=T)
#                 #for(i in 31:dim(tas.gs)[3]){
#                 #for(i in 31:(dim(tas.gs)[3]-15)){
#                 for(i in 1:dim(dT)[3]){
#                   yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],200)
#                   yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],200,irrig=T)
#                   yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],200)
#                   yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],200,irrig=T)
#                   yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,200)
#                   yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,200,irrig=T)
#                   yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],200)
#                   yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],200,irrig=T)
#                 }
#                 yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
#                 yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
#                 yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
#                 yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
#                 yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
#                 yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
#                 yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
#                 yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
#                 yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
#                 yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
#                 fn <- paste0(path.ryield,"A0_N200/",rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_movingwindow.Rdata")
#                 save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
#                      yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
#                      file=fn)
#                 opath <- paste0(path.ryield,"ncdf/A0_N200/",rcps[rcp],"/",ggcms[ggcm],"/")
#                 dir.create(opath,showWarnings=F,recursive=T)
#                 fn2 <- paste0(opath,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_movingwindow_2011_2084.nc4")
#                 fn3 <- paste0(opath,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_baseline_1980_2010_average.nc4")
#                 writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
#                             units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
#                             comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
#                             start.year=2011,end.year=2084)
#                 writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
#                             units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
#                             comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
#                             start.year=1995,end.year=1995)
#               }
#             }
#           }
#         } # else do nothing
#       }
#       if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n",job[["messages"]]))
#     }
#   }
#   
#   #yields_rf <- yields_ir <- array(NA,dim=c(dim(tas.base)[1:2],dim(tas.gs)[3]-30))
#   #for(i in 31:dim(tas.gs)[3]){
#   #  yields_rf[,,i-30] <- emulator(coeff_rf,360,tas.gs[,,i]-tas.base,pr.gs[,,i]/pr.base,200)
#   #}
#   #yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
#   #
#   #yields <- emulator(coeff_rf,C,dT,dW,N)
} #if(T)

do.standard <- TRUE
do.lowinput <- FALSE
do.highinput <- FALSE
if(T){
  if(do.cmip5){
    #fn <- paste0(path.clim,gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"/pr_mon_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",rcps[rcp],"_195001-208412.crugrid.nc")
    
    
    # read emulator parameters and compute iwd changes #### 
    
    for(cr in 1:length(crops.param)){
    #for(cr in 3){
      # read fertilizer data 
      nfert <- readmap.nc(paste0(path.fert,"agmip_",crops.fert[cr],"_apprate_fill_NPK_0.5.nc4"),"Napprate",lo="longitude",la="latitude")
      for(rcp in 1:length(rcps)){
        #for(rcp in 1){
        buf <- read.table(paste0("/p/projects/lpjml/input/scenarios/",rcps[rcp],"_CO2_1765-2200.dat"))
        co2v <- buf$V2[which(buf$V1 %in% c(2011:2084))]
        co2s <- array(NA,dim=c(720,360,74))
        for(i in 1:dim(co2s)[3]){
          co2s[,,i] <- co2v[i]
        }
        #for(gcm in gcmdo){
        #for(gcm in 1:length(gcms)){
        #for(gcm in which(gcms=="CESM1-WACCM")){
        parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
          #fn <- paste0(path.rclim,rcps[rcp],"_",gcms5[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_growingseason_climate.Rdata")
          fnc2 <- paste0(path.clim5,"../growingseason_uncapped_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
          fnc3 <- paste0(path.clim5,"../growingseason_uncapped_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_moving_31year_mean_2011-2084.nc4")
          if(file.exists(fnc2) & file.exists(fnc3)){
            dP <- readmap.nc(fnc2)
            dT <- readmap.nc(fnc3)
            dT[dT<(-1)] <- -1
            tem <- abind(co2s,dT)
            co2s <- aperm(apply(tem,c(1,2),replace.with.first.miss,2),c(2,3,1))
            dT[dT>6] <- 6
            dP[dP<0.5] <- 0.5
            dP[dP>1.3] <- 1.3
            for(ggcm in 1:length(ggcms)){
              #parloop <- foreach(ggcm=c(1:length(ggcms)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
              #for(ggcm in 6){
              if(do.standard){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A0.nc4")
                if(file.exists(fn)){
                  cat("computing iwd for",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff <- readmap.nc(fn,"K_IWD")
                  #coeff_ir <- readmap.nc(fn,"K_ir")
                  iwd_constt <- iwd_constc <- 
                    iwd <- array(NA,dim=dim(dT))
                  #iwd_base <- iwd_ir_base <- array(NA,dim=c(dim(tas.base)))
                  iwd_base <- emulator.iwd(coeff,360,0,nfert)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],nfert)
                    iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,nfert)
                    iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],nfert)
                  }
                  iwd[iwd<0 & !is.na(iwd)] <- 0
                  iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
                  iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
                  iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
                  fn <- paste0(path.riwd5,"A0/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_movingwindow_v",version,".Rdata")
                  save(iwd,iwd_constt,iwd_constc,iwd_base,
                       file=fn)
                  opath <- paste0(path.riwd5,"ncdf/A0/",rcps[rcp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for CMIP5 growing season changes",
                              comment2="emulated crop iwd, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for 1980-2010 baseline",
                              comment2="emulated crop iwd, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                }
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A1.nc4")
                if(file.exists(fn)){
                  cat("computing iwd for",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff <- readmap.nc(fn,"K_IWD")
                  iwd_constt <- iwd_constc <- 
                    iwd <- array(NA,dim=dim(dT))
                  #iwd_base <- iwd_ir_base <- array(NA,dim=c(dim(tas.base)))
                  iwd_base <- emulator.iwd(coeff,360,0,nfert)
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],nfert)
                    iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,nfert)
                    iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],nfert)
                  }
                  iwd[iwd<0 & !is.na(iwd)] <- 0
                  iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
                  iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
                  iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
                  fn <- paste0(path.riwd5,"A1/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_movingwindow_v",version,".Rdata")
                  save(iwd,iwd_constt,
                       iwd_constc,iwd_base,
                       file=fn)
                  opath <- paste0(path.riwd5,"ncdf/A1/",rcps[rcp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for CMIP5 growing season changes",
                              comment2="emulated crop iwd, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for 1980-2010 baseline",
                              comment2="emulated crop iwd, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                } # else do nothing
              }
              if(do.lowinput){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A0.nc4")
                if(file.exists(fn)){
                  cat("computing iwd for",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff <- readmap.nc(fn,"K_IWD")
                  iwd_constt <- iwd_constc <- 
                    iwd <- array(NA,dim=dim(dT))
                  #iwd_base <- iwd_ir_base <- array(NA,dim=c(dim(tas.base)))
                  iwd_base <- emulator.iwd(coeff,360.0,10)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],10)
                    iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,10)
                    iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],10)
                  }
                  iwd[iwd<0 & !is.na(iwd)] <- 0
                  iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
                  iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
                  iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
                  fn <- paste0(path.riwd5,"A0_N10/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_iwd_movingwindow_v",version,".Rdata")
                  save(iwd,iwd_constt,iwd_constc,iwd_base,
                       file=fn)
                  opath <- paste0(path.riwd5,"ncdf/A0_N10/",rcps[rcp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_iwd_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_iwd_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd",),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for CMIP5 growing season changes",
                              comment2="emulated crop iwd, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for 1980-2010 baseline",
                              comment2="emulated crop iwd, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                  
                }
              }
              if(do.highinput){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A0.nc4")
                if(file.exists(fn)){
                  cat("computing iwd for",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff <- readmap.nc(fn,"K_IWD")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  iwd_constt <- iwd_constc <- 
                    iwd <- array(NA,dim=dim(dT))
                  #iwd_base <- iwd_ir_base <- array(NA,dim=c(dim(tas.base)))
                  iwd_base <- emulator.iwd(coeff,360.0,200)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],200)
                    iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,200)
                    iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],200)
                  }
                  iwd[iwd<0 & !is.na(iwd)] <- 0
                  iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
                  iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
                  iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
                  fn <- paste0(path.riwd5,"A0_N200/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_iwd_movingwindow_v",version,".Rdata")
                  save(iwd,iwd_constt,iwd_constc,iwd_base,
                       file=fn)
                  opath <- paste0(path.riwd5,"ncdf/A0_N200/",rcps[rcp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_iwd_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_iwd_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for CMIP5 growing season changes",
                              comment2="emulated crop iwd, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for 1980-2010 baseline",
                              comment2="emulated crop iwd, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                }
              }
            }
          } # else do nothing
        }
        if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n",job[["messages"]]))
      }
    }
    

  }
  if(do.cmip6){
    for(cr in 1:length(crops.param)){
    #for(cr in 5){
      # read fertilizer data 
      nfert <- readmap.nc(paste0(path.fert,"agmip_",crops.fert[cr],"_apprate_fill_NPK_0.5.nc4"),"Napprate",lo="longitude",la="latitude")
      # skipping ssp119 for which we're still missing CO2 data
      for(ssp in 1:length(ssps)){
        #for(rcp in 1){
        buf <- read.table(paste0("/p/projects/lpjml/input/scenarios/ISIMIP3b/",rssps[ssp],"_CO2_1765_2100.dat"))
        co2v <- buf$V2[which(buf$V1 %in% c(2011:2084))]
        co2v[co2v>810] <- 810
        co2s <- array(NA,dim=c(720,360,74))
        for(i in 1:dim(co2s)[3]){
          co2s[,,i] <- co2v[i]
        }
        #for(gcm in c(10,29)){
        #parloop <- foreach(gcm=c(17,18,21,22), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
          #for(gcm in c(17,18,21,22)){
          #for(gcm in which(gcms6=="CESM1-WACCM")){
        parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
          fnc2 <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
          fnc3 <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_moving_31year_mean_2011-2084.nc4")
          if(file.exists(fnc2) & file.exists(fnc3)){
            dP <- readmap.nc(fnc2)
            dT <- readmap.nc(fnc3)
            dT[dT<(-1)] <- -1
            tem <- abind(co2s,dT)
            co2s <- aperm(apply(tem,c(1,2),replace.with.first.miss,2),c(2,3,1))
            dT[dT>6] <- 6
            dP[dP<0.5] <- 0.5
            dP[dP>1.3] <- 1.3
            for(ggcm in 1:length(ggcms)){
              #parloop <- foreach(ggcm=c(1:length(ggcms)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
              #for(ggcm in 6){
              if(do.standard){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A0.nc4")
                if(file.exists(fn)){
                  cat("computing iwd for",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff <- readmap.nc(fn,"K_IWD")
                  iwd_constt <- iwd_constc <- 
                    iwd <- array(NA,dim=dim(dT))
                  #iwd_base <- iwd_ir_base <- array(NA,dim=c(dim(tas.base)))
                  iwd_base <- emulator.iwd(coeff,360,0,nfert)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],nfert)
                    iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,nfert)
                    iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],nfert)
                  }
                  iwd[iwd<0 & !is.na(iwd)] <- 0
                  iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
                  iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
                  iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
                  fn1 <- paste0(path.riwd6,"A0/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_movingwindow_v",version,".Rdata")
                  save(iwd,iwd_constt,iwd_constc,iwd_base,
                       file=fn1)
                  opath <- paste0(path.riwd6,"ncdf/A0/",ssps[ssp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_iwd_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for CMIP5 growing season changes",
                              comment2="emulated crop iwd, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for 1980-2010 baseline",
                              comment2="emulated crop iwd, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                } else {
                  cat("cannot find",fn,"\n")
                }
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A1.nc4")
                if(file.exists(fn)){
                  cat("computing iwd for",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff <- readmap.nc(fn,"K_IWD")
                  iwd_constt <- iwd_constc <- 
                    iwd <- array(NA,dim=dim(dT))
                  #iwd_base <- iwd_ir_base <- array(NA,dim=c(dim(tas.base)))
                  iwd_base <- emulator.iwd(coeff,360,0,nfert)
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],nfert)
                    iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,nfert)
                    iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],nfert)
                  }
                  iwd[iwd<0 & !is.na(iwd)] <- 0
                  iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
                  iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
                  iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
                  fn1 <- paste0(path.riwd6,"A1/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_movingwindow_v",version,".Rdata")
                  save(iwd,iwd_constt,iwd_constc,iwd_base,
                       file=fn1)
                  opath <- paste0(path.riwd6,"ncdf/A1/",ssps[ssp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_iwd_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for CMIP5 growing season changes",
                              comment2="emulated crop iwd, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for 1980-2010 baseline",
                              comment2="emulated crop iwd, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                } 
                else {
                  cat("cannot find",fn,"\n")
                }
              }
              if(do.lowinput){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A0.nc4")
                if(file.exists(fn)){
                  cat("computing iwd for",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff <- readmap.nc(fn,"K_IWD")
                  iwd_constt <- iwd_constc <- 
                    iwd <- array(NA,dim=dim(dT))
                  #iwd_base <- iwd_ir_base <- array(NA,dim=c(dim(tas.base)))
                  iwd_base <- emulator.iwd(coeff,360.0,10)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],10)
                    iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,10)
                    iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],10)
                  }
                  iwd[iwd<0 & !is.na(iwd)] <- 0
                  iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
                  iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
                  iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
                  fn <- paste0(path.riwd6,"A0_N10/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_iwd_movingwindow_v",version,".Rdata")
                  save(iwd,iwd_constt,iwd_constc,iwd_base,
                       file=fn)
                  opath <- paste0(path.riwd6,"ncdf/A0_N10/",ssps[ssp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_iwd_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_iwd_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for CMIP5 growing season changes",
                              comment2="emulated crop iwd, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for 1980-2010 baseline",
                              comment2="emulated crop iwd, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                  
                }
              }
              if(do.highinput){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_IWD_A0.nc4")
                if(file.exists(fn)){
                  cat("computing iwd for",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff <- readmap.nc(fn,"K_IWD")
                  iwd_constt <- iwd_constc <- 
                    iwd <- array(NA,dim=dim(dT))
                  #iwd_base <- iwd_ir_base <- array(NA,dim=c(dim(tas.base)))
                  iwd_base <- emulator.iwd(coeff,360.0,200)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    iwd[,,i] <- emulator.iwd(coeff,co2s[,,i],dT[,,i],200)
                    iwd_constt[,,i] <- emulator.iwd(coeff,co2s[,,i],0,200)
                    iwd_constc[,,i] <- emulator.iwd(coeff,360,dT[,,i],200)
                  }
                  iwd[iwd<0 & !is.na(iwd)] <- 0
                  iwd_base[iwd_base<0 & !is.na(iwd_base)] <- 0
                  iwd_constt[iwd_constt<0 & !is.na(iwd_constt)] <- 0
                  iwd_constc[iwd_constc<0 & !is.na(iwd_constc)] <- 0
                  fn <- paste0(path.riwd6,"A0_N200/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_iwd_movingwindow_v",version,".Rdata")
                  save(iwd,iwd_constt,iwd_constc,iwd_base,
                       file=fn)
                  opath <- paste0(path.riwd6,"ncdf/A0_N200/",ssps[ssp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_iwd_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_iwd_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(iwd[,360:1,]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for CMIP5 growing season changes",
                              comment2="emulated crop iwd, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(iwd_base[,360:1]),var=c("iwd"),var2=c("iwd"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="mm",title="emulated crop iwd for 1980-2010 baseline",
                              comment2="emulated crop iwd, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                }
              }
            }
          } # else do nothing
        }
        if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n",job[["messages"]]))
      }
    }
  }
} #if(T)


# find matching SSP/RCP vs. GCMs matrix
if(F){
  match5 <- array("",dim=c(length(rcps),length(gcms5)))
  match6 <- array("",dim=c(length(ssps),length(gcms6)))
  # cmip5
  for(rcp in 1:length(rcps)){
    for(gcm in 1:length(gcms5)){
      fn <- paste0(path.clim5,"../growingseason_uncapped_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[1],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
      if(file.exists(fn)){
        match5[rcp,gcm] <- "X"
      }
    }
  }
  rownames(match5) <- rcps.nice
  colnames(match5) <- gcms5
  write.csv2(match5,paste0(path.figs,"cmip5_","GCM_RCP_matches.csv"))
  # cmip6
  for(ssp in 1:length(ssps)){
    for(gcm in 1:length(gcms5)){
      fn <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[1],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
      if(file.exists(fn)){
        match6[ssp,gcm] <- "X"
      }
    }
  }
  rownames(match6) <- ssps.nice
  colnames(match6) <- gcms6
  write.csv2(match6,paste0(path.figs,"cmip6_","GCM_SSP_matches.csv"))
} 

# compute distributions of changes in growing season conditions
if(F){
  minP <- 0
  maxP <- 3
  minT <- -2
  maxT <- 15
  for(cr in 1:length(crops.param)){
    # cmip5
    if(do.cmip5){
      for(rcp in 1:length(rcps)){
        vT <- vP <- NULL
        cappedP <- cappedT <- 0
        for(gcm in 1:length(gcms5)){
        #parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
          #fnc2 <- paste0(path.clim5,"../growingseason_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_pr_moving_31year_mean_2011-2084.nc4")
          #fnc3 <- paste0(path.clim5,"../growingseason_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_tas_moving_31year_mean_2011-2084.nc4")
          fnc2 <- paste0(path.clim5,"../growingseason_uncapped_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
          fnc3 <- paste0(path.clim5,"../growingseason_uncapped_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_moving_31year_mean_2011-2084.nc4")
          if(file.exists(fnc2) & file.exists(fnc3)){
            dP <- readmap.nc(fnc2)
            # only look at last time slice
            dP <- dP[,,dim(dP)[3]]
            rP <- as.vector(dP[is.finite(dP)])
            cappedP <- cappedP+length(rP[rP>maxP])
            #rP[rP>maxP] <- maxP
            #rP[rP<minP] <- minP
            dT <- readmap.nc(fnc3)
            # only look at last time slice
            dT <- dT[,,dim(dT)[3]]
            rT <- as.vector(dT[is.finite(dT)])
            #rT[rT>maxT] <- maxT
            #rT[rT<minT] <- minT
            vP <- c(vP,rP)
            vT <- c(vT,rT)
          }
        }
        cat(rcps[rcp],crops.param[cr],"capped",cappedP,"precip values at",maxP,"\n")
        fn <- paste0(path.clim5,"../cmip5_growingseaon_delta_hist_",rcps[rcp],"_",crops.param[cr],".Rdata")
        save(vP,vT,cappedP,file=fn)
      }
    }
    # cmip6
    if(do.cmip6){
      for(ssp in 1:length(ssps)){
        vT <- vP <- NULL
        cappedP <- cappedT <- 0
        for(gcm in 1:length(gcms6)){
        #parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
          fnc2 <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
          fnc3 <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_moving_31year_mean_2011-2084.nc4")
          if(file.exists(fnc2) & file.exists(fnc3)){
            dP <- readmap.nc(fnc2)
            # only look at last time slice
            dP <- dP[,,dim(dP)[3]]
            rP <- as.vector(dP[is.finite(dP)])
            cappedP <- cappedP+length(rP[rP>maxP])
            #rP[rP>maxP] <- maxP
            #rP[rP<minP] <- minP
            dT <- readmap.nc(fnc3)
            # only look at last time slice
            dT <- dT[,,dim(dT)[3]]
            rT <- as.vector(dT[is.finite(dT)])
            #rT[rT>maxT] <- maxT
            #rT[rT<minT] <- minT
            vP <- c(vP,rP)
            vT <- c(vT,rT)
          }
        }
        fn <- paste0(path.clim6,"../cmip6_growingseaon_delta_hist_",ssps[ssp],"_",crops.param[cr],".Rdata")
        cat(fn,range(dP,na.rm=T),range(dT,na.rm=T),"\n")
        save(vT,vP,cappedP,file=fn)
      }
    }
  }
}
if(F){
  #for(rcp in 1:length(rcps)){
  for(rcp in c(1,2,4)){
    fn1 <- paste0(path.clim5,"../cmip5_all_growingseaon_delta_",rcps[rcp],".Rdata")
    if(file.exists(fn1)){
      cat("reading",fn1,"\n")
      load(fn1)
    } else {
      allT5 <- allP5 <- NULL
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        fn <- paste0(path.clim5,"../cmip5_growingseaon_delta_hist_",rcps[rcp],"_",crops.param[cr],".Rdata")
        load(fn)
        allT5 <- c(allT5,vT)
        allP5 <- c(allP5,vP)
        assign(paste0(crops.param[cr],"_cmip5_",rcps[rcp],"_vT"),vT)
        assign(paste0(crops.param[cr],"_cmip5_",rcps[rcp],"_vP"),vP)
      }
      assign(paste0("all_cmip5_",rcps[rcp],"_vT"),allT5)
      assign(paste0("all_cmip5_",rcps[rcp],"_vP"),allP5)
      do.call(save,list(paste0("all_cmip5_",rcps[rcp],"_vT"),paste0("all_cmip5_",rcps[rcp],"_vP"),
                        paste0(crops.param[1],"_cmip5_",rcps[rcp],"_vT"),paste0(crops.param[1],"_cmip5_",rcps[rcp],"_vP"),
                        paste0(crops.param[2],"_cmip5_",rcps[rcp],"_vT"),paste0(crops.param[2],"_cmip5_",rcps[rcp],"_vP"),
                        paste0(crops.param[3],"_cmip5_",rcps[rcp],"_vT"),paste0(crops.param[3],"_cmip5_",rcps[rcp],"_vP"),
                        paste0(crops.param[4],"_cmip5_",rcps[rcp],"_vT"),paste0(crops.param[4],"_cmip5_",rcps[rcp],"_vP"),
                        paste0(crops.param[5],"_cmip5_",rcps[rcp],"_vT"),paste0(crops.param[5],"_cmip5_",rcps[rcp],"_vP"),
                        file=fn1))
    }
  }
  for(ssp in 1:length(ssps)){
    #for(ssp in 2){
    fn1 <- paste0(path.clim6,"../cmip6_all_growingseaon_delta_",ssps[ssp],".Rdata")
    if(file.exists(fn1)){
      cat("reading",fn1,"\n")
      load(fn1)
    } else {
      allT6 <- allP6 <- NULL
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        fn <- paste0(path.clim6,"../cmip6_growingseaon_delta_hist_",ssps[ssp],"_",crops.param[cr],".Rdata")
        load(fn)
        allT6 <- c(allT6,vT)
        allP6 <- c(allP6,vP)
        assign(paste0(crops.param[cr],"_cmip6_",ssps[ssp],"_vT"),vT)
        assign(paste0(crops.param[cr],"_cmip6_",ssps[ssp],"_vP"),vP)
      }
      assign(paste0("all_cmip6_",ssps[ssp],"_vT"),allT6)
      assign(paste0("all_cmip6_",ssps[ssp],"_vP"),allP6)
      do.call(save,list(paste0("all_cmip6_",ssps[ssp],"_vT"),paste0("all_cmip6_",ssps[ssp],"_vP"),
                        paste0(crops.param[1],"_cmip6_",ssps[ssp],"_vT"),paste0(crops.param[1],"_cmip6_",ssps[ssp],"_vP"),
                        paste0(crops.param[2],"_cmip6_",ssps[ssp],"_vT"),paste0(crops.param[2],"_cmip6_",ssps[ssp],"_vP"),
                        paste0(crops.param[3],"_cmip6_",ssps[ssp],"_vT"),paste0(crops.param[3],"_cmip6_",ssps[ssp],"_vP"),
                        paste0(crops.param[4],"_cmip6_",ssps[ssp],"_vT"),paste0(crops.param[4],"_cmip6_",ssps[ssp],"_vP"),
                        paste0(crops.param[5],"_cmip6_",ssps[ssp],"_vT"),paste0(crops.param[5],"_cmip6_",ssps[ssp],"_vP"),
                        file=fn1))
    }
  }  
  # this seems to exhaust memory too fast, split into per crop commands
  # vioplot(list(all=all_cmip6_ssp245_vT,
  #              maize=maize_cmip6_ssp245_vT,rice=rice_cmip6_ssp245_vT,soy=soy_cmip6_ssp245_vT,swheat=spring_wheat_cmip6_ssp245_vT,wwheat=winter_wheat_cmip6_ssp245_vT),
  #         col=col.rcp2[1])
  # boxplot_gs(paste0(path.figs,"growing_seaon_Tchanges_cmip5_vs_cmip6_rcp45_box.png"),
  #            all_cmip6_ssp245_vT,maize_cmip6_ssp245_vT,rice_cmip6_ssp245_vT,soy_cmip6_ssp245_vT,spring_wheat_cmip6_ssp245_vT,winter_wheat_cmip6_ssp245_vT,
  #            all_cmip5_rcp45_vT,maize_cmip5_rcp45_vT,rice_cmip5_rcp45_vT,soy_cmip5_rcp45_vT,spring_wheat_cmip5_rcp45_vT,winter_wheat_cmip5_rcp45_vT,
  #            col.rcp2[2],col.rcp2[3])
}
if(F){
  fn1 <- paste0(path.clim5,"../cmip5_all_growingseaon_delta_rcp26.Rdata")
  load(fn1)
  fn1 <- paste0(path.clim6,"../cmip6_all_growingseaon_delta_ssp126.Rdata")
  load(fn1)
  boxplot_gs2(paste0(path.figs,"growing_seaon_Tchanges_cmip5_vs_cmip6_rcp26_box_full"),
              maize_cmip6_ssp126_vT,rice_cmip6_ssp126_vT,soy_cmip6_ssp126_vT,spring_wheat_cmip6_ssp126_vT,winter_wheat_cmip6_ssp126_vT,
              maize_cmip5_rcp26_vT,rice_cmip5_rcp26_vT,soy_cmip5_rcp26_vT,spring_wheat_cmip5_rcp26_vT,winter_wheat_cmip5_rcp26_vT,
              col.rcp2[2],col.rcp2[3],reduce=F,pdf=T,xl="growing season mean T change (K)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Pchanges_cmip5_vs_cmip6_rcp26_box_full"),
              (maize_cmip6_ssp126_vP-1)*100,(rice_cmip6_ssp126_vP-1)*100,(soy_cmip6_ssp126_vP-1)*100,(spring_wheat_cmip6_ssp126_vP-1)*100,(winter_wheat_cmip6_ssp126_vP-1)*100,
              (maize_cmip5_rcp26_vP-1)*100,(rice_cmip5_rcp26_vP-1)*100,(soy_cmip5_rcp26_vP-1)*100,(spring_wheat_cmip5_rcp26_vP-1)*100,(winter_wheat_cmip5_rcp26_vP-1)*100,
              col.rcp2[2],col.rcp2[3],reduce=F,pdf=T,xl="growing season mean P change (%)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Tdiff_cmip5_vs_cmip6_rcp26_box_full"),
              maize_cmip6_ssp126_vT,rice_cmip6_ssp126_vT,soy_cmip6_ssp126_vT,spring_wheat_cmip6_ssp126_vT,winter_wheat_cmip6_ssp126_vT,
              maize_cmip5_rcp26_vT,rice_cmip5_rcp26_vT,soy_cmip5_rcp26_vT,spring_wheat_cmip5_rcp26_vT,winter_wheat_cmip5_rcp26_vT,
              col.rcp2[2],col.rcp2[3],reduce=F,diff=T,pdf=T,xl="CMIP5-CMIP6 difference in growing season T change (K)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Pdiff_cmip5_vs_cmip6_rcp26_box_full"),
              (maize_cmip6_ssp126_vP-1)*100,(rice_cmip6_ssp126_vP-1)*100,(soy_cmip6_ssp126_vP-1)*100,(spring_wheat_cmip6_ssp126_vP-1)*100,(winter_wheat_cmip6_ssp126_vP-1)*100,
              (maize_cmip5_rcp26_vP-1)*100,(rice_cmip5_rcp26_vP-1)*100,(soy_cmip5_rcp26_vP-1)*100,(spring_wheat_cmip5_rcp26_vP-1)*100,(winter_wheat_cmip5_rcp26_vP-1)*100,
              col.rcp2[2],col.rcp2[3],reduce=F,diff=T,pdf=T,xl="CMIP5-CMIP6 difference in growing season P change (%)")
  rm(maize_cmip6_ssp126_vP,rice_cmip6_ssp126_vP,soy_cmip6_ssp126_vP,spring_wheat_cmip6_ssp126_vP,winter_wheat_cmip6_ssp126_vP,
     maize_cmip5_rcp26_vP,rice_cmip5_rcp26_vP,soy_cmip5_rcp26_vP,spring_wheat_cmip5_rcp26_vP,winter_wheat_cmip5_rcp26_vP)

  fn1 <- paste0(path.clim5,"../cmip5_all_growingseaon_delta_rcp45.Rdata")
  load(fn1)
  fn1 <- paste0(path.clim6,"../cmip6_all_growingseaon_delta_ssp245.Rdata")
  load(fn1)
  boxplot_gs2(paste0(path.figs,"growing_seaon_Tchanges_cmip5_vs_cmip6_rcp45_box_full"),
              maize_cmip6_ssp245_vT,rice_cmip6_ssp245_vT,soy_cmip6_ssp245_vT,spring_wheat_cmip6_ssp245_vT,winter_wheat_cmip6_ssp245_vT,
              maize_cmip5_rcp45_vT,rice_cmip5_rcp45_vT,soy_cmip5_rcp45_vT,spring_wheat_cmip5_rcp45_vT,winter_wheat_cmip5_rcp45_vT,
              col.rcp2[2],col.rcp2[3],reduce=F,pdf=T,xl="growing season mean T change (K)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Pchanges_cmip5_vs_cmip6_rcp45_box_full"),
              (maize_cmip6_ssp245_vP-1)*100,(rice_cmip6_ssp245_vP-1)*100,(soy_cmip6_ssp245_vP-1)*100,(spring_wheat_cmip6_ssp245_vP-1)*100,(winter_wheat_cmip6_ssp245_vP-1)*100,
              (maize_cmip5_rcp45_vP-1)*100,(rice_cmip5_rcp45_vP-1)*100,(soy_cmip5_rcp45_vP-1)*100,(spring_wheat_cmip5_rcp45_vP-1)*100,(winter_wheat_cmip5_rcp45_vP-1)*100,
              col.rcp2[2],col.rcp2[3],reduce=F,pdf=T,xl="growing season mean P change (%)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Tdiff_cmip5_vs_cmip6_rcp45_box_full"),
              maize_cmip6_ssp245_vT,rice_cmip6_ssp245_vT,soy_cmip6_ssp245_vT,spring_wheat_cmip6_ssp245_vT,winter_wheat_cmip6_ssp245_vT,
              maize_cmip5_rcp45_vT,rice_cmip5_rcp45_vT,soy_cmip5_rcp45_vT,spring_wheat_cmip5_rcp45_vT,winter_wheat_cmip5_rcp45_vT,
              col.rcp2[2],col.rcp2[3],reduce=F,diff=T,pdf=T,xl="CMIP5-CMIP6 difference in growing season T change (K)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Pdiff_cmip5_vs_cmip6_rcp45_box_full"),
              (maize_cmip6_ssp245_vP-1)*100,(rice_cmip6_ssp245_vP-1)*100,(soy_cmip6_ssp245_vP-1)*100,(spring_wheat_cmip6_ssp245_vP-1)*100,(winter_wheat_cmip6_ssp245_vP-1)*100,
              (maize_cmip5_rcp45_vP-1)*100,(rice_cmip5_rcp45_vP-1)*100,(soy_cmip5_rcp45_vP-1)*100,(spring_wheat_cmip5_rcp45_vP-1)*100,(winter_wheat_cmip5_rcp45_vP-1)*100,
              col.rcp2[2],col.rcp2[3],reduce=F,diff=T,pdf=T,xl="CMIP5-CMIP6 difference in growing season P change (%)")
  rm(maize_cmip6_ssp245_vP,rice_cmip6_ssp245_vP,soy_cmip6_ssp245_vP,spring_wheat_cmip6_ssp245_vP,winter_wheat_cmip6_ssp245_vP,
     maize_cmip5_rcp45_vP,rice_cmip5_rcp45_vP,soy_cmip5_rcp45_vP,spring_wheat_cmip5_rcp45_vP,winter_wheat_cmip5_rcp45_vP)

  fn1 <- paste0(path.clim5,"../cmip5_all_growingseaon_delta_rcp85.Rdata")
  load(fn1)
  fn1 <- paste0(path.clim6,"../cmip6_all_growingseaon_delta_ssp585.Rdata")
  load(fn1)
  boxplot_gs2(paste0(path.figs,"growing_seaon_Tchanges_cmip5_vs_cmip6_rcp85_box_full"),
              maize_cmip6_ssp585_vT,rice_cmip6_ssp585_vT,soy_cmip6_ssp585_vT,spring_wheat_cmip6_ssp585_vT,winter_wheat_cmip6_ssp585_vT,
              maize_cmip5_rcp85_vT,rice_cmip5_rcp85_vT,soy_cmip5_rcp85_vT,spring_wheat_cmip5_rcp85_vT,winter_wheat_cmip5_rcp85_vT,
              col.rcp2[2],col.rcp2[3],reduce=F,pdf=T,xl="growing season mean T change (K)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Pchanges_cmip5_vs_cmip6_rcp85_box_full"),
              (maize_cmip6_ssp585_vP-1)*100,(rice_cmip6_ssp585_vP-1)*100,(soy_cmip6_ssp585_vP-1)*100,(spring_wheat_cmip6_ssp585_vP-1)*100,(winter_wheat_cmip6_ssp585_vP-1)*100,
              (maize_cmip5_rcp85_vP-1)*100,(rice_cmip5_rcp85_vP-1)*100,(soy_cmip5_rcp85_vP-1)*100,(spring_wheat_cmip5_rcp85_vP-1)*100,(winter_wheat_cmip5_rcp85_vP-1)*100,
              col.rcp2[2],col.rcp2[3],reduce=F,pdf=T,xl="growing season mean P change (%)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Tdiff_cmip5_vs_cmip6_rcp85_box_full"),
              maize_cmip6_ssp585_vT,rice_cmip6_ssp585_vT,soy_cmip6_ssp585_vT,spring_wheat_cmip6_ssp585_vT,winter_wheat_cmip6_ssp585_vT,
              maize_cmip5_rcp85_vT,rice_cmip5_rcp85_vT,soy_cmip5_rcp85_vT,spring_wheat_cmip5_rcp85_vT,winter_wheat_cmip5_rcp85_vT,
              col.rcp2[2],col.rcp2[3],reduce=F,diff=T,pdf=T,xl="CMIP5-CMIP6 difference in growing season T change (K)")
  boxplot_gs2(paste0(path.figs,"growing_seaon_Pdiff_cmip5_vs_cmip6_rcp85_box_full"),
              (maize_cmip6_ssp585_vP-1)*100,(rice_cmip6_ssp585_vP-1)*100,(soy_cmip6_ssp585_vP-1)*100,(spring_wheat_cmip6_ssp585_vP-1)*100,(winter_wheat_cmip6_ssp585_vP-1)*100,
              (maize_cmip5_rcp85_vP-1)*100,(rice_cmip5_rcp85_vP-1)*100,(soy_cmip5_rcp85_vP-1)*100,(spring_wheat_cmip5_rcp85_vP-1)*100,(winter_wheat_cmip5_rcp85_vP-1)*100,
              col.rcp2[2],col.rcp2[3],reduce=F,diff=T,pdf=T,xl="CMIP5-CMIP6 difference in growing season P change (%)")
  rm(maize_cmip6_ssp585_vP,rice_cmip6_ssp585_vP,soy_cmip6_ssp585_vP,spring_wheat_cmip6_ssp585_vP,winter_wheat_cmip6_ssp585_vP,
     maize_cmip5_rcp85_vP,rice_cmip5_rcp85_vP,soy_cmip5_rcp85_vP,spring_wheat_cmip5_rcp85_vP,winter_wheat_cmip5_rcp85_vP)
  
  # boxplot_gs2(paste0(path.figs,"growing_seaon_Tchanges_cmip5_vs_cmip6_rcp45_box_full.png"),
  #            maize_cmip6_ssp245_vT,rice_cmip6_ssp245_vT,soy_cmip6_ssp245_vT,spring_wheat_cmip6_ssp245_vT,winter_wheat_cmip6_ssp245_vT,
  #            maize_cmip5_rcp45_vT,rice_cmip5_rcp45_vT,soy_cmip5_rcp45_vT,spring_wheat_cmip5_rcp45_vT,winter_wheat_cmip5_rcp45_vT,
  #            col.rcp2[2],col.rcp2[3],reduce=F)
  # boxplot_gs(paste0(path.figs,"growing_seaon_Pchanges_cmip5_vs_cmip6_rcp45_box_full.png"),
  #            all_cmip6_ssp245_vP,maize_cmip6_ssp245_vP,rice_cmip6_ssp245_vP,soy_cmip6_ssp245_vP,spring_wheat_cmip6_ssp245_vP,winter_wheat_cmip6_ssp245_vP,
  #            all_cmip5_rcp45_vP,maize_cmip5_rcp45_vP,rice_cmip5_rcp45_vP,soy_cmip5_rcp45_vP,spring_wheat_cmip5_rcp45_vP,winter_wheat_cmip5_rcp45_vP,
  #            col.rcp2[2],col.rcp2[3],reduce=F)
  # 
  # boxplot_gs(paste0(path.figs,"growing_seaon_Tchanges_cmip5_vs_cmip6_rcp85_box_full.png"),
  #            all_cmip6_ssp585_vT,maize_cmip6_ssp585_vT,rice_cmip6_ssp585_vT,soy_cmip6_ssp585_vT,spring_wheat_cmip6_ssp585_vT,winter_wheat_cmip6_ssp585_vT,
  #            all_cmip5_rcp85_vT,maize_cmip5_rcp85_vT,rice_cmip5_rcp85_vT,soy_cmip5_rcp85_vT,spring_wheat_cmip5_rcp85_vT,winter_wheat_cmip5_rcp85_vT,
  #            col.rcp2[2],col.rcp2[3],reduce=F)
  # boxplot_gs(paste0(path.figs,"growing_seaon_Pchanges_cmip5_vs_cmip6_rcp85_box_full.png"),
  #            all_cmip6_ssp585_vP,maize_cmip6_ssp585_vP,rice_cmip6_ssp585_vP,soy_cmip6_ssp585_vP,spring_wheat_cmip6_ssp585_vP,winter_wheat_cmip6_ssp585_vP,
  #            all_cmip5_rcp85_vP,maize_cmip5_rcp85_vP,rice_cmip5_rcp85_vP,soy_cmip5_rcp85_vP,spring_wheat_cmip5_rcp85_vP,winter_wheat_cmip5_rcp85_vP,
  #            col.rcp2[2],col.rcp2[3],reduce=F)
  # 
  # boxplot_gs(paste0(path.figs,"growing_seaon_Tdiff_cmip5_vs_cmip6_rcp45_box_full.png"),
  #            all_cmip6_ssp245_vT,maize_cmip6_ssp245_vT,rice_cmip6_ssp245_vT,soy_cmip6_ssp245_vT,spring_wheat_cmip6_ssp245_vT,winter_wheat_cmip6_ssp245_vT,
  #            all_cmip5_rcp45_vT,maize_cmip5_rcp45_vT,rice_cmip5_rcp45_vT,soy_cmip5_rcp45_vT,spring_wheat_cmip5_rcp45_vT,winter_wheat_cmip5_rcp45_vT,
  #            col.rcp2[2],col.rcp2[3],reduce=F,diff=T)
  # boxplot_gs(paste0(path.figs,"growing_seaon_Tdiff_cmip5_vs_cmip6_rcp85_box_full.png"),
  #            all_cmip6_ssp585_vT,maize_cmip6_ssp585_vT,rice_cmip6_ssp585_vT,soy_cmip6_ssp585_vT,spring_wheat_cmip6_ssp585_vT,winter_wheat_cmip6_ssp585_vT,
  #            all_cmip5_rcp85_vT,maize_cmip5_rcp85_vT,rice_cmip5_rcp85_vT,soy_cmip5_rcp85_vT,spring_wheat_cmip5_rcp85_vT,winter_wheat_cmip5_rcp85_vT,
  #            col.rcp2[2],col.rcp2[3],reduce=F,diff=T)
  # boxplot_gs(paste0(path.figs,"growing_seaon_Pdiff_cmip5_vs_cmip6_rcp45_box_full.png"),
  #            all_cmip6_ssp245_vP,maize_cmip6_ssp245_vP,rice_cmip6_ssp245_vP,soy_cmip6_ssp245_vP,spring_wheat_cmip6_ssp245_vP,winter_wheat_cmip6_ssp245_vP,
  #            all_cmip5_rcp45_vP,maize_cmip5_rcp45_vP,rice_cmip5_rcp45_vP,soy_cmip5_rcp45_vP,spring_wheat_cmip5_rcp45_vP,winter_wheat_cmip5_rcp45_vP,
  #            col.rcp2[2],col.rcp2[3],reduce=F,diff=T)
  # boxplot_gs(paste0(path.figs,"growing_seaon_Pdiff_cmip5_vs_cmip6_rcp85_box_full.png"),
  #            all_cmip6_ssp585_vP,maize_cmip6_ssp585_vP,rice_cmip6_ssp585_vP,soy_cmip6_ssp585_vP,spring_wheat_cmip6_ssp585_vP,winter_wheat_cmip6_ssp585_vP,
  #            all_cmip5_rcp85_vP,maize_cmip5_rcp85_vP,rice_cmip5_rcp85_vP,soy_cmip5_rcp85_vP,spring_wheat_cmip5_rcp85_vP,winter_wheat_cmip5_rcp85_vP,
  #            col.rcp2[2],col.rcp2[3],reduce=F,diff=T)
  # vioplot_gs(paste0(path.figs,"growing_seaon_changes_cmip5_vs_cmip6_rcp45.png"),
  #            all_cmip6_ssp245_vT,maize_cmip6_ssp245_vT,rice_cmip6_ssp245_vT,soy_cmip6_ssp245_vT,spring_wheat_cmip6_ssp245_vT,winter_wheat_cmip6_ssp245_vT,
  #            all_cmip5_rcp45_vT,maize_cmip5_rcp45_vT,rice_cmip5_rcp45_vT,soy_cmip5_rcp45_vT,spring_wheat_cmip5_rcp45_vT,winter_wheat_cmip5_rcp45_vT,
  #            col.rcp2[2],col.rcp[2],col.rcp2[3],col.rcp[3])
  # vioplot_gs2(paste0(path.figs,"growing_seaon_changes_cmip5_vs_cmip6_rcp45_test.png"),
  #             all_cmip6_ssp245_vT,maize_cmip6_ssp245_vT,rice_cmip6_ssp245_vT,soy_cmip6_ssp245_vT,spring_wheat_cmip6_ssp245_vT,winter_wheat_cmip6_ssp245_vT,
  #             all_cmip5_rcp45_vT,maize_cmip5_rcp45_vT,rice_cmip5_rcp45_vT,soy_cmip5_rcp45_vT,spring_wheat_cmip5_rcp45_vT,winter_wheat_cmip5_rcp45_vT,
  #             col.rcp2[2],col.rcp[2],col.rcp2[3],col.rcp[3])
  # 
  # vioplot_gs(paste0(path.figs,"growing_seaon_changes_cmip5_vs_cmip6_rcp26.png"),
  #            all_cmip6_ssp126_vT,maize_cmip6_ssp126_vT,rice_cmip6_ssp126_vT,soy_cmip6_ssp126_vT,spring_wheat_cmip6_ssp126_vT,winter_wheat_cmip6_ssp126_vT,
  #            all_cmip5_rcp26_vT,maize_cmip5_rcp26_vT,rice_cmip5_rcp26_vT,soy_cmip5_rcp26_vT,spring_wheat_cmip5_rcp26_vT,winter_wheat_cmip5_rcp26_vT,
  #            col.rcp2[2],col.rcp[2],col.rcp2[3],col.rcp[3])
  # 
  vioplot_gs2(fn=paste0(path.figs,"growing_seaon_changes_cmip5_vs_cmip6_rcp85.png"),
             all_cmip6_ssp585_vT,maize_cmip6_ssp585_vT,rice_cmip6_ssp585_vT,soy_cmip6_ssp585_vT,spring_wheat_cmip6_ssp585_vT,winter_wheat_cmip6_ssp585_vT,
             all_cmip5_rcp85_vT,maize_cmip5_rcp85_vT,rice_cmip5_rcp85_vT,soy_cmip5_rcp85_vT,spring_wheat_cmip5_rcp85_vT,winter_wheat_cmip5_rcp85_vT,
             col.rcp2[2],col.rcp[2],col.rcp2[3],col.rcp[3])
  
}

# compute emulated yield fields
do.standard <- TRUE
do.lowinput <- FALSE
do.highinput <- FALSE
#gcmdo <- 45
#if(gcmdo>31) do.cmip6 <- F else do.cmip6 <- T
#crdo <- 3
if(F){
  if(do.cmip5){
    #fn <- paste0(path.clim,gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"/pr_mon_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",rcps[rcp],"_195001-208412.crugrid.nc")
    
    
    # read emulator parameters and compute yield changes #### 
    
    #for(cr in 1:length(crops.param)){
    for(cr in 3){
      # read fertilizer data 
      nfert <- readmap.nc(paste0(path.fert,"agmip_",crops.fert[cr],"_apprate_fill_NPK_0.5.nc4"),"Napprate",lo="longitude",la="latitude")
      for(rcp in 1:length(rcps)){
        #for(rcp in 1){
        buf <- read.table(paste0("/p/projects/lpjml/input/scenarios/",rcps[rcp],"_CO2_1765-2200.dat"))
        co2v <- buf$V2[which(buf$V1 %in% c(2011:2084))]
        co2s <- array(NA,dim=c(720,360,74))
        for(i in 1:dim(co2s)[3]){
          co2s[,,i] <- co2v[i]
        }
        #for(gcm in gcmdo){
        #for(gcm in 1:length(gcms)){
        #for(gcm in which(gcms=="CESM1-WACCM")){
        parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
          #fn <- paste0(path.rclim,rcps[rcp],"_",gcms5[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_growingseason_climate.Rdata")
          fnc2 <- paste0(path.clim5,"../growingseason_uncapped_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
          fnc3 <- paste0(path.clim5,"../growingseason_uncapped_deltas/",rcps[rcp],"/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_moving_31year_mean_2011-2084.nc4")
          if(file.exists(fnc2) & file.exists(fnc3)){
            dP <- readmap.nc(fnc2)
            dT <- readmap.nc(fnc3)
            dP <- readmap.nc(fnc2)
            dT <- readmap.nc(fnc3)
            dT[dT<(-1)] <- -1
            tem <- abind(co2s,dT)
            co2s <- aperm(apply(tem,c(1,2),replace.with.first.miss,2),c(2,3,1))
            dT[dT>6] <- 6
            dP[dP<0.5] <- 0.5
            dP[dP>1.3] <- 1.3
            for(ggcm in 1:length(ggcms)){
            #parloop <- foreach(ggcm=c(1:length(ggcms)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(ggcm in 6){
              if(do.standard){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A0.nc4")
                if(file.exists(fn)){
                  cat("computing yields for",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff_rf <- readmap.nc(fn,"K_rf")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
                    yields_rf <- yields_ir <- array(NA,dim=dim(dT))
                  #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
                  yields_rf_base <- emulator(coeff_rf,360,0,1,nfert)
                  yields_ir_base <- emulator(coeff_ir,360,0,1,nfert,irrig=T)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],nfert)
                    yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],nfert,irrig=T)
                    yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],nfert)
                    yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],nfert,irrig=T)
                    yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,nfert)
                    yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,nfert,irrig=T)
                    yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],nfert)
                    yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],nfert,irrig=T)
                  }
                  yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
                  yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
                  yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
                  yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
                  yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
                  yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
                  yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
                  yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
                  yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
                  yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
                  fn <- paste0(path.ryield5,"A0/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
                  save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
                       yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
                       file=fn)
                  opath <- paste0(path.ryield5,"ncdf/A0/",rcps[rcp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
                              comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
                              comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                }
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A1.nc4")
                if(file.exists(fn)){
                  cat("computing yields for",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff_rf <- readmap.nc(fn,"K_rf")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
                    yields_rf <- yields_ir <- array(NA,dim=dim(dT))
                  #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
                  yields_rf_base <- emulator(coeff_rf,360,0,1,nfert)
                  yields_ir_base <- emulator(coeff_ir,360,0,1,nfert,irrig=T)
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],nfert)
                    yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],nfert,irrig=T)
                    yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],nfert)
                    yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],nfert,irrig=T)
                    yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,nfert)
                    yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,nfert,irrig=T)
                    yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],nfert)
                    yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],nfert,irrig=T)
                  }
                  yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
                  yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
                  yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
                  yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
                  yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
                  yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
                  yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
                  yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
                  yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
                  yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
                  fn <- paste0(path.ryield5,"A1/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow_v",version,".Rdata")
                  save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
                       yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
                       file=fn)
                  opath <- paste0(path.ryield5,"ncdf/A1/",rcps[rcp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
                              comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
                              comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                } # else do nothing
              }
              if(do.lowinput){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A0.nc4")
                if(file.exists(fn)){
                  cat("computing yields for",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff_rf <- readmap.nc(fn,"K_rf")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
                    yields_rf <- yields_ir <- array(NA,dim=dim(dT))
                  #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
                  yields_rf_base <- emulator(coeff_rf,360,0,1,10)
                  yields_ir_base <- emulator(coeff_ir,360,0,1,10,irrig=T)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],10)
                    yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],10,irrig=T)
                    yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],10)
                    yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],10,irrig=T)
                    yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,10)
                    yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,10,irrig=T)
                    yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],10)
                    yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],10,irrig=T)
                  }
                  yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
                  yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
                  yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
                  yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
                  yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
                  yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
                  yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
                  yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
                  yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
                  yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
                  fn <- paste0(path.ryield5,"A0_N10/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_movingwindow_v",version,".Rdata")
                  save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
                       yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
                       file=fn)
                  opath <- paste0(path.ryield5,"ncdf/A0_N10/",rcps[rcp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
                              comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
                              comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                  
                }
              }
              if(do.highinput){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A0.nc4")
                if(file.exists(fn)){
                  cat("computing yields for",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff_rf <- readmap.nc(fn,"K_rf")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
                    yields_rf <- yields_ir <- array(NA,dim=dim(dT))
                  #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
                  yields_rf_base <- emulator(coeff_rf,360,0,1,200)
                  yields_ir_base <- emulator(coeff_ir,360,0,1,200,irrig=T)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],200)
                    yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],200,irrig=T)
                    yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],200)
                    yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],200,irrig=T)
                    yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,200)
                    yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,200,irrig=T)
                    yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],200)
                    yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],200,irrig=T)
                  }
                  yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
                  yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
                  yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
                  yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
                  yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
                  yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
                  yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
                  yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
                  yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
                  yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
                  fn <- paste0(path.ryield5,"A0_N200/",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_movingwindow_v",version,".Rdata")
                  save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
                       yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
                       file=fn)
                  opath <- paste0(path.ryield5,"ncdf/A0_N200/",rcps[rcp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
                              comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
                              comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                }
              }
            }
          } # else do nothing
        }
        if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n",job[["messages"]]))
      }
    }
    
    #yields_rf <- yields_ir <- array(NA,dim=c(dim(tas.base)[1:2],dim(tas.gs)[3]-30))
    #for(i in 31:dim(tas.gs)[3]){
    #  yields_rf[,,i-30] <- emulator(coeff_rf,360,tas.gs[,,i]-tas.base,pr.gs[,,i]/pr.base,200)
    #}
    #yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
    #
    #yields <- emulator(coeff_rf,C,dT,dW,N)
    
  }
  if(do.cmip6){
    for(cr in 1:length(crops.param)){
    #for(cr in 5){
      # read fertilizer data 
      nfert <- readmap.nc(paste0(path.fert,"agmip_",crops.fert[cr],"_apprate_fill_NPK_0.5.nc4"),"Napprate",lo="longitude",la="latitude")
      # skipping ssp119 for which we're still missing CO2 data
      for(ssp in 1:length(ssps)){
        #for(rcp in 1){
        buf <- read.table(paste0("/p/projects/lpjml/input/scenarios/ISIMIP3b/",rssps[ssp],"_CO2_1765_2100.dat"))
        co2v <- buf$V2[which(buf$V1 %in% c(2011:2084))]
        co2v[co2v>810] <- 810
        co2s <- array(NA,dim=c(720,360,74))
        for(i in 1:dim(co2s)[3]){
          co2s[,,i] <- co2v[i]
        }
        #for(gcm in gcmdo){
        #parloop <- foreach(gcm=c(17,18,21,22), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
        #for(gcm in c(17,18,21,22)){
        for(gcm in c(10,29)){
        #for(gcm in which(gcms6=="CESM1-WACCM")){
        #parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
          fnc2 <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_pr_moving_31year_mean_2011-2084.nc4")
          fnc3 <- paste0(path.clim6,"../growingseason_uncapped_deltas/",ssps[ssp],"/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"f",file6[gcm],"_",crops.param[cr],"_growingseason_uncapped_tas_moving_31year_mean_2011-2084.nc4")
          if(file.exists(fnc2) & file.exists(fnc3)){
            dP <- readmap.nc(fnc2)
            dT <- readmap.nc(fnc3)
            dT[dT<(-1)] <- -1
            tem <- abind(co2s,dT)
            co2s <- aperm(apply(tem,c(1,2),replace.with.first.miss,2),c(2,3,1))
            dT[dT>6] <- 6
            dP[dP<0.5] <- 0.5
            dP[dP>1.3] <- 1.3
            for(ggcm in 1:length(ggcms)){
              #parloop <- foreach(ggcm=c(1:length(ggcms)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
              #for(ggcm in 6){
              if(do.standard){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A0.nc4")
                if(file.exists(fn)){
                  cat("computing yields for",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff_rf <- readmap.nc(fn,"K_rf")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
                    yields_rf <- yields_ir <- array(NA,dim=dim(dT))
                  #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
                  yields_rf_base <- emulator(coeff_rf,360,0,1,nfert)
                  yields_ir_base <- emulator(coeff_ir,360,0,1,nfert,irrig=T)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],nfert)
                    yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],nfert,irrig=T)
                    yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],nfert)
                    yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],nfert,irrig=T)
                    yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,nfert)
                    yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,nfert,irrig=T)
                    yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],nfert)
                    yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],nfert,irrig=T)
                  }
                  yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
                  yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
                  yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
                  yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
                  yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
                  yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
                  yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
                  yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
                  yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
                  yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
                  fn <- paste0(path.ryield6,"A0/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
                  save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
                       yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
                       file=fn)
                  opath <- paste0(path.ryield6,"ncdf/A0/",ssps[ssp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
                              comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
                              comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                }
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A1.nc4")
                if(file.exists(fn)){
                  cat("computing yields for",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff_rf <- readmap.nc(fn,"K_rf")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
                    yields_rf <- yields_ir <- array(NA,dim=dim(dT))
                  #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
                  yields_rf_base <- emulator(coeff_rf,360,0,1,nfert)
                  yields_ir_base <- emulator(coeff_ir,360,0,1,nfert,irrig=T)
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],nfert)
                    yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],nfert,irrig=T)
                    yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],nfert)
                    yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],nfert,irrig=T)
                    yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,nfert)
                    yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,nfert,irrig=T)
                    yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],nfert)
                    yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],nfert,irrig=T)
                  }
                  yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
                  yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
                  yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
                  yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
                  yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
                  yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
                  yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
                  yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
                  yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
                  yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
                  fn <- paste0(path.ryield6,"A1/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow_v",version,".Rdata")
                  save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
                       yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
                       file=fn)
                  opath <- paste0(path.ryield6,"ncdf/A1/",ssps[ssp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
                              comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
                              comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                } # else do nothing
              }
              if(do.lowinput){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A0.nc4")
                if(file.exists(fn)){
                  cat("computing yields for",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff_rf <- readmap.nc(fn,"K_rf")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
                    yields_rf <- yields_ir <- array(NA,dim=dim(dT))
                  #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
                  yields_rf_base <- emulator(coeff_rf,360,0,1,10)
                  yields_ir_base <- emulator(coeff_ir,360,0,1,10,irrig=T)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],10)
                    yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],10,irrig=T)
                    yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],10)
                    yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],10,irrig=T)
                    yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,10)
                    yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,10,irrig=T)
                    yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],10)
                    yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],10,irrig=T)
                  }
                  yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
                  yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
                  yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
                  yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
                  yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
                  yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
                  yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
                  yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
                  yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
                  yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
                  fn <- paste0(path.ryield6,"A0_N10/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_movingwindow_v",version,".Rdata")
                  save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
                       yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
                       file=fn)
                  opath <- paste0(path.ryield6,"ncdf/A0_N10/",ssps[ssp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N10_emulated_yield_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
                              comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
                              comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                  
                }
              }
              if(do.highinput){
                fn <- paste0(path.coeff,ggcms[ggcm],"_",crops.param[cr],"_ggcmi_phase2_emulator_A0.nc4")
                if(file.exists(fn)){
                  cat("computing yields for",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"\n")
                  coeff_rf <- readmap.nc(fn,"K_rf")
                  coeff_ir <- readmap.nc(fn,"K_ir")
                  yields_rf_constt <- yields_ir_constt <- yields_rf_constp <- yields_ir_constp <- yields_rf_constc <- yields_ir_constc <- 
                    yields_rf <- yields_ir <- array(NA,dim=dim(dT))
                  #yields_rf_base <- yields_ir_base <- array(NA,dim=c(dim(tas.base)))
                  yields_rf_base <- emulator(coeff_rf,360,0,1,200)
                  yields_ir_base <- emulator(coeff_ir,360,0,1,200,irrig=T)
                  #for(i in 31:dim(tas.gs)[3]){
                  #for(i in 31:(dim(tas.gs)[3]-15)){
                  for(i in 1:dim(dT)[3]){
                    yields_rf[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],dP[,,i],200)
                    yields_ir[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],dP[,,i],200,irrig=T)
                    yields_rf_constt[,,i] <- emulator(coeff_rf,co2s[,,i],0,dP[,,i],200)
                    yields_ir_constt[,,i] <- emulator(coeff_ir,co2s[,,i],0,dP[,,i],200,irrig=T)
                    yields_rf_constp[,,i] <- emulator(coeff_rf,co2s[,,i],dT[,,i],1,200)
                    yields_ir_constp[,,i] <- emulator(coeff_ir,co2s[,,i],dT[,,i],1,200,irrig=T)
                    yields_rf_constc[,,i] <- emulator(coeff_rf,360,dT[,,i],dP[,,i],200)
                    yields_ir_constc[,,i] <- emulator(coeff_ir,360,dT[,,i],dP[,,i],200,irrig=T)
                  }
                  yields_rf[yields_rf<0 & !is.na(yields_rf)] <- 0
                  yields_ir[yields_ir<0 & !is.na(yields_ir)] <- 0
                  yields_rf_base[yields_rf_base<0 & !is.na(yields_rf_base)] <- 0
                  yields_ir_base[yields_ir_base<0 & !is.na(yields_ir_base)] <- 0
                  yields_rf_constt[yields_rf_constt<0 & !is.na(yields_rf_constt)] <- 0
                  yields_ir_constt[yields_ir_constt<0 & !is.na(yields_ir_constt)] <- 0
                  yields_rf_constp[yields_rf_constp<0 & !is.na(yields_rf_constp)] <- 0
                  yields_ir_constp[yields_ir_constp<0 & !is.na(yields_ir_constp)] <- 0
                  yields_rf_constc[yields_rf_constc<0 & !is.na(yields_rf_constc)] <- 0
                  yields_ir_constc[yields_ir_constc<0 & !is.na(yields_ir_constc)] <- 0
                  fn <- paste0(path.ryield6,"A0_N200/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_movingwindow_v",version,".Rdata")
                  save(yields_rf,yields_ir,yields_rf_constt,yields_ir_constt,yields_rf_constp,yields_ir_constp,
                       yields_rf_constc,yields_ir_constc,yields_rf_base,yields_ir_base,
                       file=fn)
                  opath <- paste0(path.ryield6,"ncdf/A0_N200/",ssps[ssp],"/",ggcms[ggcm],"/")
                  dir.create(opath,showWarnings=F,recursive=T)
                  fn2 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_movingwindow_2011_2084_v",version,".nc4")
                  fn3 <- paste0(opath,"cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_N200_emulated_yield_baseline_1980_2010_average_v",version,".nc4")
                  writemap.nc(fn2,list(yields_rf[,360:1,],yields_ir[,360:1,]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for CMIP5 growing season changes",
                              comment2="emulated crop yields, computed based on differences between the 1980-2010 baseline and a moving 31-year average",
                              start.year=2011,end.year=2084)
                  writemap.nc(fn3,list(yields_rf_base[,360:1],yields_ir_base[,360:1]),var=c("yield_rf","yield_ir"),var2=c("rainfed yield","irrigated yield"),crop=crops.gs2[cr],cropf=crops.nice[cr],
                              units="tDM/ha",title="emulated crop yields for 1980-2010 baseline",
                              comment2="emulated crop yields, computed based on the 1980-2010 baseline 31-year average",
                              start.year=1995,end.year=1995)
                }
              }
            }
          } # else do nothing
        }
        if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n",job[["messages"]]))
      }
    }
  }
} #if(T)

# collect last time slice per pixel for all GGCMxGCMs ####
if(F){
  rcp <- 1
  gcm <- 24
  ggcm <- 3
  cr <- 5
  do.A0 <- TRUE
  do.A1 <- TRUE
  if(do.cmip5){
    if(do.A0){
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        mask_rf <- buf_rf #mrf
        mask_ir <- buf_ir #mir
        
        
        for(rcp in 1:length(rcps)){
          global.delta <- array(0,dim=c(length(ggcms),length(gcms5),720,360))
          weighted.delta <- array(0,dim=c(length(ggcms),length(gcms5),720,360))
          #for(rcp in 4){
          for(gcm in 1:length(gcms5)){
            #parloop <- foreach(gcm=c(1:length(gcms)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 23:25){
            cat("processing",gcms5[gcm],"for",rcps[rcp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield5,"A0/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf)[3]!=74){
                  cat("WARNING, time series too short in",gcms5[gcm],ggcms[ggcm],rcps[rcp],":",dim(yields_rf),"\n")
                  next
                }
                # only keep last layer
                yields_rf_base <- yields_rf[,,1]
                yields_ir_base <- yields_ir[,,1]
                yields_rf <- yields_rf[,,dim(yields_rf)[3]]
                yields_ir <- yields_ir[,,dim(yields_ir)[3]]
                yields_rf[yields_rf_base<0.1] <- NA
                yields_ir[yields_ir_base<0.1] <- NA
                yields_rf_base[yields_rf_base<0.1] <- NA
                yields_ir_base[yields_ir_base<0.1] <- NA
                # production in kcal
                p1 <- (yields_rf_base*mask_rf+yields_ir_base*mask_ir)*ENERG_DM[cr]
                p2 <- (yields_rf*mask_rf+yields_ir*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1
                w1[!is.finite(w1)] <- NA
                prod <- p2/p1
                prod[!is.finite(prod)] <- NA
                global.delta[ggcm,gcm,,] <- prod
                wprod <- p2*w1
                wprod[!is.finite(wprod)] <- NA
                weighted.delta[ggcm,gcm,,] <- wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]-global.prod[cr,ggcm,gcm,rcp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                #fn <- paste0(path.ryield,"global_production_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms[gcm],"_movingwindow.Rdata")
                #save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }#ggcms
          }#gcms
          fn <- paste0(path.rprod5,"A0/cmip5_","global_production_map_A0_",crops.param[cr],"_",rcps[rcp],"_movingwindow_v",version,".Rdata")
          save(global.delta,weighted.delta,file=fn)
        }
      }
      
      
    }#fi(A0)
    
    if(do.A1){
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        mask_rf <- buf_rf #mrf
        mask_ir <- buf_ir #mir
        
        
        for(rcp in 1:length(rcps)){
          global.delta <- array(0,dim=c(length(ggcms),length(gcms5),720,360))
          weighted.delta <- array(0,dim=c(length(ggcms),length(gcms5),720,360))
          #for(rcp in 4){
          for(gcm in 1:length(gcms5)){
            #parloop <- foreach(gcm=c(1:length(gcms)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 23:25){
            cat("processing",gcms5[gcm],"for",rcps[rcp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              fn <- paste0(path.ryield5,"A0/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf)[3]!=74){
                  cat("WARNING, time series too short in",gcms[gcm],ggcms[ggcm],rcps[rcp],":",dim(yields_rf),"\n")
                  next
                }
                # only keep last layer
                yields_rf_base <- yields_rf[,,1]
                yields_ir_base <- yields_ir[,,1]
                yields_rf <- yields_rf[,,dim(yields_rf)[3]]
                yields_ir <- yields_ir[,,dim(yields_ir)[3]]
                yields_rf[yields_rf_base<0.1] <- NA
                yields_ir[yields_ir_base<0.1] <- NA
                yields_rf_base[yields_rf_base<0.1] <- NA
                yields_ir_base[yields_ir_base<0.1] <- NA
                # production in kcal
                p1 <- (yields_rf_base*mask_rf+yields_ir_base*mask_ir)*ENERG_DM[cr]
                p2 <- (yields_rf*mask_rf+yields_ir*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1
                w1[!is.finite(w1)] <- NA
                prod <- p2/p1
                prod[!is.finite(prod)] <- NA
                global.delta[ggcm,gcm,,] <- prod
                wprod <- p2*w1
                wprod[!is.finite(wprod)] <- NA
                weighted.delta[ggcm,gcm,,] <- wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]-global.prod[cr,ggcm,gcm,rcp,]+prod
                # make sure to not recycle things
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                #fn <- paste0(path.ryield,"global_production_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms[gcm],"_movingwindow.Rdata")
                #save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }#ggcms
          }#gcms
          fn <- paste0(path.rprod5,"A1/cmip5_","global_production_map_A1_",crops.param[cr],"_",rcps[rcp],"_movingwindow_v",version,".Rdata")
          save(global.delta,weighted.delta,file=fn)
        }
      }
    }#fi(A1)
  }
  if(do.cmip6)
  {
    if(do.A0){
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        mask_rf <- buf_rf #mrf
        mask_ir <- buf_ir #mir
        
        
        for(ssp in 1:length(ssps)){
          global.delta <- array(0,dim=c(length(ggcms),length(gcms6),720,360))
          weighted.delta <- array(0,dim=c(length(ggcms),length(gcms6),720,360))
          #for(ssp in 4){
          for(gcm in 1:length(gcms6)){
          #parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 23:25){
            cat("processing",gcms6[gcm],"for",ssps[ssp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield6,"A0/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf)[3]!=74){
                  cat("WARNING, time series too short in",gcms6[gcm],ggcms[ggcm],ssps[ssp],":",dim(yields_rf),"\n")
                  next
                }
                # only keep last layer
                yields_rf_base <- yields_rf[,,1]
                yields_ir_base <- yields_ir[,,1]
                yields_rf <- yields_rf[,,dim(yields_rf)[3]]
                yields_ir <- yields_ir[,,dim(yields_ir)[3]]
                yields_rf[yields_rf_base<0.1] <- NA
                yields_ir[yields_ir_base<0.1] <- NA
                yields_rf_base[yields_rf_base<0.1] <- NA
                yields_ir_base[yields_ir_base<0.1] <- NA
                # production in kcal
                p1 <- (yields_rf_base*mask_rf+yields_ir_base*mask_ir)*ENERG_DM[cr]
                p2 <- (yields_rf*mask_rf+yields_ir*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1
                w1[!is.finite(w1)] <- NA
                prod <- p2/p1
                prod[!is.finite(prod)] <- NA
                global.delta[ggcm,gcm,,] <- prod
                wprod <- p2*w1
                wprod[!is.finite(wprod)] <- NA
                weighted.delta[ggcm,gcm,,] <- wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]-global.prod[cr,ggcm,gcm,ssp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                #fn <- paste0(path.ryield,"global_production_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow.Rdata")
                #save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }#ggcms
          }#gcms6
          fn <- paste0(path.rprod6,"A0/","cmip6_global_production_map_A0_",crops.param[cr],"_",ssps[ssp],"_movingwindow_v",version,".Rdata")
          save(global.delta,weighted.delta,file=fn)
        }
      }
      
      
    }#fi(A0)
    
    if(do.A1){
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        mask_rf <- buf_rf #mrf
        mask_ir <- buf_ir #mir
        
        
        for(ssp in 1:length(ssps)){
          global.delta <- array(0,dim=c(length(ggcms),length(gcms6),720,360))
          weighted.delta <- array(0,dim=c(length(ggcms),length(gcms6),720,360))
          #for(ssp in 4){
          for(gcm in 1:length(gcms6)){
          #parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 23:25){
            cat("processing",gcms6[gcm],"for",ssps[ssp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              #fn <- paste0(path.ryield,"A0/",ssps[ssp],"_",gcms6[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow.Rdata")
              fn <- paste0(path.ryield6,"A0/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf)[3]!=74){
                  cat("WARNING, time series too short in",gcms6[gcm],ggcms[ggcm],ssps[ssp],":",dim(yields_rf),"\n")
                  next
                }
                # only keep last layer
                yields_rf_base <- yields_rf[,,1]
                yields_ir_base <- yields_ir[,,1]
                yields_rf <- yields_rf[,,dim(yields_rf)[3]]
                yields_ir <- yields_ir[,,dim(yields_ir)[3]]
                yields_rf[yields_rf_base<0.1] <- NA
                yields_ir[yields_ir_base<0.1] <- NA
                yields_rf_base[yields_rf_base<0.1] <- NA
                yields_ir_base[yields_ir_base<0.1] <- NA
                # production in kcal
                p1 <- (yields_rf_base*mask_rf+yields_ir_base*mask_ir)*ENERG_DM[cr]
                p2 <- (yields_rf*mask_rf+yields_ir*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1
                w1[!is.finite(w1)] <- NA
                prod <- p2/p1
                prod[!is.finite(prod)] <- NA
                global.delta[ggcm,gcm,,] <- prod
                wprod <- p2*w1
                wprod[!is.finite(wprod)] <- NA
                weighted.delta[ggcm,gcm,,] <- wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]-global.prod[cr,ggcm,gcm,ssp,]+prod
                # make sure to not recycle things
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                #fn <- paste0(path.ryield,"global_production_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow.Rdata")
                #save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }#ggcms
          }#gcms6
          fn <- paste0(path.rprod6,"A1/","cmip6_global_production_map_A1_",crops.param[cr],"_",ssps[ssp],"_movingwindow_v",version,".Rdata")
          save(global.delta,weighted.delta,file=fn)
        }
      }
    }#fi(A1)
  }
}

# don't do this??? old legacy code?
# plot variance shares as map ####
if(F){
  if(do.cmip5){
    if(do.A0){
      for(cr in 1:lenght(crops.param)){
        #for(cr in 3:5){
        for(rcp in 1:length(rcps)){
          fn <- paste0(path.ryield,"global_production_map_A0_",crops.param[cr],"_",rcps[rcp],"_movingwindow.Rdata")
          if(file.exists(fn)){
            load(fn)
            # local variance is insensitive to weighting
            # totalvar <- apply(weighted.delta,c(3,4),varall)
            # delta_meangcm <- apply(weighted.delta,c(1,3,4),mean,na.rm=T)
            # delta_meanggcm <- apply(weighted.delta,c(2,3,4),mean,na.rm=T)
            # delta_meangcm[!is.finite(delta_meangcm)] <- NA
            # delta_meanggcm[!is.finite(delta_meanggcm)] <- NA
            # meangcmvar <- apply(delta_meangcm,c(2,3),varall)
            # meanggcmvar <- apply(delta_meanggcm,c(2,3),varall)
            #dev.set(2)
            #image.plot(meanggcmvar/(meangcmvar+meanggcmvar))
            totalvar2 <- apply(global.delta,c(3,4),varall)
            delta_meangcm2 <- apply(global.delta,c(1,3,4),mean,na.rm=T)
            delta_meanggcm2 <- apply(global.delta,c(2,3,4),mean,na.rm=T)
            delta_meangcm2[!is.finite(delta_meangcm2)] <- NA
            delta_meanggcm2[!is.finite(delta_meanggcm2)] <- NA
            meangcmvar2 <- apply(delta_meangcm2,c(2,3),varall)
            meanggcmvar2 <- apply(delta_meanggcm2,c(2,3),varall)
            #dev.set(3)
            png(paste0(path.figs,"variance_shares_GGCMs_map_",crops.param2[cr],"_",rcps[rcp],".png"),width=8*600,height=3.5*600,res=600,pointsize=12)
            par(mar=c(2,3,3,1))
            image.plot(meangcmvar2/(meangcmvar2+meanggcmvar2),xlab="",ylab="",y=seq(-89.75,89.75,length.out=360),x=seq(-179.75,179.75,length.out=720),
                       main=paste("crop model share in total variance",crops.nice[cr],rcps.nice[rcp]),ylim=c(-55,65),asp=1)
            map(add=T,interior=F)
            dev.off()
            
          }
        }
      }
    }
  }
  if(do.cmip6)
  {
    if(do.A0){
      for(cr in 1:lenght(crops.param)){
        #for(cr in 3:5){
        for(rcp in 1:length(rcps)){
          fn <- paste0(path.rprod6,"A0/","cmip6_global_production_map_A0_",crops.param[cr],"_",ssps[ssp],"_movingwindow.Rdata")
          #fn <- paste0(path.ryield,"global_production_map_A0_",crops.param[cr],"_",rcps[rcp],"_movingwindow.Rdata")
          if(file.exists(fn)){
            load(fn)
            # local variance is insensitive to weighting
            # totalvar <- apply(weighted.delta,c(3,4),varall)
            # delta_meangcm <- apply(weighted.delta,c(1,3,4),mean,na.rm=T)
            # delta_meanggcm <- apply(weighted.delta,c(2,3,4),mean,na.rm=T)
            # delta_meangcm[!is.finite(delta_meangcm)] <- NA
            # delta_meanggcm[!is.finite(delta_meanggcm)] <- NA
            # meangcmvar <- apply(delta_meangcm,c(2,3),varall)
            # meanggcmvar <- apply(delta_meanggcm,c(2,3),varall)
            #dev.set(2)
            #image.plot(meanggcmvar/(meangcmvar+meanggcmvar))
            totalvar2 <- apply(global.delta,c(3,4),varall)
            delta_meangcm2 <- apply(global.delta,c(1,3,4),mean,na.rm=T)
            delta_meanggcm2 <- apply(global.delta,c(2,3,4),mean,na.rm=T)
            delta_meangcm2[!is.finite(delta_meangcm2)] <- NA
            delta_meanggcm2[!is.finite(delta_meanggcm2)] <- NA
            meangcmvar2 <- apply(delta_meangcm2,c(2,3),varall)
            meanggcmvar2 <- apply(delta_meanggcm2,c(2,3),varall)
            #dev.set(3)
            png(paste0(path.figs,"cmip6_variance_shares_GGCMs_map_",crops.param2[cr],"_",ssps[ssp],".png"),width=8*600,height=3.5*600,res=600,pointsize=12)
            par(mar=c(2,3,3,1))
            image.plot(meangcmvar2/(meangcmvar2+meanggcmvar2),xlab="",ylab="",y=seq(-89.75,89.75,length.out=360),x=seq(-179.75,179.75,length.out=720),
                       main=paste("crop model share in total variance",crops.nice[cr],rcps.nice[rcp]),ylim=c(-55,65),asp=1)
            map(add=T,interior=F)
            dev.off()
            
          }
        }
      }
    }
  }
}

# don't do this??? old legacy code?
if(F){
  for(rcp in 4){
    for(gcm in 23){
      for(cr in 1){
        for(ggcm in 1:length(ggcms)){
          fn <- paste0(path.ryield,rcps[rcp],"_",gcms[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow.Rdata")
          if(file.exists(fn)){
            load(fn)
            png(paste0(path.figs,"yieldmap_",crops.param[cr],"_rf_",rcps[rcp],"_2084_",ggcms[ggcm],".png"),width=9*600,height=4*600,res=600,pointsize=12)
            #par(mar=c(2,2,3,3))
            par(mar=c(0,0,0,0))
            image.plot(yields_rf[,,74],x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),main=paste(crops.nice[cr],"rainfed"),
                       col=terrain.colors(50)[50:1],ylim=c(-55,65),xlab="",ylab="")
            map(add=T,interior=F)
            dev.off()
            png(paste0(path.figs,"yieldmap_",crops.param[cr],"_ir_",rcps[rcp],"_2084_",ggcms[ggcm],".png"),width=9*600,height=4*600,res=600,pointsize=12)
            #par(mar=c(2,2,3,3))
            par(mar=c(0,0,0,0))
            image.plot(yields_ir[,,74],x=seq(-179.75,179.75,length.out=720),y=seq(-89.75,89.75,length.out=360),main=paste(crops.nice[cr],"rainfed"),
                       col=terrain.colors(50)[50:1],ylim=c(-55,65),xlab="",ylab="")
            map(add=T,interior=F)
            dev.off()
          }
        }
      }
    }
  }
}

# aggregate yields to production
if(F){
  rcp <- 1
  gcm <- 24
  ggcm <- 3
  cr <- 1
  do.constC <- F #TRUE
  do.constT <- do.constP <- do.A0 <- FALSE
  do.A0 <- TRUE
  do.A1 <- TRUE
  
  if(do.cmip5){
    if(do.A0){
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms5),length(rcps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(NA,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(rcp in 1:length(rcps)){
          #for(rcp in 4){
          #for(gcm in 1:length(gcms)){
          parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 23:25){
            cat("processing",gcms5[gcm],"for",rcps[rcp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield5,"A0/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              #fn <- paste0(path.ryield5,"A0/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf)[3]!=74){
                  cat("WARNING, time series too short in",gcms5[gcm],ggcms[ggcm],rcps[rcp],":",dim(yields_rf),"\n")
                  next
                }
                yields_rf[yields_rf<0.1] <- NA
                yields_ir[yields_ir<0.1] <- NA
                tabl[cr,ggcm,gcm,rcp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,rcp] <- prod.base
                global.prod.base[6,ggcm,gcm,rcp] <- global.prod.base[6,ggcm,gcm,rcp] + prod.base
                # production in kcal
                p1 <- (yields_rf*mask_rf+yields_ir*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,rcp,] <- prod
                global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,rcp,] <- wprod
                weighted.prod[6,ggcm,gcm,rcp,] <- weighted.prod[6,ggcm,gcm,rcp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]-global.prod[cr,ggcm,gcm,rcp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod5,"A0/cmip5_","global_production_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
          
        }
      }
      
      #total <- array(NA,dim=c(dim(delta_total_26),2),
      #               dimnames=list(co2=namco2,gcm=namgcm,ggcm=namggcm,rcp=namrcp))
      #total[,,,1] <- delta_total_26
      #total[,,,2] <- delta_total_85
      
      # #Using ANOVA
      # # for all GCM, GGCM, RCP but with co2 only
      # erg <- data.frame(var=c("ggcm","gcm","rcp","ggcm_x_gcm","ggcm_x_rcp","gcm_x_rcp","ggcm_x_gcm_x_rcp","SD"))
      # #global.prod <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms),length(rcps),length(2011:2084)))
      # delta <- global.prod
      # for(i in 1:dim(global.prod)[5])
      #   delta[,,,,i] <- global.prod[,,,,i]/global.prod.base
      # total <- apply(delta,c(1:4),mean,na.rm=T)
      # dimnames(total) <- list(crop=c(crops.param,"all"),ggcm=ggcms,gcm=paste(gcms,param,sep=":"),rcp=rcps)
      # for(cr in 1:(length(crops.param)+1)){
      #   mt <- melt(total[cr,,,])
      #   #ly <- aov(value~gcm+ggcm+co2+rcp,data=mt)
      #   ly <- aov(value~gcm+ggcm+rcp+gcm:rcp+ggcm:rcp+gcm:ggcm+gcm:ggcm:rcp,data=mt)
      #   summary(ly)
      #   ssq <- summary(ly)[[1]][,2]
      #   ssv <- ssq/sum(ssq)*100
      #   ssv <- c(ssv,var(as.vector(total[cr,,,]),na.rm=T)^.5)
      #   names(ssv) <- c(names(mt)[1:3],"gcmxrcp","ggcmxrcp","ggcmxgcm","gcmxggcmxrcp","SD")
      #   erg[[cr]] <- ssv
      # }
      
      
      
      
      # erg2 <- data.frame(var=c("ggcm","gcm","rcp","ggcm_x_gcm","ggcm_x_rcp","gcm_x_rcp","ggcm_x_gcm_x_rcp","SD"))
      # #global.prod <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms),length(rcps),length(2011:2084)))
      # delta <- global.prod
      # for(i in 1:dim(global.prod)[5])
      #   delta[,,,,i] <- global.prod[,,,,i]/global.prod.base
      # total <- apply(delta[,,,,60:74],c(1:4),mean,na.rm=T)
      # dimnames(total) <- list(crop=c(crops.param,"all"),ggcm=ggcms,gcm=paste(gcms,param,sep=":"),rcp=rcps)
      # for(cr in 1:(length(crops.param)+1)){
      #   mt <- melt(total[cr,,,])
      #   #ly <- aov(value~gcm+ggcm+co2+rcp,data=mt)
      #   ly <- aov(value~gcm+ggcm+rcp+gcm:rcp+ggcm:rcp+gcm:ggcm+gcm:ggcm:rcp,data=mt)
      #   summary(ly)
      #   ssq <- summary(ly)[[1]][,2]
      #   ssv <- ssq/sum(ssq)*100
      #   ssv <- c(ssv,var(as.vector(total[cr,,,]),na.rm=T)^.5)
      #   names(ssv) <- c(names(mt)[1:3],"gcmxrcp","ggcmxrcp","ggcmxgcm","gcmxggcmxrcp","SD")
      #   erg2[[cr]] <- ssv
      # }
      #fn <- paste0(path.ryield,"global_production.Rdata")
      #save(tabl,global.prod,global.prod.base,weighted.prod,file=fn)
    }#if do.A0
    if(do.A1){
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms5),length(rcps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(NA,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(rcp in 1:length(rcps)){
          #for(rcp in 4){
          #for(gcm in 1:length(gcms)){
          parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 23:25){
            cat("processing",gcms5[gcm],"for",rcps[rcp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield5,"A1/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf)[3]!=74){
                  cat("WARNING, time series too short in",gcms5[gcm],ggcms[ggcm],rcps[rcp],":",dim(yields_rf),"\n")
                  next
                }
                yields_rf[yields_rf<0.1] <- NA
                yields_ir[yields_ir<0.1] <- NA
                tabl[cr,ggcm,gcm,rcp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,rcp] <- prod.base
                global.prod.base[6,ggcm,gcm,rcp] <- global.prod.base[6,ggcm,gcm,rcp] + prod.base
                # production in kcal
                p1 <- (yields_rf*mask_rf+yields_ir*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,rcp,] <- prod
                global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,rcp,] <- wprod
                weighted.prod[6,ggcm,gcm,rcp,] <- weighted.prod[6,ggcm,gcm,rcp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]-global.prod[cr,ggcm,gcm,rcp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod5,"A1/cmip5_","global_production_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_A1_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
          
        }
      }
      
    }#if do.A1
    
    # no-co2 case ####
    if(do.constC){
      
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms5),length(rcps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # in order to add up things, NAs need to be replaced with zeros
        mrf[!is.finite(mrf)] <- 0
        mir[!is.finite(mir)] <- 0
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        prod05[!is.finite(prod05)] <- 0
        area05[!is.finite(area05)] <- 0
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(0,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        #buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_rf[(mrf+mir)<1] <- area05[(mrf+mir)<1]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        #buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        buf_ir[(mrf+mir)<1] <- 0
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(rcp in 1:length(rcps)){
          #for(rcp in 4){
          parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 1:length(gcms)){
            #for(gcm in 23:25){
            cat("processing",gcms5[gcm],"for",rcps[rcp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield5,"A0/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf_constc)[3]!=74){
                  cat("WARNING, time series too short in",gcms5[gcm],ggcms[ggcm],rcps[rcp],":",dim(yields_rf),"\n")
                  next
                }
                yields_rf_constc[yields_rf_constc<0.1] <- NA
                yields_ir_constc[yields_ir_constc<0.1] <- NA
                tabl[cr,ggcm,gcm,rcp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,rcp] <- prod.base
                global.prod.base[6,ggcm,gcm,rcp] <- global.prod.base[6,ggcm,gcm,rcp] + prod.base
                # production in kcal
                p1 <- (yields_rf_constc*mask_rf+yields_ir_constc*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,rcp,] <- prod
                global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,rcp,] <- wprod
                weighted.prod[6,ggcm,gcm,rcp,] <- weighted.prod[6,ggcm,gcm,rcp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]-global.prod[cr,ggcm,gcm,rcp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod5,"A0/cmip5_global_production_constc_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
        }
      }
      
    } #do.constC
    
    # const T case ####
    if(do.constT){
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms5),length(rcps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(NA,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(rcp in 1:length(rcps)){
          #for(rcp in 4){
          parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 1:length(gcms)){
            #for(gcm in 23:25){
            cat("processing",gcms5[gcm],"for",rcps[rcp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield5,"A0/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf_constt)[3]!=74){
                  cat("WARNING, time series too short in",gcms5[gcm],ggcms[ggcm],rcps[rcp],":",dim(yields_rf_constt),"\n")
                  next
                }
                yields_rf_constt[yields_rf_constt<0.1] <- NA
                yields_ir_constt[yields_ir_constt<0.1] <- NA
                tabl[cr,ggcm,gcm,rcp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,rcp] <- prod.base
                global.prod.base[6,ggcm,gcm,rcp] <- global.prod.base[6,ggcm,gcm,rcp] + prod.base
                # production in kcal
                p1 <- (yields_rf_constt*mask_rf+yields_ir_constt*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,rcp,] <- prod
                global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,rcp,] <- wprod
                weighted.prod[6,ggcm,gcm,rcp,] <- weighted.prod[6,ggcm,gcm,rcp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]-global.prod[cr,ggcm,gcm,rcp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod5,"A0/cmip5_global_production_constt_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
        }
      }
      
      #fn <- paste0(path.ryield,"global_production_consttas.Rdata")
      #save(tabl,global.prod,global.prod.base,weighted.prod,file=fn)  
    } # do.constT
    
    # const P case ####
    if(do.constP){
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms5),length(rcps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(NA,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(rcp in 1:length(rcps)){
          #for(rcp in 4){
          parloop <- foreach(gcm=c(1:length(gcms5)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 1:length(gcms)){
            #for(gcm in 23:25){
            cat("processing",gcms5[gcm],"for",rcps[rcp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield5,"A0/cmip5_",rcps[rcp],"_",gcms5[gcm],"_r",run5[gcm],"i1p",param5[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf_constp)[3]!=74){
                  cat("WARNING, time series too short in",gcms5[gcm],ggcms[ggcm],rcps[rcp],":",dim(yields_rf_constp),"\n")
                  next
                }
                yields_rf_constp[yields_rf_constp<0.1] <- NA
                yields_ir_constp[yields_ir_constp<0.1] <- NA
                tabl[cr,ggcm,gcm,rcp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,rcp] <- prod.base
                global.prod.base[6,ggcm,gcm,rcp] <- global.prod.base[6,ggcm,gcm,rcp] + prod.base
                # production in kcal
                p1 <- (yields_rf_constp*mask_rf+yields_ir_constp*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,rcp,] <- prod
                global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,rcp,] <- wprod
                weighted.prod[6,ggcm,gcm,rcp,] <- weighted.prod[6,ggcm,gcm,rcp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]-global.prod[cr,ggcm,gcm,rcp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod5,"A0/cmip5_global_production_constp_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
        }
      }
      
      #fn <- paste0(path.ryield,"global_production_constpr.Rdata")
      #save(tabl,global.prod,global.prod.base,weighted.prod,file=fn)  
    } # do.constP
  }  
  if(do.cmip6)
  {
    if(do.A0){
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms6),length(ssps)))
      for(cr in 1:length(crops.param)){
      #for(cr in 5){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(NA,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(ssp in 1:length(ssps)){
          #for(ssp in 4){
          #for(gcm in 1:length(gcms6)){
          #parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
          for(gcm in c(10,29)){
            cat("processing",gcms6[gcm],"for",ssps[ssp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              #fn <- paste0(path.ryield,"A0/",ssps[ssp],"_",gcms6[gcm],"_r",run[gcm],"i1p",param[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow.Rdata")
              fn <- paste0(path.ryield6,"A0/","cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf)[3]!=74){
                  cat("WARNING, time series too short in",gcms6[gcm],ggcms[ggcm],ssps[ssp],":",dim(yields_rf),"\n")
                  next
                }
                yields_rf[yields_rf<0.1] <- NA
                yields_ir[yields_ir<0.1] <- NA
                tabl[cr,ggcm,gcm,ssp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,ssp] <- prod.base
                global.prod.base[6,ggcm,gcm,ssp] <- global.prod.base[6,ggcm,gcm,ssp] + prod.base
                # production in kcal
                p1 <- (yields_rf*mask_rf+yields_ir*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,ssp,] <- prod
                global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,ssp,] <- wprod
                weighted.prod[6,ggcm,gcm,ssp,] <- weighted.prod[6,ggcm,gcm,ssp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]-global.prod[cr,ggcm,gcm,ssp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod6,"A0/","cmip6_global_production_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
          
        }
      }
      
    }#if do.A0
    if(do.A1){
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms6),length(ssps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(NA,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(ssp in 1:length(ssps)){
          #for(ssp in 4){
          #for(gcm in 1:length(gcms6)){
          parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 23:25){
            cat("processing",gcms6[gcm],"for",ssps[ssp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield6,"A1/cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A1_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf)[3]!=74){
                  cat("WARNING, time series too short in",gcms6[gcm],ggcms[ggcm],ssps[ssp],":",dim(yields_rf),"\n")
                  next
                }
                yields_rf[yields_rf<0.1] <- NA
                yields_ir[yields_ir<0.1] <- NA
                tabl[cr,ggcm,gcm,ssp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,ssp] <- prod.base
                global.prod.base[6,ggcm,gcm,ssp] <- global.prod.base[6,ggcm,gcm,ssp] + prod.base
                # production in kcal
                p1 <- (yields_rf*mask_rf+yields_ir*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,ssp,] <- prod
                global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,ssp,] <- wprod
                weighted.prod[6,ggcm,gcm,ssp,] <- weighted.prod[6,ggcm,gcm,ssp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]-global.prod[cr,ggcm,gcm,ssp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod6,"A1/","cmip6_global_production_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_A1_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
          
        }
      }
      
    }#if do.A1
    
    # no-co2 case ####
    if(do.constC){
      
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms6),length(ssps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # in order to add up things, NAs need to be replaced with zeros
        mrf[!is.finite(mrf)] <- 0
        mir[!is.finite(mir)] <- 0
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        prod05[!is.finite(prod05)] <- 0
        area05[!is.finite(area05)] <- 0
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(0,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        #buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_rf[(mrf+mir)<1] <- area05[(mrf+mir)<1]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        #buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        buf_ir[(mrf+mir)<1] <- 0
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(ssp in 1:length(ssps)){
          #for(ssp in 4){
          parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 1:length(gcms6)){
            #for(gcm in 23:25){
            cat("processing",gcms6[gcm],"for",ssps[ssp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield6,"A0/cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf_constc)[3]!=74){
                  cat("WARNING, time series too short in",gcms6[gcm],ggcms[ggcm],ssps[ssp],":",dim(yields_rf),"\n")
                  next
                }
                yields_rf_constc[yields_rf_constc<0.1] <- NA
                yields_ir_constc[yields_ir_constc<0.1] <- NA
                tabl[cr,ggcm,gcm,ssp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,ssp] <- prod.base
                global.prod.base[6,ggcm,gcm,ssp] <- global.prod.base[6,ggcm,gcm,ssp] + prod.base
                # production in kcal
                p1 <- (yields_rf_constc*mask_rf+yields_ir_constc*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,ssp,] <- prod
                global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,ssp,] <- wprod
                weighted.prod[6,ggcm,gcm,ssp,] <- weighted.prod[6,ggcm,gcm,ssp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]-global.prod[cr,ggcm,gcm,ssp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod6,"A0/","cmip6_global_production_constc_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
        }
      }
      
    } #do.constC
    
    # const T case ####
    if(do.constT){
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms6),length(ssps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(NA,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(ssp in 1:length(ssps)){
          #for(ssp in 4){
          parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 1:length(gcms6)){
            #for(gcm in 23:25){
            cat("processing",gcms6[gcm],"for",ssps[ssp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield6,"A0/cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf_constt)[3]!=74){
                  cat("WARNING, time series too short in",gcms6[gcm],ggcms[ggcm],ssps[ssp],":",dim(yields_rf_constt),"\n")
                  next
                }
                yields_rf_constt[yields_rf_constt<0.1] <- NA
                yields_ir_constt[yields_ir_constt<0.1] <- NA
                tabl[cr,ggcm,gcm,ssp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,ssp] <- prod.base
                global.prod.base[6,ggcm,gcm,ssp] <- global.prod.base[6,ggcm,gcm,ssp] + prod.base
                # production in kcal
                p1 <- (yields_rf_constt*mask_rf+yields_ir_constt*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,ssp,] <- prod
                global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,ssp,] <- wprod
                weighted.prod[6,ggcm,gcm,ssp,] <- weighted.prod[6,ggcm,gcm,ssp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]-global.prod[cr,ggcm,gcm,ssp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod6,"A0/","cmip6_global_production_constt_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
        }
      }
      
      #fn <- paste0(path.ryield,"global_production_consttas.Rdata")
      #save(tabl,global.prod,global.prod.base,weighted.prod,file=fn)  
    } # do.constT
    
    # const P case ####
    if(do.constP){
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps)))
      tabl <- array(NA,dim=c(length(crops.param),length(ggcms),length(gcms6),length(ssps)))
      for(cr in 1:length(crops.param)){
        #for(cr in 1){
        if(cr<4){
          fn <- paste0(path.lu,crops.param[cr],".nc4")
          mrf <- readmap.nc(fn,"rainfed")
          mir <- readmap.nc(fn,"irrigated")
          
        } else {
          fn <- paste0(path.gs,"phase2.masks/winter_and_spring_wheat_areas_v1_180627.nc4")  
          mrf <- readmap.nc(fn,if(cr==4) "wwh_rf_area" else if(cr==5) "swh_rf_area" else "rainfed")
          mir <- readmap.nc(fn,if(cr==4) "wwh_ir_area" else if(cr==5) "swh_ir_area" else "irrigated")
        }
        # load Monfreada prodcution for equal weighting
        load(paste0(path.prod,crops.fert[cr],"/production_harvestedarea_yield.Rdata"))
        prod05 <- prod05[,360:1] # invert latitudes to play with other maps
        area05 <- area05[,360:1] # invert latitudes to play with other maps
        
        # add 3rd dimension for easier processing later on
        mask_rf <- mask_ir <- array(NA,dim=c(dim(mrf),length(2011:2084)))
        buf_rf <- area05*mrf/(mrf+mir)
        # assume all rainfed for areas with Monfreda production but no MIRCA areas
        buf_rf[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)] <- area05[(mrf+mir)<1 & !is.na(area05) & !is.na(mrf+mir)]
        buf_ir <- area05*mir/(mrf+mir)
        # assume no irrigated production for areas with 
        buf_ir[(mrf+mir)<1 & area05>0 & !is.na(area05) & !is.na(mrf+mir)] <- 0 
        for(i in 1:length(2011:2084)){
          mask_rf[,,i] <- buf_rf #mrf
          mask_ir[,,i] <- buf_ir #mir
        }
        
        for(ssp in 1:length(ssps)){
          #for(ssp in 4){
          parloop <- foreach(gcm=c(1:length(gcms6)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
            #for(gcm in 1:length(gcms6)){
            #for(gcm in 23:25){
            cat("processing",gcms6[gcm],"for",ssps[ssp],"and",crops.param[cr],"\n")
            for(ggcm in 1:length(ggcms)){
              #for(ggcm in 3){
              # make sure to not recycle things
              fn <- paste0(path.ryield6,"A0/cmip6_",ssps[ssp],"_",gcms6[gcm],"_r",run6[gcm],"i1p",param6[gcm],"_",crops.param[cr],"_",ggcms[ggcm],"_A0_emulated_yield_movingwindow_v",version,".Rdata")
              if(file.exists(fn))
              {
                load(fn)
                if(dim(yields_rf_constp)[3]!=74){
                  cat("WARNING, time series too short in",gcms6[gcm],ggcms[ggcm],ssps[ssp],":",dim(yields_rf_constp),"\n")
                  next
                }
                yields_rf_constp[yields_rf_constp<0.1] <- NA
                yields_ir_constp[yields_ir_constp<0.1] <- NA
                tabl[cr,ggcm,gcm,ssp] <- "X"
                prod.base <- sum(yields_rf_base*mrf+yields_ir_base*mir,na.rm=T)*ENERG_DM[cr]
                global.prod.base[cr,ggcm,gcm,ssp] <- prod.base
                global.prod.base[6,ggcm,gcm,ssp] <- global.prod.base[6,ggcm,gcm,ssp] + prod.base
                # production in kcal
                p1 <- (yields_rf_constp*mask_rf+yields_ir_constp*mask_ir)*ENERG_DM[cr]
                w1 <- prod05*ENERG[cr]/p1[,,1]
                #w1[prod05<100] <- NA
                #w1[(mrf+mir)<100] <- NA
                w1[!is.finite(w1)] <- NA
                prod <- apply((p1),3,sum,na.rm=T)
                global.prod[cr,ggcm,gcm,ssp,] <- prod
                global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]+prod
                delt <- array(NA,dim=dim(yields_rf))
                for(i in 1:dim(delt)[3])
                {
                  delt[,,i] <- p1[,,i]*w1
                }
                delt[!is.finite(delt)] <- NA
                wprod <- apply(delt,3,sum,na.rm=T)
                weighted.prod[cr,ggcm,gcm,ssp,] <- wprod
                weighted.prod[6,ggcm,gcm,ssp,] <- weighted.prod[6,ggcm,gcm,ssp,] + wprod
                # for updating things: remove old value first
                #global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]-global.prod[cr,ggcm,gcm,ssp,]+prod
                rm(yields_rf,yields_ir,yields_rf_base,yields_ir_base,
                   yields_rf_constc,yields_rf_constt,yields_rf_constp,
                   yields_ir_constc,yields_ir_constt,yields_ir_constp)
                fn <- paste0(path.rprod6,"A0/","cmip6_global_production_constp_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow_v",version,".Rdata")
                save(wprod,prod,prod.base,file=fn)
                
              } else {
                next
              }
            }
          }
          if(parallel) for(job in parloop) cat(paste(str(job[["messages"]]),"\n\n",job[["messages"]]))
        }
      }
      
      #fn <- paste0(path.ryield,"global_production_constpr.Rdata")
      #save(tabl,global.prod,global.prod.base,weighted.prod,file=fn)  
    } # do.constP
  }
}

# pre-processing: combinations of 5 GCMs dropped ####
if(F){
  if(do.cmip5){
    #leave ~10% of GCMs out (equivalent to 1 of 9 GGCMs)
    index5 <- combn(c(1:9,12:46),5) # GCMs 10 and 11 without valid members, 5 of 44 is 0.11363, similar to 1 GGCM of 9 (0.11111)
    sep <- 10000
    global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
    weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
    global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps)))
    for(cr in 1:length(crops.param)){
      for(rcp in 1:length(rcps)){
        for(ggcm in 1:length(ggcms)){
          for(gcm in 1:length(gcms5)){
            fn <- paste0(path.ryield5,"global_production_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_movingwindow.Rdata")
            if(file.exists(fn)){
              load(fn)
              global.prod[cr,ggcm,gcm,rcp,] <- prod
              weighted.prod[cr,ggcm,gcm,rcp,] <- wprod
              global.prod.base[cr,ggcm,gcm,rcp] <- prod.base
              global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]+prod
              weighted.prod[6,ggcm,gcm,rcp,] <- weighted.prod[6,ggcm,gcm,rcp,]+wprod
              global.prod.base[6,ggcm,gcm,rcp] <- global.prod.base[6,ggcm,gcm,rcp]+prod.base
              #} else{
              #cat("missing",fn,"\n")
            }
          }
        }
      }
    }
    # define delta by deviding by first time step to make figures look nicer
    delta <- global.prod
    for(i in 1:dim(global.prod)[5])
      #delta[,,,,i] <- global.prod[,,,,i]/global.prod.base
      delta[,,,,i] <- global.prod[,,,,i]/global.prod[,,,,1]
    deltaw <- weighted.prod
    for(i in 1:dim(weighted.prod)[5])
      deltaw[,,,,i] <- weighted.prod[,,,,i]/weighted.prod[,,,,1]
    
    delta[!is.finite(delta)] <- NA
    deltaw[!is.finite(deltaw)] <- NA
    cat("done computing deltas",str(delta),range(delta,na.rm=T),"index",str(index5),"\n")
    #leave 5GCMs out
    parloop <- foreach(bit=c(1:ceiling(dim(index5)[2]/sep)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
      #drop5gcm <- list()
      #drop5gcm_meangcm <- list()
      #drop5gcm_meanggcm <- list()
      # do this for maize and RCP85 only
      se <- c(((bit-1)*sep+1):min(dim(index5)[2],bit*sep))
      mi <- list(min=999,tvar=NULL,meangcmvar=NULL,meanggcmvar=NULL,ind=NULL)
      ma <- list(max=-999,tvar=NULL,meangcmvar=NULL,meanggcmvar=NULL,ind=NULL)
      for(i in se){
        #if(i < 10 | i %% 1000 ==0)
        #  cat("processing",i,"of",dim(index)[2],"\n")
        delta2 <- delta[1,,-index5[,i],4,]
        if(all(is.na(delta2)))
          next
        #drop5gcm[[i]] <- delta2
        meangcm <- apply(delta2,c(1,3),mean,na.rm=T)
        #drop5gcm_meangcm[[i]] <- meangcm
        meanggcm <- apply(delta2,c(2,3),mean,na.rm=T)
        #drop5gcm_meanggcm[[i]] <- meanggcm
        tvar <- mgcmvar <- mggcmvar <- vector("numeric",74)
        for(j in 1:74){
          tvar[j]=var(as.vector(delta2[,,j]),na.rm=T)  
          mgcmvar[j]=var(as.vector(meangcm[,j]),na.rm=T)  
          mggcmvar[j]=var(as.vector(meanggcm[,j]),na.rm=T)  
        }
        if(mggcmvar[74]<mi$min){
          mi$min<-mggcmvar[74]
          mi$tvar<-tvar
          mi$meangcmvar<-mgcmvar
          mi$meanggcmvar<-mggcmvar
          mi$ind<-i
        } else if(mggcmvar[74]>ma$max){
          ma$max<-mggcmvar[74]
          ma$tvar<-tvar
          ma$meangcmvar<-mgcmvar
          ma$meanggcmvar<-mggcmvar
          ma$ind<-i
        }
        
      }
      fn <- paste0(path.ryield5,"drop5gcms_set_",bit,"-tenthousand.Rdata")
      #save(se,index,delta2,drop5gcm,drop5gcm_meangcm,drop5gcm_meanggcm,file=fn)
      save(mi,ma,file=fn)
    }
  }
  if(do.cmip6)
  {
    index6 <- combn(c(1:length(gcms6)),3) # GCMs 10 and 11 without valid members, 4 of 29 is 0.1379, 3 of 29 is 0.10344
    sep <- 10000
    global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
    weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
    global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps)))
    for(cr in 1:length(crops.param)){
      for(ssp in 1:length(ssps)){
        for(ggcm in 1:length(ggcms)){
          for(gcm in 1:length(gcms6)){
            fn <- paste0(path.rprod6,"A0/cmip6_global_production_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow.Rdata")
            if(file.exists(fn)){
              load(fn)
              global.prod[cr,ggcm,gcm,ssp,] <- prod
              weighted.prod[cr,ggcm,gcm,ssp,] <- wprod
              global.prod.base[cr,ggcm,gcm,ssp] <- prod.base
              global.prod[6,ggcm,gcm,ssp,] <- global.prod[6,ggcm,gcm,ssp,]+prod
              weighted.prod[6,ggcm,gcm,ssp,] <- weighted.prod[6,ggcm,gcm,ssp,]+wprod
              global.prod.base[6,ggcm,gcm,ssp] <- global.prod.base[6,ggcm,gcm,ssp]+prod.base
            } else{
              cat("missing",fn,"\n")
            }
          }
        }
      }
    }
    # define delta by deviding by first time step to make figures look nicer
    delta <- global.prod
    for(i in 1:dim(global.prod)[5])
      #delta[,,,,i] <- global.prod[,,,,i]/global.prod.base
      delta[,,,,i] <- global.prod[,,,,i]/global.prod[,,,,1]
    deltaw <- weighted.prod
    for(i in 1:dim(weighted.prod)[5])
      deltaw[,,,,i] <- weighted.prod[,,,,i]/weighted.prod[,,,,1]
    
    delta[!is.finite(delta)] <- NA
    deltaw[!is.finite(deltaw)] <- NA
    cat("done computing deltas",str(delta),range(delta,na.rm=T),"index",str(index6),"\n")
    #leave 5GCMs out
    #parloop <- foreach(bit=c(1:ceiling(dim(index6)[2]/sep)), .errorhandling='pass', .verbose=FALSE, .export='message') %dopar% {
    #drop5gcm <- list()
    #drop5gcm_meangcm <- list()
    #drop5gcm_meanggcm <- list()
    # do this for maize and ssp85 only
    mi <- list(min=999,tvar=NULL,meangcmvar=NULL,meanggcmvar=NULL,ind=NULL)
    ma <- list(max=-999,tvar=NULL,meangcmvar=NULL,meanggcmvar=NULL,ind=NULL)
    for(i in 1:dim(index6)[2]){
      #if(i < 10 | i %% 1000 ==0)
      #  cat("processing",i,"of",dim(index)[2],"\n")
      delta2 <- delta[1,,-index6[,i],which(ssps=="ssp585"),]
      if(all(is.na(delta2)))
        next
      #drop5gcm[[i]] <- delta2
      meangcm <- apply(delta2,c(1,3),mean,na.rm=T)
      #drop5gcm_meangcm[[i]] <- meangcm
      meanggcm <- apply(delta2,c(2,3),mean,na.rm=T)
      #drop5gcm_meanggcm[[i]] <- meanggcm
      tvar <- mgcmvar <- mggcmvar <- vector("numeric",74)
      for(j in 1:74){
        tvar[j]=var(as.vector(delta2[,,j]),na.rm=T)  
        mgcmvar[j]=var(as.vector(meangcm[,j]),na.rm=T)  
        mggcmvar[j]=var(as.vector(meanggcm[,j]),na.rm=T)  
      }
      if(mggcmvar[74]<mi$min){
        mi$min<-mggcmvar[74]
        mi$tvar<-tvar
        mi$meangcmvar<-mgcmvar
        mi$meanggcmvar<-mggcmvar
        mi$ind<-i
      } else if(mggcmvar[74]>ma$max){
        ma$max<-mggcmvar[74]
        ma$tvar<-tvar
        ma$meangcmvar<-mgcmvar
        ma$meanggcmvar<-mggcmvar
        ma$ind<-i
      }
      
    }
    fn <- paste0(path.ryield6,"drop3gcms_set.Rdata")
    #save(se,index,delta2,drop5gcm,drop5gcm_meangcm,drop5gcm_meanggcm,file=fn)
    save(mi,ma,file=fn)
  }
  #}
}

# figures and analyses ####
if(T){
  index5 <- combn(c(1:9,12:46),5) # GCMs 10 and 11 without valid members, 5 of 44 is 0.11363, similar to 1 GGCM of 9 (0.11111)
  index6 <- combn(c(1:length(gcms6)),4) # GCMs 10 and 11 without valid members, 4 of 29 is 0.1379, 3 of 29 is 0.10344
  sep <- 10000
  #rcps <- rcps[-3]
  #rcps.nice <- rcps.nice[-3]
  #rcps.nr <- rcps.nr[-3]
  cat(rcps,"\n")
  for(sets in c("","_constco2","_A1")){
    if(do.cmip5){
      
      #fn <- paste0(path.ryield,"global_production",sets,".Rdata")
      #load(fn)
      global.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      weighted.prod <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps),length(2011:2084)))
      global.prod.base <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms5),length(rcps)))
      for(cr in 1:length(crops.param)){
        tabl2 <- array("",dim=c(length(ggcms),length(gcms5)))
        for(rcp in 1:length(rcps)){
          for(ggcm in 1:length(ggcms)){
            for(gcm in 1:length(gcms5)){
              if(sets=="")
                fn <- paste0(path.rprod5,"A0/cmip5_global_production_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_movingwindow_v",version,".Rdata")
              else if(sets=="_constco2")
                fn <- paste0(path.rprod5,"A0/cmip5_global_production_constc_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_movingwindow_v",version,".Rdata")
              else if(sets=="_A1")
                fn <- paste0(path.rprod5,"A1/cmip5_global_production_",crops.param[cr],"_",rcps[rcp],"_",ggcms[ggcm],"_",gcms5[gcm],"_A1_movingwindow_v",version,".Rdata")
              if(file.exists(fn)){
                tabl2[ggcm,gcm] <- paste(tabl2[ggcm,gcm],rcps.nr[rcp],sep= if(""==tabl2[ggcm,gcm]) "" else ", ")
                load(fn)
                global.prod[cr,ggcm,gcm,rcp,] <- prod
                weighted.prod[cr,ggcm,gcm,rcp,] <- wprod
                global.prod.base[cr,ggcm,gcm,rcp] <- prod.base
                global.prod[6,ggcm,gcm,rcp,] <- global.prod[6,ggcm,gcm,rcp,]+prod
                weighted.prod[6,ggcm,gcm,rcp,] <- weighted.prod[6,ggcm,gcm,rcp,]+wprod
                global.prod.base[6,ggcm,gcm,rcp] <- global.prod.base[6,ggcm,gcm,rcp]+prod.base
                #} else{
                #cat("missing",fn,"\n")
              }
            }
          }
        }
        rownames(tabl2) <- ggcms
        colnames(tabl2) <- gcms5
        write.csv2(t(tabl2),paste0(path.figs,"cmip5_","GCM_GGCM_matches_",crops.param[cr],sets,".csv"))
      }
      # define delta by deviding by first time step to make figures look nicer
      delta <- global.prod
      for(i in 1:dim(global.prod)[5])
        #delta[,,,,i] <- global.prod[,,,,i]/global.prod.base
        delta[,,,,i] <- global.prod[,,,,i]/global.prod[,,,,1]
      deltaw <- weighted.prod
      for(i in 1:dim(weighted.prod)[5])
        deltaw[,,,,i] <- weighted.prod[,,,,i]/weighted.prod[,,,,1]
      
      delta[!is.finite(delta)] <- NA
      deltaw[!is.finite(deltaw)] <- NA
      
      # end of century mean
      delta_eoc <- apply(delta[,,,,60:74],c(1:4),mean,na.rm=T)
      delta_meangcm <- apply(delta_eoc,c(1,2,4),mean,na.rm=T)
      delta_meanggcm <- apply(delta_eoc,c(1,3,4),mean,na.rm=T)
      
      # time series
      deltats_meangcm <- apply(delta,c(1,2,4,5),mean,na.rm=T)
      deltats_meanggcm <- apply(delta,c(1,3,4,5),mean,na.rm=T)
      
      # end of century mean
      deltaw_eoc <- apply(deltaw[,,,,60:74],c(1:4),mean,na.rm=T)
      deltaw_meangcm <- apply(deltaw_eoc,c(1,2,4),mean,na.rm=T)
      deltaw_meanggcm <- apply(deltaw_eoc,c(1,3,4),mean,na.rm=T)
      
      # time series
      deltawts_meangcm <- apply(deltaw,c(1,2,4,5),mean,na.rm=T)
      deltawts_meanggcm <- apply(deltaw,c(1,3,4,5),mean,na.rm=T)
      
      #leave one-GCM out
      drop1gcm <- list()
      drop1gcm_meangcm <- list()
      drop1gcm_meanggcm <- list()
      dropw1gcm <- list()
      dropw1gcm_meangcm <- list()
      dropw1gcm_meanggcm <- list()
      for(i in 1:dim(delta)[3]){
        cat("processing",i,"\n")
        delta2 <- delta[,,-i,,]
        drop1gcm[[i]] <- delta2
        drop1gcm_meangcm[[i]] <- apply(delta2,c(1,2,4,5),mean,na.rm=T)
        drop1gcm_meanggcm[[i]] <- apply(delta2,c(1,3,4,5),mean,na.rm=T)
        delta2 <- deltaw[,,-i,,]
        dropw1gcm[[i]] <- delta2
        dropw1gcm_meangcm[[i]] <- apply(delta2,c(1,2,4,5),mean,na.rm=T)
        dropw1gcm_meanggcm[[i]] <- apply(delta2,c(1,3,4,5),mean,na.rm=T)
      }
      
      
      # for(i in 1:dim(index5)[2]){
      #   delta2 <- deltaw[,,-index5[,i],4,]
      #   dropw5gcm[[i]] <- delta2
      #   dropw5gcm_meangcm[[i]] <- apply(delta2,c(1,2,4),mean,na.rm=T)
      #   dropw5gcm_meanggcm[[i]] <- apply(delta2,c(1,3,4),mean,na.rm=T)
      # }
      
      
      #leave one-GGCM out
      drop1ggcm <- list()
      drop1ggcm_meangcm <- list()
      drop1ggcm_meanggcm <- list()
      dropw1ggcm <- list()
      dropw1ggcm_meangcm <- list()
      dropw1ggcm_meanggcm <- list()
      for(i in 1:dim(delta)[2]){
        delta2 <- delta[,-i,,,]
        drop1ggcm[[i]] <- delta2
        drop1ggcm_meangcm[[i]] <- apply(delta2,c(1,2,4,5),mean,na.rm=T)
        drop1ggcm_meanggcm[[i]] <- apply(delta2,c(1,3,4,5),mean,na.rm=T)
        delta2 <- deltaw[,-i,,,]
        dropw1ggcm[[i]] <- delta2
        dropw1ggcm_meangcm[[i]] <- apply(delta2,c(1,2,4,5),mean,na.rm=T)
        dropw1ggcm_meanggcm[[i]] <- apply(delta2,c(1,3,4,5),mean,na.rm=T)
      }
      
      for(cr in 1:length(crops.nice)){
        totalvar <- meangcmvar <- meanggcmvar <- array(NA,dim(delta)[c(4,5)])
        totalvarw <- meangcmvarw <- meanggcmvarw <- array(NA,dim(deltaw)[c(4,5)])
        for(rcp in 1:length(rcps)){
          for(i in 1:dim(totalvar)[2]){
            totalvar[rcp,i] <- var(as.vector(delta[cr,,,rcp,i]),na.rm=T)
            meangcmvar[rcp,i] <- var(as.vector(deltats_meangcm[cr,,rcp,i]),na.rm=T)
            meanggcmvar[rcp,i] <- var(as.vector(deltats_meanggcm[cr,,rcp,i]),na.rm=T)
            totalvarw[rcp,i] <- var(as.vector(deltaw[cr,,,rcp,i]),na.rm=T)
            meangcmvarw[rcp,i] <- var(as.vector(deltawts_meangcm[cr,,rcp,i]),na.rm=T)
            meanggcmvarw[rcp,i] <- var(as.vector(deltawts_meanggcm[cr,,rcp,i]),na.rm=T)
          }
        }
        
        # leave out 5GCMs for maize only
        if(cr==1 & sets==""){
          mi2 <- list(min=999,tvar=NULL,meangcmvar=NULL,meanggcmvar=NULL,ind=NULL)
          ma2 <- list(max=-999,tvar=NULL,meangcmvar=NULL,meanggcmvar=NULL,ind=NULL)
          for(bit in c(1:ceiling(dim(index5)[2]/sep))){
            #cat("processing",bit,"\n")
            fn <- paste0(path.ryield5,"drop5gcms_set_",bit,"-tenthousand.Rdata")
            load(fn)
            #drop5gcm,drop5gcm_meangcm,drop5gcm_meanggcm
            for(i in 1:length(mi)){
              if(mi$min<mi2$min){
                mi2$min<-mi$min
                mi2$tvar<- mi$tvar
                mi2$meangcmvar<-mi$meangcmvar
                mi2$meanggcmvar<-mi$meanggcmvar
                mi2$ind<-mi$ind
              } 
              if(ma$max>ma2$max){
                ma2$max<-ma$max
                ma2$tvar<- ma$tvar
                ma2$meangcmvar<-ma$meangcmvar
                ma2$meanggcmvar<-ma$meanggcmvar
                ma2$ind<-ma$ind
              }
            }
            
          }
          png(paste0(path.figs,"cmip5_","variance_shares_rcp85_",crops.param2[cr],sets,"_drop5gcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
          par(mar=c(4,4,2,2.5),xpd=NA)
          for(rcp in length(rcps)){
            plot(c(2011:2084),totalvar[rcp,],type="l",main=rcps.nice[rcp],xlab="year (AD)",ylab="variance",
                 axes=T,ylim=range(totalvar[rcp,])*1.2)
            #axis(1)
            #axis(2)
            lines(c(2011:2084),meangcmvar[rcp,],col=2)
            lines(c(2011:2084),meanggcmvar[rcp,],col=3)
            lines(c(2011:2084),meanggcmvar[rcp,]+meangcmvar[rcp,],col=4)
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
            lines(c(2011:2084),mi2$tvar,col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),mi2$meangcmvar,col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),mi2$meanggcmvar,col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),mi2$meanggcmvar+mi$meangcmvar,col=4,lty=2,lwd=0.5)
            lines(c(2011:2084),ma2$tvar,col=1,lty=3,lwd=0.5)
            lines(c(2011:2084),ma2$meangcmvar,col=2,lty=3,lwd=0.5)
            lines(c(2011:2084),ma2$meanggcmvar,col=3,lty=3,lwd=0.5)
            lines(c(2011:2084),ma2$meanggcmvar+ma$meangcmvar,col=4,lty=3,lwd=0.5)
            text(2084,mi2$meanggcmvar[74],paste(gcms5[index5[,mi2$ind]],collapse="\n"),adj=c(1,1),cex=0.5,col=1)
            text(2084,ma2$meanggcmvar[74],paste(gcms5[index5[,ma2$ind]],collapse="\n"),adj=c(1,-0.2),cex=0.5,col=1)
            
          }
          dev.off()
          
        }
        
        png(paste0(path.figs,"cmip5_","variance_shares_per_rcp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        split.screen(c(2,2))
        for(rcp in 1:length(rcps)){
          screen(rcp)
          par(mar=c(4,4,2,0.5))
          plot(c(2011:2084),totalvar[rcp,],type="l",main=rcps.nice[rcp],xlab="year (AD)",ylab="variance")
          lines(c(2011:2084),meangcmvar[rcp,],col=2)
          lines(c(2011:2084),meanggcmvar[rcp,],col=3)
          lines(c(2011:2084),meanggcmvar[rcp,]+meangcmvar[rcp,],col=4)
          if(rcp==length(rcps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
        }
        close.screen(all=T)
        dev.off()
        
        png(paste0(path.figs,"cmip5_","variance_shares_weightedprod_per_rcp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        split.screen(c(2,2))
        for(rcp in 1:length(rcps)){
          screen(rcp)
          par(mar=c(4,4,2,0.5))
          plot(c(2011:2084),totalvarw[rcp,],type="l",main=paste("production weighted",rcps.nice[rcp]),xlab="year (AD)",ylab="variance")
          lines(c(2011:2084),meangcmvarw[rcp,],col=2)
          lines(c(2011:2084),meanggcmvarw[rcp,],col=3)
          lines(c(2011:2084),meanggcmvarw[rcp,]+meangcmvarw[rcp,],col=4)
          if(rcp==length(rcps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
        }
        close.screen(all=T)
        dev.off()
        
        png(paste0(path.figs,"cmip5_","variance_shares_rcp85_",crops.param2[cr],sets,"_drop1ggcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        for(rcp in length(rcps)){
          plot(c(2011:2084),totalvar[rcp,],type="l",main=rcps.nice[rcp],xlab="year (AD)",ylab="variance",
               axes=F,ylim=range(totalvar[rcp,])*1.2)
          axis(1)
          axis(2)
          lines(c(2011:2084),meangcmvar[rcp,],col=2)
          lines(c(2011:2084),meanggcmvar[rcp,],col=3)
          lines(c(2011:2084),meanggcmvar[rcp,]+meangcmvar[rcp,],col=4)
          if(rcp==length(rcps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
          for(ii in 1:length(drop1ggcm))
          {
            delta_b <- drop1ggcm[[ii]]
            delta_b_meangcm <- drop1ggcm_meangcm[[ii]]
            delta_b_meanggcm <- drop1ggcm_meanggcm[[ii]]
            totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[rcp,i] <- var(as.vector(delta_b[cr,,,rcp,i]),na.rm=T)
              meangcmvar3[rcp,i] <- var(as.vector(delta_b_meangcm[cr,,rcp,i]),na.rm=T)
              meanggcmvar3[rcp,i] <- var(as.vector(delta_b_meanggcm[cr,,rcp,i]),na.rm=T)
            }
            lines(c(2011:2084),totalvar3[rcp,],col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),meangcmvar3[rcp,],col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[rcp,],col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[rcp,]+meangcmvar3[rcp,],col=4,lty=2,lwd=0.5)
            text(2084,totalvar3[rcp,74],ggcms[ii],adj=0,cex=0.5,col=1)
            text(2084,meangcmvar3[rcp,74],ggcms[ii],adj=0,cex=0.5,col=2)
            text(2084,meanggcmvar3[rcp,74],ggcms[ii],adj=0,cex=0.5,col=3)
            text(2084,meanggcmvar3[rcp,74]+meangcmvar3[rcp,74],ggcms[ii],adj=0,cex=0.5,col=4)
          }
        }
        dev.off()
        
        
        png(paste0(path.figs,"cmip5_","variance_shares_weightedprod_rcp85_",crops.param2[cr],sets,"_drop1ggcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        for(rcp in length(rcps)){
          plot(c(2011:2084),totalvarw[rcp,],type="l",main=paste("production weighted",rcps.nice[rcp]),xlab="year (AD)",ylab="variance",
               axes=F,ylim=range(totalvarw[rcp,],meangcmvarw[rcp,],meanggcmvarw[rcp,])*1.2)
          axis(1)
          axis(2)
          lines(c(2011:2084),meangcmvarw[rcp,],col=2)
          lines(c(2011:2084),meanggcmvarw[rcp,],col=3)
          lines(c(2011:2084),meanggcmvarw[rcp,]+meangcmvarw[rcp,],col=4)
          if(rcp==length(rcps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
          for(ii in 1:length(dropw1ggcm))
          {
            delta_b <- dropw1ggcm[[ii]]
            delta_b_meangcm <- dropw1ggcm_meangcm[[ii]]
            delta_b_meanggcm <- dropw1ggcm_meanggcm[[ii]]
            totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[rcp,i] <- var(as.vector(delta_b[cr,,,rcp,i]),na.rm=T)
              meangcmvar3[rcp,i] <- var(as.vector(delta_b_meangcm[cr,,rcp,i]),na.rm=T)
              meanggcmvar3[rcp,i] <- var(as.vector(delta_b_meanggcm[cr,,rcp,i]),na.rm=T)
            }
            lines(c(2011:2084),totalvar3[rcp,],col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),meangcmvar3[rcp,],col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[rcp,],col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[rcp,]+meangcmvar3[rcp,],col=4,lty=2,lwd=0.5)
            text(2084,totalvar3[rcp,74],ggcms[ii],adj=0,cex=0.5,col=1)
            text(2084,meangcmvar3[rcp,74],ggcms[ii],adj=0,cex=0.5,col=2)
            text(2084,meanggcmvar3[rcp,74],ggcms[ii],adj=0,cex=0.5,col=3)
            text(2084,meanggcmvar3[rcp,74]+meangcmvar3[rcp,74],ggcms[ii],adj=0,cex=0.5,col=4)
          }
        }
        dev.off()
        
        png(paste0(path.figs,"cmip5_","variance_shares_rcp85_",crops.param2[cr],sets,"_drop1gcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        for(rcp in length(rcps)){
          plot(c(2011:2084),totalvar[rcp,],type="l",main=rcps.nice[rcp],xlab="year (AD)",ylab="variance",
               axes=F,ylim=range(totalvar[rcp,])*1.2)
          axis(1)
          axis(2)
          lines(c(2011:2084),meangcmvar[rcp,],col=2)
          lines(c(2011:2084),meanggcmvar[rcp,],col=3)
          lines(c(2011:2084),meanggcmvar[rcp,]+meangcmvar[rcp,],col=4)
          if(rcp==length(rcps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
          for(ii in 1:length(drop1ggcm))
          {
            delta_b <- drop1gcm[[ii]]
            delta_b_meangcm <- drop1gcm_meangcm[[ii]]
            delta_b_meanggcm <- drop1gcm_meanggcm[[ii]]
            totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[rcp,i] <- var(as.vector(delta_b[cr,,,rcp,i]),na.rm=T)
              meangcmvar3[rcp,i] <- var(as.vector(delta_b_meangcm[cr,,rcp,i]),na.rm=T)
              meanggcmvar3[rcp,i] <- var(as.vector(delta_b_meanggcm[cr,,rcp,i]),na.rm=T)
            }
            lines(c(2011:2084),totalvar3[rcp,],col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),meangcmvar3[rcp,],col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[rcp,],col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[rcp,]+meangcmvar3[rcp,],col=4,lty=2,lwd=0.5)
            text(2084,totalvar3[rcp,74],gcms5[ii],adj=0,cex=0.5,col=1)
            text(2084,meangcmvar3[rcp,74],gcms5[ii],adj=0,cex=0.5,col=2)
            text(2084,meanggcmvar3[rcp,74],gcms5[ii],adj=0,cex=0.5,col=3)
            text(2084,meanggcmvar3[rcp,74]+meangcmvar3[rcp,74],gcms5[ii],adj=0,cex=0.5,col=4)
          }
        }
        dev.off()
        
        png(paste0(path.figs,"cmip5_","variance_shares_weightedprod_rcp85_",crops.param2[cr],sets,"_drop1gcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        for(rcp in length(rcps)){
          plot(c(2011:2084),totalvarw[rcp,],type="l",main=paste("production weighted",rcps.nice[rcp]),xlab="year (AD)",ylab="variance",
               axes=F,ylim=range(totalvarw[rcp,])*1.2)
          axis(1)
          axis(2)
          lines(c(2011:2084),meangcmvarw[rcp,],col=2)
          lines(c(2011:2084),meanggcmvarw[rcp,],col=3)
          lines(c(2011:2084),meanggcmvarw[rcp,]+meangcmvarw[rcp,],col=4)
          if(rcp==length(rcps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
          for(ii in 1:length(dropw1ggcm))
          {
            delta_b <- dropw1gcm[[ii]]
            delta_b_meangcm <- dropw1gcm_meangcm[[ii]]
            delta_b_meanggcm <- dropw1gcm_meanggcm[[ii]]
            totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[rcp,i] <- var(as.vector(delta_b[cr,,,rcp,i]),na.rm=T)
              meangcmvar3[rcp,i] <- var(as.vector(delta_b_meangcm[cr,,rcp,i]),na.rm=T)
              meanggcmvar3[rcp,i] <- var(as.vector(delta_b_meanggcm[cr,,rcp,i]),na.rm=T)
            }
            lines(c(2011:2084),totalvar3[rcp,],col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),meangcmvar3[rcp,],col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[rcp,],col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[rcp,]+meangcmvar3[rcp,],col=4,lty=2,lwd=0.5)
            text(2084,totalvar3[rcp,74],gcms5[ii],adj=0,cex=0.5,col=1)
            text(2084,meangcmvar3[rcp,74],gcms5[ii],adj=0,cex=0.5,col=2)
            text(2084,meanggcmvar3[rcp,74],gcms5[ii],adj=0,cex=0.5,col=3)
            text(2084,meanggcmvar3[rcp,74]+meangcmvar3[rcp,74],gcms5[ii],adj=0,cex=0.5,col=4)
          }
        }
        dev.off()
        
        if(length(rcps)==3){
          png(paste0(path.figs,"cmip5_","variance_relative_shares_per_rcp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=3*600,res=600,pointsize=6)
          split.screen(c(1,3))
        } else
        {
          png(paste0(path.figs,"cmip5_","variance_relative_shares_per_rcp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
          split.screen(c(2,2))
        }
        for(rcp in 1:length(rcps)){
          screen(rcp)
          par(mar=c(4,4,2,0.5))
          plot(c(2011:2084),totalvar[rcp,]/totalvar[rcp,],type="l",ylim=c(0,1.2),
               main=rcps.nice[rcp],xlab="year (AD)",ylab="relative contribution in overall variance")
          polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvar)[2]),rev(meangcmvar[rcp,]/totalvar[rcp,])),col=col.var[1],border=NA)
          lower <- upper2 <- upper <- (meanggcmvar[rcp,]+meangcmvar[rcp,])/totalvar[rcp,]
          upper[upper>1] <- 1
          polygon(c(2011:2084,2084:2011),c(meangcmvar[rcp,]/totalvar[rcp,],rev(upper)),col=col.var[2],border=NA)
          lower[upper2>1] <- 1
          upper2[upper2<1] <- 1
          polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
          
          #lines(c(2011:2084),meangcmvar[rcp,]/totalvar[rcp,],col=2)
          #lines(c(2011:2084),meanggcmvar[rcp,]/totalvar[rcp,],col=3)
          #lines(c(2011:2084),(meanggcmvar[rcp,]+meangcmvar[rcp,])/totalvar[rcp,],col=4)
          legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        }
        close.screen(all=T)
        dev.off()
        
        png(paste0(path.figs,"cmip5_","variance_relative_shares_weightedprod_per_rcp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        split.screen(c(2,2))
        for(rcp in 1:length(rcps)){
          screen(rcp)
          par(mar=c(4,4,2,0.5))
          plot(c(2011:2084),totalvarw[rcp,]/totalvarw[rcp,],type="l",ylim=c(0,1.2),
               main=paste("production weighted",rcps.nice[rcp]),xlab="year (AD)",ylab="relative contribution in overall variance")
          polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvarw)[2]),rev(meangcmvarw[rcp,]/totalvarw[rcp,])),col=col.var[1],border=NA)
          lower <- upper2 <- upper <- (meanggcmvarw[rcp,]+meangcmvarw[rcp,])/totalvarw[rcp,]
          upper[upper>1] <- 1
          polygon(c(2011:2084,2084:2011),c(meangcmvarw[rcp,]/totalvarw[rcp,],rev(upper)),col=col.var[2],border=NA)
          lower[upper2>1] <- 1
          upper2[upper2<1] <- 1
          polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
          
          #lines(c(2011:2084),meangcmvar[rcp,]/totalvar[rcp,],col=2)
          #lines(c(2011:2084),meanggcmvar[rcp,]/totalvar[rcp,],col=3)
          #lines(c(2011:2084),(meanggcmvar[rcp,]+meangcmvar[rcp,])/totalvar[rcp,],col=4)
          legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        }
        close.screen(all=T)
        dev.off()
        
        
        png(paste0(path.figs,"cmip5_","variance_relative_shares_rcp85_leave1gcmout_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        plot(c(2011:2084),totalvar[rcp,]/totalvar[rcp,],type="l",ylim=c(0,1.2),
             main=rcps.nice[rcp],xlab="year (AD)",ylab="relative contribution in overall variance")
        polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvar)[2]),rev(meangcmvar[rcp,]/totalvar[rcp,])),col=col.var[1],border=NA)
        lower <- upper2 <- upper <- (meanggcmvar[rcp,]+meangcmvar[rcp,])/totalvar[rcp,]
        upper[upper>1] <- 1
        polygon(c(2011:2084,2084:2011),c(meangcmvar[rcp,]/totalvar[rcp,],rev(upper)),col=col.var[2],border=NA)
        lower[upper2>1] <- 1
        upper2[upper2<1] <- 1
        polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
        
        for(ii in 1:length(drop1gcm))
        {
          delta_b <- drop1gcm[[ii]]
          delta_b_meangcm <- drop1gcm_meangcm[[ii]]
          delta_b_meanggcm <- drop1gcm_meanggcm[[ii]]
          totalvar2 <- meangcmvar2 <- meanggcmvar2 <- array(NA,dim(delta_b)[c(4,5)])
          for(rcp in length(rcps)){
            for(i in 1:dim(totalvar)[2]){
              totalvar2[rcp,i] <- var(as.vector(delta_b[cr,,,rcp,i]),na.rm=T)
              meangcmvar2[rcp,i] <- var(as.vector(delta_b_meangcm[cr,,rcp,i]),na.rm=T)
              meanggcmvar2[rcp,i] <- var(as.vector(delta_b_meanggcm[cr,,rcp,i]),na.rm=T)
            }
            lines(c(2011:2084),meangcmvar2[rcp,]/totalvar2[rcp,],col="darkgreen")
            lines(c(2011:2084),(meanggcmvar2[rcp,]+meangcmvar2[rcp,])/totalvar2[rcp,],col=4)
            text(2084,meangcmvar3[rcp,74]/totalvar3[rcp,74],gcms5[ii],adj=0,cex=0.5)
          }
        }
        legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        dev.off()
        
        png(paste0(path.figs,"cmip5_","variance_relative_shares_weightedprod_rcp85_leave1gcmout_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        plot(c(2011:2084),totalvarw[rcp,]/totalvarw[rcp,],type="l",ylim=c(0,1.2),
             main=paste("production weighted",rcps.nice[rcp]),xlab="year (AD)",ylab="relative contribution in overall variance")
        polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvarw)[2]),rev(meangcmvarw[rcp,]/totalvarw[rcp,])),col=col.var[1],border=NA)
        lower <- upper2 <- upper <- (meanggcmvarw[rcp,]+meangcmvarw[rcp,])/totalvarw[rcp,]
        upper[upper>1] <- 1
        polygon(c(2011:2084,2084:2011),c(meangcmvarw[rcp,]/totalvarw[rcp,],rev(upper)),col=col.var[2],border=NA)
        lower[upper2>1] <- 1
        upper2[upper2<1] <- 1
        polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
        
        for(ii in 1:length(dropw1gcm))
        {
          delta_b <- dropw1gcm[[ii]]
          delta_b_meangcm <- dropw1gcm_meangcm[[ii]]
          delta_b_meanggcm <- dropw1gcm_meanggcm[[ii]]
          totalvar2 <- meangcmvar2 <- meanggcmvar2 <- array(NA,dim(delta_b)[c(4,5)])
          for(rcp in length(rcps)){
            for(i in 1:dim(totalvar)[2]){
              totalvar2[rcp,i] <- var(as.vector(delta_b[cr,,,rcp,i]),na.rm=T)
              meangcmvar2[rcp,i] <- var(as.vector(delta_b_meangcm[cr,,rcp,i]),na.rm=T)
              meanggcmvar2[rcp,i] <- var(as.vector(delta_b_meanggcm[cr,,rcp,i]),na.rm=T)
            }
            lines(c(2011:2084),meangcmvar2[rcp,]/totalvar2[rcp,],col="darkgreen")
            lines(c(2011:2084),(meanggcmvar2[rcp,]+meangcmvar2[rcp,])/totalvar2[rcp,],col=4)
            text(2084,meangcmvar3[rcp,74]/totalvar3[rcp,74],gcms5[ii],adj=0,cex=0.5)
          }
        }
        legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        dev.off()
        
        
        png(paste0(path.figs,"cmip5_","variance_relative_shares_rcp85_leave1ggcmout_",crops.param2[cr],sets,"_v",version,".png"),
            width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        plot(c(2011:2084),totalvar[rcp,]/totalvar[rcp,],type="l",ylim=c(0,1.2),axes=F,
             main=rcps.nice[rcp],xlab="year (AD)",ylab="relative contribution in overall variance")
        axis(1)
        axis(2)
        polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvar)[2]),rev(meangcmvar[rcp,]/totalvar[rcp,])),col=col.var[1],border=NA)
        lower <- upper2 <- upper <- (meanggcmvar[rcp,]+meangcmvar[rcp,])/totalvar[rcp,]
        upper[upper>1] <- 1
        polygon(c(2011:2084,2084:2011),c(meangcmvar[rcp,]/totalvar[rcp,],rev(upper)),col=col.var[2],border=NA)
        lower[upper2>1] <- 1
        upper2[upper2<1] <- 1
        polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
        
        for(ii in 1:length(drop1ggcm))
        {
          delta_b <- drop1ggcm[[ii]]
          delta_b_meangcm <- drop1ggcm_meangcm[[ii]]
          delta_b_meanggcm <- drop1ggcm_meanggcm[[ii]]
          totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
          for(rcp in length(rcps)){
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[rcp,i] <- var(as.vector(delta_b[cr,,,rcp,i]),na.rm=T)
              meangcmvar3[rcp,i] <- var(as.vector(delta_b_meangcm[cr,,rcp,i]),na.rm=T)
              meanggcmvar3[rcp,i] <- var(as.vector(delta_b_meanggcm[cr,,rcp,i]),na.rm=T)
            }
            lines(c(2011:2084),meangcmvar3[rcp,]/totalvar3[rcp,],col="darkgreen")
            lines(c(2011:2084),(meanggcmvar3[rcp,]+meangcmvar3[rcp,])/totalvar3[rcp,],col=4)
            text(2084,meangcmvar3[rcp,74]/totalvar3[rcp,74],ggcms[ii],adj=0,cex=0.5)
            text(2084,(meanggcmvar3[rcp,74]+meangcmvar3[rcp,74])/totalvar3[rcp,74],ggcms[ii],adj=0,cex=0.5)
          }
        }
        legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        
        dev.off()
        
        png(paste0(path.figs,"cmip5_","variance_relative_shares_weightedprod_rcp85_leave1ggcmout_",crops.param2[cr],sets,"_v",version,".png"),
            width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        plot(c(2011:2084),totalvarw[rcp,]/totalvarw[rcp,],type="l",ylim=c(0,1.2),axes=F,
             main=paste("production weighted",rcps.nice[rcp]),xlab="year (AD)",
             ylab="relative contribution in overall variance")
        axis(1)
        axis(2)
        polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvarw)[2]),rev(meangcmvarw[rcp,]/totalvarw[rcp,])),col=col.var[1],border=NA)
        lower <- upper2 <- upper <- (meanggcmvarw[rcp,]+meangcmvarw[rcp,])/totalvarw[rcp,]
        upper[upper>1] <- 1
        polygon(c(2011:2084,2084:2011),c(meangcmvarw[rcp,]/totalvarw[rcp,],rev(upper)),col=col.var[2],border=NA)
        lower[upper2>1] <- 1
        upper2[upper2<1] <- 1
        polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
        
        for(ii in 1:length(dropw1ggcm))
        {
          delta_b <- dropw1ggcm[[ii]]
          delta_b_meangcm <- dropw1ggcm_meangcm[[ii]]
          delta_b_meanggcm <- dropw1ggcm_meanggcm[[ii]]
          totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
          for(rcp in length(rcps)){
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[rcp,i] <- var(as.vector(delta_b[cr,,,rcp,i]),na.rm=T)
              meangcmvar3[rcp,i] <- var(as.vector(delta_b_meangcm[cr,,rcp,i]),na.rm=T)
              meanggcmvar3[rcp,i] <- var(as.vector(delta_b_meanggcm[cr,,rcp,i]),na.rm=T)
            }
            lines(c(2011:2084),meangcmvar3[rcp,]/totalvar3[rcp,],col="darkgreen")
            lines(c(2011:2084),(meanggcmvar3[rcp,]+meangcmvar3[rcp,])/totalvar3[rcp,],col=4)
            text(2084,meangcmvar3[rcp,74]/totalvar3[rcp,74],ggcms[ii],adj=0,cex=0.5)
            text(2084,(meanggcmvar3[rcp,74]+meangcmvar3[rcp,74])/totalvar3[rcp,74],ggcms[ii],adj=0,cex=0.5)
          }
        }
        legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        
        dev.off()
        
        png(paste0(path.figs,"cmip5_","global_delta",sets,"_ts_per_rcp_",crops.param2[cr],"_v",version,".png"),width=8*600,height=5*600,res=600,pointsize=12)
        plot(c(2011:2084),delta[cr,1,1,1,],ylim=range(delta[cr,,,,],na.rm=T),type="n",main=crops.nice[cr],
             xlab="year (AD)",ylab="change in global productivity [-]")
        for(rcp in 1:length(rcps)){
          med <- apply(delta[cr,,,rcp,],3,median,na.rm=T)
          mi <- apply(delta[cr,,,rcp,],3,min,na.rm=T)
          ma <- apply(delta[cr,,,rcp,],3,max,na.rm=T)
          std <- apply(delta[cr,,,rcp,],3,sd,na.rm=T)
          polygon(c(2011:2084,2084:2011),c(med-std,rev(med+std)),col=col.rcp2[rcp],border=NA)
          lines(c(2011:2084),med,col=col.rcp[rcp],lwd=2)
          lines(c(2011:2084),mi,col=col.rcp[rcp],lwd=1,lty=2)
          lines(c(2011:2084),ma,col=col.rcp[rcp],lwd=1,lty=2)
          # adding lines for RCP uncertainty ranges
          lines(rep(2084+rcp*0.6,2),c((med-std)[length(med)],(med+std)[length(med)]),col=col.rcp[rcp])
          # for(gcm in 1:length(gcms5)){
          #   for(ggcm in 1:length(ggcms)){
          #     lines(delta[1,ggcm,gcm,rcp,],col=ggcm,lty=rcp)
          #   }
          # }
          abline(h=1,lty=2)
        }
        legend("topleft",legend=c(rcps.nice,"median","min/max","+/- 1SD",""),#x.intersp = c(1,1,1,1,1,1,-0.5,1)+1,
               col=c(col.rcp,1,1,NA,NA),lty=c(rep(1,5),2,NA,NA),lwd=c(rep(1,4),2,1,NA,NA),bty="n",
               fill=c(NA,NA,NA,NA,NA,NA,rgb(0,0,0,45/255),NA),border=NA,ncol=2)
        dev.off()
        
        png(paste0(path.figs,"cmip5_","global_delta",sets,"_ts_per_rcp_",crops.param2[cr],"_weightedprod_v",version,".png"),width=8*600,height=5*600,res=600,pointsize=12)
        plot(c(2011:2084),deltaw[cr,1,1,1,],ylim=range(deltaw[cr,,,,],na.rm=T),type="n",main=paste("production weighted",crops.nice[cr]),
             xlab="year (AD)",ylab="change in global productivity [-]")
        for(rcp in 1:length(rcps)){
          med <- apply(deltaw[cr,,,rcp,],3,median,na.rm=T)
          mi <- apply(deltaw[cr,,,rcp,],3,min,na.rm=T)
          ma <- apply(deltaw[cr,,,rcp,],3,max,na.rm=T)
          std <- apply(deltaw[cr,,,rcp,],3,sd,na.rm=T)
          polygon(c(2011:2084,2084:2011),c(med-std,rev(med+std)),col=col.rcp2[rcp],border=NA)
          lines(c(2011:2084),med,col=col.rcp[rcp],lwd=2)
          lines(c(2011:2084),mi,col=col.rcp[rcp],lwd=1,lty=2)
          lines(c(2011:2084),ma,col=col.rcp[rcp],lwd=1,lty=2)
          # adding lines for RCP uncertainty ranges
          lines(rep(2084+rcp*0.6,2),c((med-std)[length(med)],(med+std)[length(med)]),col=col.rcp[rcp])
          # for(gcm in 1:length(gcms5)){
          #   for(ggcm in 1:length(ggcms)){
          #     lines(delta[1,ggcm,gcm,rcp,],col=ggcm,lty=rcp)
          #   }
          # }
          abline(h=1,lty=2)
        }
        legend("topleft",legend=c(rcps.nice,"median","min/max","+/- 1SD",""),#x.intersp = c(1,1,1,1,1,1,-0.5,1)+1,
               col=c(col.rcp,1,1,NA,NA),lty=c(rep(1,5),2,NA,NA),lwd=c(rep(1,4),2,1,NA,NA),bty="n",
               fill=c(NA,NA,NA,NA,NA,NA,rgb(0,0,0,45/255),NA),border=NA,ncol=2)
        dev.off()
        
        for(rcp in 1:length(rcps)){
          png(paste0(path.figs,"cmip5_","global_delta",sets,"_ts_",rcps.nice[rcp],"_",crops.param2[cr],"2_v",version,".png"),width=8*600,height=5*600,res=600,pointsize=12)
          plot(c(2011:2084),delta[cr,1,1,1,],ylim=range(delta[cr,,,,],na.rm=T),type="n",main=crops.nice[cr],
               xlab="year (AD)",ylab="change in global productivity [-]")
          par(mar=c(4,4,3,4),xpd=NA)
          for(ggcm in 1:length(ggcms)){
            if(all(!is.finite(delta[cr,ggcm,,rcp,])))
              next
            med <- apply(delta[cr,ggcm,,rcp,],2,median,na.rm=T)
            mi <- apply(delta[cr,ggcm,,rcp,],2,min,na.rm=T)
            ma <- apply(delta[cr,ggcm,,rcp,],2,max,na.rm=T)
            std <- apply(delta[cr,ggcm,,rcp,],2,sd,na.rm=T)
            polygon(c(2011:2084,2084:2011),c(med-std,rev(med+std)),col=col.ggcms2[ggcm],border=NA)
            lines(c(2011:2084),med,col=col.ggcms[ggcm],lwd=2)
            lines(c(2011:2084),mi,col=col.ggcms[ggcm],lwd=1,lty=2)
            lines(c(2011:2084),ma,col=col.ggcms[ggcm],lwd=1,lty=2)
            # adding lines for RCP uncertainty ranges
            lines(rep(2084+ggcm*0.6,2),c((med-std)[length(med)],(med+std)[length(med)]),col=col.ggcms[ggcm])
            # for(gcm in 1:length(gcms5)){
            #   for(ggcm in 1:length(ggcms)){
            #     lines(delta[1,ggcm,gcm,rcp,],col=ggcm,lty=rcp)
            #   }
            # }
          }
          par(xpd=FALSE)
          abline(h=1,lty=2)
          legend("bottomleft",legend=c(ggcms,"median","min/max","+/- 1SD",""),#x.intersp = c(1,1,1,1,1,1,-0.5,1)+1,
                 col=c(col.ggcms[1:length(ggcms)],1,1,NA,NA),lty=c(rep(1,length(ggcms)+1),2,NA,NA),lwd=c(rep(1,9),2,1,NA,NA),bty="n",
                 fill=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,rgb(0,0,0,45/255),NA),border=NA,ncol=3)
          dev.off()
          
          png(paste0(path.figs,"cmip5_","global_delta",sets,"_ts_",rcps.nice[rcp],"_",crops.param2[cr],"_weightedprod_v",version,".png"),width=8*600,height=5*600,res=600,pointsize=12)
          plot(c(2011:2084),deltaw[cr,1,1,1,],ylim=range(deltaw[cr,,,,],na.rm=T),type="n",main=paste("production weighted",crops.nice[cr]),
               xlab="year (AD)",ylab="change in global productivity [-]")
          par(mar=c(4,4,3,4),xpd=NA)
          for(ggcm in 1:length(ggcms)){
            if(all(!is.finite(deltaw[cr,ggcm,,rcp,])))
              next
            med <- apply(deltaw[cr,ggcm,,rcp,],2,median,na.rm=T)
            mi <- apply(deltaw[cr,ggcm,,rcp,],2,min,na.rm=T)
            ma <- apply(deltaw[cr,ggcm,,rcp,],2,max,na.rm=T)
            std <- apply(deltaw[cr,ggcm,,rcp,],2,sd,na.rm=T)
            polygon(c(2011:2084,2084:2011),c(med-std,rev(med+std)),col=col.ggcms2[ggcm],border=NA)
            lines(c(2011:2084),med,col=col.ggcms[ggcm],lwd=2)
            lines(c(2011:2084),mi,col=col.ggcms[ggcm],lwd=1,lty=2)
            lines(c(2011:2084),ma,col=col.ggcms[ggcm],lwd=1,lty=2)
            # adding lines for RCP uncertainty ranges
            lines(rep(2084+ggcm*0.6,2),c((med-std)[length(med)],(med+std)[length(med)]),col=col.ggcms[ggcm])
            # for(gcm in 1:length(gcms5)){
            #   for(ggcm in 1:length(ggcms)){
            #     lines(delta[1,ggcm,gcm,rcp,],col=ggcm,lty=rcp)
            #   }
            # }
          }
          par(xpd=FALSE)
          abline(h=1,lty=2)
          legend("topleft",legend=c(ggcms,"median","min/max","+/- 1SD",""),#x.intersp = c(1,1,1,1,1,1,-0.5,1)+1,
                 col=c(col.ggcms[1:length(ggcms)],1,1,NA,NA),lty=c(rep(1,12),2,NA,NA),lwd=c(rep(1,9),2,1,NA,NA),bty="n",
                 fill=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,rgb(0,0,0,45/255),NA),border=NA,ncol=2)
          dev.off()
          
        }#rcp
      }#cr
    }
    if(do.cmip6){
      global.prod6 <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      weighted.prod6 <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps),length(2011:2084)))
      global.prod.base6 <- array(0,dim=c(length(crops.param)+1,length(ggcms),length(gcms6),length(ssps)))
      for(cr in 1:length(crops.param)){
        tabl2 <- array("",dim=c(length(ggcms),length(gcms6)))
        for(ssp in 1:length(ssps)){
          for(ggcm in 1:length(ggcms)){
            for(gcm in 1:length(gcms6)){
              if(sets=="")
                fn <- paste0(path.rprod6,"A0/","cmip6_global_production_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow_v",version,".Rdata")
              else if(sets=="_constco2")
                fn <- paste0(path.rprod6,"A0/","cmip6_global_production_constc_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_movingwindow_v",version,".Rdata")
              else if(sets=="_A1")
                fn <- paste0(path.rprod6,"A1/","cmip6_global_production_",crops.param[cr],"_",ssps[ssp],"_",ggcms[ggcm],"_",gcms6[gcm],"_A1_movingwindow_v",version,".Rdata")
              if(file.exists(fn)){
                tabl2[ggcm,gcm] <- paste(tabl2[ggcm,gcm],ssps.nr[ssp],sep= if(""==tabl2[ggcm,gcm]) "" else ", ")
                load(fn)
                global.prod6[cr,ggcm,gcm,ssp,] <- prod
                weighted.prod6[cr,ggcm,gcm,ssp,] <- wprod
                global.prod.base6[cr,ggcm,gcm,ssp] <- prod.base
                global.prod6[6,ggcm,gcm,ssp,] <- global.prod6[6,ggcm,gcm,ssp,]+prod
                weighted.prod6[6,ggcm,gcm,ssp,] <- weighted.prod6[6,ggcm,gcm,ssp,]+wprod
                global.prod.base6[6,ggcm,gcm,ssp] <- global.prod.base6[6,ggcm,gcm,ssp]+prod.base
              } else{
                if(ggcm==1 & gcm==1) cat("missing",fn,"\n")
              }
            }
          }
        }
        rownames(tabl2) <- ggcms
        colnames(tabl2) <- gcms6
        write.csv2(t(tabl2),paste0(path.figs,"cmip6_","GCM_GGCM_matches_",crops.param[cr],sets,".csv"))
      }
      # define delta6 by deviding by first time step to make figures look nicer
      delta6 <- global.prod6
      for(i in 1:dim(global.prod6)[5])
        #delta6[,,,,i] <- global.prod6[,,,,i]/global.prod.base6
        delta6[,,,,i] <- global.prod6[,,,,i]/global.prod6[,,,,1]
      deltaw6 <- weighted.prod6
      for(i in 1:dim(weighted.prod6)[5])
        deltaw6[,,,,i] <- weighted.prod6[,,,,i]/weighted.prod6[,,,,1]
      
      delta6[!is.finite(delta6)] <- NA
      deltaw6[!is.finite(deltaw6)] <- NA
      
      # end of century mean
      deltaw_eoc6 <- apply(deltaw6[,,,,60:74],c(1:4),mean,na.rm=T)
      deltaw_meangcm6 <- apply(deltaw_eoc6,c(1,2,4),mean,na.rm=T)
      deltaw_meanggcm6 <- apply(deltaw_eoc6,c(1,3,4),mean,na.rm=T)
      
      # time series
      deltats_meangcm6 <- apply(delta6,c(1,2,4,5),mean,na.rm=T)
      deltats_meanggcm6 <- apply(delta6,c(1,3,4,5),mean,na.rm=T)
      
      # end of century mean
      delta_eoc6 <- apply(delta6[,,,,60:74],c(1:4),mean,na.rm=T)
      delta_meangcm6 <- apply(delta_eoc6,c(1,2,4),mean,na.rm=T)
      delta_meanggcm6 <- apply(delta_eoc6,c(1,3,4),mean,na.rm=T)
      
      # time series
      deltawts_meangcm6 <- apply(deltaw6,c(1,2,4,5),mean,na.rm=T)
      deltawts_meanggcm6 <- apply(deltaw6,c(1,3,4,5),mean,na.rm=T)
      
      #leave one-GCM out
      drop1gcm6 <- list()
      drop1gcm_meangcm6 <- list()
      drop1gcm_meanggcm6 <- list()
      dropw1gcm6 <- list()
      dropw1gcm_meangcm6 <- list()
      dropw1gcm_meanggcm6 <- list()
      for(i in 1:dim(delta6)[3]){
        cat("processing",i,"\n")
        delta2 <- delta6[,,-i,,]
        drop1gcm6[[i]] <- delta2
        drop1gcm_meangcm6[[i]] <- apply(delta2,c(1,2,4,5),mean,na.rm=T)
        drop1gcm_meanggcm6[[i]] <- apply(delta2,c(1,3,4,5),mean,na.rm=T)
        delta2 <- deltaw6[,,-i,,]
        dropw1gcm6[[i]] <- delta2
        dropw1gcm_meangcm6[[i]] <- apply(delta2,c(1,2,4,5),mean,na.rm=T)
        dropw1gcm_meanggcm6[[i]] <- apply(delta2,c(1,3,4,5),mean,na.rm=T)
      }
      
      
      #leave one-GGCM out
      drop1ggcm <- list()
      drop1ggcm_meangcm <- list()
      drop1ggcm_meanggcm <- list()
      dropw1ggcm <- list()
      dropw1ggcm_meangcm <- list()
      dropw1ggcm_meanggcm <- list()
      for(i in 1:dim(delta6)[2]){
        delta2 <- delta6[,-i,,,]
        drop1ggcm[[i]] <- delta2
        drop1ggcm_meangcm[[i]] <- apply(delta2,c(1,2,4,5),mean,na.rm=T)
        drop1ggcm_meanggcm[[i]] <- apply(delta2,c(1,3,4,5),mean,na.rm=T)
        delta2 <- deltaw6[,-i,,,]
        dropw1ggcm[[i]] <- delta2
        dropw1ggcm_meangcm[[i]] <- apply(delta2,c(1,2,4,5),mean,na.rm=T)
        dropw1ggcm_meanggcm[[i]] <- apply(delta2,c(1,3,4,5),mean,na.rm=T)
      }
      
      for(cr in 1:length(crops.nice)){
        totalvar6 <- meangcmvar6 <- meanggcmvar6 <- array(NA,dim(delta6)[c(4,5)])
        totalvarw6 <- meangcmvarw6 <- meanggcmvarw6 <- array(NA,dim(deltaw6)[c(4,5)])
        for(ssp in 1:length(ssps)){
          for(i in 1:dim(totalvar6)[2]){
            totalvar6[ssp,i] <- var(as.vector(delta6[cr,,,ssp,i]),na.rm=T)
            meangcmvar6[ssp,i] <- var(as.vector(deltats_meangcm6[cr,,ssp,i]),na.rm=T)
            meanggcmvar6[ssp,i] <- var(as.vector(deltats_meanggcm6[cr,,ssp,i]),na.rm=T)
            totalvarw6[ssp,i] <- var(as.vector(deltaw6[cr,,,ssp,i]),na.rm=T)
            meangcmvarw6[ssp,i] <- var(as.vector(deltawts_meangcm6[cr,,ssp,i]),na.rm=T)
            meanggcmvarw6[ssp,i] <- var(as.vector(deltawts_meanggcm6[cr,,ssp,i]),na.rm=T)
          }
        }
        
        # leave out 5GCMs for maize only
        if(cr==1 & sets==""){
          mi6 <- list(min=999,tvar=NULL,meangcmvar6=NULL,meanggcmvar6=NULL,ind=NULL)
          ma6 <- list(max=-999,tvar=NULL,meangcmvar6=NULL,meanggcmvar6=NULL,ind=NULL)
          #for(bit in c(1:ceiling(dim(index6)[2]/sep))){
          #cat("processing",bit,"\n")
          fn <- paste0(path.ryield6,"drop3gcms_set.Rdata")
          load(fn)
          #drop5gcm,drop5gcm_meangcm,drop5gcm_meanggcm
          for(i in 1:length(mi)){
            if(mi$min<mi6$min){
              mi6$min<-mi$min
              mi6$tvar<- mi$tvar
              mi6$meangcmvar6<-mi$meangcmvar
              mi6$meanggcmvar6<-mi$meanggcmvar
              mi6$ind<-mi$ind
            } 
            if(ma$max>ma6$max){
              ma6$max<-ma$max
              ma6$tvar<- ma$tvar
              ma6$meangcmvar6<-ma$meangcmvar
              ma6$meanggcmvar6<-ma$meanggcmvar
              ma6$ind<-ma$ind
            }
          }
          
          #}
          png(paste0(path.figs,"cmip6_","variance_shares_ssp85_",crops.param2[cr],sets,"_drop5gcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
          par(mar=c(4,4,2,2.5),xpd=NA)
          for(ssp in length(ssps)){
            plot(c(2011:2084),totalvar6[ssp,],type="l",main=ssps.nice[ssp],xlab="year (AD)",ylab="variance",
                 axes=T,ylim=range(totalvar6[ssp,])*1.2)
            #axis(1)
            #axis(2)
            lines(c(2011:2084),meangcmvar6[ssp,],col=2)
            lines(c(2011:2084),meanggcmvar6[ssp,],col=3)
            lines(c(2011:2084),meanggcmvar6[ssp,]+meangcmvar6[ssp,],col=4)
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
            lines(c(2011:2084),mi6$tvar,col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),mi6$meangcmvar6,col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),mi6$meanggcmvar6,col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),mi6$meanggcmvar6+mi6$meangcmvar6,col=4,lty=2,lwd=0.5)
            lines(c(2011:2084),ma6$tvar,col=1,lty=3,lwd=0.5)
            lines(c(2011:2084),ma6$meangcmvar6,col=2,lty=3,lwd=0.5)
            lines(c(2011:2084),ma6$meanggcmvar6,col=3,lty=3,lwd=0.5)
            lines(c(2011:2084),ma6$meanggcmvar6+ma6$meangcmvar6,col=4,lty=3,lwd=0.5)
            text(2084,mi6$meanggcmvar6[74],paste(gcms6[index6[,mi6$ind]],collapse="\n"),adj=c(1,1),cex=0.5,col=1)
            text(2084,ma6$meanggcmvar6[74],paste(gcms6[index6[,ma6$ind]],collapse="\n"),adj=c(1,-0.2),cex=0.5,col=1)
            
          }
          dev.off()
          
        }
        
        png(paste0(path.figs,"cmip6_","variance_shares_per_ssp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        split.screen(c(2,2))
        for(ssp in 1:length(ssps)){
          screen(ssp)
          par(mar=c(4,4,2,0.5))
          plot(c(2011:2084),totalvar6[ssp,],type="l",main=ssps.nice[ssp],xlab="year (AD)",ylab="variance")
          lines(c(2011:2084),meangcmvar6[ssp,],col=2)
          lines(c(2011:2084),meanggcmvar6[ssp,],col=3)
          lines(c(2011:2084),meanggcmvar6[ssp,]+meangcmvar6[ssp,],col=4)
          if(ssp==4)
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
        }
        close.screen(all=T)
        dev.off()
        
        png(paste0(path.figs,"cmip6_","variance_shares_weightedprod_per_ssp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        split.screen(c(2,2))
        for(ssp in 1:length(ssps)){
          screen(ssp)
          par(mar=c(4,4,2,0.5))
          plot(c(2011:2084),totalvarw6[ssp,],type="l",main=paste("production weighted",ssps.nice[ssp]),xlab="year (AD)",ylab="variance")
          lines(c(2011:2084),meangcmvarw6[ssp,],col=2)
          lines(c(2011:2084),meanggcmvarw6[ssp,],col=3)
          lines(c(2011:2084),meanggcmvarw6[ssp,]+meangcmvarw6[ssp,],col=4)
          if(ssp==4)
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
        }
        close.screen(all=T)
        dev.off()
        
        png(paste0(path.figs,"cmip6_","variance_shares_ssp85_",crops.param2[cr],sets,"_drop1ggcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        for(ssp in length(ssps)){
          plot(c(2011:2084),totalvar6[ssp,],type="l",main=ssps.nice[ssp],xlab="year (AD)",ylab="variance",
               axes=F,ylim=range(totalvar6[ssp,])*1.2)
          axis(1)
          axis(2)
          lines(c(2011:2084),meangcmvar6[ssp,],col=2)
          lines(c(2011:2084),meanggcmvar6[ssp,],col=3)
          lines(c(2011:2084),meanggcmvar6[ssp,]+meangcmvar6[ssp,],col=4)
          if(ssp==4)
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
          for(ii in 1:length(drop1ggcm))
          {
            delta_b <- drop1ggcm[[ii]]
            delta_b_meangcm <- drop1ggcm_meangcm[[ii]]
            delta_b_meanggcm <- drop1ggcm_meanggcm[[ii]]
            totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[ssp,i] <- var(as.vector(delta_b[cr,,,ssp,i]),na.rm=T)
              meangcmvar3[ssp,i] <- var(as.vector(delta_b_meangcm[cr,,ssp,i]),na.rm=T)
              meanggcmvar3[ssp,i] <- var(as.vector(delta_b_meanggcm[cr,,ssp,i]),na.rm=T)
            }
            lines(c(2011:2084),totalvar3[ssp,],col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),meangcmvar3[ssp,],col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[ssp,],col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[ssp,]+meangcmvar3[ssp,],col=4,lty=2,lwd=0.5)
            text(2084,totalvar3[ssp,74],ggcms[ii],adj=0,cex=0.5,col=1)
            text(2084,meangcmvar3[ssp,74],ggcms[ii],adj=0,cex=0.5,col=2)
            text(2084,meanggcmvar3[ssp,74],ggcms[ii],adj=0,cex=0.5,col=3)
            text(2084,meanggcmvar3[ssp,74]+meangcmvar3[ssp,74],ggcms[ii],adj=0,cex=0.5,col=4)
          }
        }
        dev.off()
        
        
        png(paste0(path.figs,"cmip6_","variance_shares_weightedprod_ssp85_",crops.param2[cr],sets,"_drop1ggcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        for(ssp in length(ssps)){
          plot(c(2011:2084),totalvarw6[ssp,],type="l",main=paste("production weighted",ssps.nice[ssp]),xlab="year (AD)",ylab="variance",
               axes=F,ylim=range(totalvarw6[ssp,],meangcmvarw6[ssp,],meanggcmvarw6[ssp,])*1.2)
          axis(1)
          axis(2)
          lines(c(2011:2084),meangcmvarw6[ssp,],col=2)
          lines(c(2011:2084),meanggcmvarw6[ssp,],col=3)
          lines(c(2011:2084),meanggcmvarw6[ssp,]+meangcmvarw6[ssp,],col=4)
          if(ssp==length(ssps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
          for(ii in 1:length(dropw1ggcm))
          {
            delta_b <- dropw1ggcm[[ii]]
            delta_b_meangcm <- dropw1ggcm_meangcm[[ii]]
            delta_b_meanggcm <- dropw1ggcm_meanggcm[[ii]]
            totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[ssp,i] <- var(as.vector(delta_b[cr,,,ssp,i]),na.rm=T)
              meangcmvar3[ssp,i] <- var(as.vector(delta_b_meangcm[cr,,ssp,i]),na.rm=T)
              meanggcmvar3[ssp,i] <- var(as.vector(delta_b_meanggcm[cr,,ssp,i]),na.rm=T)
            }
            lines(c(2011:2084),totalvar3[ssp,],col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),meangcmvar3[ssp,],col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[ssp,],col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[ssp,]+meangcmvar3[ssp,],col=4,lty=2,lwd=0.5)
            text(2084,totalvar3[ssp,74],ggcms[ii],adj=0,cex=0.5,col=1)
            text(2084,meangcmvar3[ssp,74],ggcms[ii],adj=0,cex=0.5,col=2)
            text(2084,meanggcmvar3[ssp,74],ggcms[ii],adj=0,cex=0.5,col=3)
            text(2084,meanggcmvar3[ssp,74]+meangcmvar3[ssp,74],ggcms[ii],adj=0,cex=0.5,col=4)
          }
        }
        dev.off()
        
        png(paste0(path.figs,"cmip6_","variance_shares_ssp85_",crops.param2[cr],sets,"_drop1gcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        for(ssp in length(ssps)){
          plot(c(2011:2084),totalvar6[ssp,],type="l",main=ssps.nice[ssp],xlab="year (AD)",ylab="variance",
               axes=F,ylim=range(totalvar6[ssp,])*1.2)
          axis(1)
          axis(2)
          lines(c(2011:2084),meangcmvar6[ssp,],col=2)
          lines(c(2011:2084),meanggcmvar6[ssp,],col=3)
          lines(c(2011:2084),meanggcmvar6[ssp,]+meangcmvar6[ssp,],col=4)
          if(ssp==length(ssps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
          for(ii in 1:length(drop1ggcm))
          {
            delta_b <- drop1gcm6[[ii]]
            delta_b_meangcm <- drop1gcm_meangcm6[[ii]]
            delta_b_meanggcm <- drop1gcm_meanggcm6[[ii]]
            totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[ssp,i] <- var(as.vector(delta_b[cr,,,ssp,i]),na.rm=T)
              meangcmvar3[ssp,i] <- var(as.vector(delta_b_meangcm[cr,,ssp,i]),na.rm=T)
              meanggcmvar3[ssp,i] <- var(as.vector(delta_b_meanggcm[cr,,ssp,i]),na.rm=T)
            }
            lines(c(2011:2084),totalvar3[ssp,],col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),meangcmvar3[ssp,],col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[ssp,],col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[ssp,]+meangcmvar3[ssp,],col=4,lty=2,lwd=0.5)
            text(2084,totalvar3[ssp,74],gcms6[ii],adj=0,cex=0.5,col=1)
            text(2084,meangcmvar3[ssp,74],gcms6[ii],adj=0,cex=0.5,col=2)
            text(2084,meanggcmvar3[ssp,74],gcms6[ii],adj=0,cex=0.5,col=3)
            text(2084,meanggcmvar3[ssp,74]+meangcmvar3[ssp,74],gcms6[ii],adj=0,cex=0.5,col=4)
          }
        }
        dev.off()
        
        png(paste0(path.figs,"cmip6_","variance_shares_weightedprod_ssp85_",crops.param2[cr],sets,"_drop1gcm_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        for(ssp in length(ssps)){
          plot(c(2011:2084),totalvarw6[ssp,],type="l",main=paste("production weighted",ssps.nice[ssp]),xlab="year (AD)",ylab="variance",
               axes=F,ylim=range(totalvarw6[ssp,])*1.2)
          axis(1)
          axis(2)
          lines(c(2011:2084),meangcmvarw6[ssp,],col=2)
          lines(c(2011:2084),meanggcmvarw6[ssp,],col=3)
          lines(c(2011:2084),meanggcmvarw6[ssp,]+meangcmvarw6[ssp,],col=4)
          if(ssp==length(ssps))
            legend("topleft",legend=c("total variance","GGCM share","GCM share","GCM+GGCM shares"),col=c(1:4),lty=1,bty="n")
          for(ii in 1:length(dropw1ggcm))
          {
            delta_b <- dropw1gcm6[[ii]]
            delta_b_meangcm <- dropw1gcm_meangcm6[[ii]]
            delta_b_meanggcm <- dropw1gcm_meanggcm6[[ii]]
            totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[ssp,i] <- var(as.vector(delta_b[cr,,,ssp,i]),na.rm=T)
              meangcmvar3[ssp,i] <- var(as.vector(delta_b_meangcm[cr,,ssp,i]),na.rm=T)
              meanggcmvar3[ssp,i] <- var(as.vector(delta_b_meanggcm[cr,,ssp,i]),na.rm=T)
            }
            lines(c(2011:2084),totalvar3[ssp,],col=1,lty=2,lwd=0.5)
            lines(c(2011:2084),meangcmvar3[ssp,],col=2,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[ssp,],col=3,lty=2,lwd=0.5)
            lines(c(2011:2084),meanggcmvar3[ssp,]+meangcmvar3[ssp,],col=4,lty=2,lwd=0.5)
            text(2084,totalvar3[ssp,74],gcms6[ii],adj=0,cex=0.5,col=1)
            text(2084,meangcmvar3[ssp,74],gcms6[ii],adj=0,cex=0.5,col=2)
            text(2084,meanggcmvar3[ssp,74],gcms6[ii],adj=0,cex=0.5,col=3)
            text(2084,meanggcmvar3[ssp,74]+meangcmvar3[ssp,74],gcms6[ii],adj=0,cex=0.5,col=4)
          }
        }
        dev.off()
        
        if(length(ssps)==3){
          png(paste0(path.figs,"cmip6_","variance_relative_shares_per_ssp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=3*600,res=600,pointsize=6)
          split.screen(c(1,3))
        } else {
          png(paste0(path.figs,"cmip6_","variance_relative_shares_per_ssp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=10)
          split.screen(c(2,2))
          
        }
        for(ssp in 1:length(ssps)){
          screen(ssp)
          par(mar=c(4,4,2,0.5))
          plot(c(2011:2084),totalvar6[ssp,]/totalvar6[ssp,],type="l",ylim=c(0,1.2),
               main=ssps.nice[ssp],xlab="year (AD)",ylab="relative contribution in overall variance")
          polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvar6)[2]),rev(meangcmvar6[ssp,]/totalvar6[ssp,])),col=col.var[1],border=NA)
          lower <- upper2 <- upper <- (meanggcmvar6[ssp,]+meangcmvar6[ssp,])/totalvar6[ssp,]
          upper[upper>1] <- 1
          polygon(c(2011:2084,2084:2011),c(meangcmvar6[ssp,]/totalvar6[ssp,],rev(upper)),col=col.var[2],border=NA)
          lower[upper2>1] <- 1
          upper2[upper2<1] <- 1
          polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
          
          #lines(c(2011:2084),meangcmvar6[ssp,]/totalvar6[ssp,],col=2)
          #lines(c(2011:2084),meanggcmvar6[ssp,]/totalvar6[ssp,],col=3)
          #lines(c(2011:2084),(meanggcmvar6[ssp,]+meangcmvar6[ssp,])/totalvar6[ssp,],col=4)
          legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        }
        close.screen(all=T)
        dev.off()
        
        png(paste0(path.figs,"cmip6_","variance_relative_shares_weightedprod_per_ssp_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        split.screen(c(2,2))
        for(ssp in 1:length(ssps)){
          screen(ssp)
          par(mar=c(4,4,2,0.5))
          plot(c(2011:2084),totalvarw6[ssp,]/totalvarw6[ssp,],type="l",ylim=c(0,1.2),
               main=paste("production weighted",ssps.nice[ssp]),xlab="year (AD)",ylab="relative contribution in overall variance")
          polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvarw6)[2]),rev(meangcmvarw6[ssp,]/totalvarw6[ssp,])),col=col.var[1],border=NA)
          lower <- upper2 <- upper <- (meanggcmvarw6[ssp,]+meangcmvarw6[ssp,])/totalvarw6[ssp,]
          upper[upper>1] <- 1
          polygon(c(2011:2084,2084:2011),c(meangcmvarw6[ssp,]/totalvarw6[ssp,],rev(upper)),col=col.var[2],border=NA)
          lower[upper2>1] <- 1
          upper2[upper2<1] <- 1
          polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
          
          #lines(c(2011:2084),meangcmvar6[ssp,]/totalvar6[ssp,],col=2)
          #lines(c(2011:2084),meanggcmvar6[ssp,]/totalvar6[ssp,],col=3)
          #lines(c(2011:2084),(meanggcmvar6[ssp,]+meangcmvar6[ssp,])/totalvar6[ssp,],col=4)
          legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        }
        close.screen(all=T)
        dev.off()
        
        
        png(paste0(path.figs,"cmip6_","variance_relative_shares_ssp85_leave1gcmout_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        plot(c(2011:2084),totalvar6[ssp,]/totalvar6[ssp,],type="l",ylim=c(0,1.2),
             main=ssps.nice[ssp],xlab="year (AD)",ylab="relative contribution in overall variance")
        polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvar6)[2]),rev(meangcmvar6[ssp,]/totalvar6[ssp,])),col=col.var[1],border=NA)
        lower <- upper2 <- upper <- (meanggcmvar6[ssp,]+meangcmvar6[ssp,])/totalvar6[ssp,]
        upper[upper>1] <- 1
        polygon(c(2011:2084,2084:2011),c(meangcmvar6[ssp,]/totalvar6[ssp,],rev(upper)),col=col.var[2],border=NA)
        lower[upper2>1] <- 1
        upper2[upper2<1] <- 1
        polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
        
        for(ii in 1:length(drop1gcm6))
        {
          delta_b <- drop1gcm6[[ii]]
          delta_b_meangcm <- drop1gcm_meangcm6[[ii]]
          delta_b_meanggcm <- drop1gcm_meanggcm6[[ii]]
          totalvar2 <- meangcmvar2 <- meanggcmvar2 <- array(NA,dim(delta_b)[c(4,5)])
          for(ssp in length(ssps)){
            for(i in 1:dim(totalvar6)[2]){
              totalvar2[ssp,i] <- var(as.vector(delta_b[cr,,,ssp,i]),na.rm=T)
              meangcmvar2[ssp,i] <- var(as.vector(delta_b_meangcm[cr,,ssp,i]),na.rm=T)
              meanggcmvar2[ssp,i] <- var(as.vector(delta_b_meanggcm[cr,,ssp,i]),na.rm=T)
            }
            lines(c(2011:2084),meangcmvar2[ssp,]/totalvar2[ssp,],col="darkgreen")
            lines(c(2011:2084),(meanggcmvar2[ssp,]+meangcmvar2[ssp,])/totalvar2[ssp,],col=4)
            text(2084,meangcmvar3[ssp,74]/totalvar3[ssp,74],gcms6[ii],adj=0,cex=0.5)
          }
        }
        legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        dev.off()
        
        png(paste0(path.figs,"cmip6_","variance_relative_shares_weightedprod_ssp85_leave1gcmout_",crops.param2[cr],sets,"_v",version,".png"),width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        plot(c(2011:2084),totalvarw6[ssp,]/totalvarw6[ssp,],type="l",ylim=c(0,1.2),
             main=paste("production weighted",ssps.nice[ssp]),xlab="year (AD)",ylab="relative contribution in overall variance")
        polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvarw6)[2]),rev(meangcmvarw6[ssp,]/totalvarw6[ssp,])),col=col.var[1],border=NA)
        lower <- upper2 <- upper <- (meanggcmvarw6[ssp,]+meangcmvarw6[ssp,])/totalvarw6[ssp,]
        upper[upper>1] <- 1
        polygon(c(2011:2084,2084:2011),c(meangcmvarw6[ssp,]/totalvarw6[ssp,],rev(upper)),col=col.var[2],border=NA)
        lower[upper2>1] <- 1
        upper2[upper2<1] <- 1
        polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
        
        for(ii in 1:length(dropw1gcm6))
        {
          delta_b <- dropw1gcm6[[ii]]
          delta_b_meangcm <- dropw1gcm_meangcm6[[ii]]
          delta_b_meanggcm <- dropw1gcm_meanggcm6[[ii]]
          totalvar2 <- meangcmvar2 <- meanggcmvar2 <- array(NA,dim(delta_b)[c(4,5)])
          for(ssp in length(ssps)){
            for(i in 1:dim(totalvar6)[2]){
              totalvar2[ssp,i] <- var(as.vector(delta_b[cr,,,ssp,i]),na.rm=T)
              meangcmvar2[ssp,i] <- var(as.vector(delta_b_meangcm[cr,,ssp,i]),na.rm=T)
              meanggcmvar2[ssp,i] <- var(as.vector(delta_b_meanggcm[cr,,ssp,i]),na.rm=T)
            }
            lines(c(2011:2084),meangcmvar2[ssp,]/totalvar2[ssp,],col="darkgreen")
            lines(c(2011:2084),(meanggcmvar2[ssp,]+meangcmvar2[ssp,])/totalvar2[ssp,],col=4)
            text(2084,meangcmvar3[ssp,74]/totalvar3[ssp,74],gcms6[ii],adj=0,cex=0.5)
          }
        }
        legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        dev.off()
        
        
        png(paste0(path.figs,"cmip6_","variance_relative_shares_ssp85_leave1ggcmout_",crops.param2[cr],sets,"_v",version,".png"),
            width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        plot(c(2011:2084),totalvar6[ssp,]/totalvar6[ssp,],type="l",ylim=c(0,1.2),axes=F,
             main=ssps.nice[ssp],xlab="year (AD)",ylab="relative contribution in overall variance")
        axis(1)
        axis(2)
        polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvar6)[2]),rev(meangcmvar6[ssp,]/totalvar6[ssp,])),col=col.var[1],border=NA)
        lower <- upper2 <- upper <- (meanggcmvar6[ssp,]+meangcmvar6[ssp,])/totalvar6[ssp,]
        upper[upper>1] <- 1
        polygon(c(2011:2084,2084:2011),c(meangcmvar6[ssp,]/totalvar6[ssp,],rev(upper)),col=col.var[2],border=NA)
        lower[upper2>1] <- 1
        upper2[upper2<1] <- 1
        polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
        
        for(ii in 1:length(drop1ggcm))
        {
          delta_b <- drop1ggcm[[ii]]
          delta_b_meangcm <- drop1ggcm_meangcm[[ii]]
          delta_b_meanggcm <- drop1ggcm_meanggcm[[ii]]
          totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
          for(ssp in length(ssps)){
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[ssp,i] <- var(as.vector(delta_b[cr,,,ssp,i]),na.rm=T)
              meangcmvar3[ssp,i] <- var(as.vector(delta_b_meangcm[cr,,ssp,i]),na.rm=T)
              meanggcmvar3[ssp,i] <- var(as.vector(delta_b_meanggcm[cr,,ssp,i]),na.rm=T)
            }
            lines(c(2011:2084),meangcmvar3[ssp,]/totalvar3[ssp,],col="darkgreen")
            lines(c(2011:2084),(meanggcmvar3[ssp,]+meangcmvar3[ssp,])/totalvar3[ssp,],col=4)
            text(2084,meangcmvar3[ssp,74]/totalvar3[ssp,74],ggcms[ii],adj=0,cex=0.5)
            text(2084,(meanggcmvar3[ssp,74]+meangcmvar3[ssp,74])/totalvar3[ssp,74],ggcms[ii],adj=0,cex=0.5)
          }
        }
        legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        
        dev.off()
        
        png(paste0(path.figs,"cmip6_","variance_relative_shares_weightedprod_ssp85_leave1ggcmout_",crops.param2[cr],sets,"_v",version,".png"),
            width=8*600,height=8*600,res=600,pointsize=12)
        par(mar=c(4,4,2,2.5),xpd=NA)
        plot(c(2011:2084),totalvarw6[ssp,]/totalvarw6[ssp,],type="l",ylim=c(0,1.2),axes=F,
             main=paste("production weighted",ssps.nice[ssp]),xlab="year (AD)",
             ylab="relative contribution in overall variance")
        axis(1)
        axis(2)
        polygon(c(2011:2084,2084:2011),c(rep(0,dim(totalvarw6)[2]),rev(meangcmvarw6[ssp,]/totalvarw6[ssp,])),col=col.var[1],border=NA)
        lower <- upper2 <- upper <- (meanggcmvarw6[ssp,]+meangcmvarw6[ssp,])/totalvarw6[ssp,]
        upper[upper>1] <- 1
        polygon(c(2011:2084,2084:2011),c(meangcmvarw6[ssp,]/totalvarw6[ssp,],rev(upper)),col=col.var[2],border=NA)
        lower[upper2>1] <- 1
        upper2[upper2<1] <- 1
        polygon(c(2011:2084,2084:2011),c(lower,rev(upper2)),col=col.var[3],border=NA)
        
        for(ii in 1:length(dropw1ggcm))
        {
          delta_b <- dropw1ggcm[[ii]]
          delta_b_meangcm <- dropw1ggcm_meangcm[[ii]]
          delta_b_meanggcm <- dropw1ggcm_meanggcm[[ii]]
          totalvar3 <- meangcmvar3 <- meanggcmvar3 <- array(NA,dim(delta_b)[c(4,5)])
          for(ssp in length(ssps)){
            for(i in 1:dim(totalvar3)[2]){
              totalvar3[ssp,i] <- var(as.vector(delta_b[cr,,,ssp,i]),na.rm=T)
              meangcmvar3[ssp,i] <- var(as.vector(delta_b_meangcm[cr,,ssp,i]),na.rm=T)
              meanggcmvar3[ssp,i] <- var(as.vector(delta_b_meanggcm[cr,,ssp,i]),na.rm=T)
            }
            lines(c(2011:2084),meangcmvar3[ssp,]/totalvar3[ssp,],col="darkgreen")
            lines(c(2011:2084),(meanggcmvar3[ssp,]+meangcmvar3[ssp,])/totalvar3[ssp,],col=4)
            text(2084,meangcmvar3[ssp,74]/totalvar3[ssp,74],ggcms[ii],adj=0,cex=0.5)
            text(2084,(meanggcmvar3[ssp,74]+meangcmvar3[ssp,74])/totalvar3[ssp,74],ggcms[ii],adj=0,cex=0.5)
          }
        }
        legend("topleft",legend=c("GGCM share","GCM share","cross-terms"),fill=col.var,bty="n",ncol=2)
        
        dev.off()
        
        png(paste0(path.figs,"cmip6_","global_delta",sets,"_ts_per_ssp_",crops.param2[cr],"_v",version,".png"),width=8*600,height=5*600,res=600,pointsize=12)
        plot(c(2011:2084),delta6[cr,1,1,1,],ylim=range(delta6[cr,,,,],na.rm=T),type="n",main=crops.nice[cr],
             xlab="year (AD)",ylab="change in global productivity [-]")
        for(ssp in 1:length(ssps)){
          med <- apply(delta6[cr,,,ssp,],3,median,na.rm=T)
          mi <- apply(delta6[cr,,,ssp,],3,min,na.rm=T)
          ma <- apply(delta6[cr,,,ssp,],3,max,na.rm=T)
          std <- apply(delta6[cr,,,ssp,],3,sd,na.rm=T)
          polygon(c(2011:2084,2084:2011),c(med-std,rev(med+std)),col=col.rcp2[ssp],border=NA)
          lines(c(2011:2084),med,col=col.rcp[ssp],lwd=2)
          lines(c(2011:2084),mi,col=col.rcp[ssp],lwd=1,lty=2)
          lines(c(2011:2084),ma,col=col.rcp[ssp],lwd=1,lty=2)
          # adding lines for ssp uncertainty ranges
          lines(rep(2084+ssp*0.6,2),c((med-std)[length(med)],(med+std)[length(med)]),col=col.rcp[ssp])
          # for(gcm in 1:length(gcms6)){
          #   for(ggcm in 1:length(ggcms)){
          #     lines(delta6[1,ggcm,gcm,ssp,],col=ggcm,lty=ssp)
          #   }
          # }
          abline(h=1,lty=2)
        }
        legend("topleft",legend=c(ssps.nice,"median","min/max","+/- 1SD",""),#x.intersp = c(1,1,1,1,1,1,-0.5,1)+1,
               col=c(col.rcp,1,1,NA,NA),lty=c(rep(1,5),2,NA,NA),lwd=c(rep(1,4),2,1,NA,NA),bty="n",
               fill=c(NA,NA,NA,NA,NA,NA,rgb(0,0,0,45/255),NA),border=NA,ncol=2)
        dev.off()
        
        png(paste0(path.figs,"cmip6_","global_delta",sets,"_ts_per_ssp_",crops.param2[cr],"_weightedprod_v",version,".png"),width=8*600,height=5*600,res=600,pointsize=12)
        plot(c(2011:2084),deltaw6[cr,1,1,1,],ylim=range(deltaw6[cr,,,,],na.rm=T),type="n",main=paste("production weighted",crops.nice[cr]),
             xlab="year (AD)",ylab="change in global productivity [-]")
        for(ssp in 1:length(ssps)){
          med <- apply(deltaw6[cr,,,ssp,],3,median,na.rm=T)
          mi <- apply(deltaw6[cr,,,ssp,],3,min,na.rm=T)
          ma <- apply(deltaw6[cr,,,ssp,],3,max,na.rm=T)
          std <- apply(deltaw6[cr,,,ssp,],3,sd,na.rm=T)
          polygon(c(2011:2084,2084:2011),c(med-std,rev(med+std)),col=col.rcp2[ssp],border=NA)
          lines(c(2011:2084),med,col=col.rcp[ssp],lwd=2)
          lines(c(2011:2084),mi,col=col.rcp[ssp],lwd=1,lty=2)
          lines(c(2011:2084),ma,col=col.rcp[ssp],lwd=1,lty=2)
          # adding lines for ssp uncertainty ranges
          lines(rep(2084+ssp*0.6,2),c((med-std)[length(med)],(med+std)[length(med)]),col=col.rcp[ssp])
          # for(gcm in 1:length(gcms6)){
          #   for(ggcm in 1:length(ggcms)){
          #     lines(delta6[1,ggcm,gcm,ssp,],col=ggcm,lty=ssp)
          #   }
          # }
          abline(h=1,lty=2)
        }
        legend("topleft",legend=c(ssps.nice,"median","min/max","+/- 1SD",""),#x.intersp = c(1,1,1,1,1,1,-0.5,1)+1,
               col=c(col.rcp,1,1,NA,NA),lty=c(rep(1,5),2,NA,NA),lwd=c(rep(1,4),2,1,NA,NA),bty="n",
               fill=c(NA,NA,NA,NA,NA,NA,rgb(0,0,0,45/255),NA),border=NA,ncol=2)
        dev.off()
        
        for(ssp in 1:length(ssps)){
          png(paste0(path.figs,"cmip6_","global_delta",sets,"_ts_",ssps.nice[ssp],"_",crops.param2[cr],"2_v",version,".png"),width=8*600,height=5*600,res=600,pointsize=12)
          plot(c(2011:2084),delta6[cr,1,1,1,],ylim=range(delta6[cr,,,,],na.rm=T),type="n",main=crops.nice[cr],
               xlab="year (AD)",ylab="change in global productivity [-]")
          par(mar=c(4,4,3,4),xpd=NA)
          for(ggcm in 1:length(ggcms)){
            if(all(!is.finite(delta6[cr,ggcm,,ssp,])))
              next
            med <- apply(delta6[cr,ggcm,,ssp,],2,median,na.rm=T)
            mi <- apply(delta6[cr,ggcm,,ssp,],2,min,na.rm=T)
            ma <- apply(delta6[cr,ggcm,,ssp,],2,max,na.rm=T)
            std <- apply(delta6[cr,ggcm,,ssp,],2,sd,na.rm=T)
            polygon(c(2011:2084,2084:2011),c(med-std,rev(med+std)),col=col.ggcms2[ggcm],border=NA)
            lines(c(2011:2084),med,col=col.ggcms[ggcm],lwd=2)
            lines(c(2011:2084),mi,col=col.ggcms[ggcm],lwd=1,lty=2)
            lines(c(2011:2084),ma,col=col.ggcms[ggcm],lwd=1,lty=2)
            # adding lines for ssp uncertainty ranges
            lines(rep(2084+ggcm*0.6,2),c((med-std)[length(med)],(med+std)[length(med)]),col=col.ggcms[ggcm])
            # for(gcm in 1:length(gcms6)){
            #   for(ggcm in 1:length(ggcms)){
            #     lines(delta6[1,ggcm,gcm,ssp,],col=ggcm,lty=ssp)
            #   }
            # }
          }
          par(xpd=FALSE)
          abline(h=1,lty=2)
          legend("bottomleft",legend=c(ggcms,"median","min/max","+/- 1SD",""),#x.intersp = c(1,1,1,1,1,1,-0.5,1)+1,
                 col=c(col.ggcms[1:length(ggcms)],1,1,NA,NA),lty=c(rep(1,12),2,NA,NA),lwd=c(rep(1,9),2,1,NA,NA),bty="n",
                 fill=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,rgb(0,0,0,45/255),NA),border=NA,ncol=3)
          dev.off()
          
          png(paste0(path.figs,"cmip6_","global_delta",sets,"_ts_",ssps.nice[ssp],"_",crops.param2[cr],"_weightedprod_v",version,".png"),width=8*600,height=5*600,res=600,pointsize=12)
          plot(c(2011:2084),deltaw6[cr,1,1,1,],ylim=range(deltaw6[cr,,,,],na.rm=T),type="n",main=paste("production weighted",crops.nice[cr]),
               xlab="year (AD)",ylab="change in global productivity [-]")
          par(mar=c(4,4,3,4),xpd=NA)
          for(ggcm in 1:length(ggcms)){
            if(all(!is.finite(deltaw6[cr,ggcm,,ssp,])))
              next
            med <- apply(deltaw6[cr,ggcm,,ssp,],2,median,na.rm=T)
            mi <- apply(deltaw6[cr,ggcm,,ssp,],2,min,na.rm=T)
            ma <- apply(deltaw6[cr,ggcm,,ssp,],2,max,na.rm=T)
            std <- apply(deltaw6[cr,ggcm,,ssp,],2,sd,na.rm=T)
            polygon(c(2011:2084,2084:2011),c(med-std,rev(med+std)),col=col.ggcms2[ggcm],border=NA)
            lines(c(2011:2084),med,col=col.ggcms[ggcm],lwd=2)
            lines(c(2011:2084),mi,col=col.ggcms[ggcm],lwd=1,lty=2)
            lines(c(2011:2084),ma,col=col.ggcms[ggcm],lwd=1,lty=2)
            # adding lines for ssp uncertainty ranges
            lines(rep(2084+ggcm*0.6,2),c((med-std)[length(med)],(med+std)[length(med)]),col=col.ggcms[ggcm])
            # for(gcm in 1:length(gcms6)){
            #   for(ggcm in 1:length(ggcms)){
            #     lines(delta6[1,ggcm,gcm,ssp,],col=ggcm,lty=ssp)
            #   }
            # }
          }
          par(xpd=FALSE)
          abline(h=1,lty=2)
          legend("topleft",legend=c(ggcms,"median","min/max","+/- 1SD",""),#x.intersp = c(1,1,1,1,1,1,-0.5,1)+1,
                 col=c(col.ggcms[1:length(ggcms)],1,1,NA,NA),lty=c(rep(1,12),2,NA,NA),lwd=c(rep(1,9),2,1,NA,NA),bty="n",
                 fill=c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,rgb(0,0,0,45/255),NA),border=NA,ncol=2)
          dev.off()
          
        }#ssp
      }#cr
    }
  }#sets
  
} #if(F)

## end of script, always close worker nodes if run in parallel mode
if(parallel==T){
  closeCluster(clu)
  ## MPI-save version of R quit 
  mpi.quit()
}

