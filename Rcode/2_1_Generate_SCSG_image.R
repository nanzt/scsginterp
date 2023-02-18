# A stepwise framework for interpolating land surface temperature under cloudy
# conditions based on the solar-cloud-satellite geometry
# preprocessing: filter MODIS cloud height data 
# authors: Yuhong Chen, Zetao Cao
# last modification date: 2023-02-13
# group website: https://permalab.science/


library(raster)
library(sp)
library(tiff)
library(ggplot2)
library(MASS)
library(spData)
library(wql)
library(Metrics)
library(gstat)
library(lubridate)
library(stringi)
library(stringr)
library(wavelets)
library(waveslim)
library(cluster)
library(clusterSim)
library(splines)
library(NISTunits) # Unit conversion
library(configr)

# read config
configFile <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/Rcode/config.ini"
configList <- read.config(file = configFile)
#platform
platform <- "Aqua"


# input variables
## land surface temperature
input.LSTDir <- configList$SCSG$LST
## could top height folder
input.CTHDir <- configList$SCSG$CouldTopHeight
## sensor azimuth folder
input.SenAziDir <- configList$SCSG$SenAzi
## sensor zenith folder
input.SenZenDir <- configList$SCSG$SenZen
## solar azimuth folder
input.SolAziDir <- configList$SCSG$SolAzi
## solar zenith folder
input.SolZenDir <- configList$SCSG$SolZen
## Altitude folder
input.AltiDir <- configList$SCSG$DEM

#output data folder
output.SCSGdir= configList$SCSG$SCSGOutputDir

# Check output path.
if (file.exists(output.SCSGdir)==FALSE){
  dir.create(output.SCSGdir)
}


# read data files
LSTFiles <- list.files(input.LSTDir, "*.tif$")
CTHFiles <- list.files(input.CTHDir,"*.tif$")
SenAziFiles <- list.files(input.SenAziDir,"*.tif$")
SenZenFiles <- list.files(input.SenZenDir,"*.tif$")
SolAziFiles <- list.files(input.SolAziDir,"*.tif$")
SolZenFiles <- list.files(input.SolZenDir,"*.tif$")
AltiFile <- list.files(input.AltiDir,"*.tif$")
#common date
LSTdate <- str_extract(LSTFiles, '\\d{7}')
MOD06date <- str_extract(CTHFiles, '\\d{7}')
common_date <- Reduce(intersect,list(LSTdate,MOD06date))


# global variables
{
  # Altitude (DEM) Raster
  AltitudeRas <- raster(paste0(input.AltiDir, "/", AltiFile[1]))
  # Build index by LST raster
  IndexRas <- raster(paste0(input.LSTDir,"/",LSTFiles[1]))
  IndexRas[1:length(IndexRas)] <- 1:length(IndexRas)
}

# function: Calculate SCSG
GetSCSGRaster <- function(LST, SolZen, SolAzi, SenZen, SenAzi, CTH, Altitude, IndexRas, dateStr, savePath){
  
  # LST raster for subsequent calculation
  LSTRas <- LST
  ShadowRas <- LST
  
  # get x,y location from raster
  XY <- rasterToPoints(Altitude)
  X <- XY[,1]
  Y <- XY[,2]
  
  # convert raster to vector
  LST <- as.vector(LST)
  CTH <- as.vector(CTH)
  
  SenAzi <- as.vector(SenAzi)
  SenZen <- as.vector(SenZen)
  SolAzi <- as.vector(SolAzi)
  SolZen <- as.vector(SolZen)
  Altitude <- as.vector(Altitude)
  
  # drop pixels altitude unknown
  index <- 1:length(LST)
  DF <- data.frame (index, LST, Altitude)
  DF <- DF[!is.na(Altitude),]
  DF <- cbind(DF,X,Y)
  
  # missing and known LST
  LSTMiss <- DF[is.na(DF$LST),]
  LSTKnown <- DF[!is.na(DF$LST),]
  ## index of missing LST
  indexMiss <- LSTMiss$index
  
  # all the data
  DFAll <- data.frame(index,LST,CTH,SenAzi,SenZen,SolAzi,SolZen,Altitude)
  # data at locations with missing LST
  DFAllMiss <- subset(DFAll, DFAll$index %in% indexMiss)
  DFAllMiss <- within(DFAllMiss, {
    Y <- LSTMiss$Y
    X <- LSTMiss$X
  })
  DFAllMiss <- DFAllMiss[!is.na(DFAllMiss$SolZen),]
  
  #scale degree and convert it to radian
  DFAllMiss <- within(DFAllMiss,{
    SolZen <- NISTdegTOradian(SolZen)
    SolAzi <- NISTdegTOradian(SolAzi)
    SenZen <- NISTdegTOradian(SenZen)
    SenAzi <- NISTdegTOradian(SenAzi)
  })
  
  # missing LST but cloud top height is lower than altitude
  DFAllMiss_couldLower <- subset(DFAllMiss, DFAllMiss$CTH <= DFAllMiss$Altitude )
  
  # missing LST due to cloud
  DFAllMiss <- subset(DFAllMiss, DFAllMiss$CTH >= DFAllMiss$Altitude)
  IndexCloudLower <- DFAllMiss_couldLower$index
  
  # calculate shadow pixels
  ## cloud position in the original image
  XOri <- DFAllMiss$X
  YOri <- DFAllMiss$Y
  ## relative cloud height to the surface
  HRe <- DFAllMiss$CTH - DFAllMiss$Altitude
  ## angle information
  theta_v <- DFAllMiss$SenZen
  phi_v   <- DFAllMiss$SenAzi
  theta_s <- DFAllMiss$SolZen
  phi_s   <- DFAllMiss$SolAzi
  ## cloud orthographic projection on the surface
  XCloud <- XOri + HRe * tan(theta_v) * sin(phi_v)
  YCloud <- YOri + HRe * tan(theta_v) * cos(phi_v)
  ## cloud shadow
  XShadow <- XCloud - HRe * tan(theta_s) * sin(phi_s)
  YShadow <- YCloud - HRe * tan(theta_s) * cos(phi_s)
  # print (abs(XOri - XShadow))
  # print (abs(YOri - YShadow))
  # get LST and index of shadow pixels
  ## coordinates of shadow pixels
  XYShadow <- cbind(XShadow, YShadow)
  ## extract LST and index of shadow pixels
  IndexShadow <- extract(IndexRas,XYShadow)
  LSTShadow <- extract(LSTRas,XYShadow)
  
  # identify ABCD regions
  '%notin%' <- Negate('%in%')
  
  
  ## data frame for shadow pixels
  ShadowPixels <- data.frame(IndexShadow,LSTShadow)
  
  
  # region D
  D1 <- subset(ShadowPixels, !is.na(ShadowPixels$IndexShadow)
               & !is.na(ShadowPixels$LSTShadow))
  D <- subset(LSTKnown, LSTKnown$index %in% D1$IndexShadow)
  
  # region A
  A <- subset(LSTKnown, LSTKnown$index %notin% D$index)
  
  # region C
  C1 <- subset(ShadowPixels, !is.na(ShadowPixels$IndexShadow)
               & is.na(ShadowPixels$LSTShadow))
  C <- subset (LSTMiss, LSTMiss$index %in% C1$IndexShadow & LSTMiss$index %notin% IndexCloudLower)
  
  # region B
  B <- subset(LSTMiss, LSTMiss$index %notin% C$index)
  
  # write to raster
  LSTRas[A$index] <- 1
  LSTRas[B$index] <- 2
  LSTRas[C$index] <- 3
  LSTRas[D$index] <- 4
  LSTRas[IndexCloudLower] <- 3
  writeRaster(LSTRas, filename = 
                paste0(savePath,"/",platform, "_SCSG_ABCD_", dateStr, ".tif"), overwrite = TRUE )
  
}



# apply to all images
for ( i in 1:length(common_date)){
  
  # get date
  dateStr <- str_extract(common_date[i], '\\d{7}')
  print (dateStr)
  # read file as raster
  LSTRas    <- raster(paste0(input.LSTDir, "/", LSTFiles[i]))
  SolZenRas <- raster(paste0(input.SolZenDir, "/", SolZenFiles[grep(dateStr,SolZenFiles)]))
  SolAziRas <- raster(paste0(input.SolAziDir ,"/", SolAziFiles[grep(dateStr,SolAziFiles)]))
  SenZenRas <- raster(paste0(input.SenZenDir, "/", SenZenFiles[grep(dateStr,SenZenFiles)]))
  SenAziRas <- raster(paste0(input.SenAziDir,"/", SenAziFiles[grep(dateStr,SenAziFiles)]))
  CTHRas    <- raster(paste0(input.CTHDir, "/", CTHFiles[grep(dateStr,CTHFiles)]))
  
  #SCSG for each day
  GetSCSGRaster(LST = LSTRas, SolZen = SolZenRas, SolAzi = SolAziRas, SenZen = SenZenRas, SenAzi = SenAziRas
                , CTH = CTHRas, Altitude = AltitudeRas, IndexRas = IndexRas
                , dateStr = dateStr, savePath = output.SCSGdir)
}





