# A stepwise framework for interpolating land surface temperature under cloudy
# conditions based on the solar-cloud-satellite geometry
# preprocessing: filter MODIS cloud height data 
# authors: Yuhong Chen, Zetao Cao
# last modification date: 2023-02-13
# group website: https://permalab.science/


library(sp)
library(raster)
library(stringr)
library(configr)


# read config
configFile <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/Rcode/config.ini"
configList <- read.config(file = configFile)


# Input variables
## MODIS view time data
input.ViewTimeDir <- configList$DSRViewTime$ViewTime
## Surface downward shortwave radiation data
input.DSRDir <- configList$DSRViewTime$DSR
## Time difference between UTC and local time
input.TimeDiff <- as.numeric(configList$DSRViewTime$TimeDiff)


# Output Variables
## DSR with passtime results
output.DSR_PasstimeDir <- configList$DSRViewTime$DSRPasstime

# Check output path.
if (file.exists(output.DSR_PasstimeDir)==FALSE){
  dir.create(output.DSR_PasstimeDir)
}

# Read data (UTC)
ViewTimeFiles <- list.files(path=input.ViewTimeDir, pattern = "*.tif$", full.names = FALSE)
# Read DSR data (030/06 represent UTC 03:00/06:00)
DSR03Files <- list.files(path = input.DSRDir, pattern = "*DSR03.tif$", full.names = FALSE)
DSR06Files <- list.files(path = input.DSRDir, pattern = "*DSR06.tif$", full.names = FALSE)

# Function:: Get DSR with MODIS passtime
GetDSR_MODISPasstime <- function(ViewTimeRaster, DSR03Raster, DSR06Raster, OutPath, date){
  #Get Raster template
  Raster_Template <- DSR03Raster
  
  # convert Raster data to vectors
  ViewTime_vec <- as.vector(ViewTimeRaster)
  DSR03_vec <- as.vector(DSR03Raster)
  DSR06_vec <- as.vector(DSR06Raster)
  
  #convert vectors above to data frame matrix
  ## set index
  index1 <- 1:length(ViewTime_vec)
  ## create a data frame matrix (scale time)
  M <- data.frame(index1, ViewTime_vec*0.1, DSR03_vec, DSR06_vec)
  ## get nonempty subset
  M <- subset(M, !is.na(ViewTime_vec))
  M <- subset(M, !is.na(DSR03_vec))
  M <- subset(M, !is.na(DSR06_vec))

  ViewTime_mat=as.matrix(M$ViewTime_vec)
  DSR06_mat=as.matrix(M$DSR06_vec)
  DSR03_mat=as.matrix(M$DSR03_vec)
  TimeStep <- ViewTime_mat-9
  Slope <- (DSR06_mat-DSR03_mat)/3
  DSMatPassT_mat <-TimeStep*Slope+DSR03_mat
  DSMatPassT_vec=as.vector(DSMatPassT_mat)
  Raster_Template[M$index1] <- DSMatPassT_vec
  writeRaster(Raster_Template, filename = paste0(output.DSR_PasstimeDir,"/", date, "_DSR_MODIS_overpasstime.tif"))
}


# apply to all images
for (i in 1:length(ViewTimeFiles)){
  #read a view time file
  ViewTimeFile <- ViewTimeFiles[i]
  #get date
  dateStr <- str_extract(ViewTimeFile, pattern = "\\d{7}")
  print (dateStr)
  #get corresponding DSR data for the date
  DSR03File<- grep(dateStr, DSR03Files, value=TRUE)
  DSR06File<- grep(dateStr, DSR06Files, value=TRUE)
  
  #read as raster
  ViewTimeRaster <- raster(paste0(input.ViewTimeDir,"/",ViewTimeFile))
  DSR03Raster <- raster(paste0(input.DSRDir,"/",DSR03File))
  DSR06Raster <- raster(paste0(input.DSRDir,"/",DSR06File))
  GetDSR_MODISPasstime(ViewTimeRaster, DSR03Raster, DSR06Raster, output.DSR_PasstimeDir, dateStr)
  
}

