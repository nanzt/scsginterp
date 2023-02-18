# A stepwise framework for interpolating land surface temperature under cloudy
# conditions based on the solar-cloud-satellite geometry
# preprocessing: filter MODIS cloud height data 
# authors: Yuhong Chen, Zetao Cao
# last modification date: 2023-02-13
# group website: https://permalab.science/

# Step1: deal with missing data
# Step2: deal with small patches 
# Step3: deal with gaps caused by orbits


library(raster)
library(sp)
library(tiff)
library(igraph)
library(ggplot2)
library(MASS)
library(spData)
library(spdep)
library(wql)
library(Metrics)
require(gstat)
library(lubridate)
library(stringi)
library(wavelets)
library(waveslim)
library(cluster)
library(clusterSim)
library(splines)
library(NISTunits)
library(stringr)
library(modeest)
library(configr)


# read config
configFile <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/Rcode/config.ini"
configList <- read.config(file = configFile)


#setpath
setwd(configList$SCSGPostProcess$PathDir)
# output folder
output.ABCDdir <- configList$SCSGPostProcess$SCSGPPdir
output.oribitDir <- "1_outputOrbit"
output.missingDir <- "2_outputMissing"
output.patchDir <- "3_outputPatch"

# input variables
input.SCSGDir <- configList$SCSGPostProcess$SCSGDir
input.CTHDir <- configList$SCSGPostProcess$CouldTopHeight
input.SenAziDir <- configList$SCSGPostProcess$SenAzi
input.windowSize <- as.numeric(configList$SCSGPostProcess$windowSize)

# Check output path.
{
  if (file.exists(output.ABCDdir)==FALSE){
  dir.create(output.ABCDdir)
  }
  if (file.exists(output.missingDir)==FALSE){
    dir.create(output.missingDir)
  }
  if (file.exists(output.patchDir)==FALSE){
    dir.create(output.patchDir)
  }
  if (file.exists(output.oribitDir)==FALSE){
    dir.create(output.oribitDir)
  }
}

# get data files
SCSGFiles <- list.files(path= input.SCSGDir, pattern = "*.tif$", full.names = F)
CTHFiles <- list.files(path = input.CTHDir, pattern = "*.tif$", full.names = F)
SenAziFiles <- list.files(path = input.SenAziDir, pattern = "*.tif$", full.names = F)


### FUNCTIONS
{
  # Function: handle B-pixels by missing values of cloud and observation angle (represented by sensor azimuth)
  ProcessByMissingData <- function(SCSGRaster, CTHRaster, angleRaster, SavePath){
    
    #convert raster to raster
    SCSGVec <- as.vector(SCSGRaster)
    CTHVec <- as.vector(CTHRaster)
    angleVec <- as.vector(angleRaster)
    #convert to data frame 
    index <- 1:length(SCSGVec)
    DF <- data.frame(index, SCSGVec,CTHVec, angleVec)
    # get data gaps of cloud and angle
    DF <- subset(DF, is.na(CTHVec) | is.na(angleVec))
    # get B-pixels
    DF <- subset(DF, SCSGVec == 2)
    # convert these B-pixels to C
    SCSGRaster[DF$index] <- 3
    
    # save
    writeRaster(SCSGRaster, filename = paste0(SavePath, "/", names(SCSGRaster), "_BProcessedBymissing.tif"), overwrite = T)
    return (SCSGRaster)
  }
  
  # function: handle small patches
  ProcessSmallPatch <- function(SCSG_raster, windowSize,save_path){
    
    #function: judge whether the pixel is patch
    getMode <- function(mat){
      if (mat[13]!=2 ){
        return (mat[13])
      }
      else{
        y <- table(mat)
        B_count<-as.numeric(y["2"])
        C_count <- as.numeric(y["3"])
        if (B_count <= C_count & !is.na(C_count)){
          return (3)
        }
        else{
          return (2)
        }
      }
    }
    
    winMatrix <- matrix(1,nrow = windowSize, ncol = windowSize)
    
    ModeFocalRas <- focal(SCSG_raster, winMatrix, fun = getMode, pad = TRUE, padValue = 0)
  
    #write raster
    writeRaster(ModeFocalRas, filename = paste0(save_path, "/", names(SCSG_raster), "_PostProcess.tif"), overwrite = T)
    
    return(SCSG_raster)
  }
  
  # function: find orbit gap index
  FindOrbitGapIndex <- function(SCSGRas, angleRas){
    
    
    #build index matrix 
    IndexMat <- matrix(1:length(angleRas), nrow = nrow(angleRas), byrow = T)
    
    #angle data
    angleMat <- as.matrix(angleRas)
    angleRightMat <- cbind(angleMat[,2:ncol(angleMat)], angleMat[,1])
    diffMat <- abs(angleMat-angleRightMat)
    diffMat[,ncol(diffMat)] <- NA
    
    # find gap pixels
    FindGapPixels <- function(eachRow){
      if(length(eachRow[is.na(eachRow)]) == length(eachRow)){
        return (NA)
      }
      else{
        maxAngleDiff <- max(eachRow, na.rm = T)
        if (!is.na(maxAngleDiff) & maxAngleDiff > 100){
          return (which(eachRow == maxAngleDiff))
        }
        else{
          return (NA)
        }
      }
    }
    colIndex <- unlist(apply(diffMat, 1, FindGapPixels))
    
    # find NA pixels
    FindNAPixels <- function(eachRow){
      indexVex <- c(1:length(eachRow))
      if(length(eachRow[is.na(eachRow)]) == length(eachRow)){
        return (NA)
      }
      else if (length(eachRow[is.na(eachRow)]) < 3){
        return (NA)
      }
      else{
        return ( indexVex[is.na(eachRow)][1])
      }
    }
    
    
    colIndex <- unlist(apply(diffMat, 1, FindGapPixels))
    colIndexNa <- unlist(apply(diffMat, 1, FindNAPixels))
    
    # find corresponding index in the raster ï¼?!!!!!!RASTER 
    orbitGapLocations <- cbind(1:nrow(diffMat), colIndex)
    orbitNALocations <- cbind(1:nrow(diffMat), colIndexNa)
    gapSize=10
    orbitNALocations2 <- cbind(1:(nrow(diffMat)-gapSize), colIndexNa[(gapSize+1):length(colIndexNa)])
    
    orbitGapLocations <- rbind(orbitGapLocations,orbitNALocations)
    orbitGapLocations <- rbind(orbitGapLocations,orbitNALocations2)
    
    orbitGapIndex <- c()
    for ( i in 1:(length(colIndex)+length(colIndexNa)+length(colIndexNa)-gapSize)){
      orbitGapIndex[i] <- IndexMat[orbitGapLocations[i,1],orbitGapLocations[i,2]]
    }
    
    
    # drop nodata
    orbitGapIndex <- orbitGapIndex[!is.na(orbitGapIndex)]
    return (orbitGapIndex)
  }
  
  # funtion: Fill Gap With C regions
  FillGapWithCregion <- function(orbitIndex, SCSGRas, SavePath){
    maxLength <- length(SCSGRas)
    for (i in 1:length(orbitIndex)){
      if (SCSGRas[orbitIndex[i]] == 2){
        leftGapBoundary <- orbitIndex[i]-1
        rightGapBoundary <- orbitIndex[i]+1
        #find gap boundary
        while (SCSGRas[leftGapBoundary] == 2 & leftGapBoundary>=0){
          leftGapBoundary <- leftGapBoundary - 1
        }
        while (SCSGRas[rightGapBoundary] == 2 & rightGapBoundary <= maxLength){
          rightGapBoundary <- rightGapBoundary + 1
        }
        # set gap C region
        SCSGRas[leftGapBoundary:rightGapBoundary] <- 3
      }
      else{
        next
      }
    }
    writeRaster(SCSGRas, filename = paste0(SavePath, "/", names(SCSGRas), "_fillOrbitGap.tiff"), overwrite = T)
    return (SCSGRas)
  }
}

# Batch process SCSG files
for (i in 1:length(SCSGFiles)){
  
  # read raster
  SCSGRas <- raster(paste0(input.SCSGDir, "/", SCSGFiles[i]))
  DateStr <- str_extract(names(SCSGRas), '\\d{7}')
  print (DateStr)
  # get angle and cloud data
  angleRas <- raster(paste0(input.SenAziDir, "/", SenAziFiles[grep(DateStr, SenAziFiles)]))
  CTHRas <- raster(paste0(input.CTHDir, "/", CTHFiles[grep(DateStr, CTHFiles)]))
  #shadowRas <- raster(paste0(input.ShadowDir, "/", ShadowFiles[grep(DateStr, ShadowFiles)]))
  
  ##Steps
  
  #Step1: deal with gaps caused by orbits
  #get orbit gap index
  print( "Step1: deal with gaps caused by orbits. ")
  orbitIndex <- FindOrbitGapIndex(SCSGRas = SCSGRas, angleRas = angleRas)
  
  if (length(orbitIndex) == 0){
    print("No orbit gap.")
  }
  else{
    SCSGRas <- FillGapWithCregion(orbitIndex = orbitIndex, SCSGRas = SCSGRas, SavePath = output.oribitDir)
  }
  
  #Step2: deal with missing data
  print( "Step2: deal with missing data.")
  SCSGRas <- ProcessByMissingData(SCSGRaster = SCSGRas, CTHRaster = CTHRas, angleRaster = angleRas, output.missingDir)
  
  #Step3: deal with small patches 
  print( "Step3: deal with small patches.")
  SCSGRas <- ProcessSmallPatch(SCSG_raster = SCSGRas, windowSize = input.windowSize,save_path = output.ABCDdir)
  
}






