# A stepwise framework for interpolating land surface temperature under cloudy
# conditions based on the solar-cloud-satellite geometry
# preprocessing: filter MODIS cloud height data 
# authors: Yuhong Chen, Zetao Cao
# last modification date: 2023-02-13
# group website: https://permalab.science/


library(raster)
library(sp)
library(tiff)
library(igraph)
library(ggplot2)
library(MASS)
library(spData)
library(spdep)
# library(wql)
library(Metrics)
require(gstat)
library(lubridate)
library(stringi)
library(wavelets)
library(waveslim)
library(cluster)
library(class)
library(clusterSim)
library(stringr)
library(parallel)
library(configr)

# read config
configFile <- "D:/cyh/ISPRS_LST interpolation_cyh/Rcode/config.ini"
configList <- read.config(file = configFile)

# itpLST
itpLSTDir <-  configList$clearSkyInterpolation$originalLST
itpLSTFiles <- list.files(path = itpLSTDir,full.names = T,pattern = "*.tif$")
refLSTDir <- configList$clearSkyInterpolation$referenceLST

# output dir 
csvDir <- configList$clearSkyInterpolation$outCSV

# ndvi data
NDVIDir <- configList$clearSkyInterpolation$ndvi
NDVIFiles <- list.files(path = NDVIDir,full.names = T,pattern = "*.tif$")
NDVIDates <- as.Date(str_extract(NDVIFiles, '\\d{7}'), format = "%Y%j")

# altitude raster
DEMRas <- raster(configList$clearSkyInterpolation$DEMRas)

# SCSG directory
SCSGDir <-configList$clearSkyInterpolation$SCSGDir

# some threshold values
# common pixels with valid values should be more than this proportion
refproportion <- as.numeric(configList$clearSkyInterpolation$refproportion)
# search days
timeWinLen <- as.numeric(configList$clearSkyInterpolation$timeWinLen)
# define the distance smaller as similar pixels
threshold_similar <- as.numeric(configList$clearSkyInterpolation$similarPixelNum)

# cup core number
clnum<-detectCores() 
cpuCNum <- clnum-2


# funtion: find LST files on the near days 
GetNearLSTs <- function(dateStr, LSTDir,timeWinLen){
  
  L <- (timeWinLen-1)%/%2L
  date <- as.Date(dateStr,format = "%Y%j")
  # start date
  date1 <- date-L
  nearLSTFiles <- c()
  for (i in c((0:(L-1)),((L+1):(timeWinLen-1)))){
    nearDateStr <- as.character(date1+i,format='%Y%j')
    print (nearDateStr)
    nearLSTFile <- intersect(list.files(path = LSTDir, pattern=nearDateStr, full.names = T)
                             ,list.files(path = LSTDir, pattern=".tif$" , full.names = T))
    
    
    if (length(nearLSTFile)!=0){
      print ("Got LST file!")
      nearLSTFiles <- c(nearLSTFiles,nearLSTFile)
    }
    else {
      print ("No LST file this day.")
    }
  }
  return (nearLSTFiles)
}

# funtion: find the closest NDVI file 
GetClosestNDVI<-function(NDVIDates, dateStr){
  date <- as.Date(dateStr,format = "%Y%j")
  dateDiffs <- abs(NDVIDates-date)
  return (which.min(dateDiffs))
  
}

# function: EOF interpolation for a LST file based on a reference images
EOFforImage <- function(itpLSTRas,refLSTRas, DEMRas,NDVIRas,SCSGRas,dateStr,refDateStr){
  
  # set raster template
  RasTemplate <- itpLSTRas
  
  #raster to vector
  itpLST <<- as.vector(itpLSTRas)
  refLST <<- as.vector(refLSTRas)
  DEM <- as.vector(DEMRas)
  NDVI <- as.vector(NDVIRas)
  SCSG <- as.vector(SCSGRas)
  
  index <- 1:length(itpLST)
  
  refLSTnor <- refLST # will be normalized later
  
  refNum <- refproportion*length(itpLST)
  
  # build data frame
  df.all <- data.frame(index,itpLST,refLST, refLSTnor, DEM ,NDVI, SCSG)
  # to make sure there are enough data pairs used to interpolate
  if (length(subset(df.all, !is.na(refLST) & !is.na(itpLST))$index) < refNum){
    print(paste("Common valid pixels of",itpDateStr, "and", refDateStr, "are not enough."))
    return(NA)
  }
  else{
    # drop no data
    df.all <- subset(df.all, 
                     !is.na(refLST) 
                     & !is.na(DEM)
                     & !is.na(NDVI)
                     & !is.na(SCSG))
    # normalization 
    df.all[,4:6] <- data.Normalization(df.all[,4:6], type = "n4", normalization = "column")
    
    # pixels to be interpolated 
    df.itp <- subset(df.all, is.na(itpLST)|SCSG==4)
    
    # reference pixels with valid LST values in both images
    df.ref <- subset(df.all, !is.na(itpLST) & SCSG == 1) 
    
    # reference data for calculate similarity
    ref <<- df.ref[,4:6]
    refIndex <<- df.ref$index
    #print("ref done")
    
    # function: single pixel interpolation
    EOFforPixel <- function(oneDf, ref){
      
      # calculate Euclidean distance between pixels
      eDist <- sqrt( (ref$refLSTnor - oneDf[4])^2
                     + (ref$DEM- oneDf[5])^2
                     + (ref$NDVI - oneDf[6])^2
      )
      
      # select similar pixels: first 1000 points
      Inter_value <- quantile(eDist, threshold_similar/length(eDist))
      conditon <- (eDist < Inter_value)
      similar_index <- refIndex[conditon]
      
      #eof interpolation
      eof_data <- cbind(itpLST[similar_index], refLST[similar_index])
      colnames(eof_data) <- c("itpLST","refLST")
      eof_data <- rbind(oneDf[2:3], eof_data)
      
      mean_int <- mean(eof_data[,1],na.rm=TRUE)
      mean_ref <- mean(eof_data[,2],na.rm=TRUE)
      
      #data processing
      eof_data[,1] <- eof_data[,1] - mean_int
      eof_data[,2] <- eof_data[,2] - mean_ref
      eof_data[1,1] <- 0
      
      #svd
      s <- svd(eof_data)
      U <- s$u
      D <- s$d
      V <- s$v
      s_recon <- (U[,1] * D[1]) %*% t(V[,1])
      s_recon[,1] <- s_recon[,1] + mean_int
      s_recon[,2] <- s_recon[,2] + mean_ref
      #
      eof_result <- s_recon[1,1] 
      
      #get var
      s_recon <- s_recon[-1,]
      ref_var <- var(refLST[similar_index]-s_recon[,2])
      
      #result
      pixel_result <- c(eof_result,ref_var)
      return(pixel_result)
    }
    
    #cpu core number
    cl <- makeCluster(cpuCNum)
    clusterExport(cl,"ref")
    clusterExport(cl,"refIndex")
    clusterExport(cl,"itpLST")
    clusterExport(cl,"refLST")
    clusterExport(cl,"threshold_similar")
    
    # apply to all missing data
    #oneDf <- as.vector(as.matrix(df.itp[1,]))
    PredictValues <- parApply(cl, df.itp,1, EOFforPixel, ref)
    PredictValues <- t(PredictValues)
    stopCluster(cl)
    
    
    df.itp <- cbind(df.itp$index, PredictValues)
    
    colnames(df.itp) <- c("index","predict","ref_var")
    
    outDir <- paste0(csvDir,"/",dateStr)
    if (!dir.exists(outDir)){
      dir.create(outDir)
    }
    
    
    write.csv(df.itp,file = paste0(outDir,"/",
                                   dateStr,"_eof_var_",refDateStr,"_ref.csv"))
    print(paste0("csv completed: ", dateStr,"_eof_var_",refDateStr,"_ref.csv"))
  }
}


#batch interpolation
for ( i in 1:length(itpLSTFiles)){
  
  #a single interpolated image
  itpLSTFile <- itpLSTFiles[i]
  dateStr <- str_extract(itpLSTFile, '\\d{7}')
  
  print (paste("BEGIN*************",dateStr,"*************"))
  
  #reference image list
  itpLSTRas <- raster(itpLSTFile)
  refLSTFiles <- GetNearLSTs(dateStr,refLSTDir,timeWinLen)
  
  #SCSG file
  SCSGFile <- intersect(list.files(path = SCSGDir, pattern=dateStr , full.names = T)
                        , list.files(path = SCSGDir, pattern=".tif$" , full.names = T))
  print  (paste("SCSG:",str_extract(SCSGFile, '\\d{7}')))
  SCSGRas <- raster(SCSGFile)
  
  #NDVI
  NDVIIndex <- GetClosestNDVI(NDVIDates,dateStr)
  NDVIFile <- NDVIFiles[NDVIIndex]
  print ((paste("NDVI:", as.character(NDVIDates[NDVIIndex]))))
  NDVIRas <- raster(NDVIFile)
  
  #EOF interpolation based on multiple reference images
  for (i in 1:length(refLSTFiles)){
    
    #a reference image
    print("--------------")
    refDateStr <- str_extract(refLSTFiles[i], '\\d{7}')
    print (paste("reference:", refDateStr))
    refLSTRas <- raster(refLSTFiles[i])
    
    #EOF interpolation
    EOFforImage(itpLSTRas,refLSTRas, DEMRas,NDVIRas,SCSGRas,dateStr,refDateStr)
    
  }
  print (paste("END*************",dateStr,"*************"))
  print (paste(" "))
  
}

