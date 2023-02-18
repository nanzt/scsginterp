# A stepwise framework for interpolating land surface temperature under cloudy
# conditions based on the solar-cloud-satellite geometry
# preprocessing: filter MODIS cloud height data 
# authors: Yuhong Chen, Zetao Cao
# last modification date: 2023-02-13
# group website: https://permalab.science/


# libraries
{ library(raster)
  library(sp)
  library(MASS)
  library(spData)
  library(wql)
  library(Metrics)
  require(gstat)
  library(lubridate)
  library(wavelets)
  library(waveslim)
  library(cluster)
  library(class)
  library(clusterSim)
  library(stringr)
  library(parallel)
  library(dplyr)
  library(configr)
  library(mda)}

# cup core number
clnum<-detectCores() 
cpuCNum <- 2

# read config
configFile <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/Rcode/config.ini"
configList <- read.config(file = configFile)

# K-mean cluster number
KClusterNum <- configList$cloudInterpolation$KClusterNum


#parameters for correcting LST values in region D
#MODIS LST <- 1.01*(in-situ LST) - 4.56
correctDLSTlm_a <- 1.01
correctDLSTlm_b <- 4.56


# set data path
{
  # data path
  NDVIDir         <- configList$cloudInterpolation$NDVIDir
  AlbedoDir       <- configList$cloudInterpolation$AlbedoDir
  DSRDir          <- configList$cloudInterpolation$DSRDir
  NSSRDir         <- configList$cloudInterpolation$NSSRDir
  SCSGDir         <- configList$cloudInterpolation$SCSGDir
  LSTClearSkyDir  <- configList$cloudInterpolation$LSTClearSkyDir
  LSTOriDir       <- configList$cloudInterpolation$LSTOriDir
  LSTCloudDir     <- configList$cloudInterpolation$LSTCloudDir
  
  # global data
  DEMFile       <- configList$cloudInterpolation$DEMFile
  sinAspectFile <- configList$cloudInterpolation$sinAspectFile
  cosAspectFile <- configList$cloudInterpolation$cosAspectFile
  sinSlopeFile  <- configList$cloudInterpolation$sinSlopeFile
  cosSlopeFile  <- configList$cloudInterpolation$cosSlopeFile
  MODORMYD      <- configList$cloudInterpolation$MODORMYD
  # read global data
  DEMRas       <- raster(DEMFile)
  sinAspectRas <- raster(sinAspectFile)
  cosAspectRas <- raster(cosAspectFile)
  sinSlopeRas  <- raster(sinSlopeFile)
  cosSlopeRas  <- raster(cosSlopeFile)
  
  # list data files
  NDVIFiles        <- list.files(path = NDVIDir,               pattern = "*.tif$", full.names = T)
  AlbedoFiles      <- list.files(path = AlbedoDir,             pattern = "*.tif$", full.names = T)
  DSRFiles         <- list.files(path = DSRDir,                pattern = "*.tif$", full.names = T)
  NSSRFiles        <- list.files(path = NSSRDir,               pattern = "*.tif$", full.names = T)
  SCSGFiles        <- list.files(path = SCSGDir,               pattern = "*.tif$", full.names = T)
  LSTClearSkyFiles <- list.files(path = LSTClearSkyDir,        pattern = "*.tif$", full.names = T)
  LSTOriFiles      <- list.files(path = LSTOriDir,             pattern = "*.tif$", full.names = T)
  
  # get common date of LST, DSR, NSSR, SCSG
  LSTClearSkyDates <- str_extract(LSTClearSkyFiles, '\\d{7}')
  LSTOriDates      <- str_extract(LSTOriFiles,'\\d{7}')
  DSRDates         <- str_extract(DSRFiles,'\\d{7}')
  NSSRDates        <- str_extract(NSSRFiles, '\\d{7}')
  SCSGDates        <- str_extract(SCSGFiles, '\\d{7}')
  commonDates      <- Reduce(intersect, list(LSTClearSkyDates, LSTOriDates, DSRDates, NSSRDates, SCSGDates))
  
  print(paste0("There are ",length(commonDates), " files"))
  print(commonDates)
}

#funtion
{
  # funtion: find the closest file index
  GetClosestFileIndex <-function(FileDates, commonDate){
    commonDate <- as.Date(commonDate,format = "%Y%j")
    dateDiffs <- abs(FileDates-commonDate)
    return (which.min(dateDiffs))
    
  }
  
  ItpLSTCloud_K <- function(dateStr, DEMRas, sinAspectRas, cosAspectRas, sinSlopeRas, cosSlopeRas, 
                            LSTClearSkyRas, LSTOriRas, NDVIRas, AlbedoRas, DSRRas, NSSRRas, SCSGRas){
    
    print (paste("BEGIN*************",dateStr,"*************"))
    # raster to vector
    DEM                <- as.vector(DEMRas)
    sinAspect          <- as.vector(sinAspectRas)
    cosAspect          <- as.vector(cosAspectRas)
    sinSlope           <- as.vector(sinSlopeRas)
    cosSlope           <- as.vector(cosSlopeRas)
    NDVI               <- as.vector(NDVIRas)
    Albedo             <- as.vector(AlbedoRas)
    DSR                <- as.vector(DSRRas)
    NSSR               <- as.vector(NSSRRas)
    SCSG               <- as.vector(SCSGRas)
    LSTOri             <- as.vector(LSTOriRas)
    LSTClearSky        <- as.vector(LSTClearSkyRas)
    LSTClearSkyNor     <- LSTClearSky
    # index vector
    index <- 1:length(DEM)
    
    # build variable matrix 
    M <- data.frame(index,            
                    DEM,              
                    sinAspect,         
                    cosAspect,       
                    sinSlope,          
                    cosSlope,         
                    NDVI,               
                    Albedo,       
                    DSR,             
                    NSSR, 
                    LSTClearSkyNor,         
                    LSTClearSky, 
                    LSTOri,
                    SCSG)
    M <- subset(M, (DSR > 0 & NSSR > 0))
    # extrct C D region
    M <- subset(M, ((SCSG == 3 | SCSG == 4) & !is.na(LSTClearSky)))
    
    #normalization
    M[, 2:11] <- data.Normalization(M[, 2:11], type="n4", normalization="column")
    print("Normalized.")
    M_known <<- subset(M, (!is.na(LSTOri)))
    M_itp <- subset(M, is.na(LSTOri))
    
    #correct LST values in region D which were biased
    #MODIS LST <- 1.01*(in-situ LST) - 4.56
    M_known$LSTori <- (M_known$LSTOri+correctDLSTlm_b)/correctDLSTlm_a

    
    # K-means
    M_to_cluster <- cbind(M_itp$DEM, M_itp$NDVI,  M_itp$LSTClearSkyNor)
    k_cluster <- kmeans(M_to_cluster, KClusterNum)
    cluster <- k_cluster$cluster
    M_itp <- cbind(M_itp,cluster)
    M_itp <<- M_itp
    
    centers <- k_cluster$centers
    centers <<- data.frame(centers)
    clusterIndex <- 1:KClusterNum
    clusterIndex <- data.frame(clusterIndex)
    
    
    # predict
    {
      SimilarMarsPredict <- function(clusterId){
        
        ref_group <- M_known
        center <- centers[clusterId,1:3]
        
        #calculate distance
        dist_p <- c()
        dist_p <- sqrt( ((ref_group[,2]-center[1])^2) + ((ref_group[,7]-center[2])^2) +((ref_group[,11]-center[3])^2))
        
        #select similarity pixels index:first 1000points
        ref_group <- cbind(ref_group, dist_p)
        ref_group <- ref_group[order(ref_group[,'dist_p']),]
        ref_group <- ref_group[1:1000,]
        
        #mars fit
        mars.fit <- mars(ref_group[,2:11], ref_group[,13])
        
        itps <- subset(M_itp, cluster==clusterId )
        
        r <- predict(mars.fit, itps[,2:11])
        #r <- predict(mars.fit, int_x[2:11])
        p_res <- cbind(itps$index,r )
        return(p_res)
      }
      
      
      if (nrow(M_known) < 1000){
        print ("Reference are not enough.")
        return (NA)
      }
      else if (nrow(M_itp) == 0){
        print ("No C region.")
        return (NA)
      }
      else{
        
        # parllel
        {
          print("Parallelling.")
          
          # cpu core number
          cl <- makeCluster(cpuCNum)
          clusterExport(cl,"M_itp")
          clusterExport(cl,"M_known")
          clusterExport(cl,"centers")
          clusterExport(cl,"mars")
          # apply to all missing data
          PredictValues <- parApply(cl, clusterIndex ,1, SimilarMarsPredict)
          stopCluster(cl)
          
        }
        #PredictValues <- apply(clusterIndex,1,SimilarMarsPredict)
        PredictResults <- data.frame(PredictValues[1])
        for (i in 2:KClusterNum){
          PredictResults <- rbind(PredictResults,data.frame(PredictValues[i]))
        }
        colnames(PredictResults) <- c("index", "LST_Cloud")
        
        
        print("Predicted.")
        LSTClearSkyRas[PredictResults$index] = PredictResults$LST_Cloud
        writeRaster(LSTClearSkyRas, filename = paste0(LSTCloudDir,"/",MODORMYD,"_",dateStr,"_itp_Cloud_LST.tif"), overwrite = T)
        
      }
    }
    
    print (paste("END*************",dateStr,"*************"))
    print (paste(" "))
  }
  
  
}

# batch process

for ( i in 1:length(commonDates)){
  
  commonDate <- commonDates[i]
  print(commonDate)
  
  print("Read files this day:")
  # find closest NDVI
  NDVIDates <- as.Date(str_extract(NDVIFiles,'\\d{7}'), format = "%Y%j")
  NDVIFile <- NDVIFiles[GetClosestFileIndex(NDVIDates, commonDate)]
  print(NDVIFile)
  NDVIRas <- raster(NDVIFile)
  
  # find closest Albedo
  AlbedoDates <- as.Date(str_extract(AlbedoFiles,'\\d{7}'), format = "%Y%j")
  AlbedoFile <- AlbedoFiles[GetClosestFileIndex(AlbedoDates, commonDate)]
  print(AlbedoFile)
  AlbedoRas <- raster(AlbedoFile)
  
  # get variable on the date
  LSTClearSkyFile <- LSTClearSkyFiles[grep(commonDate,LSTClearSkyFiles)]
  print(LSTClearSkyFile)
  LSTClearSkyRas <- raster(LSTClearSkyFile)
  
  LSTOriFile <- LSTOriFiles[grep(commonDate,LSTOriFiles)]
  print(LSTOriFile)
  LSTOriRas <- raster(LSTOriFile)
  
  DSRFile <- DSRFiles[grep(commonDate,DSRFiles)]
  print(DSRFile)
  DSRRas <- raster(DSRFile)
  
  NSSRFile <- NSSRFiles[grep(commonDate,NSSRFiles)]
  print(NSSRFile)
  NSSRRas <- raster(NSSRFile)
  
  SCSGFile <- SCSGFiles[grep(commonDate,SCSGFiles)]
  print(SCSGFile)
  SCSGRas <- raster(SCSGFile)
  
  if (length(NDVIRas) != length(AlbedoRas) 
      | length(NDVIRas) != length(LSTClearSkyRas) 
      | length(NDVIRas) != length(LSTOriRas) 
      | length(NDVIRas) != length(DSRRas)
      | length(NDVIRas) != length(NSSRRas)
      | length(NDVIRas) != length(SCSGRas)){
    print("!!! Raster don't match!!! SKIP!!!")
    next
  }
  print ("Raster match, Ckecked!")
  
  ItpLSTCloud_K(dateStr = commonDate,
                DEMRas = DEMRas,
                sinAspectRas = sinAspectRas,
                cosAspectRas = cosAspectRas,
                sinSlopeRas = sinSlopeRas,
                cosSlopeRas = cosSlopeRas,
                LSTClearSkyRas = LSTClearSkyRas,
                LSTOriRas = LSTOriRas,
                NDVIRas = NDVIRas,
                AlbedoRas = AlbedoRas,
                DSRRas = DSRRas,
                NSSRRas = NSSRRas,
                SCSGRas = SCSGRas)
}





