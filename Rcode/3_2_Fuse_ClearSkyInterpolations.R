# A stepwise framework for interpolating land surface temperature under cloudy
# conditions based on the solar-cloud-satellite geometry
# preprocessing: filter MODIS cloud height data 
# authors: Yuhong Chen, Zetao Cao
# last modification date: 2023-02-13
# group website: https://permalab.science/


library(sp)
library(raster)
library(stringr)
library(ggplot2)
library(configr)

# read config
configFile <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/Rcode/config.ini"
configList <- read.config(file = configFile)


LSTDir <- configList$bayesFuse$originalLST
LSTFiles <- list.files(path = LSTDir,full.names = T,pattern = "*.tif$")

eofItpDir <- configList$bayesFuse$eofInterpolationCSV

fusingResultsDir <- configList$bayesFuse$fusingResults

# function: get eof interpolations of this date
GetEOFInterpolations <- function (Dir){
  csvFiles <- list.files(path = Dir,pattern = ".csv$",full.names = T)
  # inite date frame
  df <- data.frame(NA)
  colnames(df) <- c("index")
  
  for ( i in 1:length(csvFiles)){
    csvFile <- csvFiles[i]
    refDateStr <- str_extract(csvFile,pattern = '\\d{7}_ref')
    print (refDateStr)
    refDateStr<-substr(refDateStr,1,7)
    dfRef <- read.csv(csvFile)
    dfRef <- dfRef[,c("index","predict","ref_var")]
    # unique col names
    colnames(dfRef) <- c("index",paste0("predict_",refDateStr),paste0("ref_var_",refDateStr))
    df <- merge(df,dfRef,"index",all=T)
  }
  df <- df[!is.na(df$index),]
  return (df)
}

# Bayesian
BayesianFusion <- function (df){
  ncols <- ncol(df)
  nrows <- nrow(df)
  # initiate 
  sumVar <- rep(0,nrows)
  lstVar <- rep(0,nrows)

  
  i <- 3L
  while (i <= ncols){
    df[,i] <- 1/df[,i]
    NAIndex <-which(is.na(df[,i-1]))
    # LST prediction
    df[NAIndex,i-1] <- 0L
    # reference variance
    df[NAIndex,i] <- 0L
    
    sumVar <- sumVar+ df[,i]
    lstVar <- lstVar + df[,i-1]*df[,i]
    i <- i+2
  }
  lstByes <- lstVar/sumVar
  
  return (data.frame(lstByes))
}

for (i in 1:length(LSTFiles)){
  dateStr <- str_extract(LSTFiles[i],pattern = '\\d{7}')
  print (dateStr)
  itpDir <- paste0(eofItpDir,"/",dateStr)
  df <- GetEOFInterpolations(itpDir)
  LSTByes <- BayesianFusion(df)
  colnames(LSTByes) <- c("LSTByes")
  
  rasIndex <- df$index
  LSTByes <- cbind(rasIndex, LSTByes)
  originLSTRas <- raster(LSTFiles[i])
  itpLSTRas <- originLSTRas
  itpLSTRas[LSTByes$rasIndex] <- LSTByes$LSTByes
  writeRaster(itpLSTRas, filename = paste0(fusingResultsDir,"/LSTitpClearSky_",dateStr,".tif"), overwrite = T)
}


