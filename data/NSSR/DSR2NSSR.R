#NSSR = DSR*(1-Albedo)

library(sp)
library(raster)

#set path
DSRDir <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/data/DSR"
AlbedoDir <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/data/Albedo"
NSSRDir <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/data/NSSR"

#read raster files
DSRfiles <- list.files(path = DSRDir, pattern = "*.tif$", full.names = F)
Albedofiles <- list.files(path = AlbedoDir, pattern = "*.tif$", full.names = F)
AlbedoDates <- as.Date(str_extract(Albedofiles, '\\d{7}'), format = "%Y%j")


#function: find cloest Albedo data
GetClosestAlbedo <- function(Ablbedodates, dateStr){
  date <- as.Date(dateStr,format = "%Y%j")
  dateDiffs <- abs(AlbedoDates-date)
  return (which.min(dateDiffs))
}


#calculate NSSR
for (i in 1:length(DSRfiles)) {
  
  #DSR 
  DSRRas <- raster(paste0(DSRDir, "/", DSRfiles[i]))
  dateStr <- str_extract(DSRfiles[i], '\\d{7}')
  
  #Albedo  
  AlbedoIndex <- GetClosestAlbedo(AlbedoDates,dateStr)
  AlbedoFile <- Albedofiles[AlbedoIndex]
  print ((paste("Albedo:", as.character(AlbedoDates[AlbedoIndex]))))
  AlbedoRas <- raster(paste0(AlbedoDir, "/", AlbedoFile))
  
  #NSSR
  NSSRRas <- DSRRas*(1-AlbedoRas) 
  
  writeRaster(NSSRRas, filename = paste0(NSSRDir, "/", dateStr, "_NSSR.tif"))
}







# date <- common_date[i]
# print(paste0("Interpolation for ", date))
# Albedo_num <- which.min(abs(as.numeric(date) - as.numeric(str_extract(Albedo_files,'\\d{7}'))))
# FiAlbedo <- Albedo_files[Albedo_num]
# RasAlbedo <- raster(paste0(AlbedoDir, '/', FiAlbedo))
