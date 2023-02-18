# A stepwise framework for interpolating land surface temperature under cloudy
# conditions based on the solar-cloud-satellite geometry
# preprocessing: filter MODIS cloud height data 
# authors: Yuhong Chen, Zetao Cao
# last modification date: 2023-02-13
# group website: https://permalab.science/


library(configr)
library (sp)
library (spatial)
library(raster)

# read config
configFile <- "C:/Users/cr/Desktop/ISPRS_LST interpolation_cyh/Rcode/config.ini"
configList <- read.config(file = configFile)

# Input variables
## Could top height data.
input.CouldTopHeight_DataPath <- configList$CouldTopHeight$originalCouldTopHeight
## Filter window size.
input.WindowSize=as.numeric(configList$CouldTopHeight$cloudFilterWindowSize)

# Output Variables
## Path to save filtered could top height data.
output.FilteredCouldTopHeight_DataPath <- configList$CouldTopHeight$filteredCouldTopHeight


# Data files path.
DataPath <- input.CouldTopHeight_DataPath
OutPath <- output.FilteredCouldTopHeight_DataPath


# Get could top height raster data.
CTHRasters <- list.files(path=DataPath, pattern = '*.tif$', full.names = FALSE)

# Set filter window size.
w <- matrix(1, ncol = input.WindowSize, nrow = input.WindowSize)

# Check output path.
if (file.exists(OutPath)==FALSE){
  dir.create(OutPath)
}

# Filter could top height.
for (i in 1:length(CTHRasters)){
  # read as raster
  CTHRas <- raster(paste0(DataPath,'/',CTHRasters[i]))
  print (CTHRasters[i])
  focal (CTHRas, w, mean, filename = paste0(OutPath, "/Filtered_", CTHRasters[i]), 
         pad = TRUE, padValue = 0, na.rm = TRUE)
}

