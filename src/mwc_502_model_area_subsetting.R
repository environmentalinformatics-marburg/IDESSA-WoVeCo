# Prepare model dataset (Sentinel and aerial images)

# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "F:/analysis/theunis/"
} else {
  filepath_base <- "/media/TOSHIBA\ EXT//GFO/BushEncroachment/"
}

path_figures <- paste0(filepath_base, "figures/")
path_raster <- paste0(filepath_base, "raster/")
path_rdata <- paste0(filepath_base, "rdata/")
path_temp <- paste0(filepath_base, "temp/")
path_vector <- paste0(filepath_base, "vector/")
path_model <- paste0(filepath_base, "model/")
path_source <- paste0(filepath_base, "scripts/molopo_wc/src/")

# Libraries --------------------------------------------------------------------
library(doParallel)
library(foreach)
library(raster)
library(rgdal)
library(maptools)
library(mapview)
library(rgeos)

source(paste0(path_source, "mwc_functions.R"))

# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)


# Create combined Sentinel 1/2 data stack --------------------------------------
sen <- stack_sen(path_raster, path_temp)


# Pre-process high resolution aerial raster data -------------------------------
aerial_prj <- CRS("+proj=tmerc +lat_0=0 +lon_0=23 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Reproject aerial files to satellite projection (if necessary)
# aerial_files <- list.files(paste0(path_raster, "aerial_images/"),
#                            recursive = TRUE, full.names = TRUE, 
#                            pattern = glob2rx("*_RECT.tif"))
# 
# warp(files = aerial_files, source_prj = aerial_prj, target_prj = "EPSG:32734", 
#      outpath = paste0(path_raster, "aerial_images_utm34s/"),
#      resampling = "near")

# Read aerial files (must be same projection as low resolution datasets)
aerial_files <- list.files(paste0(path_raster, "aerial_images_utm34s/"),
                           recursive = TRUE, full.names = TRUE, 
                           pattern = glob2rx("*_RECT.tif"))

# Compute extent of all aerial files
polyg_highres <- extentRasterFiles(aerial_files)

# Compute extent of low resolution file
polyg_lowres <- raster2Polygon(sen)

# Crop lowres raster to individual rasters based on polygons and set raster
# values to pixel ID within original lowres raster
lowres_crops <- rasterCrops(sen[[1]], polyg_highres)
names(lowres_crops) <- aerial_files
lowres_crops <- lowres_crops[grep("NULL", 
                                  sapply(lowres_crops, class), invert = TRUE)]

# Sample n pixel ids from lowres data
lowres_sample_ids <- rasterSample(lowres_crops, n = 444000)
names(lowres_sample_ids) <- names(lowres_crops)

# Extract pixels from highres data based on sample locations in lowres data
aerial_files_extract_sample <- highResExtractSample(
  lowres_raster = sen[[1]], 
  sample_ids = lowres_sample_ids,
  path_highres_results = path_rdata)

# Extract pixels from lowres data based on sample locations in lowres data
lowres_sample_values <- sen[unlist(lowres_sample_ids, recursive = TRUE,
                                   use.names = FALSE)]

saveRDS(lowres_sample_values, file = paste0(path_rdata, 
                                      sprintf("lowres_sample_values.rds", s)))

