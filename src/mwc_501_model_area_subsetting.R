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


# Create combined aerial images file list --------------------------------------
aerial_prj <- CRS("+proj=tmerc +lat_0=0 +lon_0=23 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

aerial_files <- list.files(paste0(path_raster, "aerial_images/"),
                           recursive = TRUE, full.names = TRUE, 
                           pattern = glob2rx("*_RECT.tif"))


modelling_samples <- sample_xres (lowres = sen, highres_files = aerial_files, highres_prj = aerial_prj,
                                  path_rdata = path_rdata, path_temp = path_temp,
                                  path_highres_results = paste0(path_raster, "aerial_images_samples/"),
                                  n = 500000)

saveRDS(modelling_samples, file = paste0(path_rdata, "modelling_samples.rds"))
# modelling_samples <- readRDS(paste0(path_rdata, "modelling_samples.rds"))


