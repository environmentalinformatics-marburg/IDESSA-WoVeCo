# Subset Sentinel-1 based on pixels sampled for model preparation

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


# Sentinel 1 stack -------------------------------------------------------------
sen1_files <- list.files(paste0(path_raster, "S1A_IW_GRDH_1SDV/"),
                   recursive = TRUE, full.names = TRUE, 
                   pattern = glob2rx("*.img"))
sen1 <- stack(sen1_files)


# Extract sentinel-1 pixel id from previously sampled aerial snippets ----------
highres_extract_sample_files <- list.files(path_rdata,
                                           full.names = TRUE,
                                           pattern = glob2rx("highres_values_*.rds"))

lowres_sample_pixel_ids <- getLowResIDsfromHighResFiles(lowres_raster = sen1[[1]], 
                             highres_files = highres_extract_sample_files, 
                             cores = NULL)
  
saveRDS(lowres_sample_pixel_ids, 
        file = paste0(path_rdata, "lowres_sample_pixel_ids.rds"))

# Extract pixels from lowres data based on sample locations in lowres data
lowres_sample_values_sen1 <- sen1[unlist(lowres_sample_pixel_ids, 
                                         recursive = TRUE,
                                         use.names = FALSE)]

saveRDS(lowres_sample_values_sen1, 
        file = paste0(path_rdata, "lowres_sample_values_sen1.rds"))

