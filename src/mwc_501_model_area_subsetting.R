# Prepare model dataset (Sentinel and aerial images)

# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "E:/analysis/theunis/"
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


sample_xres (lowres = sen, highres_files = aerial_files, highres_prj = aerial_prj,
             path_rdata = path_rdata, path_temp = path_temp,
             n = 100)



mask_d[1]

  class(s[[1]])

class(p[[1]])
# Get classified aerial images and extract data based on Sentinel geometry -----
aerial_files <- list.files(paste0(path_raster, "aerial_images_classified/"),
                           full.names = TRUE, pattern = glob2rx("*classified.tif"))

modelling_data <- lapply(aerial_files, function(f){
  aerial <- raster(f)
  aerial <- projectRaster(aerial, crs = projection(sen[[1]]), method="ngb")
  sen_crop <- crop(sen[[1]], aerial)
  sen_crop_sp <- rasterToPolygons(sen_crop, fun=NULL, n=4, na.rm=TRUE, 
                                  digits=12, dissolve=FALSE)
  aerial_reclass <- reclassify(aerial, c(1, 3.5, 0, 3.5, 4.5, 1, 4.5, 10, 0))
  aerial_extract <- extract(aerial_reclass, sen_crop_sp, fun = mean)
})

saveRDS(modelling_data, file = paste0(path_rdata, "modelling_data_interm.rds"))

# Prepare final modelling dataset ----------------------------------------------
# ogrListLayers(paste0(path_vector, "plots.shp"))
plots <-readOGR(dsn = paste0(path_vector, "plots.shp"),
                "plots")

# proj4string(plots)
# proj4string(sen)

#Extract raster values using training plots
modelling_data <- extract(sen, plots, df=TRUE, na.rm=TRUE)
save(modelling_data, file = paste0(path_rdata, "modelling_data_bu.RData"))
# head(modelling_data)
# str(modelling_data)
# str(plots)

# Merge extracted data to plot data alternative method
# id <- as.integer(rownames(plots@data)) + 1 # Add id
# plots$id <- id
# plots@data <- merge(plots@data, plots_exstract, by = "id", all = TRUE)

modelling_data$class<-plots$CLASS[modelling_data$ID]
head(modelling_data)

# Save extracted data to disk
save(modelling_data, file = paste0(path_model, "modelling_dataset.RData"))
