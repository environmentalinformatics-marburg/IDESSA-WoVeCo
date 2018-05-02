# Prepare combined stack of sentinel 1 and 2 data

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

source(paste0(path_source, "mwc_functions.R"))

# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)


# Create combined Sentinel 1/2 data stack --------------------------------------
sen <- stack_sen(path_raster, path_temp)


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
