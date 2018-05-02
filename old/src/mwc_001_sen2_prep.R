# Check Sentinel-2 data and compute artifical images

# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "F:/GFO/BushEncroachment/"
} else {
  filepath_base <- "/media/TOSHIBA\ EXT//GFO/BushEncroachment/"
}

path_figures <- paste0(filepath_base, "figures/")
path_raster <- paste0(filepath_base, "raster/")
path_rdata <- paste0(filepath_base, "rdata/")
path_temp <- paste0(filepath_base, "temp/")
path_vector <- paste0(filepath_base, "vector/")


# Libraries --------------------------------------------------------------------
library(doParallel)
library(foreach)
library(mapview)
library(raster)
library(RStoolbox)


# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)


# Get Sentinel-2 dataset names -------------------------------------------------
sen2_files <- list.files(path_raster, recursive = TRUE, full.names = TRUE,
                         pattern = glob2rx("*.img"))

sen2_files <- split(sen2_files, as.factor(dirname(sen2_files)))
names(sen2_files) <- paste0(substr(basename(names(sen2_files)), 1, 4), 
                            substr(basename(names(sen2_files)), 26, 33))
names(sen2_files)



# Compute some pixel-based artifical images ------------------------------------

sen2_data <- lapply(sen2_files, function(ds){
  return(stack(ds))
})
names(sen2_data) <- names(sen2_files)

sen2_data
# Compute NDVI
# for(ds in seq(length(sen2_data))){
#   ndvi <- spectralIndices(sen2_data[[ds]], red = "B4", nir = "B8", indices = "NDVI")
#   filename <- paste0(dirname(sen2_data[[ds]][[1]]@file@name), "/ndvi.img")
#   writeRaster(ndvi, filename = filename)
# }


#Compute PCA
#cl <- makeCluster(detectCores()-1)
#registerDoParallel(cl)

#pca <- function(ds, sen2_data, path_temp){
#  rasterOptions(tmpdir = path_temp)
#  subst <- sen2_data[[ds]][[which(names(sen2_data[[ds]])%in%paste0("B",2:8)=="TRUE")]]
#  pca_model <- rasterPCA(subst, nSamples = 0.05*ncell(subst), nCOMP=3)
#  filename <- paste0(dirname(sen2_data[[ds]][[1]]@file@name), "/pca_model.RData")
#  save(pca_model, file = filename)
#  filename <- paste0(dirname(sen2_data[[ds]][[1]]@file@name), "/pca.img")
#  writeRaster(pca_model$map, filename = filename, bylayer = TRUE)
#}

#foreach(ds = seq(length(sen2_data)), .packages = c("raster", "RStoolbox")) %dopar% pca(ds, sen2_data, path_temp)
#stopCluster(cl)

#unlink(path_temp, recursive=TRUE)
