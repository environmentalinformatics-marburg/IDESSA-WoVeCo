# Check Sentinel-2 data and compute artifical images

# Set path ---------------------------------------------------------------------
if(Sys.info()["sysname"] == "Windows"){
  filepath_base <- "L:/GFO/BushEncroachment/"
} else {
  filepath_base <- "/media/TOSHIBA\ EXT//GFO/BushEncroachment/"
}

path_figures <- paste0(filepath_base, "figures/")
path_raster <- paste0(filepath_base, "raster/")
path_rdata <- paste0(filepath_base, "rdata/")
path_temp <- paste0(filepath_base, "temp/")
path_vector <- paste0(filepath_base, "vector/")


# Libraries --------------------------------------------------------------------
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
  stack(ds)
})
names(sen2_data) <- names(sen2_files)

for(ds in seq(length(sen2_data))){
  # nir <- ds[[which("B8" == names(ds))]]
  # red <- ds[[which("B4" == names(ds))]]
  ndvi <- spectralIndices(sen2_data[[ds]], red = "B4", nir = "B8", indices = "NDVI")
  filename <- paste0(dirname(sen2_data[[ds]][[1]]@file@name), "/ndvi.img")
  writeRaster(ndvi, filename = filename)
}

# ndvi <- lapply(seq(length(sen2_data)), function(ds ){
#   # nir <- ds[[which("B8" == names(ds))]]
#   # red <- ds[[which("B4" == names(ds))]]
#   ndvi <- spectralIndices(sen2_data[[ds]], red = "B4", nir = "B8", indices = "NDVI")
#   filename <- paste0(dirname(sen2_data[[ds]][[1]]@file@name), "/ndvi.img")
#   writeRaster(ndvi, filename = filename)
#   return(ndvi)
# })
# 
#                                                                                                                     
# 
# names(ndvi) <- names(sen2_files)


# 
# writeRaster(ndvi, filename = )

# for(i in a){
#   print(i)
# }
# 
# for(i in seq(length(a))){
#   print(names(a[i]))
# }
# 
# myreturn <- lapply(a, function(i){
#   return(i)
# })
# 
# myreturn <- lapply(seq(length(a)), function(i){
#   return(i)
# })

