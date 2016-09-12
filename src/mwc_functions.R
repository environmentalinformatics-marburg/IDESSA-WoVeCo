# Extract samples of high resolution pixels based on low resolution raster overlay
sample_xres <- function(lowres, highres_files, highres_prj = NULL,
                        path_rdata, path_temp, path_highres_results,
                        n = 10000){
  library(raster)
  library(rgdal)
  library(sp)
  
  rasterOptions(tmpdir = path_temp)
  
  
  # Get overall extent of all highres images -----------------------------------
  print("Computing extent of high resolution datasets...")
  highres_extent <- SpatialPolygonsDataFrame(
    as(extent(raster(highres_files[[1]])), "SpatialPolygons"),
    data = data.frame(NAME = highres_files[[1]]))
  
  for(f in highres_files[2:length(highres_files)]){
    highres_extent <- bind(highres_extent, SpatialPolygonsDataFrame(
      as(extent(raster(f)), "SpatialPolygons"),
      data = data.frame(NAME = f)))
  }
  
  if(!is.null(highres_prj)){
    projection(highres_extent) <- highres_prj
  }
  highres_extent <- spTransform(highres_extent, CRS(projection(lowres)))
  saveRDS(highres_extent, file = paste0(path_rdata, "highres_extent.rds"))
  # highres_extent <- readRDS(paste0(path_rdata, "highres_extent.rds"))
  
  
  # Intersect low resolution and high resolution images ------------------------
  print("Computing extent of low resolution datasets...")
  lowres_extent <- as(extent(lowres), "SpatialPolygons")
  projection(lowres_extent) <- projection(lowres)
  
  mask <- intersect(lowres_extent, highres_extent)
  mask <- disaggregate(mask)

  n_sample <- round(n / length(mask), 0)
  
  for(m in seq(length(mask))){
    print(paste0("Processing mask subset ", m, " of ", length(mask)))
    
    highres_filename <- as.character(mask[m, ]@data[1,1])
    
    act_lowres_subset <- crop(lowres[[1]], mask[m, ], cellnumbers = TRUE)
    act_lowres_subset_vec <- rasterToPolygons(act_lowres_subset, fun=NULL, n=4, na.rm=TRUE, 
                                       digits=12, dissolve=FALSE)
    names(act_lowres_subset_vec) <- highres_filename
    
    # Define samples for actual area
    set.seed(m)
    lowres_samples_nbr <- sample(length(act_lowres_subset_vec), n_sample)
    if(m == 1){
      lowres_samples <- act_lowres_subset_vec[lowres_samples_nbr,]
      names(lowres_samples) <- highres_filename
      act_lowres_subset_samples <- lowres_samples
    } else {
      act_lowres_subset_samples <- act_lowres_subset_vec[lowres_samples_nbr,]
      names(act_lowres_subset_samples) <- highres_filename
      lowres_samples <- append(lowres_samples, act_lowres_subset_samples)
    }
    saveRDS(lowres_samples, file = paste0(path_rdata, "lowres_subset_samples_", 
                                             sprintf("%08d", m), ".rds"))
    
    
    # Extract information from high resolution images for each sampled pixel
    act_highres <- raster(highres_files[grep(highres_filename, highres_files)])
    act_lowres_pixels <- spTransform(act_lowres_subset_samples, CRS(projection(act_highres)))
    
    for(p in seq(length(act_lowres_pixels))){
      act_pix <- tryCatch(crop(act_highres, act_lowres_pixels[p, ], snap = "in"),
                          error = function(e)e)
      if(!inherits(act_pix, "error")){
        filepath <- paste0(path_highres_results, 
                           "/subarea_",
                           sprintf("%08d", m),
                           "/highres_info_", 
                           sprintf("%08d", m),
                           sprintf("_%08d", p),
                           ".tif")
        dir.create(dirname(filepath), showWarnings = FALSE)
        writeRaster(act_pix, filename = filepath, overwrite = TRUE)
      }
    }

    # Extract information from low resolution images for each sampled pixel
    if(m == 1){
      lowres_pixel_data <- extract(lowres, act_lowres_subset_samples, cellnumbers = TRUE)
      lowres_pixel_data <- list(do.call("rbind", lowres_pixel_data))
      names(lowres_pixel_data) <- highres_filename
    } else {
      act_lowres_pixel_data <- extract(lowres, act_lowres_subset_samples, cellnumbers = TRUE)
      act_lowres_pixel_data <- list(do.call("rbind", act_lowres_pixel_data))
      names(act_lowres_pixel_data) <- highres_filename
      
      lowres_pixel_data <- append(lowres_pixel_data, act_lowres_pixel_data)
    }
    saveRDS(lowres_pixel_data, file = paste0(path_rdata, "lowres_pixel_data_", 
                                             sprintf("%08d", m), ".rds"))
    
  }
  
  # Combine return information -------------------------------------------------
  names(lowres_samples) <- sapply(seq(length(mask)), function(m) as.character(mask[m, ]@data[1,1]))

  highres_pixel_samples <-   lapply(seq(length(mask)), function(m){
    list.files(path_highres_results, recursive = TRUE, full.names = TRUE,
               pattern = glob2rx(sprintf("highres_info_%08d*.tif", m)))
  })
  names(highres_pixel_samples) <- sapply(seq(length(mask)), function(m) as.character(mask[m, ]@data[1,1]))
  

  info <- list(LowRes_Pixel_Data = lowres_pixel_data, HighRes_Pixel_Samples = highres_pixel_samples,
               Sample_Mask = mask, LowRes_Samples = lowres_samples)
}




#Prepare combined stack of sentinel 1 and 2 data
stack_sen <- function(path_raster, path_temp){
  library(raster)
  library(rgdal)
  
  
  # Additional settings ----------------------------------------------------------
  rasterOptions(tmpdir = path_temp)
  
  
  # Process Sentinel-2 dataset names -------------------------------------------------
  drs_S2<-dir(path = path_raster, 
              pattern = "S2A_USER_MTD_SAFL2A_PDMC.*data", full.names = TRUE)
  
  sen2<- foreach(i = drs_S2, season = c("summer", "autumn", "winter")) %do% {
    fls <- list.files(i, pattern = ".img$", full.names = TRUE)
    bands <- paste0("B", c(2:8, 11, 12))
    
    # band layers
    bands_id <- sapply(bands, function(i) grep(i, fls))
    bands_id <- unlist(bands_id)
    fls_bands <- fls[bands_id]
    rst_bands <- stack(fls_bands)
    
    # ndvi layers
    ndvi_id <- grep("ndvi", fls)
    fls_ndvi <- fls[ndvi_id]
    rst_ndvi <- stack(fls_ndvi)
    
    # pca layers
    pca_id <- grep("pca", fls)
    fls_pca <- fls[pca_id]
    rst_pca <- stack(fls_pca)
    
    # stack all layers
    rst_all <- stack(rst_bands, rst_ndvi, rst_pca)
    names(rst_all) <- paste(names(rst_all), season, sep = "_")
    
    return(rst_all)
  }
  
  sen2 <- stack(sen2)
  
  # Process Sentinel-1 dataset names ---------------------------------------------
  fls <- list.files(path = paste0(path_raster, "S1A_IW_GRDH_1SDV"), 
                    pattern = ".img$", full.names = TRUE)
  sen1 <- stack(fls)
  
  
  # Create combined Sentinel 1/2 data stack --------------------------------------
  sen <- stack(sen1, sen2)
}

