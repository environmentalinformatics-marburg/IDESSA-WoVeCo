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
path_model <- paste0(filepath_base, "model/")

# Libraries --------------------------------------------------------------------
library(gpm)
library(raster)
library(ggplot2)


# Additional settings ----------------------------------------------------------
rasterOptions(tmpdir = path_temp)

# Load modelling data and create gpm object ------------------------------------
extracteddata<-get(load(paste0(path_model, "modelling_dataset.RData")))

meta <- createGPMMeta(dataset = extracteddata, 
                      selector = 1, response = c(1, 53), 
                      independent = c(2:52), meta = NULL)
lc <- gpm(extracteddata, meta)

# lc_resamples <- resamplingsByVariable(x = lc@data$input, 
#                                        selector =  as.character(unlist(lc@data$input["sample"])), grabs = 1, 
#                                        resample = 2)

lc_resamples <- lapply(seq(5), function(x){
  seq(1:nrow(lc@data$input))
})
  
lc_trte <- splitMultResp(x = lc@data$input, 
                          p = 0.75,
                          response = "class",
                          resamples = lc_resamples)

response <- "class"
independent <- lc@meta$input$INDEPENDENT
n_vars <- c(seq(length(independent)))
n_vars <- c(51) # c(10, 51)

models <- trainModel(x = lc, 
                     response = response, independent = independent,
                     resamples = lc_trte, n_var = n_vars,
                     mthd = "rf", seed_nbr = 11, cv_nbr = 5,
                     var_selection = "sd",
                     filepath_tmp = path_temp)



save(models, file = paste0(path_model, "models_rf_2016-08-25.RData"))

# models[[1]][[1]]$testing



# var_imp_scale <- compVarImp(models, scale = TRUE)
# str(var_imp_scale)
# 
# 
# var_imp_plot <- plotVarImp(var_imp)
# 
# var_imp_heat <- plotVarImpHeatmap(var_imp_scale, xlab = "Class", ylab = "Band")
# 
# tstat <- compContTests(models, mean = TRUE)
# 
# tstat_mean <- merge(tstat[[1]], prevalence, by.x = "Response", by.y = "RESPONSE")
# 
# tstat_mean[order(tstat_mean$Kappa_mean, decreasing = TRUE),]
# 
# ggplot(data = tstat_mean, aes(x = OCCURENCE, y = Kappa_mean)) + geom_point() + geom_smooth()
# 
# 
# variables <- lapply(models, function(m){
#   samples <- lapply(m, function(s){
#     if(inherits(s$model, "try-error")){
#       NULL
#     } else {
#       var <- as.data.frame(s$model$finalModel$importance)
#       var$variables <- rownames(var)
#     }
#     return(var)
#   })
#   samples <- do.call("rbind", samples)
#   return(samples)
# })
# variables <- do.call("rbind", variables)
# unique(variables$variables)
