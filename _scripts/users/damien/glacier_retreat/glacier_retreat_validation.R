##' ---
##' title : Mont blanc glacier reatreat validation
##' author : damien g.
##' ---

##' @title: Mont blanc glacier reatreat validation
##' @description: Here we will validate our glacial retreat models
##'   with Simon's glacier classification derived from aerian images.
##'   We used the thresholds optimizing ROC scores to transform models 
##'   continuous projections into binary ones.
##' 
##' @note Anorther script that evaluate in a different way our models with Simon's photos
##'   is available in the script glacier_orthophoto_extent.R
##'   
##' @date 2016-05-11
##' @author damien g.
##' @license GPL-2

##+ , echo = FALSE
fig.path <- 'figs_grv/'
dir.create(fig.path, recursive = TRUE, showWarnings = FALSE)
knitr::opts_chunk$set(fig.width=12, fig.height=8, fig.path=fig.path,
                      echo=TRUE, warning=FALSE, message=FALSE,
                      cache  = TRUE)

##' ## script initialization
rm(list = ls())
# setwd("~/Work/MB_ATLAS_2016/_workdir/")

library(raster)

##' ## load models scores

load("~/Work/MB_ATLAS_2016/_workdir/mod.eval.RData")

##' ## produce models binary maps according to ROC threshold

mods <- paste0('mod', c(1:3, 5:7))
for (mod_ in mods){
  cat("\n>", mod_)
  ice.masks.outdir <- paste0("~/Work/MB_ATLAS_2016/_outputs/MB_glacier_masks_", mod_)
  dir.create(ice.masks.outdir, showWarnings = FALSE, recursive = TRUE)
  
  cont.ice.masks.files <- list.files(mod_, pattern = "^pred_[0-9]{4}.tif$", full.names = TRUE)
  bin.ice.masks.files <- paste0("MB_glacier_mask_", sub("pred_", "", basename(cont.ice.masks.files)))
  names(bin.ice.masks.files) <- basename(cont.ice.masks.files)
  
  for(f_ in cont.ice.masks.files){
    r_ <- raster(f_, RAT = FALSE)
    writeRaster(r_ > as.numeric(mod.eval$cutoff[mod.eval$mod == mod_ & grepl("ROC", rownames(mod.eval))]),
                filename = file.path(ice.masks.outdir, bin.ice.masks.files[basename(f_)]),
                overwrite = TRUE)
  }
}

##' ## rpoduce simple graphs to visualize model projections accuracy

for (mod_ in mods){
  cat("\n-----------------------------------------------------------------------------")
  cat("\n", mod_)
  cat("\n-----------------------------------------------------------------------------")
  
  ## define couple of variables
  path.to.shape <- "/media/georgeda/equipes/emabio/MONT_BLANC_2016/_data/Glaciers_MontBlanc/Glacier_shapefiles/Simon_valid"
  ref.proj <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  
  ## laad al needed rasters
  ice.pred.1949 <- raster(paste0("~/Work/MB_ATLAS_2016/_outputs/MB_glacier_masks_", mod_, "/MB_glacier_mask_1949.tif"))
  ice.pred.1967 <- raster(paste0("~/Work/MB_ATLAS_2016/_outputs/MB_glacier_masks_", mod_, "/MB_glacier_mask_1967.tif"))
  ice.pred.1989 <- raster(paste0("~/Work/MB_ATLAS_2016/_outputs/MB_glacier_masks_", mod_, "/MB_glacier_mask_1989.tif"))
  ice.pred.1952 <- raster(paste0("~/Work/MB_ATLAS_2016/_outputs/MB_glacier_masks_", mod_, "/MB_glacier_mask_1952.tif"))
  ice.pred.1988 <- raster(paste0("~/Work/MB_ATLAS_2016/_outputs/MB_glacier_masks_", mod_, "/MB_glacier_mask_1988.tif"))
  ice.pred.1980 <- raster(paste0("~/Work/MB_ATLAS_2016/_outputs/MB_glacier_masks_", mod_, "/MB_glacier_mask_1980.tif"))
  
  ## get the glacier shapefiles
  shp.list.files <- list.files(path.to.shape, pattern = ".shp$", recursive = TRUE, full.names = TRUE)
  shp.list <- lapply(shp.list.files, function(f_) shapefile(f_))
  names(shp.list) <- sub(".shp$", "", basename(shp.list.files))
  # ras.list <- sub(".shp$", ".grd", shp.list) ## rasterized version of the shapefiles
  
  par(mfrow = c(3,3))
  for(shp.n_ in names(shp.list)){
    r_ <- get(paste0("ice.pred.", sub("^.*_", "", shp.n_)))
    shp_ <- shp.list[[shp.n_]]
    shp_ <- spTransform(shp_, crs(r_))
    shp.ext_ <- extent(shp_)
    plot(crop(r_, shp.ext_), main = shp.n_)
    plot(shp_, add = TRUE)
  }
}
