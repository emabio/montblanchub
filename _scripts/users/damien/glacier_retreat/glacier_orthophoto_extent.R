##' ---
##' title: rasterize orthophotos studied areas and produce visual outpouts of glacial validations
##' ---


##' script init
rm(list = ls())
setwd("~/Work/MB_ATLAS_2016/_workdir/")

wd.dir <- "~/Work/MB_ATLAS_2016/_workdir/"

library(raster)
library(dplyr)
library(tidyr)

## define couple of variables
path.to.shape <- "/media/georgeda/equipes/emabio/MONT_BLANC_2016/_data/Glaciers_MontBlanc/Glacier_shapefiles/Simon_valid"
path.to.ortho <- "/media/georgeda/equipes/emabio/MONT_BLANC_2016/_data/Glaciers_MontBlanc/Images/Simon_valid"

ref.ras <- raster(file.path(wd.dir, "mod1/pred_1850.tif"))

## get the glacier shapefiles
shp.list <- list.files(path.to.shape, pattern = ".shp$", recursive = TRUE, full.names = TRUE)
ortho.list <- list.files(path.to.ortho, pattern = ".tif$", recursive = TRUE, full.names = TRUE)

shp.df <- data.frame(shape = shp.list, stringsAsFactors = FALSE) %>% 
  mutate(date = as.numeric(sub(".shp", "", sub(".*_", "", shape))),
         glacier = sub("_.*$", "", basename(shape)))

ortho.df <- data.frame(ortho = ortho.list, stringsAsFactors = FALSE) %>% 
  mutate(date = as.numeric(sub(".tif", "", sub(".*_", "", ortho))),
         glacier = sub("_.*$", "", basename(ortho)))

shp.ortho.df <- inner_join(shp.df, ortho.df)
shp.ortho.df$mask <- ""

##+

for(shp.id_ in 1:nrow(shp.ortho.df)){
  cat("\n>", shp.id_, "/", nrow(shp.ortho.df))
  ## get the glacial area
  shp.f_ <- shp.ortho.df$shape[shp.id_]
  shp_ <- shapefile(shp.f_)
  shp_ <- spTransform(shp_, crs(ref.ras))
  shp.ext_ <- extent(shp_)
  
  ## get the background area (nonglace)
  ortho.f_ <- shp.ortho.df$ortho[shp.id_]
  ortho_ <- raster(ortho.f_)
  ## find the optimal aggregatiopn factor
  agg.fact <- floor(min(res(ref.ras) / max(res(ortho_))))
  ortho.mask_ <- aggregate(ortho_, agg.fact)
  ortho.mask_[ortho.mask_ == 0] <- NA
  ortho.mask_ <- projectRaster(ortho.mask_, ref.ras)
  ortho.mask_ <- crop(ortho.mask_, shp.ext_)
  ortho.mask_[!is.na(ortho.mask_[])] <- 0
  
  ## save the rasterized orthophoto based glaciacial classif
  shp.ras <- rasterize(shp_, ortho.mask_, getCover = TRUE)
  shp.ras <- shp.ras + ortho.mask_
  shp.ras.f_ <- sub(".shp", ".grd", shp.f_)
  writeRaster(shp.ras, filename = shp.ras.f_, overwrite = TRUE)
  
  ## keep a link to the created file
  shp.ortho.df$mask[shp.id_] <- shp.ras.f_
}

##+

## do the model evaluation
bin.ref.cover.thresh <- 50 ## the pecentage of cover of the shapefile over wich we consider a ref pix to belong to a glacier
tab.scores <- NULL
for(shp.id_ in 1:nrow(shp.ortho.df)){
  cat("\n>", shp.id_, "/", nrow(shp.ortho.df))
  ras.ref_ <- raster(shp.ortho.df$mask[shp.id_])
  ## transform the ref raster (% ice cover) into binary one
  ras.ref_ <- reclassify(ras.ref_, c(-Inf, bin.ref.cover.thresh, 0, bin.ref.cover.thresh, Inf, 1))
  for(mod_ in paste0("mod", c(1,5:7))){
    ras.pred_ <- raster(file.path(wd.dir, mod_, paste0("pred_", shp.ortho.df$date[shp.id_], ".tif")))
    ras.pred_ <- crop(ras.pred_, ras.ref_)
    ras.pred_ <- mask(ras.pred_, ras.ref_)
    fit.obs.df <- na.omit(data.frame(fit = ras.pred_[], obs = ras.ref_[]))
    tss.score <- biomod2::Find.Optim.Stat(Stat = 'TSS', Fit = fit.obs.df$fit, Obs = fit.obs.df$obs)
    acc.score <- biomod2::Find.Optim.Stat(Stat = 'ACCURACY', Fit = fit.obs.df$fit, Obs = fit.obs.df$obs)
    pod.score <- biomod2::Find.Optim.Stat(Stat = 'POD', Fit = fit.obs.df$fit, Obs = fit.obs.df$obs)
    far.score <- biomod2::Find.Optim.Stat(Stat = 'FAR', Fit = fit.obs.df$fit, Obs = fit.obs.df$obs)
    colnames(tss.score) <- paste0("TSS.", colnames(tss.score))
    colnames(acc.score) <- paste0("ACCURACY.", colnames(acc.score))
    colnames(pod.score) <- paste0("POD.", colnames(pod.score))
    colnames(far.score) <- paste0("FAR.", colnames(far.score))
    tab.scores <- rbind(tab.scores, as.data.frame(cbind(shp.ortho.df[shp.id_,], data.frame(model = mod_), tss.score, acc.score, pod.score, far.score)))
  }
}

##+

## produce some graphical outputs
tab.scores

tab.scores %>% group_by(model) %>% summarize( TSS.mean = mean(TSS.best.stat), 
                                              TSS.sens.mean = mean(TSS.sensitivity),
                                              TSS.spe.mean = mean(TSS.specificity),
                                              ACCURACY.mean = mean(ACCURACY.best.stat),
                                              POD.mean = mean(POD.best.stat),
                                              FAR.mean = mean(FAR.best.stat))

tab.scores %>% group_by(glacier) %>% summarize( TSS.mean = mean(TSS.best.stat),
                                                TSS.sens.mean = mean(TSS.sensitivity),
                                                TSS.spe.mean = mean(TSS.specificity),
                                                ACCURACY.mean = mean(ACCURACY.best.stat),
                                                POD.mean = mean(POD.best.stat),
                                                FAR.mean = mean(FAR.best.stat))

tab.scores %>% group_by(glacier, date) %>% summarize( TSS.mean = mean(TSS.best.stat), 
                                                      TSS.sens.mean = mean(TSS.sensitivity),
                                                      TSS.spe.mean = mean(TSS.specificity),
                                                      ACCURACY.mean = mean(ACCURACY.best.stat),
                                                      POD.mean = mean(POD.best.stat),
                                                      FAR.mean = mean(FAR.best.stat))

##+ 

## ggplots
library(ggplot2)
gg.dat <- tab.scores %>% mutate(id = paste0(glacier, "_", date)) %>% group_by(id, model) %>% select(contains("best.stat")) %>% gather(stat, score, -c(id,model)) %>% mutate(stat = sub(".best.stat", "", as.character(stat)))
gg <- ggplot(gg.dat, aes(id, score, color = model, shape = stat)) + geom_point(size = 4) + facet_grid(~stat) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
gg


##+ 

## rasters plots
library(rasterVis)
ras.plot.dir <- file.path(wd.dir, "validation_rasters_plots")
dir.create(ras.plot.dir, recursive = TRUE, showWarnings = FALSE)

## visual checking of data
my.at <- seq(0,100, 5)
my.at.div <- seq(-100,100,5)
blue.cols <- c('#f7fbff', '#deebf7', '#c6dbef', '#9ecae1',  '#6baed6', '#4292c6',  '#2171b5', '#08519c',  '#08306b')
div.cols <- c('#67001f', '#b2182b', '#d6604d', '#f4a582', '#fddbc7', '#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061')
myThemeCont <- rasterTheme(region=colorRampPalette(blue.cols)(10))
myThemeDiv <- rasterTheme(region=colorRampPalette(div.cols)(10))

##+ 

for(shp.id_ in 1:nrow(shp.ortho.df)){
  cat("\n>", shp.id_, "/", nrow(shp.ortho.df))
  glacier_ <- shp.ortho.df$glacier[shp.id_]
  date_ <- shp.ortho.df$date[shp.id_]
  f.out_ <- file.path(ras.plot.dir, paste0("validation_", glacier_, "_", date_, ".png"))
  ras.ref_ <- raster(shp.ortho.df$mask[shp.id_])
  mods <- unique(tab.scores$model)
  stk.pred_ <- raster::stack(file.path(wd.dir, mods, paste0("pred_", date_, ".tif")))
  names(stk.pred_) <- mods
  stk.pred_ <- crop(stk.pred_, ras.ref_)
  stk.pred_ <- mask(stk.pred_, ras.ref_)
  shp_ <- shapefile(shp.ortho.df$shape[shp.id_])
  shp_ <- spTransform(shp_, crs(ras.pred_))
  nl <- nlayers(stk.pred_)
  m <- matrix(1:nl, nrow=2)
  
  ## cont models
  # png(f.out_, width = 2*480, height = 2*480)
  for (i in 1:nl){
    p <-  levelplot(subset(stk.pred_, i), par.settings = myThemeCont, margin = FALSE, main = paste0(mods[i], "\t", tab.scores$glacier[shp.id_], "\t", tab.scores$date[shp.id_]), cex = 2) + 
      latticeExtra::layer(sp.lines(shp_, lwd=2, col='orange'), data=list(shp_ = shp_))
    print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
  }
  # dev.off()
  
  ## bin models
  # png(sub(".png", "_bin.png", f.out_), width = 2*480, height = 2*480)
  for (i in 1:nl){
    p <-  levelplot(subset(stk.pred_, i) > tab.scores$TSS.cutoff[tab.scores$model == mods[i] & tab.scores$glacier == glacier_ & tab.scores$date == date_], par.settings = myThemeCont, margin = FALSE, main = paste0(mods[i], "\t", tab.scores$glacier[shp.id_], "\t", tab.scores$date[shp.id_]), cex = 2) + 
      latticeExtra::layer(sp.lines(shp_, lwd=2, col='orange'), data=list(shp_ = shp_))
    print(p, split=c(col(m)[i], row(m)[i], ncol(m), nrow(m)), more=(i<nl))
  }
  # dev.off()
  
#   levelplot(stk.pred_, par.settings = myThemeCont, margin = FALSE, main = paste0(tab.scores$glacier[shp.id_], "\t", tab.scores$date[shp.id_])) + 
#     latticeExtra::layer(sp.lines(shp_, lwd=1.5, col='black'), data=list(shp_ = shp_))
  
}

