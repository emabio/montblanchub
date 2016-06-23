##' ---
##' title : Mont Blanc glacier retreat modelling
##' author : damien g.
##' ---


##' ## description 
##' Here we will test several modelling procedure to model MB glacier
##' retreat from LGM (~1850-1880), 1952, 1970-1975, 2003-2008. 
##' 
##' The 2 main modelling ideas tha will be developped here are:
##' 
##'   - linear mixed model with random effect due to pixel => Carlson et al.
##'   - linear mixed model with random effect due to glacier
##'   - (non-linear mixed models ?)
##'   

##+ , echo = FALSE
fig.path <- 'figs_grm/'
dir.create(fig.path, recursive = TRUE, showWarnings = FALSE)
knitr::opts_chunk$set(fig.width=12, fig.height=8, fig.path=fig.path,
                      echo=TRUE, warning=FALSE, message=FALSE,
                      cache  = TRUE)

##' ## script inititialisation -------------------------------------------------

##+
rm(list = ls())
# setwd("~/Work/MB_ATLAS_2016/_workdir")

## load a bench of packages
library(raster)
library(foreach)
library(parallel)
library(doParallel)
library(dplyr)
library(tidyr)
library(lme4)
library(animation)
library(rasterVis)
library(ggplot2)
library(MuMIn)

## refer the path to the input and output directories on the network
path.to.dat <- "/media/georgeda/equipes/emabio/MONT_BLANC_2016/_data"
path.to.res <- "/media/georgeda/equipes/emabio/MONT_BLANC_2016/_results"
path.to.fig <- "/media/georgeda/equipes/emabio/MONT_BLANC_2016/_figures"


##' ## load input data ---------------------------------------------------------

## the reference 100m grid of the MB
mb.grid.ref <- raster(file.path(path.to.dat, "MB_Atlas_2016_masks/MB_Atlas_2016_mask_100m_ETRS89.img"))
mb.grid.shp <- shapefile(file.path(path.to.dat, "MB_Atlas_2016_masks/MB_Atlas_2016_mask_ETRS89.shp"))

## the dem of the area
mb.dem <- raster(file.path(path.to.dat, "DEM25_Atlas_2016_ETRS89.tif"))

## calculate the slope and the aspect of the area
mb.slope <- terrain(mb.dem, opt = 'slope', filename = "SLOPE25_Atlas_2016_ETRS89.tif", overwrite = TRUE)
mb.aspect <- terrain(mb.dem, opt = 'aspect', filename = "ASPECT25_Atlas_2016_ETRS89.tif", overwrite = TRUE)

## resample the dem, slope, and aspect on our ref grid
mb.dem <- resample(mb.dem, mb.grid.ref, filename = "DEM100_Atlas_2016_ETRS89.tif", overwrite = TRUE)
mb.slope <- resample(mb.slope, mb.grid.ref, filename = "SLOPE100_Atlas_2016_ETRS89.tif", overwrite = TRUE)
mb.aspect <- resample(mb.aspect, mb.grid.ref, filename = "ASPECT100_Atlas_2016_ETRS89.tif", overwrite = TRUE)

## the dem of the area
mb.dem <- raster(file.path(path.to.dat, "DEM25_Atlas_2016_ETRS89.tif"))
mb.dem <- resample(mb.dem, mb.grid.ref, filename = "DEM100_Atlas_2016_ETRS89.tif", overwrite = TRUE)

## the glacier shapefiles
mb.glac.shp.PAG <- shapefile(file.path(path.to.dat, "Glaciers_MontBlanc/Glacier_shapefiles/MAJ_JR/PAGMergeETRS.shp"))
mb.glac.shp.1952 <- shapefile(file.path(path.to.dat, "Glaciers_MontBlanc/Glacier_shapefiles/MAJ_JR/20160412_1952MergeETRS_nunatak.shp"))
mb.glac.shp.1975 <- shapefile(file.path(path.to.dat, "Glaciers_MontBlanc/Glacier_shapefiles/MAJ_JR/_1975MergeETRS.shp"))
mb.glac.shp.2008 <- shapefile(file.path(path.to.dat, "Glaciers_MontBlanc/Glacier_shapefiles/MAJ_JR/_2008MergeETRS.shp"))
## remove all the newly appeared glacier
mb.glac.shp.2008 <- mb.glac.shp.2008[mb.glac.shp.2008@data$OBJECTID_1 <  150, ]

## we will use a dynamic mask to remove dynamic in pixs where we are sure that pixs were never unglaced
## we will use the altitude 3100m as the current level and an elevation rate of 170m over 1984 - 2010 for
## the ELA (equilibrium-line altitude) ( ref: Rabatel et al. 2010 )
ela.inc.rate <- 170/(2010-1984)
ice.lim.2008 <- 3100
ice.lim.df <- data.frame(year = 1850:2100) %>% 
  mutate(ele = ice.lim.2008 + (year - 2008) * ela.inc.rate)
ice.lim.df$ele[ice.lim.df$year < 2008] <- ice.lim.2008

# ## the past ice mask is defined at 3100m from PAG (1850) to 2008
# path.to.ice.mask <- "mb_ice_masks"
# dir.create(path.to.ice.mask, showWarnings = FALSE, recursive = TRUE)
# mb.dem.on.grid.ref <- mb.dem  * mb.grid.ref
# mb.past.ice.mask <- reclassify(mb.dem.on.grid.ref, c(-Inf, ice.lim.2008, 0, ice.lim.2008, Inf, 1), 
#                                filename = file.path(path.to.ice.mask, "ice_mask_past.tif"),
#                                overwrite = TRUE)
# ## the future ice masks
# for(yr_ in 2008:2100){
#   ela.inc_ <- ela.inc.rate * (yr_ - 2008)
#   ice.lim_ <- ice.lim.2008 + ela.inc_
#   mb.ice.mask_ <- reclassify(mb.dem.on.grid.ref, c(-Inf, ice.lim_, 0, ice.lim_, Inf, 1), 
#                              filename = file.path(path.to.ice.mask, paste0("ice_mask_", yr_, ".tif")),
#                              overwrite = TRUE)
# }

## keep only the "mer de glace area"
# plot(mb.grid.ref)
# plot(mb.glac.shp.PAG, add = TRUE)
# ext.mdg <- drawExtent()
# save(ext.mdg, file = "ext.mdg.RData")
# load("ext.mdg.RData")
# mb.grid.ref <- crop(mb.grid.ref, ext.mdg)
# mb.dem <- crop(mb.dem, ext.mdg)
# mb.slope <- crop(mb.slope, ext.mdg)
# mb.aspect <- crop(mb.aspect, ext.mdg)

## define the period names and the associated date
period.names <- c("PAG", "1952", "1975", "2008")
period.year <- c(1850, 1952, 1975, 2008)

##' ## rasterize the shapefiles to get percentage of cover by cell -------------

cover.dir.out <- "cover_raster"
dir.create(cover.dir.out, showWarnings = FALSE, recursive = TRUE)
mb.glac.cover.list <- mclapply(period.names, function(d_){
  filename_ <- file.path(cover.dir.out, paste0("mb_glacier_cover_", d_,".img"))
  mb.glac.cover_ <- rasterize(get(paste0("mb.glac.shp.", d_)), mb.grid.ref, getCover = TRUE, 
                              filename = filename_,
                              overwrite = TRUE)
  writeRaster(mb.glac.cover_, filename = filename_, overwrite = TRUE)
  return(filename_)
}, mc.cores = 4)

##' **note** because of a bug in rasterize function we are not able to produce the accurate version
##' of 1975 raster file so we will use the ones manually done by Julien.
file.copy(from = file.path(path.to.dat, "Glaciers_MontBlanc/Glacier_shapefiles/MAJ_JR/_1975prct.img"),
          to = mb.glac.cover.list[[which(period.year == 1975)]], 
          overwrite = TRUE)

mb.glac.cover.stk <- raster::stack(unlist(mb.glac.cover.list))

## visual checking of data
my.at <- seq(0,100, 5)
my.at.div <- seq(-100,100,5)
blue.cols <- c('#f7fbff', '#deebf7', '#c6dbef', '#9ecae1',  '#6baed6', '#4292c6',  '#2171b5', '#08519c',  '#08306b')
div.cols <- c('#67001f', '#b2182b', '#d6604d', '#f4a582', '#fddbc7', '#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061')
myThemeCont <- rasterTheme(region=colorRampPalette(blue.cols)(10))
myThemeDiv <- rasterTheme(region=colorRampPalette(div.cols)(10))

mb.glac.cover.change.PAG <- subset(mb.glac.cover.stk, 2:4) - subset(mb.glac.cover.stk, 1)
names(mb.glac.cover.change.PAG) <- paste0("from_PAG_to_", period.names[2:4])

mb.glac.cover.change.1952 <- subset(mb.glac.cover.stk, 3:4) - subset(mb.glac.cover.stk, 2)
names(mb.glac.cover.change.1952) <- paste0("from_1952_to_", period.names[3:4])

mb.glac.cover.change.1975 <- subset(mb.glac.cover.stk, 4) - subset(mb.glac.cover.stk, 3)
names(mb.glac.cover.change.1975) <- paste0("from_1975_to_", period.names[4])

##+ glacier_cover_input
levelplot(mb.glac.cover.stk, at = my.at, par.settings = myThemeCont)
##+ glacier_cover_change_PAG_input
levelplot(mb.glac.cover.change.PAG, at = my.at.div, par.settings = myThemeDiv)
##+ glacier_cover_change_1952_input
levelplot(mb.glac.cover.change.1952, at = my.at.div, par.settings = myThemeDiv)
##+ glacier_cover_change_1975_input
levelplot(mb.glac.cover.change.1975, at = my.at.div, par.settings = myThemeDiv)

##' rasterize the shapefiles to get the glacial mass id ------------------------

glacID.dir.out <- "glacID_raster"
dir.create(glacID.dir.out, showWarnings = FALSE, recursive = TRUE)
mb.glac.id.list <- mclapply(period.names, function(d_){
  filename_ <- file.path(glacID.dir.out, paste0("mb_glacier_id_", d_,".img"))
  mb.glac.cover_ <- rasterize(get(paste0("mb.glac.shp.", d_)), mb.grid.ref, field = "OBJECTID_1", 
                              filename = filename_,
                              overwrite = TRUE)
  return(filename_)
}, mc.cores = 4)

##' **note** because of a bug in rasterize function we are not able to produce the accurate version
##' of 1975 raster file so we will use the ones manually done by Julien.
file.copy(from = file.path(path.to.dat, "Glaciers_MontBlanc/Glacier_shapefiles/MAJ_JR/_1975MergeETRS.img"),
          to = mb.glac.id.list[[which(period.year == 1975)]], 
          overwrite = TRUE)

mb.glac.id.stk <- stack(unlist(mb.glac.id.list))

##+ glacier_ids_input
qual.cols <- c('#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f')
levelplot(mb.glac.id.stk, col.regions=qual.cols)

##' ## construct the distance to the ice hedge variable ------------------------

dist.to.border <- function(r_, thresh_ = 10){
  r.inv <- r_
  r.inv[] <- 1
  r.inv[r_[] > thresh_] <- NA
  r.dist <- distance(r.inv)
  return(r.dist)
}

mb.glac.dist.stk <- stack(lapply(1:nlayers(mb.glac.cover.stk), function(l_) dist.to.border(subset(mb.glac.cover.stk, l_))))
names(mb.glac.dist.stk) <- sub('cover', 'dist', names(mb.glac.cover.stk))
levelplot(mb.glac.dist.stk)

##' ## construct the dataset for model analysis --------------------------------

mb.grid.ref.ids <- data.frame(cellID = 1:ncell(mb.glac.id.stk[[1]]),
                              alti = mb.dem[],
                              slope = mb.slope[],
                              aspect = mb.aspect[],
                              glacialID = as.factor(mb.glac.id.stk[[1]][])) %>% na.omit

mb.glac.dat.list <- mclapply(1:length(period.names), function(d_){
  df_ <- bind_cols( mb.grid.ref.ids,
                    data.frame(
                      year = period.year[d_],
                      year_rel = period.year[d_] - min(period.year),
                      cover = mb.glac.cover.stk[[d_]][mb.grid.ref.ids$cellID],
                      dist.to.border = mb.glac.dist.stk[[d_]][mb.grid.ref.ids$cellID]))
}, mc.cores = 4 ) 
mb.glac.dat.df <- bind_rows(mb.glac.dat.list)

##' We consider that theire is no change in the pixels upper than 3100m since LGM.. we will
##' take the 2008 classification as a reference or we concider 100% of ice everywhere

mb.glac.dat.df_ <- mb.glac.dat.df %>% 
  group_by(cellID) %>%
  filter(alti > 3100, year == 2008) %>% 
  select(cellID, cover, alti) 

##' summarize the state of pixs over 3100m in 2008

summary(mb.glac.dat.df_)
ggplot(mb.glac.dat.df_, aes(cover)) + geom_histogram()
ggplot(mb.glac.dat.df_, aes(x = alti, y = cover)) + geom_point()
ggplot(mb.glac.dat.df_ %>% 
  mutate(alti_class = cut(alti, breaks = c(3000, 3500, 4000, 4500, 5000),
                          labels = c('3100-3500 m', '3500-4000 m', '4000-4500 m', 'over 4500 m'))),
  aes(alti_class, cover)) + geom_boxplot(varwidth = TRUE)

## spatially represent the pc cover in 2008 over 3100m
mb.glac.cover.over.3100 <- subset(mb.glac.cover.stk, "mb_glacier_cover_2008")
mb.glac.cover.over.3100[mb.dem[] < 3100] <- NA 
plot(mb.glac.cover.over.3100); plot(mb.grid.shp, add = TRUE)

mb.glac.dat.df__ <- mb.glac.dat.df %>% filter(alti > 3100) %>% select(-cover) %>% left_join(mb.glac.dat.df_)

mb.glac.dat.df <-  mb.glac.dat.df %>% filter(alti <= 3100) %>% bind_rows(mb.glac.dat.df__)



##' ### number of glacial pixels along time

gg.dat <- mb.glac.dat.df %>% group_by(year) %>% summarize(nb.pix.glac.100 = sum(cover >= 100),
                                                          nb.pix.glac.75 = sum(cover >= 75),
                                                          nb.pix.glac.50 = sum(cover >= 50),
                                                          nb.pix.glac.25 = sum(cover >= 25)) %>%
  gather("glac.thresh.str", "nb.pix.glac", starts_with("nb.pix.glac"))

##+ nb_glacial_pix_fct_year
gg <- ggplot(gg.dat, aes(year, nb.pix.glac, fill = glac.thresh.str)) +
  geom_bar(stat = 'identity', position = "dodge")
gg

gg <- ggplot(gg.dat, aes(glac.thresh.str, nb.pix.glac, fill = as.factor(year))) +
  geom_bar(stat = 'identity', position = "dodge")
gg

##' ### trend in number of glacial pixels by glacier

gg.dat <- mb.glac.dat.df %>% group_by(year, glacialID) %>% summarize(nb.pix.glac.25 = sum(cover >= 25)) %>% ungroup
gg.dat.ref <- gg.dat %>% filter(year == 1850) %>% mutate(nb.pix.glac.ref = nb.pix.glac.25) %>% select(glacialID, nb.pix.glac.ref)
gg.dat <- gg.dat %>% left_join(gg.dat.ref) %>% mutate(pc.glace.change = 100 * (nb.pix.glac.25 - nb.pix.glac.ref) / nb.pix.glac.ref)

##+ trend__glacial_pix_by_glacier
gg <- ggplot(gg.dat , aes(x = year, y = pc.glace.change, group = glacialID)) + geom_smooth(method = 'lm', se = FALSE)
gg


##' ### create prediction dataset ----------------------------------------------
mb.glac.pred.dat.list <- mclapply(unique(sort(c(seq(1850, 2100, 5), seq(1990, 2020,1), c(1949, 1952, 1967, 1980, 1988, 1989)))), function(yr_){
  data.frame(mb.grid.ref.ids, year = yr_, year_rel = yr_ - min(period.year))
}, mc.cores = 8)
mb.glac.pred.dat.df <- bind_rows(mb.glac.pred.dat.list)

##' ### build a function to reconstruct raster stack from prediction -----------
build_stack_from_pred <- function(pred_, df_ , ras_ = mb.grid.ref, output.dir_ = "mod1", mc.cores_ = 8){
  dir.create(output.dir_, showWarnings = FALSE, recursive = TRUE)
  df__ <- cbind(df_, pred = pred_)
  years_ <- unique(df_$year)
  pred.ras.list <- mclapply(years_, function(yr_){
    df___ <- df__ %>% filter(year == yr_)
    ras__ <- ras_
    ras__[df___$cellID] <- df___$pred
    writeRaster(ras__, filename = file.path(output.dir_, paste0("pred_",yr_ , ".tif")), overwrite = TRUE)
  }, mc.cores = mc.cores_)
}

##' ## produce models ----------------------------------------------------------

##' ### define a function to keep 2008 glacial cover values for pixs over ela

# pred.df <- pred.fut.mod1
# obs.df <- mb.glac.dat.df
update.pred.cover <- function(pred.df, obs.df, ice.lim.df, stat = 'cover'){
  new.pred.df_ <- data.frame()
  ## keep only the 2008 cover values as ref
  obs.df <- obs.df %>% filter(year == 2008)
  for(yr_ in unique(pred.df$year)){ ## loop over years
    cat("\n yr:", yr_)
    ## get the value of the ela for the year of interest
    ela_ <- ice.lim.df$ele[ice.lim.df$year == yr_]
    ## get unchanged pixels (below ela)
    pred.no.change_ <- pred.df %>% filter(year == yr_, alti <= ela_)
    ## put the 2008 cover value on pixs over ela
    if(stat == 'cover'){
      pred.updated_ <- pred.df %>% filter(year == yr_, alti > ela_) %>% 
        left_join(select(obs.df, cellID, cover)) %>% 
        mutate(pred.cover = cover) %>% select(- cover)
    } else if(stat == 'dist.to.border'){ 
      ## a bit contrintuitive in this case because the distance to the border should change beyond ela
      pred.updated_ <- pred.df %>% filter(year == yr_, alti > ela_) %>% 
        left_join(select(obs.df, cellID, dist.to.border)) %>% 
        mutate(pred.cover = dist.to.border) %>% select(- dist.to.border)
    }
    ## merge yearly updated predictions
    new.pred_ <- bind_rows(pred.no.change_, pred.updated_) %>% arrange(cellID)
    ## add predictions to the output data.frame
    new.pred.df_ <- bind_rows(new.pred.df_, new.pred_)
  }
  return(new.pred.df_)
}

##' ### mod1: simple linear mixed models with random effect on cellID

# mb.glac.dat.df <- mb.glac.dat.df %>% ungroup %>% na.omit %>% mutate(year_rel = as.numeric(year_rel), cellID.fact = as.factor(cellID))

##+ , eval = TRUE
mod1 <- lmer(cover ~ year_rel * alti + slope + (1+year_rel|cellID), 
             data = mb.glac.dat.df)
save(mod1, file = "mod1.RData")

##+ , echo = FALSE
load("mod1.RData")

##+
mod1
r.squaredGLMM(mod1)

## do data.frame prediction
pred.fut.mod1 <- cbind(mb.glac.pred.dat.df, pred.cover = predict(mod1, mb.glac.pred.dat.df))
## keep 2008 values for points over ELA
pred.fut.mod1 <- update.pred.cover(pred.fut.mod1, obs.df = mb.glac.dat.df, ice.lim.df = ice.lim.df)
## remove unrelalitc values
pred.fut.mod1$pred.cover[pred.fut.mod1$pred.cover < 0]  <- 0 
pred.fut.mod1$pred.cover[pred.fut.mod1$pred.cover > 100]  <- 100


pred.dat.mod1 <- cbind(mb.glac.dat.df, pred.cover = fitted(mod1))
## keep 2008 values for points over ELA
pred.dat.mod1 <- update.pred.cover(pred.dat.mod1, obs.df = mb.glac.dat.df, ice.lim.df = ice.lim.df)
## remove unrelalitc values
pred.dat.mod1$pred.cover[pred.dat.mod1$pred.cover < 0]  <- 0 
pred.dat.mod1$pred.cover[pred.dat.mod1$pred.cover > 100]  <- 100

pred.dat.mod1.list <- build_stack_from_pred(pred_ = pred.dat.mod1$pred.cover, df_ = mb.glac.dat.df, output.dir_ = "mod1")
pred.dat.mod1.stk <- stack(pred.dat.mod1.list)

##+ mod1_pred_calib_dates
levelplot(pred.dat.mod1.stk, at = my.at, par.settings = myThemeCont)

## convert into rasters
pred.fut.mod1.list <- build_stack_from_pred(pred_ = pred.fut.mod1$pred.cover, df_ = pred.fut.mod1, output.dir_ = "mod1")
pred.fut.mod1.stk <- stack(pred.fut.mod1.list)

## produce a gif
years_ <- as.numeric(unique(sub("^.*_", "", names(pred.fut.mod1.stk))))

##+ , eval = TRUE
saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod1.list[[t_]], main = list(years_[t_]), 
                    at = my.at, par.settings = myThemeCont))
  }
}, movie.name = "mod1_1850_2100_MB.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)


##' ![mod1.gif](mod1_1850_2100_MB.gif)
##' 


##' ### mod2: simple linear mixed models with random effect on glacial ID

# mod2.1 <- lmer(cover ~ year_rel * alti + slope + (1+year_rel+alti+slope|glacialID), 
#                data = mb.glac.dat.df)
# mod2.2 <- lmer(cover ~ year_rel * alti + slope + (1+year_rel+slope|glacialID), 
#                data = mb.glac.dat.df)
# mod2.3 <- lmer(cover ~ year_rel * alti + slope + (1+year_rel|glacialID), 
#                data = mb.glac.dat.df)
# 
# r.squaredGLMM(mod2.1)
# r.squaredGLMM(mod2.2)
# r.squaredGLMM(mod2.3)

##+ , eval = TRUE
do.drege = FALSE
if(do.drege){
  clusterType <- "PSOCK"
  clust <- try(makeCluster(getOption("cl.cores", 8), type = clusterType))
  clusterExport(clust, "mb.glac.dat.df")
  invisible(clusterCall(clust, "library", "lme4", character.only = TRUE))
  # mod2 <- pdredge(lmer(cover ~ year_rel + alti + slope + year_rel:alti + year_rel:slope + 
  #                        (year_rel|glacialID) + (slope|glacialID) + (alti|glacialID) + 
  #                        (year_rel:alti|glacialID) + (year_rel:slope|glacialID), 
  #                      data = na.omit(mb.glac.dat.df)), cluster = clust, fixed = c("year_rel", "alti", "slope"))
  mod2 <- pdredge(lmer(cover ~ year_rel + alti + slope + year_rel:alti + year_rel:slope + 
                         (year_rel|glacialID) + (slope|glacialID) + (alti|glacialID), 
                       data = mb.glac.dat.df), cluster = clust, fixed = c("year_rel", "alti", "slope"))
  stopCluster(clust)
} else {
  mod2 <- lmer(cover ~ year_rel + alti + slope + year_rel:alti + year_rel:slope + 
                         (year_rel|glacialID) + (slope|glacialID) + (alti|glacialID), 
                       data = mb.glac.dat.df)
}

save(mod2, file = "mod2.RData")

##+ , echo = FALSE
load("mod2.RData")

##+ , 
mod2

##+ , eval = TRUE
if(do.drege){
  mod2.best <- get.models(mod2, 1)[[1]]
} else {
  mod2.best <- mod2
}

save(mod2.best, file = "mod2.best.RData")

##+ , echo = FALSE
load("mod2.best.RData")

##+ 
r.squaredGLMM(mod2.best)

## do data.frame prediction
pred.fut.mod2 <- cbind(mb.glac.pred.dat.df, pred.cover = predict(mod2.best, mb.glac.pred.dat.df))
## keep 2008 values for points over ELA
pred.fut.mod2 <- update.pred.cover(pred.fut.mod2, obs.df = mb.glac.dat.df, ice.lim.df = ice.lim.df)
## remove unrelalitc values
pred.fut.mod2$pred.cover[pred.fut.mod2$pred.cover < 0]  <- 0 
pred.fut.mod2$pred.cover[pred.fut.mod2$pred.cover > 100]  <- 100

## convert into rasters
pred.fut.mod2.list <- build_stack_from_pred(pred_ = pred.fut.mod2$pred.cover, df_ = pred.fut.mod2, output.dir_ = "mod2")
pred.fut.mod2.stk <- stack(pred.fut.mod2.list)

## produce a gif
years_ <- as.numeric(unique(sub("^.*_", "", names(pred.fut.mod2.stk))))

##+ , eval = TRUE
saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod2.list[[t_]], main = list(years_[t_]), 
                    at = my.at, par.settings = myThemeCont))
  }
}, movie.name = "mod2_1850_2100_MB.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)

##' ![mod2.gif](mod2_1850_2100_MB.gif)

##+ , eval = TRUE
saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod1.list[[t_]] > 25, main = list(years_[t_]), 
                    at = seq(0,1, by = 0.5), par.settings = myThemeCont))
  }
}, movie.name = "mod1_1850_2100_MB_bin.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)


#' ### mod3: quadratic linear mixed models with random effect on glacial ID

# mod3.1 <- lmer(cover ~ year_rel + alti + slope + I(year_rel^2) + I(alti^2) + I(slope^2) + (1+year_rel+alti+slope+I(year_rel^2)|glacialID), 
#                data = mb.glac.dat.df)
# mod3.2 <- lmer(cover ~ year_rel + alti + slope + I(year_rel^2) + I(alti^2) + I(slope^2) + (1+year_rel+alti+slope+I(year_rel^2)+I(slope^2)||glacialID), 
#                data = mb.glac.dat.df)
# # mod3.3 <- lmer(cover ~ year_rel + alti + slope + I(year_rel^2) + I(alti^2) + I(slope^2) + (1+year_rel+alti+slope+I(year_rel^2)+I(slope^2)+I(alti^2)||glacialID), 
# #                data = mb.glac.dat.df) ## failing

##+ , eval = TRUE
if(do.drege){
  # Set up the cluster
  clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
  clust <- try(makeCluster(getOption("cl.cores", 16), type = clusterType))
  clusterExport(clust, "mb.glac.dat.df")
  invisible(clusterCall(clust, "library", "lme4", character.only = TRUE))
  # mod3 <- pdredge(lmer(cover ~ year_rel + alti + slope + I(year_rel^2) + I(alti^2) + I(slope^2) + year_rel:alti + year_rel:slope + 
  #                (year_rel|glacialID) + (alti|glacialID) + (slope|glacialID) + (I(year_rel^2)|glacialID) + (I(alti^2)|glacialID) + 
  #                (I(slope^2)|glacialID), 
  #              data = na.omit(mb.glac.dat.df)), cluster = clust, fixed = c("year_rel", "alti", "slope"))
  # mod3 <- pdredge(lmer(cover ~ year_rel + alti + slope + I(year_rel^2) + I(alti^2) + I(slope^2) + year_rel:alti + year_rel:slope + 
  #                        (year_rel|glacialID) + (alti|glacialID) + (slope|glacialID) + (I(year_rel^2)|glacialID) + 
  #                        (I(slope^2)|glacialID) , 
  #                      data = na.omit(mb.glac.dat.df)), cluster = clust, fixed = c("year_rel", "alti", "slope"))
  mod3 <- pdredge(lmer(cover ~ year_rel + alti + slope + I(alti^2) + I(slope^2) + year_rel:alti + year_rel:slope + 
                         (year_rel|glacialID) + (alti|glacialID) + (slope|glacialID), 
                       data = mb.glac.dat.df), cluster = clust, fixed = c("year_rel", "alti", "slope"))
  stopCluster(clust)
} else {
  mod3 <- lmer(cover ~ year_rel + alti + slope + I(alti^2) + I(slope^2) + year_rel:alti + year_rel:slope + 
         (year_rel|glacialID) + (alti|glacialID) + (slope|glacialID), 
       data = mb.glac.dat.df)
}
save(mod3, file = "mod3.RData")

##+ , echo = FALSE
load("mod3.RData")

##+ 
mod3

##+ , eval = TRUE
if(do.drege){
  mod3.best <- get.models(mod3, 1)[[1]]
} else {
  mod3.best <- mod3
}

save(mod3.best, file = "mod3.best.RData")

##+ , echo = FALSE
load("mod3.best.RData")

##+ 
r.squaredGLMM(mod3.best)



## do data.frame prediction
pred.fut.mod3 <- cbind(mb.glac.pred.dat.df, pred.cover = predict(mod3.best, mb.glac.pred.dat.df))
## keep 2008 values for points over ELA
pred.fut.mod3 <- update.pred.cover(pred.fut.mod3, obs.df = mb.glac.dat.df, ice.lim.df = ice.lim.df)
## remove unrelalitc values
pred.fut.mod3$pred.cover[pred.fut.mod3$pred.cover < 0]  <- 0 
pred.fut.mod3$pred.cover[pred.fut.mod3$pred.cover > 100]  <- 100

## convert into rasters
pred.fut.mod3.list <- build_stack_from_pred(pred_ = pred.fut.mod3$pred.cover, df_ = pred.fut.mod3, output.dir_ = "mod3")
pred.fut.mod3.stk <- stack(pred.fut.mod3.list)

## produce a gif
years_ <- as.numeric(sub("^.*_", "", names(pred.fut.mod3.stk)))

##+ , eval = TRUE
saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod3.list[[t_]], main = list(years_[t_]), 
                    at = my.at, par.settings = myThemeCont))
  }
}, movie.name = "mod3_1850_2100_MB.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)

##' ![mod3.gif](mod3_1850_2100_MB.gif)

# ##' ### mod4: non-linear mixed models with random effect on glacial ID
# 
# # startvec <- c(Asym.y = 100, xmid.y = 150, scal.y = 100)
# # (nm1 <- nlmer(cover ~ SSlogis(year_rel, Asym.y, xmid.y, scal.y) ~ xmid.y|glacialID,
# #               mb.glac.dat.df, start = startvec))
# # 
# # lmer(cover ~ year_rel + alti + slope + I(alti^2) + I(slope^2) + year_rel:alti + year_rel:slope + 
# #        (year_rel|glacialID) + (alti|glacialID) + (slope|glacialID), 
# #      data = na.omit(mb.glac.dat.df)), cluster = clust, fixed = c("year_rel", "alti", "slope")

##' ### mod5: linear mixed models with hierarchical random effect

mod5 <- lmer(cover ~ year_rel + alti + slope + I(alti^2) + I(slope^2) + year_rel:alti + year_rel:slope + 
               + (1 + year_rel|cellID/glacialID) + (alti + slope|glacialID), 
             data = mb.glac.dat.df)
mod5
r.squaredGLMM(mod5)
save(mod5, file = 'mod5.RData')

mod5.best <- mod5

## do data.frame prediction
pred.fut.mod5 <- cbind(mb.glac.pred.dat.df, pred.cover = predict(mod5.best, mb.glac.pred.dat.df))
## keep 2008 values for points over ELA
pred.fut.mod5 <- update.pred.cover(pred.fut.mod5, obs.df = mb.glac.dat.df, ice.lim.df = ice.lim.df)
## remove unrelalitc values
pred.fut.mod5$pred.cover[pred.fut.mod5$pred.cover < 0]  <- 0 
pred.fut.mod5$pred.cover[pred.fut.mod5$pred.cover > 100]  <- 100

## convert into rasters
pred.fut.mod5.list <- build_stack_from_pred(pred_ = pred.fut.mod5$pred.cover, df_ = pred.fut.mod5, output.dir_ = "mod5")
pred.fut.mod5.stk <- stack(pred.fut.mod5.list)

## produce a gif
years_ <- as.numeric(sub("^.*_", "", names(pred.fut.mod5.stk)))

##+ , eval = TRUE
saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod5.list[[t_]], main = list(years_[t_]), 
                    at = my.at, par.settings = myThemeCont))
  }
}, movie.name = "mod5_1850_2100_MB.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)

saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod5.list[[t_]] > 25, main = list(years_[t_]), 
                    at = seq(0,1, by = 0.5), par.settings = myThemeCont))
  }
}, movie.name = "mod5_1850_2100_MB_bin.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)

##' ## model with the distance to the border as explanatory

mod6 <- lmer(dist.to.border ~ year_rel * alti + slope  + (1+year_rel|cellID),
             data = mb.glac.dat.df)

mod6
r.squaredGLMM(mod6)
save(mod6, file = 'mod6.RData')

mod6.best <- mod6

## do data.frame prediction
pred.fut.mod6 <- cbind(mb.glac.pred.dat.df, pred.cover = predict(mod6.best, mb.glac.pred.dat.df))
## keep 2008 values for points over ELA
pred.fut.mod6 <- update.pred.cover(pred.fut.mod6, obs.df = mb.glac.dat.df, ice.lim.df = ice.lim.df, stat = 'dist.to.border')
## remove unrelalitc values
pred.fut.mod6$pred.cover[pred.fut.mod6$pred.cover < 0]  <- 0 
pred.fut.mod6$pred.cover[pred.fut.mod6$pred.cover > 1500]  <- 1500

## convert into rasters
pred.fut.mod6.list <- build_stack_from_pred(pred_ = pred.fut.mod6$pred.cover, df_ = pred.fut.mod6, output.dir_ = "mod6")
pred.fut.mod6.stk <- stack(pred.fut.mod6.list)

## produce a gif
years_ <- as.numeric(sub("^.*_", "", names(pred.fut.mod6.stk)))

##+ , eval = TRUE
saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod6.list[[t_]], main = list(years_[t_]), 
                    at = seq(0,1500, by = 50), par.settings = myThemeCont))
  }
}, movie.name = "mod6_1850_2100_MB.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)

saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod6.list[[t_]] > 100, main = list(years_[t_]), 
                    at = seq(0,1, by = 0.5), par.settings = myThemeCont))
  }
}, movie.name = "mod6_1850_2100_MB_bin.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)


##'  distance to the ice border with hierarchical random effect
mod7 <- lmer(dist.to.border ~ year_rel + alti + slope + I(alti^2) + I(slope^2) + year_rel:alti + year_rel:slope + 
               + (1 + year_rel|cellID/glacialID) + (alti + slope|glacialID), 
             data = mb.glac.dat.df)
mod7
r.squaredGLMM(mod7)
save(mod7, file = 'mod7.RData')

mod7.best <- mod7

## do data.frame prediction
pred.fut.mod7 <- cbind(mb.glac.pred.dat.df, pred.cover = predict(mod7.best, mb.glac.pred.dat.df))
## keep 2008 values for points over ELA
pred.fut.mod7 <- update.pred.cover(pred.fut.mod7, obs.df = mb.glac.dat.df, ice.lim.df = ice.lim.df, stat = 'dist.to.border')
## remove unrelalitc values
pred.fut.mod7$pred.cover[pred.fut.mod7$pred.cover < 0]  <- 0 
max.dist <- max(pred.fut.mod7$pred.cover)
pred.fut.mod7$pred.cover[pred.fut.mod7$pred.cover > 1500]  <- 1500


## convert into rasters
pred.fut.mod7.list <- build_stack_from_pred(pred_ = pred.fut.mod7$pred.cover, df_ = pred.fut.mod7, output.dir_ = "mod7")
pred.fut.mod7.stk <- stack(pred.fut.mod7.list)

## produce a gif
years_ <- as.numeric(sub("^.*_", "", names(pred.fut.mod7.stk)))

##+ , eval = TRUE
saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod7.list[[t_]], main = list(years_[t_]), 
                    at = seq(0,1500, by = 50), par.settings = myThemeCont))
  }
}, movie.name = "mod7_1850_2100_MB.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)

saveGIF({
  for(t_ in 1:length(years_)){
    print(levelplot(pred.fut.mod7.list[[t_]] > 100, main = list(years_[t_]), 
                    at = seq(0,1, by = 0.5), par.settings = myThemeCont))
  }
}, movie.name = "mod7_1850_2100_MB_bin.gif", interval = 0.33, nmax = 200, ani.width = 600,
ani.height = 600)



# 
# ##' ### test with a simple linear model
# 
# mb.glac.dat.df_ <- mb.glac.dat.df %>% mutate(glacialID = as.factor(glacialID))
# mod5 <- lm(cover ~ year_rel + alti + slope + I(alti^2) + I(slope^2) + year_rel:alti + 
#                year_rel:slope + glacialID + glacialID:year + glacialID:alti + glacialID:slope, 
#               data = mb.glac.dat.df_)
# 
# mod5.best <- step(mod.lm)
# summary(mod5.best)
# 
# ## do data.frame prediction
# mb.glac.pred.dat.df_ <- mb.glac.pred.dat.df %>% mutate(glacialID = as.factor(glacialID)) 
# pred.fut.mod5 <- cbind(mb.glac.pred.dat.df_, pred.cover = predict(mod5.best, mb.glac.pred.dat.df_))
# pred.fut.mod5$pred.cover[pred.fut.mod5$pred.cover < 0]  <- 0 
# pred.fut.mod5$pred.cover[pred.fut.mod5$pred.cover > 100]  <- 100
# 
# ## convert into rasters
# pred.fut.mod5.list <- build_stack_from_pred(pred_ = pred.fut.mod5$pred.cover, df_ = pred.fut.mod5, output.dir_ = "mod5")
# pred.fut.mod5.stk <- stack(pred.fut.mod5.list)
# 
# ## produce a gif
# years_ <- as.numeric(sub("^.*_", "", names(pred.fut.mod5.stk)))
# 
# ##+ , eval = FALSE
# saveGIF({
#   for(t_ in 1:length(years_)){
#     print(levelplot(pred.fut.mod5.list[[t_]], main = list(years_[t_]), 
#                     at = my.at, par.settings = myThemeCont))
#   }
# }, movie.name = "mod5_1850_2100_MB.gif", interval = 0.33, nmax = 200, ani.width = 600,
# ani.height = 600)




# ##' ## Bonus graphs
# 
# cover.thresh_ <- 25
# glacialID_ <- 103
# var.y_ <- "alti_fact"
# var.x_ <- "year"
# dots.grp <- sapply(c(var.y_, var.x_), . %>% {as.formula(paste0('~', .))})
# 
# gg.dat <- mb.glac.dat.df %>% 
#   mutate(alti_fact = cut(alti, seq(1000, 5000, 250), dig.lab = 5),
#          slope_fact = cut(slope, seq(0, 1.25, .10), dig.lab = 3)) %>%
#   filter(glacialID == glacialID_, cover > cover.thresh_) %>%
#   group_by_(.dots = dots.grp) %>% summarize(nb.glac.pix = n())
# 
# gg <- ggplot()

##' ## Evaluate models perfs

##' overall evaluation on the 4 dates mixed up

valid.cover.dat <- rbind(left_join(mb.glac.dat.df %>% mutate(obs.cover = cover) %>% select(cellID, year, obs.cover), pred.fut.mod1 %>% mutate(mod = "mod1")),
                   left_join(mb.glac.dat.df %>% mutate(obs.cover = cover) %>% select(cellID, year, obs.cover), pred.fut.mod2 %>% mutate(mod = "mod2")),
                   left_join(mb.glac.dat.df %>% mutate(obs.cover = cover) %>% select(cellID, year, obs.cover), pred.fut.mod3 %>% mutate(mod = "mod3"))) 

valid.dist.to.border.dat <- rbind(left_join(mb.glac.dat.df %>% mutate(obs.dist.to.border = dist.to.border) %>% select(cellID, year, obs.dist.to.border), pred.fut.mod5 %>% mutate(pred.dist.to.border = pred.cover, mod = "mod5") %>% select(-pred.cover)),
                                  left_join(mb.glac.dat.df %>% mutate(obs.dist.to.border = dist.to.border) %>% select(cellID, year, obs.dist.to.border), pred.fut.mod6 %>% mutate(pred.dist.to.border = pred.cover, mod = "mod6") %>% select(-pred.cover)),
                                  left_join(mb.glac.dat.df %>% mutate(obs.dist.to.border = dist.to.border) %>% select(cellID, year, obs.dist.to.border), pred.fut.mod7 %>% mutate(pred.dist.to.border = pred.cover, mod = "mod7") %>% select(-pred.cover))) 

##+

mod.eval <- data.frame()
for(mod_ in paste0("mod", 1:3)){
  cat("\n", mod_)
  mod.eval <- rbind(mod.eval, 
                    cbind(biomod2::Find.Optim.Stat(Stat = 'TSS', 
                                             Obs = as.numeric(valid.cover.dat %>% filter(mod == mod_) %>% .[['obs.cover']] > 0), 
                                             Fit = valid.cover.dat %>% filter(mod == mod_) %>% .[['pred.cover']]), mod = mod_),
                    cbind(biomod2::Find.Optim.Stat(Stat = 'ROC', 
                                             Obs = as.numeric(valid.cover.dat %>% filter(mod == mod_) %>% .[['obs.cover']] > 0), 
                                             Fit = valid.cover.dat %>% filter(mod == mod_) %>% .[['pred.cover']]), mod = mod_))
}

##+

for(mod_ in paste0("mod", 5:7)){
  cat("\n", mod_)
  mod.eval <- rbind(mod.eval,
                    cbind(biomod2::Find.Optim.Stat(Stat = 'TSS',
                                             Obs = as.numeric(valid.dist.to.border.dat %>% filter(mod == mod_) %>% .[['obs.dist.to.border']] > 0),
                                             Fit = valid.dist.to.border.dat %>% filter(mod == mod_) %>% .[['pred.dist.to.border']]), mod = mod_),
                    cbind(biomod2::Find.Optim.Stat(Stat = 'ROC', 
                                             Obs = as.numeric(valid.dist.to.border.dat %>% filter(mod == mod_) %>% .[['obs.dist.to.border']] > 0), 
                                             Fit = valid.dist.to.border.dat %>% filter(mod == mod_) %>% .[['pred.dist.to.border']]), mod = mod_))
}

save(mod.eval, file = "mod.eval.RData")
knitr::kable(mod.eval)

##' ## Plot the glacial area predoctons in 2100 to spot best models

r.m1 <- raster("~/Work/MB_ATLAS_2016/_workdir/mod1/pred_2100.tif")
r.m2 <- raster("~/Work/MB_ATLAS_2016/_workdir/mod2/pred_2100.tif")
r.m3 <- raster("~/Work/MB_ATLAS_2016/_workdir/mod3/pred_2100.tif")

r.m5 <- raster("~/Work/MB_ATLAS_2016/_workdir/mod5/pred_2100.tif")
r.m6 <- raster("~/Work/MB_ATLAS_2016/_workdir/mod6/pred_2100.tif")
r.m7 <- raster("~/Work/MB_ATLAS_2016/_workdir/mod7/pred_2100.tif")

r.m1.bin <- r.m1 > as.numeric(mod.eval$cutoff[mod.eval$mod == 'mod1' & grepl("TSS", rownames(mod.eval))])
r.m2.bin <- r.m2 > as.numeric(mod.eval$cutoff[mod.eval$mod == 'mod2' & grepl("TSS", rownames(mod.eval))])
r.m3.bin <- r.m3 > as.numeric(mod.eval$cutoff[mod.eval$mod == 'mod3' & grepl("TSS", rownames(mod.eval))])

r.m5.bin <- r.m5 > as.numeric(mod.eval$cutoff[mod.eval$mod == 'mod5' & grepl("TSS", rownames(mod.eval))])
r.m6.bin <- r.m6 > as.numeric(mod.eval$cutoff[mod.eval$mod == 'mod6' & grepl("TSS", rownames(mod.eval))])
r.m7.bin <- r.m7 > as.numeric(mod.eval$cutoff[mod.eval$mod == 'mod7' & grepl("TSS", rownames(mod.eval))])


stk <- stack(r.m1.bin, r.m2.bin, r.m3.bin, r.m5.bin, r.m6.bin, r.m7.bin)
names(stk) <- paste0("mod", c(1:3, 5:7))
levelplot(stk)

##' ## Summary
##' 
##' As we see in this analysis, models that based on percentage of ice cover (mod1, mod2, mod3 and mod5) strugled to predict 
##' ice removal in the future (not enough melting after 2008). Models based on distance to the border (mod6, mod7) of the glacier seems to 
##' work much better on this points and keep in the mean time a high predictability power according to TSS and ROC scores.
##' As both mod6 and mod7 have comparable evaluation scores, we decided to keep the more parcimonious one (i.e. mod6). At the 
##' end we will use the model:
##' 
##' `mod6 <- lmer(dist.to.border ~ year_rel * alti + slope  + (1+year_rel|cellID), data = mb.glac.dat.df)`
##' 
##' as reference model to produce ice cover future masks. This model is very similar to Carlson et al. one exept that we focus now
##' on distance to glacier border instead of percentage of cover.
##' 

##' ## Production of ice mask for FATEHD

ice.masks.outdir <- "~/Work/MB_ATLAS_2016/_outputs/MB_glacier_masks"
dir.create(ice.masks.outdir, showWarnings = FALSE, recursive = TRUE)

cont.ice.masks.files <- list.files("mod6", pattern = "^pred_[0-9]{4}.tif$", full.names = TRUE)
bin.ice.masks.files <- paste0("MB_glacier_mask_", sub("pred_", "", basename(cont.ice.masks.files)))
names(bin.ice.masks.files) <- basename(cont.ice.masks.files)

for(f_ in cont.ice.masks.files){
  r_ <- raster(f_, RAT = FALSE)
  writeRaster(r_ > as.numeric(mod.eval$cutoff[mod.eval$mod == 'mod6' & grepl("TSS", rownames(mod.eval))]),
              filename = file.path(ice.masks.outdir, bin.ice.masks.files[basename(f_)]),
              overwrite = TRUE)
}


##' -- end of script -----------------------------------------------------------