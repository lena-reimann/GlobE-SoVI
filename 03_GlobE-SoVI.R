######################################
#### produce final GlobE-SoVI map ####
######################################
# by Lena Reimann
# October 31, 2023

# This script corresponds to the data processing described in '3) GlobE-SoVI' of
# Reimann et al. "An empirical social vulnerability map for flood risk assessment at global scale (‘GlobE-SoVI’)".
# The final output is a global social vulnerability map calculated as the weighted sum of all relevant 
# vulnerability variables, scaled from 1-10, both in raster format (i.e. GeoTIFF) and vector format (i.e. shapefile) per administrative unit.

rm(list=ls())

path = "./path_to_data/"
lib = "./path_to_R_libraries"

# packages
library(sp, lib.loc = lib)
library(raster, lib.loc = lib)
library(dplyr, lib.loc = lib)
library(sf, lib.loc = lib)
library(rgdal, lib.loc = lib)
library(iterators, lib.loc = lib)
library(scales, lib.loc = lib)
library(igraph, lib.loc = lib)


####----------- 1.) Load data -------------####

## regression coefficients ##
coef = read.csv(paste0(path, "regression_coefficients.csv"))

## base data ##
# water mask (https://sedac.ciesin.columbia.edu/data/set/gpw-v4-land-water-area-rev11/data-download)
name = "gpw_v4_data_quality_indicators_rev11_watermask_30_sec.tif"
name = paste(path, name, sep = "/")
wat = raster(name)

# GADM administrative unit polygons (https://gadm.org/download_world.html)
adm = st_read(paste0(path, "gadm_410.gdb"), layer = "gadm")


## relevant vulnerability characteristics ##
# age groups
ag1 <- sprintf("ag%s", seq(0,60, by= 5))
ag1[length(ag1)+1] <- paste0(sprintf("ag%s", length(ag1) * 5), "+")
ag1 = c(paste0(ag1, "f"), paste0(ag1, "m"))

name = paste(ag1, "prop.tif", sep = "_")
name = paste(path, name, sep = "/")
prop_ag = stack(name)

# proportions f and m
name = c("f_prop.tif", "m_prop.tif")
name = paste(path, name, sep = "/")
prop = stack(name)

# education
name = c("mschf2010.tif", "mschm2010.tif")
name = paste(path, name, sep = "/")
mys = stack(name)

# income
name = c("gnicf2010.tif", "gnicm2010.tif")
name = paste(path, name, sep = "/")
inc = stack(name)

# settlements
name = "GHS_SMOD2015_WGS.tif" # year 2015 as the year closest to 2010
name = paste(path, name, sep = "/")
settl = raster(name)

# walking time to health facility
name = "2020_walking_only_travel_time_to_healthcare.geotiff"
name = paste(path, name, sep = "/")  
tvl = raster(name)


####------------------ 2.) Functions --------------------####

setNA <- function(x) {x[x == 0] <- NA; return(x)}
setNA_edu <- function(x) {x[x > 20] <- NA; return(x)}


# scale social vulnerability raster from 1-10
vul_sc = function(vul, p) { #vul = vulnerability raster with positive and negative values; 
  #p = quantiles to calculate: provide vector of two if outliers should be removed, otherwise vector of one
  require(raster)
  
  # 1. remove outliers
  if (length(p) == 2) {
    # calculate quantiles
    qu = quantile(vul, probs = p, na.rm = T, ncells = (vul@ncols * vul@nrows / 100))

    # establish min and max for reclassify
    max = vul@data@max
    min = vul@data@min
    
    # reclassify raster based on quantiles
    rcl = c((min-1), qu[1], qu[1], qu[2], max+1, qu[2])
    rcl = matrix(rcl, ncol=3, byrow=TRUE)
    vul = reclassify(vul, rcl)
  }
  
  # 2. make all values start from 0
  max = vul@data@max
  min = vul@data@min
  
  if (min < 0) {
    # make raster to transform data into positive domain
    rcl = c(min-1, max+1, (min*-1)) 
    rcl = matrix(rcl, ncol=3, byrow=TRUE)
    fct = reclassify(vul, rcl)
    
    vul = overlay(vul, fct, fun = function(x, y) {return(x + y)}) # values range from 0 to max
    
  } else if (min > 0) {
    
    # make raster to transform data starting from 0
    rcl = c(min-1, max+1, min) 
    rcl = matrix(rcl, ncol=3, byrow=TRUE)
    fct = reclassify(vul, rcl)
    
    vul = overlay(vul, fct, fun = function(x, y) {return(x - y)}) # values range from 0 to max-min
  }
  
  # 3. scale values
  #establish max
  min = vul@data@min
  max = vul@data@max
  
  #calculate scaling factor
  fct = 9 / max
  
  rcl = c(min-1, max+1, fct)
  rcl = matrix(rcl, ncol=3, byrow=TRUE)
  fct = reclassify(vul, rcl)
  
  #scale vul values
  vul_sc = overlay(vul, fct, fun = function(x, y) {return(x * y + 1)})
  
  return(vul_sc)
}


####------------------ 3.) Analysis --------------------####

## 1. produce variable rasters (and write) ##

# produce land mask from water mask
rcl = c(-0.5, 2.5, 1,  2.5, 3.5, NA)
rcl = matrix(rcl, ncol=3, byrow=TRUE)

name = paste(path, "gpw_v4_land.tif", sep = "/")
lan = reclassify(wat, rcl, filename = name, overwrite = T)

# education
mys.gen = overlay(prop, mys, fun = function(x, y) {return(x * y / 10)})
mys.gen = calc(mys.gen, sum, na.rm = T)
mys.gen = calc(mys.gen, fun = setNA)
mys.gen = calc(mys.gen, fun = setNA_edu)
writeRaster(mys.gen, paste0(path, "MYS_gender.tif"), overwrite = T)

# income
inc.gap = overlay(inc[[1]], inc[[2]], fun = function(x,y) {return((y - x) / y * 100)})
inc.gap = calc(inc.gap, fun = setNA) # set income gap == 0 to NA because in those instances we do not have income data per gender
writeRaster(inc.gap, paste0(path, "income_gap.tif"), overwrite = T)

# settlements (SMOD class 13) 
rcl = c(0, 12, 0,  12.5, 13.5, 1,  14, 50, 0)
rcl = matrix(rcl, ncol=3, byrow=TRUE)
settl = reclassify(settl, rcl)
settl = mask(settl, lan)
writeRaster(settl, paste0(path, "rural_settlements.tif"), overwrite = T)

# elderly
eld.tot = stack(prop_ag[[14]], prop_ag[[28]])
eld.tot = calc(eld.tot, fun = function(x) {x = sum(x, na.rm = T) * 100; return(x)})
eld.tot = mask(eld.tot, lan)
writeRaster(eld.tot, paste0(path, "percent_elderly.tif"), overwrite = T)

# walking time
tvl = extend(tvl, eld.wgt) # adjust extents
tvl = calc(tvl,  fun = function(x) {x / 60})
writeRaster(tvl, paste0(path, "walk_healthcare_hours.tif"), overwrite = T)


## 2. calculate GlobE-SoVI ##

# raster weighting based on regression coefficients
settl.wgt = calc(settl,  fun = function(x) {x * 0.3958673}) #percent exposed transformed to presence multiplied by mean exposure (=0.0830145816717685) * (mean(smod_rel)== 4.771297
edu.wgt = calc(mys.gen,  fun = function(x) {x * -0.209976450721468})
eld.wgt = calc(eld.tot,  fun = function(x) {x * 0.0540403921327649})
wlk.wgt = calc(tvl,  fun = function(x) {x * 0.00514173490948157})
inc.wgt = calc(inc.gap,  fun = function(x) {x * 0.0180278658712546})

# vulnerability sum
vul = sum(edu.wgt, inc.wgt, settl.wgt, eld.wgt, wlk.wgt)

# scale vulnerability map to 1-10
p = c(0.01, 0.99)
vul = vul_sc(vul, p)
writeRaster(vul, paste0(path, "GlobE-SoVI.tif"), overwrite = T)


## 3. aggregate to administrative unit level

# create raster stack
char = stack(mys.gen, inc.gap, settl, eld.tot, tvl, vul)

# zonal statistics
adm = adm[,1]
char = zonal(char, adm, fun = mean, na.rm = T, df = T)

# change colnames
colnames(char) <- c("ID", "MYS_gen", "inc_gap", "rur_settl", "eld_per", "wlk_hrs", "GlobE-SoVI")

# make settlement percentage
char$rur_settl <- char$rur_settl * 100

# make spatial
adm.char = merge(adm, char, by.x = "UID", by.y = "ID")

# write csv and shp
write.csv(char,  paste0(path, "GlobE-SoVI.csv"))
write_sf(adm.char, paste0(path, "GlobE-SoVI.shp"))

