######################################
#### produce final GlobE-SoVI map ####
######################################
# by Lena Reimann
# May 31, 2023

# This script corresponds to the data processing described in '3) GlobE-SoVI' of
# Reimann et al. "An empirical social vulnerability map (‘GlobE-SoVI’) for flood risk assessment at global scale".
# The final output is a global social vulnerability raster calculated as the weighted sum of all relevant 
# vulnerability variables scaled from 1-10 based on the 98th percent confidence interval (i.e. 1st and 99th percentiles removed).

rm(list=ls())

path = "./path_to_data/"
lib = "./path_to_R_libraries"

# packages
library(sp, lib.loc = lib)
library(raster, lib.loc = lib)
library(terra, lib.loc = lib)
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

set0 <- function(x) {x[is.na(x)] <- 0; return(x)}
setNA <- function(x) {x[x == 0] <- NA; return(x)}
setNA_edu <- function(x) {x[x > 20] <- NA; return(x)}

# scale social vulnerability raster
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
    
    vul = overlay(vul, fct, fun = function(x, y) {return(x + y)}) ## values range from 0 to max
    
  } else if (min > 0) {
    
    # make raster to transform data starting from 0
    rcl = c(min-1, max+1, min) 
    rcl = matrix(rcl, ncol=3, byrow=TRUE)
    fct = reclassify(vul, rcl)
    
    vul = overlay(vul, fct, fun = function(x, y) {return(x - y)}) ## values range from 0 to max-min
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

#### produce V index ####

# produce land mask from water mask
rcl = c(-0.5, 2.5, 1,  2.5, 3.5, NA)
rcl = matrix(rcl, ncol=3, byrow=TRUE)

name = "gpw_v4_land.tif"
name = paste(path, name, sep = "/")
lan = reclassify(wat, rcl, filename = name, overwrite = T)

# edu 
mys.gen = overlay(prop, mys, fun = function(x, y) {return(x * y / 10)})
mys.gen = calc(mys.gen, sum, na.rm = T)
mys.gen = calc(mys.gen, fun = setNA)
mys.gen = calc(mys.gen, fun = setNA_edu)
writeRaster(mys.gen, paste0(path_vul, "/mys_gen.tif"), overwrite = T)

edu.wgt = calc(mys.gen,  fun = function(x) {x * -0.199541185387587})
writeRaster(edu.wgt, paste0(path_vul, "/edu_wgt.tif"), overwrite = T)

# income
inc.gap = overlay(inc[[1]], inc[[2]], fun = function(x,y) {return((y - x) / y * 100)})

# set income gap == 0 to NA because in those instances we do not have income data per gender
inc.gap = calc(inc.gap, fun = setNA)
writeRaster(inc.gap, paste0(path_vul, "/inc_gap.tif"), overwrite = T)

inc.wgt = calc(inc.gap,  fun = function(x) {x * 0.0209104692643852})
#inc.wgt = mask(inc.wgt, lan)
writeRaster(inc.wgt, paste0(path_vul, "/inc_wgt.tif"), overwrite = T)
#inc.wgt = raster(paste0(path_vul, "/inc_wgt.tif"))


# settl13 --> percent exposed transformed to presence multiplied by mean exposure
#settl = crop(settl, e)
rcl = c(0, 12, 0,  12.5, 13.5, 1,  14, 50, 0)
rcl = matrix(rcl, ncol=3, byrow=TRUE)
settl = reclassify(settl, rcl)
settl = mask(settl, lan)

writeRaster(settl, paste0(path_vul, "/settl_smod13.tif"), overwrite = T)

settl.wgt = calc(settl,  fun = function(x) {x * 0.396087224486764}) # coef from lm (=0.0830145816717685) * (mean(smod_rel)== 4.771297
writeRaster(settl.wgt, paste0(path_vul, "/settl_wgt.tif"), overwrite = T)
#settl.wgt = raster(paste0(path_vul, "/settl_wgt.tif"))


# elderly tot proportion --> percent
#prop_ag = crop(prop_ag, e)
eld.tot = stack(prop_ag[[14]], prop_ag[[28]])
eld.tot = calc(eld.tot, fun = function(x) {
  x = sum(x, na.rm = T) * 100; return(x)})
eld.tot = mask(eld.tot, lan)
writeRaster(eld.tot, paste0(path_vul, "/eld_tot_perc.tif"), overwrite = T)

eld.wgt = calc(eld.tot,  fun = function(x) {x * 0.0469823888996626})
writeRaster(eld.wgt, paste0(path_vul, "/eld_wgt.tif"), overwrite = T)
#eld.wgt = raster(paste0(path_vul, "/eld_wgt.tif"))


# walk --> from minutes to hours
# adjust extents
#tvl = crop(tvl, e)
tvl = extend(tvl, eld.wgt)
name = paste0(path_vul, "/wlk_glo_hrs.tif")

tvl = calc(tvl,  fun = function(x) {x / 60})
writeRaster(tvl,filename = name, overwrite = T)

wlk.wgt = calc(tvl,  fun = function(x) {x * 0.00520586310110901})
writeRaster(wlk.wgt, paste0(path_vul, "/wlk_wgt.tif"), overwrite = T)
#wlk.wgt = raster(paste0(path_vul, "/wlk_wgt.tif"))



## vulnerability sum ##
vul = sum(edu.wgt, inc.wgt, settl.wgt, eld.wgt, wlk.wgt)#, na.rm = T) #
writeRaster(vul, paste0(path_vul, "/vul_fat.log_sum_ngov.tif"), overwrite = T)
#vul = raster(paste0(path_vul, "/vul_fat.log_sum_ngov.tif"))


# scale vulnerability map to 0-1
# a) original values
p = 0
vul_ori = vul_sc(vul, p)
writeRaster(vul_ori, paste0(path_vul, "/vul_sum_ori_ngov_sc.tif"), overwrite = T)

# b) 1st and 99th qu removed
p = c(0.01, 0.99)
vul_out = vul_sc(vul, p)
writeRaster(vul_out, paste0(path_vul, "/vul_sum_out_ngov_sc_p1_NA.tif"), overwrite = T)


## vulnerability exp ## --> final vul
vul = calc(vul, fun = function(x) {exp(x)})
vul = mask(vul, lan)
writeRaster(vul, paste0(path_vul, "/vul_fat.log_exp_ngov.tif"), overwrite = T)


# scale vulnerability map to 0-1
# a) original values
p = 0
vul_ori = vul_sc(vul, p)
writeRaster(vul_ori, paste0(path_vul, "/vul_exp_ori_ngov_sc.tif"), overwrite = T)

# b) 1st and 99th qu removed
p = c(0.01, 0.99)
vul_out = vul_sc(vul, p)
writeRaster(vul_out, paste0(path_vul, "/vul_exp_out_ngov_sc_p1_NA.tif"), overwrite = T)





