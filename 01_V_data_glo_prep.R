###################################################
#### analyze spatial patterns of observed data ####
###################################################
# by Lena Reimann
# Oct 29, 2022

# Goal: process all V datasets (i.e. age + sex; edu; inc; urb/rur) for the entire globe
# NOTES: now also includes dep ratios (Mar 21, 2022)


## define basic job characterstics ##

rm(list=ls())

# setwd
run = "cluster" # change here if run on PC

if (run == "cluster") {
  path_vul = "./vulnerability/global"
  
  #.libPaths("~/r_libs")
  lib = "./r_libs"
  
} else {
  path <- "U:/Data/"
  #path <- "C:/Users/admin/Documents/Lena/process_vul_data_GLOBAL/" #remote desktop
  setwd(path)
  
  path_vul = "Population/vulnerability/global"
  
  #.libPaths()
  lib = "C:/Users/lrn238/AppData/Local/RLIB"  
}

# packages
library(sp, lib.loc = lib)
library(raster, lib.loc = lib)
library(dplyr, lib.loc = lib)
library(sf, lib.loc = lib)
library(rgdal, lib.loc = lib)



####----------- 1.) load data -------------####
 
# admin units
#iso codes
#file = "data/admin_units/ISO3166_codes.csv"
#isos = read.csv(file)

#shp
#file = "data/admin_units/gadm36_0.shp" #edited: Code for Kosovo changed to KOS
#adm0_shp <- read_sf(file)

#tif
name = "gpw/gpw_v4_national_identifier_grid_rev11_30_sec.tif"
file = paste(path_vul, name, sep = "/")
adm0_tif <- raster(file)

print(paste(name, "loaded"))

#eucl allocation raster of admin ids
#name = "gpw/gpw_v4_nat_eucl.tif"
#file = paste(path_vul, name, sep = "/")
#adm0_ecl <- raster(file)

#buffered admin units
name = "gpw/gpw_v4_nat_id_buf_mask.tif"
file = paste(path_vul, name, sep = "/")
adm0_cst <- raster(file)

print(paste(name, "loaded"))


# edu + inc 
name = c("gdl/mschf2010_gfil_EK.tif", "gdl/mschm2010_gfil_EK.tif") # years of schooling multiplied by 10 for int values and gaps filled with hdi (Elco's approach)
file = paste(path_vul, name, sep = "/")
edu = stack(file)

print(paste(name, "loaded"))

name = c("gdl/gnicf2010_gfil_EK.tif", "gdl/gnicm2010_gfil_EK.tif") # gaps filled with hdi (Elco's approach)
file = paste(path_vul, name, sep = "/")
inc = stack(file)

print(paste(name, "loaded"))


# age + sex
pop = "gpw/gpw_v4_basic_demographic_characteristics_rev11_"
data_type <- "cntm_30_sec"
t  <- 2010

# define ag to be loaded
ag1 <- sprintf("ag%s", seq(0,60, by= 5))
ag1[length(ag1)+1] <- paste0(sprintf("ag%s", length(ag1) * 5), "+")
ag1 = c(paste0(ag1, "f"), paste0(ag1, "m"))

ag2 <- sprintf("a0%s", seq(0,60, by= 5))
ag2[1:2] <- c("a000", "a005")
ag2[length(ag2)+1] <- paste0(sprintf("a0%s", length(ag2) * 5))

ag3 <- sprintf("0%s", seq(4,64, by= 5))
ag3[c(1:2)] <- c("004", "009")
ag3[length(ag3)+1] <- "plus"

ag_f <- stack(paste(paste(paste0(paste(path_vul, pop, sep = "/"), ag2, "_", ag3, "ft"), t, data_type, sep = "_"), 
                    "tif", sep = "."))
print(paste(pop, "female loaded"))

ag_m <- stack(paste(paste(paste0(paste(path_vul, pop, sep = "/"), ag2, "_", ag3, "mt"), t, data_type, sep = "_"), 
                    "tif", sep = "."))
print(paste(pop, "male loaded"))

# all age groups in one stack
ag <- stack(ag_f, ag_m)
# name layers
names(ag) = ag1

remove(ag_f, ag_m)


# settlements
name = c("SMOD/GHS_SMOD2000_WGS.tif", "SMOD/GHS_SMOD2015_WGS.tif") #smod projected to wgs coordinates
file = paste(path_vul, name, sep = "/")
settl = raster(file)


####------------------ 2.) Functions --------------------####


## admin units ##

# produce raster with admin id for the respective cntr
adm0_id <- function(adm0_shp, adm0_tif, id) { # input: global adm id polygons and tif, cntr id of respective cntr
  require(raster, sf)
  
  adm0_cntr = adm0_shp[adm0_shp$ID == id, ]
  adm0_cntr = crop(adm0_tif, adm0_cntr)

  rcl <- rbind(c(-Inf, id-0.5, NA), c(id+0.5, Inf, NA))
  adm0_cntr <- reclassify(adm0_cntr, rcl)

  return(adm0_cntr)
}


#  make admin unit mask
adm_m = function(adm0_ecl, adm0_cst) { # input: admin unit raster with eucl allo of cntr ids; buffered admin unit raster
  require(raster)

  # set all values in coastal mask to 1
  adm0_cst = calc(adm0_cst, fun = function(x) { x[!is.na(x)] <- 1; return(x)})
  
  # mask eucl allo raster with coastal mask
  adm0_ecl = mask(adm0_ecl, adm0_cst)
  
  return(stack(adm0_ecl, adm0_cst))
}



## Age and sex ##

# calculate AG or sex proportions (raster-based)
prop <- function(stack, aim)  { # input: rasters female and male combined in one stack; aim ("age" or "sex")

  # if aim == sex make stack with sexes summed up
  if (aim == "sex") {
    # sum up all AGs per sex and stack again
    f = calc(stack[[1 : (nlayers(stack) / 2)]], fun = function(x) sum(x))
    m = calc(stack[[(nlayers(stack) / 2 + 1) : nlayers(stack)]], fun = function(x) sum(x))
    stack = stack(f, m)
  }

  # sum up all pop (cross-check for t2010)
  tot <- calc(stack, fun = function(x) sum(x))
  
  # calculate AG proportions per cell --> sum of all AG == 1
  prop <- overlay(stack, tot, fun = function(x,y) x / y)
  
  if (aim == "sex") {
    return(list(prop, stack))
  } else {
    return(list(prop, tot))
  }
}

# calculate dependency ratios
dep_ratio = function(stack)  { # input: stack of age group proportions (f and m combined)
  
  # a. separate female versus male AGs 
  prop_f = stack[[1 : (nlayers(stack)/2)]]
  prop_m = stack[[(nlayers(stack)/2 + 1) : nlayers(stack)]]
  
  # b. combine them to come to total prop per AG --> nice to have!
  # define name
  # ags
  ag_tot <- sprintf("ag%s", seq(0,60, by= 5))
  ag_tot[length(ag_tot)+1] <- paste0(sprintf("ag%s", length(ag_tot) * 5), "+")
  
  #overlay
  prop_tot = overlay(prop_f, prop_m, fun = sum)#, filename = name, bylayer = T) 
  
  # c. sum up prop per AGs
  # establish layers for the AGs needed
  chi = c(1:3)
  wrk = c(4: (nlayers(prop_tot)-1))
  eld = c(nlayers(prop_tot))
  
  #calculate broad AG props
  chi = calc(prop_tot[[chi]], fun = function(x) {sum(x, na.rm = T)}) #proportion of children (0-14)
  wrk = calc(prop_tot[[wrk]], fun = function(x) {sum(x, na.rm = T)}) #proportion of working age pop (15-64)
  eld = calc(prop_tot[[eld]], fun = function(x) {sum(x, na.rm = T)}) #proportion of elderly (65+)
  
  dep_ags = stack(chi, wrk, eld)
  
  # d. calc dep ratios
  cdr = overlay(chi, eld, wrk, fun = function(x, y, z) {(x + y) / z})
  ydr = overlay(chi, wrk, fun = function(x, y) {x / y})
  odr = overlay(eld, wrk, fun = function(x, y) {x / y})  
  
  dep_rat = stack(cdr, ydr, odr)
  
  return(list(prop_tot, dep_ags, dep_rat)) # return list with total AG prop stack, 
  #stack of summed up AGs for dep ratios, 
  #stack of dep ratios
}


## education ##

# scale MYS for different AGs
sc_edu = function(edu, ag, mask, ag_name, path) { #input: stack of MYS and pop numbers per AG, vector with ag names, path to write to

  # ags to do this for 
  #ags = c("ag0", "ag5", "ag10", "ag15", "ag20", "ag25+") 
  
  for (i in 1:(nlayers(ag)/2)) {
    # extract f and m layers of the respective AG
    sex = stack(ag[[i]], ag[[nlayers(ag)/2 + i ]])
    
    # set edu to 0 where ag 0
    edu_ag = overlay(edu, sex, fun = function(x,y) {x[y==0] <- NA; return(x)})
    
    # multiplication factor
    fac = ifelse(i == 1, 0, (i - 0.5) / 5) # 5 = AG interval

    # scale for the first 5 AGs
    if (i <= 5) {
      edu_ag = edu_ag * fac
    } else {
      edu_ag = edu_ag
    }

    # divide by 10 to get to YS and write
    edu_ag = calc(edu_ag, fun = function(x) {x / 10})
    
    # mask data to remove very high values instead of NAs in water cells
    #edu_ag = mask(edu_ag, mask)
    edu_ag = calc(edu_ag, fun = function(x) {x[x>20] <- NA; return(x)})
    
    # write raster
    name = paste(paste("edu", c(ag_name[[i]], ag_name[[length(ag_name)/2 + i ]]), sep = "_"), "tif", sep = ".")
    name = paste(path, name, sep = "/")  
    writeRaster(edu_ag, filename = name, format = "GTiff", bylayer = T, overwrite = T)
  }
}


# calculate MYS per cell 
edu_mean = function(edu_sc, ag) { #input: raster stack with MYS scaled per AG and stack with AG numbers
  
  edu_sc = overlay(edu_sc, ag, fun = function(x,y) {x * y})
  
  edu_tot = calc(edu_sc, fun = function(x) { sum(x, na.rm = T)})
  pop_tot = calc(ag, fun = function(x) { sum(x, na.rm = T)})
  
  edu_mean = overlay(edu_tot, pop_tot, fun = function(x,y) {x / y})
  
  return(edu_mean)
}


## income ##

# produce income per sex raster
inc_sex = function(inc, sex) { #input: stack of income and pop per sex
  
  # set inc to 0 where ag 0
  inc_sex = overlay(inc, sex, fun = function(x,y) {x[y==0] <- NA; return(x)})
  
  return(inc_sex)
}


## settlements ##
settl_cl = function(settl) { #input: raster with codes for each settlement type
  # values to reclass
  val = unique(settl)
  
  #create empty stack
  settl_st = stack()
  
  for (i in 1:length(val)) {
    #reclassify
    rcl = rbind(c(-Inf, val[i]-0.5, NA), c(val[i]+0.5, Inf, NA))
    rcl = reclassify(settl, rcl)
    
    #stack
    settl_st = stack(settl_st, rcl)
  }
  return(settl_st)
}


####-------------- 2.) data preprocessing ---------------####

start <- Sys.time()

# merge admin units with iso code IDs
#adm0_shp = merge(adm0_shp, isos, by.x = "GID_0", by.y = "ISO", all = T)


# mask eucl rasters of edu and inc ## only relevant when global analysis
#edu = mask(edu, adm0_tif)
#inc = mask(inc, adm0_tif)


####-------------- 3.) test on example cntr ---------------####

id <- 0 

# example France
#id    <- 250
#iso   <- "FRA"
#iso1  <- "fra"

# make cntr admin raster
if (id == 0) {
  adm0_cntr = adm0_tif
  remove(adm0_tif)
} else {
  adm0_cntr <- adm0_id(adm0_shp, adm0_tif, id)
}


# make admin mask (only needed for country-level analysis)
#crop to cntr
#adm0_ecl_cntr = crop(adm0_ecl, adm0_cntr)
#adm0_cst_cntr = crop(adm0_cst, adm0_cntr)

#apply function
#adm_mask = adm_m(adm0_ecl_cntr, adm0_cst_cntr)

#extract rasters
#adm0_id_buf = adm_mask[[1]]
#adm0_buf_mask = adm_mask[[2]]

#path = "Administrative_units/gpw_v4/gpw_v4_national-identifier-grid/"
#name = c("gpw_v4_id_buf.tif", "gpw_v4_buf_mask.tif")
#name = paste(path, name, sep = "/")
#writeRaster(adm_mask, filename = name, format = "GTiff", bylayer = T, overwrite = T)



#### ---- 3a. preprocess data ---- ####

## age and sex ##

# calculate proportions of AGs and total pop
# crop data to cntr
#ag_cntr = mask(crop(ag, adm0_cntr), adm0_cntr)
ag_cntr = ag
remove(ag)

prop_cntr <- prop(ag_cntr, "age")

prop_ag = prop_cntr[[1]]
#pop_tot = prop_cntr[[2]]

# write rasters
name = paste(ag1, "prop.tif", sep = "_")
name = paste(path_vul, name, sep = "/")
writeRaster(prop_cntr[[1]], filename = name, format = "GTiff", bylayer = T, overwrite = T)
#prop_ag = stack(paste0("Population/", name))


# calculate dep ratios
dep_ratio_stacks = dep_ratio(prop_ag)

# write rasters
# total proportions
ag_tot <- sprintf("ag%s", seq(0,60, by= 5))
ag_tot[length(ag_tot)+1] <- paste0(sprintf("ag%s", length(ag_tot) * 5), "+")

name = paste(ag_tot, "prop.tif", sep = "_")
name = paste(path_vul, name, sep = "/")
writeRaster(dep_ratio_stacks[[1]], filename = name, format = "GTiff", bylayer = T, overwrite = T)

# age groups for dep ratios
name = c(paste0("chi_prop.tif"), paste0("wrk_prop.tif"), paste0("eld_prop.tif"))
name = paste(path_vul, name, sep = "/")
writeRaster(dep_ratio_stacks[[2]], filename = name, format = "GTiff", bylayer = T, overwrite = T)

# dep ratios
name = c(paste0("cdr.tif"), paste0("ydr.tif"), paste0("odr.tif"))
name = paste(path_vul, name, sep = "/")
writeRaster(dep_ratio_stacks[[3]], filename = name, format = "GTiff", bylayer = T, overwrite = T)


# calculate proportions of sexes
prop_cntr = prop(ag_cntr, "sex")

#prop_sex = prop_cntr[[1]]
#pop_sex = prop_cntr[[2]]

# check whether proportions add up --> mean should be ~1
#mean(values(calc(prop_ag, fun = function(x) sum(x))), na.rm = T)
#mean(values(calc(prop_sex, fun = function(x) sum(x))), na.rm = T)

# write rasters
# proportions
name = c("f_prop.tif", "m_prop.tif")
name = paste(path_vul, name, sep = "/")
writeRaster(prop_cntr[[1]], filename = name, format = "GTiff", bylayer = T, overwrite = T)

# tot pop sex
name = c("f_tot.tif", "m_tot.tif")
name = paste(path_vul, name, sep = "/")
writeRaster(prop_cntr[[2]], filename = name, format = "GTiff", bylayer = T, overwrite = T)


## edu ##

# crop to country
#edu_cntr = mask(crop(edu, adm0_cntr), adm0_cntr)
edu_cntr = edu
remove(edu)
#ag_cntr = crop(ag, adm0_cntr)

gc()

# scale edu for each AG
sc_edu(edu_cntr, ag_cntr, adm0_tif, ag1, path_vul) #scaled education per AG (about 1h 15m per AG --> total of 18h)


# MYS per cell weighted by AG 

# load edu_sc
name = paste(paste("edu", c(ag1, ag1[length(ag1)/2]), sep = "_"), "tif", sep = ".")
name = paste(path_vul, name, sep = "/")  

edu_sc = stack(name)

# run function
edu_mn = edu_mean(edu_sc, ag_cntr)

# write mean edu raster
name = paste("edu_mean", "tif", sep = ".")
name = paste(path_vul, name, sep = "/")
writeRaster(edu_mn, filename = name, format = "GTiff", bylayer = T, overwrite = T)


## inc ##

# crop to country
#inc_cntr = mask(crop(inc, adm0_cntr), adm0_cntr)
inc_cntr = inc
remove(inc)

# load pop tot per sex
name = c("f_tot.tif", "m_tot.tif")
name = paste(path_vul, name, sep = "/")

pop_sex = stack(name)

# calculate income per sex
inc_fm = inc_sex(inc_cntr, pop_sex)

# write rasters
name = c("gnicf2010_gpw.tif", "gnicm2010_gpw.tif")
name = paste(path_vul, name, sep = "/")
writeRaster(inc_fm, filename = name, format = "GTiff", bylayer = T, overwrite = T)


## settlement type ##
settl_cntr = settl
remove(settl)
val = unique(settl_cntr)
# run function
settl_st = settl_cl(settl_cntr)

# write rasters
name = paste(paste("smod15", paste0("cl", val), sep = "_"), "tif", sep = ".")
name = paste(path_vul, name, sep = "/")
writeRaster(settl_st, filename = name, format = "GTiff", bylayer = T, overwrite = T)


