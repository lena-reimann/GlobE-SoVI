##########################################
#### preprocess admin unit level data ####
##########################################
# by Lena Reimann
# Oct 21, 2022

# Goal: 1) preprocess gdl data based on V6 SHDI database: fill gaps, select variables, merge with admin units, produce rasters (done in python)
#       2) preprocess national-level data needed to fill data gaps in subnational SHDI data and fill the gaps

# Comment: one way to run the entire script without relying on python would be to produce one gdl admin unit raster 
# with unique ids and use that as a basis for reclassifying the shdi values per variable

rm(list=ls())

path = "./vulnerability/global"
lib = "./r_libs"

# packages
library(sp, lib.loc = lib)
#library(raster, lib.loc = lib)
library(terra, lib.loc = lib)
library(dplyr, lib.loc = lib)
library(sf, lib.loc = lib)
library(rgdal, lib.loc = lib)
library(igraph, lib.loc = lib)
library(iterators, lib.loc = lib)

#for parallel processing
library("parallel", lib.loc = lib)
library("foreach", lib.loc = lib)
library("doParallel", lib.loc = lib)
library("igraph", lib.loc = lib)
library("snow", lib.loc = lib)
library("doSNOW", lib.loc = lib)

cores = 2

####----------- 1.) load data -------------####

# iso codes
iso = read.csv(paste(path, "gpw/gpw_national_identifier_un_wb.csv", sep = "/"))
iso = iso[, c("Value", "ISOCODE")]

# ISO code raster
gpw = rast(paste(path, "gpw/gpw_v4_nat_eucl.tif", sep = "/"))
gpw_ori = rast(paste(path, "gpw/gpw_v4_national_identifier_grid_rev11_30_sec.tif", sep = "/"))

# gdl data
gdl <- read_sf(paste(path, "gdl/shdi2022_World_large.shp", sep = "/"))
shdi <- read.csv(paste(path, "gdl/SHDI-SGDI-Total 5.0.csv", sep = "/"), stringsAsFactors = F)

gdl_tif_ori = rast(paste(path, "gdl/gdl_admin_ids.tif", sep = "/")) # load raster of GDL admin units based on unique ids
gdl_tif_euc = rast(paste(path, "gdl/gdl_admin_eucl.tif", sep = "/")) # load raster of GDL admin units based on unique ids, extended with eucl allo

# load national-level
hdi = read.csv(paste(path, "HDR21-22_edu+inc_2010.csv", sep = "/"), stringsAsFactors = F)
gov_all = read.csv(paste(path, "observed_yr.csv", sep = "/"), stringsAsFactors = F) # Andrijevic 



####------------------ 2.) Functions --------------------####

## 2010 gap filling gdl data ##
gap_fill_yrs <- function(shdi_db, shdi_yr, years, iso, var) {
  
  ## fill gaps in years
  for (j in 1:length(years)) { 
    
    sub = subset(shdi_db, year == years[j] & iso_code == iso[i])
    
    if (sum(sub[,var[1]], na.rm = T) == 0 & j != length(years)) { # if data cannot be filled
      next
    } else if (nrow(sub) == 0) {
      shdi_yr[shdi_yr$iso_code == iso[i], var[3]] <- 0
      
      print(paste0("no replacement year"))
      break
    } else if (sum(sub[,var[1]], na.rm = T) == 0 & j == length(years)) {
      shdi_yr[shdi_yr$iso_code == iso[i], var[3]] <- rep(0, nrow(sub))
      
      print(paste0("no replacement year"))
      break
    } else {
      shdi_yr[shdi_yr$iso_code == iso[i], var[1]] <- sub[, var[1]]
      shdi_yr[shdi_yr$iso_code == iso[i], var[2]] <- sub[, var[2]]
      shdi_yr[shdi_yr$iso_code == iso[i], var[3]] <- rep(years[j], nrow(sub))
      
      print(paste0("taken ", years[j]))
      break
    }
  } # end of j loop

    return(shdi_yr)
}

## fill gaps in shdi data with national-level data from hdi and harmonize extents
fill_gap_r = function(subn_ori, nat, subn_euc) { #subnational raster (original and eucl allo) & national raster
  require(terra)
  
  # fill gaps with national-level data
  subn_gfil = merge(subn_ori, nat) # use as mask
  subn_gfil_eu = merge(subn_euc, nat)
  
  # mask gap-filled raster based on eucl with gap-filled raster without eucl
  subn_gfil = mask(subn_gfil_eu, subn_gfil)
  
  return(subn_gfil)
}


####------------------ 3.) Analysis --------------------####

#### 1. process GDL data ####

## Step 1a. gap filling ##

# select 2010
shdi2010 = shdi[which(shdi$year == 2010),]

# prep gap filling
iso_gdl = unique(shdi2010$iso_code)
yrs = c(rbind(seq(2011, 2019), seq(2009, 2001)))

# baseline year (will be replaced when using different year)
shdi2010$yr_edu <- 2010
shdi2010$yr_inc <- 2010

# define variables (column names)
edu = c("mschf", "mschm", "yr_edu")
inc = c("gnicf", "gnicm", "yr_inc")

# fill gaps
for (i in 1:length(iso_gdl)) { # check each cntr individually
  sub = subset(shdi2010, shdi2010$iso_code == iso_gdl[i])
  
  #education
  if (sum(sub[,edu[1]], na.rm = T) == 0) { # if 2010 data missing
    shdi2010 <- gap_fill_yrs(shdi, shdi2010, yrs, iso_gdl, edu)
  } else {
    print(paste0("2010 available"))
  }
  
  #income
  if (sum(sub[,inc[1]], na.rm = T) == 0) { # if 2010 data missing
    shdi2010 <- gap_fill_yrs(shdi, shdi2010, yrs, iso_gdl, inc)
  } else {
    print(paste0("2010 available"))
  }
} # end of i loop

# select columns
shdi2010 = data.frame(shdi2010[,4:5], shdi2010$continent, 
                      round(shdi2010$msch * 10, 0), round(shdi2010$mschf * 10, 0), round(shdi2010$mschm * 10,0), # multiply by 10 and round values to be able to create integer values in the rasters
                      round(shdi2010$gnic, 0), round(shdi2010$gnicf, 0), round(shdi2010$gnicm, 0), 
                      shdi2010$yr_edu, shdi2010$yr_inc)
colnames(shdi2010) = c("gdlcode", "level", "continent",
                       "msch", "mschf", "mschm", "gnic", "gnicf", "gnicm", "yr_edu", "yr_inc")



## Step 1b. replace NA for female versus male with national totals ##

shdi2010$edu_sex <- 1
shdi2010$inc_sex <- 1

# education
shdi2010[which(is.na(shdi2010$mschf)), "edu_sex"] <- 0
shdi2010[which(is.na(shdi2010$mschf)), "mschf"] <- shdi2010[which(is.na(shdi2010$mschf)), "msch"] 
shdi2010[which(is.na(shdi2010$mschm)), "mschm"] <- shdi2010[which(is.na(shdi2010$mschm)), "msch"] 

# income
shdi2010[which(is.na(shdi2010$gnicf)), "inc_sex"] <- 0
shdi2010[which(is.na(shdi2010$gnicf)), "gnicf"] <- shdi2010[which(is.na(shdi2010$gnicf)), "gnic"] 
shdi2010[which(is.na(shdi2010$gnicm)), "gnicm"] <- shdi2010[which(is.na(shdi2010$gnicm)), "gnic"] 

#write.csv(shdi2010, "Population/global_data_lab/SHDI-SGDI-2010.csv", row.names = F)
#shdi2010 = read.csv("Population/global_data_lab/SHDI-SGDI-2010.csv")


## Step 1c. merge with admin units ##
# add shape ids
id = seq(0, nrow(gdl)-1)
gdl = cbind(id, gdl)

# merge to polygons
glo_sf = merge(gdl[,-c(3)], shdi2010, by = ("gdlcode"), all = T)

# check for duplicates (due to all = T) and remove
dupl = glo_sf[which(glo_sf$gdlcode %in% gdl$gdlcode == F), ] #don't contain any geometries

glo_sf = glo_sf[-c(which(glo_sf$gdlcode %in% gdl$gdlcode == F)), ] #all geometries 

#write_sf(glo_sf, "Population/global_data_lab/SHDI-SGDI-2010.shp")
#glo_sf = read_sf("Population/global_data_lab/SHDI-SGDI-2010.shp")


# aggregate for national values ## revise - reason for not just selected national totals? ## 
#glo_nat = aggregate(x = glo_sf, 
#                    by = list(glo_sf$iso_code), #[("index")],
#                    FUN = mean)
#write_sf(glo_nat[,c(1,5:ncol(glo_nat))], "Population/global_data_lab/SHDI-SGDI-2010_nat.shp")


## Step 1d. rasterize admin units with the variables of interest ##

# extract relevant variables for reclassification (make sure that ids are aligned correctly)
vars = c("mschf", "mschm", "gnicf", "gnicm") #define variables
#vars = c("gnicm") #for testing

shdi = glo_sf[, c("id", vars)]
st_geometry(shdi) <- NULL

for (i in 1:length(vars)) {
  rcl = c(shdi$id, shdi[,vars[i]]) # determine reclass values
  rcl = matrix(rcl, ncol = 2) # make matrix
  
  name = paste0(path, "/gdl/", vars[i], "2010.tif")
  classify(gdl_tif_ori, rcl, filename = name, overwrite = T, datatype = "INT4U")
  
  name = paste0(path, "/gdl/", vars[i], "2010_eucl.tif")
  classify(gdl_tif_euc, rcl, filename = name, overwrite = T, datatype = "INT4U")
}



#### 2. process HDI data ####

# merge HDI with isos
#hdi = merge(hdi, iso, by.x = "iso3", by.y = "ISOCODE", all = T)
hdi = merge(iso, hdi, by.x = "ISOCODE", by.y = "iso3", all.x = T)

# replace values in GRL and PRI with values from DK and US
#hdi[which(hdi$ISOCODE == "GRL"), 4:11] <- hdi[which(hdi$ISOCODE == "DNK"), 4:11]
#hdi[which(hdi$ISOCODE == "PRI"), 4:11] <- hdi[which(hdi$ISOCODE == "USA"), 4:11]

# fill gaps in gender-based edu and inc
hdi$edu_sex <- 1
hdi$inc_sex <- 1

#education
hdi[which(is.na(hdi$mschf)), "edu_sex"] <- 0
hdi[which(is.na(hdi$mschf)), "mschf"] <- hdi[which(is.na(hdi$mschf)), "msch"] 
hdi[which(is.na(hdi$mschm)), "mschm"] <- hdi[which(is.na(hdi$mschm)), "msch"] 

#income
hdi[which(is.na(hdi$gnicf)), "inc_sex"] <- 0
hdi[which(is.na(hdi$gnicf)), "gnicf"] <- hdi[which(is.na(hdi$gnicf)), "gnic"] 
hdi[which(is.na(hdi$gnicm)), "gnicm"] <- hdi[which(is.na(hdi$gnicm)), "gnic"] 


# multiply edu*10 and round
hdi$msch = round(hdi$msch * 10, 0)
hdi$mschf = round(hdi$mschf * 10, 0)
hdi$mschm = round(hdi$mschm * 10, 0)

# round inc
hdi$gnic = round(hdi$gnic, 0)
hdi$gnicf = round(hdi$gnicf, 0)
hdi$gnicm = round(hdi$gnicm, 0)

# reclassify gpw_ori raster
vars = c("mschf", "mschm", "gnicf", "gnicm") #define variables

for (i in 1:length(vars)) {
  rcl = c(hdi$Value, hdi[,vars[i]]) # determine reclass values
  rcl = matrix(rcl, ncol = 2) # make matrix
  
  name = paste0(path, "/hdi_", vars[i], "2010.tif")
  classify(gpw_ori, rcl, filename = name, overwrite = T, datatype = "INT4U")
}



#### 3. fill GDL gaps ####

cl <- makeCluster(cores)
registerDoParallel(cl, cores = cores)

# fill gaps with HDI data
foreach (i = 1:length(vars),
         .packages = c("terra", "igraph", "dplyr")) %dopar%  {
           # load shdi rasters         
           name = paste0(path, "/gdl/", vars[i], "2010.tif")
           shdi = rast(name)
           
           name = paste0(path, "/gdl/", vars[i], "2010_eucl.tif")
           shdi_euc = rast(name)
           
           # load HDI rasters
           name = paste0(path, "/hdi_", vars[i], "2010.tif")
           hdi = rast(name)
           
           # fill shdi gaps
           shdi_gfil = fill_gap_r(shdi, hdi, shdi_euc)
           
           # path and vars need to be included here because otherwise not handed to each worker
           #path = "./vulnerability/global"
           #vars = c("mschf", "mschm", "gnicf", "gnicm")
           name = paste0(path, "/gdl/", vars[i], "2010_gfil_EK.tif")
           writeRaster(shdi_gfil, filename = name, overwrite = T, datatype = "INT4U")
         }
stopCluster(cl)  

  
#### 4. process governance data ####

gov_all$countrycode = recode(gov_all$countrycode, 
                             "ADO" = "AND",
                             "KSV" = "KOS",
                             "TMP" = "TLS",
                             "WBG" = "PSE")

# subset to 2010
gov = subset(gov_all, year == 2010, select = c(countrycode,governance))
gov = merge(gov, iso, by.x = "countrycode", by.y = "ISOCODE", all = T)

# use ANT gov for BES, CUW, SXM
gov[gov$countrycode == "BES", "governance"] <- gov[gov$countrycode == "ANT", "governance"] 
gov[gov$countrycode == "CUW", "governance"] <- gov[gov$countrycode == "ANT", "governance"] 
gov[gov$countrycode == "SXM", "governance"] <- gov[gov$countrycode == "ANT", "governance"] 

# remove ANT
gov = subset(gov, countrycode != "ANT")

### replace NAs with other values??
# e.g. SJM - Norway; SMR - ITA; MCO - FRA? ### BUT: would have to be done for each country!


# fill in NAs from other time steps
gov2011 = subset(gov_all, year == 2011, select = c(countrycode, governance))
gov[which(gov$countrycode == "JEY"), "governance"] <- gov2011[which(gov2011$countrycode == "JEY"), "governance"]
gov[which(gov$countrycode == "SSD"), "governance"] <- gov2011[which(gov2011$countrycode == "SSD"), "governance"]

# replace NAs by -999
#gov = gov %>%
#  mutate(governance = replace(governance, is.na(governance), -999 ))


#reclass gpw raster
rcl = c(gov$Value, gov$governance) # determine reclass values
rcl = matrix(rcl, ncol = 2) # make matrix
#rcl = rcl[-c(which(is.na(rcl[,1]))), ] # remove NAs

name = paste(path, "gov2010.tif", sep = "/")
classify(gpw_ori, rcl, filename = name, overwrite = T)
#gov = rast(name)





#### Data correlations for edu and inc ####

# establish correlation between edu and inc
#cor(shdi2010$msch, log(shdi2010$gnic), method = c("pearson"))

# Creating the plot
#png(file = "Population/global_data_lab/cor_edu_inc_log_Europe.png")
#yrs_school = shdi2010[which(shdi2010$continent == "Europe"), "msch"]
#GNI_cap_log = log(shdi2010[which(shdi2010$continent == "Europe"), "gnic"] / 1000)
#plot(yrs_school, GNI_cap_log, pch = 19, col = "lightblue")

# Regression line
#abline(lm(GNI_cap_log ~ yrs_school), col = "red", lwd = 3)

# Pearson correlation
#text(paste("Correlation: ", round(cor(yrs_school, GNI_cap_log), 2)), x = 8, y = 5)

#dev.off()


# Creating the plot entire globe
#png(file = "Population/global_data_lab/cor_edu_inc_log.png")
#yrs_school = shdi2010[, "msch"]
#GNI_cap_log = log(shdi2010[, "gnic"] / 1000)
#plot(yrs_school, GNI_cap_log, pch = 19, col = "lightblue")

# Regression line
#abline(lm(GNI_cap_log ~ yrs_school), col = "red", lwd = 3)

# Pearson correlation
#text(paste("Correlation: ", round(cor(yrs_school, GNI_cap_log), 2)), x = 2, y = 5)

#dev.off()



