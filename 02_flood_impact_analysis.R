##########################################################################################
#### Estimate regression model based on vulnerability characteristics per flood event ####
##########################################################################################
# by Lena Reimann
# October 31, 2023

# This script corresponds to the data processing described in '2) Flood Impact Analysis' of
# Reimann et al. "An empirical social vulnerability map for flood risk assessment at global scale (‘GlobE-SoVI’)".
# After calculating vulnerability characteristics per flood event, different regression model are tested before the final model 
# model configuration is found along with the vulnerability variables that significantly contribute to predicting flood fatalities.
# Technical note: the part 'I - establish vulnerability characteristics per flood event' was run with 8 cores and 255G of memory

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
library(corrplot, lib.loc = lib)
library(stats, lib.loc = lib)
library(car, lib.loc = lib)


#for parallel processing
library("parallel", lib.loc = lib)
library("foreach", lib.loc = lib)
library("doParallel", lib.loc = lib)
library("igraph", lib.loc = lib)
library("snow", lib.loc = lib)
library("doSNOW", lib.loc = lib)

cores <- 8


####----------- 1.) Load data -------------####

## exposure ##
# load pop raster for harmonization
file = c("GHS_POP_E2000_GLOBE_R2019A_4326_30ss_V1_0.tif", 
         "GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif")  
file = paste0(path, file)
pop = stack(file)


## hazard ##
# load events from Global Floods Database (GFD) (https://global-flood-database.cloudtostreet.ai/)
ev = read.csv(paste0(path, "gfd_event_stats_20215_13_iso.csv"), encoding = "UTF-8") # flood events with iso country codes added as a column "iso"
 
# extract event id
id = ev$system.index


## vulnerability ##
# age (produced in 01_vulnerability_drivers)
ag1 <- sprintf("ag%s", seq(0,60, by= 5))
ag1[length(ag1)+1] <- paste0(sprintf("ag%s", length(ag1) * 5), "+")
ag1 = c(paste0(ag1, "f"), paste0(ag1, "m"))

name = paste(ag1, "prop.tif", sep = "_")
name = paste(path, name, sep = "/")
prop_ag = stack(name)

# education (produced in 01_vulnerability_drivers)
name = paste(paste("edu", ag1, sep = "_"), "tif", sep = ".")
name = paste(path, name, sep = "/")  
mys = stack(name)

# income (produced in 01_vulnerability_drivers)
name = c("gnicf2010.tif", "gnicm2010.tif")
name = paste(path, name, sep = "/")
inc = stack(name)

# settlements
# GHS-SMOD (https://ghsl.jrc.ec.europa.eu/ghs_smod2019.php)
name = c("GHS_SMOD2000_WGS.tif", "GHS_SMOD2015_WGS.tif") #smod settlement types of 2000 and 2015 projected to WGS coordinates
name = paste(path, name, sep = "/")
settl = raster(name)

# stack all v 
vul = stack(prop_ag, mys, inc, settl)

# other potential vul indicators (added later)
# travel to time health facility (walking & motorized) (https://data.malariaatlas.org/maps)
name = c("2020_walking_only_travel_time_to_healthcare.geotiff", "2020_motorized_travel_time_to_healthcare.geotiff")
name = paste(path, name, sep = "/")  
tvl = stack(name)

# protection standards coastal and riverine (https://nhess.copernicus.org/articles/16/1049/2016/nhess-16-1049-2016-supplement.zip)
name = c("Flopros_prot_stand_cst.tif", "Flopros_prot_stand_riv.tif") # coastal and riverine flood protection standards per administrative unit converted into a global raster of 30 arc seconds spatial resolution
name = paste(path, name, sep = "/")  
prot = stack(name)

# poverty levels (https://datacatalog.worldbank.org/search/dataset/0042041)
name = c("GSAP2_tot_190.tif", "GSAP2_tot_320.tif", "GSAP2_tot_550.tif", "GSAP2_pop_tot.tif", "GSAP2_cns_d_md.tif", "GSAP2_cns_d_mn.tif") # poverty indicators per administrative unit converted into a global raster of 30 arc seconds spatial resolution
name = paste(path, name, sep = "/")  
pov = stack(name)

# national-level datasets
#governance index (https://github.com/marina-andrijevic/governance2019)
gov = read.csv(paste(path, "governance2019-master/data/observed_yr.csv", sep = "/"))

#adaptation readiness index (https://gain.nd.edu/our-work/country-index/download-data/)
nd.gain = read.csv(paste(path, "nd-gain/resources/readiness/readiness.csv", sep = "/"))
  

####------------------ 2.) Functions --------------------####

## hazard ##
# aggregate flood plains to 30 arc seconds
event = function(id, pop) { # input: vector with event ids; 
  #stack of population data with the required time steps; 
  require(raster)
  
  # unzip data 
  test = dir.exists(paste(path, id, sep = "/"))
  
  if (test == F) { # only if not done yet
    unzip(paste0(path, id, ".zip"), exdir = paste(path, id, sep = "/"))
  }
  
  # check whether 'original' file exists (some floodplains are so large that the data come in several tiffs)
  name = paste0(path, id, "/", id, ".tif")
  test = file.exists(name)
  
  # if not: mosaic existing files
  if (test == F) {
    files = list.files(paste0(path, id, "/"), pattern = "\\.tif$")
    
    # make list
    flo = list()
    dur = list()
    for (i in 1:length(files)) {
      flo[[i]] = raster(paste0(path, id, "/", files[i]), band = 1)
      dur[[i]] = raster(paste0(path, id, "/", files[i]), band = 2)
    }
    
    # make mosaic (but dont write)
    flo$fun <- max
    dur$fun <- max
    
    flo$na.rm <- T
    dur$na.rm <- T
    
    flo = do.call(mosaic, flo)
    dur = do.call(mosaic, dur)

  } else { # load rasters
    
    flo = raster(paste0(path, id, "/", id, ".tif"), band = 1)
    dur = raster(paste0(path, id, "/", id, ".tif"), band = 2)
  }
  
  # determine extent
  e = extent(flo)
  
  # crop pop to flood extent
  pop_fp = crop(pop, e)

  # check whether aggregated files exist (here: duration)
  name = paste0(path, "/", id, "_dur.tif")
  test = file.exists(name)
  
  if (test == F) { 
    # aggregate flooded land to resolution of pop
    flo4 = aggregate(flo, fact = 4, fun = mean)
    dur4 = aggregate(dur, fact = 4, fun = mean)
    
    # resample to pop
    name = paste0(path, "/", id, "_flo.tif")
    flo4 = resample(flo4, pop_fp, method = "ngb", filename = name, overwrite = T)

    name = paste0(path, "/", id, "_dur.tif")
    dur4 = resample(dur4, pop_fp, method = "ngb", filename = name, overwrite = T)
  }
}

## Exposure & vulnerability ##
# produce df of flood event with E and V characteristics per cell
vul_char_rast = function(id, ag, pop, vul) { # input: vector with event ids; 
  #vector with AGs by f/m
  #stack of population data with the required time steps; 
  #stack of vul characteristics
  require(raster)
  
  ## 1. conversion of vul char from raster to table ##
  # load aggregated floodplain rasters
  name = paste0(path, "/", id, "_flo.tif")
  flo4 = raster(name)
  
  name = paste0(path, "/", id, "_dur.tif")
  dur4 = raster(name)
  
  # determine extent
  e = extent(flo4)  
  
  # crop pop and vul to flood extent
  pop_fp = crop(pop, e)
  vul_fp = crop(vul, e)
  
  # set to NAs
  flo4 = calc(flo4, fun = function(x) {x[x == 0] <- NA; return(x)})
  dur4 = calc(dur4, fun = function(x) {x[x == 0] <- NA; return(x)})
  
  # determine area flooded
  flo4_sqkm = overlay(area(flo4), flo4, fun = function(x,y) { x * y })
  
  # convert raster to points
  flo4_df = as.data.frame(rasterToPoints(flo4))
  sqkm_df = as.data.frame(rasterToPoints(flo4_sqkm))
  dur4_df = as.data.frame(rasterToPoints(dur4))
  
  # merge both df
  ev_df = left_join(flo4_df, sqkm_df, by = c("x", "y"))
  ev_df = left_join(ev_df, dur4_df, by = c("x", "y"))
  
  # change colnames
  colnames(ev_df) = c("x", "y", "flo", "sqk", "dur")
  
  # round values
  ev_df$x <- round(ev_df$x, digits = 7)
  ev_df$y <- round(ev_df$y, digits = 7)
  
  
  ## 2. add population to df (not exposed pop! calculate E later) ##
  
  # select pop and smod year based on event year
  ev_yr = as.numeric(substr(strsplit(id, "_")[[1]][6], 1, 4))
  
  if(ev_yr < 2008) {
    pop_fp = pop_fp[[1]]
    vul_fp = vul_fp[[1:(nlayers(vul)-1)]]
  } else {
    pop_fp = pop_fp[[2]]
    vul_fp = vul_fp[[-c((nlayers(vul)-1))]]
  }
  
  # set to NAs and convert to points (reduces size)
  pop_fp_df = as.data.frame(rasterToPoints(overlay(pop_fp, flo4, fun = function(x, y) {x[y == 0] <- NA; return(x)})))
  
  # change colnames
  colnames(pop_fp_df) = c("x", "y", "pop")
  
  # round values
  pop_fp_df$x <- round(pop_fp_df$x, digits = 7)
  pop_fp_df$y <- round(pop_fp_df$y, digits = 7)
  
  # merge hazard and exposure by coord
  ev_df = left_join(ev_df, pop_fp_df, by = c("x", "y"))
  
  
  ## 3. add vulnerability characteristics to df ##

  # set to NAs and convert to points
  vul_fp_df = overlay(vul_fp, flo4, fun = function(x, y) {x[y == 0] <- NA; return(x)})
  vul_fp_df = rasterToPoints(vul_fp_df)
  vul_fp_df = as.data.frame(vul_fp_df)
  
  # change colnames
  colnames(vul_fp_df) = c("x", "y", ag1, paste("edu", ag1, sep = "_"), "gnicf", "gnicm", "smod")
  
  # round values
  vul_fp_df$x <- round(vul_fp_df$x, digits = 7)
  vul_fp_df$y <- round(vul_fp_df$y, digits = 7)
  
  # merge hazard and vulnerability data by coord
  ev_df = left_join(ev_df, vul_fp_df, by = c("x", "y"))
  
  # write basic vulnerability characteristics df
  name = paste0(path, id, "_vul_char.csv")
  write.csv(ev_df, name, row.names = F)
}

# produce df of flood event with E and V characteristics per event
vul_char_event = function(ev_df, ag, id) { # input: table with vulnerability characteristics per raster cell;
  #vector with AGs by f/m
  #id of the event analyzed (for writing the files)
  require(raster)
  
  ## 1. calculate event char per cell and event ##
  
  # a. event characteristic totals
  flo_tot = sum(ev_df$flo, na.rm = T)
  
  # duration
  dur_mn = mean(ev_df$dur, na.rm = T)
  dur_md = median(ev_df$dur, na.rm = T)
  
  # b. pop exposed
  pop_exp = ev_df$flo * ev_df$pop
  pop_day = pop_exp * ev_df$dur
  
  exp_tot = sum(pop_exp, na.rm = T)
  pop_day_tot = exp_tot * dur_mn
  
  ## 2. age characteristics ##
  
  # a. broad AG proportions & dependency ratios
  
  #children f/m/t
  ags = ag[1:3]
  ags = which(colnames(ev_df) %in% ags) #select the respective AGs
  chi_f = rowSums(ev_df[,ags], na.rm = T) * pop_exp
  
  ags = ag[(length(ag)/2+1) : (length(ag)/2+3)]
  ags = which(colnames(ev_df) %in% ags) 
  chi_m = rowSums(ev_df[,ags], na.rm = T) * pop_exp
  
  chi_t = chi_f + chi_m
  
  #working f/m/t
  ags = c(ag[4 : (length(ag)/2-1)])
  ags = which(colnames(ev_df) %in% ags)
  wrk_f = rowSums(ev_df[,ags], na.rm = T) * pop_exp
  
  ags = ag[(length(ag)/2+4) : (length(ag)-1)]
  ags = which(colnames(ev_df) %in% ags)
  wrk_m = rowSums(ev_df[,ags], na.rm = T) * pop_exp
  
  wrk_t = wrk_f + wrk_m
  
  #elderly
  ags = c(ag[(length(ag)/2)])
  ags = which(colnames(ev_df) %in% ags)
  eld_f = ev_df[,ags] * pop_exp # different approach because only one AG
  
  ags = ag[length(ag)]
  ags = which(colnames(ev_df) %in% ags)
  eld_m = ev_df[,ags] * pop_exp
  
  eld_t = eld_f + eld_m
  
  #broad AG totals
  chi_f_tot = sum(chi_f, na.rm = T)
  wrk_f_tot = sum(wrk_f, na.rm = T)
  eld_f_tot = sum(eld_f, na.rm = T)
  
  chi_m_tot = sum(chi_m, na.rm = T)
  wrk_m_tot = sum(wrk_m, na.rm = T)
  eld_m_tot = sum(eld_m, na.rm = T)
  
  chi_tot = sum(chi_t, na.rm = T)
  wrk_tot = sum(wrk_t, na.rm = T)
  eld_tot = sum(eld_t, na.rm = T)
  
  #dependency ratios
  ydr = chi_t / wrk_t
  odr = eld_t / wrk_t
  cdr = (chi_t + eld_t) / wrk_t
  
  ydr_tot = chi_tot / wrk_tot
  odr_tot = eld_tot / wrk_tot
  cdr_tot = (chi_tot + eld_tot) / wrk_tot
  
  
  # b. mean and median age per event
  ags = which(colnames(ev_df) %in% ag)
  
  #tot pop by ag
  ags_tot = colSums(pop_exp * ev_df[,ags], na.rm = T)
  ags_tot = ags_tot[1 : (length(ags_tot) / 2)] + ags_tot[((length(ags_tot) / 2)+1) : length(ags_tot)]

  #weights on the AGs
  ags_wgt = seq(2, 67, by = 5)
  ags_wgt[length(ags_wgt)] <- 72 # replace last AG weight because larger AG than the rest
  
  #mean
  mean_age = sum(ags_tot * ags_wgt) / exp_tot
  
  # median 
  #define 'middle person'
  per = round(sum(ags_tot) / 2)
  
  #calc cumulative pop across AGs
  cum_ag = cumsum(ags_tot)
  
  #numeric AGs vector
  ags_num = seq(0, 65, by = 5)
  
  #establish position of median AG and actual AG
  pos = which(cum_ag > per)[1]
  med_ag = ags_num[pos]
  
  # account for 0s and NAs
  med_age = ifelse(med_ag == 0, med_ag + cum_ag[pos+1] * 5, med_ag + (per - cum_ag[pos-1]) / (cum_ag[pos] - cum_ag[pos-1]) * 5)
  med_age = ifelse(is.na(med_ag), NA, med_age)
  
  
  # c. female proportion + total 
  f_prop = rowSums(ev_df[,ags[1:(length(ags)/2)]], na.rm = T)
  f_exp = f_prop * pop_exp
  
  f_exp_tot = sum(f_exp, na.rm = T)
  f_prop_tot = f_exp_tot / exp_tot
  
  
  ## 3. education: mean MYS ##
  edu = which(colnames(ev_df) %in% paste("edu", ag, sep = "_")) #establish edu cols
  # mean years exposed weighted by each AG individually
  mys_exp = rowSums(ev_df[, ags] * pop_exp * ev_df[, edu], na.rm = T) / pop_exp
  # weighted by gender
  mys_gen = (f_prop * ev_df$edu_ag25f) + ((1 - f_prop) * ev_df$edu_ag25m)
  
  # MYS of exposed pop for entire event (weighted by AGs)
  mys_exp_tot = sum(colSums(pop_exp * ev_df[,ags], na.rm = T) * colMeans(ev_df[, edu], na.rm = T)) / exp_tot
  mys_gen_tot = mean(mys_gen, na.rm = T)
  
  ## 4. exposed income ##
  incf_exp = ev_df$gnicf * f_exp
  incm_exp = ev_df$gnicm * (pop_exp - f_exp)
  inct_exp = incf_exp + incm_exp # total exposed income
  
  inc_gap = (ev_df$gnicm - ev_df$gnicf) / ev_df$gnicm
  inc_gap_exp = (incm_exp - incf_exp) / incm_exp
  minc_exp = inct_exp / pop_exp # mean exposed income
  
  #event totals
  incf_exp_tot = sum(incf_exp, na.rm = T)
  incm_exp_tot = sum(incm_exp, na.rm = T)
  inct_exp_tot = incf_exp_tot + incm_exp_tot
  
  inc_gap_tot = (sum(ev_df$gnicm, na.rm = T) - sum(ev_df$gnicf, na.rm = T)) / sum(ev_df$gnicm, na.rm = T)
  inc_gap_exp_tot = (incm_exp_tot - incf_exp_tot) / incm_exp_tot
  
  minc_exp_tot = inct_exp_tot / exp_tot # mean total exposed income
  
  ## 5. exposed settlements relative to all non-water flooded cells ##
  smod11_rel = sum(ev_df$smod == 11) / sum(ev_df$smod != 10)
  smod12_rel = sum(ev_df$smod == 12) / sum(ev_df$smod != 10)
  smod13_rel = sum(ev_df$smod == 13) / sum(ev_df$smod != 10)
  smod21_rel = sum(ev_df$smod == 21) / sum(ev_df$smod != 10)
  smod22_rel = sum(ev_df$smod == 22) / sum(ev_df$smod != 10)
  smod23_rel = sum(ev_df$smod == 23) / sum(ev_df$smod != 10)
  smod30_rel = sum(ev_df$smod == 30) / sum(ev_df$smod != 10)
  
  rur_rel = (sum(ev_df$smod == 11) + sum(ev_df$smod == 12) + sum(ev_df$smod == 13)) / sum(ev_df$smod != 10)
  urb_rel = (sum(ev_df$smod == 22) + sum(ev_df$smod == 23) + sum(ev_df$smod == 30)) / sum(ev_df$smod != 10)
  sub_rel = sum(ev_df$smod == 21) / sum(ev_df$smod != 10)
  
  rur_rel+urb_rel+sub_rel
  
  ## 6. assemble df of cell-based and event-based characteristics ##
  
  # create new df and retain hazard info
  ev_char = data.frame(ev_df[,c("x", "y", "flo", "sqk", "dur")], 
                       pop_exp, pop_day, 
                       chi_f, chi_m, chi_t, wrk_f, wrk_m, wrk_t, eld_f, eld_m, eld_t, ydr, odr, cdr,
                       f_prop, f_exp,
                       mys_exp, mys_gen,
                       incf_exp, incm_exp, inct_exp, inc_gap, inc_gap_exp, minc_exp
  )
  
  # write exposed vulnerability characteristics df
  name_df = paste0(path, id, "_vul_exp.csv")
  write.csv(ev_char, name_df, row.names = F)
  
  # create total event vector
  return(c(flo_tot, dur_mn, 
           exp_tot, pop_day_tot, 
           chi_f_tot, chi_m_tot, chi_tot, wrk_f_tot, wrk_m_tot, wrk_tot, eld_f_tot, eld_m_tot, eld_tot, ydr_tot, odr_tot, cdr_tot,
           f_prop_tot, f_exp_tot,
           mys_exp_tot, mys_gen_tot,
           incf_exp_tot, incm_exp_tot, inct_exp_tot, inc_gap_tot, inc_gap_exp_tot, minc_exp_tot,
           smod11_rel, smod12_rel, smod13_rel, smod21_rel, smod22_rel, smod23_rel, smod30_rel, rur_rel, urb_rel, sub_rel,
           mean_age, med_age
  ))
}

# produce df of flood event with additional vulnerability characteristics added afterwards (can be done with any additional V data)
attr = function(id, tvl, prot, pov) { # input: vector with event ids; 
  #stack of travel times; flopros protection standards; poverty data
  require(raster)
  
  ## 1. conversion from raster to table ##

  # load aggregated floodplain rasters
  name = paste0(path, "/", id, "_flo.tif")
  flo4 = raster(name)
  
  max = max(values(flo4), na.rm = T)
  
  if (max != 0) {
    # determine extent
    e = extent(flo4)  
    
    # set 0s to NA
    flo4 = calc(flo4, fun = function(x) {x[x == 0] <- NA; return(x)})
    # set all cells with some flooding to 1
    flo4 = calc(flo4, fun = function(x) {x[x > 0] <- 1; return(x)})
    
    # crop vul data to flood extent
    tvl_fp = crop(tvl, e)
    prot_fp = crop(prot, e)
    pov_fp = crop(pov, e)
    
    # a. deal with pop-independent data first
    
    # tvl time
    zon.mn = zonal(tvl_fp, flo4, fun = 'mean')
    zon.ma = zonal(tvl_fp, flo4, fun = 'max')
    
    wlk.mn = zon.mn[,2]
    wlk.ma = zon.ma[,2]  
    
    mot.mn = zon.mn[,3]
    mot.ma = zon.ma[,3]
    
    # protection levels
    zon.mn = zonal(prot_fp, flo4, fun = 'mean')
    
    prot.cst.mn = zon.mn[,2]
    prot.riv.mn = zon.mn[,3]
    prot.tot.mn = sum(prot.cst.mn, prot.riv.mn) / 2
    
    # b. deal with pop-dependent data
    
    # poverty
    #use only data located in the fp
    pov_fp = overlay(pov_fp, flo4, fun = function(x, y) {x * y})
    
    #count cells per value in layer with pop_tot
    num = freq(pov_fp[[4]])
    
    #calculate new values for pop_tot
    num[,2] = num[,1] / num[,2]
    
    # use as reclassification matrix and write classified raster into stack
    pov_fp[[nlayers(pov_fp)]] = reclassify(pov_fp[[nlayers(pov_fp)]], rcl = num)
    
    # zonals
    zon.su = zonal(pov_fp, flo4, fun = 'sum')
    zon.md = zonal(pov_fp, flo4, fun = 'median')
    zon.mn = zonal(pov_fp, flo4, fun = 'mean')
    
    # extract relevant variables
    rel.190 = zon.su[,2] / zon.su[,5]
    rel.320 = zon.su[,3] / zon.su[,5]
    rel.550 = zon.su[,4] / zon.su[,5]
    
    cns_md = zon.md[,6] 
    cns_mn = zon.mn[,7]
    
    # combine all vars in one vector
    c(wlk.mn, wlk.ma, mot.mn, mot.ma, 
      prot.cst.mn, prot.riv.mn, prot.tot.mn,
      rel.190, rel.320, rel.550, cns_md, cns_mn
    )
  } else { # if no flooding in the floodplain data (did not show in the previous analysis because we would convert the raster cells to points, which would then be 0s)
    rep(0, 12) #12 = number of variables calculated
  }
}


####------------------ 3.) Analysis --------------------####

#### I - establish vulnerability characteristics per flood event ####

## 1. loop through events ##
cl <- makeCluster(cores)
registerDoParallel(cl, cores = cores)

t_start <- Sys.time()

sum_ev = 
foreach(i = 1:length(id), .combine = 'rbind',
          .packages = c("raster", "igraph", "dplyr")) %dopar% {
            
  # a. run raster aggregation function
  event(id[i], pop)
            
  # b. run function for producing a table with V char per raster cell
  vul_char_rast(id[i], ag1, pop, vul)   
  
  # c. run function for producing a table with exposed V char per event
  #load df
  name = paste0(path, id[i], "_vul_char.csv")
  ev_df = read.csv(name)
  
  #establish column names
  ag1 = colnames(ev_df %>% dplyr:: select(starts_with("ag")))

  #run function
  ev_char = vul_char_event(ev_df, ag1, id[i])
  
  #return vector per event
  c(id[i], ev_char)
}

stopCluster(cl) 
t_stop <- Sys.time()
print(t_stop - t_start)


## 2. make a df from the matrix ##

#load one event df and use colnames
cols = colnames(read.csv(paste0(path, id[1], "_vul_exp.csv"), encoding = "UTF-8"))

#add additional colnames
cols = c("id" , "flo", "dur_mean",
         cols[6 : length(cols)], 
         "smod11_rel", "smod12_rel", "smod13_rel", "smod21_rel", "smod22_rel", "smod23_rel", "smod30_rel", 
         "rur_rel", "urb_rel", "sub_rel",
         "mean_age", "med_age")

sum_ev = as.data.frame(sum_ev)
colnames(sum_ev) <- cols

## 3. merge with selected event characteristics ##

#determine start and end date
start = as.Date(substr(id, 15, 23), "%Y%m%d")
n = 7
end = as.Date(substr(id, nchar(id)-n, nchar(id)), "%Y%m%d")

#select characteristics
ev_sel = data.frame(ev$system.index, ev$iso, ev$dfo_dead)

# add total duration in days
ev_sel$dur_tot = as.integer(round(difftime(end, start, units = "days")) +1) # inclusive definition

colnames(ev_sel) = c("system.index", "iso", "dfo_dead", "dur_tot")

# merge with events
ev_char = merge(ev_sel, sum_ev, by.x = "system.index", by.y = "id")


## 4. run function for adding the other potential vul characteristics ##
cl <- makeCluster(cores)
registerDoParallel(cl, cores = cores)

sum_ev = 
  foreach(i = 1:length(id), .combine = 'rbind',
          .packages = c("raster", "igraph", "dplyr")) %dopar% {

            attr.ev = attr(id[i], tvl, prot, pov)
            
            #return vector per event
            c(id[i], attr.ev)
          }

stopCluster(cl) 

# make a df from the matrix
# colnames
cols = c("id", "wlk.mn", "wlk.ma", "mot.mn", "mot.ma",
         "prot.cst.mn", "prot.riv.mn", "prot.tot.mn",
         "rel.190", "rel.320", "rel.550", "cns_md", "cns_mn"
)

sum_ev = as.data.frame(sum_ev)
colnames(sum_ev) <- cols

# merge with events
ev.tot = merge(ev_char, sum_ev, by.x = "system.index", by.y = "id")

# write csv
write.csv(ev.tot, paste0(path, "event_vul_char.csv"), row.names = F)


#### II - regression analysis ####

## 1. prep data ##

# extract event year
ev$year = extr_year(ev$dfo_began)

# merge calculated variables with event year
ev.tot = merge(ev[,c("system.index", "year")], ev.tot, by = "system.index")


# data cleaning
#protection levels
ev.tot[which(ev.tot$prot.cst.mn == "NaN"),"prot.cst.mn"] <- NA
ev.tot[which(ev.tot$prot.riv.mn == "NaN"),"prot.riv.mn"] <- NA
ev.tot[which(ev.tot$prot.tot.mn == "NaN"),"prot.tot.mn"] <- NA

ev.tot$prot.tot.mn = ifelse(is.na(ev.tot$prot.tot.mn), ev.tot$prot.riv.mn, ev.tot$prot.tot.mn)

#poverty levels
ev.tot[which(is.infinite(ev.tot$rel.190)),"rel.190"] <- NA
ev.tot[which(is.infinite(ev.tot$rel.320)),"rel.320"] <- NA
ev.tot[which(is.infinite(ev.tot$rel.550)),"rel.550"] <- NA

ev.tot[which(ev.tot$rel.190 == "NaN"),"rel.190"] <- NA
ev.tot[which(ev.tot$rel.320 == "NaN"),"rel.320"] <- NA
ev.tot[which(ev.tot$rel.550 == "NaN"),"rel.550"] <- NA

ev.tot[which(ev.tot$cns_md == "NaN"),"cns_md"] <- NA
ev.tot[which(ev.tot$cns_md == "NaN"),"cns_mn"] <- NA

# calculate relative share of age and income characteristics
#age
ev.tot$chi_f = ev.tot$chi_f / ev.tot$pop_exp
ev.tot$chi_m = ev.tot$chi_m / ev.tot$pop_exp
ev.tot$chi_t = ev.tot$chi_t / ev.tot$pop_exp

ev.tot$wrk_f = ev.tot$wrk_f / ev.tot$pop_exp
ev.tot$wrk_m = ev.tot$wrk_m / ev.tot$pop_exp
ev.tot$wrk_t = ev.tot$wrk_t / ev.tot$pop_exp

ev.tot$eld_f = ev.tot$eld_f / ev.tot$pop_exp
ev.tot$eld_m = ev.tot$eld_m / ev.tot$pop_exp
ev.tot$eld_t = ev.tot$eld_t / ev.tot$pop_exp

# income
ev.tot$inc_gap = (ev.tot$incm_exp - ev.tot$incf_exp) / ev.tot$incm_exp

# settlements: rural from smod 11 and 12
ev.tot$rur11_12 = ev.tot$smod11_rel + ev.tot$smod12_rel


# calculate additional vars potentially of interest
ev.tot$fat.log = log(ev.tot$dfo_dead)

#replace infinites with 0s
ev.tot[which(is.infinite(ev.tot$fat.log)),"fat.log"] <- 0


# add national-level variables
# adaptation readiness index
ev.tot = merge(ev.tot, nd.gain[,c("ISO3", "X2010")], by.x = "iso", by.y = "ISO3", all.x = T)
colnames(ev.tot) [colnames(ev.tot) == "X2010"] <- "nd_gain"

# governance index
#subset
gov = subset(gov, year == 2010, select = c(countrycode, governance))
#merge
ev.tot = merge(ev.tot, gov, by.x = "iso", by.y = "countrycode", all.x = T)
colnames(ev.tot) [colnames(ev.tot) == "governance"] <- "gov"

# subset data for analysis: remove all instances where fat > E
ev.char = ev.tot[-c(which(ev.tot$dfo_dead > ev.tot$pop_exp)), ] # 911 events
ev.char = subset(ev.char, select = -c(flo, dur_mean, pop_day, f_exp, inc_gap_exp)) # drop calculated variables not needed for the analysis

# rearrange dataframe to have all variables derived from the same dataset next to each other
#determine column order
order <- c(1:4, 52, 5:18, 37:38, 19:26, 46:50, 27:36, 51, 39:45, 53:54)

#Rearrange columns using dplyr
ev.char <- ev.char %>% dplyr::select(order)

#transform variables
ev.char$pop_exp = ev.char$pop_exp / 1000000

ev.char$chi_f = ev.char$chi_f*100
ev.char$chi_m = ev.char$chi_m*100
ev.char$chi_t = ev.char$chi_t*100
ev.char$wrk_f = ev.char$wrk_f*100
ev.char$wrk_m = ev.char$wrk_m*100
ev.char$wrk_t = ev.char$wrk_t*100
ev.char$eld_f = ev.char$eld_f*100
ev.char$eld_m = ev.char$eld_m*100
ev.char$eld_t = ev.char$eld_t*100
ev.char$ydr = ev.char$ydr*100
ev.char$odr = ev.char$odr*100
ev.char$cdr = ev.char$cdr*100

ev.char$f_prop = ev.char$f_prop*100

ev.char$incf_exp = ev.char$incf_exp/ 1000000
ev.char$incm_exp = ev.char$incm_exp/ 1000000
ev.char$inct_exp = ev.char$inct_exp/ 1000000
ev.char$inc_gap = ev.char$inc_gap*100

ev.char$rel.190 = ev.char$rel.190*100
ev.char$rel.320 = ev.char$rel.320*100
ev.char$rel.550 = ev.char$rel.550*100

ev.char$smod11_rel = ev.char$smod11_rel*100
ev.char$smod12_rel = ev.char$smod12_rel*100
ev.char$smod13_rel = ev.char$smod13_rel*100
ev.char$smod21_rel = ev.char$smod21_rel*100
ev.char$smod22_rel = ev.char$smod22_rel*100
ev.char$smod23_rel = ev.char$smod23_rel*100
ev.char$smod30_rel = ev.char$smod30_rel*100
ev.char$rur_rel = ev.char$rur_rel*100
ev.char$urb_rel = ev.char$urb_rel*100
ev.char$sub_rel = ev.char$sub_rel*100
ev.char$rur11_12 = ev.char$rur11_12*100

ev.char$wlk.mn = ev.char$wlk.mn/60
ev.char$wlk.ma = ev.char$wlk.ma/60
ev.char$mot.mn = ev.char$mot.mn/60
ev.char$mot.ma = ev.char$mot.ma/60


## 2. trend analysis ##
trend = lm(fat.log ~ year, data = ev.char) 
summary(trend)


## 3. explore linear regression model ##

# MODEL 1: basic model without V
model1 = lm(fat.log ~ pop_exp + dur_tot, data = ev.char)
summary(model1)


# Semi-automated stepwise approach
#basic model for stepwise approach based on variable selection (i.e. significant variables per variable group)
model_basic = lm(fat.log ~ pop_exp + dur_tot + 
                chi_m + chi_t + wrk_t + eld_f + odr + ydr + mean_age + med_age + 
                inc_gap + minc_exp + incm_exp + mys_exp + mys_gen + 
                smod11_rel + smod12_rel + smod13_rel + smod23_rel + smod30_rel + #rur11_12 + #rur_rel + urb_rel +
                prot.riv.mn + cns_mn + cns_md + rel.190 + rel.320 + rel.550 + 
                wlk.mn + wlk.ma + gov + nd_gain , 
              data = ev.char)

summary(model_basic) #Table S4

#run stepwise model: MODEL 2
model2 = step(model_basic, direction = "both")
summary(model2)

#significance (until all vars p< 0.05)
model3 = update(model2, . ~. - med_age)
model3 = update(model3, . ~. - chi_t)
model3 = update(model3, . ~. - ydr)
model3 = update(model3, . ~. - cns_md)

#test for multicollinearity
vif(model3)

#remove collinear variables step by step
model3 = update(model3, . ~. - mys_exp)

#remove insignificant vars
model3 = update(model3, . ~. - rel.550)
model3 = update(model3, . ~. - wrk_t)
model3 = update(model3, . ~. - odr)

vif(model3)

#remove collinear vars
model3 = update(model3, . ~. - smod11_rel)

#remove insignificant vars
model3 = update(model3, . ~. - smod30_rel)
model3 = update(model3, . ~. - smod12_rel)

vif(model3)

#remove collinear vars
model3 = update(model3, . ~. - minc_exp)

#remove insignificant vars
model3 = update(model3, . ~. - prot.riv.mn)

vif(model3)

#remove collinear vars
model3 = update(model3, . ~. - cns_mn)

#remove insignificant vars
model3 = update(model3, . ~. - gov)

#remove double-counting vars
model3 = update(model3, . ~. - wlk.mn)
summary(model3) #.272 

#add other variables again based on significance and R2, ensuring that VIF <5 (test all possible variables)
test = update(model3, . ~. + eld_t) #.2779
test = update(model3, . ~. + prot.riv.mn) #.2765

vif(test)
summary(test)

# FINAL MODEL CONFIGURATION: Model 3
model3 = update(model3, . ~. + eld_t) 
summary(model3)
#model3 = lm(fat.log ~ pop_exp + dur_tot + smod13_rel + mys_gen + eld_t + 
#                 wlk.ma + inc_gap, data = ev.char) #.2779

coef = model3$coefficients
write.csv(coef, paste0(path, "regression_coefficients.csv"))

