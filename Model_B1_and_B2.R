# A new, mechanistic model of weather-dependent Septoria tritici blotch disease risk.
# Thomas M. Chaloner*, Helen N. Fones*, Varun Varma, Daniel P Bebber and Sarah J Gurr.
# *These authors contributed equally to this work

############################################################################################################################################################
# Supplimentary code
# Model B1 (wetness included) / B2 (wentess excluded)
# Sections that begin with "!!!" must be included to run model B1
############################################################################################################################################################

############################################################################################################################################################
# Load libraries
library(ggplot2); library(raster); library(ncdf4); library(R.utils); library(accelerometry); library(abind)

############################################################################################################################################################
# Define Z. tritici thermal performance curve according to Bernard et al., (2013) (https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.12134)
TPC <- function(t, LPmin, Curv, topt) {
  r <- (1 / (LPmin + Curv * (t - topt) ^ 2)) / 0.001447178 # 0.001447178 is the rate at optimum temperature (maximum rate). Dividing by 0.001447178 rescales rate to be bounded between 0 and 1, enabling relative rates to be calculated - Fig. 3.
  r
}
LPmin = 691; Curv = (0.85 * 24); topt = 18.4

# LPmin (minimum latent period, hours), Curv (scaling parameter, hoC-2), topt (optimum temperature, oC)
# Note - Bernard et al., (2013) define LPmin was in dpi (28.8) and Curv was in dpioC-2 (0.85). Model B, in contrast runs in hours, fitting the temporal scale of the climate dataset, and so the equation was adjusted accordingly.
############################################################################################################################################################

############################################################################################################################################################
# Define α and γ
# Bernard et al., (2013) (https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.12134) definitions:
# SPOmax = "maximum percentage of the inoculated area covered by spores".
# LP = the time needed for a generation of the pathogen, defined as the time elapsed from inoculation to 37% of SPOmax, assessed from the sporulation area fitted curve (Suffert et al., 2013). The value of 37% corresponds to the ordinate at the point of inflection of Gompertz curve (Winsor, 1932)".

# For model B, we estimated the Weibull parameters as α = 848 and γ = 3.65, matching the infection rate of ~0.37 (0.377) at LPmin at Topt

############################################################################################################################################################
# Ascospores (sexual spores)
############################################################################################################################################################

# Load climate data
load("CT hourly.Rdata"); load("CM hourly.Rdata")
# Objects loaded into R environment explained below
# ctx = temperature, cmx = wetness (3d climate arrays) - climate data obtained at 3h estimates from JRA55 and interpolated (wetness and temperature) to hourly estimates.

# Create hour, year and month vectors (timestamp labels)
times <- as.POSIXlt(3600*as.numeric(dimnames(cmx)$hour), origin = "1800-01-01")
hourly.num <- head(dimnames(cmx)$hour,1):tail(dimnames(cmx)$hour,1)
hourly <- as.POSIXlt(seq(head(times,1), tail(times,1), by = 3600)) # a running number of hours (all values are unique)
years <- hourly$year+1900 # what year it was - a values of year repeat for as many times as there are hours in a year
months <- hourly$mon + 1 # repeated seq of 1:12 same length as number of hours in the whole times series.

##########################################################################################################################################################################################################################################################################################
# Create an (ascospore) infection function
asco_infection_function <- function(canopy.wetness, tmp, hour, LPmin, Curv, topt, a, g){
  # For a single pixel, for a single growth season
  # canopy.wetness = canopy wetness data
  # tmp = air temperature data
  # hour = list all the hours for a pixel for a given
  # LPmin, Curv, topt = Z. tritici thermal performance curve parameters - see above
  # a, g (α and γ) = scale and shape parameter for Weibull distribution
  
    if(all(is.na(canopy.wetness))) return(NA) else # exit the function if all values of canopy.wetness are NA, i.e. because pixel is in the sea
    #!!! canopy.wetness <- as.numeric(canopy.wetness) #Convert the values into a string
    #!!! canopy.wetness[canopy.wetness > 0] <- 1 # convert wetness to binary (is it wet or not)
    #!!! length(canopy.wetness) #this is the length of the time period (length of the growth season)
      
    tmp <- as.numeric(tmp) # convert the values of temp into a string
    h <- hour[1:length(tmp)] # these are the hours of the period of interest
    ct <- tmp
    #!!! wt <- canopy.wetness
    
    # Temperature matrix
    gs_length <- length(ct)-1 # gs_length = length of growth season - 1
    ti_temp <- matrix(ct[1:gs_length]+diff(ct)/2, nc = gs_length, nr = gs_length, byrow = F) # mean temperatures between hours (av temp experienced)
    rate <- TPC(ti_temp, LPmin, Curv, topt) # multiply all these temperatures in the matrix by the TPC
    im_temp <- matrix(1:gs_length, nc = gs_length, nr =  gs_length) - rep((1:gs_length)-1,each=gs_length) # wet period end time
    im0 <- im_temp-1 # matrix to remove top right hand corner (repeated part of time) (im0 contains pos + neg series of numbers) - these spores have not yet landed on the leaf surface
    
    # Wetness matrix
    #!!! ti_wet <- matrix(wt[1:gs_length], nc = gs_length, nr = gs_length, byrow = F) #create a matrix for wetness
    
    # Hazard function
    Hm_temp <- ((im_temp*rate/a)^g - (im0*rate/a)^g) # matrix giving instantaneous hazard for each cohort at interval, given hours that have passed
    dimnames(Hm_temp) <- list("i" = 1:gs_length, "j" = 1:gs_length)
    #!!! Hm_temp[ti_wet == 0] <- 0  # replace anywhere that is dry (0) with a hazard of 0 - this prevents any ascospores from successfully infecting and hence sporulating
    Hm_temp[im0 < 0] <- 0 # replace values with zeros when j < i (removes the top right side of matrix)
    Hm <- Hm_temp # instantaneous hazard at a given hour, given all hours that have passed (change name of obeject)
    Hmc <- apply(Hm, 2, cumsum) # cumulative hazard
    Fm <- (1-exp(-Hmc)) # cumulative infection at a given hour, given cumulative hazard
    
    # !!! Ascospore delivery - ensure ascospore cohort size is seasonally variable
    season_length_Fm <- dim(Fm)[2]  # length of the growing season (N - 1 climate datae, due to the interpolation of temperature data)
    n <- 1
    t <- seq(0, season_length_Fm, by = 1)
    h <- (0.0035 / 2.2) # paramteres  chosen to provide an approximate fit of peak during winter and tail to summer
    ascospore_burden <- (n*t^2*exp(-h*t)) / max(n*t^2*exp(-h*t)) # season ascospore budren determined from equation provided by Kitchen et al., (2016), with modeified parameters (see above) (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161887)
    ascospore_burden_matrix <- matrix(ascospore_burden[1:season_length_Fm], nc = season_length_Fm, nr = season_length_Fm, byrow = TRUE)
    Fm <- Fm * ascospore_burden_matrix # multiply germination and ascospore cohort size matrices together. I.e. at peak winter cohort size = 1, so germination is * 1. 
    
    #
    Fd <- diff(c(0,rowSums(Fm))) # total infection, at a specified hour
    names(Fd) <- names(ct)[-1] 
    Fd
}

##########################################################################################################################################################################################################################################################################################
# Run infection function - model B1/2 - Ascospores

# define α (a) and γ (g)
LPmin = 691; Curv = (0.85 * 24); topt = 18.4 # define parameters for Z. tritici thermal performannce curve - see above

X <- length(unique(years)); years_for_loop <- unique(years)    # number of years + vector of years
ind <- expand.grid(1:7, 1:9) #ID for pixels # expansion for locations in a grid

# Set working directory

for(j in 1:X){
  cms_winter <- cmx[,,years == years_for_loop[j] & c(months == 10 | months == 11 | months == 12)]
  cms_sprsum <- cmx[,,years == years_for_loop[j + 1] & c(months == 1 | months == 2 | 
                                          months == 3 | months == 4 |
                                          months == 5 | months == 6 |
                                            months == 7)] 
  
  cts_winter <- ctx[,,years == years_for_loop[j] & c(months == 10 | months == 11 | months == 12)]
  cts_sprsum <- ctx[,,years == years_for_loop[j + 1] & c(months == 1 | months == 2 | 
                                          months == 3 | months == 4 |
                                          months == 5 | months == 6 |
                                            months == 7)] 
  
  # For the last run it only includes winter 2016, as there is no summer 2017 in the dataset
  # Combine arrays to create single growth season array
  cms_season <- abind(cms_winter, cms_sprsum, along = 3)
  cts_season <- abind(cts_winter, cts_sprsum, along = 3)
  hour_season <- as.numeric(dimnames(cms_season)[[3]])
  inf_risk <- list()
for(i in 1:63){ # for each pixel in the dataset
  inf_risk[[i]] <- asco_infection_function(canopy.wetness = cms_season[ind[i,1], ind[i,2],], # all wetness values for a particular pixel, for a particular season
                                  tmp = cts_season[ind[i,1], ind[i,2],], # all temps for a particular pixel, for a particular season
                                  hour = hour_season, # all the hours for a pixel for a given
                                  LPmin = LPmin, Curv = Curv, topt = topt, a = a, g = g) # model parameters - see above

} 

  # Model output for 63 pixels sepertely held as seperate objects in list - change data structure to 3d array (same as climate data)
  array_infection_risk <- array(NA, dim = c(7,9,length(hour_season) - 1)) # same lat and lon as cms, and - 1 less than the #hours in the season (due to interpolation). We do not include an NA for h1 because output is not used to drive another process, like "germination" driving "growth" in models A, A1 and A2.
  for(i in 1:63){
    array_infection_risk[ind[i,1],ind[i,2],] <- inf_risk[[i]]
  }
  cms_season_dim <- cms_season[,,2:dim(cms_season)[[3]]]
  dimnames(array_infection_risk) <- dimnames(cms_season_dim)
  
  # Save model output indiviudally for each growth season
  save(array_infection_risk, file = paste("ASCO Infection Risk", years_for_loop[j], "to", years_for_loop[j + 1], "R.data"))
  
  cms_winter <- NULL; cms_sprsum <- NULL; cts_winter <- NULL; cts_sprsum <- NULL; cms_season <- NULL; cts_season <- NULL; hour_season <- NULL #clear objects

  }

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

############################################################################################################################################################
# Pycnidiospores (asexual spores)
############################################################################################################################################################

rm(list=ls()) # clear environment

############################################################################################################################################################
# Load libraries
library(ggplot2); library(raster); library(ncdf4); library(R.utils); library(accelerometry); library(abind)

############################################################################################################################################################

############################################################################################################################################################
# Redefine Z. tritici thermal performance curve according to Bernard et al., (2013) (https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.12134)

TPC <- function(t, LPmin, Curv, topt) {
  r <- (1 / (LPmin + Curv * (t - topt) ^ 2)) / 0.001447178
  r
}
LPmin = 691; Curv = (0.85 * 24); topt = 18.4

############################################################################################################################################################

# Load climate data
load("CT hourly.Rdata"); load("CM hourly.Rdata"); load("Rain_hourly.Rdata")
# Objects loaded into R environment explained below
# ctx = temperature, cmx = wetness, final_rain = precipitation (3d climate arrays) - climate data obtained at 3h estimates from JRA55 and interpolated (wetness and temperature) or divided (precip) to hourly estimates.

############################################################################################################################################################

# Create hour, year and month vectors (timestamp labels)
times <- as.POSIXlt(3600*as.numeric(dimnames(cmx)$hour), origin = "1800-01-01")
hourly.num <- head(dimnames(cmx)$hour,1):tail(dimnames(cmx)$hour,1)
hourly <- as.POSIXlt(seq(head(times,1), tail(times,1), by = 3600)) # a running number of hours (all values are unique)
years <- hourly$year+1900 # what year it was - a values of year repeat for as many times as there are hours in a year
months <- hourly$mon + 1 # repeated seq of 1:12 same length as number of hours in the whole times series.

times_rain <- as.POSIXlt(3600*as.numeric(dimnames(final_rain)$hour), origin = "1800-01-01")
years_rain <- times_rain$year+1900
months_rain <- times_rain$mon + 1

##########################################################################################################################################################################################################################################################################################
# Create an (pycnidiospore) infection function

pyc_infection_function <- function(canopy.wetness, tmp, hour, LPmin, Curv, topt, a, g, precip_combined_interest, cts_october_nov_interest){
  # For a single pixel, for a single growth season
  # canopy.wetness = canopy wetness data
  # tmp = air temperature data
  # hour = list all the hours for a pixel for a given
  # LPmin, Curv, topt = Z. tritici thermal performance curve parameters
  # a, g (α and γ) = scale and shape parameter for Weibull distribution
  # precip_combined_interest = precip data 
  # cts_october_nov_interest - subset of hours for october and november
  
  if(all(is.na(canopy.wetness))) return(NA) else # exit the function if all values of canopy.wetness are NA, i.e. because pixel is in the sea
  #!!! canopy.wetness <- as.numeric(canopy.wetness) # confirm a vector of numbers
  #!!! canopy.wetness[canopy.wetness > 0] <- 1 # convert wetness to binary (is it wet or not)
  
  tmp <- as.numeric(tmp) # convert the values of temp into a string
  h <- hour[1:length(tmp)] # these are the hours of the period of interest
  ct <- tmp
  #!!! wt <- canopy.wetness
  
  # Temperature matrix
  gs_length <- length(ct)-1 # gs_length = length of growth season - 1
  ti_temp <- matrix(ct[1:gs_length]+diff(ct)/2, nc = gs_length, nr = gs_length, byrow = F) # mean temperatures between hours (av temp experienced)
  rate <- TPC(ti_temp, LPmin, Curv, topt) # multiply all these temperatures in the matrix by the TPC
  im_temp <- matrix(1:gs_length, nc = gs_length, nr =  gs_length) - rep((1:gs_length)-1,each=gs_length) # wet period end time
  im0 <- im_temp-1 # matrix to remove top right hand corner (repeated part of time) (im0 contains pos + neg series of numbers) - these spores have not yet landed on the leaf surface
  
  # Wetness matrix
  #!!! ti_wet <- matrix(wt[1:gs_length], nc = gs_length, nr = gs_length, byrow = F) #create a matrix for wetness 
  
  # Hazard function
  Hm_temp <- ((im_temp*rate/a)^g - (im0*rate/a)^g) # matrix giving instantaneous hazard for each cohort at interval, given hours that have passed
  dimnames(Hm_temp) <- list("i" = 1:gs_length, "j" = 1:gs_length)
  #!!! Hm_temp[ti_wet == 0] <- 0  # replace anywhere that is dry (0) with a hazard of 0 - this prevents any pycnidiospores from successfully infecting and hence sporulating
  Hm_temp[im0 < 0] <- 0 # replace values with zeros when j < i (removes the top right side of matrix)
  Hm <- Hm_temp # instantaneous hazard at a given hour, given all hours that have passed (change name of obeject)
  Hmc <- apply(Hm, 2, cumsum) # cumulative hazard
  Fm <- (1-exp(-Hmc)) # cumulative infection at a given hour, given cumulative hazard
  
  # Pycnidiospore infection relies on preciptation (spores are released from pycnidia and land on available leaves only when it rains)
  precip_combined_interest[precip_combined_interest > 0] <- 1 # convert precip to binary
  precip_combined_interest <- precip_combined_interest[1:(length(precip_combined_interest) - 1)] # vectorise
  precip_status_cohort <- matrix(precip_combined_interest[1:length(precip_combined_interest)], nc = length(precip_combined_interest), nr = length(precip_combined_interest), byrow = TRUE) # create precip matrix. byrow = T, so hours that did not recieve precip, the whole column is converted to 0. If there was precip, mutiplied by 1. 
  Fm <- Fm * precip_status_cohort # Germination only occurs during periods of precip

  # Remove pycnidia landing on leaf in oct or nov (we assume that there is a 2 month delay, the approximate latent period for infection from Bernard et al., (2013) for average temperature experienced during october).
  hours_to_remove_pyc <- length(cts_october_nov_interest) # length of period
  PYC_REMOVE <- rep(0, hours_to_remove_pyc) # cohort size of 0
  PYC_INCLUE <- rep(1, length(tmp) - hours_to_remove_pyc - 1) # relative cohort size of 1 for 1st December onwards
  pyc_combine <- c(PYC_REMOVE, PYC_INCLUE)
  
  # Create a matrix by row, that removes spore cohorts created during october or november
  delete_octnov_pyc_cohorts <- matrix(pyc_combine[1:(length(tmp) - 1)], nc = (length(tmp) - 1), nr = (length(tmp) - 1), byrow = TRUE) # by row = T
  Fm <- Fm * delete_octnov_pyc_cohorts # periods before the 1st December are recalculated as 0 (*0), periods on or after 1st Decemver are unchanged (*1)
  
  #
  Fd <- diff(c(0,rowSums(Fm))) # total infection, at a specified hour
  names(Fd) <- names(ct)[-1] 
  Fd
}

##########################################################################################################################################################################################################################################################################################
# Pycnidiospores (asexual)

# define α (a) and γ (g)
LPmin = 691; Curv = (0.85 * 24); topt = 18.4 # define parameters for Z. tritici thermal performannce curve

X <- length(unique(years)); years_for_loop <- unique(years)    # number of years + vector of years
ind <- expand.grid(1:7, 1:9) # ID for pixels # expansion for locations in a grid

# Set working directory

for(j in 1:X){
  cms_winter <- cmx[,,years == years_for_loop[j] & c(months == 10 | months == 11 | months == 12)]
  cms_sprsum <- cmx[,,years == years_for_loop[j + 1] & c(months == 1 | months == 2 | 
                                                           months == 3 | months == 4 |
                                                           months == 5 | months == 6 |
                                                           months == 7)] # ensure year is following year
  
  cts_winter <- ctx[,,years == years_for_loop[j] & c(months == 10 | months == 11 | months == 12)]
  cts_sprsum <- ctx[,,years == years_for_loop[j + 1] & c(months == 1 | months == 2 | 
                                                           months == 3 | months == 4 |
                                                           months == 5 | months == 6 |
                                                           months == 7)] # ensure year is following year
  # For the last run it only includes winter 2016, as there is no summer 2017 in the dataset
  # Combine arrays to create single growth season array
  cms_season <- abind(cms_winter, cms_sprsum, along = 3)
  cts_season <- abind(cts_winter, cts_sprsum, along = 3)
  hour_season <- as.numeric(dimnames(cms_season)[[3]])
  
  # Pull out precipitation data for 1 growing season which is last 3 months of year j and first 7 months od year j+1 (total of 10 months)
  precipx_winter <- final_rain[,,years_rain == years_for_loop[j] & c(months_rain == 10 | months_rain == 11 | months_rain == 12)]
  precipx_sprsum <- final_rain[,,years_rain == years_for_loop[j + 1] & c(months_rain == 1 | months_rain == 2 | 
                                                                             months_rain == 3 | months_rain == 4 |
                                                                             months_rain == 5 | months_rain ==6 |
                                                                             months_rain == 7)] #ensure year is following year
  # Winter and sprsum variables ofr cms and cts = 3d Arrays
  # Combine the winter precipitation values (3 months) with spring + summer values (7 months) in the 3rd dimension array of all locations (X, Y) and Z length=10
  precip_combined <- abind(precipx_winter, precipx_sprsum, along = 3)
  hour_interest_rain <- as.numeric(dimnames(precip_combined)[[3]])
  
  # Subset hours during october or novemeber of the growing season - we restrict pycnidiospores so that they cannot land during this period, but we include these periods in model so that array dimensions of ascospores and pycnidiospores are the same
  cts_october_nov <- ctx[,, years == years_for_loop[j] & c(months == 10 | months == 11)] # cts_october_nov will enable pycnidiospores landing on the leaf surface to be restricted to 1st December onwards
  ###############################################################################################################
  
  inf_risk <- list()
  for(i in 1:63){ # for each pixel in the dataset
    inf_risk[[i]] <- pyc_infection_function(canopy.wetness = cms_season[ind[i,1], ind[i,2],], # all wetness values for a particular pixel, for a particular season
                          tmp = cts_season[ind[i,1], ind[i,2],], # all temps for a particular pixel, for a particular season
                          hour = hour_season, # all the hours for a pixel for a given
                          LPmin = LPmin, Curv = Curv, topt = topt, a = a, g = g, # model parameters - see above
                          precip_combined_interest = precip_combined[ind[i,1], ind[i,2],], # all precip values for a particular pixel, for a particular season
                          cts_october_nov_interest = cts_october_nov[ind[i,1], ind[i,2],]) # subset of hours during october or novemeber of the growing season  
    
  } 
  
  # Model output for 63 pixels sepertely held as seperate objects in list - change data structure to 3d array (same as climate data)
  array_infection_risk <- array(NA, dim = c(7,9,length(hour_season) - 1)) # same lat and lon as cms, and - 1 less than the #hours in the season (due to interpolation). We do not include an NA for h1 because output is not used to drive another process, like "germination" driving "growth" in models A, A1 and A2.
  for(i in 1:63){
    array_infection_risk[ind[i,1],ind[i,2],] <- inf_risk[[i]]
  }
  cms_season_dim <- cms_season[,,2:dim(cms_season)[[3]]]
  dimnames(array_infection_risk) <- dimnames(cms_season_dim)
  
  
  # Save model output indiviudally for each growth season
  save(array_infection_risk, file = paste("PYC Infection Risk", years_for_loop[j], "to", years_for_loop[j + 1], "R.data"))
  
  cms_winter <- NULL; cms_sprsum <- NULL; cts_winter <- NULL; cts_sprsum <- NULL; cms_season <- NULL; cts_season <- NULL; hour_season <- NULL
  precipx_winter <- NULL; precipx_sprsum <- NULL; precip_combined <- NULL; hour_interest_rain <- NULL; cts_october_nov <- NULL #clear objects
  
}
