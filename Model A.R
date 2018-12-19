# A new, mechanistic model of weather-dependent Septoria tritici blotch disease risk.
# Thomas M. Chaloner*, Helen N. Fones*, Varun Varma, Daniel P Bebber and Sarah J Gurr.
# *These authors contributed equally to this work

############################################################################################################################################################
# Supplimentary code
# Model A
# Model A1, A2, and A3 code are identical, except for tmin, topt and tmax specification for "growth".
# Sections that begin with "!!!" require selection based on spore type (ascospore vs. pycnidiospore) and biological process being modelled (germination vs. growth)
############################################################################################################################################################

# Load libraries
library(ggplot2); library(ncdf4); library(R.utils); library(raster); library(accelerometry); library(abind)

# Load climate data
load("CT hourly.Rdata"); load("CM hourly.Rdata"); load("Rain_hourly.Rdata")
# Objects loaded into R environment explained below
# ctx = temperature, cmx = wetness, final_rain = precipitation (3d climate arrays) - climate data obtained at 3h estimates from JRA55 and interpolated (wetness and temperature) or divided (precip) to hourly estimates.

# Thermal death survival function - the act of drying (hour 1 of dry period) - Fig. 2
survival_function_dry <- function(x, c = 1.47885, m = 0.16889){ # c and m calculated from regression analysis of percentage Z tritici thermal death
  prop_dead_1h <- ((m*x + c) / 100) # turn death into a proportion (/ 100)
  prop_dead_1h[prop_dead_1h < 0] <- 0 # if regression predicts negative death, change to 0 (cannot have negative death - this is at very cold (< -8 oC) temperatures)
  prop_dead_1h <- prop_dead_1h + 0.3612 # 0.3612 calculated as the average death due to the act of drying (as opposed to the temperature-dependent death during the dry period)
  surv <- 1 - prop_dead_1h # proportion survival is 1 - death
  surv
}

# Thermal death survival function - all hours of dry after hour 1 - Fig. 2
# Same eq. as above, excluding the 0.3612 additional death
survival_function_dry_continous <- function(x, c = 1.47885, m = 0.16889){
  prop_dead_1h <- ((m*x + c) / 100)
  prop_dead_1h[prop_dead_1h < 0] <- 0 
  prop_dead_1h <- prop_dead_1h
  surv <- 1 - prop_dead_1h
  surv
}

# No death under wet conditions - Fig. 2

# Beta function to calculate temperature dependent rate (0, 1) Magarey et al., (2005) (https://apsjournals.apsnet.org/doi/abs/10.1094/PHYTO-95-0092)
rt <- function(x, tmin, topt, tmax){
  r <- ((tmax - x)/(tmax - topt))*((x - tmin)/(topt - tmin))^((topt-tmin)/(tmax - topt))
  r[r < 0] <- 0; r[is.na(r)] <- 0 # Where rates are calculated as <0, set to 0 (this occurs <tmin and >tmax).
  r
}

# Create hour, year and month vectors (timestamp labels)
times <- as.POSIXlt(3600*as.numeric(dimnames(cmx)$hour), origin = "1800-01-01")
hourly.num <- head(dimnames(cmx)$hour,1):tail(dimnames(cmx)$hour,1)
hourly <- as.POSIXlt(seq(head(times,1), tail(times,1), by = 3600)) # a running number of hours (all values are unique)
years <- hourly$year+1900 # what year it was - a values of year repeat for as many times as there are hours in a year
months <- hourly$mon + 1 # repeated seq of 1:12 same length as number of hours in the whole times series.

##########################################################################################################################################################################################################################################################################################
# Create function of germination/growth

germ <- function(canopy.wetness, tmp, hour, tmin, topt, tmax, a, g, array_germ_risk_vector, precip_combined_interest, cts_october_nov_interest){
	# For a single pixel, for a single growth season
    # canopy.wetness = canopy wetness data
	# tmp = air temperature data
	# hour = list all the hours for a pixel for a given
	# !!! tmin, topt, tmax = beta parameters for temperature dependent rate
	# !!! a, g (α and γ) = scale and shape parameter for Weibull distribution
	# !!! array_germ_risk_vector = germination values (from .~GERMINATION.R) - same str as canopy.wetness and tmp  
    # !!! precip_combined_interest = precip data
    # !!! cts_october_nov_interest - subset of hours for october and november

  if(all(is.na(canopy.wetness))) return(NA) else # exit the function if all values of canopy.wetness are NA, i.e. because pixel is in the sea
  canopy.wetness <- as.numeric(canopy.wetness) # confirm a vector of numbers
  canopy.wetness[canopy.wetness > 0] <- 1 # convert wetness to binary (is it wet or not)

  tmp <- as.numeric(tmp) # convert to vector of temperature values
  h <- hour[1:length(canopy.wetness)] #a list of all hour values in a growing season
  ct <- tmp # rename for logical names temperature
  wt <- canopy.wetness # rename for logical names wetness
  
  # Temperature matrix
  gs_length <- length(ct)-1 # gs_length = length of growth season - 1
  ti_temp <- matrix(ct[1:gs_length]+diff(ct)/2, nc = gs_length, nr = gs_length, byrow = F) # matrix of temperature mid-points
  rate <- rt(ti_temp, tmin, topt, tmax) # calulates the relative rates (0,1) along the beta function of each temperature
  im_temp <- matrix(1:gs_length, nc = gs_length, nr =  gs_length) - rep((1:gs_length)-1,each=gs_length) 
  im0 <- im_temp-1
  
  # Wetness matrix
  ti_wet <- matrix(wt[1:gs_length], nc = gs_length, nr = gs_length, byrow = F) # create a matrix for wetness, where last value is removed (to ensure same dimension as ti_temp)
	# midpoints not used as with temperature
	# Similar matrix to ti_temp , remmber the values are binary (dry or wet, 0 or 1)
  
  # Hazard function
  Hm_temp <- rate*((im_temp/a)^g - (im0/a)^g) # Probability of germination/growth multiplied by temperature dependence rate
  dimnames(Hm_temp) <- list("i" = 1:gs_length, "j" = 1:gs_length)
  Hm_temp[ti_wet == 0] <- 0 # if dry then no probability of germination/growth
  Hm_temp[im0 < 0] <- 0 # make upper triangle 0, we are interested in the lower part of the triangle
  Hm <- Hm_temp # instantaneous hazard at a given hour, given all hours that have passed
  
  # You will need these matrices later in model, depending on whether the hour is wet or dry. The numbers are NOT correct, but matrix will act as a place holder and update later
  prop_survive_dry <- survival_function_dry(ti_temp) # temperature dependent survival for dry hour, when wet previous hour
  prop_dead_dry <- 1 - prop_survive_dry # temperature dependent death for dry hour, when wet previous hour

perc_survive_dry_continous <- survival_function_dry_continous(ti_temp) # temperature dependent survival for dry hour, when it was dry previous hour


  # Turn all the top right hand corner to 1s (survival) and 0s (death) - ensures cumprod equation is not effected by cohorts that have not yet landed on a leaf surface
  prop_survive_dry[im0 < 0] <- 1 # survival = 1, spores have not yet landed
  perc_survive_dry_continous[im0 < 0] <- 1 # survival = 1, spores have not yet landed
  prop_dead_dry[im0 < 0] <- 0 # death = 0, spores have not yet landed
  
  prop_survive_wet <- prop_survive_dry; prop_survive_wet[,] <- 1 # survival under wet conditions = 1 (no death under wet conditions)
  prop_dead_wet <- prop_dead_dry; prop_dead_wet[,] <- 0 # death under wet conditions = 0 (no death under wet conditions)
  
  # Determine starting matrix depending on whether it is wet or dry in hour 1
  if(ti_wet[1,1] == 1){
    av <- prop_survive_wet # number available to transition in hour 1
    perc_dead <- prop_dead_wet
  } else {
    av <- prop_survive_dry # number available to transition in hour 1
    perc_dead <- prop_dead_dry
  } 
  
  # Determine hazard and thermal death in subsuquent hours. For loop b/c dependent on previous hour
  Fm <- ti_temp * 0 # create a new empty hazard matrix, with the same dimensions as ti_temp/ti_wet
  Fm[1,] <- av[1,]*(1-exp(-Hm[1,])) # Fm[1,] = number that transition to next stage of life cycle, given that x died (in this case for h = 1)
  
  i2 <- vector()
  N <- nrow(ti_temp)
  
  # Each hour will have a hazard and thermal death dependent on whether it is wet or dry + temperature
  for(i in 2:N){ # from 2nd hour in growth season to Nth hour
    i2 <- i - 1 # id for previous hour
    if(ti_wet[i,1] == 1){ # if  wet in current hour
      av[i,] <- (av[i2,] - Fm[i2,]) * 1 # the number available to germinate in hour i = (available to germ in i2 - those that did germ in i2) * 1 (because no death under wet conditions)
    } else { # current hour is dry
      if(ti_wet[i2,1] == 0){ # if dry in previous hour
        av[i,] <- (av[i2,] - Fm[i2,]) * perc_survive_dry_continous[i,] # the number available to germinate in hour i = (available to germ in i2 - those that did germ in i2) * the fraction that survived dryness in that hour (calculated from function above) i.e. 0.6 i.e. 60%
      } else { # previous hour was wet (the act of drying has occured)
        av[i,] <- (av[i2,] - Fm[i2,]) * prop_survive_dry[i,] # same as above but uses function that includes act of drying in death
      }
    }
    perc_dead[i,] <- (av[i2,] - Fm[i2,]) - av[i,] # number that didnt germinate in previous hour - those now available to germinate
    Fm[i,] <- av[i,]*(1-exp(-Hm[i,])) # proportion that transition (germinate or grow and penetrate stomata). Calculated from the cohort size availabe * the temperature dependent instantaneous hazard
  }
  
  ############################################################################################################
  ############################################################################################################
  # !!! Ascospore delivery - ensure ascospore cohort size is seasonally variable (used for ascospore germination model)
  season_length_Fm <- dim(Fm)[2] # length of the growing season (N - 1 climate datae, due to the interpolation of temperature data)
  n <- 1
  t <- seq(0, season_length_Fm, by = 1)
  h <- (0.0035 / 2.2) # paramteres  chosen to provide an approximate fit of peak during winter and tail to summer
  ascospore_burden <- (n*t^2*exp(-h*t)) / max(n*t^2*exp(-h*t)) # season ascospore budren determined from equation provided by Kitchen et al., (2016), with modeified parameters (see above) (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161887)
  ascospore_burden_matrix <- matrix(ascospore_burden[1:season_length_Fm], nc = season_length_Fm, nr = season_length_Fm, byrow = TRUE)
  Fm <- Fm * ascospore_burden_matrix # multiply germination and ascospore cohort size matrices together. I.e. at peak winter cohort size = 1, so germination is * 1.
  ############################################################################################################
  ############################################################################################################
  
  ############################################################################################################
  ############################################################################################################
  # !!! Pycnidiospore germination relies on preciptation (spores are released from pycnidia and land on available leaves only when it rains)
  precip_combined_interest[precip_combined_interest > 0] <- 1 # convert precip to binary
  precip_combined_interest <- precip_combined_interest[1:(length(precip_combined_interest) - 1)] # vectorise
  precip_status_cohort <- matrix(precip_combined_interest[1:length(precip_combined_interest)], nc = length(precip_combined_interest), nr = length(precip_combined_interest), byrow = TRUE) # create precip matrix. byrow = T, so hours that did not recieve precip, the whole column is converted to 0. If there was precip, mutiplied by 1.
  Fm <- Fm * precip_status_cohort # Germination only occurs during periods of precip
  
  # Remove pycnidia landing on leaf in oct or nov (we assume that there is a 2 month delay, the approximate latent period for infection from Bernard et al., (2013) for average temperature experienced during october).
  hours_to_remove_pyc <- length(cts_october_nov_interest) # length of period
  PYC_REMOVE <- rep(0, hours_to_remove_pyc) # cohort size of 0
  PYC_INCLUE <- rep(1, length(tmp) - hours_to_remove_pyc - 1) # relative cohort size of 1 for 1st December onwards
  pyc_combine <- c(PYC_REMOVE, PYC_INCLUE)
  
  # create a matrix by row, that removes spore cohorts created during october or november
  delete_octnov_pyc_cohorts <- matrix(pyc_combine[1:(length(tmp) - 1)], nc = (length(tmp) - 1), nr = (length(tmp) - 1), byrow = TRUE) 
  Fm <- Fm * delete_octnov_pyc_cohorts # Periods before the 1st December are recalculated as 0 (*0), periods on or after 1st Decemver are unchanged (*1)
  ########################################################################################
  ############################################################################################################
  
  ############################################################################################################
  ############################################################################################################
  # !!! Ascospore OR pycnidiospore growth
  # Muliply by cohorts that have gemrinated - gives fraction of spores available to go on to next stage (growth)
  array_germ_risk_vector <- array_germ_risk_vector[1:(length(array_germ_risk_vector) - 1)] # create a vector of germination for hours 1 to hour N
  array_germ_risk_vector[1] <- 0 # first vector value is "NA", set to 0. this shifts it all by 1 hour, so what germinated between hour 1 and 2, can now grow between hours 2 and 3
  Fd_germination <- matrix(array_germ_risk_vector[1:length(array_germ_risk_vector)], nc = length(array_germ_risk_vector), nr = length(array_germ_risk_vector), byrow = TRUE) # array_germ_risk_vector is from germination. Create a matrix byrow = T. When this is multiplied by Fm (below) it corrects growth cohort size to be dependent on germination having to have occured
  array_germ_risk_vector <- NULL
  Fm <- Fm * Fd_germination
  ############################################################################################################
  ############################################################################################################
  
  ############################################################################################################
  Fm_cum <- apply(Fm, 2, cumsum) # cumulative germination/growth across intervals
  
  #############################################################################################################
  Fd <- diff(c(0,rowSums(Fm_cum))) # total germination / growth in a specific hour
  Fd <- c(NA, Fd) # NA is added for hour 1. This is because the number of hours in the model is N - 1 growing season, because temperature is interpolated. Adding "NA" to hour 1 corrects for this, so the model output has the same dimensions as the climate data.
  Fd 
}

  ############################################################################################################

# Run model
# !!! specify tmin, topt, tmax, α (a) and γ (g)

X <- length(unique(years)); years_for_loop <- unique(years)    # number of years + vector of years
ind <- expand.grid(1:7, 1:9) # ID for pixels # expansion for locations in a grid

# Set working directory

for(j in 1:X){ # loop runs for number of years in the data
	
	# For all locations, pull out data for 1 growing season which is last 3 months of year j and first 7 months od year j+1 (total of 10 months)
  cms_winter <- cmx[,,years == years_for_loop[j] & c(months == 10 | months == 11 | months == 12)] # extract wetness values for last 3 months of the jth year

  cms_sprsum <- cmx[,,years == years_for_loop[j + 1] & c(months == 1 | months == 2 | 
                                                           months == 3 | months == 4 |
                                                           months == 5 | months == 6 |
                                                           months == 7)]
  
	# 10 months growing season data for temperature - similar logic to wetness (above)
  cts_winter <- ctx[,,years == years_for_loop[j] & c(months == 10 | months == 11 | months == 12)]
  cts_sprsum <- ctx[,,years == years_for_loop[j + 1] & c(months == 1 | months == 2 | 
                                                           months == 3 | months == 4 |
                                                           months == 5 | months == 6 |
                                                           months == 7)]
	
	# winter and sprsum variables ofr cms and cts = 3d Arrays
	# Combine the winter values (3 months) with spring + summer values (7 months) in the 3rd dimension array of all locations (X, Y) and Z length=10
  cms_season <- abind(cms_winter, cms_sprsum, along = 3)
  cts_season <- abind(cts_winter, cts_sprsum, along = 3)

	# extract the hour labels from the climate data and store to a vector
  hour_season <- as.numeric(dimnames(cms_season)[[3]])

  ############################################################################################################  
  ############################################################################################################
  # !!! Ascospore growth
	# Loading ascospore germination data for the same growth season - external file which is generated on a growing season by growing season basis.
	# as the loop runs, this data gets overwritten by the growing season of interest for the jth iteration
  load(paste("ASCO GERMIN Risk", years_for_loop[j], "to", years_for_loop[j + 1], "R.data")) # jth will rewrite over the last array_germ_risk
  ############################################################################################################
  ############################################################################################################  

  ############################################################################################################
  ############################################################################################################  
  # !!! Pycnidiospore germination
  # FOR ALL locations
  # Pull out precipitation data for 1 growing season which is last 3 months of year j and first 7 months od year j+1 (total of 10 months)
  precipx_winter <- final_rain[,,years_rain == years_for_loop[j] & c(months_rain == 10 | months_rain == 11 | months_rain == 12)]
  precipx_sprsum <- final_rain[,,years_rain == years_for_loop[j + 1] & c(months_rain == 1 | months_rain == 2 | 
                                                                             months_rain == 3 | months_rain == 4 |
                                                                             months_rain == 5 | months_rain == 6 |
                                                                             months_rain == 7)] # ensure year is following year
  # Winter and sprsum variables ofr cms and cts = 3d Arrays
  # Combine the winter precipitation values (3 months) with spring + summer values (7 months) in the 3rd dimension array of all locations (X, Y) and Z length=10
  precip_combined <- abind(precipx_winter, precipx_sprsum, along = 3)
  hour_interest_rain <- as.numeric(dimnames(precip_combined)[[3]])
  
  # Subset hours during october or novemeber of the growing season - we restrict pycnidiospores so that they cannot land during this period, but we include these periods in model so that array dimensions of ascospores and pycnidiospores are the same
  cts_october_nov <- ctx[,, years == years_for_loop[j] & c(months == 10 | months == 11)] # cts_october_nov will enable pycnidiospore germination to be restricted to 1st December onwards
  ############################################################################################################
  ############################################################################################################
  
  ############################################################################################################
  ############################################################################################################
  # !!! Pycnidiospore growth
  # Loading pycnidiospore germination data for the same growth season - external file which is generated on a growing season by growing season basis.
  # as the loop runs, this data gets overwritten by the growing season of interest for the jth iteration
  load(paste("PYC GERMIN Risk", years_for_loop[j], "to", years_for_loop[j + 1], "R.data")) # jth will rewrite over the last array_germ_risk
  ############################################################################################################
  ############################################################################################################
 
	# Create an output storage list where each location in a list element
  All.summary_growth_risk <- list()
  
  for(i in 1:63){ # for each pixel in the dataset
    All.summary_growth_risk[[i]] <- germ(canopy.wetness = cms_season[ind[i,1], ind[i,2],], # all wetness for a particular pixel, for a particular season
                                      tmp = cts_season[ind[i,1], ind[i,2],], # all temps for a particular pixel, for a particular season
                                      hour = hour_season, # all the hours for a pixel for a given growing season
                                      tmin = tmin, topt = topt, tmax = tmax, a = a, g = g,
                                      array_germ_risk_vector = array_germ_risk[ind[i,1], ind[i,2],], # !!! Ascospore OR pycnidiospore germination - germination data for a particular pixel's growing season, calculated from germination model for particular spore
                                      precip_combined_interest = precip_combined[ind[i,1], ind[i,2],], # !!! Pycidiospore germination - all precip values for a particular pixel, for a particular season
                                      cts_october_nov_interest = cts_october_nov[ind[i,1], ind[i,2],]) # !!! Pycnidospore germination - subset of hours during october or novemeber of the growing season

  }
  
  # Model output for 63 pixels sepertely held as seperate objects in list - change data structure to 3d array (same as climate data)
  array_infection_risk <- array(NA, dim = c(7,9,length(hour_season))) # same lat and lon as cms/cts
  for(i in 1:63){
    array_infection_risk[ind[i,1],ind[i,2],] <- All.summary_growth_risk[[i]]
  }
  dimnames(array_infection_risk) <- dimnames(cms_season)
  
  # Save model output indiviudally for each growth season
  save(array_infection_risk, file = paste("[SPORE TYPE] [PROCESS] Risk", years_for_loop[j], "to", years_for_loop[j + 1], "R.data"))
  
  cms_winter <- NULL; cms_sprsum <- NULL; cts_winter <- NULL; cts_sprsum <- NULL; cms_season <- NULL; cts_season <- NULL; hour_season <- NULL # clear objects
  
}
