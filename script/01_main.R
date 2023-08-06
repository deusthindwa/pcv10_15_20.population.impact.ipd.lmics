#Authors: Deus & Dan
#Date: 01/03/2023
#Title: Rebound to normal RSV dynamics post COVID-19 suppression

#====================================================================

#load a package "pacman" used for for installing and loading other packages
if(!require(pacman)) install.packages("pacman")

#use pacman to load packages for analysis
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr", "tidyr", "broom", "rio", "scales", "boot", "magrittr",  
                        "mvtnorm", "zoo", "stringr", "patchwork", "PropCIs", "reshape2","purrr", "tsibble", "here"))

#set seed for entire session to ensure reproducibility using a task call
addTaskCallback(function(...) {set.seed(12345); TRUE})

#turn off the task call for set seed if needed
#removeTaskCallback(1)

#====================================================================

#set directory for downloaded datasets
source("script/02_fileCache.R")

#set archiving path for downloaded datasets
source("script/03_runIfExpired.R")

#load RSV case time series datasets
source("script/04_loadCases.R")

#plot RSV time series and weekly dynamics
source("script/05_plotTsWkdyn.R")

#Compute RSV onset timing
source("script/06_computeOnset.R")

#plot RSV onset timing
source("script/06_plotOnset.R")

#compute RSV peak timing 
source("script/07_computePeak.R")

#plot RSV peak timing
source("script/07_plotPeak.R")

#compute RSV growth rates
source("script/08_computeGrowth.R")

#plot RSV growth rates
source("script/08_plotGrowth.R")

#compute intensity
source("script/09_computeIntensity.R")

#plot RSV intensity
source("script/09_plotIntensity.R")

#compute and plot univariate regression models of onset, peak, growth rates and intensity
source("script/10_univarRegress.R")

#compute multivariate regression output of onset, peak, growth rates and intensity
source("script/11_multivarRegress.R")

#dynamic time warping and time series classification
source("script/12_dynTimeWarp.R")

#dynamic time warping and time series classification
source("script/13_computeDesc.R")
