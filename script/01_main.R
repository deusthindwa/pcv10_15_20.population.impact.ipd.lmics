#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines against pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#CONSIDERATIONS
# We to measure the potential benefits of newer PCVs over existing PCVs
# We evaluate the impact of SII's PCV10, Merck's PCV15 and Pfizer's PCV20 on childhood IPD in low/middle/high-income settings 
# IPD is considered endpoint because its easy to measure
# We account for the fact that newer vaccines may not be as effective as PCV13 due to serological inferiority of higher valency PCV
# Analysis assumes a country's PCV program under which carriage/IPD was observed post-PCV is matured under current PCV schedule

#====================================================================

#load a package "pacman" used for for installing and loading other packages
if(!require(pacman)) install.packages("pacman")

#load packages for analysis
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr", "tidyr", "broom", "rio", "scales", "boot", "magrittr", "MASS", "ggridges", "htmlTable", "lme4", 
                        "mvtnorm", "zoo", "stringr", "patchwork", "PropCIs", "reshape2","purrr", "tsibble", "missForest", "MLmetrics", "Metrics", "pbapply", "here"))

#set seed using a task call for entire session to ensure reproducibility
addTaskCallback(function(...) {set.seed(1988); TRUE})

#turn off the task call to reset seed if needed
#removeTaskCallback(1)

#====================================================================

#load datasets for pneumococcal carriage and IPD for different countries
source("script/02_manageData.R")

#describe pneumococcal serotype carriage and IPD isolates
source("script/03_israelData.R")

#fit pre-PCV IPD (+ invasiveness) in Malawi and South Africa to infer pre-PCV carriage
source("script/04_brazilData.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/05_southafricaData.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/06_malawiData.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/07_ipdcarrDescription.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/08_modelValidation.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/09_modelForecast.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/10_modelSensitivity.R")
