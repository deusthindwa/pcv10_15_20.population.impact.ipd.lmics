#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines against pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

##CONSIDERATION
# We to measure the potential benefits of newer PCVs over existing PCVs.
# We evaluate the impact of SII's PCV10, Merck's PCV15 and Pfizer's PCV20 on childhood IPD in low-middle income settings. 
# IPD is considered endpoint because its easy to measure. 
# We could underestimate preventable burden if PPV23 was already in use but considered countries do not have routine adult PPV23 vaccination. 
# Newer vaccines may not be as effective as PCV13 due to serological inferiority of higher valency PCVs.
# Hence, consider serological inferiority of higher valency vaccines leading to faster waning below protective thresholds
# Analysis assumes each country's PCV regimen/schedule under which carriage/IPD was observed in mature PCV era. Schedule changes may change results

#====================================================================

#load a package "pacman" used for for installing and loading other packages
if(!require(pacman)) install.packages("pacman")

#use pacman to load packages for analysis
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr", "tidyr", "broom", "rio", "scales", "boot", "magrittr", "MASS", "ggridges",  
                        "mvtnorm", "zoo", "stringr", "patchwork", "PropCIs", "reshape2","purrr", "tsibble", "missForest", "here"))

#set seed using a task call for entire session to ensure reproducibility
addTaskCallback(function(...) {set.seed(1988); TRUE})

#turn off the task call to reset seed if needed
#removeTaskCallback(1)

#====================================================================

#load datasets for pneumococcal carriage and IPD for different countries
source("script/02_manageData.R")

#fit pre-PCV IPD (+ invasiveness) in Malawi and South Africa to infer pre-PCV carriage
source("script/03_carriageInfer.R")

#describe pneumococcal serotype carriage and IPD isolates
source("script/04_descCarrIPD.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/05_predModels.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/06_evalModels.R")

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/087_vaccineImpact.R")
