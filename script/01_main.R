#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

## NOTE
# We to measure the potential benefits of newer PCVs over existing PCVs.
# We evaluate the impact of SII's PCV10, Merck's PCV15 and Pfizer's PCV20 on childhood IPD in low-middle income settings. 
# The first set of analyses is a simple set of calculations looking at the fraction of IPD that could be prevented with these new vaccines. 
# IPD is considered endpoint because its easy to measure. 
# We could underestimate preventable burden if PPV23 was already in use but considered low-middle income countries do not have no routine PPV23 programs. 
# Newer vaccines are considered as effeCtive as PCV13 but this may be incorrect if additional serotypes reduce effectiveness.

#====================================================================

#load a package "pacman" used for for installing and loading other packages
if(!require(pacman)) install.packages("pacman")

#use pacman to load packages for analysis
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr", "tidyr", "broom", "rio", "scales", "boot", "magrittr",  
                        "mvtnorm", "zoo", "stringr", "patchwork", "PropCIs", "reshape2","purrr", "tsibble", "here"))

#set seed for entire session to ensure reproducibility using a task call
addTaskCallback(function(...) {set.seed(1988); TRUE})

#turn off the task call for to reset seed if needed
#removeTaskCallback(1)

#====================================================================

#load datasets for pneumococcal carriage and IPD during PCV13 era
source("script/02_loadData.R")


