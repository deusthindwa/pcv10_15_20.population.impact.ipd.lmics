#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#import invasiveness data from Navajo et al. and store it in computer hard drive
# data_inv <- rio::import("https://raw.githubusercontent.com/weinbergerlab/Invasiveness_Navajo/main/Results/mcmc_invasive_single_stage.csv")
# data_inv %>% readr::write_csv(x = ., file = here("output", "data_inv.csv"))
data_inv <-
  rio::import(here("output", "data_inv.csv")) %>%
  dplyr::select(everything(), -V1, -log.inv.prec.age1) %>%
  dplyr::rename("log_inv" = "log.inv.age1")

#import carriage data for different settings from Croucher et al.
# data_carr <- rio::import("https://raw.githubusercontent.com/nickjcroucher/progressionEstimation/main/data-raw/S_pneumoniae_infant_serotype.csv")
# data_carr %>% readr::write_csv(x = ., file = here("output", "data_carr.csv"))
data_carr <-
  rio::import(here("output", "data_carr.csv")) %>%
  dplyr::select(study, time_interval, type, carriage_samples, carriage) %>%
  dplyr::rename("period" = "time_interval", "nsamples" = "carriage_samples", "ncarr" = "carriage", "st" = "type")

#import ipd data for different settings from Croucher et al.
# data_ipd <- rio::import("https://raw.githubusercontent.com/nickjcroucher/progressionEstimation/main/data-raw/S_pneumoniae_infant_serotype.csv")
# data_ipd %>% readr::write_csv(x = ., file = here("output", "data_ipd.csv"))
data_ipd <-
  rio::import(here("output", "data_ipd.csv")) %>%
  dplyr::select(study, time_interval, type, surveillance_population, disease) %>%
  dplyr::rename("period" = "time_interval", "npop" = "surveillance_population", "nipd" = "disease", "st" = "type")

#====================================================================

# #simulate some case fatality ratio data
# data_cfr <- data.frame(
#   st = c("1", "2", "3", "4", "5", "6A/C", "6B", "7F", "8", "9A", "9N", "9V", "10A", "11A", "12A", "12F", "13", "14", "15A", "15B", "15C", "16F", "17F", "18C", "19A", "19F", "20", "22F", "23A", "23F", "24F", "31", "33F", "34", "35B", "35F", "38"),
#   invSt = c("1", NA, "3", "4", "5", "6B", "8", NA, "9A", "9N", NA, "10A", "11A", "12F", "13", NA, NA, "14", "15A", "15B", NA, "15C", "18C",NA, "19A", "19F", "20", "23A", "23F", NA, NA, "24F", "31", NA, "33F", "34", NA), #serotypes with invasiveness data
#   inv = sample(c(20:100), size = 37, replace = TRUE, prob = rep(0.5, 81)),
#   cfr = sample(c(10:40), size = 37, replace = TRUE, prob = rep(0.5, 31))) %>%
#   dplyr::mutate(inv = ifelse(is.na(invSt), NA_integer_, inv))

#====================================================================

# #simulate data on age/serotype-specific carriage
# data_carr <- data.frame(
#   agem = rep(c("0-6m", "7-12m", "13-18m", "19-24m", "25-30m", "31-36m", "37-42m", "43-48m", "49-54m", "55-60m"), 37)) %>%
#   dplyr::group_by(agem) %>%
#   dplyr::mutate(st = c("1", "2", "3", "4", "5", "6A/C", "6B", "7F", "8", "9A", "9N", "9V", "10A", "11A", "12A", "12F", 
#                                   "13", "14", "15A", "15B", "15C", "16F", "17F", "18C", "19A", "19F", "20", "22F", "23A", "23F", 
#                                   "24F", "31", "33F", "34", "35B", "35F", "38"),
#                 carrPos = sample(c(5:30), size = 37, replace = TRUE, prob = rep(0.5, 26)),
#                 carrSamp = 200,
#                 prev = carrPos/carrSamp) %>%
#   dplyr::ungroup()

#====================================================================

# #simulate data on age/serotype-specific invasive disease
# data_ipd <- data.frame(
#   agem = rep(c("0-6m", "7-12m", "13-18m", "19-24m", "25-30m", "31-36m", "37-42m", "43-48m", "49-54m", "55-60m"), 37)) %>%
#   dplyr::group_by(agem) %>%
#   dplyr::mutate(st = c("1", "2", "3", "4", "5", "6A/C", "6B", "7F", "8", "9A", "9N", "9V", "10A", "11A", "12A", "12F", 
#                        "13", "14", "15A", "15B", "15C", "16F", "17F", "18C", "19A", "19F", "20", "22F", "23A", "23F", 
#                        "24F", "31", "33F", "34", "35B", "35F", "38"),
#                 ipd = sample(c(20:50), size = 37, replace = TRUE, prob = rep(0.5, 31))) %>%
#   dplyr::ungroup()

#====================================================================

# #simulate data on age/serotype-specific carbons per capsular polysaccharide repeat unit)
# data_caps <- data.frame(
#   st = c("1", "2", "3", "4", "5", "6A/C", "6B", "7F", "8", "9A", "9N", "9V", "10A", "11A", "12A", "12F", 
#          "13", "14", "15A", "15B", "15C", "16F", "17F", "18C", "19A", "19F", "20", "22F", "23A", "23F", 
#          "24F", "31", "33F", "34", "35B", "35F", "38"),
#   caps = sample(c(2:15), size = 37, replace = TRUE, prob = rep(0.5, 14)))
