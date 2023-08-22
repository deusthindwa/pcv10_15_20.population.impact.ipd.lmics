#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#import invasiveness data from Navajo et al. and store it in computer hard drive
#data_inv <- rio::import("https://raw.githubusercontent.com/weinbergerlab/Invasiveness_Navajo/main/Results/mcmc_invasive_single_stage.csv")
#data_inv %>% readr::write_csv(x = ., file = here("data", "data_inv.csv"))
data_inv <-
  rio::import(here("data", "data_inv.csv")) %>%
  dplyr::select(everything(), -V1, -log.inv.prec.age1) %>%
  dplyr::rename("log_inv" = "log.inv.age1")

#import carriage and ipd datasets
#data_all <- rio::import("https://raw.githubusercontent.com/nickjcroucher/progressionEstimation/main/data-raw/S_pneumoniae_infant_serotype.csv")
#data_all %>% readr::write_csv(x = ., file = here("data", "data_all.csv"))
data_all <- 
  rio::import(here("data", "data_all.csv")) %>%
  dplyr::mutate(country = word(study, 1, sep = "\\."),
                phase = if_else(str_detect(study, "pre") == TRUE, "pre-pcv", "post-pcv")) %>%
  dplyr::select(country, phase, time_interval, type, carriage_samples, carriage, surveillance_population, disease) %>%
  dplyr::rename("period" = "time_interval", "st" = "type",  "nsamples"= "carriage_samples", "ncarr" = "carriage",  "npop" = "surveillance_population", "nipd" = "disease") %>%
  dplyr::mutate(prevcarr = ncarr/nsamples,
                log_prevcarr = log(prevcarr+0.5),
                log_nipd = log(nipd+0.5)) %>%
  dplyr::filter(phase == "pre-pcv") %>%
  dplyr::select(country:ncarr, prevcarr, npop:log_nipd)

#restructure data for plotting
data_desc <-
  bind_rows(
    data_all %>%
      dplyr::select(country, phase, st, log_prevcarr) %>%
      dplyr::rename("est" = "log_prevcarr") %>%
      dplyr::mutate(cat = "carriage"),
    
    data_all %>%
      dplyr::select(country, phase, st, log_nipd) %>%
      dplyr::rename("est" = "log_nipd") %>%
      dplyr::mutate(cat = "ipd")
    )


