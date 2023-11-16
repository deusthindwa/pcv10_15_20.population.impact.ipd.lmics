#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
#MALAWI DATA MANIPULATION
#====================================================================

#notes
#PCV13 introduced in Nov 2011 under 3+0 dosing schedule
#serotyping done using...

#import pre-PCV13 pediatric invasive disease data
mw_ipdb2011 <-
  rio::import(here("data", "malawi_ipd_2011-.csv")) %>%
  dplyr::select(date_collected , serotype, age_years) %>%
  dplyr::rename("datec" = "date_collected", 
                "st" = "serotype", 
                "agey" = "age_years") %>%
  dplyr::mutate(datec = lubridate::dmy(datec),
                yearc = as.integer(lubridate::year(datec)),
                agey = as.numeric(agey),
                st = as.factor(st))

#age imputation via missForest
mw_ipdb2011 <-
  mw_ipdb2011 %>%
  mutate(agey = missForest(mw_ipdb2011 %>% dplyr::select(st, yearc, agey))$ximp$agey) %>% 
  dplyr::filter(agey <=5L, yearc >2005L, st !="NoGrowth") %>%
  dplyr::select(datec, yearc, agey, st)
                
#import post-PCV13 pediatric invasive disease data
mw_ipda2015 <-
  rio::import(here("data", "malawi_ipd_2015-18.csv")) %>%
  dplyr::select(date_collected , serotype, age_years) %>%
  dplyr::rename("datec" = "date_collected", 
                "st" = "serotype", 
                "agey" = "age_years") %>%
  dplyr::mutate(datec = lubridate::dmy(datec),
                yearc = as.integer(lubridate::year(datec)),
                agey = as.numeric(agey),
                st = as.factor(st))

#age imputation via missForest
mw_ipda2015 <-
  mw_ipda2015 %>%
  mutate(agey = missForest(mw_ipda2015 %>% dplyr::select(st, yearc, agey))$ximp$agey) %>% 
  dplyr::filter(agey <=5L, st !="NoGrowth", st != "Sample_not_located") %>%
  dplyr::select(datec, yearc, agey, st)

#import post-PCV13 pediatric carriage data
mw_cara2015 <-
  rio::import(here("data", "malawi_carriage_2015-18.csv")) %>%
  dplyr::select(collection_date, serotype_latex, serotype_pneumocat, serotype_add, age_y) %>%
  dplyr::mutate(st = serotype_latex,
                st = if_else(serotype_latex =="NVT" & serotype_pneumocat !="", serotype_pneumocat,
                             if_else(serotype_latex == "NVT" & serotype_pneumocat =="", serotype_add, serotype_latex))) %>%
  dplyr::mutate(st = if_else(st =="", serotype_latex, st)) %>%
  dplyr::mutate(st = if_else(st == "Culture_negative" | st == "Failed" | st == "NT_fail", "None", st)) %>%
  dplyr::rename("datec" = "collection_date", 
                "agey" = "age_y") %>%
  dplyr::mutate(datec = lubridate::dmy(datec),
                yearc = as.integer(lubridate::year(datec)),
                agey = as.numeric(agey),
                st = as.factor(st)) %>%
  dplyr::select(datec, yearc, agey, st) %>%
  dplyr::filter(agey <=5L) #there are no missing ages

#redefine serotypes so they are sensible
mw_ipdb2011 <- 
  mw_ipdb2011 %>%
  dplyr::mutate(st = if_else(st %in% c("6A", "6B", "6A/6B"), "6A/6B", st),
                country = "Malawi")

mw_ipda2015 <- 
  mw_ipda2015 %>%
  dplyr::mutate(st = if_else(st %in% c("6A", "6B", "6A/6B"), "6A/6B", st),
                country = "Malawi")

mw_cara2015 <- 
  mw_cara2015 %>%
  dplyr::mutate(st = if_else(st %in% c("6A", "6B", "6A/6B"), "6A/6B", st),
                country = "Malawi")

#====================================================================
#ISRAEL DATA MANIPULATION
#====================================================================

#notes
#PCV7 introduced in July 2009 and replaced by PCV13 in Nov 2010 under 2+1 dosing schedule
#serotyping done using....

#import pre- and post-PCV7 pediatric invasive disease data
is_ipdb2009 <- is_ipda2013 <-
  rio::import(here("data", "israel_ipd_2009-16.csv")) %>%
  dplyr::select(Obtaining_Year, Obtaining_Month, Ethnicity, Pnc_Serotype, agem2) %>%
  dplyr::mutate(datec = lubridate::make_date(year = Obtaining_Year, month = Obtaining_Month, day = "01")) %>%
  dplyr::rename("eth" = "Ethnicity",
                "st" = "Pnc_Serotype",
                "agey" = "agem2") %>%
  dplyr::mutate(datec = lubridate::ymd(datec),
                yearc = as.integer(lubridate::year(datec)),
                agey = as.numeric(agey)/12,
                st = as.factor(st),
                eth = as.factor(eth),
                country = "Israel") %>%
  dplyr::select(datec, yearc, agey, st, eth, country)

#filter only pre-PCV7 pediatric invasive disease data
is_ipdb2009 <- 
  is_ipdb2009 %>% 
  dplyr::filter(yearc == 2009, st != "ND", st != "OMNI NEG")

#filter only post-PCV7 pediatric invasive disease data in a mature program
is_ipda2013 <- 
  is_ipda2013 %>%
  dplyr::filter(yearc >= 2013, st != "ND", st != "OMNI NEG")

#import pre- and post-PCV13 pediatric carriage data
is_carb2009 <- is_cara2013 <-
  rio::import(here("data", "israel_carriage_2009-16.csv")) %>%
  dplyr::select(swab_year, swab_month, Ethnicity, Pnc_Serotype, agem2) %>%
  dplyr::mutate(datec = lubridate::make_date(year = swab_year, month = swab_month, day = "01")) %>%
  dplyr::rename("eth" = "Ethnicity",
                "st" = "Pnc_Serotype",
                "agey" = "agem2") %>%
  dplyr::mutate(datec = lubridate::ymd(datec),
                yearc = as.integer(lubridate::year(datec)),
                agey = as.numeric(agey)/12,
                st = as.factor(st),
                eth = as.factor(eth)) %>%
  dplyr::select(datec, yearc, agey, st, eth)

#filter only pre-PCV7 pediatric carriage data
is_carb2009 <- 
  is_carb2009 %>% 
  dplyr::filter(yearc == 2009) %>%
  dplyr::mutate(st = as.factor(if_else(st == "OMNI NEG", NA_character_, st)),
                country = "Israel")

#filter only post-PCV7 pediatric carriage data in a mature program
is_cara2013 <- 
  is_cara2013 %>%
  dplyr::filter(yearc >= 2013) %>%
  dplyr::mutate(st = as.factor(if_else(st == "OMNI NEG", NA_character_, st)),
                country = "Israel")

#====================================================================
#SOUTH AFRICA DATA MANIPULATION
#====================================================================

#notes
#PCV7 introduced in April 2009 and replaced by PCV13 in Nov 2011 under 2+1 dosing schedule
#serotyping done using...

#import pre- and post-PCV7 pediatric invasive disease data
sa_ipdb2009 <- sa_ipda2015 <-
  rio::import(here("data", "southafrica_ipd_2005-19_viable.csv")) %>%
  dplyr::rename("yearc" = "YEAR",
                "st" = "SEROTYPE",
                "agey" = "age_category",
                "hiv" = "HIV",
                "ccount" = "viable_case_count") %>%
  dplyr::mutate(yearc = as.integer(yearc),
                agey = as.factor(agey),
                st = as.factor(st),
                hiv = as.factor(hiv),
                ccount = as.integer(ccount)) %>%
  dplyr::select(yearc, agey, st, hiv, ccount) %>%
  dplyr::filter(grepl("months", agey)) %>%
  dplyr::mutate(agey = as.integer(str_sub(agey, 1, 2))/12,
                hiv = as.factor(if_else(hiv == "Negative", "Neg",
                              if_else(hiv == "Positive", "Pos", NA_character_))),
                country = "South Africa") %>%
  
#expand to make a relational dataset (each row = each case)
  dplyr::arrange(yearc, st, hiv) %>% 
  utils::type.convert(as.is = TRUE) %>% 
  tidyr::uncount(ccount) %>%
  
#check serotypes so the y are sensible
  dplyr::filter(st != "NT", st != "NEG38", st != "POOL G", st !="")

#filter only pre-PCV7 pediatric invasive disease data (2005-2008)
sa_ipdb2009 <- 
  sa_ipdb2009 %>% 
  dplyr::filter(yearc <2009)

#filter only post-PCV13 pediatric invasive disease data (2015-2019)
sa_ipda2015 <- 
  sa_ipda2015 %>%
  dplyr::filter(yearc >2014)

#====================================================================
#BRAZIL DATA MANIPULATION
#====================================================================

#notes
#PCV7 introduced in March 2010 under 3+0 dosing schedule





#create vaccine serotypes (PCV7, PCV10-GSK, PCV10-SII, PCV13, PCV15, PCV20)






#====================================================================
#INVASIVENESS DATA MANIPULATION
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
  dplyr::rename("period" = "time_interval", 
                "st" = "type",  
                "nsamples"= "carriage_samples", 
                "ncarr" = "carriage",  
                "npop" = "surveillance_population", 
                "nipd" = "disease") %>%
  dplyr::mutate(prevcarr = ncarr/nsamples,
                log_prevcarr = log(prevcarr+0.5),
                log_npop = log(npop)) %>%
  #dplyr::filter(phase == "pre-pcv") %>%
  dplyr::select(country:ncarr, prevcarr, log_prevcarr, log_npop, nipd)


#case when code
dplyr::mutate(seas = case_when(st %in% c("6A", "6B", "6A/6B") ~ "6A/6B",
                               st %in% c("11A", "11D", "11A/D") ~ "11A/D",
                               st %in% c("11A", "11D", "11A/D") ~ "11A/D",
                               st %in% c("11A", "11D", "11A/D") ~ "11A/D",
                               st %in% c("11A", "11D", "11A/D") ~ "11A/D",
                               st %in% c("11A", "11D", "11A/D") ~ "11A/D",
                               st %in% c("11A", "11D", "11A/D") ~ "11A/D",
                               st %in% c("11A", "11D", "11A/D") ~ "11A/D",
                               st %in% c("11A", "11D", "11A/D") ~ "11A/D",
                               TRUE ~ NA_character_)) 

