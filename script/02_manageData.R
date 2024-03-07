#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines on pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#Malawi data
#pcv13 introduced in Nov 2011 under 3+0 dosing schedule
#disease serotypes only done for pcv13 serotypes 
#carriage serotypes limited to all pcv13 and a few nonpcv13 serotypes

#import pre-pcv13 ipd data
mw_ipd <-
  bind_rows(
    rio::import(here("data", "malawi_ipd_2005_11.csv")) %>%
      dplyr::select(date_collected , serotype, age_years) %>%
      dplyr::rename("datec" = "date_collected", 
                    "st" = "serotype", 
                    "agey" = "age_years") %>%
      dplyr::mutate(datec = lubridate::dmy(datec),
                    yearc = as.integer(lubridate::year(datec)),
                    agey = as.numeric(agey),
                    st = as.factor(st)),
    
    rio::import(here("data", "malawi_ipd_2015_18.csv")) %>%
      dplyr::select(date_collected , serotype, age_years) %>%
      dplyr::rename("datec" = "date_collected", 
                    "st" = "serotype", 
                    "agey" = "age_years") %>%
      dplyr::mutate(datec = lubridate::dmy(datec),
                    yearc = as.integer(lubridate::year(datec)),
                    agey = as.numeric(agey),
                    st = as.factor(st))) 

#age imputation via missForest
mw_ipd <-
  mw_ipd %>%
  dplyr::mutate(agey = missForest(mw_ipd %>% dplyr::select(st, yearc, agey))$ximp$agey) %>% 
  dplyr::filter(agey <5, yearc >2005L, st != "NoGrowth", st != "Sample_not_located") %>%
  dplyr::select(datec, yearc, agey, st)

#import post-PCV13 pediatric carriage data
mw_car <-
  rio::import(here("data", "malawi_carriage_2015_18.csv")) %>%
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
  dplyr::filter(agey <5) %>% #there are no missing ages
  dplyr::mutate(st = if_else(st == "None", NA_character_, st))

#====================================================================

#Israel data
#PCV7 introduced in July 2009 and replaced by PCV13 in Nov 2010 under 2+1 dosing schedule
#complete serotypes are available

#import pre- and post-pcv7/13 ipd data
is_ipd <-
  rio::import(here("data", "israel_ipd_2009_16.csv")) %>%
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
  dplyr::select(datec, yearc, agey, st, eth, country) %>% 
  dplyr::filter(st != "ND", st != "OMNI NEG", agey <5)

#import pre- and post-pcv7/13 carriage data
is_car <- 
  rio::import(here("data", "israel_carriage_2009_16.csv")) %>%
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
  dplyr::select(datec, yearc, agey, st, eth) %>% 
  dplyr::mutate(st = as.factor(if_else(st == "OMNI NEG", NA_character_, st)), country = "Israel") %>%
  dplyr::filter(agey <5)

#====================================================================

#South Africa data
#PCV7 introduced in April 2009 and replaced by PCV13 in Nov 2011 under 2+1 dosing schedule
#only ipd data is available with all serotypes completed whereas as carriage data is not available

#import pre- and post-pcv7/13 ipd data
sa_ipd <-
  rio::import(here("data", "southafrica_ipd_2005_19.csv")) %>%
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
  
#check serotypes so they all have serotypes
  dplyr::filter(st != "NT", st != "NEG38", st != "POOL G", st !="", agey <5) %>%
  dplyr::select(yearc, agey, st, country)

#====================================================================

#Brazil data
#pcv10 introduced in March 2010 under 3+0 dosing schedule

#import pre- and post-pcv10 ipd data
br_ipd <-
  rio::import(here("data", "brazil_ipd_2005_19.xlsx")) %>%
  dplyr::rename("yearc" = "year",
                "st" = "serotype",
                "agey" = "age (months)") %>%
  dplyr::mutate(yearc = as.integer(yearc),
                agey = as.numeric(agey)/12,
                st = as.factor(st)) %>%
  dplyr::select(yearc, agey, st) %>%
  dplyr::mutate(country = "Brazil") %>%

#check serotypes so they all have serotypes
  dplyr::filter(st != "NT", agey < 5)

#import pre- and post-pcv10 carriage data
br_car <-
  rio::import(here("data", "brazil_carriage_2010_13_17.xlsx")) %>%
  dplyr::rename("yearc" = "YEAR_COLLECTION",
                "st" = "SEROTYPE",
                "agey" = "AGE_MONTHS") %>%
  dplyr::mutate(yearc = as.integer(yearc),
                agey = as.numeric(agey)/12,
                st = as.factor(st)) %>%
  dplyr::select(yearc, agey, st) %>%
  dplyr::mutate(country = "Brazil") %>%
  
#check serotypes so they all have serotypes
  dplyr::filter(agey < 5) %>%
  dplyr::mutate(st = if_else(st == "negative" | st == "NT", NA_character_, st))

#====================================================================

#import invasiveness data
st_inv <-
  rio::import(here("data", "st_invasiveness.csv")) %>%
  dplyr::filter(st != "NT") %>%
  mutate(pcv7pfz = if_else(grepl("\\b(4|6B|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV7", "NVT"),
         pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10-SII", "NVT"),
         pcv10gsk = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV10-GSK", "NVT"),
         pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
         pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "NVT"),
         pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "NVT"),
         inv = as.numeric(exp(log.inv))) %>%
  dplyr::select(st, inv, everything(), -log.inv.prec, -log.inv)

#calculate weighted mean invasiveness of VT and NVT 
inv_vtg <- st_inv %>% group_by(pcv20pfz) %>% mutate(inv_vtg = weighted.mean(inv)) %>% ungroup() %>% distinct(pcv20pfz, inv_vtg)
inv_vt <- inv_vtg$inv_vtg[1]
inv_nvt <- inv_vtg$inv_vtg[2]
rm(inv_vtg)

#====================================================================

#load data for vaccine effectiveness against each serotype
#against serotype disease 
st_VEd <-
  rio::import(here("data", "serotype_VEd.xlsx"))

#against serotype carriage 
st_VEc <-
  rio::import(here("data", "serotype_VEc.xlsx"))
