#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines on pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#Malawi data
#pcv13 introduced in Nov 2011 under 3+0 dosing schedule
#disease serotypes only done for pcv13 serotypes 
#carriage serotypes limited to all pcv13 and a few nonpcv13 serotypes

#load IPD data from a database of IPD isolates from the GPS project
mw_ipd <-
  rio::import(here("data", "monocle_gps.csv")) %>% 
  dplyr::select(Country, Year, In_silico_serotype, Phenotypic_serotype, Age_years, Age_months, Age_days) %>%
  dplyr::mutate(Age_years = as.numeric(Age_years), Age_months = as.numeric(Age_months), Age_days = as.numeric(Age_days)) %>%
  dplyr::mutate(agey = if_else(is.na(Age_years) & is.na(Age_months), Age_days/365.25, if_else(is.na(Age_years) & is.na(Age_days), Age_months/12, Age_years)),
                Year = as.integer(Year)) %>%
  dplyr::select(Country, Year, In_silico_serotype, agey) %>%
  
  dplyr::filter(Country == "MALAWI", Year >= 2009L, (agey <5)) %>% #restrict to <11 years-old pop in Malawi between 2009-2015
  dplyr::mutate(In_silico_serotype = if_else(In_silico_serotype=="UNTYPABLE" | In_silico_serotype=="SWISS_NT", "NT", In_silico_serotype)) %>% 
  
  dplyr::rename("yearc"="Year", "st"="In_silico_serotype", "country"="Country") %>%
  dplyr::mutate(country = "Malawi")

#import post-PCV13 pediatric carriage data
mw_car <-
  dplyr::bind_rows(
  rio::import(here("data", "malawi_carriage_2015_18.csv")) %>%
  dplyr::select(collection_date, yearc, serotype_latex, serotype_pneumocat, serotype_add, age_y) %>%
  dplyr::mutate(st = serotype_latex,
                st = if_else(serotype_latex =="NVT" & serotype_pneumocat !="", serotype_pneumocat,
                             if_else(serotype_latex == "NVT" & serotype_pneumocat =="", serotype_add, serotype_latex))) %>%
  dplyr::mutate(st = if_else(st =="", serotype_latex, st)) %>%
  dplyr::filter(st != "Failed", st != "NT_fail") %>%
  dplyr::mutate(st = if_else(st == "Culture_negative", "None", st)) %>%
  dplyr::rename("datec" = "collection_date", 
                "agey" = "age_y") %>%
  dplyr::mutate(datec = lubridate::dmy(datec),
                #yearc = lubridate::year(datec),
                agey = as.numeric(agey),
                st = as.factor(st)) %>%
  dplyr::select(datec, yearc, agey, st) %>%
  dplyr::filter(agey <5) %>% #there are no missing ages
  dplyr::mutate(st = if_else(st == "None", NA_character_, st)),

rio::import(here("data", "malawi_carriage_2019.xlsx")) %>%
  dplyr::mutate(datec = lubridate::ymd(collection_date),
                #yearc = lubridate::year(datec),
                yearc = as.integer(yearc),
                agey = as.numeric(age_y),
                st = if_else(Serotype == "No carriage", NA_character_, Serotype),
                st = as.factor(st)) %>%
  dplyr::select(datec, yearc, agey, st))

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
  
  #dplyr::filter(hiv != "Pos") %>%

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

# #England data
# #pcv17 introduced in September 2006 and pcv13 in April 2010 under 2+1 dosing schedule
# en_ipd <-
#   rio::import(here("data", "eng_ipd_2000_2016.xlsx")) %>%
#   dplyr::mutate(yearc = as.integer(stringr::str_sub(yearc, 1, 4)),
#                 st = stringr::str_sub(st, 3, nchar(st))) %>%
#   dplyr::filter((agegp == "a <2y" | agegp == "b 2-4y"), yearc >=2005) %>% #only focus on under 5 years and from 2007 post-pcv7
#   dplyr::select(yearc, st, Count)

#====================================================================

#PCV serotype coverage
pcv10gsk = c("1","5","6B","7F","9V","14","19F","23F","4","18C") #6A added for cross-protection
pcv10sii = c("1","5","6A","6B","7F","9V","14","19F","23F","19A")
pcv13pfz = c("1","5","6A","6B","7F","9V","14","19F","23F","19A","4","18C","3")
pcv15mek = c("1","5","6A","6B","7F","9V","14","19F","23F","19A","4","18C","3","22F","33F")
pcv20pfz = c("1","5","6A","6B","7F","9V","14","19F","23F","19A","4","18C","3","22F","33F","8","10A","11A","12F","15B")

#invasiveness data incidence is calculated as cases/(1000 people per year). This doesn't change the estimates at all, just the scaling
#so if the mean incidence of st1 is 0.008, and prevalence is 0.001 (0.1%), we get an invasiveness estimate of 8, and log(invasiveness) of 2.1.
#the scaling of invasiveness doesn't really matter, especially if trying to use it in a different population because one would have to use a scaling factor to match the incidence in that other population anyway.
sg_inv <- 
  rio::import(here("data", "invasive_mcmc_single_stage_multinomial.csv")) %>%
  dplyr::select(-V1) %>%
  dplyr::mutate(sg = gsub("([0-9]+).*$", "\\1", st), sg = if_else(st=='NT','NT', sg)) %>%
  dplyr::group_by(sg) %>%
  dplyr::rename("log.inv" = "log.inv.age1") %>%
  dplyr::summarize(log.inv = sum(log.inv * log.inv.prec.age1/sum(log.inv.prec.age1))) %>% #inverse variance-weighted mean of invasiveness within serogroup
  dplyr::mutate(inv = as.numeric(exp(log.inv))) %>%
  dplyr::select(everything(), -log.inv)

#import invasiveness data
st_inv <-
  rio::import(here("data", "invasive_st.csv")) %>%
  dplyr::filter(st != "NT") %>%
  mutate(inv = as.numeric(exp(log.inv))) %>%
  dplyr::select(st, inv, everything(), -log.inv.prec, -log.inv)

#calculate weighted mean invasiveness of VT and NVT
inv_vt <- st_inv %>% dplyr::filter(st %in% pcv20pfz) %>% mutate(x = weighted.mean(inv)) %>% dplyr::select(x)
inv_nvt <- st_inv %>% dplyr::filter(!(st %in% pcv20pfz)) %>% mutate(x = weighted.mean(inv)) %>% dplyr::select(x)
inv_vt <- inv_vt$x[1]
inv_nvt <- inv_nvt$x[1]

#====================================================================

#load data for vaccine effectiveness against each serotype
#against serotype disease and carriage
st_VEd <- rio::import(here("data", "serotype_VEd.xlsx"))
st_VEc <- rio::import(here("data", "serotype_VEc.xlsx"))

#against serotype disease and carriage (sensitivty)
st_VEdX <- rio::import(here("data", "serotype_VEdX.xlsx"))
st_VEcX <- rio::import(here("data", "serotype_VEcX.xlsx"))

#colors for sampling year, and vt and nvt groups
palx <- c("2005" = "#8B0019", "2006" = "#D16A5F", "2007" = "#E358A0", "2008" = "#8872D9", "2009" = "#C04AE1" ,"2010" = "#D598D9", "2011" = "#E0B1BC", "2012" = "#E3B45B", "2013" = "#D5E1E2", "2014" = "#84967C", "2015" = "#86AADA", "2016" = "#77DAD9", "2017" = "#84E4A6", "2018" = "#75EA5C", "2019" = "#C8DE5A")
paly <- c("PCV10G" = "#C04AE1", "PCV10S" = "#84967C", "PCV13" = "#D16A5F", "PCV15" = "#77DAD9", "PCV20" = "#E3B45B")
palz <- c("PCV10G_NVT" = "#C04AE1", "PCV10S_NVT" = "#84967C", "PCV13_NVT" = "#D16A5F", "PCV15_NVT" = "#77DAD9", "PCV20_NVT" = "#E3B45B")
