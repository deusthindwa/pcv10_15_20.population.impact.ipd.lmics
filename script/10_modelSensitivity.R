#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries


#====================================================================
# BASELINE
#====================================================================

#Impact dataset
impactDS <-
  dplyr::bind_rows(sa_val, is_val, br_val, mw_val) %>%
  dplyr::filter(period == 2) %>%
  dplyr::rename("prev0"="carr_est", "ipd0"="N_cases") %>%
  dplyr::mutate(newVT_10g = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_13 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_15 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L),
                newVT_20 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L))

#check VTs versus NVTs levels
impactDS %>% group_by(country, newVT_10g) %>% tally(prev0)
impactDS %>% group_by(country, newVT_10) %>% tally(prev0)
impactDS %>% group_by(country, newVT_13) %>% tally(prev0)
impactDS %>% group_by(country, newVT_15) %>% tally(prev0)
impactDS %>% group_by(country, newVT_20) %>% tally(prev0)

#match vaccine effectiveness against each serotype disease and carriage dataset
impactDS <- 
  impactDS %>%
  left_join(st_VEd %>% dplyr::select(st, VEd_pcv10g, VEd_pcv10s, VEd_pcv13, VEd_pcv15, VEd_pcv20)) %>%
  replace(is.na(.), 0) %>%
  left_join(st_VEc %>% dplyr::select(st, VEc_pcv10g, VEc_pcv10s, VEc_pcv13, VEc_pcv15, VEc_pcv20)) %>%
  replace(is.na(.), 0)

#configure proportion of NVT that needs to fill the hole in the niche (gamma)
#set the relative change in the proportion of new vaccine coverage (vax_uptake)
gamma = 1
vax_uptake = 0.85

#determined whether there will be serotype expansion or not
#set 'a' to be true if NVT and false if VT
impactDS <- 
  impactDS %>% 
  dplyr::group_by(country) %>%
  dplyr::mutate(a10g = (((1-VEc_pcv10g) >1)|((1-VEc_pcv10g) ==1 & newVT_10g !=1)),
                a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
                a13 = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT_13 !=1)),
                a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
                a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1))) %>%
  dplyr::ungroup()

#computed change in colonization post-PCV switch
impactDS <-
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(col10g = dplyr::if_else(a10g == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10g*(1-a10g)))/sum(a10g*prev0))), sum(prev0*VEc_pcv10g*(1-a10g))),
                col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*(1-a10))),
                col13 = dplyr::if_else(a13 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*(1-a13)))/sum(a13*prev0))), sum(prev0*VEc_pcv13*(1-a13))),
                col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*(1-a15))),
                col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*(1-a20)))) %>%
  dplyr::ungroup()

#compute the expected IPD incidence of each serotype in the new PCV era after switch
impactDS <- 
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(ipd10g = ipd0 * (1-VEd_pcv10g * vax_uptake) * (col10g * vax_uptake),
                ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10 * vax_uptake),
                ipd13 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (col13 * vax_uptake),
                ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15 * vax_uptake),
                ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20 * vax_uptake)) %>%
  dplyr::ungroup()

#check total ipd (all serotypes) after PCV switch
bs_samples = 10000 #number of bootstrap samples

#ISRAEL
is_ipdf <- impactDS %>% dplyr::filter(country == "Israel") 
is_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd10))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd15))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd20))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Israel")

is_impact <-
  bind_rows(
    is_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    is_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    is_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#SOUTH AFRICA
sa_ipdf <- impactDS %>% dplyr::filter(country == "South Africa") 
sa_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd10))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd15))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd20))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "South Africa")

sa_impact <-
  bind_rows(
    sa_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    sa_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    sa_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#MALAWI
mw_ipdf <- impactDS %>% dplyr::filter(country == "Malawi") 
mw_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd10))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd15))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd20))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Malawi")

mw_impact <-
  bind_rows(
    mw_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    mw_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    mw_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#BRAZIL
br_ipdf <- impactDS %>% dplyr::filter(country == "Brazil") 
br_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd15))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd20))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Brazil")

br_impact <-
  bind_rows(
    br_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    br_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    br_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#Combine all datasets
sens_baseline <-
  bind_rows(is_impact, sa_impact, mw_impact, br_impact) %>%
  dplyr::mutate(cat = "Baseline scenario")


#====================================================================
# 20% SEROTYPE REPLACEMENT
#====================================================================

#Impact dataset
impactDS <-
  dplyr::bind_rows(sa_val, is_val, br_val, mw_val) %>%
  dplyr::filter(period == 2) %>%
  dplyr::rename("prev0"="carr_est", "ipd0"="N_cases") %>%
  dplyr::mutate(newVT_10g = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_13 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_15 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L),
                newVT_20 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L))

#check VTs versus NVTs levels
impactDS %>% group_by(country, newVT_10g) %>% tally(prev0)
impactDS %>% group_by(country, newVT_10) %>% tally(prev0)
impactDS %>% group_by(country, newVT_13) %>% tally(prev0)
impactDS %>% group_by(country, newVT_15) %>% tally(prev0)
impactDS %>% group_by(country, newVT_20) %>% tally(prev0)

#match vaccine effectiveness against each serotype disease and carriage dataset
impactDS <- 
  impactDS %>%
  left_join(st_VEd %>% dplyr::select(st, VEd_pcv10g, VEd_pcv10s, VEd_pcv13, VEd_pcv15, VEd_pcv20)) %>%
  replace(is.na(.), 0) %>%
  left_join(st_VEc %>% dplyr::select(st, VEc_pcv10g, VEc_pcv10s, VEc_pcv13, VEc_pcv15, VEc_pcv20)) %>%
  replace(is.na(.), 0)

#configure proportion of NVT that needs to fill the hole in the niche (gamma)
#set the relative change in the proportion of new vaccine coverage (vax_uptake)
gamma = 0.20
vax_uptake = 0.85

#determined whether there will be serotype expansion or not
#set 'a' to be true if NVT and false if VT
impactDS <- 
  impactDS %>% 
  dplyr::group_by(country) %>%
  dplyr::mutate(a10g = (((1-VEc_pcv10g) >1)|((1-VEc_pcv10g) ==1 & newVT_10g !=1)),
                a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
                a13 = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT_13 !=1)),
                a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
                a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1))) %>%
  dplyr::ungroup()

#computed change in colonization post-PCV switch
impactDS <-
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(col10g = dplyr::if_else(a10g == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10g*(1-a10g)))/sum(a10g*prev0))), sum(prev0*VEc_pcv10g*(1-a10g))),
                col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*(1-a10))),
                col13 = dplyr::if_else(a13 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*(1-a13)))/sum(a13*prev0))), sum(prev0*VEc_pcv13*(1-a13))),
                col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*(1-a15))),
                col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*(1-a20)))) %>%
  dplyr::ungroup()

#compute the expected IPD incidence of each serotype in the new PCV era after switch
impactDS <- 
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(ipd10g = ipd0 * (1-VEd_pcv10g * vax_uptake) * (col10g * vax_uptake),
                ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10 * vax_uptake),
                ipd13 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (col13 * vax_uptake),
                ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15 * vax_uptake),
                ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20 * vax_uptake)) %>%
  dplyr::ungroup()

#check total ipd (all serotypes) after PCV switch
bs_samples = 10000 #number of bootstrap samples

#ISRAEL
is_ipdf <- impactDS %>% dplyr::filter(country == "Israel") 
is_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd10))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd15))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd20))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Israel")

is_impact <-
  bind_rows(
    is_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    is_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    is_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#SOUTH AFRICA
sa_ipdf <- impactDS %>% dplyr::filter(country == "South Africa") 
sa_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd10))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd15))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd20))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "South Africa")

sa_impact <-
  bind_rows(
    sa_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    sa_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    sa_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#MALAWI
mw_ipdf <- impactDS %>% dplyr::filter(country == "Malawi") 
mw_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd10))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd15))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd20))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Malawi")

mw_impact <-
  bind_rows(
    mw_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    mw_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    mw_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#BRAZIL
br_ipdf <- impactDS %>% dplyr::filter(country == "Brazil") 
br_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd15))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd20))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Brazil")

br_impact <-
  bind_rows(
    br_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    br_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    br_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#Combine all datasets
sens_replacement <-
  bind_rows(is_impact, sa_impact, mw_impact, br_impact) %>%
  dplyr::mutate(cat = "20% ST replacement")
  

#====================================================================
# 20% VACCINATION COVERAGE
#====================================================================

#Impact dataset
impactDS <-
  dplyr::bind_rows(sa_val, is_val, br_val, mw_val) %>%
  dplyr::filter(period == 2) %>%
  dplyr::rename("prev0"="carr_est", "ipd0"="N_cases") %>%
  dplyr::mutate(newVT_10g = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_13 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_15 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L),
                newVT_20 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L))

#check VTs versus NVTs levels
impactDS %>% group_by(country, newVT_10g) %>% tally(prev0)
impactDS %>% group_by(country, newVT_10) %>% tally(prev0)
impactDS %>% group_by(country, newVT_13) %>% tally(prev0)
impactDS %>% group_by(country, newVT_15) %>% tally(prev0)
impactDS %>% group_by(country, newVT_20) %>% tally(prev0)

#match vaccine effectiveness against each serotype disease and carriage dataset
impactDS <- 
  impactDS %>%
  left_join(st_VEd %>% dplyr::select(st, VEd_pcv10g, VEd_pcv10s, VEd_pcv13, VEd_pcv15, VEd_pcv20)) %>%
  replace(is.na(.), 0) %>%
  left_join(st_VEc %>% dplyr::select(st, VEc_pcv10g, VEc_pcv10s, VEc_pcv13, VEc_pcv15, VEc_pcv20)) %>%
  replace(is.na(.), 0)

#configure proportion of NVT that needs to fill the hole in the niche (gamma)
#set the relative change in the proportion of new vaccine coverage (vax_uptake)
gamma = 1
vax_uptake = 0.20

#determined whether there will be serotype expansion or not
#set 'a' to be true if NVT and false if VT
impactDS <- 
  impactDS %>% 
  dplyr::group_by(country) %>%
  dplyr::mutate(a10g = (((1-VEc_pcv10g) >1)|((1-VEc_pcv10g) ==1 & newVT_10g !=1)),
                a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
                a13 = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT_13 !=1)),
                a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
                a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1))) %>%
  dplyr::ungroup()

#computed change in colonization post-PCV switch
impactDS <-
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(col10g = dplyr::if_else(a10g == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10g*(1-a10g)))/sum(a10g*prev0))), sum(prev0*VEc_pcv10g*(1-a10g))),
                col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*(1-a10))),
                col13 = dplyr::if_else(a13 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*(1-a13)))/sum(a13*prev0))), sum(prev0*VEc_pcv13*(1-a13))),
                col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*(1-a15))),
                col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*(1-a20)))) %>%
  dplyr::ungroup()

#compute the expected IPD incidence of each serotype in the new PCV era after switch
impactDS <- 
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(ipd10g = ipd0 * (1-VEd_pcv10g * vax_uptake) * (col10g * vax_uptake),
                ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10 * vax_uptake),
                ipd13 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (col13 * vax_uptake),
                ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15 * vax_uptake),
                ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20 * vax_uptake)) %>%
  dplyr::ungroup()

#check total ipd (all serotypes) after PCV switch
bs_samples = 10000 #number of bootstrap samples

#ISRAEL
is_ipdf <- impactDS %>% dplyr::filter(country == "Israel") 
is_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd10))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd15))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd20))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Israel")

is_impact <-
  bind_rows(
    is_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    is_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    is_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#SOUTH AFRICA
sa_ipdf <- impactDS %>% dplyr::filter(country == "South Africa") 
sa_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd10))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd15))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd20))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "South Africa")

sa_impact <-
  bind_rows(
    sa_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    sa_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    sa_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#MALAWI
mw_ipdf <- impactDS %>% dplyr::filter(country == "Malawi") 
mw_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd10))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd15))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd20))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Malawi")

mw_impact <-
  bind_rows(
    mw_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    mw_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    mw_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#BRAZIL
br_ipdf <- impactDS %>% dplyr::filter(country == "Brazil") 
br_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd15))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd20))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Brazil")

br_impact <-
  bind_rows(
    br_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    br_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    br_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#Combine all datasets
sens_coverage <-
  bind_rows(is_impact, sa_impact, mw_impact, br_impact) %>%
  dplyr::mutate(cat = "20% vaccine coverage")


#====================================================================
# SEROTYPE 3 AS NON-VACCINE TYPE
#====================================================================

#Impact dataset
impactDS <-
  dplyr::bind_rows(sa_val, is_val, br_val, mw_val) %>%
  dplyr::filter(period == 2) %>%
  dplyr::rename("prev0"="carr_est", "ipd0"="N_cases") %>%
  dplyr::mutate(newVT_10g = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_13 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L), #consider ST3 as VT since PCV13 is the reference
                newVT_15 = if_else(grepl("\\b(1|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L), #consider ST3 as NVT since PCV15 is the comparator
                newVT_20 = if_else(grepl("\\b(1|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L)) #consider ST3 as NVT since PCV13 is the comparator

#check VTs versus NVTs levels
impactDS %>% group_by(country, newVT_10g) %>% tally(prev0)
impactDS %>% group_by(country, newVT_10) %>% tally(prev0)
impactDS %>% group_by(country, newVT_13) %>% tally(prev0)
impactDS %>% group_by(country, newVT_15) %>% tally(prev0)
impactDS %>% group_by(country, newVT_20) %>% tally(prev0)

#match vaccine effectiveness against each serotype disease and carriage dataset
impactDS <- 
  impactDS %>%
  left_join(st_VEdX %>% dplyr::select(st, VEd_pcv10g, VEd_pcv10s, VEd_pcv13, VEd_pcv15, VEd_pcv20)) %>%
  replace(is.na(.), 0) %>%
  left_join(st_VEcX %>% dplyr::select(st, VEc_pcv10g, VEc_pcv10s, VEc_pcv13, VEc_pcv15, VEc_pcv20)) %>%
  replace(is.na(.), 0)

#configure proportion of NVT that needs to fill the hole in the niche (gamma)
#set the relative change in the proportion of new vaccine coverage (vax_uptake)
gamma = 1
vax_uptake = 0.85

#determined whether there will be serotype expansion or not
#set 'a' to be true if NVT and false if VT
impactDS <- 
  impactDS %>% 
  dplyr::group_by(country) %>%
  dplyr::mutate(a10g = (((1-VEc_pcv10g) >1)|((1-VEc_pcv10g) ==1 & newVT_10g !=1)),
                a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
                a13 = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT_13 !=1)),
                a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
                a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1))) %>%
  dplyr::ungroup()

#computed change in colonization post-PCV switch
impactDS <-
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(col10g = dplyr::if_else(a10g == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10g*(1-a10g)))/sum(a10g*prev0))), sum(prev0*VEc_pcv10g*(1-a10g))),
                col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*(1-a10))),
                col13 = dplyr::if_else(a13 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*(1-a13)))/sum(a13*prev0))), sum(prev0*VEc_pcv13*(1-a13))),
                col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*(1-a15))),
                col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*(1-a20)))) %>%
  dplyr::ungroup()

#compute the expected IPD incidence of each serotype in the new PCV era after switch
impactDS <- 
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(ipd10g = ipd0 * (1-VEd_pcv10g * vax_uptake) * (col10g * vax_uptake),
                ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10 * vax_uptake),
                ipd13 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (col13 * vax_uptake),
                ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15 * vax_uptake),
                ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20 * vax_uptake)) %>%
  dplyr::ungroup()

#check total ipd (all serotypes) after PCV switch
bs_samples = 10000 #number of bootstrap samples

#ISRAEL
is_ipdf <- impactDS %>% dplyr::filter(country == "Israel") 
is_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd10))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd15))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd20))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Israel")

is_impact <-
  bind_rows(
    is_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    is_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    is_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#SOUTH AFRICA
sa_ipdf <- impactDS %>% dplyr::filter(country == "South Africa") 
sa_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd10))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd15))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd20))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "South Africa")

sa_impact <-
  bind_rows(
    sa_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    sa_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    sa_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#MALAWI
mw_ipdf <- impactDS %>% dplyr::filter(country == "Malawi") 
mw_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd10))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd15))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd20))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Malawi")

mw_impact <-
  bind_rows(
    mw_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    mw_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    mw_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#BRAZIL
br_ipdf <- impactDS %>% dplyr::filter(country == "Brazil") 
br_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd15))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd20))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Brazil")

br_impact <-
  bind_rows(
    br_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    br_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    br_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#Combine all datasets
sens_serotype3 <-
  bind_rows(is_impact, sa_impact, mw_impact, br_impact) %>%
  dplyr::mutate(cat = "Serotype 3 as NVT")


#====================================================================
#====================================================================

#plot of vaccine preventable IPD distributions
A <-
  bind_rows(sens_baseline, sens_replacement, sens_coverage, sens_serotype3) %>%
  dplyr::mutate(cat = factor(cat, levels = c("Baseline scenario", "Serotype 3 as NVT", "20% ST replacement", "20% vaccine coverage"))) %>%
  ggplot() +
  geom_density(aes(x = impact, group = cat, color = cat), size = 1, alpha = 0.5, position = position_dodge(width = 0.00)) +
  theme_bw(base_size = 20, base_family = "American typewriter") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  scale_y_continuous(limit = c(0, 16), breaks = seq(0, 15.5, 3)) +
  scale_x_continuous(limit = c(-0.3, 0.70), breaks = seq(-0.3, 0.70, 0.2)) +
  labs(title = "", x = "vaccine net impact (additional preventable IPD)", y = "probability density") + 
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 22), strip.background = element_rect(fill = "gray90")) +
  geom_text(data = dplyr::filter(mw_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = -0.15, y = 14, label = 'favors\nexisting PCV', size = 6, family = "American typewriter") +
  geom_curve(data = dplyr::filter(mw_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = -0.05, y = 9, yend = 9, xend = -0.20, linewidth = 1.5, curvature = 0, arrow = arrow(length = unit(0.7, 'cm'))) +
  geom_text(data = dplyr::filter(mw_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = 0.15, y = 14, label = 'favors\nnew PCV', size = 6, family = "American typewriter") +
  geom_curve(data = dplyr::filter(mw_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = 0.05, y = 9, yend = 9, xend = 0.20, linewidth = 1.5, curvature = 0, arrow = arrow(length = unit(0.7, 'cm'))) +
  facet_grid(factor(country, levels = c("Israel", "South Africa", "Brazil", "Malawi")) ~ pcv, scales = "free_x") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 3)) +
  guides(color = guide_legend(title = "Epidemiologic scenarios")) +
  theme(legend.text = element_text(size = 14), theme(legend.position = "right"), legend.title = element_text(size = 14), legend.key.size = unit(1.2,"cm"))

#save combined plots
ggsave(here("output", "Parameter_sensitivity.png"),
       plot = (A), 
       width = 24, height = 13, unit = "in", dpi = 300)

#====================================================================

#summarizing absolute uncertainty for scenarios
B <-
  bind_rows(sens_baseline, sens_replacement, sens_coverage, sens_serotype3) %>%
  dplyr::mutate(country = factor(country, levels = c("Israel", "South Africa", "Brazil", "Malawi")),
                cat = factor(cat, levels = c("Baseline scenario", "Serotype 3 as NVT", "20% ST replacement", "20% vaccine coverage"))) %>%
  group_by(country, pcv, cat) %>%
  summarise(imp1M = round(quantile(impact, 0.500), digits = 3),
            imp1L = round(quantile(impact, 0.025), digits = 3),
            imp1U = round(quantile(impact, 0.975), digits = 3))

#store estimates
write_csv(B, here("output", "Parameter_sensitivity.csv"))

#====================================================================

#summarizing relative uncertainty to base scenario
C <-
  dplyr::bind_rows(
    bind_cols(sens_baseline, sens_serotype3 %>% dplyr::rename("impactx"="impact", "countryx"="country", "pcvx"="pcv", "catx"="cat" )),
    bind_cols(sens_baseline, sens_replacement %>% dplyr::rename("impactx"="impact", "countryx"="country", "pcvx"="pcv", "catx"="cat" )),
    bind_cols(sens_baseline, sens_coverage %>% dplyr::rename("impactx"="impact", "countryx"="country", "pcvx"="pcv", "catx"="cat" ))) %>%
  dplyr::mutate(difImpact = (impactx-impact)) %>%
  dplyr::select(country, pcv, catx, difImpact) %>%
  dplyr::filter(!is.na(difImpact)) %>%
  dplyr::mutate(country = factor(country, levels = c("Israel", "South Africa", "Brazil", "Malawi")),
                catx = factor(catx, levels = c("Serotype 3 as NVT", "20% ST replacement", "20% vaccine coverage"))) %>%
  group_by(country, pcv, catx) %>%
  summarise(imp1M = round(quantile(difImpact, 0.500), digits = 3),
            imp1L = round(quantile(difImpact, 0.025), digits = 3),
            imp1U = round(quantile(difImpact, 0.975), digits = 3))

#store estimates
write_csv(C, here("output", "Parameter_sensitivityRel.csv"))



