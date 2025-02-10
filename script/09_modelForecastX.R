# #====================================================================
# #model forecasts for Israel
# #====================================================================
# 
# #observed post-pcv13 ipd combined with observed pre-pcv ipd and predicted post-pcv ipd
# is_ipdf <-
#   is_ipd %>%
#   dplyr::mutate(yearc = str_c("y", yearc)) %>%
#   base::split(list(.$yearc))
# 
# is_ipdf <-
#   dplyr::bind_rows(bind_rows(is_ipdf$y2014, is_ipdf$y2015, is_ipdf$y2016), .id = "yearc") %>%
#   dplyr::mutate(yearg = "2014-16") %>%
#   dplyr::group_by(yearg, st) %>%
#   dplyr::tally() %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", st)) %>%
#   dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
#   tidyr::separate_rows(., st) %>%
#   dplyr::group_by(yearg, st) %>%
#   dplyr::tally(n/3) %>%#annualise
#   dplyr::ungroup() %>%
#   dplyr::rename("ipd0" = "n")
# 
# #using carriage prevalence during pcv13 era
# is_carf <-
#   is_car %>%
#   dplyr::filter(datec >= date('2014-01-01')) %>%
#   dplyr::group_by(st) %>%
#   dplyr::tally() %>%
#   dplyr::ungroup() %>%
#   dplyr::add_row(st = c("5", "9A", "27", "24B", "36"), n = rep(0.5, 5)) %>% #either isolated in post-PCV13 IPD
#   dplyr::mutate(st = if_else(is.na(st), "NA", st), st = if_else(st == "15BC", "15B/15C", st)) %>%
#   dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
#   tidyr::separate_rows(., st) %>%
#   dplyr::group_by(st) %>%
#   dplyr::tally(n/3) %>% #annualise from 8mos
#   dplyr::ungroup() %>%
#   dplyr::mutate(prev0 = n/sum(n)) %>%
#   dplyr::filter(st !="NA") %>%
#   dplyr::rename("ncarr0" = "n")
# 
# #using post-pcv13 carriage prevalence
# is_carf <-
#   is_carf %>%
#   mutate(newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
#          newVT_13 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
#          newVT_15 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L),
#          newVT_20 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L))
# 
# #check VTs versus NVTs levels
# is_carf %>% group_by(newVT_10) %>% tally(prev0)
# is_carf %>% group_by(newVT_13) %>% tally(prev0)
# is_carf %>% group_by(newVT_15) %>% tally(prev0)
# is_carf %>% group_by(newVT_20) %>% tally(prev0)
# 
# #match vaccine effectiveness against each serotype disease and carriage dataset
# is_carf <-
#   is_carf %>%
#   left_join(st_VEd %>% dplyr::select(st, VEd_pcv10s, VEd_pcv13, VEd_pcv15, VEd_pcv20)) %>%
#   replace(is.na(.), 0)
# 
# is_carf <-
#   is_carf %>%
#   left_join(st_VEc %>% dplyr::select(st, VEc_pcv10s, VEc_pcv13, VEc_pcv15, VEc_pcv20)) %>%
#   replace(is.na(.), 0)
# 
# #configure proportion of NVT that needs to fill the hole in the niche
# #set the relative change in the proportion of new vaccine uptake
# gamma = 1
# vax_uptake = 1
# 
# #determined whether there will be serotype expansion or not
# #set 'a' to be true if NVT and false if VT
# is_carf <- is_carf %>% dplyr::mutate(a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
#                                      a13 = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT_13 !=1)),
#                                      a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
#                                      a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1)))
# 
# #computed change in colonization post-PCV switch
# is_carf <-
#   is_carf %>%
#   dplyr::mutate(col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*vax_uptake*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*vax_uptake*(1-a10))),
#                 col13 = dplyr::if_else(a13 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*vax_uptake*(1-a13)))/sum(a13*prev0))), sum(prev0*VEc_pcv13*vax_uptake*(1-a13))),
#                 col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*vax_uptake*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*vax_uptake*(1-a15))),
#                 col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*vax_uptake*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*vax_uptake*(1-a20))))
# 
# #join IPD data over the same time period 2013-16
# is_ipdf <-  is_ipdf %>% left_join(is_carf)
# 
# #compute the expected IPD incidence of each serotype in the new PCV era after switch
# is_ipdf <-
#   is_ipdf %>%
#   dplyr::mutate(ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10),
#                 ipd13 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (col13),
#                 ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15),
#                 ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20))
# 
# #check total ipd (all serotypes) after PCV switch
# bs_samples = 10000 #number of bootstrap samples
# 
# is_impact <-
#   dplyr::bind_cols(
#     pcv10 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd10))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
#     pcv15 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd15))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
#     pcv20 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd20))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1)))) %>%
#   dplyr::mutate(country = "Israel")
# 
# is_impact <-
#   bind_rows(
#     is_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10-SII"),
#     is_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
#     is_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))
# 
# #====================================================================
# #model forecasts for South Africa
# #====================================================================
# 
# #observed post-pcv13 ipd combined with observed pre-pcv ipd and predicted post-pcv ipd
# sa_ipdf <-
#   sa_ipd %>%
#   dplyr::mutate(yearc = str_c("y", yearc)) %>%
#   base::split(list(.$yearc))
# 
# sa_ipdf <-
#   dplyr::bind_rows(bind_rows(sa_ipdf$y2015, sa_ipdf$y2016, sa_ipdf$y2017, sa_ipdf$y2018, sa_ipdf$y2019), .id = "yearc") %>%
#   dplyr::mutate(yearg = "2015-19") %>%
#   dplyr::group_by(yearg, st) %>%
#   dplyr::tally() %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
#   dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
#   tidyr::separate_rows(., st) %>%
#   dplyr::group_by(yearg, st) %>%
#   dplyr::tally(n/5) %>%#annualise
#   dplyr::ungroup() %>%
#   dplyr::rename("ipd0" = "n")
# 
# #using carriage prevalence during pcv13 era
# sa_carf <-
#   sa_carInf %>%
#   dplyr::filter(yearc >= 2015) %>%
#   dplyr::select(yearc, st, prev) %>%
#   dplyr::group_by(st) %>%
#   summarise(prev0 = sum(prev)/5) %>% #annualise
#   dplyr::ungroup()
# 
# #using post-pcv13 carriage prevalence
# sa_carf <-
#   sa_carf %>%
#   mutate(newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
#          newVT_13 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
#          newVT_15 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L),
#          newVT_20 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L))
# 
# #check VTs versus NVTs levels
# sa_carf %>% group_by(newVT_10) %>% tally(prev0)
# sa_carf %>% group_by(newVT_13) %>% tally(prev0)
# sa_carf %>% group_by(newVT_15) %>% tally(prev0)
# sa_carf %>% group_by(newVT_20) %>% tally(prev0)
# 
# #match vaccine effectiveness against each serotype disease and carriage dataset
# sa_carf <-
#   sa_carf %>%
#   left_join(st_VEd %>% dplyr::select(st, VEd_pcv10s, VEd_pcv13, VEd_pcv15, VEd_pcv20)) %>%
#   replace(is.na(.), 0)
# 
# sa_carf <-
#   sa_carf %>%
#   left_join(st_VEc %>% dplyr::select(st, VEc_pcv10s, VEc_pcv13, VEc_pcv15, VEc_pcv20)) %>%
#   replace(is.na(.), 0)
# 
# #configure proportion of NVT that needs to fill the hole in the niche
# #set the relative change in the proportion of new vaccine uptake
# gamma = 1
# vax_uptake = 1
# 
# #determined whether there will be serotype expansion or not
# #set 'a' to be true if NVT and false if VT
# sa_carf <- sa_carf %>% dplyr::mutate(a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
#                                      a13 = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT_13 !=1)),
#                                      a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
#                                      a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1)))
# 
# #computed change in colonization post-PCV switch
# sa_carf <-
#   sa_carf %>%
#   dplyr::mutate(col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*vax_uptake*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*vax_uptake*(1-a10))),
#                 col13 = dplyr::if_else(a13 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*vax_uptake*(1-a13)))/sum(a13*prev0))), sum(prev0*VEc_pcv13*vax_uptake*(1-a13))),
#                 col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*vax_uptake*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*vax_uptake*(1-a15))),
#                 col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*vax_uptake*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*vax_uptake*(1-a20))))
# 
# #join IPD data over the same time period 2013-16
# sa_ipdf <-  sa_ipdf %>% left_join(sa_carf)
# 
# #compute the expected IPD incidence of each serotype in the new PCV era after switch
# sa_ipdf <-
#   sa_ipdf %>%
#   dplyr::mutate(ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10),
#                 ipd13 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (col13),
#                 ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15),
#                 ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20))
# 
# #check total ipd (all serotypes) after PCV switch
# bs_samples = 10000 #number of bootstrap samples
# 
# sa_impact <-
#   dplyr::bind_cols(
#     pcv10 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd10))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
#     pcv15 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd15))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
#     pcv20 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd20))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1)))) %>%
#   dplyr::mutate(country = "South Africa")
# 
# sa_impact <-
#   bind_rows(
#     sa_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10-SII"),
#     sa_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
#     sa_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))
# 
# #====================================================================
# #model forecasts for Malawi
# #====================================================================
# 
# #observed post-pcv13 ipd combined with observed pre-pcv ipd and predicted post-pcv ipd
# mw_ipdf <-
#   mw_ipdInf %>%
#   dplyr::group_by(st) %>%
#   dplyr::tally(ipd/4) %>%#annualise
#   dplyr::ungroup() %>%
#   dplyr::rename("ipd0" = "n")
# 
# #using carriage prevalence during pcv13 era
# mw_carf <-
#   mw_ipdInf %>%
#   dplyr::select(yearc, st, prev) %>%
#   dplyr::group_by(st) %>%
#   summarise(prev0 = sum(prev)/4) %>% #annualise
#   dplyr::ungroup()
# 
# #using post-pcv13 carriage prevalence
# mw_carf <-
#   mw_carf %>%
#   mutate(newVT_13 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
#          newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
#          newVT_15 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L),
#          newVT_20 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L))
# 
# #check VTs versus NVTs levels
# mw_carf %>% group_by(newVT_10) %>% tally(prev0)
# mw_carf %>% group_by(newVT_13) %>% tally(prev0)
# mw_carf %>% group_by(newVT_15) %>% tally(prev0)
# mw_carf %>% group_by(newVT_20) %>% tally(prev0)
# 
# #match vaccine effectiveness against each serotype disease and carriage dataset
# mw_carf <-
#   mw_carf %>%
#   left_join(st_VEd %>% dplyr::select(st, VEd_pcv10s, VEd_pcv13, VEd_pcv15, VEd_pcv20)) %>%
#   replace(is.na(.), 0)
# 
# mw_carf <-
#   mw_carf %>%
#   left_join(st_VEc %>% dplyr::select(st, VEc_pcv10s, VEc_pcv13, VEc_pcv15, VEc_pcv20)) %>%
#   replace(is.na(.), 0)
# 
# #configure proportion of NVT that needs to fill the hole in the niche
# #set the relative change in the proportion of new vaccine uptake
# gamma = 1
# vax_uptake = 1
# 
# #determined whether there will be serotype expansion or not
# #set 'a' to be true if NVT and false if VT
# mw_carf <- mw_carf %>% dplyr::mutate(a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
#                                      a13 = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT_13 !=1)),
#                                      a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
#                                      a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1)))
# 
# #computed change in colonization post-PCV switch
# mw_carf <-
#   mw_carf %>%
#   dplyr::mutate(col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*vax_uptake*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*vax_uptake*(1-a10))),
#                 col13 = dplyr::if_else(a13 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*vax_uptake*(1-a13)))/sum(a13*prev0))), sum(prev0*VEc_pcv13*vax_uptake*(1-a13))),
#                 col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*vax_uptake*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*vax_uptake*(1-a15))),
#                 col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*vax_uptake*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*vax_uptake*(1-a20))))
# 
# #join IPD data over the same time period 2013-16
# mw_ipdf <-  mw_ipdf %>% left_join(mw_carf)
# 
# #compute the expected IPD incidence of each serotype in the new PCV era after switch
# mw_ipdf <-
#   mw_ipdf %>%
#   dplyr::mutate(ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10),
#                 ipd13 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (col13),
#                 ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15),
#                 ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20))
# 
# #check total ipd (all serotypes) after PCV switch
# bs_samples = 10000 #number of bootstrap samples
# 
# mw_impact <-
#   dplyr::bind_cols(
#     pcv10 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd10))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
#     pcv15 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd15))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
#     pcv20 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd20))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1)))) %>%
#   dplyr::mutate(country = "Malawi")
# 
# mw_impact <-
#   bind_rows(
#     mw_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10-SII"),
#     mw_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
#     mw_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))
# 
# #====================================================================
# #model forecasts for Brazil
# #====================================================================
# 
# br_ipdf <-
#   br_ipd %>%
#   dplyr::group_by(yearc, st) %>%
#   dplyr::tally() %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(n = if_else(str_length(st)>3, n/3, n)) %>%
#   tidyr::separate_rows(., st) %>%
#   dplyr::filter(yearc >=2013) %>% #7yr-mature PCV10 program
#   dplyr::group_by(st) %>%
#   dplyr::tally(n/7) %>% #annualise
#   dplyr::rename("ipd0" = "n") %>%
#   dplyr::ungroup()
# 
# br_carf <-
#   br_car %>%
#   dplyr::group_by(yearc, st) %>%
#   dplyr::tally() %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(n = if_else(is.na(st), n, n/(stringr::str_count(st, fixed(' ')) + 1))) %>%
#   tidyr::separate_rows(., st) %>%
#   dplyr::filter(yearc >=2013) %>% #2y-mature PCV10 program
#   dplyr::group_by(st) %>%
#   dplyr::tally(n/2) %>% #annualise
#   dplyr::add_row(st = c("1","10F","11B","18B","24B","25A","25F","33F","36","4","5","7F","9N"), n = rep(0.5, 13)) %>% #isolated in IPD not carriage 
#   dplyr::mutate(n = n/sum(n)) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(st !="NA") %>%
#   dplyr::rename("prev0" = "n")
# 
# #using post-pcv13 carriage prevalence
# br_carf <-
#   br_carf %>%
#   mutate(newVT_10g = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, 1L, 0L),
#          newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
#          newVT_15 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L),
#          newVT_20 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, 1L, 0L))
# 
# #check VTs versus NVTs levels
# br_carf %>% group_by(newVT_10) %>% tally(prev0)
# br_carf %>% group_by(newVT_10g) %>% tally(prev0)
# br_carf %>% group_by(newVT_15) %>% tally(prev0)
# br_carf %>% group_by(newVT_20) %>% tally(prev0)
# 
# #match vaccine effectiveness against each serotype disease and carriage dataset
# br_carf <-
#   br_carf %>%
#   left_join(st_VEd %>% dplyr::select(st, VEd_pcv10s, VEd_pcv10g, VEd_pcv15, VEd_pcv20)) %>%
#   replace(is.na(.), 0)
# 
# br_carf <-
#   br_carf %>%
#   left_join(st_VEc %>% dplyr::select(st, VEc_pcv10s, VEc_pcv10g, VEc_pcv15, VEc_pcv20)) %>%
#   replace(is.na(.), 0)
# 
# #configure proportion of NVT that needs to fill the hole in the niche
# #set the relative change in the proportion of new vaccine uptake
# gamma = 1
# vax_uptake = 1
# 
# #determined whether there will be serotype expansion or not
# #set 'a' to be true if NVT and false if VT
# br_carf <- br_carf %>% dplyr::mutate(a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
#                                      a10g = (((1-VEc_pcv10g) >1)|((1-VEc_pcv10g) ==1 & newVT_10g !=1)),
#                                      a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
#                                      a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1)))
# 
# #computed change in colonization post-PCV switch
# br_carf <-
#   br_carf %>%
#   dplyr::mutate(col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*vax_uptake*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*vax_uptake*(1-a10))),
#                 col10g = dplyr::if_else(a10g == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10g*vax_uptake*(1-a10g)))/sum(a10g*prev0))), sum(prev0*VEc_pcv10g*vax_uptake*(1-a10g))),
#                 col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*vax_uptake*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*vax_uptake*(1-a15))),
#                 col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*vax_uptake*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*vax_uptake*(1-a20))))
# 
# #join IPD data over the same time period 2013-16
# br_ipdf <-  br_ipdf %>% left_join(br_carf)
# 
# #compute the expected IPD incidence of each serotype in the new PCV era after switch
# br_ipdf <- 
#   br_ipdf %>%
#   dplyr::mutate(ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10),
#                 ipd10g = ipd0 * (1-VEd_pcv10g * vax_uptake) * (col10g),
#                 ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15),
#                 ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20))
# 
# #check total ipd (all serotypes) after PCV switch
# bs_samples = 10000 #number of bootstrap samples
# 
# br_impact <-
#   dplyr::bind_cols(
#     pcv10 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
#     pcv15 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd15))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
#     pcv20 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd20))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1)))) %>%
#   dplyr::mutate(country = "Brazil")
# 
# br_impact <-
#   bind_rows(
#     br_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10-SII"),
#     br_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
#     br_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))
# 
# #====================================================================
# #summarise estimates
# #====================================================================
# 
# #summarising uncertainty
# all_impactSum <-
#   bind_rows(is_impact, sa_impact, mw_impact, br_impact) %>%
#   group_by(country, pcv) %>%
#   summarise(imp1M = round(quantile(impact, 0.500), digits = 3),
#             imp1L = round(quantile(impact, 0.025), digits = 3),
#             imp1U = round(quantile(impact, 0.975), digits = 3))
# 
# #store estimates
# write_csv(all_impactSum, here("output", "impact_summary.csv"))
# 
# #====================================================================
# #combined impact plots
# #====================================================================
# 
# #plot of vaccine preventable IPD distributions
# 
# A <-
#   ggplot() +
#   geom_density(data = bind_rows(is_impact, sa_impact, mw_impact, br_impact), aes(x = impact, group = pcv, fill = country), size = 1, alpha = 0.3) +
#   theme_bw(base_size = 16, base_family = "American typewriter") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_y_continuous(limit = c(0, 15.5), breaks = seq(0, 15.5, 3)) +
#   facet_grid(country ~ pcv, scales = "free_x") +
#   labs(title = "", x = "vaccine impact (proportion of preventable IPD)", y = "probability density") +
#   theme(legend.text = element_text(size = 14), legend.position = "none", legend.title = element_text(size = 14), legend.key.size = unit(1.2,"cm")) +
#   theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30), strip.background = element_rect(fill = "gray90")) +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
# 
#   geom_text(data = dplyr::filter(br_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = -0.30, y = 14, label = 'favors existing PCV', size = 7, family = "American typewriter") +
#   geom_curve(data = dplyr::filter(br_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = -0.05, y = 12, yend = 12, xend = -0.30, linewidth = 1.5, curvature = 0, arrow = arrow(length = unit(0.7, 'cm'))) +
#   geom_text(data = dplyr::filter(br_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = 0.25, y = 14, label = 'favors new PCV', size = 7, family = "American typewriter") +
#   geom_curve(data = dplyr::filter(br_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = 0.05, y = 12, yend = 12, xend = 0.30, linewidth = 1.5, curvature = 0, arrow = arrow(length = unit(0.7, 'cm')))
# 
# #save combined plots
# ggsave(here("output", "pcvImpact.png"),
#        plot = (A),
#        width = 20, height = 14, unit = "in", dpi = 300)
