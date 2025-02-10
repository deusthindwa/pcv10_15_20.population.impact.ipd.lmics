#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries


# #fit a poisson mixed model for israel
# XX <-
#   mw_val %>%
#   dplyr::left_join(sa_val %>% dplyr::rename('sa' = 'N_cases') %>% dplyr::select(st, sa)) %>% dplyr::mutate(sa = if_else(is.na(sa), 0, sa)) %>%
#   dplyr::left_join(br_val %>% dplyr::rename('br' = 'N_cases') %>% dplyr::select(st, br)) %>% dplyr::mutate(br = if_else(is.na(br), 0, br)) %>%
#   dplyr::left_join(is_val %>% dplyr::rename('is' = 'N_cases') %>% dplyr::select(st, is)) %>% dplyr::mutate(is = if_else(is.na(is), 0, is))
# 
# ipd_model <- glmer(N_cases ~ log(sa+1) + log(br+1) + log(is+1) + (1|st),
#                    family = poisson(link = "log"),
#                    na.action = na.omit,
#                    data = XX,
#                    verbose = TRUE)
# 
# mw_val <-
#   XX %>%
#   dplyr::mutate(N_casesx = predict(ipd_model, type = "response")) %>%
#   dplyr::mutate(N_cases = if_else(N_cases==0, N_casesx, N_cases)) %>%
#   dplyr::select(everything(), -sa, -br, -is, -N_casesx)

#====================================================================
#model validation in England
#====================================================================

#estimate carriage prevalence using ipd and invasiveness
en_val <- 
  en_ipd %>%
  dplyr::mutate(period = if_else(yearc >=2007 & yearc <=2009, 1, #prePCV13, postPCV7
                                 if_else(yearc >=2012, 2, NA_integer_))) %>% #postPCV13
  dplyr::filter(!is.na(period)) %>%
  #tidyr::separate_rows(., st) %>%
  dplyr::group_by(period, st) %>%
  dplyr::summarize(N_cases = sum(Count)) %>%
  dplyr::ungroup() %>%
  tidyr::complete(st, period, fill = list(N_cases = 0)) %>%
  dplyr::arrange(period, -N_cases) %>%
  dplyr::mutate(sg= gsub("([0-9]+).*$", "\\1", st),
                sg = if_else(st == 'NT','NT', sg),
                country = "England",
                yrsn = if_else(period == 1, 2,
                               if_else(period == 2, 3, 
                                       if_else(period == 3, 5, NA_integer_))),
                N_cases = N_cases/yrsn) %>%
  dplyr::select(country, st, sg, period, N_cases) %>%
  
  dplyr::left_join(st_inv, by = 'st') %>%
  dplyr::mutate(inv = if_else(is.na(inv), inv_nvt, inv)) %>%
  dplyr::group_by(country, period) %>%
  dplyr::mutate(carr_est = (N_cases+1)/inv,
                carr_est = carr_est/sum(carr_est, na.rm=T),
                period = as.factor(period)) %>%
  dplyr::ungroup() %>%
  dplyr::select(everything(), -sg, -inv)

#====================================================================
#model validation in Israel
#====================================================================

is_val <-
  dplyr::full_join(
    is_ipd %>%
      dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
      dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, 1), country = "Israel", yearc = yearc) %>%
      tidyr::separate_rows(., st) %>%
      group_by(country, yearc, st) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(period = if_else(yearc == 2010, 1L, if_else(yearc >= 2014 & yearc <= 2016, 2L, NA_integer_))) %>%
      dplyr::filter(!is.na(period)) %>%
      dplyr::group_by(country, period, st) %>%
      dplyr::summarize(N_cases = sum(n)) %>%
      dplyr::ungroup() %>%
      tidyr::complete(st, period, country, fill = list(N_cases = 1)) %>%
      dplyr::mutate(yrsn = if_else(period == 1, 1, if_else(period == 2, 3, NA_integer_)), N_cases = N_cases/yrsn) %>%
      dplyr::select(country, st, period, N_cases),
    
    is_car %>%
      dplyr::mutate(st = if_else(st == "15BC", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
      dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, 1), country = "Israel", yearc = yearc) %>%
      group_by(country, yearc, st) %>%
      tally() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(period = if_else(yearc == 2010, 1L, if_else(yearc >= 2014 & yearc <= 2016, 2L, NA_integer_))) %>%
      dplyr::filter(!is.na(period)) %>%
      dplyr::group_by(country, period, st) %>%
      dplyr::summarize(N = sum(n)) %>%
      dplyr::mutate(NN = sum(N), carr_est = N/NN) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(st)) %>%
      tidyr::complete(st, period, country, fill = list(carr_est = 0.0001)) %>%
      dplyr::mutate(yrsn = if_else(period == 1, 1, if_else(period == 2, 3, NA_integer_)), carr_est = carr_est/yrsn) %>%
      dplyr::select(country, st, period, carr_est)) %>%
  
    dplyr::mutate(N_cases = if_else(is.na(N_cases), min(N_cases, na.rm = T), N_cases),
                carr_est = if_else(is.na(carr_est), min(carr_est, na.rm = T), carr_est),
                period = factor(period))

#====================================================================
#model validation in South Africa
#====================================================================

sa_val <- 
  sa_cari %>%
  dplyr::mutate(period = if_else(yearc >=2009 & yearc <=2011, 1L,
                                 if_else(yearc >=2015 & yearc <=2019, 2L, NA_integer_))) %>%
  dplyr::filter(!is.na(period)) %>%
  dplyr::group_by(period, st) %>%
  dplyr::summarise(N_cases = sum(n)) %>%
  dplyr::mutate(yrsn = if_else(period == 1, 3, if_else(period == 2, 5, NA_integer_)), N_cases = N_cases/yrsn) %>%
  dplyr::ungroup() %>%
  tidyr::complete(st, period, fill = list(N_cases = 0)) %>%
  dplyr::left_join(st_inv, by = 'st') %>%
  dplyr::mutate(inv = if_else(is.na(inv), inv_nvt, inv)) %>% #replace mean invasiveness for all NVTs
  dplyr::group_by(period) %>%
  dplyr::mutate(carr_est = (N_cases+1)/inv, 
                carr_est = carr_est/sum(carr_est, na.rm=T),
                period = as.factor(period),
                country = "South Africa") %>%  #inferred carriage
  dplyr::ungroup() %>%
  dplyr::select(country, st, period, N_cases, carr_est)

#====================================================================
#combining data
#====================================================================

#using post-pcv7 (pre-pcv13) carriage and IPD
effVE <-
  bind_rows(en_val, sa_val, is_val) %>%
  dplyr::filter(period==1) %>%
  dplyr::rename("prev1"="carr_est", "N1_cases"="N_cases") %>%
  dplyr::mutate(newVT = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L))

#check VTs versus NVTs
effVE %>% group_by(country, newVT) %>% tally(prev1)

#match vaccine effectiveness against each serotype disease and carriage dataset
effVE <-
  effVE %>%
  left_join(st_VEd %>% dplyr::select(st, VEd_pcv13)) %>%
  replace(is.na(.), 0) %>%
  left_join(st_VEc %>% dplyr::select(st, VEc_pcv13)) %>%
  replace(is.na(.), 0)

#configure proportion of NVT that needs to fill the hole in the niche
#set the relative change in the proportion of new vaccine coverage
gamma = 1
vax_uptake = 0.90

#determined whether there will be serotype expansion or not
#set 'a' to be true if NVT and false if VT
effVE <-
  effVE %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(a =((1-VEc_pcv13) >1) | ((1-VEc_pcv13) ==1 & newVT !=1)) %>%
  dplyr::ungroup()

#computed change in colonization post-PCV switch (#NVTs expand and VT shrink or remain constant)
effVE <-
  effVE %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(RR_col = if_else(a == TRUE, (1+((gamma*sum(prev1*VEc_pcv13*(1-a)))/sum(a*prev1))), sum(prev1*VEc_pcv13*(1-a)))) %>%
  dplyr::ungroup()

#compute the expected IPD incidence of each serotype in the new PCV era after switch
effVE <- 
  effVE %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(P_cases = N1_cases*(1-VEd_pcv13*vax_uptake)*(RR_col*vax_uptake)) %>%
  dplyr::ungroup() %>%
  dplyr::select(country, st, N1_cases, P_cases)

#check total ipd (all serotypes) before and after PCV switch
effVE %>% dplyr::group_by(country) %>% dplyr::summarise(sum(N1_cases)) %>% dplyr::ungroup()
effVE %>% dplyr::group_by(country) %>% dplyr::summarise(sum(P_cases)) %>% dplyr::ungroup()

#match the prediction dataset with observed dataset during post PCV switch
effVE <-
  bind_rows(en_val, sa_val, is_val) %>%
  dplyr::filter(period==2) %>%
  dplyr::rename("prev2"="carr_est", "N2_cases"="N_cases") %>%
  dplyr::select(country, st, N2_cases) %>%
  dplyr::left_join(effVE)

#====================================================================
#plotting the predictions and observed data data
#====================================================================

#combine the data for plotting
effVE <-
  bind_rows(
    effVE %>%
      dplyr::select(country, st, N2_cases, N1_cases) %>%
      dplyr::rename("cases"="N1_cases") %>%
      dplyr::mutate(cat = "Observed pre-PCV13 IPD"),

    effVE %>%
      dplyr::select(country, st, N2_cases, P_cases) %>%
      dplyr::rename("cases"="P_cases") %>%
      dplyr::mutate(cat = "Predicted post-PCV13 IPD")) %>%

  #plot only NVTs
  dplyr::mutate(newVT = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L))

A <- 
  effVE %>%
  dplyr::mutate(country = factor(country, levels = c("England", "Israel", "South Africa"))) %>%
  dplyr::group_by(country, cat) %>%
  dplyr::mutate(r = round((stats::cor(N2_cases, cases, method = "pearson", use = 'pairwise.complete.obs'))[1], digits = 2),
                w = round((stats::cor(N2_cases, cases, method = "spearman", use = 'pairwise.complete.obs'))[1], digits = 2),
                e = round(ie2misc::mae(N2_cases, cases)/ie2misc::madstat(N2_cases), digits = 3)) %>%
  dplyr::ungroup() %>%
  ggplot(aes(x = log(N2_cases+1), y = log(cases+1), label = st, color = st)) + 
  geom_point(aes(shape = factor(newVT)), size = 4, stroke = 1.5) + 
  facet_grid(cat ~ country, scales = "free_y") + 
  scale_x_continuous(breaks = seq(0, 6, 1), limits = c(0, 6)) + 
  scale_y_continuous(breaks = seq(0, 6, 1), limits = c(0, 6)) + 
  geom_abline(linetype = "dashed", size = 0.5) + 
  scale_shape_manual(values=c(1, 4)) + 
  geom_text(aes(x = log(N2_cases+1), y = log(cases+1), label = st, color = st),  size = 3.5, fontface = "bold", vjust = 0, hjust = -0.4) + 
  geom_text(aes(x = 0, y = 6, label = paste0("r = ", r)), hjust = 0, vjust = 0.5, size = 4, family = "serif", color = "gray20") +
  geom_text(aes(x = 0, y = 5.75, label = paste0("ω = ", w)), hjust = 0, vjust = 0.5, size = 4, family = "serif", color = "gray20") +
  geom_text(aes(x = 0, y = 5.5, label = paste0("ε = ", e)), hjust = 0, vjust = 0.5, size = 4, family = "serif", color = "gray20") +
  theme_bw(base_size = 20, base_family = "American typewriter") +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 22), strip.background=element_rect(fill="gray90")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 3)) +
  labs(title = "", x = "Observed post-PCV13 IPD log_cases", y = "Invasive pneumococcal disease (IPD) log_cases") +
  theme(legend.position = "none")

ggsave(here("output", "model_validation.png"),
       plot = A,
       width = 16, height = 11, unit = "in", dpi = 300) #width = 16, height = 13 | width = 24, height = 16
