#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
#model validation in Israel
#====================================================================

#input external datasets (USA and RSA from monocle GPS project)
gps_rsa <- 
  rio::import(here("data", "monocle_rsa.csv")) %>%
  dplyr::mutate(agey = if_else(is.na(Age_years), Age_months/12, Age_years)) %>%
  dplyr::filter(agey <=5, Phenotypic_serotype != "NT") %>%
  dplyr::rename("year" = "Year", "st" = "Phenotypic_serotype", "pcvera" = "Vaccine_period") %>%
  dplyr::select(year, st, pcvera) %>%
  dplyr::group_by(st) %>%
  dplyr::tally() %>% 
  dplyr::rename("n_rsa" = "n")#already #annual

gps_usa <- 
  rio::import(here("data", "monocle_usa.csv")) %>%
  dplyr::mutate(agey = if_else(is.na(Age_years), Age_months/12, Age_years)) %>%
  dplyr::filter(agey <=5, Phenotypic_serotype != "NT") %>%
  dplyr::rename("year" = "Year", "st" = "Phenotypic_serotype", "pcvera" = "Vaccine_period") %>%
  dplyr::select(year, st, pcvera)  %>%
  dplyr::group_by(st) %>%
  tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", st)) %>%
  dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
  tidyr::separate_rows(., st) %>% #tidy serotypes to atomic values
  dplyr::group_by(st) %>%
  dplyr::tally(n) %>%
  dplyr::ungroup() %>%
  dplyr::rename("n_usa" = "n") %>% 
  dplyr::mutate(n_usa = n_usa/3) #annualise

#using carriage prevalence
gps_car <-
  is_car %>% 
  dplyr::filter(datec <= date('2010-06-01')) %>%
  dplyr::group_by(st) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::add_row(st = c("2", "11B", "46", "36", "7B", "29", "24B", "31", "33A", "35F", "9A", "6D"), n = rep(0.5, 12)) %>% #either isolated in post-PCV13 IPD or pre-PCV13 IPD and not sampled in pre-PCV13 carriage
  dplyr::mutate(st = if_else(is.na(st), "NA", st), st = if_else(st == "15BC", "15B/15C", st)) %>%
  dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
  tidyr::separate_rows(., st) %>%
  dplyr::group_by(st) %>%
  dplyr::tally(n*(12/8)) %>% #annualise from 8mos
  dplyr::ungroup() %>%
  dplyr::mutate(prev0 = n/sum(n)) %>%
  dplyr::filter(st !="NA") %>%
  dplyr::rename("ncarr0" = "n") 

#join IPD data over the same time period 2013-16
gps_ipd <-
  is_ipd %>% 
  dplyr::filter(datec <= date('2010-06-01')) %>%
  dplyr::group_by(st) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", st)) %>%
  dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
  tidyr::separate_rows(., st) %>%
  dplyr::group_by(st) %>%
  dplyr::tally(n*(12/8)) %>%#annualise from 8mos
  dplyr::ungroup() %>%
  dplyr::rename("ipd0" = "n") %>%
  dplyr::add_row(st = c("10A", "11A", "13", "16F", "17F", "18A", "23B", "24B", "24F", "29", "31", "33A", "34", "35F", "36", "6C", "7B", "9N"), ipd0 = rep(1, 18)) #add serotypes found in post-PCV13 IPD

#combine all datasets
gps_fit <-
  gps_ipd %>% 
  left_join(st_inv %>% dplyr::select(st, inv)) %>%
  mutate(inv = if_else(is.na(inv), inv_nvt, inv)) %>%
  left_join(gps_car) %>%
  left_join(gps_usa) %>%
  left_join(gps_rsa) %>%
  mutate(ncarr0 = if_else(is.na(ncarr0), 0L, ncarr0), 
         n_usa = if_else(is.na(n_usa), 0L, n_usa),
         n_rsa = if_else(is.na(n_rsa), 0L, n_rsa))

#fit a poisson mixed model to infer post-PCV7 IPD in Israel
gps_model <- glmer(ipd0 ~ log(inv*(ncarr0+0.5)) + log(n_usa+0.5) + log(n_rsa+0.5) + (1|st), 
                family = poisson(link = "log"), 
                na.action = na.omit, 
                data = gps_fit, 
                verbose = TRUE)
coef(summary(gps_model))

#append IPD predictions to the dataset
gps_fit <- 
  gps_fit %>%
  dplyr::select(st) %>%
  dplyr::mutate(ipd0 = predict(gps_model, type = "response"))

#====================================================================

#using post-pcv7 (pre-pcv13) carriage prevalence
is_effVE <-
  gps_car %>%
  mutate(newVT = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L))

#check VTs versus NVTs
is_effVE %>% group_by(newVT) %>% tally(prev0)

#match vaccine effectiveness against each serotype disease and carriage dataset
is_effVE <-
  is_effVE %>%
  left_join(st_VEd %>% dplyr::select(st, VEd_pcv13)) %>%
  replace(is.na(.), 0)

is_effVE <-
  is_effVE %>%
  left_join(st_VEc %>% dplyr::select(st, VEc_pcv13)) %>%
  replace(is.na(.), 0)

#configure proportion of NVT that needs to fill the hole in the niche
#set the relative change in the proportion of new vaccine uptake
gamma = 1
vax_uptake = 1

#determined whether there will be serotype expansion or not
#set 'a' to be true if NVT and false if VT
is_effVE <-
  is_effVE %>%
  dplyr::mutate(a = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT !=1)))

#computed change in colonization post-PCV switch
is_effVE <-
  is_effVE %>%
  dplyr::mutate(RR_col = dplyr::if_else(a == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*vax_uptake*(1-a)))/sum(a*prev0))), #NVTs will expand to fill the hole in the NP niche
                                        sum(prev0*VEc_pcv13*vax_uptake*(1-a)))) #VTs expand if relative VE is lower for newer PCV than old PCV (otherwise remain same or decline)

#join IPD data over the same time period 2013-16
is_effDI <-
  gps_fit %>%
  left_join(is_effVE)

#compute the expected IPD incidence of each serotype in the new PCV era after switch
is_effDI <- 
  is_effDI %>%
  dplyr::mutate(ipd1 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (RR_col))

#check total ipd (all serotypes) before and after PCV switch
sum(is_effDI$ipd1)
sum(is_effDI$ipd0)

#====================================================================

#observed post-pcv13 ipd combined with observed pre-pcv ipd and predicted post-pcv ipd
XX <- 
  is_ipd %>%
  dplyr::mutate(yearc = str_c("y", yearc)) %>%
  base::split(list(.$yearc))

XX <- 
  dplyr::bind_rows(bind_rows(XX$y2014, XX$y2015, XX$y2016), .id = "yearc") %>%
  dplyr::mutate(yearg = "2014-16") %>%
  dplyr::group_by(yearg, st) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", st)) %>%
  dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
  tidyr::separate_rows(., st) %>% 
  dplyr::group_by(yearg, st) %>%
  dplyr::tally(n/3) %>%#annualise
  dplyr::ungroup() %>%
  dplyr::rename("ipdObs" = "n") %>%
  left_join(is_effDI) %>%
  mutate(vtg = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "VT", "NVT"))

#====================================================================
#model validation in South Africa
#====================================================================

#create ipd serotype dataset to match those of invasiveness
#infer carriage data (carriage  <- ipd / invasiveness)
sa_effVE <-
  sa_ipd %>%
  dplyr::filter(yearc == 2010 | yearc == 2011) %>%
  dplyr::select(yearc, st) %>%
  dplyr::group_by(st) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::rename("ipd" = "n") %>%
  dplyr::add_row(st = c("12B", "20", "23A", "23B", "24", "27", "28F", "29", "6C", "6D", "15F", "18F", "28A", "35A"), ipd = rep(1,14)) %>% #add serotypes found in post-PCV13 IPD
  #dplyr::add_row(st = c("10F", "18B", "18F", "19B", "20", "24", "25A", "28A", "28F", "33D", "35A", "37", "37C", "12B", "15F", "45", "10B", "11B", "27", "7B", "9A", "27", "6D"), ipd = rep(1,23)) %>% #add serotypes found in post-PCV13 IPD
  mutate(ipd = if_else(str_length(st)>3, ipd/2, ipd)) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd) %>%
  
  #split multiple serotypes with "/" into new rows
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "C", "15C", st)) %>%
  group_by(st) %>%
  summarise(ipd0 = sum(ipd)/2) %>%#annualise from 2yrs
  ungroup() %>%
  
  #join IPD serotypes with invasiveness
  left_join(st_inv %>% dplyr::select(st, inv)) %>%
  mutate(inv = if_else(is.na(inv), inv_nvt, inv)) %>%
  
#fill the NAs on serotype group
  mutate(prev0 = ipd0/inv, scalex = 0.5/sum(prev0), prev0 = scalex*prev0) %>%
  dplyr::select(st, ipd0, prev0) %>%
  mutate(newVT = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L))

#check VTs versus NVTs
sa_effVE %>% group_by(newVT) %>% tally(prev0)

#match vaccine effectiveness against each serotype diasese and carriage dataset
sa_effVE <-
  sa_effVE %>%
  left_join(st_VEd %>% dplyr::select(st, VEd_pcv13)) %>%
  replace(is.na(.), 0)

sa_effVE <-
  sa_effVE %>%
  left_join(st_VEc %>% dplyr::select(st, VEc_pcv13)) %>%
  replace(is.na(.), 0)

#configure proportion of NVT that needs to fill the hole in the niche
#set the relative change in the proportion of new vaccine uptake
gamma = 1
vax_uptake = 1

#determined whether there will be serotype expansion or not
#set 'a' to be true if NVT and false if VT
sa_effVE <-
  sa_effVE %>%
  dplyr::mutate(a = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT !=1)))

#computed change in colonization post-PCV switch
#estimated_x = 1 + ((gamma *sum(prev0*VEc_pcv13(1-a)))/(sum(a*prev0)))
sa_effVE <-
  sa_effVE %>%
  dplyr::mutate(RR_col = dplyr::if_else(a == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*vax_uptake*(1-a)))/sum(a*prev0))), #NVTs will expand to fill the hole in the NP niche
                                        sum(prev0*VEc_pcv13*vax_uptake*(1-a)))) #VTs expand if relative VE is lower for newer PCV than old PCV (otherwise remain same or decline)

#compute the expected IPD incidence of each serotype in the new PCV era after switch
sa_effDI <- 
  sa_effVE %>%
  dplyr::mutate(ipd1 = ipd0 * (1-VEd_pcv13) * (RR_col) * vax_uptake)

#check total ipd (all serotypes) before and after PCV switch
sum(sa_effDI$ipd1)
sum(sa_effDI$ipd0)

#observed post-pcv13 ipd combined with observed pre-pcv ipd and predicted post-pcv ipd
YY <- 
  sa_ipd %>%
  dplyr::select(yearc, st) %>%
  dplyr::mutate(yearc = str_c("y", yearc)) %>%
  base::split(list(.$yearc))

YY <- 
  dplyr::bind_rows(bind_rows(YY$y2015, YY$y2016, YY$y2017, YY$y2018, YY$y2019), .id = "yearc") %>%
  dplyr::mutate(yearg = "2015-19") %>%
  
  dplyr::group_by(yearg, st) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", st), st = if_else(st == "12B/F", "12B/12F", st)) %>%
  dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
  tidyr::separate_rows(., st) %>% 
  dplyr::group_by(yearg, st) %>%
  dplyr::tally(n/5) %>%#annualise from 5yrs
  dplyr::ungroup() %>%
  dplyr::rename("ipdObs" = "n") %>%
  left_join(sa_effDI) %>%
  mutate(vtg = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "VT", "NVT")) %>%
  dplyr::ungroup()

#====================================================================
#combined plots for Israel and South Africa
#====================================================================

XX_a <- 
  XX %>%
  dplyr::mutate(s = round((stats::cor(ipdObs, ipd0, method = c("spearman")))[1], digits = 3),
                p = round((stats::cor(ipdObs, ipd0, method = c("pearson")))[1], digits = 3),
                smape = round(Metrics::smape(ipdObs, ipd0)/2, digits = 3)) %>%
  
  ggplot() +
  geom_point(aes(x = log(ipd0), y = log(ipdObs), color = vtg), stroke = 2, size = 1, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = -2.7, y = 4.6, label = paste0("Spearman = ", s)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = -2.7, y = 4, label = paste0("Pearson = ", p)), size = 7, family = "American typewriter", color = "black") +
  #geom_text(aes(x = -2.7, y = 3.4, label = paste0("sMAPE = ", smape)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = log(ipd0), y = log(ipdObs), label = st, color = vtg),  size = 3, fontface = "bold", vjust = 0, hjust = -0.4) +
  scale_x_continuous(limit = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_y_continuous(limit = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) + 
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "ISRAEL", x = "observed pre-PCV13 log_ipd (2010)", y = "observed post-PCV13 log_ipd (2014-16)") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.9, 0.1), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20), strip.background = element_rect(fill = "gray80"))

XX_b <- 
  XX %>%
  dplyr::mutate(s = round((stats::cor(ipdObs, ipd1, method = c("spearman")))[1], digits = 3),
                p = round((stats::cor(ipdObs, ipd1, method = c("pearson")))[1], digits = 3),
                smape = round(Metrics::smape(ipdObs, ipd1)/2, digits = 3)) %>%
  
  ggplot() +
  geom_point(aes(x = log(ipd1), y = log(ipdObs), color = vtg), stroke = 2, size = 1, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = -2.7, y = 4.6, label = paste0("Spearman = ", s)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = -2.7, y = 4, label = paste0("Pearson = ", p)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = -2.7, y = 3.4, label = paste0("sMAPE = ", smape)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = log(ipd1), y = log(ipdObs), label = st, color = vtg),  size = 3, fontface = "bold", vjust = 0, hjust = -0.4) +
  scale_x_continuous(limit = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_y_continuous(limit = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) + 
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "predicted post-PCV13 log_ipd (2014-16)", y = "observed post-PCV13 log_ipd (2014-16)") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.9, 0.1), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20), strip.background = element_rect(fill = "gray80"))

YY_a <- 
  YY %>%
  dplyr::mutate(s = round((stats::cor(ipdObs, ipd0, method = c("spearman")))[1], digits = 3),
                p = round((stats::cor(ipdObs, ipd0, method = c("pearson")))[1], digits = 3),
                smape = round(Metrics::smape(ipdObs, ipd0)/2, digits = 3)) %>%
  
  ggplot() +
  geom_point(aes(x = log(ipd0), y = log(ipdObs), color = vtg), stroke = 2, size = 1, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = -2.7, y = 4.6, label = paste0("Spearman = ", s)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = -2.7, y = 4, label = paste0("Pearson = ", p)), size = 7, family = "American typewriter", color = "black") +
  #geom_text(aes(x = -2.7, y = 3.4, label = paste0("sMAPE = ", smape)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = log(ipd0), y = log(ipdObs), label = st, color = vtg),  size = 3, fontface = "bold", vjust = 0, hjust = -0.4) +
  scale_x_continuous(limit = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_y_continuous(limit = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) + 
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "SOUTH AFRICA", x = "observed pre-PCV13 log_ipd (2010-11)", y = "observed post-PCV13 log_ipd (2015-19)") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.9, 0.1), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20), strip.background = element_rect(fill = "gray80"))

YY_b <- 
  YY %>%
  dplyr::mutate(s = round((stats::cor(ipdObs, ipd1, method = c("spearman")))[1], digits = 3),
                p = round((stats::cor(ipdObs, ipd1, method = c("pearson")))[1], digits = 3),
                smape = round(Metrics::smape(ipdObs, ipd1)/2, digits = 3)) %>%

  ggplot() +
  geom_point(aes(x = log(ipd1), y = log(ipdObs), color = vtg), stroke = 2, size = 1, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = -2.7, y = 4.6, label = paste0("Spearman = ", s)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = -2.7, y = 4, label = paste0("Pearson = ", p)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = -2.7, y = 3.4, label = paste0("sMAPE = ", smape)), size = 7, family = "American typewriter", color = "black") +
  geom_text(aes(x = log(ipd1), y = log(ipdObs), label = st, color = vtg),  size = 3, fontface = "bold", vjust = 0, hjust = -0.4) +
  scale_x_continuous(limit = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_y_continuous(limit = c(-5, 5), breaks = seq(-5, 5, 1)) +
  scale_color_manual(values = c("#F8766D", "#00BFC4")) + 
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "predicted post-PCV13 log_ipd (2015-19)", y = "observed post-PCV13 log_ipd (2015-19)") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.9, 0.1), legend.title = element_blank()) +
  theme(strip.text.x = element_text(size = 20), strip.text.y = element_text(size = 20), strip.background = element_rect(fill = "gray80"))


#save combined plots
ggsave(here("output", "modelValidity.png"),
       plot = (XX_a/XX_b | YY_a/YY_b), 
       width = 12, height = 14, unit = "in", dpi = 300)
