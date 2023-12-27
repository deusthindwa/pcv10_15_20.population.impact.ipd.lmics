#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries
 
#====================================================================
#INVASIVENESS DATA
#====================================================================

#import invasiveness data
invasivenes <-
  rio::import(here("data", "invasiveness_global2.csv")) %>%
  dplyr::filter(st != "NT") %>%
  mutate(pcv7pfz = if_else(grepl("\\b(4|6A|6B|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV7", "NVT"), #add 6A cross-protection
         pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10-sii", "NVT"),
         pcv10gsk = if_else(grepl("\\b(1|4|5|6A|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV10-gsk", "NVT"), #add 6A cross-protection
         pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
         pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "NVT"),
         pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "NVT"),
         exp.inv = as.numeric(exp(log.inv)),
         log.inv = as.numeric(log.inv)) %>%
  dplyr::select(st, log.inv, exp.inv, everything(), -log.inv.prec)

#calculate weighted mean invasiveness of VT and NVT for serotypes with unknown invasiveness
inv_pcv20pfz <- invasivenes %>% group_by(pcv20pfz) %>% mutate(wgt.inv = weighted.mean(exp.inv)) %>% ungroup() %>% distinct(pcv20pfz, wgt.inv)
inv_pcv15mek <- invasivenes %>% group_by(pcv15mek) %>% mutate(wgt.inv = weighted.mean(exp.inv)) %>% ungroup() %>% distinct(pcv15mek, wgt.inv)
inv_pcv13pfz <- invasivenes %>% group_by(pcv13pfz) %>% mutate(wgt.inv = weighted.mean(exp.inv)) %>% ungroup() %>% distinct(pcv13pfz, wgt.inv)
inv_pcv10sii <- invasivenes %>% group_by(pcv10sii) %>% mutate(wgt.inv = weighted.mean(exp.inv)) %>% ungroup() %>% distinct(pcv10sii, wgt.inv)
inv_pcv10gsk <- invasivenes %>% group_by(pcv10gsk) %>% mutate(wgt.inv = weighted.mean(exp.inv)) %>% ungroup() %>% distinct(pcv10gsk, wgt.inv)
inv_pcv7pfz  <- invasivenes %>% group_by(pcv7pfz) %>% mutate(wgt.inv = weighted.mean(exp.inv)) %>% ungroup() %>% distinct(pcv7pfz, wgt.inv)
inv_vt = (inv_pcv20pfz$wgt.inv[1] + inv_pcv15mek$wgt.inv[1] + inv_pcv13pfz$wgt.inv[1] + inv_pcv10sii$wgt.inv[1] + inv_pcv10gsk$wgt.inv[1] + inv_pcv7pfz$wgt.inv[2])/6 #average VT invasiveness
inv_nvt = (inv_pcv20pfz$wgt.inv[2] + inv_pcv15mek$wgt.inv[2] + inv_pcv13pfz$wgt.inv[2] + inv_pcv10sii$wgt.inv[2] + inv_pcv10gsk$wgt.inv[2] + inv_pcv7pfz$wgt.inv[1])/6 #average NVT invasiveness

#====================================================================
#PREDICTIONS FORM SOUTH AFRICA
#====================================================================

#create ipd serotype dataset to match those of invasiveness
#infer carriage data pre-pcv13 introduction in South Africa (carriage  <- ipd / invasiveness)
sa_carb2009 <-
  sa_ipdb2009 %>%
  dplyr::select(yearc, st) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd" = "n") %>%
  mutate(ipd = if_else(str_length(st)>3, ipd/2, ipd)) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd) %>%
  
  #split multiple serotypes with "/" into new rows
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "C", "15C", st)) %>%
  group_by(st) %>%
  summarise(ipd = sum(ipd)) %>% 
  ungroup() %>%
  
  #join IPD serotypes with invasiveness
  left_join(invasivenes) %>%
  mutate(exp.inv = if_else(is.na(exp.inv), inv_nvt, exp.inv)) %>%
  mutate(log.inv = log(exp.inv)) %>%
  
  #fill the NAs on serotype group
  mutate(pcv7pfz = if_else(grepl("\\b(4|6A|6B|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV7", "NVT"), #add 6A cross-protection
         pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10-sii", "NVT"),
         pcv10gsk = if_else(grepl("\\b(1|4|5|6A|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV10-gsk", "NVT"), #add 6A cross-protection
         pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
         pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "NVT"),
         pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "NVT"),
         prev1 = ipd/exp.inv,
         scalex = 0.8/sum(prev1), #assume total prevalence is up to 80% Nzeze et al.
         prev2 = scalex*prev1) %>% #scale prevalence according
  dplyr::select(st, ipd, exp.inv, log.inv, prev1, prev2, everything())

#====================================================================

#create carriage serotype dataset for pre PCV13 introduction
sa_carb2011 <-
  sa_ipd %>%
  dplyr::filter(yearc >=2010, yearc <=2011) %>%
  dplyr::select(yearc, st) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd" = "n") %>%
  mutate(ipd = if_else(str_length(st)>3, ipd/2, ipd)) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd) %>%
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "C", "15C", st)) %>%
  group_by(st) %>%
  summarise(ipd = sum(ipd)) %>% 
  ungroup() %>%
  
  #join IPD serotypes with invasiveness
  left_join(invasivenes) %>%
  mutate(exp.inv = if_else(is.na(exp.inv), inv_nvt, exp.inv)) %>%
  mutate(log.inv = log(exp.inv)) %>%
  
  #fill the NAs on serotype group
  mutate(pcv7pfz = if_else(grepl("\\b(4|6A|6B|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV7", "NVT"), #add 6A cross-protection
         pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10-sii", "NVT"),
         pcv10gsk = if_else(grepl("\\b(1|4|5|6A|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV10-gsk", "NVT"), #add 6A cross-protection
         pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
         pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "NVT"),
         pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "NVT"),
         prev1 = ipd/exp.inv,
         scalex = 0.70/sum(prev1), #assume total prevalence is up to 80% Nzeze et al.
         prev2 = scalex*prev1) %>% #scale prevalence according
  dplyr::select(st, ipd, exp.inv, log.inv, prev1, prev2, everything())

#====================================================================

#create ipd serotype dataset to match those of invasiveness
#infer carriage data post-pcv13 introduction in south africa (carriage <- ipd / invasiveness)
sa_cara2015 <-
  sa_ipda2015 %>%
  dplyr::select(yearc, st) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd" = "n") %>%
  mutate(ipd = if_else(str_length(st)>3, ipd/2, ipd)) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd) %>%

  #split multiple serotypes with "/" into new rows
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "C", "15C", st), st = if_else(st == "F", "12F", st)) %>%
  group_by(st) %>%
  summarise(ipd = sum(ipd)) %>% 
  ungroup() %>%

  #join IPD serotypes with invasiveness
  left_join(invasivenes) %>%
  mutate(exp.inv = if_else(is.na(exp.inv), inv_nvt, exp.inv)) %>%
  mutate(log.inv = log(exp.inv)) %>%
  
  #fill the NAs on serotype group
  mutate(pcv7pfz = if_else(grepl("\\b(4|6B|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV7", "NVT"),
         pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10-sii", "NVT"),
         pcv10gsk = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV10-gsk", "NVT"),
         pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
         pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "NVT"),
         pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "NVT"),
         prev1 = ipd/exp.inv,
         scalex = 0.5/sum(prev1), #assume total prevalence is up to 50% Olwagen et al.
         prev2 = scalex*prev1) %>% #scale prevalence according
  dplyr::select(st, ipd, exp.inv, log.inv, prev1, prev2, everything())

#====================================================================

#plots of the relationship between serotype-specific carriage prevalence and ipd
A1 <-
ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv7pfz == "PCV7"), aes(x = log(ipd), y = prev2), color = "#8B0019", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv7pfz == "PCV7"), aes(x = log(ipd), y = prev2), color = "#8B0019", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV7 VT", x = "observed log_ipd", y = "inferred VT carriage prevalence") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16))

A2 <-
ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv7pfz == "NVT"), aes(x = log(ipd), y = prev2), color = "#8B0019", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv7pfz == "NVT"), aes(x = log(ipd), y = prev2), color = "#8B0019", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV7 NVT", x = "observed log_ipd", y = "inferred NVT carriage prevalence") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16))

B1 <-
ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv10sii == "PCV10-sii"), aes(x = log(ipd), y = prev2), color = "#D16A5F", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv10sii == "PCV10-sii"), aes(x = log(ipd), y = prev2), color = "#D16A5F", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV10-sii VT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

B2 <-
  ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv10sii == "NVT"), aes(x = log(ipd), y = prev2), color = "#D16A5F", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv10sii == "NVT"), aes(x = log(ipd), y = prev2), color = "#D16A5F", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV10-sii NVT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

C1 <-
ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv10gsk == "PCV10-gsk"), aes(x = log(ipd), y = prev2), color = "#E358A0", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv10gsk == "PCV10-gsk"), aes(x = log(ipd), y = prev2), color = "#E358A0", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV10-gsk VT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

C2 <-
  ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv10gsk == "NVT"), aes(x = log(ipd), y = prev2), color = "#E358A0", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv10gsk == "NVT"), aes(x = log(ipd), y = prev2), color = "#E358A0", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV10-gsk NVT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

D1 <-
ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv13pfz == "PCV13"), aes(x = log(ipd), y = prev2), color = "#8872D9", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv13pfz == "PCV13"), aes(x = log(ipd), y = prev2), color = "#8872D9", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV13 VT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

D2 <-
  ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv13pfz == "NVT"), aes(x = log(ipd), y = prev2), color = "#8872D9", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv13pfz == "NVT"), aes(x = log(ipd), y = prev2), color = "#8872D9", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV13 NVT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

E1 <-
ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv15mek == "PCV15"), aes(x = log(ipd), y = prev2), color = "#C04AE1", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv15mek == "PCV15"), aes(x = log(ipd), y = prev2), color = "#C04AE1", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV15 VT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

E2 <-
  ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv15mek == "NVT"), aes(x = log(ipd), y = prev2), color = "#C04AE1", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv15mek == "NVT"), aes(x = log(ipd), y = prev2), color = "#C04AE1", size = 2, stroke = 2, shape = 4) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV15 NVT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

F1 <-
ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv20pfz == "PCV20"), aes(x = log(ipd), y = prev2, shape = "pre-PCV13"), color = "#E3B45B", size = 2, stroke = 2) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv20pfz == "PCV20"), aes(x = log(ipd), y = prev2, shape = "post-PCV13"), color = "#E3B45B", size = 2, stroke = 2) +
  geom_text(data = dplyr::filter(sa_carb2009, pcv20pfz == "PCV20"), aes(x = log(ipd), y = prev2, label = st), color = "black", size = 3, fontface = "bold", vjust = 0, hjust = 0) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  scale_shape_manual(name='Period', breaks = c("pre-PCV13", "post-PCV13"), values=c("pre-PCV13" = 1, "post-PCV13" = 4))+
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV20 VT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

F2 <-
  ggplot() + 
  geom_point(data = dplyr::filter(sa_carb2009, pcv20pfz == "NVT"), aes(x = log(ipd), y = prev2), color = "#E3B45B", size = 2, stroke = 2, shape = 1) +
  geom_point(data = dplyr::filter(sa_cara2015, pcv20pfz == "NVT"), aes(x = log(ipd), y = prev2), color = "#E3B45B", size = 2, stroke = 2, shape = 4) +
  geom_text(data = dplyr::filter(sa_carb2009, pcv20pfz == "NVT"), aes(x = log(ipd), y = prev2, label = st), color = "black", size = 3, fontface = "bold", vjust = 0, hjust = 0) +
  scale_x_continuous(limit = c(1, 7), breaks = seq(1, 7, 1)) +
  scale_y_continuous(limit = c(0, 0.08), breaks = seq(0, 0.08, 0.02), labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "PCV20 NVT", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 0))

#save combined plots
ggsave(here("output", "sfig9_carstDist.png"),
       plot = (A1 | B1 | C1 | D1 | E1 | F1) / (A2 | B2 | C2 | D2 | E2 | F2), 
       width = 26, height = 14, unit = "in", dpi = 300)

#====================================================================

#plots of the relationship between serotype group-specific carriage prevalence and ipd
sa_pcv <-
  bind_rows(
    sa_carb2009 %>% group_by(pcv7pfz) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv7pfz, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV7") %>% mutate(reg = "PCV7", period = "pre-PCV"),
    
    sa_cara2015 %>% group_by(pcv7pfz) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv7pfz, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV7") %>% mutate(reg = "PCV7", period = "post-PCV"),
    
    sa_carb2009 %>% group_by(pcv10sii) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv10sii, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV10-sii") %>% mutate(reg = "PCV10-sii", period = "pre-PCV"),
    
    sa_cara2015 %>% group_by(pcv10sii) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv10sii, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV10-sii") %>% mutate(reg = "PCV10-sii", period = "post-PCV"),
    
    sa_carb2009 %>% group_by(pcv10gsk) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv10gsk, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV10-gsk") %>% mutate(reg = "PCV10-gsk", period = "pre-PCV"),
    
    sa_cara2015 %>% group_by(pcv10gsk) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv10gsk, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV10-gsk") %>% mutate(reg = "PCV10-gsk", period = "post-PCV"),
    
    sa_carb2009 %>% group_by(pcv13pfz) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv13pfz, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV13") %>% mutate(reg = "PCV13", period = "pre-PCV"),
    
    sa_cara2015 %>% group_by(pcv13pfz) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv13pfz, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV13") %>% mutate(reg = "PCV13", period = "post-PCV"),
    
    sa_carb2009 %>% group_by(pcv15mek) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv15mek, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV15") %>% mutate(reg = "PCV15", period = "pre-PCV"),
    
    sa_cara2015 %>% group_by(pcv15mek) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv15mek, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV15") %>% mutate(reg = "PCV15", period = "post-PCV"),
    
    sa_carb2009 %>% group_by(pcv20pfz) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv20pfz, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV20") %>% mutate(reg = "PCV20", period = "pre-PCV"),
    
    sa_cara2015 %>% group_by(pcv20pfz) %>% tally(prev2) %>% 
      pivot_wider(names_from = pcv20pfz, values_from = n) %>%
      ungroup() %>% rename("VT" = "PCV20") %>% mutate(reg = "PCV20", period = "post-PCV"))

A <-
  bind_rows(
    sa_pcv %>% dplyr::select(NVT, reg, period) %>% rename("car" = "NVT") %>% mutate(`serotype group` = "NVT"),
    sa_pcv %>% dplyr::select(VT, reg, period)%>% rename("car" = "VT") %>% mutate(`serotype group` = "VT")) %>%
  
  ggplot() +
  geom_col(aes(x = factor(reg,levels(factor(reg))[c(6,1,2,3,4,5)]), y=car, fill = `serotype group`), position = position_dodge()) +
  scale_y_continuous(limit = c(0, 0.6), breaks = seq(0, 0.6, 0.2), labels = scales::percent_format(accuracy = 1)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "", x = "PCV regimen", y = "inferred serotype group carriage prevalence") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  facet_grid(.~ factor(period, levels = c("pre-PCV", "post-PCV")))

#save combined plots
ggsave(here("output", "sfig10_carsgDist.png"),
       plot = (A), 
       width = 20, height = 7, unit = "in", dpi = 300)

#====================================================================

#model validation of inferred carriage data using pre- and post-PCV7 data
sa_ipdPred1 <- 
  sa_ipd %>% 
  dplyr::filter(yearc <2009) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd1" = "n") %>%
  mutate(ipd1 = if_else(str_length(st)>3, (ipd1/2+0.5), ipd1)) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd1) %>%
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "C", "15C", st)) %>%
  group_by(st) %>%
  summarise(ipd1 = sum(ipd1)/4) %>%  #4y surveillance from 2005-2008
  ungroup()

sa_ipdPred2 <- 
  sa_ipd %>% 
  dplyr::filter(yearc >=2010, yearc <=2011) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd2" = "n") %>%
  mutate(ipd2 = if_else(str_length(st)>3, (ipd2/2+0.5), ipd2)) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd2) %>%
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "C", "15C", st)) %>%
  group_by(st) %>%
  summarise(ipd2 = sum(ipd2)/2) %>%  #2y surveillance from 2010-2011
  ungroup()

sa_ipdPred3 <- 
  sa_ipd %>% 
  dplyr::filter(yearc >=2015) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd3" = "n") %>%
  mutate(ipd3 = if_else(str_length(st)>3, (ipd3/2+0.5), ipd3)) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd3) %>%
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "C", "15C", st), st = if_else(st == "F", "12F", st)) %>%
  group_by(st) %>%
  summarise(ipd3 = sum(ipd3)/5) %>%  #5y surveillance from 2015-2019
  ungroup()

sa_ipdPred <-
  sa_ipdPred1 %>%
  full_join(sa_ipdPred2) %>%
  full_join(sa_ipdPred3) %>%
  full_join(sa_carb2009 %>% dplyr::select(st, prev2) %>% rename("prev1" = "prev2")) %>%
  full_join(sa_carb2011 %>% dplyr::select(st, prev2)) %>%
  full_join(sa_cara2015 %>% dplyr::select(st, prev2) %>% rename("prev3" = "prev2")) %>%
  mutate(cond = if_else(is.na(ipd1) | is.na(ipd2) | is.na(ipd3) | is.na(prev1) | is.na(prev2) | is.na(prev3), NA_integer_, 1L)) %>%
  dplyr::filter(!is.na(cond)) %>%
  dplyr::select(everything(), -cond) %>%
  mutate(ipd1 = log(ipd1+1),
         ipd2 = log(ipd2+1),
         ipd3 = log(ipd3+1),
         prev1 = log(prev1),
         prev2 = log(prev2),
         prev3 = log(prev3),
         prevfR = prev2/prev1,
         prevfR = log(prevfR))

#fit a negative-binomial model to predict post-PCV7 IPD
#ipd_2010-2011 = ipd_2005-2008 * (carr_prev_2010-11/carr_prev_2005-2008)
sa_model1 <- tidy(MASS::glm.nb(ipd2 ~ ipd1 + prevfR, data = sa_ipdPred), 
                   exponentiate = TRUE, 
                   conf.int = TRUE, 
                   conf.level = 0.95)
sa_ipdPred <- 
  sa_ipdPred %>% 
  dplyr::mutate(ipdfit1 = sa_model1$estimate[1] + sa_model1$estimate[2]*exp(ipd1-1) + sa_model1$estimate[3]*exp(prevfR)) #exp(ipd1-1) to revert to original IPD values from log(ipd1+1)

rsq <- function (x, y) cor(x, y) ^ 2

#visualize relationship between carriage prevalence and IPD
A <- 
  sa_ipdPred %>%
  mutate(r = abs(round((stats::cor(ipd2, log(ipdfit1), method = c("pearson")))[1], digits = 3)),
         p = scientific(stats::cor.test(ipd2, log(ipdfit1))[3]$p.value),
         r2 = abs(round(rsq(ipd2, log(ipdfit1)), digits = 3))) %>%
  ggplot() +
  geom_point(aes(x = (ipd2), y = log(ipdfit1)), stroke = 2, size = 0.5, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = 1, y = 5.8, label = paste0("r = ", r[1])), size = 4, family = "American typewriter") +
  geom_text(aes(x = 1, y = 5.5, label = paste0("p = ", p[1])), size = 4, family = "American typewriter") +
  scale_x_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_y_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "A", x = "observed log_ipd", y = "predicted log_ipd") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = "none")

B <- 
  sa_ipdPred %>%
  mutate(pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10-sii", "NVT")) %>%
  group_by(pcv10sii) %>%
  mutate(r = abs(round((stats::cor(ipd2, log(ipdfit1), method = c("pearson")))[1], digits = 3)),
         p = scientific(stats::cor.test(ipd2, log(ipdfit1))[3]$p.value),
         r2 = abs(round(rsq(ipd2, log(ipdfit1)), digits = 3))) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x = (ipd2), y = log(ipdfit1), color = pcv10sii), stroke = 2, size = 0.5, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = 1, y = 5.8, label = paste0("r = ", r[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  geom_text(aes(x = 1, y = 5.5, label = paste0("p = ", p[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  scale_x_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_y_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "B", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.7,0.1), legend.title = element_blank())

C <- 
  sa_ipdPred %>%
  mutate(pcv10gsk = if_else(grepl("\\b(1|4|5|6A|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV10-gsk", "NVT")) %>% #add 6A cross-protection 
  group_by(pcv10gsk) %>%
  mutate(r = abs(round((stats::cor(ipd2, log(ipdfit1), method = c("pearson")))[1], digits = 3)),
         p = scientific(stats::cor.test(ipd2, log(ipdfit1))[3]$p.value),
         r2 = abs(round(rsq(ipd2, log(ipdfit1)), digits = 3))) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x = (ipd2), y = log(ipdfit1), color = pcv10gsk), stroke = 2, size = 0.5, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = 1, y = 5.8, label = paste0("r = ", r[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  geom_text(aes(x = 1, y = 5.5, label = paste0("p = ", p[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  scale_x_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_y_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "C", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.7,0.1), legend.title = element_blank())

D <- 
  sa_ipdPred %>%
  mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT")) %>% 
  group_by(pcv13pfz) %>%
  mutate(r = abs(round((stats::cor(ipd2, log(ipdfit1), method = c("pearson")))[1], digits = 3)),
         p = scientific(stats::cor.test(ipd2, log(ipdfit1))[3]$p.value),
         r2 = abs(round(rsq(ipd2, log(ipdfit1)), digits = 3))) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x = (ipd2), y = log(ipdfit1), color = pcv13pfz), stroke = 2, size = 0.5, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = 1, y = 5.8, label = paste0("r = ", r[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  geom_text(aes(x = 1, y = 5.5, label = paste0("p = ", p[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  scale_x_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_y_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "D", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.7,0.1), legend.title = element_blank())

E <- 
  sa_ipdPred %>%
  mutate(pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "NVT")) %>% 
  group_by() %>%
  mutate(r = abs(round((stats::cor(ipd2, log(ipdfit1), method = c("pearson")))[1], digits = 3)),
         p = scientific(stats::cor.test(ipd2, log(ipdfit1))[3]$p.value),
         r2 = abs(round(rsq(ipd2, log(ipdfit1)), digits = 3))) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x = (ipd2), y = log(ipdfit1), color = pcv15mek), stroke = 2, size = 0.5, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = 1, y = 5.8, label = paste0("r = ", r[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  geom_text(aes(x = 1, y = 5.5, label = paste0("p = ", p[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  scale_x_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_y_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "E", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.7,0.1), legend.title = element_blank())

F <- 
  sa_ipdPred %>%
  mutate(pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "NVT")) %>% 
  group_by() %>%
  mutate(r = abs(round((stats::cor(ipd2, log(ipdfit1), method = c("pearson")))[1], digits = 3)),
         p = scientific(stats::cor.test(ipd2, log(ipdfit1))[3]$p.value),
         r2 = abs(round(rsq(ipd2, log(ipdfit1)), digits = 3))) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(x = (ipd2), y = log(ipdfit1), color = pcv20pfz), stroke = 2, size = 0.5, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = 1, y = 5.8, label = paste0("r = ", r[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  geom_text(aes(x = 1, y = 5.5, label = paste0("p = ", p[1])), size = 4, family = "American typewriter", color = "#00BFC4") +
  scale_x_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  scale_y_continuous(limit = c(0, 6), breaks = seq(0, 6, 1)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "F", x = "observed log_ipd", y = "") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.7,0.1), legend.title = element_blank())

#save combined plots
ggsave(here("output", "sfig12_sa_ipdPredictions.png"),
       plot = (A | B | C | D | E | F), 
       width = 24, height = 7, unit = "in", dpi = 300)

#====================================================================
#PREDICTIONS FROM MALAWI
#====================================================================

mw_ipdb2011_obs <-
  mw_ipdb2011 %>%
  dplyr::select(yearc, st) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd" = "n") %>%
  mutate(ipd = if_else(str_length(st)>3 & str_length(st)<=7, ipd/2, 
                       if_else(str_length(st)>7, ipd/3, ipd))) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd, ipd) %>%
  
  #split multiple serotypes with "/" into new rows
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "23B1", "23B", st), 
         st = if_else(st == "12FAB", "12F", st),
         ipd = if_else(st == "12F", 1, ipd)) %>%
  group_by(st) %>%
  summarise(ipd = sum(ipd)) %>% 
  ungroup()

#====================================================================

mw_ipda2015_obs <-
  mw_ipda2015 %>%
  dplyr::select(yearc, st) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd" = "n") %>%
  mutate(ipd = if_else(str_length(st)>3 & str_length(st)<=7, ipd/2, 
                       if_else(str_length(st)>7, ipd/3, ipd))) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, ipd) %>%
  
  #split multiple serotypes with "/" into new rows
  tidyr::separate_rows(., st) %>%
  group_by(st) %>%
  summarise(ipd = sum(ipd)) %>% 
  ungroup()

#====================================================================

mw_cara2015_obs <-
  mw_cara2015 %>%
  dplyr::select(yearc, st) %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  rename("prev" = "n") %>%
  #dplyr::filter(st != "None") %>%
  mutate(prev = if_else(str_length(st)>3 & str_length(st)<=7, prev/2, 
                        if_else(str_length(st)>7, prev/3, prev))) %>% #multiple serotypes in a sample, split into half
  dplyr::select(st, prev) %>%
  
  #split multiple serotypes with "/" into new rows
  tidyr::separate_rows(., st) %>%
  mutate(st = if_else(st == "23B1", "23B", st)) %>%
  group_by(st) %>%
  summarise(prev = sum(prev)) %>% 
  ungroup() %>%
  mutate(prev = if_else(st == "None", 2*prev, prev),#None was wrongly halved
         prev = prev/sum(prev))

#====================================================================

#2015-2019 samples summary
mw_none = 1004 #no samples
mw_nvt = 1110 #other non-vaccine samples

#create ipd serotype dataset to match those of invasiveness
#infer carriage data pre-pcv13 introduction in south africa (carriage  <- ipd / invasiveness)
mw_ipda2015_pred <-
  mw_cara2015_obs %>%

  #join IPD serotypes with invasiveness
  left_join(invasivenes) %>%
  mutate(exp.inv = if_else(is.na(exp.inv), inv_nvt, exp.inv),
         exp.inv = if_else(st == "None", 0, exp.inv)) %>%
  mutate(log.inv = log(exp.inv)) %>%
  
  #fill the NAs on serotype group
  mutate(pcv7pfz = if_else(grepl("\\b(4|6A|6B|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV7", "NVT"), #add 6A cross-protection
         pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10-sii", "NVT"),
         pcv10gsk = if_else(grepl("\\b(1|4|5|6A|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV10-gsk", "NVT"), #add 6A cross-protection
         pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
         pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "NVT"),
         pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "NVT"),
         ipd = prev*exp.inv,
         scalex = sum(mw_ipda2015_obs$ipd)/sum(ipd), # sum(mw_ipda2015_obs$ipd) = assume total ipd 97 based on reported IPD
         ipd = scalex*ipd) %>% #scale ipd according
  dplyr::select(st, prev, exp.inv, ipd, everything(), -log.inv)

#====================================================================

#use cleaned observed pre-pcv13 IPD dataset intro in Malawi
mw_carb2011_pred <-
  mw_ipdb2011_obs %>%

  #join IPD serotypes with invasiveness
  left_join(invasivenes) %>%
  mutate(exp.inv = if_else(is.na(exp.inv), inv_pcv13pfz$wgt.inv[2], exp.inv)) %>%
  mutate(log.inv = log(exp.inv),
         pcv13pfz = if_else(is.na(pcv13pfz), "NVT", pcv13pfz)) %>%

  #fill the PCV13 serotype group since only group with complete serotypes
  mutate(prev = ipd/exp.inv,
         scalex = 0.45/sum(prev), #assume total prevalence is up to 45%, Ellen et al.
         prev = scalex*prev) %>% #scale prev according
  dplyr::select(st, prev, exp.inv, ipd, pcv13pfz)

#====================================================================

#generate ppotential pre-PCV carriage data
#postpcvIPD = prepcvIPD * postpcvCarr / prepcvCarr
#prepcvCarr = prepcvIPD * postpcvCarr / postpcvIPD
# 
# bind_cols(
#   mw_ipdb2011_obs %>%
#     mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT")) %>%
#     group_by(pcv13pfz) %>%
#     tally(ipd) %>%
#     ungroup() %>%
#     rename("ipd2011" = "n"),
#   
#   mw_ipda2015_obs %>%
#     mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT")) %>%
#     group_by(pcv13pfz) %>%
#     tally(ipd) %>%
#     ungroup() %>%
#     rename("ipd2015" = "n") %>%
#     dplyr::select(ipd2015),
#   
#   mw_cara2015_obs %>%
#     mutate(pcv13pfz = if_else(st == "None", "None", 
#                               if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"))) %>%
#     group_by(pcv13pfz) %>%
#     tally(prev) %>%
#     ungroup() %>%
#     rename("car2015" = "n") %>%
#     dplyr::filter(pcv13pfz != "None") %>%
#     dplyr::select(car2015)) %>%
#   
#   mutate(car2011 = ipd2011 * car2015 / ipd2015) %>%
#   
#   mutate(scalex = 0.70/sum(car2011), #assume total prevalence is up to 45%, Ellen et al.
#          car2011 = scalex*car2011) #scale prev according

A <-
  left_join(
    mw_ipda2015_obs %>%
      mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
             st = if_else(pcv13pfz == "NVT", "NVT", st)) %>%
      group_by(st) %>%
      tally(ipd) %>% 
      rename("n1"="n"),
    
    mw_ipda2015_pred %>%
      mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
             st = if_else(pcv13pfz == "NVT", "NVT", st)) %>%
      group_by(st) %>%
      tally(ipd) %>% 
      rename("n2"="n")) %>%
  mutate(r = abs(round((stats::cor(n1, n2, method = c("pearson")))[1], digits = 3)),
         p = scientific(stats::cor.test(n1, n2)[3]$p.value),
         r2 = abs(round(rsq(n1, n2), digits = 3))) %>%
  
  ggplot() +
  geom_point(aes(x = log(n1), y = log(n2), color = st), stroke = 2, size = 0.5, shape = 4) +
  geom_abline(linetype = "dashed") +
  geom_text(aes(x = 1, y = 4.8, label = paste0("r = ", r[1])), size = 4, family = "American typewriter") +
  geom_text(aes(x = 1, y = 4.5, label = paste0("p = ", p[1])), size = 4, family = "American typewriter") +
  scale_x_continuous(limit = c(0, 5), breaks = seq(0, 6, 1)) +
  scale_y_continuous(limit = c(0, 5), breaks = seq(0, 6, 1)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "observed log_ip, 2015-2019", y = "predicted log_ipd, 2015-2019") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(legend.position = c(0.9,0.4), legend.title = element_blank()) +
  theme(legend.key.size = unit(0.4, "cm"))


#save combined plots
ggsave(here("output", "sfig12_mw_ipdPredictions.png"),
       plot = (A ), 
       width = 8, height = 6, unit = "in", dpi = 300)
