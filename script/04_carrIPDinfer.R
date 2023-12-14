#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries
 
#====================================================================
#INVASIVENESS DATA FROM SOUTH AFRICA 
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
inv_nvt = (inv_pcv20pfz$wgt.inv[2] + inv_pcv15mek$wgt.inv[2] + inv_pcv13pfz$wgt.inv[2] + inv_pcv10sii$wgt.inv[2] + inv_pcv10gsk$wgt.inv[2] + inv_pcv7pfz$wgt.inv[2])/6 #average NVT invasiveness
inv_vt = (inv_pcv20pfz$wgt.inv[1] + inv_pcv15mek$wgt.inv[1] + inv_pcv13pfz$wgt.inv[1] + inv_pcv10sii$wgt.inv[1] + inv_pcv10gsk$wgt.inv[1] + inv_pcv7pfz$wgt.inv[1])/6 #average NVT invasiveness

#create ipd serotype dataset to match those of invasiveness
#infer carriage data pre-pcv13 introduction in south africa (carriage  <- ipd / invasiveness)
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


#plots of the relationship between serotype group-specific carriage prevalence and ipd and ivasiveness
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
  labs(title = "", x = "PCV regimen", y = "serotype group carriage prevalence") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
  theme(plot.title = element_text(hjust = 0.1, vjust = -10), axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  facet_grid(.~ factor(period, levels = c("pre-PCV", "post-PCV")))

#save combined plots
ggsave(here("output", "sfig10_carsgDist.png"),
       plot = (A), 
       width = 20, height = 7, unit = "in", dpi = 300)

#====================================================================
#INVASIVENESS DATA FROM MALAWI
#====================================================================



