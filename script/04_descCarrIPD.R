#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
#IPD SEROTYPE-SPECIFIC DISTRIBUTION PLOT
#====================================================================
#set the colors for each year
palx <- c("2005" = "#8B0019", "2006" = "#D16A5F", "2007" = "#E358A0", "2008" = "#8872D9", "2009" = "#C04AE1",
          "2010" = "#D598D9", "2011" = "#E0B1BC", "2012" = "#E3B45B", "2013" = "#D5E1E2", "2014" = "#84967C",
          "2015" = "#86AADA", "2016" = "#77DAD9", "2017" = "#84E4A6", "2018" = "#75EA5C", "2019" = "#C8DE5A")

#Malawi
A <-
  bind_rows(mw_ipdb2011, mw_ipda2015) %>%
  group_by(country, yearc, st) %>%
  tally() %>%
  mutate(N = sum(n)) %>%
  ungroup() %>% 
  filter(!is.na(st)) %>%
  group_by(country, st) %>%
  mutate(NN = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = log(n+0.5), y = reorder(st, NN), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack") +
  theme_bw(base_size = 16, base_family = "Lato") +
  scale_fill_manual(values = palx) +
  labs(title = "", x = "log_number of IPD isolates", y = "pneumococcal serotype") + 
  facet_wrap(.~country) +
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.4), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#Israel
B <-
  is_ipd %>%
  group_by(country, yearc, st) %>%
  tally() %>%
  mutate(N = sum(n)) %>%
  ungroup() %>% 
  filter(!is.na(st)) %>%
  group_by(country, st) %>%
  mutate(NN = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = log(n+0.5), y = reorder(st, NN), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack") +
  theme_bw(base_size = 16, base_family = "Lato") +
  scale_fill_manual(values = palx) +
  labs(title = "", x = "log_number of IPD isolates", y = "pneumococcal serotype") + 
  facet_wrap(.~country) +
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.4), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#South Africa
C <-
  #bind_rows(sa_ipdb2009, sa_ipda2015) %>%
  sa_ipd %>%
  group_by(country, yearc, st) %>%
  tally() %>%
  mutate(N = sum(n)) %>%
  ungroup() %>%filter(!is.na(st)) %>%
  group_by(country, st) %>%
  mutate(NN = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = log(n+0.5), y = reorder(st, NN), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack") +
  theme_bw(base_size = 16, base_family = "Lato") +
  scale_fill_manual(values = palx) +
  labs(title = "", x = "log_number of IPD isolates", y = "pneumococcal serotype") + 
  facet_wrap(.~country) +
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.4), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig1_ipdstDist.png"),
       plot = ((A/B) | C), 
       width = 18, height = 20, unit = "in", dpi = 300)

#====================================================================
#IPD VACCINE SEROTYPE-GROUP DISTRIBUTION PLOT
#====================================================================

paly = c("PCV7" = "#C04AE1", "PCV10-gsk" = "#D598D9", "PCV10-sii" = "#D16A5F", "PCV13" = "#E3B45B", "PCV15" = "#C8DE5A", "PCV20" = "#84967C")

#combine all IPD datasets and generate vaccine serotype
sg_ipd <-
bind_rows(
  mw_ipdb2011 %>% dplyr::select(yearc, st, country), 
  mw_ipda2015 %>% dplyr::select(yearc, st, country),
  is_ipd %>% dplyr::select(yearc, st, country),
  sa_ipd %>% dplyr::select(yearc, st, country)) %>%
  dplyr::mutate(pcv20pfz = if_else(grepl("1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F", st) == TRUE, "PCV20", "NVT"),
                pcv15mek = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F", st) == TRUE, "PCV15", "NVT"),
                pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", "NVT"),
                pcv10sii = if_else(grepl("1|5|6A|6B|7F|9V|14|19A|19F|23F", st) == TRUE, "PCV10-sii", "NVT"),
                pcv10gsk = if_else(grepl("1|4|5|6B|7F|9V|14|18C|19F|23F", st) == TRUE, "PCV10-gsk", "NVT"),
                pcv7pfz = if_else(grepl("4|6B|9V|14|18C|19F|23F", st) == TRUE, "PCV7", "NVT"))

D <-
bind_rows(
  sg_ipd %>%
  group_by(country, yearc, pcv20pfz) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N) %>%
  ungroup() %>% 
  rename("sg" = "pcv20pfz"),
  
  sg_ipd %>%
  group_by(country, yearc, pcv15mek) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N) %>%
  ungroup() %>% 
  rename("sg" = "pcv15mek"),
  
  sg_ipd %>%
  group_by(country, yearc, pcv13pfz) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N) %>%
  ungroup() %>% 
  rename("sg" = "pcv13pfz"),
  
  sg_ipd %>%
  group_by(country, yearc, pcv10sii) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N) %>%
  ungroup() %>% 
  rename("sg" = "pcv10sii"),
  
  sg_ipd %>%
  group_by(country, yearc, pcv10gsk) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N) %>%
  ungroup() %>% 
  rename("sg" = "pcv10gsk"),
  
  sg_ipd %>%
  group_by(country, yearc, pcv7pfz) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N) %>%
  ungroup() %>% 
  rename("sg" = "pcv7pfz")) %>%
  filter(!is.na(sg)) %>%
  
  dplyr::filter(sg != "NVT") %>%
  
  ggplot() +
  geom_point(aes(x = yearc, y = p, group = sg), size = 1.5, position = position_dodge(width = 1.5), stroke = 2, shape = 1) +
  geom_line(aes(x = yearc, y = p, color = fct_rev(fct_relevel(factor(sg), "PCV7", after = 0)), group = sg), size = 1.5, position = position_dodge(width = 1.5)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  scale_color_manual(values = paly) +
  scale_x_continuous(breaks = c(2005, 2009, 2013, 2017), limits = c(2005, 2019)) +
  labs(title = "", x = "sampling year", y = "proportion of vaccine serotype IPD") + 
  facet_wrap(.~country) +
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(color = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig2_ipdsgDist.png"),
       plot = (D), 
       width = 12, height = 9, unit = "in", dpi = 300)

#====================================================================
#CARRIAGE SEROTYPE-SPECIFIC DISTRIBUTION PLOT
#====================================================================

e <-
  mw_cara2015 %>%
  mutate(st = if_else(st == "None", NA_character_, st)) %>%
  group_by(country, yearc, st) %>%
  tally() %>%
  mutate(N = sum(n)) %>%
  ungroup() %>% 
  filter(!is.na(st), st == "NVT") %>%
  mutate(p = n/N) %>%
  ggplot() +
  geom_col(aes(x = p, y = reorder(st, p), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack") +
  theme_bw(base_size = 16, base_family = "Lato") +
  scale_fill_manual(values = palx) +
  labs(title = "", x = "Other NVT carriage prevalence", y = "") + 
  scale_x_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.2), labels = scales::percent_format(accuracy = 1)) + 
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

E <-
  mw_cara2015 %>%
  mutate(st = if_else(st == "None", NA_character_, st)) %>%
  group_by(country, yearc, st) %>%
  tally() %>%
  mutate(N = sum(n)) %>%
  ungroup() %>% 
  filter(!is.na(st), st != "NVT") %>%
  mutate(p = n/N) %>%
  ggplot() +
  geom_col(aes(x = p, y = reorder(st, p), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack") +
  theme_bw(base_size = 16, base_family = "Lato") +
  scale_fill_manual(values = palx) +
  labs(title = "", x = "serotype carriage prevalence", y = "pneumococcal serotype") + 
  facet_wrap(.~country) +
  scale_x_continuous(limit = c(0, 0.2), breaks = seq(0, 0.2, 0.04), labels = scales::percent_format(accuracy = 1)) + 
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.4), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

F <-
  is_car %>%
  group_by(country, yearc, st) %>%
  tally() %>%
  mutate(N = sum(n)) %>%
  ungroup() %>% 
  filter(!is.na(st)) %>%
  group_by(country, st) %>%
  mutate(NN = sum(n)) %>%
  mutate(p = n/N) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = p, y = reorder(st, p), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack") +
  theme_bw(base_size = 16, base_family = "Lato") +
  scale_fill_manual(values = palx) +
  labs(title = "", x = "serotype carriage prevalence", y = "pneumococcal serotype") + 
  facet_wrap(.~country) +
  scale_x_continuous(limit = c(0, 0.2), breaks = seq(0, 0.2, 0.04), labels = scales::percent_format(accuracy = 1)) + 
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.4), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig3_carstDist.png"),
       plot = (E | inset_element(e, right = 0.98, left = 0.3, bottom = 0.6, top = 0.8) | F), 
       width = 18, height = 13, unit = "in", dpi = 300)

#====================================================================
#CARRIAGE VACCINE SEROTYPE-GROUP DISTRIBUTION PLOT
#====================================================================

paly = c("PCV7" = "#C04AE1", "PCV10-gsk" = "#D598D9", "PCV10-sii" = "#D16A5F", "PCV13" = "#E3B45B", "PCV15" = "#C8DE5A", "PCV20" = "#84967C")

#combine all IPD datasts and generate vaccine serotype
sg_car <-
  bind_rows(
    mw_cara2015 %>% dplyr::select(yearc, st, country), 
    is_car %>% dplyr::select(yearc, st, country)) %>%
  mutate(st = if_else(st == "None", NA_character_, st)) %>%
  dplyr::mutate(pcv20pfz = if_else(grepl("1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F", st) == TRUE, "PCV20", 
                                   if_else(is.na(st), NA_character_, "NVT")),
                pcv15mek = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F", st) == TRUE, "PCV15", 
                                   if_else(is.na(st), NA_character_, "NVT")),
                pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", 
                                   if_else(is.na(st), NA_character_, "NVT")),
                pcv10sii = if_else(grepl("1|5|6A|6B|7F|9V|14|19A|19F|23F", st) == TRUE, "PCV10-sii", 
                                   if_else(is.na(st), NA_character_, "NVT")),
                pcv10gsk = if_else(grepl("1|4|5|6B|7F|9V|14|18C|19F|23F", st) == TRUE, "PCV10-gsk", 
                                   if_else(is.na(st), NA_character_, "NVT")),
                pcv7pfz = if_else(grepl("4|6B|9V|14|18C|19F|23F", st) == TRUE, "PCV7", 
                                  if_else(is.na(st), NA_character_, "NVT")))

G <-
  bind_rows(
    sg_car %>%
      group_by(country, yearc, pcv20pfz) %>%
      tally() %>%
      mutate(N = sum(n), p = n/N) %>%
      ungroup() %>% 
      rename("sg" = "pcv20pfz"),
    
    sg_car %>%
      group_by(country, yearc, pcv15mek) %>%
      tally() %>%
      mutate(N = sum(n), p = n/N) %>%
      ungroup() %>% 
      rename("sg" = "pcv15mek"),
    
    sg_car %>%
      group_by(country, yearc, pcv13pfz) %>%
      tally() %>%
      mutate(N = sum(n), p = n/N) %>%
      ungroup() %>% 
      rename("sg" = "pcv13pfz"),
    
    sg_car %>%
      group_by(country, yearc, pcv10sii) %>%
      tally() %>%
      mutate(N = sum(n), p = n/N) %>%
      ungroup() %>% 
      rename("sg" = "pcv10sii"),
    
    sg_car %>%
      group_by(country, yearc, pcv10gsk) %>%
      tally() %>%
      mutate(N = sum(n), p = n/N) %>%
      ungroup() %>% 
      rename("sg" = "pcv10gsk"),
    
    sg_car %>%
      group_by(country, yearc, pcv7pfz) %>%
      tally() %>%
      mutate(N = sum(n), p = n/N) %>%
      ungroup() %>% 
      rename("sg" = "pcv7pfz")) %>%
  filter(!is.na(sg)) %>%
  
  dplyr::filter(sg != "NVT") %>%
  
  ggplot() +
  geom_point(aes(x = factor(yearc), y = p, group = sg), size = 1.5, position = position_dodge(width = 1.5), stroke = 2, shape = 1) +
  geom_line(aes(x = factor(yearc), y = p, color = fct_rev(fct_relevel(factor(sg), "PCV7", after = 0)), group = sg), size = 1.5, position = position_dodge(width = 1.5)) +
  theme_bw(base_size = 16, base_family = "Lato") +
  scale_color_manual(values = paly) +
  labs(title = "", x = "sampling year", y = "vaccine serotype carriage prevalence") + 
  facet_wrap(.~country, scales = "free_x") +
  scale_y_continuous(limit = c(0, 0.65), breaks = seq(0, 0.65, 0.15), labels = scales::percent_format(accuracy = 1)) + 
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(color = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig4_carsgDist.png"),
       plot = (G), 
       width = 12, height = 9, unit = "in", dpi = 300)

#====================================================================
#IPD AGE DISTRIBUTION
#====================================================================

H <-
bind_rows(
  mw_ipdb2011 %>% mutate(era = "pre-PCV13"),
  mw_ipda2015 %>% mutate(era = "post-PCV13"),
  is_ipdb2009 %>% mutate(era = "pre-PCV13"),
  is_ipda2013 %>% mutate(era = "post-PCV13"),
  sa_ipdb2009 %>% mutate(era = "pre-PCV13"),
  sa_ipda2015 %>% mutate(era = "post-PCV13")) %>%
  
  ggplot(aes(x = agey, color = era, fill = era)) +
  geom_density(size = 1.5, alpha = 0.3) +
  theme_bw(base_size = 16, base_family = "Lato") +
  labs(title = "", x = "age (years)", y = "density of IPD isolates") + 
  facet_wrap(.~country, scales = "free_x") +
  scale_x_continuous(limit = c(1, 5), breaks = seq(1, 5, 1)) + 
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  guides(fill = guide_legend(title = "", label = FALSE), color = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig5_ipdAgeDist.png"),
       plot = (H), 
       width = 10, height = 7, unit = "in", dpi = 300)









#set color for all plots
cols <- c("IPD"="darkviolet", "carriage" = "chartreuse4")

#invasiveness descriptive plot
A <-
  data_inv %>% 
  ggplot() +
  geom_point(aes(x = reorder(st, log_inv), y = log_inv, color = log_inv), size = 2.5, shape = 5, stroke = 2) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(A)", x = "pneumococcal serotype", y = "log_invasiveness") + 
  scale_color_distiller(palette = "Reds", direction = 1) +
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Brazil
B <-
data_all %>%
  dplyr::filter(country == "Bogota") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1))) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Bogota", values = cols) + 
  facet_wrap(.~ factor(phase, levels = c("pre-pcv", "post-pcv")), scales = "free_x") +
  labs(title = "(B)", x = "invasive disease isolates", y = "pneumococcal serotype") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Kenya
C <-
  data_all %>%
  dplyr::filter(country == "Netherlands") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1))) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Netherlands", values = cols) + 
  labs(title = "(C)", x = "invasive disease isolates", y = "") + 
  facet_wrap(.~ factor(phase, levels = c("pre-pcv", "post-pcv")), scales = "free_x") +
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Malawi
D <-
  data_all %>%
  dplyr::filter(country == "Caracas") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1))) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Caracas", values = cols) + 
  facet_wrap(.~ factor(phase, levels = c("pre-pcv", "post-pcv")), scales = "free_x") +
  labs(title = "(D)", x = "invasive disease isolates", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for South Africa
E <-
  data_all %>%
  dplyr::filter(country == "Czech") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1))) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Czech", values = cols) + 
  facet_wrap(.~ factor(phase, levels = c("pre-pcv", "post-pcv")), scales = "free_x") +
  labs(title = "(E)", x = "invasive disease isolates", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#save combined plots
ggsave(here("output", "fig1_carripdDesc.png"),
       plot = (A / (B | C | D | E | plot_layout(widths = c(2,2,1,1))) + plot_layout(heights = c(1,2.3))), 
       width = 18, height = 12, unit = "in", dpi = 300)

#delete all plot objects
rm(list = grep("data_", ls(), value = TRUE, invert = TRUE))


