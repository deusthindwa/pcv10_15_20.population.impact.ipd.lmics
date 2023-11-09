#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#set the colors for each year
palx <- c("2006" = "#D16A5D", "2007" = "#E358A0", "2008" = "#8872D9", "2009" = "#C04AE1", "2010" = "#D598D9",
          "2011" = "#E0B1BC", "2012" = "#E3B45B", "2013" = "#D5E1E2", "2014" = "#84967C", "2015" = "#86AADA",
          "2016" = "#77DAD9", "2017" = "#84E4A6", "2018" = "#75EA5C", "2019" = "#C8DE5A")


#description of serotype-specific IPD isolates
#select limited vars for plotting from all datasets
X <-
bind_rows(
  mw_ipdb2011 %>% dplyr::select(yearc, st, country), 
  mw_ipda2015 %>% dplyr::select(yearc, st, country),
  is_ipdb2009 %>% dplyr::select(yearc, st, country),
  is_ipda2013 %>% dplyr::select(yearc, st, country),
  sa_ipdb2009 %>% dplyr::select(yearc, st, country),
  sa_ipda2015 %>% dplyr::select(yearc, st, country)) %>%
  dplyr::mutate(stx = if_else(grepl("4|6B|9V|14|18C|19F|23F|1|5|7F|3|6A|19A|22F|33F|15B|8|10A|11A|12F", st) == TRUE, st, "NVT")) %>%

bind_rows() %>%
  group_by(country, yearc, st) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N*100) %>%
  ungroup() %>%filter(!is.na(st)) %>%
  group_by(country, st) %>%
  mutate(NN = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = log(n+0.5), y = reorder(st, NN), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack", stat = "identity") +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_fill_manual(values = palx) +
  labs(title = "Malawi", x = "log_number of IPD isolates", y = "pneumococcal serotype") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  facet_wrap(.~country) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.4), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))


#save combined plots
ggsave(here("output", "fig2_serotypeRank.png"),
       plot = (A | B), 
       width = 12, height = 8, unit = "in", dpi = 300)


#Israel
bind_rows(is_ipdb2009, is_ipda2013) %>%
  group_by(yearc, st) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N*100) %>%
  ungroup() %>%filter(!is.na(st)) %>%
  group_by(st) %>%
  mutate(NN = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = log(n+0.5), y = reorder(st, NN), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack", stat = "identity") +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_fill_manual(values = palx) +
  labs(title = "Israel", x = "log_number of IPD isolates", y = "pneumococcal serotype") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.4), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#South Africa
X <-
bind_rows(sa_ipdb2009, sa_ipda2014) %>%
  group_by(yearc, st) %>%
  tally() %>%
  mutate(N = sum(n), p = n/N*100) %>%
  ungroup() %>%filter(!is.na(st)) %>%
  group_by(st) %>%
  mutate(NN = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_bar(aes(x = log(n+0.5), y = reorder(st, NN), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack", stat = "identity") +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_fill_manual(values = palx) +
  labs(title = "Israel", x = "log_number of IPD isolates", y = "pneumococcal serotype") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 12), legend.position = c(0.7, 0.4), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))



#set color for all plots
cols <- c("IPD"="darkviolet", "carriage" = "chartreuse4")

#====================================================================

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


