#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#description of the invasiveness, carriage and ipd datasets
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
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Brazil
B <-
data_all %>%
  dplyr::filter(country == "Bogota") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1)), breaks = c(0, 50, 100, 150)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Bogota", values = cols) + 
  labs(title = "(B)", x = "invasive disease isolates", y = "pneumococcal serotype") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.8,0.1), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Kenya
C <-
  data_all %>%
  dplyr::filter(country == "Alabama") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1)), breaks = c(0, 10, 20, 30)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Alabama", values = cols) + 
  labs(title = "(C)", x = "invasive disease isolates", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.8,0.1), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Malawi
D <-
  data_all %>%
  dplyr::filter(country == "Caracas") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1)), breaks = c(0, 10, 20, 30, 40, 50)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Caracas", values = cols) + 
  labs(title = "(D)", x = "invasive disease isolates", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.8,0.1), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for South Africa
E <-
  data_all %>%
  dplyr::filter(country == "Czech") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1)), breaks = c(0, 10, 20, 30, 40, 50)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Czech", values = cols) + 
  labs(title = "(E)", x = "invasive disease isolates", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.8,0.1), legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

ggsave(here("output", "fig1_carripdDesc.png"),
       plot = (A / (B | C | D | E) + plot_layout(heights = c(1,2.3))), 
       width = 18, height = 10, unit = "in", dpi = 300)
