#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#description of the invasiveness, carriage and ipd datasets

#invasiveness descriptive plot
A <-
  data_inv %>% 
  ggplot() +
  geom_point(aes(x = reorder(st, log_inv), y = log_inv), size = 2.5, shape = 4, stroke = 2) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(A), invasiveness", x = "Log_invasiveness", y = "Pneumococcal serotype") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Malawi
B <-
  data_desc %>%
  dplyr::filter(country == "Bogota") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, est), x = est, color = cat), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(values=c("darkgreen", "darkred")) +
  scale_x_continuous(sec.axis = sec_axis(~ . - 20/10)) +
  labs(title = "(B)", x = "Log_carriage prevalence/IPD", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.7,0.2), legend.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Bogota")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Malawi
C <-
  data_desc %>%
  dplyr::filter(country == "Alabama") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, est), x = est, color = cat), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(values=c("darkgreen", "darkred")) +
  labs(title = "(C)", x = "Log_carriage prevalence/IPD", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.7,0.2), legend.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Alabama")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Malawi
D <-
  data_desc %>%
  dplyr::filter(country == "Caracas") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, est), x = est, color = cat), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(values=c("darkgreen", "darkred")) +
  labs(title = "(D)", x = "Log_carriage prevalence/IPD", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.7,0.2), legend.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Caracas")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Malawi
E <-
  data_desc %>%
  dplyr::filter(country == "Czech") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, est), x = est, color = cat), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(values=c("darkgreen", "darkred")) +
  labs(title = "(E)", x = "Log_carriage prevalence/IPD", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(legend.text = element_text(size = 12), legend.position = c(0.7,0.2), legend.title = element_text(size = 12)) +
  guides(color = guide_legend(title = "Czech")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

ggsave(here("output", "fig1_carripdDesc.png"),
       plot = (A / (B | C | D | E) + plot_layout(heights = c(1,2.5))), 
       width = 20, height = 10, unit = "in", dpi = 300)
