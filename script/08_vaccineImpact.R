#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
#====================================================================

#plot of serotype replacement scenarios (post-VT eliminated)
A <-
  bind_rows(is_err, sa_err) %>%
  ggplot() +
  geom_line(aes(x = sr, y = err_diff, color = country), lty = "twodash", size = 1) +
  #geom_point(aes(x = sr[which.min(err_diff)], y = min(err_diff), color = country), size = 2, stroke = 1, shape = 8) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "proportion of NVT replacing VT carriage", y = "predicted vs observed median IRR error") + 
  scale_x_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.25)) + 
  theme(strip.text.x = element_text(size = 0), strip.background = element_rect(fill = "gray90")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#plot of model validation (obs vs pred distributions)
B <-
  bind_rows(is_validate, sa_validate) %>%
  ggplot() +
  geom_density_ridges(aes(x = log(irr), y = sr, fill = sr), size = 0.5, alpha = 0.5) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "log_incidence rate ratio", y = "Density") +
  facet_grid(.~country, scales = "free_y") +
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 0)) +
  guides(fill = guide_legend(title = "Model validation and scenarios\nof serotype replacement (SR)")) +
  theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig6_obsvspred.png"),
       plot = (A | B),
       width = 24, height = 9, unit = "in", dpi = 300)

#====================================================================

#plot of vaccine preventable IPD distributions
C <-
  bind_rows(is_pcvsamples, sa_pcvsamples) %>%
  ggplot() +
  geom_density(aes(x = 1-irr1, group = pcv, fill = "no SR"), size = 0.6, alpha = 0.3) +
  #geom_density(aes(x = 1-irr2, group = pcv, fill = "baseline SR"), size = 0.6, alpha = 0.3) +
  geom_density(aes(x = 1-irr3, group = pcv, fill = "complete SR"), size = 0.6, alpha = 0.3) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(country~pcv, scales = "free_y") +
  scale_fill_manual(name = "scenarios of serotype\nreplacement (SR)", values = c("no SR" = "blue", "baseline SR" = "red", "complete SR" = "green")) +
  labs(title = "", x = "vaccine impact (proportion of preventable IPD)", y = "density") + 
  theme(legend.text = element_text(size = 14), legend.position = "right", legend.title = element_text(size = 14), legend.key.size = unit(1.2,"cm")) +
  theme(strip.text.x = element_text(size = 26), strip.text.y = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig8_vaximpact.png"),
       plot = (C), 
       width = 24, height = 9, unit = "in", dpi = 300)

#====================================================================


