#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
#====================================================================

#plot of model validation (obs vs pred)
predictions <-
  bind_rows(is_validate, sa_validate, mw_validate) %>%
  mutate(sr = if_else(sr == "observed IRR", "obs", 
                      if_else(sr == "predicted IRR, baseline SR", "predba", 
                              if_else(sr == "predicted IRR, complete SR", "predco", "predno")))) %>%
  group_by(country, sr) %>%
  summarise(irr_mean = quantile(irr, 0.500),
            irr_low = quantile(irr, 0.025),
            irr_high = quantile(irr, 0.975)) %>%
  pivot_wider(names_from = sr, values_from = c(irr_mean, irr_low, irr_high)) %>%
  mutate(srp = if_else(country == "Israel", is_err$sr_min[1],
                       if_else(country == "South Africa", sa_err$sr_min[1], 
                               if_else(country == "Malawi", mw_err$sr_min[1], NA_real_))))

A <-
  predictions %>%
  ggplot(aes(x = irr_mean_obs, y = irr_mean_predno, color = country)) +
  geom_point(shape = 4, stroke = 2, size = 3) + 
  geom_errorbar(aes(xmin = irr_low_obs, xmax = irr_high_obs), width = 0, size = 1, alpha = 0.5) +
  geom_errorbar(aes(ymin = irr_low_predno, ymax = irr_high_predno), width = 0, size = 1, alpha = 0.4) +
  geom_abline(linetype = "dashed") +
  #coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(limit = c(0, 1.5), breaks = seq(0, 1.5, 0.4)) + 
  scale_y_continuous(limit = c(0, 1.5), breaks = seq(0, 1.5, 0.4)) + 
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "observed IRR", y = "Predicted IRR (SR = 0)") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), legend.position = "none")

B <-
  predictions %>%
  ggplot(aes(x = irr_mean_obs, y = irr_mean_predco, color = country)) +
  geom_point(shape = 4, stroke = 2, size = 3) + 
  geom_errorbar(aes(xmin = irr_low_obs, xmax = irr_high_obs), width = 0, size = 1, alpha = 0.5) +
  geom_errorbar(aes(ymin = irr_low_predco, ymax = irr_high_predco), width = 0, size = 1, alpha = 0.5) +
  geom_abline(linetype = "dashed") +
  #coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(limit = c(0, 1.5), breaks = seq(0, 1.5, 0.4)) + 
  scale_y_continuous(limit = c(0, 1.5), breaks = seq(0, 1.5, 0.4)) + 
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "observed IRR", y = "Predicted IRR (SR = 1)") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), legend.position = "none")

C <-
  predictions %>%
  ggplot(aes(x = irr_mean_obs, y = irr_mean_predba, color = country)) +
  geom_point(shape = 4, stroke = 2, size = 3) + 
  geom_errorbar(aes(xmin = irr_low_obs, xmax = irr_high_obs), width = 0, size = 1, alpha = 0.5) +
  geom_errorbar(aes(ymin = irr_low_predba, ymax = irr_high_predba), width = 0, size = 1, alpha = 0.5) +
  geom_text(aes(x = 0.8, y = 0.25, label = paste0("SR = ", srp[1])), color = "#F8766D", size = 5, fontface = "bold", family = "American typewriter") +
  geom_text(aes(x = 0.8, y = 0.2, label = paste0("SR = ", srp[2])), color = "#7CAE00", size = 5, fontface = "bold", family = "American typewriter") +
  geom_text(aes(x = 0.8, y = 0.15, label = paste0("SR = ", srp[3])), color = "#619CFF", size = 5, fontface = "bold", family = "American typewriter") +
  geom_text(aes(x = 0.8, y = 0.1, label = paste0("SR = ", 0)), color = "#C77CFF", size = 5, fontface = "bold", family = "American typewriter") +
  geom_abline(linetype = "dashed") +
  #coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(limit = c(0, 1.5), breaks = seq(0, 1.5, 0.4)) + 
  scale_y_continuous(limit = c(0, 1.5), breaks = seq(0, 1.5, 0.4)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "observed IRR", y = "Predicted IRR") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2), legend.position = "right")

#save combined plots
ggsave(here("output", "sfig11_obspred.png"),
       plot = (A | B |C), 
       width = 20, height = 6, unit = "in", dpi = 300)
  
#====================================================================

#plot of vaccine preventable IPD distributions
C <-
  bind_rows(is_pcvsamples, sa_pcvsamples, mw_pcvsamples) %>%
  ggplot() +
  geom_density(aes(x = 1-irr1, group = pcv, fill = "no SR"), size = 0.6, alpha = 0.3) +
  geom_density(aes(x = 1-irr2, group = pcv, fill = "baseline SR"), size = 0.6, alpha = 0.3) +
  geom_density(aes(x = 1-irr3, group = pcv, fill = "complete SR"), size = 0.6, alpha = 0.3) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(country ~ pcv, scales = "free_x") +
  scale_fill_manual(name = "scenarios of serotype\nreplacement (SR)", values = c("no SR" = "blue", "baseline SR" = "red", "complete SR" = "green")) +
  labs(title = "", x = "vaccine impact (proportion of preventable IPD)", y = "density") + 
  theme(legend.text = element_text(size = 14), legend.position = "right", legend.title = element_text(size = 14), legend.key.size = unit(1.2,"cm")) +
  theme(strip.text.x = element_text(size = 26), strip.text.y = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig8_vaximpact.png"),
       plot = (C), 
       width = 24, height = 12, unit = "in", dpi = 300)

#====================================================================



