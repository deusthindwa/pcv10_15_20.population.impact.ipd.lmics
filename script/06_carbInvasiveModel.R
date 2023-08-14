#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#model
#Expected_IPD = constant + log_carbon/repeat_unit + log_est_invasiveness

#combine polysaccharide capsule, ipd and invasiveness datasets
data_modC <- 
  data_modP %>% 
  dplyr::left_join(data_caps) %>%
  dplyr::mutate(log_caps = log(caps),
                log_fit_inv = log(fit_inv)) %>%
  dplyr::select(st, fit_inv, log_fit_inv, caps, log_caps, ipd, log_ipd)

#fit a negative-binomial model of carriage prevalence and estimated invasiveness
model_caps <- tidy(MASS::glm.nb(log_ipd ~ log_caps + log_fit_inv, data = data_modC), 
                   exponentiate = TRUE, 
                   conf.int = TRUE, 
                   conf.level = 0.95)

data_modC <- data_modC %>% 
  dplyr::mutate(fit_ipd = model_caps$estimate[1] + model_caps$estimate[2]*log_caps + model_caps$estimate[3]*log_fit_inv)

#visualize relationship between polysaccharide capsule and fitted IPD
data_modC %>%
  ggplot() +
  geom_point(aes(x = caps, y = fit_ipd, color = st)) +
  theme_bw() + 
  theme(legend.position = "none")
