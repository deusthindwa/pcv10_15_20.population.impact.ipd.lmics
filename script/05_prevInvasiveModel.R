#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#model
#expected_IPD = constant + log_carriage_prevalence + log_est_invasiveness
#use this model to predict expected IPD for each country's carriage data (?using estimated par for one country)

#combine carriage, ipd and invasiveness datasets
data_modP <- 
  data_cfr %>% 
  dplyr::left_join(data_carr) %>%
  dplyr::left_join(data_ipd) %>%
  dplyr::mutate(log_prev = log(prev),
                log_fit_inv = log(fit_inv),
                log_ipd = log(ipd)) %>%
  dplyr::select(st, fit_inv, log_fit_inv, prev, log_prev, ipd, log_ipd)

#fit a negative-binomial model of carriage prevalence and estimated invasiveness
model_carr <- tidy(MASS::glm.nb(log_ipd ~ log_prev + log_fit_inv, data = data_modP), 
                   exponentiate = TRUE, 
                   conf.int = TRUE, 
                   conf.level = 0.95)

data_modP <- data_modP %>% 
  dplyr::mutate(fit_ipd = model_carr$estimate[1] + model_carr$estimate[2]*log_prev + model_carr$estimate[3]*log_fit_inv)

#visualize relationship between carriage prevalence and fitted IPD
data_modP %>%
  ggplot() +
  geom_point(aes(x = prev, y = fit_ipd, color = st)) +
  theme_bw() + 
  theme(legend.position = "none")
