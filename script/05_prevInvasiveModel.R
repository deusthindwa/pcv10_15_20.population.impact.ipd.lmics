#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#model
#expected_IPD = constant + log_carriage_prevalence + log_est_invasiveness
#use this model to predict expected IPD for each country's carriage data (?using estimated par for one country)

#combine carriage, ipd and invasiveness datasets
data_Bogota <- 
  data_all %>% dplyr::filter(country == "Bogota") %>%
  left_join(data_inv)

#fit a negative-binomial model of carriage prevalence and estimated invasiveness
model1 <- tidy(MASS::glm.nb(log_nipd ~ log_prevcarr + log_inv, data = data_Bogota), 
                   exponentiate = TRUE, 
                   conf.int = TRUE, 
                   conf.level = 0.95)

bogota <- bogota %>% 
  dplyr::mutate(fit_ipd = model_carr$estimate[1] + model_carr$estimate[2]*log_ncarr + model_carr$estimate[3]*log_inv)

#visualize relationship between predicted and observed IPD
bogota %>%
  ggplot() +
  geom_point(aes(x = fit_ipd, y = log_nipd, color = st)) +
  theme_bw() + 
  theme(legend.position = "right")
