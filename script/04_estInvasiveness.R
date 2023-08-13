#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#model
#log_est_invasiveness = constant + log_case_fatality_ratio
#invasiveness are likely available for limited # of serotypes hence CFR
#invasiveness estimates likely unstable hence CFR

#fit a linear model of invasiveness and case fatality ratio
data_cfr <- 
  data_cfr %>%  
  mutate(log_inv = log(inv),
         log_cfr = log(cfr))

model_cfr <- tidy(lm(log_inv ~ log_cfr, data = data_cfr), 
                  exponentiate = TRUE, 
                  conf.int = TRUE, 
                  conf.level = 0.95)

data_cfr <- data_cfr %>% 
  dplyr::mutate(fit_inv = model_cfr$estimate[1] + model_cfr$estimate[2]*log_cfr)

#visualize relationship between case fatality rate and fitted IPD
data_cfr %>%
  ggplot() +
  geom_point(aes(x = cfr, y = inv, color = st)) +
  geom_line(aes(x = cfr, y = fit_inv)) +
  theme_bw() + 
  theme(legend.position = "none")
