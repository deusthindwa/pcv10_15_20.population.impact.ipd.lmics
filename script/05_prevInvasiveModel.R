#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#model
#expected_IPD = constant + log_carriage_prevalence + log_est_invasiveness
#use this model to predict expected IPD for each country's carriage data (?using estimated par for one country)

#combine carriage, ipd and invasiveness datasets
# data_carr_ipd <- rio::import("https://raw.githubusercontent.com/nickjcroucher/progressionEstimation/main/data-raw/S_pneumoniae_infant_serotype.csv")
# data_carr_ipd %>% readr::write_csv(x = ., file = here("output", "data_carr_ipd.csv"))
data_carr_ipd <- 
  rio::import(here("output", "data_carr_ipd.csv")) %>%
  dplyr::mutate(country = word(study, 1, sep = "\\."),
                phase = str_detect(study, "pre"),
                phase = if_else(str_detect(study, "pre") == TRUE, "pre-pcv carriage", "post-pcv carriage")) %>%
  dplyr::select(country, phase, time_interval, type, carriage_samples, carriage, surveillance_population, disease) %>%
  dplyr::rename("period" = "time_interval", "st" = "type",  "nsamples"= "carriage_samples", "ncarr" = "carriage",  "npop" = "surveillance_population", "nipd" = "disease") %>%
  dplyr::mutate(log_ncarr = log(ncarr+1),
                log_nipd = log(nipd+1))

bogota <-
  data_inv %>%
  left_join(data_carr_ipd) %>%
  dplyr::filter(country == "Bogota", phase == "pre-pcv carriage")

#fit a negative-binomial model of carriage prevalence and estimated invasiveness
model_carr <- tidy(MASS::glm.nb(log_nipd ~ log_ncarr + log_inv, data = bogota), 
                   exponentiate = TRUE, 
                   conf.int = TRUE, 
                   conf.level = 0.95)

bogota <- bogota %>% 
  dplyr::mutate(fit_ipd = model_carr$estimate[1] + model_carr$estimate[2]*log_ncarr + model_carr$estimate[3]*log_inv)

#visualize relationship between predicted and observed IPD
bogota %>%
  ggplot() +
  geom_point(aes(x = fit_ipd, y = log_nipd)) +
  theme_bw() + 
  theme(legend.position = "right")



bogota1 %>%
  ggplot() +
  geom_point(aes(x = prop_fit_ipd, y = prop_nipd)) +
  theme_bw() + 
  theme(legend.position = "right")
