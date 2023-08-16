#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#data description of carriage
data_carr <-
data_carr %>%
  dplyr::mutate(country = word(study, 1, sep = "\\."),
                phase = str_detect(study, "pre"),
                phase = if_else(str_detect(study, "pre") == TRUE, "prePCV", "postPCV")) %>%
  dplyr::select(country, phase, period, everything(), -study)

#carriage descriptive plot
data_carr1 <-
data_carr %>% 
  dplyr::filter(country == "Bogota")

data_carr1 %>%
  dplyr::group_by(type, phase) %>%
  dplyr::mutate(prev = mean(prevcarr)) %>%
  dplyr::ungroup() %>%
  ggplot() +
  geom_point(aes(y = type, x = prev)) +
  theme_bw() +
  facet_grid(.~phase)

#====================================================================

#data description of ipd
data_ipd <-
  data_ipd %>%
  dplyr::mutate(country = word(study, 1, sep = "\\."),
                phase = str_detect(study, "pre"),
                phase = if_else(str_detect(study, "pre") == TRUE, "prePCV", "postPCV")) %>%
  dplyr::select(country, phase, period, everything(), -study)

#carriage descriptive plot
data_ipd1 <-
  data_ipd %>% 
  dplyr::filter(country == "Bogota")

data_ipd1 %>%
  dplyr::group_by(type, phase) %>%
  ggplot() +
  geom_point(aes(y = type, x = log(nipd))) +
  theme_bw() +
  facet_grid(.~phase)
