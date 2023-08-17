#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#manipulation of carriage and ipd data
data_carr1 <-
data_carr %>%
  dplyr::mutate(country = word(study, 1, sep = "\\."),
                phase = str_detect(study, "pre"),
                phase = if_else(str_detect(study, "pre") == TRUE, "pre-pcv carriage", "post-pcv carriage"),
                ncarr = ncarr+1) %>% #add 1 positive carriage sample to avoid math errors
  dplyr::select(country, phase, period, everything(), -study)

data_ipd1 <-
  data_ipd %>%
  dplyr::mutate(country = word(study, 1, sep = "\\."),
                phase = str_detect(study, "pre"),
                phase = if_else(str_detect(study, "pre") == TRUE, "pre-pcv ipd", "post-pcv ipd"),
                nipd = nipd+1) %>% #add 1 positive ipd isolate to avoid math errors
  dplyr::select(country, phase, period, everything(), -study)

#====================================================================

#save descriptive plot for carriage and ipd
ggsave(here("output", "fig1_carripdDesc.png"),
plot = (

#carriage descriptive plot
data_carr1 %>% 
  dplyr::filter(country == "Bogota") %>%
  ggplot() +
  geom_point(aes(y = st, x = log(ncarr), color = phase), size = 2.5) +
  theme_bw() +
  facet_grid(.~factor(phase, levels = c("pre-pcv carriage", "post-pcv carriage"))) +
  theme(legend.position = "none") |

#ipd description plot
data_ipd1 %>% 
  dplyr::filter(country == "Bogota") %>%
  ggplot() +
  geom_point(aes(y = st, x = log(nipd), color = phase), size = 2.5) +
  theme_bw() +
  facet_grid(.~ factor(phase, levels = c("pre-pcv ipd", "post-pcv ipd"))) + 
  theme(legend.position = "none")

), width = 10, height = 8, unit = "in", dpi = 300)
  
