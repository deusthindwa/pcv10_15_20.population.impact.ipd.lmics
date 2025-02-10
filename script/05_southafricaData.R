#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#South Africa IPD distribution
A <-
  sa_ipd %>%
  dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
  dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, 1), country = "South Africa", yearc = factor(yearc)) %>%
  tidyr::separate_rows(., st) %>%
  group_by(country, yearc, st) %>%
  tally() %>%
  mutate(N = sum(n)) %>%
  ungroup() %>% 
  filter(!is.na(st)) %>%
  group_by(country, st) %>%
  mutate(NN = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = log(n+0.5), y = reorder(st, NN), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack") +
  theme_bw(base_size = 20, base_family = "American Typewriter") +
  scale_fill_manual(values = palx) +
  scale_x_continuous(limit = c(0, 53), breaks = seq(0, 53, 4)) +
  labs(title = "", x = "Reported log_number of IPD isolates in South Africa", y = "Pneumococcal serotype") + 
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 14), legend.position = c(0.7, 0.4), legend.title = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

ggsave(here("output", "Southafrica_ipd_distribution.png"), plot = A, width = 25, height = 15, unit = "in", dpi = 300)

#====================================================================

#South Africa carriage distribution
B <-
  sa_ipd %>%
  dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
  dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, 1), country = "South Africa", yearc = factor(yearc)) %>%
  tidyr::separate_rows(., st) %>%
  group_by(country, yearc, st) %>%
  dplyr::tally() %>%
  dplyr::ungroup() %>%
  dplyr::left_join(st_inv, by = 'st') %>%
  dplyr::mutate(inv = if_else(is.na(inv), inv_nvt, inv)) %>%
  dplyr::group_by(country, yearc) %>%
  dplyr::mutate(p = (n+1)/inv, p = p/sum(p, na.rm=T)*0.2) %>%  #inferred carriage
  dplyr::ungroup() %>%
  group_by(st) %>%
  mutate(NN = sum(n)) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = p, y = reorder(st, p), fill = fct_rev(factor(yearc))), size = 0.3, color = "black", position = "stack") +
  theme_bw(base_size = 20, base_family = "American Typewriter") +
  scale_fill_manual(values = palx) +
  labs(title = "", x = "Inferred serotype carriage prevalence in South Africa", y = "Pneumococcal serotype") + 
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16)) +
  scale_x_continuous(limit = c(0, 0.37), breaks = seq(0, 0.37, 0.03), labels = scales::percent_format(accuracy = 1)) + 
  guides(fill = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 14), legend.position = c(0.7, 0.4), legend.title = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 3))

#save combined plots
ggsave(here("output", "Southafrica_carriage_distribution.png"), plot = B, width = 25, height = 15, unit = "in", dpi = 300)

#====================================================================

#infer carriage data (carriage  <- ipd / invasiveness)
sa_cari <-
  sa_ipd %>%
  dplyr::mutate(st = if_else(st == "12B/F", "12B/12F", if_else(st == "15B/C", "15B/15C", st))) %>%
  dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, if_else(str_length(st) > 7, 1/3, 1))) %>%
  tidyr::separate_rows(., st) %>%
  group_by(yearc, st) %>%
  summarise(n = sum(n)) %>% 
  ungroup() %>%
  tidyr::complete(st, yearc, fill = list(n = 0)) %>%
  dplyr::left_join(st_inv, by = 'st') %>%
  dplyr::mutate(inv = if_else(is.na(inv), inv_nvt, inv)) %>% #replace mean invasiveness for all NVTs
  dplyr::group_by(yearc) %>%
  dplyr::mutate(p = (n+1)/inv, 
                p = p/sum(p, na.rm=T)*0.2,
                country = "South Africa") %>%  #inferred carriage
  dplyr::ungroup()
