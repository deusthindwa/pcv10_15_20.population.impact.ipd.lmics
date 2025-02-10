#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================

#Impact dataset
impactDS <-
  dplyr::bind_rows(sa_val, is_val, br_val, mw_val) %>%
  dplyr::filter(period == 2) %>%
  dplyr::rename("prev0"="carr_est", "ipd0"="N_cases") %>%
  dplyr::mutate(newVT_10g = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_10 = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_13 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, 1L, 0L),
                newVT_15 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F|22F|33F)\\b", st) == TRUE, 1L, 0L),
                newVT_20 = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F|22F|33F|8|10A|11A|12F|15B)\\b", st) == TRUE, 1L, 0L))

#check VTs versus NVTs levels
impactDS %>% group_by(country, newVT_10g) %>% tally(prev0)
impactDS %>% group_by(country, newVT_10) %>% tally(prev0)
impactDS %>% group_by(country, newVT_13) %>% tally(prev0)
impactDS %>% group_by(country, newVT_15) %>% tally(prev0)
impactDS %>% group_by(country, newVT_20) %>% tally(prev0)

#match vaccine effectiveness against each serotype disease and carriage dataset
impactDS <- 
  impactDS %>%
  left_join(st_VEd %>% dplyr::select(st, VEd_pcv10g, VEd_pcv10s, VEd_pcv13, VEd_pcv15, VEd_pcv20)) %>%
  replace(is.na(.), 0) %>%
  left_join(st_VEc %>% dplyr::select(st, VEc_pcv10g, VEc_pcv10s, VEc_pcv13, VEc_pcv15, VEc_pcv20)) %>%
  replace(is.na(.), 0)

#configure proportion of NVT that needs to fill the hole in the niche (gamma)
#set the relative change in the proportion of new vaccine coverage (vax_uptake)
gamma = 1
vax_uptake = 0.85

#determined whether there will be serotype expansion or not
#set 'a' to be true if NVT and false if VT
impactDS <- 
  impactDS %>% 
  dplyr::group_by(country) %>%
  dplyr::mutate(a10g = (((1-VEc_pcv10g) >1)|((1-VEc_pcv10g) ==1 & newVT_10g !=1)),
                a10 = (((1-VEc_pcv10s) >1)|((1-VEc_pcv10s) ==1 & newVT_10 !=1)),
                a13 = (((1-VEc_pcv13) >1)|((1-VEc_pcv13) ==1 & newVT_13 !=1)),
                a15 = (((1-VEc_pcv15) >1)|((1-VEc_pcv15) ==1 & newVT_15 !=1)),
                a20 = (((1-VEc_pcv20) >1)|((1-VEc_pcv20) ==1 & newVT_20 !=1))) %>%
  dplyr::ungroup()

#computed change in colonization post-PCV switch
impactDS <-
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(col10g = dplyr::if_else(a10g == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10g*(1-a10g)))/sum(a10g*prev0))), sum(prev0*VEc_pcv10g*(1-a10g))),
                col10 = dplyr::if_else(a10 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv10s*(1-a10)))/sum(a10*prev0))), sum(prev0*VEc_pcv10s*(1-a10))),
                col13 = dplyr::if_else(a13 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv13*(1-a13)))/sum(a13*prev0))), sum(prev0*VEc_pcv13*(1-a13))),
                col15 = dplyr::if_else(a15 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv15*(1-a15)))/sum(a15*prev0))), sum(prev0*VEc_pcv15*(1-a15))),
                col20 = dplyr::if_else(a20 == TRUE, (1 + ((gamma*sum(prev0*VEc_pcv20*(1-a20)))/sum(a20*prev0))), sum(prev0*VEc_pcv20*(1-a20)))) %>%
  dplyr::ungroup()

#compute the expected IPD incidence of each serotype in the new PCV era after switch
impactDS <- 
  impactDS %>%
  dplyr::group_by(country) %>%
  dplyr::mutate(ipd10g = ipd0 * (1-VEd_pcv10g * vax_uptake) * (col10g * vax_uptake),
                ipd10 = ipd0 * (1-VEd_pcv10s * vax_uptake) * (col10 * vax_uptake),
                ipd13 = ipd0 * (1-VEd_pcv13 * vax_uptake) * (col13 * vax_uptake),
                ipd15 = ipd0 * (1-VEd_pcv15 * vax_uptake) * (col15 * vax_uptake),
                ipd20 = ipd0 * (1-VEd_pcv20 * vax_uptake) * (col20 * vax_uptake)) %>%
  dplyr::ungroup()

#====================================================================

#check total ipd (all serotypes) after PCV switch
bs_samples = 10000 #number of bootstrap samples

#ISRAEL
is_ipdf <- impactDS %>% dplyr::filter(country == "Israel") 
is_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd10))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd15))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd20))/sum(is_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(is_ipdf$ipd0), p = (sum(is_ipdf$ipd0)-sum(is_ipdf$ipd13))/sum(is_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Israel")

is_impact <-
  bind_rows(
    is_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    is_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    is_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#SOUTH AFRICA
sa_ipdf <- impactDS %>% dplyr::filter(country == "South Africa") 
sa_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd10))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd15))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd20))/sum(sa_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(sa_ipdf$ipd0), p = (sum(sa_ipdf$ipd0)-sum(sa_ipdf$ipd13))/sum(sa_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "South Africa")

sa_impact <-
  bind_rows(
    sa_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    sa_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    sa_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#MALAWI
mw_ipdf <- impactDS %>% dplyr::filter(country == "Malawi") 
mw_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd10))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd15))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd20))/sum(mw_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(mw_ipdf$ipd0), p = (sum(mw_ipdf$ipd0)-sum(mw_ipdf$ipd13))/sum(mw_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Malawi")

mw_impact <-
  bind_rows(
    mw_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    mw_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    mw_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#BRAZIL
br_ipdf <- impactDS %>% dplyr::filter(country == "Brazil") 
br_impact <-
  dplyr::bind_cols(
    pcv10 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv15 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd15))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1))),
    pcv20 = replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd20))/sum(br_ipdf$ipd0), size = 1))) - replicate(bs_samples, mean(rbinom(n = sum(br_ipdf$ipd0), p = (sum(br_ipdf$ipd0)-sum(br_ipdf$ipd10g))/sum(br_ipdf$ipd0), size = 1)))) %>%
  dplyr::mutate(country = "Brazil")

br_impact <-
  bind_rows(
    br_impact %>% dplyr::select(pcv10, country) %>% dplyr::rename("impact" = "pcv10") %>% dplyr::mutate(pcv = "PCV10S"),
    br_impact %>% dplyr::select(pcv15, country) %>% dplyr::rename("impact" = "pcv15") %>% dplyr::mutate(pcv = "PCV15"),
    br_impact %>% dplyr::select(pcv20, country) %>% dplyr::rename("impact" = "pcv20") %>% dplyr::mutate(pcv = "PCV20"))

#====================================================================

#plot of vaccine preventable IPD distributions
A <-
  bind_rows(is_impact, sa_impact, mw_impact, br_impact) %>%
  ggplot() +
  geom_density(aes(x = impact, group = pcv, fill = pcv), size = 1, alpha = 0.5) +
  theme_bw(base_size = 20, base_family = "American typewriter") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = paly) +
  scale_y_continuous(limit = c(0, 16), breaks = seq(0, 16, 3)) +
  scale_x_continuous(limit = c(-0.3, 0.6)) +
  labs(title = "", x = "vaccine net impact (additional preventable IPD)", y = "probability density") + 
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 22), strip.background = element_rect(fill = "gray90")) +
  geom_text(data = dplyr::filter(mw_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = -0.15, y = 14, label = 'favors\nexisting PCV', size = 6, family = "American typewriter") +
  geom_curve(data = dplyr::filter(mw_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = -0.05, y = 9, yend = 9, xend = -0.20, linewidth = 1.5, curvature = 0, arrow = arrow(length = unit(0.7, 'cm'))) +
  geom_text(data = dplyr::filter(mw_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = 0.15, y = 14, label = 'favors\nnew PCV', size = 6, family = "American typewriter") +
  geom_curve(data = dplyr::filter(mw_impact, pcv == "PCV15"), aes(x = impact, group = pcv), x = 0.05, y = 9, yend = 9, xend = 0.20, linewidth = 1.5, curvature = 0, arrow = arrow(length = unit(0.7, 'cm'))) +
  facet_grid(factor(country, levels = c("Israel", "South Africa", "Brazil", "Malawi")) ~ pcv, scales = "free_x") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 3)) +
  theme(legend.position = "none")

#save combined plots
ggsave(here("output", "PCV_net_benefit.png"),
       plot = (A), 
       width = 20, height = 13, unit = "in", dpi = 300)

#====================================================================

#summarizing uncertainty
B <-
  bind_rows(is_impact, sa_impact, mw_impact, br_impact) %>%
  dplyr::mutate(country = factor(country, levels = c("Israel", "South Africa", "Brazil", "Malawi"))) %>%
  group_by(country, pcv) %>%
  summarise(imp1M = round(quantile(impact, 0.500), digits = 3),
            imp1L = round(quantile(impact, 0.025), digits = 3),
            imp1U = round(quantile(impact, 0.975), digits = 3))

#store estimates
write_csv(B, here("output", "PCV_net_benefit.csv"))
