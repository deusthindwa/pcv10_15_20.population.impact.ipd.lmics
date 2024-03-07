#Authors: Deus & Dan
#Date: from 03/03/2024
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries
 
#====================================================================

#infer carriage in pre-pv13 era (carriage  <- ipd / invasiveness)
mw_carInf <-
  mw_ipd %>%
  dplyr::select(yearc, st) %>%
  group_by(yearc, st) %>%
  tally() %>%
  ungroup() %>%
  rename("ipd" = "n") %>%
  mutate(st = if_else(st == "12FAB_44_46", "NVT", st),
         ipd = if_else(str_length(st)>3 & str_length(st)<=7, ipd/2, 
                       if_else(str_length(st)>7, ipd/3, ipd))) %>% 
  dplyr::select(yearc, st, ipd) %>%
  
#split multiple serotypes with "/" into new rows
  tidyr::separate_rows(., st) %>%
  group_by(yearc, st) %>%
  summarise(ipd = sum(ipd)) %>% 
  ungroup() %>%
  
  #join with invasiveness
  dplyr::left_join(st_inv) %>%
  dplyr::mutate(inv = if_else(is.na(st), 0, if_else(is.na(inv), inv_nvt, inv))) %>%
  
  #infer prevalence
  mutate(prev = ipd/inv,
         scaley = if_else(yearc >=2005 & yearc <=2011, 0.8/sum(prev), 0.5/sum(prev)),
         prev = scaley*prev) %>% 
  dplyr::select(yearc, st, ipd, inv, prev, everything(), -scaley) 

#replace missing serotype group with NVT
mw_carInf[is.na(mw_carInf)] <- "NVT"

#====================================================================

#2015-2018 samples summary
mw_none = 1004 #no samples
mw_nvt = 1109 #other non-vaccine samples

#infer IPD data (ipd  <- carriage * invasiveness)
mw_ipdInf <-
  mw_car %>%
  dplyr::select(yearc, st) %>%
  group_by(yearc, st) %>%
  tally() %>%
  ungroup() %>%
  rename("prev" = "n") %>%
  mutate(prev = if_else(is.na(st), prev, 
                        if_else(str_length(st)>3 & str_length(st)<=7, prev/2, 
                                if_else(str_length(st)>7, prev/3, prev)))) %>%
  dplyr::select(yearc, st, prev) %>%
  
#split multiple serotypes with "/" into new rows
  tidyr::separate_rows(., st) %>%
  group_by(yearc, st) %>%
  mutate(prev = sum(prev)) %>% 
  ungroup() %>%
  group_by(yearc) %>%
  mutate(prev = prev/sum(prev)) %>% 
  ungroup() %>%
  dplyr::filter(!is.na(st)) %>%

#join with invasiveness
  dplyr::left_join(st_inv) %>%
  dplyr::mutate(inv = if_else(is.na(inv), inv_nvt, inv)) %>%
  dplyr::mutate(ipd = prev*inv,
                scaley = 89/sum(ipd),
                ipd = scaley*ipd) %>%
  dplyr::select(yearc, st, prev, inv, ipd, everything(), -scaley)  
  
#replace missing serotype group with NVT
mw_ipdInf[is.na(mw_ipdInf)] <- "NVT"

#====================================================================


#generate ppotential pre-PCV carriage data
#postpcvIPD = prepcvIPD * postpcvCarr / prepcvCarr
#prepcvCarr = prepcvIPD * postpcvCarr / postpcvIPD
# 
# bind_cols(
#   mw_ipdb2011_obs %>%
#     mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT")) %>%
#     group_by(pcv13pfz) %>%
#     tally(ipd) %>%
#     ungroup() %>%
#     rename("ipd2011" = "n"),
#   
#   mw_ipda2015_obs %>%
#     mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT")) %>%
#     group_by(pcv13pfz) %>%
#     tally(ipd) %>%
#     ungroup() %>%
#     rename("ipd2015" = "n") %>%
#     dplyr::select(ipd2015),
#   
#   mw_cara2015_obs %>%
#     mutate(pcv13pfz = if_else(st == "None", "None", 
#                               if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"))) %>%
#     group_by(pcv13pfz) %>%
#     tally(prev) %>%
#     ungroup() %>%
#     rename("car2015" = "n") %>%
#     dplyr::filter(pcv13pfz != "None") %>%
#     dplyr::select(car2015)) %>%
#   
#   mutate(car2011 = ipd2011 * car2015 / ipd2015) %>%
#   
#   mutate(scalex = 0.70/sum(car2011), #assume total prevalence is up to 45%, Ellen et al.
#          car2011 = scalex*car2011) #scale prev according

# A <-
#   left_join(
#     mw_ipda2015_obs %>%
#       mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
#              st = if_else(pcv13pfz == "NVT", "NVT", st)) %>%
#       group_by(st) %>%
#       tally(ipd) %>% 
#       rename("n1"="n"),
#     
#     mw_ipda2015_pred %>%
#       mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
#              st = if_else(pcv13pfz == "NVT", "NVT", st)) %>%
#       group_by(st) %>%
#       tally(ipd) %>% 
#       rename("n2"="n")) %>%
#   mutate(r = abs(round((stats::cor(n1, n2, method = c("pearson")))[1], digits = 3)),
#          p = scientific(stats::cor.test(n1, n2)[3]$p.value),
#          r2 = abs(round(rsq(n1, n2), digits = 3))) %>%
#   
#   ggplot() +
#   geom_point(aes(x = log(n1), y = log(n2), color = st), stroke = 2, size = 0.5, shape = 4) +
#   geom_abline(linetype = "dashed") +
#   geom_text(aes(x = 1, y = 4.8, label = paste0("r = ", r[1])), size = 4, family = "American typewriter") +
#   geom_text(aes(x = 1, y = 4.5, label = paste0("p = ", p[1])), size = 4, family = "American typewriter") +
#   scale_x_continuous(limit = c(0, 5), breaks = seq(0, 6, 1)) +
#   scale_y_continuous(limit = c(0, 5), breaks = seq(0, 6, 1)) +
#   theme_bw(base_size = 16, base_family = "American typewriter") +
#   labs(title = "", x = "observed PCV13 VT log_ipd, 2015-2019", y = "predicted PCV13 VT log_ipd, 2015-2019") +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 2)) +
#   theme(legend.position = c(0.9,0.4), legend.title = element_blank()) +
#   theme(legend.key.size = unit(0.4, "cm"))
# 
# 
# #save combined plots
# ggsave(here("output", "sfig6_ipd_predValidation.png"),
#        plot = (A ), 
#        width = 8, height = 6, unit = "in", dpi = 300)
# 
