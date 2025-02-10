#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
# ipd and carriage data restructuring for plotting
#====================================================================

#description of carriage targeted by each PCV
desc_all <-
  
  bind_rows(
    
    #Israel IPD
    is_ipd %>% 
      dplyr::group_by(yearc, st) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
      dplyr::mutate(n = if_else(str_length(st) > 3, n/2, n)) %>%
      
      tidyr::separate_rows(., st) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally(n) %>%
      dplyr::mutate(p = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(source = "Proportion of IPD", country = "Israel")%>%
      dplyr::select(yearc,st,n,p,source,country),
    
    #Israel carriage
    is_car %>% 
      dplyr::group_by(yearc, st) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(st = if_else(is.na(st), "NA", st), st = if_else(st == "15BC", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
      dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
      tidyr::separate_rows(., st) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally(n) %>%
      dplyr::mutate(p = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(st !="NA") %>%
      dplyr::mutate(source = "Prevalence of carriage", country = "Israel")%>%
      dplyr::select(yearc,st,n,p,source,country),
    
    #South Africa IPD  
    sa_ipd %>%
      dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
      dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, 1)) %>%
      tidyr::separate_rows(., st) %>%
      group_by(yearc, st) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally(n) %>%
      dplyr::mutate(p = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(source = "Proportion of IPD", country = "South Africa")%>%
      dplyr::select(yearc,st,n,p,source,country),
    
    #South Africa carriage
    sa_cari %>%
      dplyr::select(yearc, st, n, p) %>%
      dplyr::mutate(p=p*3.5) %>%
      dplyr::mutate(source = "Prevalence of carriage", country = "South Africa") %>%
      dplyr::select(yearc,st,n,p,source,country),
    
    #Brazil IPD
    br_ipd %>%
      dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, 1)) %>%
      tidyr::separate_rows(., st) %>%
      group_by(country, yearc, st) %>%
      tally() %>%
      dplyr::mutate(p = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(source = "Proportion of IPD", country = "Brazil")%>%
      dplyr::select(yearc,st,n,p,source,country),
    
    #Brazil carriage
    br_car %>%
      dplyr::mutate(st = if_else(st == "15b", "15B", if_else(st == "15c", "15C", if_else(st == "16f", "16F", st)))) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(n = if_else(is.na(st), n, n/(stringr::str_count(st, fixed(' ')) + 1))) %>%
      tidyr::separate_rows(., st) %>%
      group_by(yearc, st) %>%
      tally(n) %>%
      dplyr::mutate(p = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(st)) %>%
      dplyr::mutate(source = "Prevalence of carriage", country = "Brazil") %>%
      dplyr::select(yearc,st,n,p,source,country),
    
    #Malawi IPD (based on actual reported IPD from monocle)
    mw_ipd %>%
      dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
      dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, 1)) %>%
      tidyr::separate_rows(., st) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally() %>%
      dplyr::mutate(p = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(source = "Proportion of IPD", country = "Malawi")%>%
      dplyr::select(yearc,st,n,p,source,country),
      #dplyr::filter(yearc !=2015),
    
    #Malawi IPD (based on inferred IPD from PCVPA surveys)
    mw_ipdi %>%
      dplyr::select(yearc, st, n, p) %>%
      dplyr::mutate(p=p*5) %>%
      dplyr::mutate(source = "Proportion of IPD", country = "Malawi") %>%
      dplyr::filter(!is.na(st)) %>%
      dplyr::select(yearc,st,n,p,source,country) %>%
      dplyr::filter(yearc !=2015),
    
    #Malawi carriage
    mw_car %>%
      dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", if_else(st == "12B/F", "12B/12F", st))) %>%
      dplyr::mutate(n = if_else(str_length(st) > 3, 1/2, if_else(str_length(st) > 7, 1/3, 1))) %>%
      tidyr::separate_rows(., st) %>%
      group_by(yearc, st) %>%
      tally() %>%
      dplyr::mutate(p = n/sum(n), yearc=as.integer(yearc)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(st)) %>%
      dplyr::mutate(source = "Prevalence of carriage", country = "Malawi")%>%
      dplyr::select(yearc,st,n,p,source,country)
    
  ) %>%

#assign vaccine serotype groups
  dplyr::mutate(pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "PCV20_NVT"),
                pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "PCV15_NVT"),
                pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "PCV13_NVT"),
                pcv10gsk = if_else(grepl("\\b(1|4|5|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "PCV10G", "PCV10G_NVT"), #6A excluded
                pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10S", "PCV10S_NVT"))

#====================================================================

desc_all <- 
  bind_rows(
    desc_all %>%
      dplyr::select(everything(), -pcv20pfz, -pcv15mek, -pcv13pfz, -pcv10sii) %>%
      dplyr::group_by(country, yearc, source, pcv10gsk) %>%
      dplyr::summarise(n = sum(n), p = sum(p)) %>%
      dplyr::ungroup() %>%
      dplyr::rename("pcvg" = "pcv10gsk"),
    
    desc_all %>%
      dplyr::select(everything(), -pcv20pfz, -pcv15mek, -pcv13pfz, -pcv10gsk) %>%
      dplyr::group_by(country, yearc, source, pcv10sii) %>%
      dplyr::summarise(n = sum(n), p = sum(p)) %>%
      dplyr::ungroup() %>%
      dplyr::rename("pcvg" = "pcv10sii"),
    
    desc_all %>%
      dplyr::select(everything(), -pcv20pfz, -pcv15mek, -pcv10sii, -pcv10gsk) %>%
      dplyr::group_by(country, yearc, source, pcv13pfz) %>%
      dplyr::summarise(n = sum(n), p = sum(p)) %>%
      dplyr::ungroup() %>%
      dplyr::rename("pcvg" = "pcv13pfz"),
    
    desc_all %>%
      dplyr::select(everything(), -pcv20pfz, -pcv13pfz, -pcv10sii, -pcv10gsk) %>%
      dplyr::group_by(country, yearc, source, pcv15mek) %>%
      dplyr::summarise(n = sum(n), p = sum(p)) %>%
      dplyr::ungroup() %>%
      dplyr::rename("pcvg" = "pcv15mek"),
    
    desc_all %>%
      dplyr::select(everything(), -pcv15mek, -pcv13pfz, -pcv10sii, -pcv10gsk) %>%
      dplyr::group_by(country, yearc, source, pcv20pfz) %>%
      dplyr::summarise(n = sum(n), p = sum(p)) %>%
      dplyr::ungroup() %>%
      dplyr::rename("pcvg" = "pcv20pfz"))

desc_vt <-
  desc_all %>%
  dplyr::filter(pcvg == "PCV10G" | pcvg == "PCV10S" | pcvg == "PCV13" | pcvg == "PCV15" | pcvg == "PCV20") %>%
  dplyr::mutate(country = factor(country, levels = c("Israel", "South Africa", "Brazil", "Malawi"))) %>%
  dplyr::filter(!((country=="Israel" | country=="South Africa" | country=="Malawi") & pcvg=="PCV10G"))

desc_nvt <-
  desc_all %>%
  dplyr::filter(pcvg == "PCV10G_NVT" | pcvg == "PCV10S_NVT" | pcvg == "PCV13_NVT" | pcvg == "PCV15_NVT" | pcvg == "PCV20_NVT") %>%
  dplyr::mutate(country = factor(country, levels = c("Israel", "South Africa", "Brazil", "Malawi"))) %>%
  dplyr::filter(!((country=="Israel" | country=="South Africa" | country=="Malawi") & pcvg=="PCV10G_NVT"))

#====================================================================

#ggplot of IPD frequency and carriage prevalence of vaccine-serotypes
A <-
  desc_vt %>%
  ggplot() +
  geom_line(aes(x = yearc, y = p, color = pcvg, group = pcvg), size = 1, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = yearc, y = p, color = pcvg, size = n), stroke = 2.5, shape = 1, position = position_dodge(width = 0.2)) +
  geom_vline(data=desc_vt %>% dplyr::filter(country=="Israel"), aes(xintercept = 2009L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_vt %>% dplyr::filter(country=="Israel"), aes(xintercept = 2010L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_vt %>% dplyr::filter(country=="South Africa"), aes(xintercept = 2009L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_vt %>% dplyr::filter(country=="South Africa"), aes(xintercept = 2011L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_vt %>% dplyr::filter(country=="Brazil"), aes(xintercept = 2010L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_vt %>% dplyr::filter(country=="Brazil"), aes(xintercept = 2016L), linetype = "dotted", linewidth = 1.5, color = 'blue') +
  geom_vline(data=desc_vt %>% dplyr::filter(country=="Malawi"), aes(xintercept = 2011L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  theme_bw(base_size = 20, base_family = "American typewriter") +
  scale_color_manual(values = paly) +
  facet_grid(source ~ country, scales = "free") +
  guides(color = guide_legend(title = "PCV type: "), size = guide_legend(title = "Sample size: ")) +
  scale_x_continuous(breaks = seq(2005, 2019, 2), limits = c(2004, 2020)) +
  labs(title = "", x = "year of sample collection", y = "") + 
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16), legend.position = "right", legend.title = element_text(size = 16, face="bold"), legend.key.size = unit(1.2,"cm")) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30), strip.background = element_rect(fill = "gray90")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 3))

#save combined plots
ggsave(here("output", "ipd_carprevVT.png"),
       plot = A, 
       width = 27, height = 13, unit = "in", dpi = 300)

#store estimates
write_csv(desc_vt, here("output", "ipd_carprevVT.csv"))


#ggplot of IPD frequency and carriage prevalence of non-vaccine serotypes
B <-
  desc_nvt %>%
  ggplot() +
  geom_line(aes(x = yearc, y = p, color = pcvg, group = pcvg), size = 1, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = yearc, y = p, color = pcvg, size = n), stroke = 2.5, shape = 1, position = position_dodge(width = 0.2)) +
  geom_vline(data=desc_nvt %>% dplyr::filter(country=="Israel"), aes(xintercept = 2009L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_nvt %>% dplyr::filter(country=="Israel"), aes(xintercept = 2010L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_nvt %>% dplyr::filter(country=="South Africa"), aes(xintercept = 2009L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_nvt %>% dplyr::filter(country=="South Africa"), aes(xintercept = 2011L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_nvt %>% dplyr::filter(country=="Brazil"), aes(xintercept = 2010L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  geom_vline(data=desc_nvt %>% dplyr::filter(country=="Brazil"), aes(xintercept = 2016L), linetype = "dotted", linewidth = 1.5, color = 'blue') +
  geom_vline(data=desc_nvt %>% dplyr::filter(country=="Malawi"), aes(xintercept = 2011L), linetype = "dotted", linewidth = 1.5, color = 'black') +
  theme_bw(base_size = 20, base_family = "American typewriter") +
  scale_color_manual(values = palz) +
  facet_grid(source ~ country, scales = "free") +
  guides(color = guide_legend(title = "PCV type: "), size = guide_legend(title = "Sample size: ")) +
  scale_x_continuous(breaks = seq(2005, 2019, 2), limits = c(2004, 2020)) +
  labs(title = "", x = "year of sample collection", y = "") + 
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16), legend.position = "right", legend.title = element_text(size = 16, face="bold"), legend.key.size = unit(1.2,"cm")) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30), strip.background = element_rect(fill = "gray90")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 3))

#save combined plots
ggsave(here("output", "ipd_carprevNVT.png"),
       plot = B, 
       width = 27, height = 13, unit = "in", dpi = 300)

#store estimates
write_csv(desc_nvt, here("output", "ipd_carprevNVT.csv"))





XX <-
  desc_all %>%
  dplyr::filter(source=='Prevalence of carriage') %>%
  
  ggplot(mapping = aes(x = factor(yearc), y = p, color = pcvg, fill = pcvg)) + 
  geom_bar(stat = "identity", size = 0.7) + #color = "black"
  #geom_text(aes(label = n), color = "black", position = position_stack(vjust = 0.5)) +
  #scale_fill_brewer() +
  theme_bw() +
  labs(title = "B", x = "Contact frequency", y = "") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(face = "bold", size = 12, angle = 90, vjust = 0.5, hjust = 0.3), axis.text.y = element_text(face = "bold", size = 12)) +
  theme(plot.title = element_text(size = 22), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) +
  theme(legend.position = "none")









