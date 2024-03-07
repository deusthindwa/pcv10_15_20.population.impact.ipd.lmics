#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on pediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
# ipd and carriage data restructuring for plotting
#====================================================================

all_desc <-
  bind_rows(
    is_ipd %>% 
      dplyr::group_by(yearc, st) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(st = if_else(st == "15B/C", "15B/15C", st)) %>%
      dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
      tidyr::separate_rows(., st) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally(n) %>%
      dplyr::mutate(n = n/sum(n), source = "invasive disease", country = "Israel") %>%
      dplyr::ungroup(),
    
    is_car %>% 
      dplyr::group_by(yearc, st) %>%
      dplyr::tally() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(st = if_else(is.na(st), "NA", st), st = if_else(st == "15BC", "15B/15C", st)) %>%
      dplyr::mutate(n = if_else(str_length(st)>3, n/2, n)) %>%
      tidyr::separate_rows(., st) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally(n) %>%
      dplyr::mutate(n = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(st !="NA") %>%
      dplyr::mutate(source = "carriage", country = "Israel"),
    
    #====================================================================
    
    sa_carInf %>%
      dplyr::select(yearc, st, ipd) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally(ipd) %>%
      dplyr::mutate(n = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(source = "invasive disease", country = "South Africa"),
    
    sa_carInf %>%
      dplyr::select(yearc, st, prev) %>%
      dplyr::rename("n" = "prev") %>%
      dplyr::mutate(source = "carriage", country = "South Africa"),
    
    #====================================================================
    
    mw_ipdInf %>%
      dplyr::select(yearc, st, ipd) %>%
      dplyr::group_by(yearc, st) %>%
      dplyr::tally(ipd) %>%
      dplyr::mutate(n = n/sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(source = "invasive disease", country = "Malawi"),
    
    mw_ipdInf %>%
      dplyr::select(yearc, st, prev) %>%
      dplyr::rename("n" = "prev") %>%
      dplyr::mutate(source = "carriage", country = "Malawi")) %>%
  
  #====================================================================

#assign serotype groups
  dplyr::mutate(pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV20", "NVT"),
                pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "PCV15", "NVT"),
                pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "PCV13", "NVT"),
                pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "PCV10", "NVT"))

#====================================================================

all_desc <- 
  bind_rows(
    all_desc %>%
      dplyr::select(everything(), -pcv20pfz, -pcv15mek, -pcv13pfz) %>%
      dplyr::group_by(country, yearc, source, pcv10sii) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::rename("pcvg" = "pcv10sii"),
    
    all_desc %>%
      dplyr::select(everything(), -pcv20pfz, -pcv15mek, -pcv10sii) %>%
      dplyr::group_by(country, yearc, source, pcv13pfz) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::rename("pcvg" = "pcv13pfz"),
    
    all_desc %>%
      dplyr::select(everything(), -pcv20pfz, -pcv13pfz, -pcv10sii) %>%
      dplyr::group_by(country, yearc, source, pcv15mek) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::rename("pcvg" = "pcv15mek"),
    
    all_desc %>%
      dplyr::select(everything(), -pcv15mek, -pcv13pfz, -pcv10sii) %>%
      dplyr::group_by(country, yearc, source, pcv20pfz) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::rename("pcvg" = "pcv20pfz")) %>%
  dplyr::filter(pcvg != "NVT")

#====================================================================

#set the colors
palx <- c("PCV10" = "#84967C", "PCV13" = "#D16A5F", "PCV15" = "#77DAD9", "PCV20" = "#E3B45B")


#ggplot of IPD frequency and carriage prevalence
A <-
  ggplot() +
  geom_line(data = all_desc, aes(x = yearc, y = n, color = pcvg), size = 2) +
  geom_point(data = all_desc, aes(x = yearc, y = n, color = pcvg), size = 2, stroke = 3, shape = 4) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  scale_color_manual(values = palx) +
  facet_wrap(source ~ country, scales = "free") +
  labs(title = "", x = "year of sampling or sample isolation", y = "frequency of vaccine-serotypegroup (or prevalence of vaccine-serotype group)") + 
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16), legend.position = "right", legend.title = element_text(size = 0), legend.key.size = unit(1.2,"cm")) +
  theme(strip.text.x = element_text(size = 30), strip.text.y = element_text(size = 30), strip.background = element_rect(fill = "gray90")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "ipd_carprev.png"),
       plot = (A), 
       width = 20, height = 15, unit = "in", dpi = 300)
