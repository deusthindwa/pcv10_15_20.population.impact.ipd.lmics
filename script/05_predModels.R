#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
#POST-PCV PREDICTION MODEL
#====================================================================

#assumptions
#-assumes perfectly monitored homogeneous population
#-also work even if IPD and carriage surveillance are imperfectly sensitive as long as don't change post-vax
#-assumes VT are completely eliminated in post-PCV era
#-assumes invasiveness ratio remains constant before and after PCV introduction
#-assumes estimates only apply in mature PCV program, has been in use for a long time
#-assumes IPD and carriage are a representative of the population from which samples are generated
#-assumes the mean carriage duration of VT and NVT is similar allowing calculation of odds VT (VT/NVT)

#required equations
#Q = D/C, Qvt = Dvt/Cvt, Qnvt = Dnvt/Cnvt, invasiveness
#D = Cvt * Qvt + Cnvt * Qnvt, IPD rate in pre-PCV era}
#Cvt' = 0 assumes vaccine type is eliminated after PCV introduction
#y is the prop of pre-PCV VT being replaced by post-PCV NVT
#Cnvt' = Cnvt + y * Cvt
#Qnvt' = Qnvt, NVT invasiveness doesn't change after PCV introduction
#D' = Cvt' * Qvt' + Cnvt' * Qnvt', IPD rate in post-PCV era
#D' = 0 + (Cnvt + y * Cvt) * Qnvt = (Cnvt + y * Cvt) * Qnvt
#IRR = D'/D = [(Cnvt + y * Cvt) * Qnvt]/[Cvt * Qvt + Cnvt * Qnvt], divide by Cnvt, then Qnvt to get below equation
#IRR = [y*c+1]/[d+1], where c = Cvt/Cnvt, d = Dvt/Dnvt

#====================================================================

#Israel
#calculate observed IPD incidence rate ratio (IRR)
bind_cols(
  is_ipdb2009 %>%
    mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", "NVT")) %>%
    group_by(pcv13pfz) %>%
    tally() %>%
    ungroup() %>%
    rename("n1" = "n") %>%
    mutate(fup1 = 1,
           N1 = 7500000,
           incid1 = n1/(N1*fup1)),

  is_ipda2013 %>%
    mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", "NVT")) %>%
    group_by(pcv13pfz) %>%
    tally() %>%
    ungroup() %>%
    rename("n2" = "n") %>%
    mutate(fup2 = 4,
           N2 = 9000000,
           incid2 = n2/(N2*fup2)) %>%
    dplyr::select(everything(), -pcv13pfz)) %>%
  mutate(irr = incid2/incid1)

#calculate expected IPD incidence rate ratio (IRR)
y = 0.68
bind_cols(
  is_carb2009 %>% 
    mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", 
                              if_else(is.na(st), NA_character_, "NVT"))) %>%
    group_by(pcv13pfz) %>%
    tally() %>%
    mutate(p = n/sum(n)) %>%
    dplyr::select(everything(), -n) %>%
    pivot_wider(names_from = pcv13pfz, values_from = p) %>%
    ungroup() %>%
    mutate(c = PCV13/NVT) %>%
    rename("cNVT" = "NVT",  "cPCV13"= "PCV13"),

  is_ipdb2009 %>%
    mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", "NVT")) %>%
    group_by(pcv13pfz) %>%
    tally() %>%
    ungroup() %>%
    pivot_wider(names_from = pcv13pfz, values_from = n) %>%
    ungroup() %>%
    mutate(d = PCV13/NVT)) %>%
  mutate(irr = (y*c+1)/(d+1))

#====================================================================

#estimate invasiveness of NVT (non-PCV13)
data_inv %>%
  mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", 
                            if_else(st == "NT", NA_character_, "NVT"))) %>%
  group_by(pcv13pfz) %>%
  summarise(inv = mean(log_inv)) %>%
  ungroup()

#expected IPD rates post-PCV introduction for VT and NVT
is_ipda2013 %>%
  mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", "NVT")) %>%
  group_by(pcv13pfz) %>%
  tally() %>%
  ungroup() %>%
  rename("n2" = "n") %>%
  mutate(fup2 = 4,
         N2 = 9000000,
         incid2 = n2/(fup2)) %>%
  dplyr::select(everything(), -pcv13pfz)

#serotype carriage data before PCV introduction
is_carb2009 %>% 
  mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "PCV13", 
                            if_else(is.na(st), NA_character_, "NVT"))) %>%
  group_by(pcv13pfz) %>%
  tally() %>%
  mutate(p = n/sum(n)) %>%
  dplyr::select(everything(), -n) %>%
  pivot_wider(names_from = pcv13pfz, values_from = p) %>%
  ungroup() %>%
  mutate(c = PCV13/NVT) %>%
  rename("cNVT" = "NVT",  "cPCV13"= "PCV13")

#estimate post-IPD rates using serotype carriage data pre-PCV




































