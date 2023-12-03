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
#====================================================================

#define constants and scenarios
fuy1 = 0.5; popn1 = 7500000 #2009 follow up time and population size
fuy2 = 4; popn2 = 9000000 #2013+ follow up time and population size
bs_samples = 10000 #bootstrap sampling

#observed IPD incidence rate ratio for PCV13 serotypes (postPCV13 vs prePCV13)
x1 <-
  is_ipdb2009 %>%
    mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "dVT1", "dNVT1")) %>%
    group_by(pcv13pfz) %>%
    tally() %>%
    pivot_wider(names_from = pcv13pfz, values_from = n) %>%
    ungroup() %>%
    mutate(dTot1 = sum(dNVT1, dVT1))
    
x2 <-
  is_ipda2013 %>%
  mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "dVT2", "dNVT2")) %>%
  group_by(pcv13pfz) %>%
  tally() %>%
  pivot_wider(names_from = pcv13pfz, values_from = n) %>%
  ungroup() %>%
  mutate(dTot2 = sum(dNVT2, dVT2))

x <-
bind_cols(
  tibble(dVT1 = replicate(bs_samples, mean(rbinom(n = x1$dVT1, p = x1$dVT1/x1$dTot1, size = 1)*x1$dTot1)), dTot1 = x1$dTot1),
  tibble(dVT2 = replicate(bs_samples, mean(rbinom(n = x2$dVT2, p = x2$dVT2/x2$dTot2, size = 1)*x2$dTot2)), dTot2 = x2$dTot2)) %>%
  mutate(incidVT1 = dVT1/(popn1*fuy1), incidVT2 = dVT2/(popn2*fuy2), irr0 = incidVT2/incidVT1)

#model-based IPD incidence rate ratio (using only prePCV13 carriage and IPD data)
y1 <-
  is_carb2009 %>% 
    mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "cVT", if_else(is.na(st), "None", "cNVT"))) %>%
    group_by(pcv13pfz) %>%
    tally() %>%
    pivot_wider(names_from = pcv13pfz, values_from = n) %>%
    ungroup() %>%
    mutate(cTot = sum(cNVT, cVT, None),
           None = 1000*None, #scale up numbers by 1000 to ensure carriage prevalence simulation >0
           cNVT = 1000*cNVT,
           cVT = 1000*cVT,
           cTot = 1000*cTot)

y2 <-
  is_ipdb2009 %>%
  mutate(pcv13pfz = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F", st) == TRUE, "dVT", "dNVT")) %>%
  group_by(pcv13pfz) %>%
  tally() %>%
  pivot_wider(names_from = pcv13pfz, values_from = n) %>%
  ungroup() %>%
  mutate(dTot = sum(dNVT, dVT),
         dNVT = 1000*dNVT, #scale up numbers by 1000 to ensure disease prevalence simulation >0
         dVT = 1000*dVT,
         dTot = 1000*dTot)

y <-
  bind_cols(
    tibble(cVT = replicate(bs_samples, mean(rbinom(n = y1$cVT, p = y1$cVT/y1$cTot, size = 1))),
           cNVT = replicate(bs_samples, mean(rbinom(n = y1$cNVT, p = y1$cNVT/y1$cTot, size = 1)))),
    tibble(dVT = replicate(bs_samples, mean(rbinom(n = y2$dVT, p = y2$dVT/y2$dTot, size = 1))), 
           dNVT = replicate(bs_samples, mean(rbinom(n = y2$dNVT, p = y2$dNVT/y2$dTot, size = 1))))) %>%
  mutate(cNVT = if_else(cNVT == 0, 0.05, if_else(cNVT == 1, 0.999, cNVT)),
         dNVT = if_else(dNVT == 0, 0.05*y2$dTot, if_else(dNVT == 1, 0.999*y2$dTot, dNVT)))

#explore the optimal serotype replacement parameter value over a 3 digit decimal place grid (assume VT elimination)
j = 1
sr = seq(0, 1, by = 0.001) 
err_DS <- tibble(sr = sr, err_diff = rep(NA, 1001), err_diffL = rep(NA, 1001), err_diffU = rep(NA, 1001))
for (i in sr) {
  err_DS$err_diff[j] = abs(median(x$irr0) - median((i*(y$cVT/y$cNVT)+1)/((y$dVT/y$dNVT)+1)))
  err_DS$err_diffL[j] = abs(quantile(x$irr0, 0.025) - quantile((i*(y$cVT/y$cNVT)+1)/((y$dVT/y$dNVT)+1), 0.025))
  err_DS$err_diffU[j] = abs(quantile(x$irr0, 0.975) - quantile((i*(y$cVT/y$cNVT)+1)/((y$dVT/y$dNVT)+1), 0.975))
  j = j+1
}
err_DS <- err_DS %>% mutate(sr_min = round(sr[which.min(err_diff)], digits = 2),
                            sr_minL = round(sr[which.min(err_diffL)], digits = 2),
                            sr_minU = round(sr[which.min(err_diffU)], digits = 2),
                            country = "Israel")

#plot serotype replacement scenarios (post-VT eliminated)
A <-
err_DS %>%
  ggplot() +
  geom_line(aes(x = sr, y = err_diff), lty = "twodash", size = 1) +
  geom_point(aes(x = sr[which.min(err_diff)], y = min(err_diff)), size = 2, stroke = 1, shape = 4, color = "black") +
  geom_text(aes(x = sr[which.min(err_diff)], y = min(err_diff), label = paste0(sr_min, " (", sr_minU,"-", sr_minL, ")")), size = 4, angle = "0", vjust = 0.3, hjust = -0.1, fontface = "bold") +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "proportion of NVT replacing VT carriage", y = "predicted vs observed median IRR error") + 
  facet_grid(.~country, scales = "free_y") +
  scale_x_continuous(limit = c(0, 1), breaks = seq(0, 1, 0.25)) + 
  theme(strip.text.x = element_text(size = 0), strip.background = element_rect(fill = "gray90")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#calculate the incidence rate ratio based on predictive function
y <-
  y %>% mutate(irr1 = (0*(cVT/cNVT)+1)/((dVT/dNVT)+1), 
               irr2 = (err_DS1$sr_min[1]*(cVT/cNVT)+1)/((dVT/dNVT)+1), 
               irr3 = (1*(cVT/cNVT)+1)/((dVT/dNVT)+1))

B <-
bind_rows(
x %>% dplyr::select(irr0) %>% mutate(sr = "observed IRR", country = "Israel") %>% rename("irr" = "irr0"),
y %>% dplyr::select(irr1) %>% mutate(sr = "predicted IRR, no SR", country = "Israel") %>% rename("irr" = "irr1"),
y %>% dplyr::select(irr2) %>% mutate(sr = "predicted IRR, estimated SR", country = "Israel", country = "Israel") %>% rename("irr" = "irr2"),
y %>% dplyr::select(irr3) %>% mutate(sr = "predicted IRR, complete SR", country = "Israel") %>% rename("irr" = "irr3")) %>%

ggplot() +
  geom_density_ridges(aes(x = log(irr), y = sr, fill = sr), size = 0.5, alpha = 0.5) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "log_incidence rate ratio", y = "Density") + 
  facet_grid(.~country, scales = "free_y") +
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 0)) +
  guides(fill = guide_legend(title = "Model validation and scenarios\nof serotype replacement (SR)")) +
  theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig6_obsvspred.png"),
       plot = (A | B), 
       width = 16, height = 9, unit = "in", dpi = 300)


#testing conditions that VT are still in circulation
i = 1
stR = seq(0, 1, by = 0.01)
rVT = seq(0, 1, by = 0.01) 
res_DS <- data_frame(stR = rep(NA, 10201), rVT = rep(NA, 10201), irr = rep(NA, 10201))

for (j in stR) {
  
  for (k in rVT) {
    res_DS$stR[i] = j
    res_DS$rVT[i] = k
    res_DS$irr[i] = ((k*(y$dVT/y$dNVT)) + j*(y$cVT/y$cNVT)+1) / ((y$dVT/y$dNVT)+1)
    i = i+1
  }
}

C <-
  res_DS %>%
  ggplot() +
  geom_tile(aes(x = rVT, y = stR, fill = (1-res_DS$irr)), linejoin = "bevel") +
  scale_fill_gradientn("z", colours = terrain.colors(100, rev=TRUE, alpha = 0.8), limits = c(-1.6,1)) +
  geom_point(aes(x = 0, y = err_DS$sr_min[1]), shape = 4, size = 4, stroke = 2) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  labs(title = "", x = "residual VT carriage in mature PCV era", y = "proportion of NVT replacing VT carriage") +
  theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  guides(fill = guide_legend(title = "proportion of\npreventable IPD")) +
  theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig7_rVTstRsens.png"),
       plot = (C), 
       width = 10, height = 10, unit = "in", dpi = 300)

#====================================================================
#====================================================================

#compute expected pcv impact with a prediction model
pcv_carr <-
  bind_cols(
    is_cara2013 %>% 
      mutate(pcv10sii = if_else(grepl("1|5|6A|6B|7F|9V|14|19A|19F|23F", st) == TRUE, "cVT", if_else(is.na(st), "None", "cNVT")),
             pcv10gsk = if_else(grepl("1|4|5|6B|7F|9V|14|18C|19F|23F", st) == TRUE, "cVT", if_else(is.na(st), "None", "cNVT")),
             pcv15mek = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F", st) == TRUE, "cVT", if_else(is.na(st), "None", "cNVT")),
             pcv20pfz = if_else(grepl("1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F", st) == TRUE, "cVT", if_else(is.na(st), "None", "cNVT"))) %>%
      pivot_longer(names_to = "pcv", cols = c(pcv10sii, pcv10gsk, pcv15mek, pcv20pfz)) %>%
      group_by(pcv, value) %>%
      tally() %>%
      pivot_wider(names_from = value, values_from = n) %>%
      ungroup() %>%
      mutate(cNVT = cNVT/fuy2, cVT = cVT/fuy2, None = None/fuy2, cTot = sum(None, cVT, cNVT)/fuy2) %>% #annualise disease cases
      mutate(pcNVT = cNVT/cTot, pcVT = cVT/cTot), 
    
    is_ipda2013 %>%
      mutate(pcv10sii = if_else(grepl("1|5|6A|6B|7F|9V|14|19A|19F|23F", st) == TRUE, "dVT", "dNVT"),
             pcv10gsk = if_else(grepl("1|4|5|6B|7F|9V|14|18C|19F|23F", st) == TRUE, "dVT", "dNVT"),
             pcv15mek = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F", st) == TRUE, "dVT", "dNVT"),
             pcv20pfz = if_else(grepl("1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F", st) == TRUE, "dVT", "dNVT")) %>%
      pivot_longer(names_to = "pcv", cols = c(pcv10sii, pcv10gsk, pcv15mek, pcv20pfz)) %>%
      ungroup() %>%
      group_by(pcv, value) %>%
      tally() %>%
      pivot_wider(names_from = value, values_from = n) %>%
      ungroup() %>%
      mutate(dNVT = dNVT/fuy2, dVT = dVT/fuy2, dTot = sum(dNVT, dVT)/fuy2) %>% #annualise carriage prevalence
      mutate(pdNVT = dNVT/dTot, pdVT = dVT/dTot) %>%
      dplyr::select(everything(), -pcv))

#generate random samples around carriage and disease mean values
pcv_samples <- tibble(pcv = c(rep("pcv10gsk", bs_samples), 
                              rep("pcv10sii", bs_samples), 
                              rep("pcv15mek", bs_samples), 
                              rep("pcv20pfz", bs_samples)),
                      cVT = NA, cNVT = NA, dVT = NA, dNVT = NA)
k = 1
l = bs_samples
for (i in 1:4) {
  for (j in k:l) {
    pcv_samples$cVT[j] = mean(rbinom(sum(pcv_carr$cTot[i]), p = pcv_carr$pcVT[i], size = 1))
    pcv_samples$cNVT[j] = mean(rbinom(sum(pcv_carr$cTot[i]), p = pcv_carr$pcNVT[i], size = 1))
    pcv_samples$dVT[j] = mean(rbinom(sum(pcv_carr$dTot[i]), p = pcv_carr$pdVT[i], size = 1)*dTot)
    pcv_samples$dNVT[j] = mean(rbinom(sum(pcv_carr$dTot[i]), p = pcv_carr$pdNVT[i], size = 1)*dTot)
  }
k = l + 1
l = l + bs_samples
}

pcv_samples <-
  pcv_samples %>%
  mutate(irr1 = round((0*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4),
         irr2 = round((err_DS$sr_min[1]*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4),
         irr3 = round((1*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4),
         country = "Israel")

#summary preventable disease estimates
pcv_samples %>% 
  group_by(pcv) %>%
  summarise(irr1M = 1-quantile(irr1, 0.500),#flip the 95%CI koz of subtracting from 1
            irr1L = 1-quantile(irr1, 0.975),
            irr1U = 1-quantile(irr1, 0.025),
            irr2M = 1-quantile(irr2, 0.500),
            irr2L = 1-quantile(irr2, 0.975),
            irr2U = 1-quantile(irr2, 0.025),
            irr3M = 1-quantile(irr3, 0.500),
            irr3L = 1-quantile(irr3, 0.975),
            irr3U = 1-quantile(irr3, 0.025))

#plot preventable disease distributions
C <-
pcv_samples %>%
  ggplot() +
  geom_density(aes(x = 1-irr1, group = pcv, fill = "no SR"), size = 0.6, alpha = 0.4) +
  geom_density(aes(x = 1-irr2, group = pcv, fill = "estimated SR"), size = 0.6, alpha = 0.4) +
  geom_density(aes(x = 1-irr3, group = pcv, fill = "complete SR"), size = 0.6, alpha = 0.4) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(country~pcv, scales = "free_y") +
  scale_fill_manual(name = "scenarios of serotype\nreplacement (SR)", values = c("no SR" = "blue", "estimated SR" = "red", "complete SR" = "green")) +
  labs(title = "", x = "proportion of preventable IPD", y = "density") + 
  scale_x_continuous(limit = c(-1, 1), breaks = seq(-1, 1, 0.2)) + 
  theme(legend.text = element_text(size = 14), legend.position = "right", legend.title = element_text(size = 14), legend.key.size = unit(1.2,"cm")) +
  theme(strip.text.x = element_text(size = 26), strip.text.y = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "sfig8_vaximpact.png"),
       plot = (C), 
       width = 22, height = 7, unit = "in", dpi = 300)

