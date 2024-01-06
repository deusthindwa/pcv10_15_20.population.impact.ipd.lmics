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
#MODEL VALIDATION
#====================================================================

#define constants and scenarios
fuy1 = 4 #ipd isolates from 2005-2008
fuy2 = 5 #ipd isolates from 2015-2019
bs_samples = 10000 #number of bootstrap samples

#observed IPD incidence rate ratio for PCV13 serotypes (postPCV13 vs prePCV13)
x1 <-
  sa_ipdb2009 %>%
  mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "dVT1", "dNVT1")) %>%
  group_by(pcv13pfz) %>%
  tally() %>%
  pivot_wider(names_from = pcv13pfz, values_from = n) %>%
  ungroup() %>%
  mutate(dTot1 = sum(dNVT1, dVT1)) %>%
  mutate(dTot1 = dTot1/fuy1, dNVT1 = dNVT1/fuy1, dVT1 = dVT1/fuy1)
   
x2 <-
  sa_ipda2015 %>%
  mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "dVT2", "dNVT2")) %>%
  group_by(pcv13pfz) %>%
  tally() %>%
  pivot_wider(names_from = pcv13pfz, values_from = n) %>%
  ungroup() %>%
  mutate(dTot2 = sum(dNVT2, dVT2)) %>%
  mutate(dTot2 = dTot2/fuy2, dNVT2 = dNVT2/fuy2, dVT2 = dVT2/fuy2)

x <-
  bind_cols(
    tibble(dVT1 = replicate(bs_samples, mean(rbinom(n = x1$dTot1, p = x1$dVT1/x1$dTot1, size = 1)*x1$dTot1)), dTot1 = x1$dTot1),
    tibble(dVT2 = replicate(bs_samples, mean(rbinom(n = x2$dTot2, p = x2$dVT2/x2$dTot2, size = 1)*x2$dTot2)), dTot2 = x2$dTot2)) %>%
  mutate(incidVT1 = dVT1/dTot1, incidVT2 = dVT2/dTot2, irr0 = incidVT2/incidVT1)

#model-based IPD incidence rate ratio (using only prePCV13 inferred-carriage and observed IPD data)
y1 <-
sa_pcv %>%
  dplyr::filter(reg == "PCV13", period == "pre-PCV") %>%
  dplyr::mutate(cNVT = NVT*1000, cVT = VT*1000, None = (1000-cNVT-cVT), cTot = sum(cNVT, cVT, None)) %>% #imaginery populaton of 1000 individuals
  dplyr::select(None, cNVT, cVT, cTot) %>%
  dplyr::mutate(None = None/fuy1, cNVT = cNVT/fuy1, cVT = cVT/fuy1, cTot = cTot/fuy1) 

y2 <-
  sa_ipdb2009 %>%
  mutate(pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "dVT", "dNVT")) %>%
  group_by(pcv13pfz) %>%
  tally() %>%
  pivot_wider(names_from = pcv13pfz, values_from = n) %>%
  ungroup() %>%
  mutate(dTot = sum(dNVT, dVT)) %>%
  mutate(dTot = dTot/fuy1, dNVT = dNVT/fuy1, dVT = dVT/fuy1)

y <-
  bind_cols(
    tibble(cVT = replicate(bs_samples, mean(rbinom(n = y1$cTot, p = y1$cVT/y1$cTot, size = 1))),
           cNVT = replicate(bs_samples, mean(rbinom(n = y1$cTot, p = y1$cNVT/y1$cTot, size = 1)))),

    tibble(dVT = replicate(bs_samples, mean(rbinom(n = y2$dTot, p = y2$dVT/y2$dTot, size = 1))), 
           dNVT = replicate(bs_samples, mean(rbinom(n = y2$dTot, p = y2$dNVT/y2$dTot, size = 1)))))

#explore the optimal serotype replacement parameter value over a 3 digit decimal place grid (assume VT elimination)
j = 1
sr = seq(0, 1, by = 0.001) 
sa_err <- tibble(sr = sr, err_diff = rep(NA, 1001))
for (i in sr) {
  sa_err$err_diff[j] = abs(median(x$irr0) - median((i*(y$cVT/y$cNVT)+1)/((y$dVT/y$dNVT)+1)))
  j = j+1
}
sa_err <- 
  sa_err %>% 
  mutate(sr_min = round(sr[which.min(err_diff)], digits = 2),
         country = "South Africa")

#calculate the incidence rate ratio based on predictive function
y <-
  y %>% 
  mutate(irr1 = (0*(cVT/cNVT)+1)/((dVT/dNVT)+1), 
         irr2 = (sa_err$sr_min[1]*(cVT/cNVT)+1)/((dVT/dNVT)+1), 
         irr3 = (1*(cVT/cNVT)+1)/((dVT/dNVT)+1))

sa_validate <-
  bind_rows(
    x %>% dplyr::select(irr0) %>% mutate(sr = "observed IRR", country = "South Africa") %>% rename("irr" = "irr0"),
    y %>% dplyr::select(irr1) %>% mutate(sr = "predicted IRR, no SR", country = "South Africa") %>% rename("irr" = "irr1"),
    y %>% dplyr::select(irr2) %>% mutate(sr = "predicted IRR, baseline SR", country = "South Africa") %>% rename("irr" = "irr2"),
    y %>% dplyr::select(irr3) %>% mutate(sr = "predicted IRR, complete SR", country = "South Africa") %>% rename("irr" = "irr3"))

#====================================================================
#MODEL PREDICTIONS
#====================================================================

#compute expected pcv impact with a prediction model
pcv_carr <-
  bind_cols(
    bind_rows(
      sa_pcv %>%
        dplyr::filter(reg == "PCV10-gsk", period == "post-PCV") %>%
        dplyr::mutate(pcv = "pcv10gsk", cNVT = NVT*2000, cVT = VT*2000, None = (2000-cNVT-cVT), cTot = sum(cNVT, cVT, None)) %>% #imaginery populaton of 2000 individuals
        dplyr::select(pcv, None, cNVT, cVT, cTot) %>%
        dplyr::mutate(None = None/fuy1, cNVT = cNVT/fuy1, cVT = cVT/fuy1, cTot = cTot/fuy1, pcNVT = cNVT/cTot, pcVT = cVT/cTot),
      
      sa_pcv %>%
        dplyr::filter(reg == "PCV10-sii", period == "post-PCV") %>%
        dplyr::mutate(pcv = "pcv10sii", cNVT = NVT*2000, cVT = VT*2000, None = (2000-cNVT-cVT), cTot = sum(cNVT, cVT, None)) %>% #imaginery populaton of 2000 individuals
        dplyr::select(pcv, None, cNVT, cVT, cTot) %>%
        dplyr::mutate(None = None/fuy1, cNVT = cNVT/fuy1, cVT = cVT/fuy1, cTot = cTot/fuy1, pcNVT = cNVT/cTot, pcVT = cVT/cTot),
      
      sa_pcv %>%
        dplyr::filter(reg == "PCV13", period == "post-PCV") %>%
        dplyr::mutate(pcv = "pcv13pfz", cNVT = NVT*2000, cVT = VT*2000, None = (2000-cNVT-cVT), cTot = sum(cNVT, cVT, None)) %>% #imaginery populaton of 2000 individuals
        dplyr::select(pcv, None, cNVT, cVT, cTot) %>%
        dplyr::mutate(None = None/fuy1, cNVT = cNVT/fuy1, cVT = cVT/fuy1, cTot = cTot/fuy1, pcNVT = cNVT/cTot, pcVT = cVT/cTot),
      
      sa_pcv %>%
        dplyr::filter(reg == "PCV15", period == "post-PCV") %>%
        dplyr::mutate(pcv = "pcv15mek", cNVT = NVT*2000, cVT = VT*2000, None = (2000-cNVT-cVT), cTot = sum(cNVT, cVT, None)) %>% #imaginery populaton of 2000 individuals
        dplyr::select(pcv, None, cNVT, cVT, cTot) %>%
        dplyr::mutate(None = None/fuy1, cNVT = cNVT/fuy1, cVT = cVT/fuy1, cTot = cTot/fuy1, pcNVT = cNVT/cTot, pcVT = cVT/cTot),
      
      sa_pcv %>%
        dplyr::filter(reg == "PCV20", period == "post-PCV") %>%
        dplyr::mutate(pcv = "pcv20pfz", cNVT = NVT*2000, cVT = VT*2000, None = (2000-cNVT-cVT), cTot = sum(cNVT, cVT, None)) %>% #imaginery populaton of 2000 individuals
        dplyr::select(pcv, None, cNVT, cVT, cTot) %>%
        dplyr::mutate(None = None/fuy1, cNVT = cNVT/fuy1, cVT = cVT/fuy1, cTot = cTot/fuy1, pcNVT = cNVT/cTot, pcVT = cVT/cTot)),
    
    sa_ipda2015 %>%
      mutate(pcv10sii = if_else(grepl("\\b(1|5|6A|6B|7F|9V|14|19A|19F|23F)\\b", st) == TRUE, "dVT", "dNVT"),
             pcv10gsk = if_else(grepl("\\b(1|4|5|6A|6B|7F|9V|14|18C|19F|23F)\\b", st) == TRUE, "dVT", "dNVT"), #add 6A for cross-protection
             pcv13pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|23F)\\b", st) == TRUE, "dVT", "dNVT"),
             pcv15mek = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "dVT", "dNVT"),
             pcv20pfz = if_else(grepl("\\b(1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F)\\b", st) == TRUE, "dVT", "dNVT")) %>%
      pivot_longer(names_to = "pcv", cols = c(pcv10sii, pcv10gsk, pcv13pfz, pcv15mek, pcv20pfz)) %>%
      ungroup() %>%
      group_by(pcv, value) %>%
      tally() %>%
      pivot_wider(names_from = value, values_from = n) %>%
      mutate(dTot = sum(dNVT, dVT)) %>%
      mutate(dNVT = dNVT/fuy2, dVT = dVT/fuy2, dTot = dTot/fuy2) %>% 
      mutate(pdNVT = dNVT/dTot, pdVT = dVT/dTot) %>%
      ungroup() %>%
      dplyr::select(everything(), -pcv))

#generate random samples around carriage and disease mean values
sa_pcvsamples <- tibble(pcv = c(rep("pcv10gsk", bs_samples), 
                                rep("pcv10sii", bs_samples), 
                                rep("pcv13pfz", bs_samples),
                                rep("pcv15mek", bs_samples), 
                                rep("pcv20pfz", bs_samples)),
                        cVT = NA, cNVT = NA, dVT = NA, dNVT = NA)
k = 1
l = bs_samples
for (i in 1:nrow(pcv_carr)) {
  for (j in k:l) {
    sa_pcvsamples$cVT[j] = mean(rbinom(sum(pcv_carr$cTot[i]), p = pcv_carr$pcVT[i], size = 1))
    sa_pcvsamples$cNVT[j] = mean(rbinom(sum(pcv_carr$cTot[i]), p = pcv_carr$pcNVT[i], size = 1))
    sa_pcvsamples$dVT[j] = mean(rbinom(sum(pcv_carr$dTot[i]), p = pcv_carr$pdVT[i], size = 1)*pcv_carr$dTot[i])
    sa_pcvsamples$dNVT[j] = mean(rbinom(sum(pcv_carr$dTot[i]), p = pcv_carr$pdNVT[i], size = 1)*pcv_carr$dTot[i])
  }
k = l + 1
l = l + bs_samples
}

sa_pcvsamples <-
  sa_pcvsamples %>%
  mutate(irr1 = round((0*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4),
         irr2 = round((sa_err$sr_min[1]*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4),
         irr3 = round((1*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4),
         country = "South Africa")

#summary preventable disease estimates
sa_impactEst <-
  sa_pcvsamples %>% 
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

#next benefit of next genberation PCV over pcv13
sa_netImpact <-
  bind_rows(
    bind_cols(
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv20pfz") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3),
      
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv13pfz") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3) %>%
        dplyr::select(irr1, irr2, irr3) %>%
        dplyr::rename("irr1x" = "irr1", "irr2x" = "irr2", "irr3x" = "irr3")),
    
    bind_cols(
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv15mek") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3),
      
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv13pfz") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3) %>%
        dplyr::select(irr1, irr2, irr3) %>%
        dplyr::rename("irr1x" = "irr1", "irr2x" = "irr2", "irr3x" = "irr3")),
    
    bind_cols(
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv13pfz") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3),
      
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv13pfz") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3) %>%
        dplyr::select(irr1, irr2, irr3) %>%
        dplyr::rename("irr1x" = "irr1", "irr2x" = "irr2", "irr3x" = "irr3")),
    
    bind_cols(
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv10sii") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3),
      
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv13pfz") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3) %>%
        dplyr::select(irr1, irr2, irr3) %>%
        dplyr::rename("irr1x" = "irr1", "irr2x" = "irr2", "irr3x" = "irr3")),
    
    bind_cols(
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv10gsk") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3),
      
      sa_pcvsamples %>%
        dplyr::filter(pcv == "pcv13pfz") %>%
        dplyr::mutate(irr1 = 1-irr1, irr2 = 1-irr2, irr3 = 1-irr3) %>%
        dplyr::select(irr1, irr2, irr3) %>%
        dplyr::rename("irr1x" = "irr1", "irr2x" = "irr2", "irr3x" = "irr3"))) %>%
  
  dplyr::mutate(imp1 = irr1-irr1x, imp2 = irr2-irr2x, imp3 = irr3-irr3x)

#summary of net preventable disease estimates
sa_impactEst <-
  sa_netImpact %>% 
  group_by(pcv) %>%
  summarise(imp1M = quantile(imp1, 0.500),
            imp1L = quantile(imp1, 0.025),
            imp1U = quantile(imp1, 0.975),
            imp2M = quantile(imp2, 0.500),
            imp2L = quantile(imp2, 0.025),
            imp2U = quantile(imp2, 0.975),
            imp3M = quantile(imp3, 0.500),
            imp3L = quantile(imp3, 0.025),
            imp3U = quantile(imp3, 0.975))
