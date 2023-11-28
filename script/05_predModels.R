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

#Israel pcv13 impact
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
  
  mutate(irr1 = (0*c+1)/(d+1),
         irr2 = (0.679*c+1)/(d+1),
         irr3 = (1*c+1)/(d+1))

#====================================================================

#Israel pcv10sii impact
#calculate expected IPD incidence rate ratio (IRR)
  bind_cols(
    is_cara2013 %>% 
      mutate(pcv10sii = if_else(grepl("1|5|6A|6B|7F|9V|14|19A|19F|23F", st) == TRUE, "PCV10sii", 
                                if_else(is.na(st), NA_character_, "NVT"))) %>%
      group_by(pcv10sii) %>%
      tally() %>%
      mutate(p = n/sum(n)) %>%
      dplyr::select(everything(), -n) %>%
      pivot_wider(names_from = pcv10sii, values_from = p) %>%
      ungroup() %>%
      mutate(c = PCV10sii/NVT) %>%
      rename("cNVT" = "NVT",  "cPCV10sii"= "PCV10sii"),
    
    is_ipda2013 %>%
      mutate(pcv10sii = if_else(grepl("1|5|6A|6B|7F|9V|14|19A|19F|23F", st) == TRUE, "PCV10sii", "NVT")) %>%
      group_by(pcv10sii) %>%
      tally() %>%
      ungroup() %>%
      pivot_wider(names_from = pcv10sii, values_from = n) %>%
      ungroup() %>%
      mutate(d = PCV10sii/NVT)) %>%
  
  mutate(irr1 = (0*c+1)/(d+1),
         irr2 = (0.679*c+1)/(d+1),
         irr3 = (1*c+1)/(d+1))

#====================================================================

#Israel pcv10gsk impact
#calculate expected IPD incidence rate ratio (IRR)
  bind_cols(
    is_cara2013 %>% 
      mutate(pcv10gsk = if_else(grepl("1|4|5|6B|7F|9V|14|18C|19F|23F", st) == TRUE, "PCV10gsk", 
                                if_else(is.na(st), NA_character_, "NVT"))) %>%
      group_by(pcv10gsk) %>%
      tally() %>%
      mutate(p = n/sum(n)) %>%
      dplyr::select(everything(), -n) %>%
      pivot_wider(names_from = pcv10gsk, values_from = p) %>%
      ungroup() %>%
      mutate(c = PCV10gsk/NVT) %>%
      rename("cNVT" = "NVT",  "cPCV10gsk"= "PCV10gsk"),
    
    is_ipda2013 %>%
      mutate(pcv10gsk = if_else(grepl("1|4|5|6B|7F|9V|14|18C|19F|23F", st) == TRUE, "PCV10gsk", "NVT")) %>%
      group_by(pcv10gsk) %>%
      tally() %>%
      ungroup() %>%
      pivot_wider(names_from = pcv10gsk, values_from = n) %>%
      ungroup() %>%
      mutate(d = PCV10gsk/NVT)) %>%
  
  mutate(irr1 = (0*c+1)/(d+1),
         irr2 = (0.679*c+1)/(d+1),
         irr3 = (1*c+1)/(d+1))

#====================================================================

#Israel pcv15 impact
#calculate expected IPD incidence rate ratio (IRR)
  bind_cols(
    is_cara2013 %>% 
      mutate(pcv15mek = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F", st) == TRUE, "PCV15", 
                                if_else(is.na(st), NA_character_, "NVT"))) %>%
      group_by(pcv15mek) %>%
      tally() %>%
      mutate(p = n/sum(n)) %>%
      dplyr::select(everything(), -n) %>%
      pivot_wider(names_from = pcv15mek, values_from = p) %>%
      ungroup() %>%
      mutate(c = PCV15/NVT) %>%
      rename("cNVT" = "NVT",  "cPCV15"= "PCV15"),
    
    is_ipda2013 %>%
      mutate(pcv15mek = if_else(grepl("1|3|4|5|6A|6B|7F|9V|14|18C|19A|19F|22F|23F|33F", st) == TRUE, "PCV15", "NVT")) %>%
      group_by(pcv15mek) %>%
      tally() %>%
      ungroup() %>%
      pivot_wider(names_from = pcv15mek, values_from = n) %>%
      ungroup() %>%
      mutate(d = PCV15/NVT)) %>%
  
  mutate(irr1 = (0*c+1)/(d+1),
         irr2 = (0.679*c+1)/(d+1),
         irr3 = (1*c+1)/(d+1))

#====================================================================

#Israel pcv20 impact
#calculate expected IPD incidence rate ratio (IRR)
  bind_cols(
    is_cara2013 %>% 
      mutate(pcv20pfz = if_else(grepl("1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F", st) == TRUE, "PCV20", 
                                if_else(is.na(st), NA_character_, "NVT"))) %>%
      group_by(pcv20pfz) %>%
      tally() %>%
      mutate(p = n/sum(n)) %>%
      dplyr::select(everything(), -n) %>%
      pivot_wider(names_from = pcv20pfz, values_from = p) %>%
      ungroup() %>%
      mutate(c = PCV20/NVT) %>%
      rename("cNVT" = "NVT",  "cPCV20"= "PCV20"),
    
    is_ipda2013 %>%
      mutate(pcv20pfz = if_else(grepl("1|3|4|5|6A|6B|7F|8|9V|10A|11A|12F|14|15B|18C|19A|19F|22F|23F|33F", st) == TRUE, "PCV20", "NVT")) %>%
      group_by(pcv20pfz) %>%
      tally() %>%
      ungroup() %>%
      pivot_wider(names_from = pcv20pfz, values_from = n) %>%
      ungroup() %>%
      mutate(d = PCV20/NVT)) %>%
  
  mutate(irr1 = (0*c+1)/(d+1),
         irr2 = (0.679*c+1)/(d+1),
         irr3 = (1*c+1)/(d+1))

#====================================================================
#====================================================================
fuy = 4
is_pcv10sii <-
bind_cols(
is_cara2013 %>% 
  mutate(pcv10sii = if_else(grepl("1|5|6A|6B|7F|9V|14|19A|19F|23F", st) == TRUE, "cVT", if_else(is.na(st), "None", "cNVT"))) %>%
  group_by(pcv10sii) %>%
  tally() %>%
  pivot_wider(names_from = pcv10sii, values_from = n) %>%
  ungroup() %>%
  mutate(cNVT = cNVT/fuy, cVT = cVT/fuy, None = None/fuy, cTot = sum(None, cVT, cNVT)), #annualise carriage prevalence

is_ipda2013 %>%
  mutate(pcv10sii = if_else(grepl("1|5|6A|6B|7F|9V|14|19A|19F|23F", st) == TRUE, "dVT", "dNVT")) %>%
  group_by(pcv10sii) %>%
  tally() %>%
  pivot_wider(names_from = pcv10sii, values_from = n) %>%
  ungroup() %>%
  mutate(dNVT = dNVT/fuy, dVT = dVT/fuy, dTot = sum(dVT, dNVT)) #annualise disease cases
)

#set the mean carriage prevalence and disease incidence
cVT_mean = is_pcv10sii$cVT/(is_pcv10sii$cTot)
cNVT_mean = is_pcv10sii$cNVT/(is_pcv10sii$cTot)
#c = cVT_mean/cNVT_mean
dVT_mean = is_pcv10sii$dVT
dNVT_mean = is_pcv10sii$dNVT
dTot = is_pcv10sii$dTot
#d = dVT_mean/dNVT_mean

#serotype replacement status
y1 = 0     #no replacement 
y2 = 0.679 #estimated/baseline replacement
y3 = 1     #complete replacement 

#bootstrap sampling
bs_samples = 1000

simul_ps <-  
  tibble(cVT = replicate(bs_samples, mean(rbinom(sum(is_pcv10sii[,3]), p = cVT_mean, size = 1))),
         cNVT = replicate(bs_samples, mean(rbinom(sum(is_pcv10sii[,2]), p = cNVT_mean, size = 1))),
         dVT = replicate(bs_samples, mean(rbinom(sum(is_pcv10sii[,6]), p = dVT_mean/(dTot), size = 1)*dTot)),
         dNVT = replicate(bs_samples, mean(rbinom(sum(is_pcv10sii[,5]), p = dNVT_mean/(dTot), size = 1)*dTot))) %>%
  mutate(irr1 = round((y1*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4),
         irr2 = round((y2*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4),
         irr3 = round((y3*(cVT/cNVT) + 1) / ((dVT/dNVT) + 1), 4))

simul_ps %>%
ggplot(aes(x = cVT, y = dVT/dTot, fill = irr1)) +
  geom_tile()


#bootstrap sampling
for(b in 1:bs_samples){
  
  #carriage
  cVT = dVT*0+(rbinom(n = length(cN_pos), size = cN_pos, prob = cVT_mean)/cN_pos)
  for(c in names(cVT_posterior)){
    cVT[c] = sample(cVT_posterior[,c],1)
  }
  cVT[which(cVT ==0)] = 0.001
  cVT[which(cVT ==1)] = 0.999
  cNVT = 1-cVT
  cVTNVT = cVT/cNVT
  
  cVT = dVT*0+(rbinom(n = length(cN_pos), size = cN_pos, prob = cVT_mean)/cN_pos)
  for(c in names(cVT_posterior)){
    cVT[c] = sample(cVT_posterior[,c],1)
  }
  cVT[which(cVT ==0)] = 0.001
  cVT[which(cVT ==1)] = 0.999
  cNVT = 1-cVT
  cVTNVT = cVT/cNVT
  
  #disease
  dVT_temp = dVT*0+(rbinom(n = length(dVT), size = dVT_sample, prob = dVT/100)/dVT_sample)
  dVT_temp[which(dVT_temp==0)] = 0.001
  dVT_temp[which(dVT_temp==1)] = 0.999
  dNVT_temp = 1-dVT_temp
  dVTNVT_temp = dVT_temp/dNVT_temp
  
  #make prediction
  my_pred_vacc_eff_disease_boots[b,] = predicted_vacc_eff_Disease(cVTNVT, dVTNVT_temp, lambda)
  noreplace_pred_vacc_eff_disease_boots[b,] = predicted_vacc_eff_Disease(cVTNVT, dVTNVT_temp, 0)				
}

my_pred_IRR_disease = round(1-predicted_vacc_eff_Disease(cVTNVT_mean, dVTNVT, lambda), 2)
my_pred_IRR_disease_low = round(apply(1-my_pred_vacc_eff_disease_boots, 2, quantile,probs = 0.025,na.rm = T), 2)
my_pred_IRR_disease_high = round(apply(1-my_pred_vacc_eff_disease_boots, 2, quantile,probs = 0.975, na.rm = T), 2)

noreplace_pred_IRR_disease = round(1-predicted_vacc_eff_Disease(cVTNVT_mean, dVTNVT, 0), 2)
noreplace_pred_IRR_disease_low = round(apply(1-noreplace_pred_vacc_eff_disease_boots, 2, quantile, probs = 0.025, na.rm = T), 2)
noreplace_pred_IRR_disease_high = round(apply(1-noreplace_pred_vacc_eff_disease_boots, 2, quantile,probs = 0.975, na.rm = T), 2)

