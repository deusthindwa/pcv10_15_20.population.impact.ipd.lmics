#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries
 
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
  dplyr::select(everything(), -pcv13pfz)
) %>%
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



#fit a negative-binomial model of carriage prevalence and estimated invasiveness
model1 <- MASS::glm.nb(nipd ~ log_prevcarr + log_inv, 
                       data = data_Bogota0,
                       control = glm.control(maxit = 25, trace = T), 
                       link = log)

#add model estimates to the dataset
data_Bogota0 <- 
  data_Bogota0 %>% 
  dplyr::mutate(fit_ipd = model1$coefficients[1] + model1$coefficients[2]*log_prevcarr + model1$coefficients[3]*log_inv)



#====================================================================

#model
#expected_IPD = constant + log_carriage_prevalence + log_est_invasiveness
#use this model to predict expected IPD for each country's carriage data (?using estimated par for one country)

#combine carriage, ipd and invasiveness datasets
data_Bogota0 <- 
  data_all %>% dplyr::filter(country == "Bogota", phase == "pre-pcv") %>%
  left_join(data_inv)

#say we are counting the # of IPD isolates; if we increase the the total population under surveillance 
#would this increase in the number of IPD isolates observed? If yes, the total population should be 
#considered as offset, otherwise include log_total_population as covariate in the model.
#with offset, the change in the rate of the outcome per the offset variable 
#for changes in the independent variables. Use incident rate ratios ....IRR.

#fit a negative-binomial model of carriage prevalence and estimated invasiveness
model1 <- MASS::glm.nb(nipd ~ log_prevcarr + log_inv, 
                       data = data_Bogota0,
                       control = glm.control(maxit = 25, trace = T), 
                       link = log)

#add model estimates to the dataset
data_Bogota0 <- 
  data_Bogota0 %>% 
  dplyr::mutate(fit_ipd = model1$coefficients[1] + model1$coefficients[2]*log_prevcarr + model1$coefficients[3]*log_inv)

#====================================================================

#visualize relationship between predicted and observed IPD
spearmanCI <- function(x, y, alpha = 0.05){ #function for spearman rank correlation coefficient 95%CIs
  rs <- cor(x, y, method = "spearman", use = "complete.obs")
  n <- sum(complete.cases(x, y))
  round(sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2))), digits = 2)
}

A <-
  data_Bogota0 %>%
  dplyr::filter(!is.na(fit_ipd)) %>%
  mutate(spcoef = round(cor(fit_ipd, nipd,  method = "spearman"), digits = 2)) %>%
  
  ggplot(aes(x = fit_ipd, y = log(nipd), color = st)) +
  geom_point(size = 2.5, shape = 1, stroke = 2) +
  geom_text(aes(label = st), size = 3, angle = 45, fontface = "bold", vjust = 2, hjust = 0.5) +
  geom_text(aes(x = 3, y = 6, label = paste0("Spearman, p = ", spcoef, ", \n[95%CI = ", spearmanCI(fit_ipd, nipd)[1], "-", spearmanCI(fit_ipd, nipd)[2], "]")), color = "black", size = 4) +
  expand_limits(x = c(-1,7), y = c(-1,7)) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(A)", x = "invasive disease (model estimate)", y = "log_invasive disease (observed)") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

 #====================================================================

#visualize relationship between serotype rank order and Obs-Pred cumulative disease
 data_Bogota0 <- 
   data_Bogota0 %>%
   dplyr::filter(!is.na(fit_ipd)) %>%
   dplyr::mutate(nipdR = dplyr::min_rank(nipd),
                 fit_ipdR = dplyr::min_rank(fit_ipd)) %>%
   dplyr::arrange(fit_ipdR) %>%
   dplyr::mutate(fit_ipdCum = cumsum(fit_ipd)/sum(fit_ipd),
                 st_contrib = fit_ipd/sum(fit_ipd)) %>%
   dplyr::arrange(nipdR) %>%
   dplyr::mutate(nipdCum = cumsum(nipd)/sum(nipd))

B <-
  data_Bogota0 %>%
  ggplot() +
  geom_point(aes(x = fit_ipdCum, y = fit_ipdR), color = "darkblue", size = 2.5, shape = 1, stroke = 2) +
  geom_line(aes(x = nipdCum, y = nipdR ), color = "darkred", size = 1, linetype = 1) +
  geom_line(aes(x = fit_ipdCum, y = fit_ipdR), color = "darkblue", size = 1, linetype = 1) +
  geom_abline(intercept = 0, slope = 26, color = "black", linetype = 2, size = 1) + 
  geom_text(aes(x = fit_ipdCum, y = fit_ipdR, label = st), size = 3, angle = 0, fontface = "bold", vjust = 2, hjust = -0.5) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(B)", x = "Cumulative proportion of  IPD\n(Observed and Predicted)", y = "Serotype rank") + 
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  #theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))

#====================================================================

#save combined plots
ggsave(here("output", "fig2_serotypeRank.png"),
       plot = (A | B), 
       width = 12, height = 8, unit = "in", dpi = 300)

#delete all plot objects
rm(list = grep("data_|model", ls(), value = TRUE, invert = TRUE))
