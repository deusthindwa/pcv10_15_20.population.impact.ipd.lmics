#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries
 
#====================================================================
#INVASIVENESS DATA FROM SOUTH AFRICA 
#====================================================================

#import invasiveness data
#invasivenes <- rio::import("https://raw.githubusercontent.com/weinbergerlab/Invasiveness_Navajo/main/Results/mcmc_invasive_single_stage.csv")
invasivenes <-
  rio::import(here("data", "invasiveness_global.csv")) %>%
  dplyr::select(everything(), -V1, -st.index) %>%
  dplyr::rename("log_inv" = "log.inv.age1", "log_var" = "log.inv.prec.age1") %>%
  dplyr::filter(st != "NT") %>%
  dplyr::mutate(pcv20pfz = if_else(grepl(pcv20pfz1, st) == TRUE, "PCV20", "NVT"),
                pcv15mek = if_else(grepl(pcv15mek1, st) == TRUE, "PCV15", "NVT"),
                pcv13pfz = if_else(grepl(pcv13pfz1, st) == TRUE, "PCV13", "NVT"),
                pcv10sii = if_else(grepl(pcv10sii1, st) == TRUE, "PCV10-sii", "NVT"),
                pcv10gsk = if_else(grepl(pcv10gsk1, st) == TRUE, "PCV10-gsk", "NVT"),
                pcv7pfz = if_else(grepl(pcv7pfz1, st) == TRUE, "PCV7", "NVT"))

#calculate weighted mean invasiveness of VT and NVT for serotypes with unknown invasiveness
inv_pcv20pfz <- invasivenes %>% group_by(pcv20pfz) %>% mutate(inv = weighted.mean(log_inv, log_var)) %>% ungroup() %>% distinct(pcv20pfz, inv)
inv_pcv15mek <- invasivenes %>% group_by(pcv15mek) %>% mutate(inv = weighted.mean(log_inv, log_var)) %>% ungroup() %>% distinct(pcv15mek, inv)
inv_pcv13pfz <- invasivenes %>% group_by(pcv13pfz) %>% mutate(inv = weighted.mean(log_inv, log_var)) %>% ungroup() %>% distinct(pcv13pfz, inv)
inv_pcv10sii <- invasivenes %>% group_by(pcv10sii) %>% mutate(inv = weighted.mean(log_inv, log_var)) %>% ungroup() %>% distinct(pcv10sii, inv)
inv_pcv10gsk <- invasivenes %>% group_by(pcv10gsk) %>% mutate(inv = weighted.mean(log_inv, log_var)) %>% ungroup() %>% distinct(pcv10gsk, inv)
inv_pcv7pfz  <- invasivenes %>% group_by(pcv7pfz) %>% mutate(inv = weighted.mean(log_inv, log_var)) %>% ungroup() %>% distinct(pcv7pfz, inv)

#create ipd serotype dataset to match those of invasiveness
#infer carriage data pre-pcv13 introduction in south africa (ipd <- carriage * invasiveness)
sa_carb2009 <-
  sa_ipdb2009 %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  mutate(p = n/sum(n), ipd = n*p, log_ipd = log(ipd)) %>%
  dplyr::select(st, log_ipd)

sa_carb2009 <-
  left_join(sa_carb2009, invasivenes %>% 
              dplyr::select(st, log_inv, log_var, pcv13pfz)) %>%
  
  mutate(pcv13pfz = if_else(grepl(pcv13pfz1, st) == TRUE, "PCV13", "NVT"),
         log_inv = if_else(pcv13pfz == "PCV13" & is.na(log_inv), inv_pcv13pfz$inv[1],
                           if_else(pcv13pfz == "NVT" & is.na(log_inv), inv_pcv13pfz$inv[2], log_inv))) %>%
  
  dplyr::select(pcv13pfz, everything(), -log_var)

sa_carb2009 <- 
  sa_carb2009 %>%
  mutate(log_carr = log_ipd - log_inv, carr = exp(log_carr))

#create ipd serotype dataset to match those of invasiveness
#infer carriage data ppst-pcv13 introduction in south africa (ipd <- carriage * invasiveness)
sa_cara2015 <-
  sa_ipda2015 %>%
  group_by(st) %>%
  tally() %>%
  ungroup() %>%
  mutate(p = n/sum(n), ipd = n*p, log_ipd = log(ipd)) %>%
  dplyr::select(st, log_ipd)

sa_cara2015 <-
  left_join(sa_cara2015, invasivenes %>% 
              dplyr::select(st, log_inv, log_var, pcv13pfz)) %>%
  
  mutate(pcv13pfz = if_else(grepl(pcv13pfz1, st) == TRUE, "PCV13", "NVT"),
         log_inv = if_else(pcv13pfz == "PCV13" & is.na(log_inv), inv_pcv13pfz$inv[1],
                           if_else(pcv13pfz == "NVT" & is.na(log_inv), inv_pcv13pfz$inv[2], log_inv))) %>%
  
  dplyr::select(pcv13pfz, everything(), -log_var)

sa_cara2015 <- 
  sa_cara2015 %>%
  mutate(log_carr = log_ipd - log_inv, carr = exp(log_carr))

#====================================================================
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








#set color for all plots
cols <- c("IPD"="darkviolet", "carriage" = "chartreuse4")

#invasiveness descriptive plot
A <-
  data_inv %>% 
  ggplot() +
  geom_point(aes(x = reorder(st, log_inv), y = log_inv, color = log_inv), size = 2.5, shape = 5, stroke = 2) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  labs(title = "(A)", x = "pneumococcal serotype", y = "log_invasiveness") + 
  scale_color_distiller(palette = "Reds", direction = 1) +
  theme(axis.text.y = element_text(face = "bold", size = 10)) + 
  theme(legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 0)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Brazil
B <-
  data_all %>%
  dplyr::filter(country == "Bogota") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1))) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Bogota", values = cols) + 
  facet_wrap(.~ factor(phase, levels = c("pre-pcv", "post-pcv")), scales = "free_x") +
  labs(title = "(B)", x = "invasive disease isolates", y = "pneumococcal serotype") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Kenya
C <-
  data_all %>%
  dplyr::filter(country == "Netherlands") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1))) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Netherlands", values = cols) + 
  labs(title = "(C)", x = "invasive disease isolates", y = "") + 
  facet_wrap(.~ factor(phase, levels = c("pre-pcv", "post-pcv")), scales = "free_x") +
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for Malawi
D <-
  data_all %>%
  dplyr::filter(country == "Caracas") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1))) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Caracas", values = cols) + 
  facet_wrap(.~ factor(phase, levels = c("pre-pcv", "post-pcv")), scales = "free_x") +
  labs(title = "(D)", x = "invasive disease isolates", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

#carriage description for South Africa
E <-
  data_all %>%
  dplyr::filter(country == "Czech") %>%
  ggplot() +
  geom_point(aes(y = reorder(st, nipd), x = nipd, color = "IPD"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  geom_point(aes(y = reorder(st, prevcarr), x = prevcarr*800, color = "carriage"), size = 2.5, shape = 4, stroke = 2, position = position_dodge2(width = 0.5), stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis( ~ . /800, name = "carriage prevalence", labels = scales::percent_format(accuracy = 1))) +
  theme_bw(base_size = 14, base_family = "American Typewriter") +
  scale_color_manual(name = "Czech", values = cols) + 
  facet_wrap(.~ factor(phase, levels = c("pre-pcv", "post-pcv")), scales = "free_x") +
  labs(title = "(E)", x = "invasive disease isolates", y = "") + 
  theme(axis.text.y = element_text(face = "bold", size = 10), axis.title.x = element_text(colour = "darkviolet"), axis.title.x.top = element_text(colour = "chartreuse4")) + 
  theme(legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#====================================================================

# #save combined plots
# ggsave(here("output", "fig1_carripdDesc.png"),
#        plot = (A / (B | C | D | E | plot_layout(widths = c(2,2,1,1))) + plot_layout(heights = c(1,2.3))), 
#        width = 18, height = 12, unit = "in", dpi = 300)
# 
# #delete all plot objects
# rm(list = grep("data_", ls(), value = TRUE, invert = TRUE))
# 
# 
# 
# model1 <- MASS::glm.nb(nipd ~ log_prevcarr + log_inv,
#                        data = mw_carb2011,
#                        control = glm.control(maxit = 25, trace = T),
#                        link = log)
# 
# #add model estimates to the dataset
# mw_carb2011 <-
#   mw_carb2011 %>%
#   dplyr::mutate(fit_ipd = model1$coefficients[1] + model1$coefficients[2]*log_prevcarr + model1$coefficients[3]*log_inv)
# 
# #visualize relationship between predicted and observed IPD
# spearmanCI <- function(x, y, alpha = 0.05){ #function for spearman rank correlation coefficient 95%CIs
#   rs <- cor(x, y, method = "spearman", use = "complete.obs")
#   n <- sum(complete.cases(x, y))
#   round(sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2))), digits = 2)
# }
# 
# data_Bogota0 %>%
#   dplyr::filter(!is.na(fit_ipd)) %>%
#   mutate(spcoef = round(cor(fit_ipd, nipd,  method = "spearman"), digits = 2)) %>%
#   
#   ggplot(aes(x = fit_ipd, y = log(nipd), color = st)) +
#   geom_point(size = 2.5, shape = 1, stroke = 2) +
#   geom_text(aes(label = st), size = 3, angle = 45, fontface = "bold", vjust = 2, hjust = 0.5) +
#   geom_text(aes(x = 3, y = 6, label = paste0("Spearman, p = ", spcoef, ", \n[95%CI = ", spearmanCI(fit_ipd, nipd)[1], "-", spearmanCI(fit_ipd, nipd)[2], "]")), color = "black", size = 4) +
#   expand_limits(x = c(-1,7), y = c(-1,7)) +
#   theme_bw(base_size = 14, base_family = "American Typewriter") +
#   labs(title = "(A)", x = "invasive disease (model estimate)", y = "log_invasive disease (observed)") +
#   theme(axis.text.y = element_text(face = "bold", size = 10)) +
#   theme(legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 12)) +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))
# 
# 
