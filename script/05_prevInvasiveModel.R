#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

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
