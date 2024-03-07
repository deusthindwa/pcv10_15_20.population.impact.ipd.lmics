#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
#MODEL SENSITIVITY TO RESIDUAL CARRIAGE
#====================================================================

#split datasets by their PCV
is_pcvirr10gsk <- is_pcvsamples %>% filter(pcv == "pcv10gsk")
is_pcvirr10sii <- is_pcvsamples %>% filter(pcv == "pcv10sii")
is_pcvirr15mek <- is_pcvsamples %>% filter(pcv == "pcv15mek")
is_pcvirr20pfz <- is_pcvsamples %>% filter(pcv == "pcv20pfz")

#testing conditions that VT are still in circulation in presence of replacement
i = 1
stR = seq(0, 1, by = 0.01)
rVT = seq(0, 1, by = 0.01)
is_resds <- data_frame(stR = rep(NA, 10201), rVT = rep(NA, 10201), pcv10gsk = rep(NA, 10201), pcv10sii = rep(NA, 10201), pcv15mek = rep(NA, 10201), pcv20pfz = rep(NA, 10201))

for (j in stR) {
  
  for (k in rVT) {
    is_resds$stR[i] = j
    is_resds$rVT[i] = k
    is_resds$pcv10gsk[i] = (k*(is_pcvirr10gsk$dVT/is_pcvirr10gsk$dNVT) + (j*(1-k)*(is_pcvirr10gsk$cVT/is_pcvirr10gsk$cNVT) + 1)) / ((is_pcvirr10gsk$dVT/is_pcvirr10gsk$dNVT) + 1)
    is_resds$pcv10sii[i] = (k*(is_pcvirr10sii$dVT/is_pcvirr10sii$dNVT) + (j*(1-k)*(is_pcvirr10sii$cVT/is_pcvirr10sii$cNVT) + 1)) / ((is_pcvirr10sii$dVT/is_pcvirr10sii$dNVT) + 1)
    is_resds$pcv15mek[i] = (k*(is_pcvirr15mek$dVT/is_pcvirr15mek$dNVT) + (j*(1-k)*(is_pcvirr15mek$cVT/is_pcvirr15mek$cNVT) + 1)) / ((is_pcvirr15mek$dVT/is_pcvirr15mek$dNVT) + 1)
    is_resds$pcv20pfz[i] = (k*(is_pcvirr20pfz$dVT/is_pcvirr20pfz$dNVT) + (j*(1-k)*(is_pcvirr20pfz$cVT/is_pcvirr20pfz$cNVT) + 1)) / ((is_pcvirr20pfz$dVT/is_pcvirr20pfz$dNVT) + 1)
    i = i+1
  }
}
is_resds <- is_resds %>% mutate(country = "Israel")

#====================================================================

#split datasets by their PCV
sa_pcvirr10gsk <- sa_pcvsamples %>% filter(pcv == "pcv10gsk")
sa_pcvirr10sii <- sa_pcvsamples %>% filter(pcv == "pcv10sii")
sa_pcvirr15mek <- sa_pcvsamples %>% filter(pcv == "pcv15mek")
sa_pcvirr20pfz <- sa_pcvsamples %>% filter(pcv == "pcv20pfz")

#testing conditions that VT are still in circulation in presence of replacement
i = 1
stR = seq(0, 1, by = 0.01)
rVT = seq(0, 1, by = 0.01)
sa_resds <- data_frame(stR = rep(NA, 10201), rVT = rep(NA, 10201), pcv10gsk = rep(NA, 10201), pcv10sii = rep(NA, 10201), pcv15mek = rep(NA, 10201), pcv20pfz = rep(NA, 10201))

for (j in stR) {
  
  for (k in rVT) {
    sa_resds$stR[i] = j
    sa_resds$rVT[i] = k
    sa_resds$pcv10gsk[i] = (k*(sa_pcvirr10gsk$dVT/sa_pcvirr10gsk$dNVT) + (j*(1-k)*(sa_pcvirr10gsk$cVT/sa_pcvirr10gsk$cNVT) + 1)) / ((sa_pcvirr10gsk$dVT/sa_pcvirr10gsk$dNVT) + 1)
    sa_resds$pcv10sii[i] = (k*(sa_pcvirr10sii$dVT/sa_pcvirr10sii$dNVT) + (j*(1-k)*(sa_pcvirr10sii$cVT/sa_pcvirr10sii$cNVT) + 1)) / ((sa_pcvirr10sii$dVT/sa_pcvirr10sii$dNVT) + 1)
    sa_resds$pcv15mek[i] = (k*(sa_pcvirr15mek$dVT/sa_pcvirr15mek$dNVT) + (j*(1-k)*(sa_pcvirr15mek$cVT/sa_pcvirr15mek$cNVT) + 1)) / ((sa_pcvirr15mek$dVT/sa_pcvirr15mek$dNVT) + 1)
    sa_resds$pcv20pfz[i] = (k*(sa_pcvirr20pfz$dVT/sa_pcvirr20pfz$dNVT) + (j*(1-k)*(sa_pcvirr20pfz$cVT/sa_pcvirr20pfz$cNVT) + 1)) / ((sa_pcvirr20pfz$dVT/sa_pcvirr20pfz$dNVT) + 1)
    i = i+1
  }
}
sa_resds <- sa_resds %>% mutate(country = "South Africa")

#====================================================================

#split datasets by their PCV
mw_pcvirr10gsk <- mw_pcvsamples %>% filter(pcv == "pcv10gsk")
mw_pcvirr10sii <- mw_pcvsamples %>% filter(pcv == "pcv10sii")
mw_pcvirr15mek <- mw_pcvsamples %>% filter(pcv == "pcv15mek")
mw_pcvirr20pfz <- mw_pcvsamples %>% filter(pcv == "pcv20pfz")

#testing conditions that VT are still in circulation in presence of replacement
i = 1
stR = seq(0, 1, by = 0.01)
rVT = seq(0, 1, by = 0.01)
mw_resds <- data_frame(stR = rep(NA, 10201), rVT = rep(NA, 10201), pcv10gsk = rep(NA, 10201), pcv10sii = rep(NA, 10201), pcv15mek = rep(NA, 10201), pcv20pfz = rep(NA, 10201))

for (j in stR) {
  
  for (k in rVT) {
    mw_resds$stR[i] = j
    mw_resds$rVT[i] = k
    mw_resds$pcv10gsk[i] = (k*(mw_pcvirr10gsk$dVT/mw_pcvirr10gsk$dNVT) + (j*(1-k)*(mw_pcvirr10gsk$cVT/mw_pcvirr10gsk$cNVT) + 1)) / ((mw_pcvirr10gsk$dVT/mw_pcvirr10gsk$dNVT) + 1)
    mw_resds$pcv10sii[i] = (k*(mw_pcvirr10sii$dVT/mw_pcvirr10sii$dNVT) + (j*(1-k)*(mw_pcvirr10sii$cVT/mw_pcvirr10sii$cNVT) + 1)) / ((mw_pcvirr10sii$dVT/mw_pcvirr10sii$dNVT) + 1)
    mw_resds$pcv15mek[i] = (k*(mw_pcvirr15mek$dVT/mw_pcvirr15mek$dNVT) + (j*(1-k)*(mw_pcvirr15mek$cVT/mw_pcvirr15mek$cNVT) + 1)) / ((mw_pcvirr15mek$dVT/mw_pcvirr15mek$dNVT) + 1)
    mw_resds$pcv20pfz[i] = (k*(mw_pcvirr20pfz$dVT/mw_pcvirr20pfz$dNVT) + (j*(1-k)*(mw_pcvirr20pfz$cVT/mw_pcvirr20pfz$cNVT) + 1)) / ((mw_pcvirr20pfz$dVT/mw_pcvirr20pfz$dNVT) + 1)
    i = i+1
  }
}
mw_resds <- mw_resds %>% mutate(country = "Malawi")

#====================================================================

#combine all datasets
A <-
  bind_rows(is_resds, mw_resds, sa_resds) %>%
  pivot_longer(c("pcv10gsk", "pcv10sii", "pcv15mek", "pcv20pfz"), names_to = "pcv", values_to = "irr") %>% #combine all sensitivity datasets
  ggplot() +
  geom_tile(aes(x = rVT, y = stR, fill = (1-irr)), linejoin = "bevel") +
  scale_fill_gradientn("z", colours = terrain.colors(100, rev = TRUE, alpha = 0.8)) +
  theme_bw(base_size = 16, base_family = "American typewriter") +
  facet_grid(country ~ pcv) +
  labs(title = "", x = "proportion of residual VT carriage in a mature PCV era", y = "proportion of NVT carriage replacement") +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
  guides(fill = guide_legend(title = "preventable \n    VT IPD")) +
  theme(legend.text = element_text(size = 14), legend.position = "right", legend.title = element_text(size = 14), legend.key.size = unit(1.2,"cm")) +
  theme(strip.text.x = element_text(size = 26), strip.text.y = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

#save combined plots
ggsave(here("output", "fig4_sensResidReplace.png"),
       plot = (A), 
       width = 24, height = 12, unit = "in", dpi = 300)
