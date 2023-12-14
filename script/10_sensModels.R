#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================



# #testing conditions that VT are still in circulation
# i = 1
# stR = seq(0, 1, by = 0.01)
# rVT = seq(0, 1, by = 0.01)
# res_DS <- data_frame(stR = rep(NA, 10201), rVT = rep(NA, 10201), irr = rep(NA, 10201))
# 
# for (j in stR) {
#   
#   for (k in rVT) {
#     res_DS$stR[i] = j
#     res_DS$rVT[i] = k
#     res_DS$irr[i] = (k*(y$dVT/y$dNVT) + j*(y$cVT/y$cNVT)+1) / ((y$dVT/y$dNVT)+1)
#     i = i+1
#   }
# }
# 
# C <-
#   res_DS %>%
#   ggplot() +
#   geom_tile(aes(x = rVT, y = stR, fill = (1-res_DS$irr)), linejoin = "bevel") +
#   scale_fill_gradientn("z", colours = terrain.colors(100, rev = TRUE, alpha = 0.8)) +
#   geom_point(aes(x = 0, y = err_DS$sr_min[1]), shape = 4, size = 4, stroke = 2) +
#   theme_bw(base_size = 16, base_family = "American typewriter") +
#   labs(title = "", x = "residual VT carriage in mature PCV era", y = "proportion of NVT replacing VT carriage") +
#   theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
#   guides(fill = guide_legend(title = "proportion of\npreventable IPD")) +
#   theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))
# 
# #save combined plots
# ggsave(here("output", "sfig7_rVTstRsens.png"),
#        plot = (C), 
#        width = 10, height = 10, unit = "in", dpi = 300)





# B <-
#   is_validate %>%
#   ggplot() +
#   geom_density_ridges(aes(x = log(irr), y = sr, fill = sr), size = 0.5, alpha = 0.5) +
#   theme_bw(base_size = 16, base_family = "American typewriter") +
#   labs(title = "", x = "log_incidence rate ratio", y = "Density") + 
#   facet_grid(.~country, scales = "free_y") +
#   theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 0)) +
#   guides(fill = guide_legend(title = "Model validation and scenarios\nof serotype replacement (SR)")) +
#   theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))

# #save combined plots
# ggsave(here("output", "sfig6_obsvspred.png"),
#        plot = (A | B), 
#        width = 16, height = 9, unit = "in", dpi = 300)

# #testing conditions that VT are still in circulation
# i = 1
# stR = seq(0, 1, by = 0.01)
# rVT = seq(0, 1, by = 0.01)
# res_DS <- data_frame(stR = rep(NA, 10201), rVT = rep(NA, 10201), irr = rep(NA, 10201))
# 
# for (j in stR) {
#   
#   for (k in rVT) {
#     res_DS$stR[i] = j
#     res_DS$rVT[i] = k
#     res_DS$irr[i] = (k*(y$dVT/y$dNVT) + j*(y$cVT/y$cNVT)+1) / ((y$dVT/y$dNVT)+1)
#     i = i+1
#   }
# }
# 
# C <-
#   res_DS %>%
#   ggplot() +
#   geom_tile(aes(x = rVT, y = stR, fill = (1-res_DS$irr)), linejoin = "bevel") +
#   scale_fill_gradientn("z", colours = terrain.colors(100, rev = TRUE, alpha = 0.8)) +
#   geom_point(aes(x = 0, y = err_DS$sr_min[1]), shape = 4, size = 4, stroke = 2) +
#   theme_bw(base_size = 16, base_family = "American typewriter") +
#   labs(title = "", x = "residual VT carriage in mature PCV era", y = "proportion of NVT replacing VT carriage") +
#   theme(strip.text.x = element_text(size = 26), strip.background = element_rect(fill = "gray90")) +
#   theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14)) +
#   guides(fill = guide_legend(title = "proportion of\npreventable IPD")) +
#   theme(legend.text = element_text(size = 12), legend.position = "right", legend.title = element_text(size = 12)) +
#   theme(panel.border = element_rect(colour = "black", fill = NA, size = 2))
# 
# #save combined plots
# ggsave(here("output", "sfig7_rVTstRsens.png"),
#        plot = (C), 
#        width = 10, height = 10, unit = "in", dpi = 300)
