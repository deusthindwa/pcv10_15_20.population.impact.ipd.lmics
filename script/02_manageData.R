#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
  
#simulate some case fatality ratio data
data_cfr <- data.frame(
  id = 1:37, #id of each serotype
  st = c("1", "2", "3", "4", "5", "6A/C", "6B", "7F", "8", "9A", "9N", "9V", "10A", "11A", "12A", "12F", "13", "14", "15A", "15B", "15C", "16F", "17F", "18C", "19A", "19F", "20", "22F", "23A", "23F", "24F", "31", "33F", "34", "35B", "35F", "38"),
  invSt = c("1", NA, "3", "4", "5", "6B", "8", NA, "9A", "9N", NA, "10A", "11A", "12F", "13", NA, NA, "14", "15A", "15B", NA, "15C", "18C",NA, "19A", "19F", "20", "23A", "23F", NA, NA, "24F", "31", NA, "33F", "34", NA), #serotypes with invasiveness data
  inv = sample(c(20:100), size = 37, replace = TRUE, prob = rep(0.5, 81)),
  cfr = sample(c(10:40), size = 37, replace = TRUE, prob = rep(0.5, 31))) %>%
  dplyr::mutate(inv = ifelse(is.na(invSt), NA_integer_, inv))

#fit a linear model of invasiveness and case fatality ratio
data_cfr <- 
  data_cfr %>%  
  mutate(log_inv = log(inv),
         log_cfr = log(cfr))

model_cfr <- tidy(lm(log_inv ~ log_cfr, data = data_cfr), 
               exponentiate = TRUE, 
               conf.int = TRUE, 
               conf.level = 0.95)

data_cfr <- data_cfr %>% dplyr::mutate(fit_inv = model_cfr$estimate[1] + model_cfr$estimate[2]*log_cfr)

#visualize model fit
data_cfr %>%
ggplot() +
  geom_point(aes(x = cfr, y = inv, color = st)) +
  geom_line(aes(x = cfr, y = fit_inv)) +
  theme_bw() + 
  theme(legend.position = "right")

#====================================================================

#simulate some data on age/serotype-specific carriage and IPD
data_prev <- data.frame(
  id = 1:2368,
  agem = rep(c("0-6m", "7-12m", "13-18m", "19-24m", "25-30m", "31-36m", "37-42m", "43-48m", "49-54m", "55-60m"), 2368)) %>%
  dplyr::group_by(agem) %>%
  dplyr::mutate(carrSt = sample(c("1", "2", "3", "4", "5", "6A/C", "6B", "7F", "8", "9A", "9N", "9V", "10A", "11A", "12A", "12F", 
                                  "13", "14", "15A", "15B", "15C", "16F", "17F", "18C", "19A", "19F", "20", "22F", "23A", "23F", 
                                  "24F", "31", "33F", "34", "35B", "35F", "38"), size = 2368, replace = TRUE, prob = rep(c(0.5), 37)),
                carrPos = sample(c(100:2000), size = 37, replace = TRUE, prob = rep(0.5, 81)),
                carrSamp = 2368) %>%
  dplyr::ungroup()





















