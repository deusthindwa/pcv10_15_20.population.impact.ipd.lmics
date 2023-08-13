#Authors: Deus & Dan
#Date: from 01/08/2023
#Title: Potential benefits of newer pneumococcal vaccines on paediatric invasive pneumococcal disease in low- and middle-countries

#====================================================================
  
#simulate some case fatality ratio data
cfr_data <- data.frame(
  id = 1:37,
  st = c("1", "2", "3", "4", "5", "6A/C", "6B", "7F", "8", "9A", "9N", "9V", "10A", "11A", "12A", "12F", "13", "14", "15A", "15B", "15C", "16F", "17F", "18C", "19A", "19F", "20", "22F", "23A", "23F", "24F", "31", "33F", "34", "35B", "35F", "38"),
  invSt = c("1", "3", "4", "5", "6B", "8", "9A", "9N", "10A", "11A", "12F", "13", "14", "15A", "15B", "15C", "18C", "19A", "19F", "20", "23A", "23F", "24F", "31", "33F", "34", rep(NA_character_, 11)),
  inv = sample(c(1:100), size = 37, replace = TRUE, prob = rep(0.5, 100)),
  cfr = sample(c(1:40), size = 37, replace = TRUE, prob = rep(0.5, 40))
)

#fit a linear model


#simulate some IPD data
ipd_data <- data.frame(
  id = 1:640,
  agem = rep(c("0-6m", "7-12m", "13-18m", "19-24m", "25-30m", "31-36m", "37-42m", "43-48m", "49-54m", "55-60m"), 64),
  survyr = sample(c(2015, 2016, 2017, 2018, 2019), size = 640, replace = TRUE, prob = rep(c(0.2), 5)),
  st = sample(c("1", "3", "4", "5", "6A", "6B", "7F", "8", "9V", "10A", "11A", "12F", "14", "15B", "18C", "19A", "19F", "22F", "23F", "33F"), size = 640, replace = TRUE, prob = rep(c(0.05), 20)),
  ipdno = sample(c(1:200), size = 640, replace = TRUE, prob = rep(c(0.01), 200)),
  poprisk = 100000
)
