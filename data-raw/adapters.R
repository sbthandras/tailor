rm(list = ls())

library(tailor)
data(rbps)

adapters <- find_all_adapters(
  ids = rbps$Core_ORF,
  data = rbps,
  method = "plateau",
  type = "xbar.one"
)

save(adapters, file = "data/adapters.rda")
