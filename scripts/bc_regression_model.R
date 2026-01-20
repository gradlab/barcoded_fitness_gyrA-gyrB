library(brms)
library(dplyr)
library(tidyverse)

data <- read.csv('model_input.csv')

data$log.prop <- log(data$prop_to_ref)

fit2 <- brm(log.prop ~ -1 + gyra:time + (-1 + gyrb:time|gyra) + unique_id + (-1 + time|barcode), data=data, cores = 4, family=gaussian())


#add gyra and gyrb effects

gyras <- c('S.D','S.A','S.G','S.N','F.D','F.A','F.G','F.N')
col1start = 'b_gyra'
col1end = ':time'
col2start = 'r_gyra['
col2end = ',gyrb:time]'
modeloutput = as.data.frame(fit2)
combined <- data.frame(matrix(ncol = 0, nrow = nrow(modeloutput)))

for (gyra in gyras) {
  col1 = paste(col1start, gyra, col1end, sep="")
  col2 = paste(col2start, gyra, col2end, sep="")
  newname = paste(gyra, 'gyrb', sep='_')
  newname2 = paste('gyrb_given',gyra, sep='_')
  combined[[newname]] = modeloutput[[col1]] + modeloutput[[col2]]
  combined[[newname2]] = modeloutput[[col2]]
  combined[[gyra]] = modeloutput[[col1]]
  #if (gyra != 'S.D') {
    #combined[[gyra]] = modeloutput[[col1]]
  #}
}

#get median and ci
summary <- combined %>%
  pivot_longer(everything(), 
               names_to = "variable", 
               values_to = "value") %>%
  group_by(variable) %>%
  summarise(median = median(value),
            lower_CI = quantile(value, 0.025),
            upper_CI = quantile(value, 0.975))

summary <- summary %>% mutate(has_gyrb = str_detect(variable, "_gyrb"),
                   gyra = str_extract(variable, "[A-Z]\\.[A-Z]")) %>%
                   filter(!str_detect(variable, "gyrb_.*"))


write.csv(summary, file = "gyra+gyrb_effects.csv", row.names = FALSE)