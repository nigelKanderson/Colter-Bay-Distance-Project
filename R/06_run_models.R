library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(performance)

#data <- readRDS("data/data_out.rds")

#glimpse(data)

run_models <- function(data) {
  library(tidyverse)
  library(glmmTMB)
  library(DHARMa)
  library(performance)
  
  model_data <- data %>%
    filter(
      !is.na(detections),
      !is.na(jd),
      !is.na(pct_forest)
    )
  
  m0 <- glmmTMB(
    detections ~ 1 + 
      (1|site) +
      (1|year),
    data = model_data,
    family = nbinom2()
  )
  
  m1 <- glmmTMB(
    detections ~
      jd +
      I(jd^2) +
      mean_moonlight +
      max_moonlight +
      mean_phase +
      (1|site) +
      (1|year),
    data = model_data,
    family = nbinom2()
  )
  
  m2 <- glmmTMB(
    detections ~
      jd +
      I(jd^2) +
      mean_moonlight +
      max_moonlight +
      mean_phase +
      pct_forest +
      color * intensity +
      (1|site) +
      (1|year),
    data = model_data,
    family = nbinom2()
  )
  
  simple_model1 <- glmmTMB(
    detections ~
      color * intensity +
      mean_phase +
      pct_nonforest +
      I(jd^2) +
      (1|site),
    data = model_data,
    ziformula = ~1,
    family = nbinom2()
  )
  
  list(
    m0 = m0,
    m1 = m1,
    m2 = m2,
    sm = simple_model1
  )
  
}




