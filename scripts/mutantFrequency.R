# mutant frequency figures & stats

# packages
library(here)
library(tidyverse)
library(multcomp)
library(multcompView)
library(knitr)
library(kableExtra)
# common plot parameters
source(here("scripts/plotParameter.R"))

# output directory for plots & tables
outDir <- "/Users/daniel/Documents/work/projects/TreatmentStrategies/Paper"

# data
dataFile <- "20180424-sensitiveCommunity.csv"

phenotypeData <- read_csv(
    file = here(dataFile),
    col_types = cols(
      plate = col_double(),
      well = col_double(),
      row = col_character(),
      col = col_double(),
      rep = col_double(),
      transfer = col_double(),
      turnoverStrain = col_character(),
      infectionToWell = col_integer(),
      treatment = col_character(),
      OD = col_double(),
      p256 = col_double(),
      p = col_character()
    ))

phenotypeData %>%
  filter(turnoverStrain == "WT") %>%
  dplyr::select(turnoverStrain, treatment, OD) %>%
  mutate(growth = OD > 0.1) %>%
  group_by(treatment) %>%
  dplyr::summarize(
    nRes = sum(growth),
    nTot = n(),
    f = nRes / nTot
  )


plotData <- phenotypeData %>%
  filter(turnoverStrain == "WT") %>%
  dplyr::select(turnoverStrain, treatment, OD)

ggplot(plotData, aes(factor(treatment), OD)) +
  geom_jitter() +
  geom_violin()
