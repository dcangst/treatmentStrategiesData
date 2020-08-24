library(tidyverse)
library(flan)
library(here)

# output directory for plots & tables
outDir <- "/Users/daniel/Documents/work/projects/TreatmentStrategies/Paper"

countDataRaw <- read_csv(
  file = here("20200820-fluctuationTest.csv"),
  col_types = cols(
    n = col_double(),
    strain = col_character(),
    drug = col_character(),
    mut = col_double(),
    tot = col_double()))

countData <- countDataRaw %>%
  group_by(strain, drug) %>%
  group_split()

flanTest <- list()
estimate <- list()

for (i in seq_along(countData)){
  flanTest[[i]] <- flan.test(
    mc = countData[[i]]$mut,
    fn = countData[[i]]$tot)
  
  estimate[[i]] <- tibble(
    strain = countData[[i]]$strain[1],
    drug = countData[[i]]$drug[1],
    mutProb = flanTest[[i]]$estimate[1],
    mutProb95low = flanTest[[i]]$conf.int[1, 1],
    mutProb95high = flanTest[[i]]$conf.int[2, 1],
    fitness = flanTest[[i]]$estimate[2],
    fitness95low = flanTest[[i]]$conf.int[1, 2],
    fitness95high = flanTest[[i]]$conf.int[2, 2] 
  )
}

estimates <- bind_rows(estimate) %>%
  mutate(
    strain_drug = factor(
      str_c(strain, "_", drug),
      levels = c("mutS_Nal", "mutSSmR_Nal", "mutS_Sm", "mutSNalR_Sm")))

estimates

ggplot(
  data = estimates,
  mapping = aes(strain_drug, mutProb)) +
  geom_pointrange(aes(ymin = mutProb95low, ymax = mutProb95high)) +
  scale_y_log10(name = "mutation probability (95% conf. interval)")
