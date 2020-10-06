# Prepare Data for Plotting and Analysis
library(tidyverse)
library(here)

dataFiles <- c(
  "20180424-raw-sensitiveCommunity.csv",
  "20190627-raw-doubleResistantCommunity.csv",
  "20190816-raw-singleResistantCommunity.csv",
  "20191124-raw-combinationTreatments.csv"
)

# add Phenotypes
## Phenotype definition
p <- expand.grid(
  agarN = c(FALSE, TRUE),
  agarA = c(FALSE, TRUE),
  agarB = c(FALSE, TRUE),
  agarAB = c(FALSE, TRUE),
  growth = c(FALSE, TRUE),
  growthA = c(FALSE, TRUE),
  growthB = c(FALSE, TRUE),
  growthAB = c(FALSE, TRUE)
)
p$p256 <- 0
for (i in 1:dim(p)[1]) {
  p$p256[i] <- sum(p[i, 1:8] * 2^(0:7))
}
# E: Error (default)
p$p <- "E"
# phenotype U (uninfected): no growth, empty well
eSel <- !p$agarN & !p$agarA & !p$agarB & !p$agarAB & !p$growth & !p$growthA & !p$growthB & !p$growthAB
p[eSel, ]$p <- "U"
# phenotype WT: growth only on agarN and in liquid without treatment
wtSel <- (p$agarN | p$growth) & !p$agarA & !p$agarB & !p$agarAB & !p$growthA & !p$growthB & !p$growthAB
p[wtSel, ]$p <- "S"
# phenotype A: growth in liquid with A or on plate with A, not in other conditions
aSel <- (p$agarA | p$growthA) & (!p$agarB & !p$agarAB & !p$growthB & !p$growthAB)
p[aSel, ]$p <- "A"
# phenotype B: growth in liquid with B or on plate with B, not in other conditions
bSel <- (p$agarB | p$growthB) & (!p$agarA & !p$agarAB & !p$growthA & !p$growthAB)
p[bSel, ]$p <- "B"
# phenotype A/B, a mixed culture of A and B, no growth on AB but growth in A and B
abmixSel <- !(p$agarAB | p$growthAB) & ((p$agarA | p$growthA) & (p$agarB | p$growthB))
p[abmixSel, ]$p <- "A/B"
# phenotype AB, double Resistance, growth on A and B
abSel <- p$agarAB | p$growthAB
p[abSel, ]$p <- "AB"
# set impossible combinations to E (growth in liquid is always only measured in one condition during the experiment)
impossibleSel <- !p$growth & (p$growthA | p$growthB | p$growthAB)
p[impossibleSel, ]$p <- "E"

# save csv for reference
write_csv(p, here("1-phenotypeDefinition.csv"))


## add phenotype to data & save
# set threshold for OD
odThreshold <- 0.1
phenotypes <- list()
frequencies <- list()
meanFrequencies <- list()

for (i in seq_along(dataFiles)) {
  rawData <- read_csv(
    file = here("raw", dataFiles[i]),
    col_types = "iiciiicccdllllc"
  )
  phenotypes[[i]] <- rawData %>%
    mutate(
      growth = OD > odThreshold,
      growthA = growth & treatment == "A",
      growthB = growth & treatment == "B",
      growthAB = growth & (treatment == "AB" | treatment == "AB22")
    ) %>%
    left_join(p,
      by = c(
        "agarA", "agarAB", "agarB", "agarN",
        "growth", "growthA", "growthB", "growthAB"
      )
    ) %>%
    # remove blanks
    filter(type == "exp") %>%
    # remove unneeded variable
    select(-(agarA:growthAB))

  frequencies[[i]] <- phenotypes[[i]] %>%
    group_by(plate, transfer, rep) %>%
    count(p) %>%
    mutate(
      p = factor(p, levels = c("U", "S", "A", "B", "A/B", "AB", "E")),
      f = n / sum(n)
    ) %>%
    # complete the table with all possible phenotypes
    complete(
      plate, transfer, rep, p,
      fill = list(n = 0, f = 0)
    )

  # correct NA phenotypes (complete() added 0 for other phenotypes)
  naPhenotypes <- filter(frequencies[[i]], is.na(p))
  for (j in seq_along(naPhenotypes[[1]])) {
    frequencies[[i]][
      frequencies[[i]]$plate == naPhenotypes$plate[j] &
        frequencies[[i]]$transfer == naPhenotypes$transfer[j] &
        frequencies[[i]]$rep == naPhenotypes$rep[j], 5:6
    ] <- NA
  }
  # remove NA 'phenotype'
  frequencies[[i]] <- filter(frequencies[[i]], !is.na(p))

  meanFrequencies[[i]] <- frequencies[[i]] %>%
    group_by(plate, transfer, p) %>%
    summarise(
      mean = mean(f),
      n = n(),
      se = sd(f) / sqrt(n),
      max = max(f),
      min = min(f)
    )
  # save in file
  write_csv(phenotypes[[i]], path = here(str_replace(dataFiles[i], "-raw-", "-")))
  write_csv(frequencies[[i]], path = here(str_replace(dataFiles[i], "-raw-", "-F-")))
  write_csv(meanFrequencies[[i]], path = here(str_replace(dataFiles[i], "-raw-", "-meanF-")))
}
