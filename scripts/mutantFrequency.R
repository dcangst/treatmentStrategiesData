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
outDir <- "/Users/daniel/polybox/work/projects/TreatmentStrategies/paper"

# data
dataFile <- "20180424-sensitiveCommunity.csv"
phenotypeDefFile <- "1-phenotypeDefinition.csv"

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
  )
)

phenotypeDefinition <- read_csv(phenotypeDefFile)

filter(phenotypeDefinition, p == "A" & agarA == "TRUE")$p256

growthAgar <- function(treatment, p256) {
  out <- list()
  for (i in seq_along(p256)) {
    if (treatment[[i]] == "none") {
      out[[i]] <- p256[[i]] != 0
    } else if (treatment[[i]] == "A") {
      out[[i]] <- p256[[i]] %in% filter(phenotypeDefinition, p %in% c("A", "A/B", "AB") & agarA == "TRUE")$p256
    } else if (treatment[[i]] == "B") {
      out[[i]] <- p256[[i]] %in% filter(phenotypeDefinition, p %in% c("B", "A/B", "AB") & agarB == "TRUE")$p256
    } else if (treatment[[i]] == "AB") {
      out[[i]] <- p256[[i]] %in% filter(phenotypeDefinition, p %in% c("AB") & agarAB == "TRUE")$p256
    }
  }

  return(unlist(out))
}

odPlotData <- phenotypeData %>%
  filter(turnoverStrain == "WT") %>%
  dplyr::select(turnoverStrain, treatment, OD)

odPlot <- ggplot(odPlotData, aes(factor(treatment), OD)) +
  geom_jitter() +
  geom_violin()

# odPlot

growthAgarTable <- expand.grid(treatment = unique(phenotypeData$treatment), p256 = 1:256) %>%
  as_tibble() %>%
  mutate(
    growthAgar = growthAgar(treatment, p256)
  )

mutFrequency <- phenotypeData %>%
  filter(turnoverStrain == "WT") %>%
  dplyr::select(turnoverStrain, treatment, OD, p256) %>%
  left_join(growthAgarTable, by = c("treatment", "p256")) %>%
  mutate(
    growth = OD > 0.1
  ) %>%
  group_by(treatment) %>%
  dplyr::summarize(
    nRes = sum(growth),
    nTot = n(),
    f = nRes / nTot,
    nResAgar = sum(growthAgar, na.rm = TRUE),
    nTotAgar = n(),
    fAgar = nResAgar / nTotAgar,
    .groups = "keep"
  )

labelNudge <- 0.01

mutFrequencyPlotData <- mutFrequency %>%
  filter(treatment != "none") %>%
  mutate(
    fLabel = if_else(f + labelNudge > 0, f + labelNudge, f + labelNudge)
  )

mutFrequencyPlot <- ggplot(
  data = mutFrequencyPlotData,
  mapping = aes(
    x = treatment, y = f,
    label = str_c(nRes, "/", nTot)
  )
) +
  geom_col() +
  geom_text(mapping = aes(y = fLabel), vjust = 0, size = 1.5) +
  scale_y_continuous("mutant frequency", limits = c(0, 1)) +
  scale_x_discrete("") +
  plotTheme +
  theme(axis.title.x = element_blank())

plotWidth <- 50
plotHeight <- 90

ggsave(
  filename = file.path(outDir, "figures", "temp", "mutantFrequency.pdf"),
  # device = "tiff",
  # compression = "lzw", type = "cairo",
  # dpi = 600,
  plot = mutFrequencyPlot,
  width = plotWidth, height = plotHeight, units = "mm"
)
