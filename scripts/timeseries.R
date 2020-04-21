# Timeseries plot for experiment with sensitive community
# and timeseries plots for all the experiments with varying community for SI

# Packages
library(here)
library(tidyverse)
# common plot parameters
source(here("figures/scripts/plotParameter.R"))

# data
dataFiles <- c(
  "20180424-meanF-sensitiveCommunity.csv",
  "20190627-meanF-doubleResistantCommunity.csv",
  "20190816-meanF-singleResistantCommunity.csv")

meanFrequencies <- list()
for (i in dataFiles) {
  meanFrequencies[[i]] <- read_csv(
    file = here("data", i),
    col_types = "iicdiddd") %>%
    # drop phenotype E if it doesn't show up
    filter(!(p == "E" & mean == 0)) %>%
    # order of phenotypes (for legend presentation)
    mutate(
      p = fct_drop(factor(p, levels = names(strainColors)))
    )
}

simDataFile <- "SIM-F-sensitiveCommunity.csv"
simFrequencies <- read_csv(
  file = here("data", "simulation", simDataFile),
  col_types = "iicd"
  ) %>%
  # add A/B for plotting purposes
  bind_rows(
    tibble(plate = 1, transfer = 0, p = "A/B", f = NA)
  ) %>%
  mutate(
    p = fct_drop(factor(p, levels = names(strainColors)))
  )


# plot

plotAllPlatesRibbon <- ggplot(
  data = filter(meanFrequencies[[1]], transfer <= 40) %>%
    mutate(p = fct_relevel(p, levels = c("U", "A", "A/B"))),
  mapping = aes(x = transfer, y = mean, color = p)) +
  facet_wrap(vars(plate), ncol = 2, nrow = 3, labeller = labeller(plate = plateLabels)) +
  geom_ribbon(aes(x = transfer, ymin = min, ymax = max, fill = p), alpha = 0.25, inherit.aes = FALSE)

plotAllPlatesSensitiveSim <- plotAllPlatesRibbon +
  # add simulation data
  geom_line(data = filter(simFrequencies, transfer <= 40) %>%
    # change factor for legend sorting
    mutate(p = fct_relevel(p, levels = c("U", "A", "A/B"))),
    mapping = aes(y = f),
    size = 0.3) +
  scale_colour_manual(
    name = "Phenotype",
    values = strainColors,
    guide = guide_legend(nrow = 3)) +
  scale_fill_manual(
    name = "Phenotype",
    values = strainColors)  +
  scale_y_continuous(
      name = "population phenotype frequency",
      limits = c(0, 1), breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
  plotTheme +
  theme(
    legend.position = c(0.75, 0.45),
    legend.box.background = element_rect(size = 0.3),
    legend.key.size = unit(7, "pt"),
    legend.key.width = unit(7, "pt"),
    legend.key.height = unit(7, "pt"),
    legend.title = element_blank(),
    legend.margin = margin(0, 4, 4, 4),
    legend.justification = c(0, 0),
    legend.direction = "vertical")

plotWidth <- 110
plotHeight <- 100
ggsave(
  filename = here("figures", "fig1partA-timeseries-sensitiveCommunity.pdf"),
  plot = plotAllPlatesSensitiveSim,
  dpi = 600,
  width = plotWidth, height = plotHeight, units = "mm")

# SI plots
# complete timeseries sensitiveCommunity
# make bigger than main plot
plotWidth <- 110
plotHeight <- 150

plotAllPlatesSI <- plotAllPlatesRibbon +
  # mean line
  geom_line(size = 0.3) +
  scale_colour_manual(
    name = "Phenotype",
    values = strainColors,
    guide = guide_legend(nrow = 1)) +
  scale_fill_manual(
    name = "Phenotype",
    values = strainColors,
    guide = "none")  +
  scale_y_continuous(
      name = "mean frequency and range (n = 4)",
      limits = c(0, 1), breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
  plotTheme +
  theme(legend.position = "bottom")

# complete timeseries sensitiveCommunity
plateLabelsT53 <- plateLabels
plateLabelsT53[2] <- "b) mono A / B (40+)"
plateLabelsT53[3] <- "c) mono B / A (40+)"
plotAllPlatesSISensitive <- plotAllPlatesSI +
    facet_wrap(vars(plate), ncol = 2, nrow = 3, labeller = labeller(plate = plateLabelsT53))
ggsave(
  filename = here("SI", "S3_Fig.pdf"),
  plot = plotAllPlatesSISensitive %+% meanFrequencies[[1]],
  width = plotWidth, height = plotHeight, units = "mm")

# complete timeseries singleResistantCommunity
ggsave(
  filename = here("SI", "S4_Fig.pdf"),
  plot = plotAllPlatesSI %+% meanFrequencies[[3]],
  width = plotWidth, height = plotHeight, units = "mm")

# complete timeseries doubleResistantCommunity
ggsave(
  filename = here("SI", "S5_Fig.pdf"),
  plot = plotAllPlatesSI %+% meanFrequencies[[2]],
  width = plotWidth, height = plotHeight, units = "mm")
