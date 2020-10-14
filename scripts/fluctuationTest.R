library(tidyverse)
library(flan)
library(here)
library(scales)

source(here("scripts/plotParameter.R"))

# output directory for plots & tables
outDir <- "/Users/daniel/polybox/work/projects/TreatmentStrategies/paper"

countDataRaw <- read_csv(
  file = here("20200820-fluctuationTest.csv"),
  col_types = cols(
    n = col_double(),
    strain = col_character(),
    drug = col_character(),
    mut = col_double(),
    tot = col_double()
  )
)

countData <- countDataRaw %>%
  group_by(strain, drug) %>%
  group_split()

flanTest <- list()
estimate <- list()

for (i in seq_along(countData)) {
  flanTest[[i]] <- flan.test(
    mc = countData[[i]]$mut,
    fn = countData[[i]]$tot
  )

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
      str_c(strain, " -> ", drug, "R"),
      levels = c("mutS -> NalR", "mutSSmR -> NalR", "mutS -> SmR", "mutSNalR -> SmR")
    )
  ) 

estimates

strainLabels <- c(
  "mutS" = "strain S\n(ΔmutS)",
  "mutSNalR" = "strain A\n(ΔmutS NalR)",
  "mutSSmR" = "strain B\n(ΔmutS SmR)")

ftPlot <- ggplot(
    data = estimates,
    mapping = aes(drug, mutProb)) +
  facet_grid(cols = vars(strain),
    scales = "free_x", space = "free_x", 
    switch = "x", labeller = labeller(strain = strainLabels)) +
  geom_pointrange(aes(ymin = mutProb95low, ymax = mutProb95high)) +
  scale_y_log10(
    "mutation probability (95% conf. interval)",
    labels = label_scientific()) +
  scale_x_discrete(NULL, labels = c("Nal" = "Nal 40μg/ml", "Sm" = "Sm 100μg/ml")) +
  plotTheme +
  theme(
    # facet labels
     strip.background = element_rect(fill = NA, color = NA),
     strip.text.x = element_text(hjust = 0.5),
     strip.placement = "outside",
    # axis.line.x  = element_line(color = "black", size = 0.25)
  )
ftPlot


plotWidth <- 4 * 25
plotHeight <- 90

ggsave(
  filename = file.path(outDir, "SI", "FigS2_partB.pdf"),
  device = cairo_pdf,
  # device = "tiff",
  # compression = "lzw", type = "cairo",
  # dpi = 600,
  plot = ftPlot,
  width = plotWidth, height = plotHeight, units = "mm"
)
