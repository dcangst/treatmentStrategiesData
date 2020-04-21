# phenotype frequency figures & stats

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
dataFiles <- c(
  "20180424-F-sensitiveCommunity.csv",
  "20190816-F-singleResistantCommunity.csv",
  "20190627-F-doubleResistantCommunity.csv")

experimentNames <- c(
  "Scenario 0",
  "Scenario I",
  "Scenario II")

experimentNamesPrint <- c(
  "textstyle(Scenario)~symbol(\"\\306\")",
  "textstyle(Scenario)~I",
  "textstyle(Scenario)~II")

frequenciesList <- list()
for (i in seq_along(dataFiles)) {
  frequenciesList[[i]] <- read_csv(
    file = here(dataFiles[i]),
    col_types = "iiicid") %>%
    filter(p != "E") %>%
    mutate(
      community = experimentNames[i],
      communityPrint = experimentNamesPrint[i]
    )
}

# Summarize data
meanFrequencies <- bind_rows(frequenciesList) %>%
mutate(
  p = factor(p, levels = names(strainColors)[1:6]),
  community = factor(community, levels = experimentNames),
  communityPrint = factor(communityPrint, levels = experimentNamesPrint),
  ) %>%
filter(transfer %in% 9:12, p != "E") %>%
group_by(community, communityPrint, plate, rep, p) %>%
summarise(
  meanF = mean(f),
  n = n(),
  se = sd(f) / sqrt(n),
  max = max(f),
  min = min(f))

frequenciesAllRes <- bind_rows(frequenciesList) %>%
  mutate(
    p = factor(p, levels = names(strainColors)[1:6]),
    community = factor(community, levels = experimentNames),
    communityPrint = factor(communityPrint, levels = experimentNamesPrint),
    ) %>%
  filter(transfer %in% 9:12, p != "E") %>%
  # combine 'resistant phenotypes' i.e. all resistant to at least one of the used drugs
  mutate(
    pRes = case_when(
      p %in% c("A", "B", "AB", "A/B") ~ "R",
      p == "U" ~ "U",
      p == "S" ~ "S")
  ) %>%
  group_by(community, communityPrint, plate, transfer, rep, pRes) %>%
  summarize(
    f2 = sum(f)) %>%
  mutate(
    pRes = factor(pRes, levels = c("U", "S", "R"))
  )
  
meanFrequenciesAllRes <- frequenciesAllRes %>%
  group_by(community, communityPrint, plate, rep, pRes) %>%
  summarise(
    meanF = mean(f2),
    n = n(),
    se = sd(f2) / sqrt(n),
    max = max(f2),
    min = min(f2))

# ANOVA & generalized linear hypothesis test

## functions
modelFunction <- function(df) {
  # independent var: phenotype & phenotype/plate interaction
  lm(meanF ~ p + p * plate, data = df)
}

compLetters <- function(test) {
  # calculate "letter-based representation of pairwise comparisons"
  out <- tibble(
      comp = names(test$tstat),
      pValue = test$pvalues) %>%
    separate(comp, into = c("p", "comparison"), sep = ":") %>%
    mutate(
      comparison = stringi::stri_reverse(str_replace_all(comparison," ",""))
    ) %>%
    group_by(p) %>%
    nest() %>%
    mutate(
      compVector = map(data, function(x) {
        out <- x$pValue
        names(out) <- x$comparison
        return(out)}),
      compLettersList = map(compVector, multcompLetters),
      compLetters = map(compLettersList, ~enframe(.x$Letters, name = c("plate"), value = "compLetter"))
    ) %>%
    dplyr::select(p, compLetters) %>%
    unnest(cols = c(compLetters))
  return(out)
}

meanFrequenciesRSU <- meanFrequenciesAllRes %>%
  ungroup() %>%
  group_by(community) %>%
  mutate(
    p = fct_drop(pRes),
    plate = as_factor(plate))

## prepare matrices for multiple comparisons
tmp <- expand.grid(plate = unique(meanFrequenciesRSU$plate),
                   p = unique(meanFrequenciesRSU$p))
modelMatrix <- model.matrix(~ p * plate, data = tmp)

tukey <- contrMat(table(meanFrequenciesRSU$plate), "Tukey")
tukeyMatrix1 <- cbind(
  tukey,
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)),
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)))
rownames(tukeyMatrix1) <- paste(levels(meanFrequenciesRSU$p)[1], rownames(tukeyMatrix1), sep = ":")

tukeyMatrix2 <- cbind(
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)),
  tukey,
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)))
rownames(tukeyMatrix2) <- paste(levels(meanFrequenciesRSU$p)[2], rownames(tukeyMatrix2), sep = ":")

tukeyMatrix3 <- cbind(
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)),
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)),
  tukey)
rownames(tukeyMatrix3) <- paste(levels(meanFrequenciesRSU$p)[3], rownames(tukeyMatrix3), sep = ":")

tukeyMatrix <- rbind(tukeyMatrix1, tukeyMatrix2, tukeyMatrix3)
colnames(tukeyMatrix) <- rep(colnames(tukey), 3)

# compute for all communities (takes a while)
testRSU <- meanFrequenciesRSU %>%
  nest() %>%
  mutate(
    model = map(data, modelFunction),
    anova = map(model, anova),
    mcomp = map(model, ~glht(., linfct = tukeyMatrix %*% modelMatrix)),
    mcompSum = map(mcomp, summary),
    plotLetters = map(mcompSum, ~compLetters(.$test)))

names(testRSU$anova) <- levels(meanFrequenciesRSU$community)
names(testRSU$mcompSum) <- levels(meanFrequenciesRSU$community)

plotLettersRAS <- testRSU %>%
  dplyr::select(community, plotLetters) %>%
  unnest(cols = c(plotLetters)) %>%
  mutate(
    plate = as.integer(plate),
    p = factor(p, levels = c("U", "S", "R"))
  )

# save statistical tables
options(knitr.kable.NA = '')
for (i in seq_along(testRSU$mcompSum)) {
  aovTable <- data.frame(testRSU$anova[[i]]) %>%
    rename("Sum Sq" = Sum.Sq, "Mean Sq" = Mean.Sq, "F" = F.value, "$\\Pr(>F)$" = Pr..F.) %>%
    mutate("$\\Pr(>F)$" = if_else(
      "$\\Pr(>F)$" < 0.001,
      "$< 0.001$", as.character("$\\Pr(>F)$")))
  row.names(aovTable) <- row.names(testRSU$anova[[i]])
  texTableAov <- kable(aovTable,
    format = "latex",
    caption = str_c(
      "\\textbf{", levels(meanFrequenciesRSU$community)[i], ": ",
      "Effect of treatment strategy on the frequency of uninfected and resistant populations (ANOVA).}"),
    longtable = TRUE, booktabs = TRUE,
    label = str_c("SI-StatPhenotypeAnovaRSU", i),
    format.args = list(digits = 3, drop0trailing = TRUE),
    escape = FALSE
  )

  statPrintRSU <- tibble(
      p = names(testRSU$mcompSum[[i]]$test$coefficients),
      Estimate = signif(testRSU$mcompSum[[i]]$test$coefficients, 3),
      "Std. Error" = signif(testRSU$mcompSum[[i]]$test$sigma, 3),
      t = signif(testRSU$mcompSum[[i]]$test$tstat, 3),
      "$\\Pr(>|t|)$" = if_else(
        testRSU$mcompSum[[i]]$test$pvalue < 0.001,
        "$< 0.001$", as.character(signif(testRSU$mcompSum[[i]]$test$pvalue, 3))
        )
      ) %>%
    separate(p, into = c("Phenotype", "linfct"), sep = ":") %>%
    mutate(
      "Linear Hypothesis" = str_replace(str_c(linfct, " == 0"), "-", "--")) %>%
    dplyr::select(Phenotype, "Linear Hypothesis", Estimate, "Std. Error", "$\\Pr(>|t|)$")

  texTableMcomp <- kable(statPrintRSU, format = "latex",
    caption = str_c(
      "\\textbf{", levels(meanFrequenciesRSU$community)[i], ": ",
      "Multiple comparison of phenotype frequencies between treatment strategies.}"),
    longtable = TRUE, booktabs = TRUE,
    label = str_c("SI-StatPhenotypeCompRSU-", i),
    format.args = list(drop0trailing = TRUE),
    align = "lrrrr",
    escape = FALSE
  ) %>%
  kable_styling(latex_options = c("repeat_header"))
  
  write(str_c(texTableAov, "\n\n", texTableMcomp),
    file = file.path(outDir, "SI", "tables", str_c("tableSI-phenotypeFrequencyComparisonRSU-raw-",i,".tex"))) 
}

# Figures

resColors <- c(
  "U" = "#b3b3b3",
  "S" = "#0cb7eb",
  "R" = "#eb400c")

dataRUletters <- meanFrequenciesAllRes %>%
  group_by(community, communityPrint, plate, pRes) %>%
  summarize(
    meanF = max(meanF) + 0.1
  ) %>%
  left_join(plotLettersRAS, by = c("community", "plate", "pRes" = "p")) %>%
  filter(pRes %in% c("U", "R"))

freqPlotAllRes <- ggplot(
  data = filter(meanFrequenciesAllRes %>% filter(pRes %in% c("U", "R"))),
  aes(x = as.character(plate), y = meanF, colour = pRes)) +
  facet_wrap(vars(communityPrint), ncol = 1, scales = "free", labeller = label_parsed) +
  stat_summary(
    mapping = aes(x = as.character(plate), y = meanF, group = pRes, colour = pRes), inherit.aes = FALSE,
    fun.y = median, fun.ymin = median, fun.ymax = median,
    geom = "errorbar", width = 2 / 3, size = 0.5,
    position = position_dodge(width = 0.75), show.legend = TRUE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
    size = 1.5, alpha = 0.5, shape = 20, stroke = 0, show.legend = FALSE) +
  geom_text(
    data = dataRUletters,
    mapping = aes(label = compLetter),
    position = position_dodge(width = 0.75),
    vjust = 0.5, hjust = 0.5,
    size = 5 / .pt,
    alpha = 0.5,
    show.legend = FALSE
    ) +
  scale_y_continuous(
    name = "population phenotype frequency",
    limits = c(0, 1.1), breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
  scale_x_discrete(NULL,
    limits = 1:6,
    labels = str_sub(plateLabels, 4)
  ) +
  scale_colour_manual(
      name = "Phenotype",
      values = resColors) +
  geom_vline(xintercept = seq(1.5, 5.5, 1), size = 0.3, colour = "grey92") +
  plotTheme +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = c(0.75, 0.85),
    legend.box.background = element_rect(size = 0.3),
    legend.key.size = unit(7, "pt"),
    legend.key.width = unit(7, "pt"),
    legend.key.height = unit(7, "pt"),
    legend.title = element_blank(),
    legend.margin = margin(2, 2, 2, 2),
    legend.justification = c(0, 0),
    legend.direction = "horizontal",
    strip.text.x = element_text(size = 8, hjust = 0, face = "bold"))

plotWidth <- 100
plotHeight <- 90

ggsave(
  filename = file.path(outDir, "figures", "fig3.tiff"),
  device = "tiff",
  compression = "lzw", type = "cairo",
  dpi = 600,
  plot = freqPlotAllRes,
  width = plotWidth, height = plotHeight, units = "mm")

# compare high-level and low-level resistance

## data

dataFile <- "20191124-raw-combinationTreatments.csv"

rawData <- read_csv(
  file = here("raw", dataFile),
  col_types = "iiciiicicdllllc")

odThreshold <- 0.1
# phenotypes
p <- read_csv(
  file = here("1-phenotypeDefinition.csv"),
  col_types = "llllllllic")
# treatment Info
comboTreatments <- tribble(
  ~plate, ~cA, ~cB,
      1,   0,   0,
      2,   2,   2,
      3,   1,   2,
      4,   2,   1,
      5,   1,   1,
      6, 0.5,   1,
      7,   1, 0.5,
      8, 0.5, 0.5,
  )


# calculate phenotypes not sensitive to high levels of drugs (i.e. treatment with 2x MIC and/or growth on plates)
  phenotypesHL <- rawData %>%
    mutate(
      growth = OD > odThreshold,
      growthA = growth & treatment == "A",
      growthB = growth & treatment == "B",
      growthAB = growth & (treatment == "AB" | treatment == "AB22")) %>%
    left_join(p,
      by = c("agarA", "agarAB", "agarB", "agarN",
            "growth", "growthA", "growthB", "growthAB")) %>%
    # remove blanks
    filter(type == "exp") %>%
    # remove unneeded variable
    dplyr::select(- (agarA:growthAB)) %>%
    mutate(
      pType = "HL"
    )

# calculate phenotypes not sensitive to current level of drugs (i.e. not sensitive to the current treatment)
  phenotypesLL <- rawData %>%
    mutate(
      growth = OD > odThreshold,
      growthA = growth & treatment == "A",
      growthB = growth & treatment == "B",
      growthAB = growth & (treatment %in% c("AB22", "AB12", "AB21", "AB11", "AB051", "AB105", "AB0505"))) %>%
    left_join(p,
      by = c("agarA", "agarAB", "agarB", "agarN",
            "growth", "growthA", "growthB", "growthAB")) %>%
    # remove blanks
    filter(type == "exp") %>%
    # remove unneeded variable
    dplyr::select(- (agarA:growthAB)) %>%
    mutate(
      pType = "LL"
    )

# combine & summarize
  phenotypes <- bind_rows(phenotypesHL, phenotypesLL)

  frequenciesCombo <- phenotypes %>%
    group_by(plate, transfer, rep, pType) %>%
    count(p) %>%
    mutate(
      p = factor(p, levels = c("U", "S", "A", "B", "A/B", "AB", "E")),
      f = n / sum(n)
      ) %>%
    # complete the table with all possible phenotypes
    complete(
      plate, transfer, rep, pType, p,
      fill = list(n = 0, f = 0))

  meanFrequenciesCombo <- frequenciesCombo %>%
    group_by(plate, transfer, p, pType) %>%
    summarise(
      mean = mean(f),
      n = n(),
      se = sd(f) / sqrt(n),
      max = max(f),
      min = min(f)) %>%
    # drop phenotype E if it doesn't show up
    ungroup() %>%
    filter(!(p == "E" & mean == 0)) %>%
    # order of phenotypes (for legend presentation)
    mutate(
      p = fct_drop(factor(p, levels = names(strainColors)))
    ) %>%
    # add treatment info
    left_join(comboTreatments, by = "plate")

## Fig: fDoubleResistance Combination
  meanFrequenciesComboT5 <- meanFrequenciesCombo %>%
    filter(transfer == 5, p == "AB") %>%
    mutate(
      id = str_c(pType, cA, cB),
      plotA = if_else(cA == 0, 0, log(cA, 2) + 2), plotB = if_else(cB == 0, 0, log(cB, 2) + 2),
      x1 = plotA, y1 = plotB,
      x2 = plotA + 1, y2 = plotB + 1,
      x3 = if_else(pType == "LL", plotA + 1, plotA), y3 = if_else(pType == "LL", plotB, plotB + 1)) %>%
    pivot_longer(x1:y3, names_to = c(".value", "set"), names_pattern = "(.)(.)") %>%
    filter(!(id %in% c("HL00", "LL00")))

  fABcomparison <- ggplot(meanFrequenciesComboT5, aes(x = x, y = y)) +
    geom_polygon(aes(fill = mean, group = id), colour = "white", size = 1) +
    scale_fill_viridis_c("mean\nfrequency of AB",
      breaks = seq(0, 0.18, 0.06)) +
    scale_y_continuous("Nalidixic acid concentration",
      breaks = c(0.5, 1.5, 2.5, 3.5),
      labels = c("0", expression(0.5 %*% MIC), expression(1 %*% MIC), expression(2 %*% MIC))) +
    scale_x_continuous("Streptomycin concentration",
      breaks = c(0.5, 1.5, 2.5, 3.5),
      labels = c("0", expression(0.5 %*% MIC), expression(1 %*% MIC), expression(2 %*% MIC))) +
    plotTheme +
    theme(
      panel.grid.major = element_blank(),
      # legend.position = c(.05, .95),
      legend.justification = c("left", "top"),
      legend.box.just = "right",
      legend.margin = margin(0, 0, 0, 0))

  plotWidth <- 100
  plotHeight <- 75

  ggsave(
    filename = file.path(outDir, "figures", "fig4.tiff"),
    device = "tiff",
    compression = "lzw", type = "cairo",
    dpi = 600,
    plot = fABcomparison,
    width = plotWidth, height = plotHeight, units = "mm")

# Timeseries for SI

  # labels for missing combinations
  noDataLabels <- expand.grid(cA = c(0, 0.5, 1, 2), cB = c(0, 0.5, 1, 2)) %>%
    as_tibble() %>%
    left_join(comboTreatments, by = c("cA", "cB")) %>%
    mutate(
      transfer = 2.5,
      mean = 0.5,
      description = if_else(is.na(plate), "no data", ""))

  plotAllPlatesComboTx <- ggplot(
    data = meanFrequenciesCombo,
    mapping = aes(x = transfer, y = mean, color = p)) +
    facet_grid(
      cols = vars(cA), rows = vars(cB),
      labeller = label_bquote(
        cols = .(cA) %*% MIC,
        rows = .(cB) %*% MIC)) +
    geom_line(size = 0.3) +
    geom_ribbon(aes(x = transfer, ymin = min, ymax = max, fill = p), alpha = 0.25, inherit.aes = FALSE) +
    geom_text(
      data = noDataLabels,
      aes(label = description), color = "grey", fontface = "italic",
      show.legend = FALSE) +
    scale_y_continuous(
        name = "mean phenotype frequency and range (n = 4)",
        limits = c(0, 1), breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
    scale_colour_manual(
      name = "Phenotype",
      values = strainColors,
      guide = guide_legend(nrow = 3)) +
    scale_fill_manual(
      name = "Phenotype",
      values = strainColors,
      guide = "none") +
    plotTheme +
    theme(
      legend.position = c(0.85, 0.75),
      legend.box.background = element_rect(size = 0.3),
      legend.key.size = unit(7, "pt"),
      legend.key.width = unit(7, "pt"),
      legend.key.height = unit(7, "pt"),
      legend.title = element_blank(),
      legend.margin = margin(0, 4, 4, 4),
      legend.justification = c(0, 0),
      legend.direction = "vertical",
      strip.text.x = element_text(size = 7, hjust = 0.5, face = "bold"),
      strip.text.y = element_text(size = 7, hjust = 0.5, face = "bold"),
      plot.margin = unit(c(1.5, 1.5, 1, 1), "lines"))

  plotWidth <- 183
  plotHeight <- 80

  # Timeseries with phenotypes not sensitive to high levels of drugs (i.e. treatment with 2x MIC and/or plates)
  plotDataHL <- filter(meanFrequenciesCombo, pType == "HL") %>%
      mutate(p = fct_relevel(p, levels = c("U", "A", "A/B")))

  ggsave(
    filename = file.path(outDir, "SI", "S6_FigPartA-timeseries-combinationTreatments-HLphenotype.pdf"),
    plot = plotAllPlatesComboTx %+% plotDataHL,
    width = plotWidth, height = plotHeight, units = "mm")

  # Timeseries with phenotypes not sensitive to current level of drugs (i.e. not sensitive to the current treatment)
  plotDataLL <- filter(meanFrequenciesCombo, pType == "LL") %>%
      mutate(p = fct_relevel(p, levels = c("U", "A", "A/B")))
  
  ggsave(
    filename = file.path(outDir,"SI", "S6_FigPartB-timeseries-combinationTreatments-LLphenotype.pdf"),
    plot = plotAllPlatesComboTx %+% plotDataLL,
    width = plotWidth, height = plotHeight, units = "mm")

## Additional Stats - phenotype frequency for diff. combination therapies
  
  frequenciesAllResHL <- frequenciesCombo %>%
    mutate(
      p = factor(p, levels = names(strainColors)[1:6])
      ) %>%
    filter(transfer == 5, p != "E", pType == "LL") %>% 
    mutate(
      pRes = case_when(
        p %in% c("A", "B", "AB", "A/B") ~ "R",
        p == "U" ~ "U",
        p == "S" ~ "S")
    ) %>%
    group_by(plate, transfer, rep, pRes) %>%
    summarize(
      f2 = sum(f)) %>%
    mutate(
      pRes = factor(pRes, levels = c("U", "S", "R")),
      communityPrint = "textstyle(Scenario)~symbol(\"\\306\")",
      community = "Scenario 0"
    )
    
  meanFrequenciesAllResHL <- frequenciesAllResHL %>%
    group_by(community, communityPrint, plate, rep, pRes) %>%
    summarise(
      meanF = mean(f2),
      n = n(),
      se = sd(f2) / sqrt(n),
      max = max(f2),
      min = min(f2))

# interaction
modelFunction <- function(df) {
  lm(meanF ~ p + p * plate, data = df)
}

meanFrequenciesHLRSU <- meanFrequenciesAllResHL %>%
  ungroup() %>%
  group_by(community) %>%
  mutate(
    p = fct_drop(pRes),
    plate = as_factor(plate))

# interaction matrix
tmp <- expand.grid(plate = unique(meanFrequenciesHLRSU$plate),
                   p = unique(meanFrequenciesHLRSU$p))
modelMatrix <- model.matrix(~ p * plate, data = tmp)

tukey <- contrMat(table(meanFrequenciesHLRSU$plate), "Tukey")
tukeyMatrix1 <- cbind(
  tukey,
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)),
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)))
rownames(tukeyMatrix1) <- paste(levels(meanFrequenciesHLRSU$p)[1], rownames(tukeyMatrix1), sep = ":")

tukeyMatrix2 <- cbind(
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)),
  tukey,
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)))
rownames(tukeyMatrix2) <- paste(levels(meanFrequenciesHLRSU$p)[2], rownames(tukeyMatrix2), sep = ":")

tukeyMatrix3 <- cbind(
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)),
  matrix(0, nrow = nrow(tukey), ncol = ncol(tukey)),
  tukey)
rownames(tukeyMatrix3) <- paste(levels(meanFrequenciesHLRSU$p)[3], rownames(tukeyMatrix3), sep = ":")

tukeyMatrix <- rbind(tukeyMatrix1, tukeyMatrix2, tukeyMatrix3)
colnames(tukeyMatrix) <- rep(colnames(tukey), 3)

testHLRSU <- meanFrequenciesHLRSU %>%
  nest() %>%
  mutate(
    model = map(data, modelFunction),
    anova = map(model, anova),
    mcomp = map(model, ~glht(., linfct = tukeyMatrix %*% modelMatrix)),
    mcompSum = map(mcomp, summary),
    plotLetters = map(mcompSum, ~compLetters(.$test)))

names(testHLRSU$anova) <- levels(meanFrequenciesHLRSU$community)
names(testHLRSU$mcompSum) <- levels(meanFrequenciesHLRSU$community)

plotLettersHLRAS <- testHLRSU %>%
  dplyr::select(community, plotLetters) %>%
  unnest(cols = c(plotLetters)) %>%
  mutate(
    plate = as.integer(plate),
    p = factor(p, levels = c("U", "S", "R"))
  )

testHLRSU$anova
testHLRSU$mcompSum

# plot

resColors <- c(
  "U" = "#b3b3b3",
  "S" = "#0cb7eb",
  "R" = "#eb400c")

dataHLRUletters <- meanFrequenciesAllResHL %>%
  group_by(community, communityPrint, plate, pRes) %>%
  summarize(
    meanF = max(meanF) + 0.1
  ) %>%
  left_join(plotLettersHLRAS, by = c("community", "plate", "pRes" = "p")) %>%
  filter(pRes %in% c("U", "R"))

plateLabelsCombo <- c(
  "no treatment",
  "2 MIC A\n2 MIC B",
  "1 MIC A\n2 MIC B",
  "2 MIC A\n1 MIC B",
  "1 MIC A\n1 MIC B",
  "0.5 MIC A\n1 MIC B",
  "1 MIC A\n0.5 MIC B",
  "0.5 MIC A\n0.5 MIC B")

plotData <- filter(meanFrequenciesAllResHL %>% filter(pRes %in% c("U", "R"))) %>%
  left_join(comboTreatments, by = "plate")

freqPlotAllResHL <- ggplot(
  data = plotData,
  aes(x = as.character(plate), y = meanF, colour = pRes)) +
  facet_wrap(vars(communityPrint), ncol = 1, scales = "free", labeller = label_parsed) +
  stat_summary(
    mapping = aes(x = as.character(plate), y = meanF, group = pRes, colour = pRes), inherit.aes = FALSE,
    fun.y = median, fun.ymin = median, fun.ymax = median,
    geom = "errorbar", width = 2 / 3, size = 0.5,
    position = position_dodge(width = 0.75), show.legend = TRUE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1),
    size = 1.5, alpha = 0.5, shape = 20, stroke = 0, show.legend = FALSE) +
  geom_text(
    data = dataHLRUletters,
    mapping = aes(label = compLetter),
    position = position_dodge(width = 0.75),
    vjust = 0.5, hjust = 0.5,
    size = 5 / .pt,
    alpha = 0.5,
    show.legend = FALSE
    ) +
  scale_y_continuous(
    name = "population phenotype frequency",
    limits = c(0, 1.1), breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
  scale_x_discrete(NULL,
    limits = 1:8,
    labels = plateLabelsCombo
  ) +
  scale_colour_manual(
      name = "Phenotype",
      values = resColors) +
  geom_vline(xintercept = seq(1.5, 5.5, 1), size = 0.3, colour = "grey92") +
  plotTheme +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = c(0.75, 0.85),
    legend.box.background = element_rect(size = 0.3),
    legend.key.size = unit(7, "pt"),
    legend.key.width = unit(7, "pt"),
    legend.key.height = unit(7, "pt"),
    legend.title = element_blank(),
    legend.margin = margin(2, 2, 2, 2),
    legend.justification = c(0, 0),
    legend.direction = "horizontal",
    strip.text.x = element_text(size = 8, hjust = 0, face = "bold"))


plotWidth <- 100
plotHeight <- 50

# ggsave(
#   filename = file.path(outDir, "figures", "fig-phenotypeFrequenciesComboLL-R-U.pdf"),
#   plot = freqPlotAllResHL,
#   width = plotWidth, height = plotHeight, units = "mm",
#   device = cairo_pdf)

freqPlotAllResHL
