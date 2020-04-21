# common paramaters

strainColors <- c(
  "U" = "#b3b3b3",
  "S" = "#0cb7eb",
  "A" = "#ebe507",
  "B" = "#ffa10a",
  "A/B" = "#662f89",
  "AB" = "#acdf0c",
  "E" = "#e53b6a")

plateLabels <- str_c(
  letters[1:6], ") ",
  c("no treatment", "mono A", "mono B", "combination", "cycling", "mixing"))
names(plateLabels) <- 1:6

plotTheme <- theme_bw() +
  theme(
    # overall text size
    text = element_text(size = 7),
    # legend
    legend.position = "right",
    legend.key = element_rect(color = NA),
    # facet labels
    strip.background = element_rect(fill = NA, color = NA),
    strip.text.x = element_text(size = 8, hjust = 0, face = "bold"),
    strip.text.y = element_text(size = 7, hjust = 0.5),
    # panel
    panel.border = element_rect(color = NA),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_line(size = 0.2),
    # axis
    axis.line = element_line(size = 0.2),
    axis.ticks = element_line(size = 0.2),
    axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
    axis.title.y = element_text(margin = margin(0, 10, 0, 0))
  )
