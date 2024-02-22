library(tidyverse)
library(rstatix)
library(ggpubr)

# Methylation probes data
Methyl_Probes <- read.csv("Methyl_Probes_Transpose_LIHC_Legacy_KDM8_sample_type.csv", sep = ",", header = TRUE, check.names = FALSE)

# List of probes to analyze
probe_names <- c("cg00632019", "cg01755539", "cg02411582", "cg02871891", "cg03101936", "cg07361448")

# Loop through each probe
for (probe in probe_names) {
  
  # Print the current probe
  print(paste("****", "", probe))
  
  # Compare mean DNA methylation by sample type
  mean_comparison <- compare_means(as.formula(paste(probe, "~ sample_type")), data = Methyl_Probes)
  
  # Generate the boxplot
  jpeg(file = paste0("plot_", probe, ".jpeg"), width = 5, height = 5, units = "in", res = 600)
  
  p <- ggboxplot(
    data = Methyl_Probes,
    x = "sample_type",
    y = probe,
    color = "sample_type",
    palette = c("red", "deepskyblue3"),
    ylab = "Mean DNA methylation (Beta-value)",
    xlab = "Sample Type",
    add = "jitter",
    ylim = c(0, 1),
    title = probe
  )
  
  pp<-(p + theme(plot.title = element_text(hjust = 0.5, size = 10)) +
    stat_compare_means())
  
  print(pp)
  dev.off()
  
  # Perform unpaired Wilcoxon test
  #wilcox.test(Methyl_Probes$sample_type, Methyl_Probes[[probe]], paired = FALSE)
}
