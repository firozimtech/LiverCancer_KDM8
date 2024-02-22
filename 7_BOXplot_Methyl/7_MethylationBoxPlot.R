library(tidyverse)
library(dplyr)


# methylation probes data
Methyl_Probes<-read.csv("Methyl_Probes_Transpose_LIHC_Legacy_KDM8_sample_type.csv", sep = ",", header = TRUE, check.names=FALSE)

View(head(Methyl_Probes))



# Running the test on unpaired samples would run the Mann-Whitney U Test. 
# This is also called the Mann-Whitney-Wilcoxon test, which tests differences in the magnitude between groups.
# In R, you would use wilcox.test(x, y, paired=FALSE) to run this test.

library(rstatix)
library(ggpubr)

Methyl_Probes[,c("cg02871891","sample_type")]

# compare mean by wilcox Test
compare_means(cg02871891 ~ sample_type, data = Methyl_Probes)

#cg02871891, cg03101936, and cg02411582
jpeg(file="plot_cg02411582.jpeg",
     width=5, height=5, units="in", res=600)

p<-ggboxplot(Methyl_Probes, x = "sample_type", y = "cg02411582", 
          color = "sample_type", palette = c("red", "deepskyblue3"),
          ylab = "Mean DNA methylation (B-value)", xlab = "sample_type", add = "jitter", ylim=c(0,1),
          #legend title, labels and position: legend = "right"
          title = "cg02411582")
#p + (plot.title = element_text(hjust = 0.5))  + stat_compare_means()
p  + stat_compare_means()
dev.off()
# unpaired test
wilcox.test( patients.exp$sample_type, patients.exp$JMJD5, paired=FALSE)


