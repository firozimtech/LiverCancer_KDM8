library(tidyverse)
library(dplyr)

# upload normalized gene expression data
Gene.exp<-read.table("Eexpression_Normalized_Transpose_LIHC_LEGACY.csv", sep=",", header =  T, check.names=FALSE)
# change colnames
colnames(Gene.exp)[1] <-"cases"
View(head(Gene.exp))

# upload patient information of gene expression 
patients<-read.table("Expression_patient_information.csv", sep=",", header =  T, check.names=FALSE)
class(patients)
dim(patients)
View(head(patients))

# merge two files based upon common "patients"
#patients.exp<-inner_join(patients, Gene.exp, by = "cases", all = FALSE)
patients.exp<-inner_join(patients, Gene.exp, by = "cases")

View(head(patients.exp))

dim(patients)
dim(Gene.exp)
dim(patients.exp)



# Running the test on unpaired samples would run the Mann-Whitney U Test. 
# This is also called the Mann-Whitney-Wilcoxon test, which tests differences in the magnitude between groups.
# In R, you would use wilcox.test(x, y, paired=FALSE) to run this test.

library(rstatix)
library(ggpubr)

patients.exp[,c("JMJD5","sample_type")]

# compare mean by wilcox Test
compare_means(JMJD5 ~ sample_type, data = patients.exp)


jpeg(file="plot_JMJD5.jpeg",
     width=5, height=5, units="in", res=600)

p<-ggboxplot(patients.exp, x = "sample_type", y = "JMJD5", 
          color = "sample_type", palette = c("red", "deepskyblue3"),
          ylab = "KDM8 (Normalized)", xlab = "sample_type", add = "jitter", 
          #legend title, labels and position: legend = "right"
          title = "Gene expression")
#p + (plot.title = element_text(hjust = 0.5))  + stat_compare_means()
p  + stat_compare_means()
dev.off()
# unpaired test
wilcox.test( patients.exp$sample_type, patients.exp$JMJD5, paired=FALSE)


