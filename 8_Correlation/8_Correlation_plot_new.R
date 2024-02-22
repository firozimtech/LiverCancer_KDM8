# differential correlation analysis using DGCA
# https://www.rdocumentation.org/packages/DGCA/versions/1.0.2
library(dbplyr)
library(tidyverse)
library(ggplot2)
library(DGCA)
library(plotrix)
library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)

#-------------------- prepare the clean input dataset -------------
# there are two final datasets
# (1) dataset1 - matrix of gene expression: row as gene, column as samples
# (2) dataset2 - design matrix showing samples catagery(vancer,noncancer) as column, rowshing sample belongs to 


# read the expression and methylation values
# check.names=FALSE (avoid converting "-" in column name to "." eg "TCGA-12-34" TCGA.12.34")

#Patient_exprs_methyl <-Patient_exprs_log_methyl
Patient_exprs_methyl  <-read.table("Final_ExpressNormLog2_methyl_data.csv", sep=",", header =  T, check.names=FALSE)
class(Patient_exprs_methyl)
dim(Patient_exprs_methyl)
View(head(Patient_exprs_methyl))

Patient_exprs_methyl<-Patient_exprs_methyl[,-1]
# change the col name
colnames(Patient_exprs_methyl)[which(names(Patient_exprs_methyl) == "INS-IGF2")] <- "INSIGF2"

#change row names
rownames(Patient_exprs_methyl)<-Patient_exprs_methyl[,2]
Patient_exprs_methyl<-Patient_exprs_methyl[,-c(1,2,3)]

# Select the only "Primary Tumor" gene expression data
Tumor_exp<-Patient_exprs_methyl %>% filter(sample_type_2 == "Primary Tumor")
dim(Tumor_exp)
View(head(Tumor_exp))

# transpose the data
t_Tumor_exp<-t(Tumor_exp)
#t_Tumor_exp<-t_Tumor_exp[-c(1),]
View(head(t_Tumor_exp))
dim(t_Tumor_exp)

# Select the only "Solid Tissue Normal" gene expression data
Normal_exp<-Patient_exprs_methyl %>% filter(sample_type_2 == "Solid Tissue Normal")
dim(Normal_exp)
View(head(Normal_exp))

# transpose the data
t_Normal_exp<-t(Normal_exp)
#t_Normal_exp<-t_Normal_exp[-c(1),]
View(head(t_Normal_exp))
dim(t_Normal_exp)


#------ correlation analysis for multiple probes for both tumor and normal -------
library(ggpubr)

# List of probe names
probe_names <- c("cg00632019", "cg01755539", "cg02411582", "cg02871891", "cg03101936", "cg07361448")

# Loop through each probe
for (probe in probe_names) {
  # Define the output file names for tumor and normal samples
  tumor_output_file <- paste0("LIHC-Tumor-", probe, ".jpeg")
  normal_output_file <- paste0("LIHC-Normal-", probe, ".jpeg")
  
  #print(paste("****", "", tumor_output_file))
  #print(paste("****", "", normal_output_file))
  
  ##### Generate the plot for tumor samples #####
  jpeg(file = tumor_output_file, width = 5, height = 5, units = "in", res = 600)
  
  p_tumor <- ggscatter(
    Tumor_exp, x = "JMJD5", y = probe,
    title = "LIHC-Primary Tumor",
    size = 2, # point size
    add = "reg.line", # Add regression line
    conf.int = TRUE,  # Add confidence interval
    xlab = "KDM8 log2(normalized+1)",
    ylab = paste0(probe, " (Beta-value)")
    )
  
  pp_tumor<- (p_tumor + theme(axis.title = element_text(size = 10)) +
    theme(title = element_text(size = 10)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_cor(method = "pearson", label.x = 10, label.y = 1) # Add correlation coefficient, print p-value t specific position
    )
  
    print(pp_tumor)
    dev.off()

  
  
  #### Generate the plot for normal samples ####
  jpeg(file = normal_output_file, width = 5, height = 5, units = "in", res = 600)
  
  p_normal <- ggscatter(
    Normal_exp, x = "JMJD5", y = probe,
    title = "LIHC-Solid Tissue Normal",
    size = 2, # point size
    add = "reg.line", # Add regression line
    conf.int = TRUE, # Add confidence interval
    xlab = "KDM8 log2(normalized+1)",
    ylab = paste0(probe, " (Beta-value)")
  )
  
  pp_normal<-(p_normal + theme(axis.title = element_text(size = 10)) +
    theme(title = element_text(size = 10)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_cor(method = "pearson", label.x = 10, label.y = 1) # Add correlation coefficient, print p-value t specific position
    )
  
    print(pp_normal)
    dev.off()
  
}


