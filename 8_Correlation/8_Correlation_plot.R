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


#------ correlation analysis -------
library(ggpubr)
# List of probe names
#probe_names <- c("cg03101936", "cg02871891", "cg07361448", "cg00632019")
probe_names <- "cg03101936"

#jpeg(file="LIHC-Tumor-cg03101936.jpeg", width = 5,   height = 5, units="in", res = 600)
jpeg(file="LIHC-Normal-cg03101936.jpeg", width = 5,   height = 5, units="in", res = 600)

p<-ggscatter( #Tumor_exp, x = "JMJD5", y = "cg03101936",
              Normal_exp, x = "JMJD5", y = "cg03101936",
              #title = ("LIHC-Primary Tumor"),
              title = ("LIHC-Solid Tissue Normal"),
              size = 2, # point size
             add = "reg.line", # Add regression line
             conf.int = TRUE,  # Add confidence interval
             #cor.coef = TRUE, cor.method = "pearson",
             xlab = "KDM8 log2(normalized+1)", ylab = "cg03101936 (Beta-value)")


p + theme(axis.title=element_text(size=10)) +
  theme(title=element_text(size=10)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(method = "pearson", label.x = 10, label.y =1)  # Add correlation coefficient, print p-value t specific position

dev.off()


# multiplr regression
# Select columns starting with 'cg' from normal
cg_cols <- grep("^cg", names(Normal_exp), value = TRUE)
cg_cols
normal_data_cg <- Normal_exp[, c("JMJD5",cg_cols)]
normal_data_cg
View(normal_data_cg)
normal_model <- lm(JMJD5 ~ ., data = normal_data_cg)
summary(normal_model)


tumor_data_cg <- Tumor_exp[, c("JMJD5",cg_cols)]
tumor_data_cg
View(tumor_data_cg)
tumor_model <- lm(JMJD5 ~ ., data = tumor_data_cg)
summary(tumor_model)

