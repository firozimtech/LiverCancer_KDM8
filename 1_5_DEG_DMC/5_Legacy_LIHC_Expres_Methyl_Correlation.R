# the progam is to calculate correlation
library(dbplyr)
library(tidyverse)

# read the expression and methylation values
# check.names=FALSE (avoid converting "-" in column name to "." eg "TCGA-12-34" TCGA.12.34")

Patient_exprs_methyl <-Patient_exprs_log_methyl
#Patient_exprs_methyl<-read.table("Final_ExpressNorm_methyl_data.csv", sep=",", header =  T, check.names=FALSE)
#Patient_exprs_log_methyl<-read.table("Final_ExpressNormLog2_methyl_data.csv", sep=",", header =  T, check.names=FALSE)
View(head(Patient_exprs_methyl))



# Select the only "Primary Tumor" gene expression data
Tumor_exp<-Patient_exprs_methyl %>% filter(sample_type_2 == "Primary Tumor")
dim(Tumor_exp)
Tumor_exp<-Tumor_exp[,c(-1,-2,-3,-4)]
View(head(Tumor_exp))
dim(Tumor_exp)

# Select the only "Solid Tissue Normal" gene expression data
Normal_exp<-Patient_exprs_methyl %>% filter(sample_type_2 == "Solid Tissue Normal")
dim(Normal_exp)
Normal_exp<-Normal_exp[,c(-1,-2,-3,-4)]
dim(Normal_exp)
View(head(Normal_exp))





library(Hmisc)
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
#Normal_exp_matrix <-as.matrix(Normal_exp)

Normal_res2<-rcorr(as.matrix(Normal_exp, type="pearson"))
#res2<-rcorr(as.matrix(mtcars[,1:7]))
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
Normal_cor_pval<-flattenCorrMatrix(Normal_res2$r, Normal_res2$P)

# filter the correlation data
Normal_cor_pval_all<-Normal_cor_pval %>% filter(cor >0.3 | cor < -0.3) %>% filter(p < 0.05) 
View(Normal_cor_pval_all)
write.csv(Normal_cor_pval_all, "Gene_Normal_cor_pval_all.csv")


Tumor_res2<-rcorr(as.matrix(Tumor_exp, type="pearson"))
Tumor_cor_pval<-flattenCorrMatrix(Tumor_res2$r, Tumor_res2$P)

# filter the correlation data
Tumor_cor_pval_all<-Tumor_cor_pval %>% filter(cor >0.3 | cor < -0.3) %>% filter(p < 0.05) #%>% filter(row == "JMJD5" | column == "JMJD5")
View(Tumor_cor_pval_all)

write.csv(Tumor_cor_pval_all, "Gene_Cancer_cor_pval_all.csv")
