# name of patinets having both expression and methylation data
comm_patients_Methyl_Expr<- read.csv("Common_Patients_Express_Methy_2.csv", sep = ",", header = TRUE, check.names=FALSE)

# methylation probes data
Methyl_Probes<-read.csv("Methyl_Probes_Transpose_LIHC_Legacy_KDM8_sample_type.csv", sep = ",", header = TRUE, check.names=FALSE)

# gene expression data (log2(Normalized+1))
Gene_expr_log<-read.csv("Expression_Normalized_Transpose_log2+1_LIHC_LEGACY.csv", sep = ",", header = TRUE, check.names=FALSE)

# gene expression data (Normalized)
Gene_expr<-read.csv("Expression_Normalized_Transpose_LIHC_LEGACY.csv", sep = ",", header = TRUE, check.names=FALSE)

DGEs_list <- read.csv("../DGE_LIHC_LEGACY_UPreg_Downreg.csv", sep = ",", header = TRUE)
GeneSymbol<-DGEs_list$GeneSymbol
View(DGEs_list)
dim(DGEs_list)


dim(comm_patients_Methyl_Expr)
dim(Methyl_Probes)
dim(Gene_expr_log)
dim(Gene_expr)

View(head(comm_patients_Methyl_Expr))
View(head(Methyl_Probes))
View(head(Gene_expr_log))
View(head(Gene_expr))


### (A)merge sample information+log2(Normalized expression+1) + Methylation data

#replace the wrong name "Sep-14" with "SEP14"
#replace the wrong name "Sep-03" with "SEP3"
names(Gene_expr_log)[names(Gene_expr_log) == "Sep-14"] <-"SEPT14"
names(Gene_expr_log)[names(Gene_expr_log) == "Sep-03"] <-"SEPT3"

Gene_expr_log[,"SEP14"]

# select only gene showing DEGs
Gene_expr_log <-Gene_expr_log %>% select("cases_express",GeneSymbol)
View(head(Gene_expr_log))
dim(Gene_expr_log)



Patient_exprs_log <-inner_join(comm_patients_Methyl_Expr, Gene_expr_log, by ="cases_express" )
View(head(Patient_exprs_log))
dim(Patient_exprs_log)

Patient_exprs_log_methyl<-inner_join (Patient_exprs_log, Methyl_Probes, by ="cases_methyl" )
dim(Patient_exprs_log_methyl)

Patient_exprs_log_methyl<-Patient_exprs_log_methyl[,c(-3,-4)]
View(head(Patient_exprs_log_methyl))
write.csv(Patient_exprs_log_methyl, "Final_ExpressNormLog2_methyl_data.csv", sep = ",")
dim(Patient_exprs_log_methyl)

### (B) merge sample information+ Normaized expression + Methylation data

#replace the wrong name "Sep-14" with "SEP14"
#replace the wrong name "Sep-03" with "SEP3"
names(Gene_expr)[names(Gene_expr) == "Sep-14"] <-"SEPT14"
names(Gene_expr)[names(Gene_expr) == "Sep-03"] <-"SEPT3"

Gene_expr[,"SEPT3"]

# select only gene showing DEGs
Gene_expr <-Gene_expr %>% select("cases_express",GeneSymbol)
View(head(Gene_expr))
dim(Gene_expr)


Patient_exprs <-inner_join(comm_patients_Methyl_Expr, Gene_expr, by ="cases_express" )
View(head(Patient_exprs))
dim(Patient_exprs)

Patient_exprs_methyl<-inner_join (Patient_exprs, Methyl_Probes, by ="cases_methyl" )
dim(Patient_exprs_methyl)

Patient_exprs_methyl<-Patient_exprs_methyl[,c(-3,-4)]
View(head(Patient_exprs_methyl))
write.csv(Patient_exprs_methyl, "Final_ExpressNorm_methyl_data.csv", sep = ",")

dim(Patient_exprs_methyl)
