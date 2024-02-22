library(tidyverse)
library(dplyr)
# Get all patients that have DNA methylation and gene expression.
patients.met<-data.frame(getResults(query.met, cols = c("cases.submitter_id","cases", "sample_type")))
patients.exp<-data.frame(getResults(query.exp, cols = c("cases.submitter_id", "cases", "sample_type")))

write.csv(patients.met, "Methyl_patient_information.csv")
write.csv(patients.exp, "Expression_patient_information.csv")
View(head(patients.exp))
dim(patients.exp)

View(head(patients.met))
dim(patients.met)

xx<-(patients.met %>% count(sample_type))
yy<-patients.exp %>% count(sample_type)

# merge two files based upon common "patients"
patients.exp.meth<-inner_join(patients.exp, patients.met, by = "cases.submitter_id", all = FALSE)
write.csv(patients.exp.meth, "Expression_methy_patient_information.csv")

write.csv(myProbes2, "Methyl_Probes_KDM8.csv")
intersect(patients.exp, patients.met)

#patients.exp %in% patients.met
#common.patients <- intersect(
#  substr(getResults(query.met, cols = "cases"), 1, 12),
#  substr(getResults(query.exp, cols = "cases"), 1, 12)
#)

Common_Patients<-read.csv("Common_Patients_Express_Methy_2.csv", sep = ",", header = TRUE)
View(Common_Patients)
