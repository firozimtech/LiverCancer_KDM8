# http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/casestudy.html#Case_study_n_3:_Integration_of_methylation_and_expression_for_ACC
library(TCGAbiolinks)
library(SummarizedExperiment)
query.met <- GDCquery(
  project = "TCGA-LIHC", 
  legacy = TRUE,
  data.category = "DNA methylation",
  platform = "Illumina Human Methylation 450",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)
GDCdownload(query.met)


lihc.met <- GDCprepare(
  query = query.met,
  save = TRUE, 
  save.filename = "lihcDNAmet.rda",
  summarizedExperiment = TRUE
)

#--------------------------------------------
# STEP 2: Processing DNA methylation
#--------------------------------------------
# To load the RDA data from current working directory
lihc.met <- get(load("../lihcDNAmet.rda"))

# # remove probes with NA (similar to na.omit)
lihc.met <- subset(lihc.met,subset = (rowSums(is.na(assay(lihc.met))) == 0))

# remove probes in chromosomes X, Y and NA
lihc.met <- subset(lihc.met, subset = !as.character(seqnames(lihc.met)) %in% c("chrNA","chrX","chrY"))

# see header of the samples matrix
colnames(SummarizedExperiment::colData(lihc.met))
#sample_type(lihc.met)

# To get the sample type information
aa<- colData(lihc.exp)[c("definition", "sample_type")]
write.csv(aa, "Sample_information_Methyl_LIHC_Legacy.csv")

# make the group of NT(Solid Tissue Normal), and TP(Primary Tumor)
group1 <- TCGAquery_SampleTypes(colnames(lihc.met), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(lihc.met), typesample = c("TP"))

# select the NT and TP, remove others (recurrent tumor)
lihc.met <- lihc.met[,c(group1, group2)]

View(head(lihc.met))
# Volcano plot and DMC file
lihc.met.DMC <- TCGAanalyze_DMC(lihc.met, 
                                groupCol = "sample_type",
                                group1 = "Primary Tumor", #Primary solid Tumor "TP"
                                group2= "Solid Tissue Normal", #Solid Tissue Normal "NT"
                                p.cut = 10^-3, # adjusted P-value
                                diffmean.cut = 0.25,
                                legend = "State",
                                plot.filename = "Meth_volcano_LIHC_TPvsNT_new.png")


# Save lihc.met.DMC
save(lihc.met.DMC, file = "lihc_met_DMC.RData")


# download the methylation probe data
meth_data <- assay(lihc.met)
View(head(meth_data))
write.csv(meth_data, "Methyl_Probes_LIHC_Legacy.csv")

#Transpose the data
meth_dataTrans <- t(meth_data)
meth_dataTrans2<-round(meth_dataTrans, digits =5)

meth_dataTrans2[1:5,1:5]

write.csv(meth_dataTrans2, "Methyl_Probes_Transpose_LIHC_Legacy.csv")

# download the methylation probe genomic information
meth_genome <-rowRanges(lihc.met)
write.csv(meth_genome, "Methyl_Probes_genome_LIHC_Legacy.csv")

#--------------------------------------------
# 2.1 - Mean methylation of samples, and Probes
# -------------------------------------------

TCGAvisualize_meanMethylation(lihc.met,
                              groupCol = "sample_type",
                              group.legend = "TCGA-LIHC",
                              filename = "LIHC_mean_methylation.pdf",
                              print.pvalue = TRUE,
                              width = 5,
                              height = 5,
                              dpi = 600,
                              y.limits = c(0,1),
                              #title = "cg00632019",
                              jitter.size = 1)


##### Figures for methylation ########
# selecting all Methylation probes with "JMJD5" gene
myProbes<-rowData(lihc.met)
#str(myProbes)

# Gene Symbol with pattern match
myProbes2 <- myProbes[grep("JMJD5", myProbes@listData$Gene_Symbol), ]

write.csv(myProbes2, "Methyl_Probes_KDM8.csv")

#-------------------------All Probes Figures in KDM8-----------------------------------------
# Making figures of methylation mean in all probes og KDM8 
myProbes3 <- rownames(myProbes2)
#myProbes3 <- "cg02871891"
colnames(met_probe)
for (probe_name in myProbes3){
  # plot the differential analysis of a probe in different sample, tumor vs normal
  
  # Change the probe_name for different probes, automatically write jpeg name according to probe
  #probe_name = "cg00636124"
  probe_output_figure <- paste("LIHC_KDM8_probe", "-",probe_name,"-Meth.pdf", sep="")
  
  met_probe<-lihc.met[probe_name,]
  TCGAvisualize_meanMethylation(met_probe,
                                groupCol = "sample_type",
                                group.legend = probe_name,
                                filename = probe_output_figure,
                                print.pvalue = TRUE,
                                width = 5,
                                height = 5,
                                dpi = 600,
                                y.limits = c(0,1),
                                legend.ncols = 3,
                                jitter.size = 1)
  
  
}


sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

