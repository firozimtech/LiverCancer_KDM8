# http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/casestudy.html#Case_study_n_3:_Integration_of_methylation_and_expression_for_ACC
library(SummarizedExperiment)
library(TCGAbiolinks)
query.exp <- GDCquery(project = "TCGA-LIHC", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))
GDCdownload(query.exp)
lihc.exp <- GDCprepare(query = query.exp, save = TRUE, save.filename = "lihcExp.rda")

# get subtype information
dataSubt <- TCGAquery_subtype(tumor = "LIHC")

# get clinical data
dataClin <- GDCquery_clinic(project = "TCGA-LIHC","clinical") 

# Which samples are Primary Tumor
dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP") 
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")


dataPrep <- TCGAanalyze_Preprocessing(
  object = lihc.exp, 
  cor.cut = 0.6
)                      

dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfo,
  method = "gcContent"
)                

dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)   

dataDEGs <- TCGAanalyze_DEA(
  mat1 = dataFilt[,dataSmNT],
  mat2 = dataFilt[,dataSmTP],
  Cond1type = "Normal",
  Cond2type = "Tumor",
  fdr.cut = 1 ,
  logFC.cut = 0.00004,
  method = "glmLRT"
)  

# save DEGS
write.csv(dataDEGs, "DGE_LIHC_LEGACY.csv")

dataDEGs["JMJD5",]
# save the expression data normalised
exp_data <- assay(dataFilt)
write.csv(dataFilt, "Eexpression_Normalized_LIHC_LEGACY.csv")
dataFiltTrans<-t(dataFilt)
View(head(dataFiltTrans))

write.csv(dataFiltTrans, "Eexpression_Normalized_Transpose_LIHC_LEGACY.csv")

dataFiltTransLog<-log2(dataFiltTrans+1)
write.csv(dataFiltTransLog, "Eexpression_Normalized_Transpose_log2+1_LIHC_LEGACY.csv")
dataFiltTransLog[1:3,1:3]

colnames(SummarizedExperiment::colData(lihc.exp))

rowRanges(lihc.exp)
GeneSymbol_information<-as.data.frame(rowRanges(lihc.exp))
write.csv(GeneSymbol_information, "GeneSymbol_information_LIHC_LEGACY.csv")

#---- TCGA volcano plot
# to increase tje maximum overlap gene name label in Volcano plot
# options(ggrepel.max.overlaps = Inf)
options(ggrepel.max.overlaps = 50)
TCGAVisualize_volcano(x = dataDEGs$logFC,
                      y = dataDEGs$FDR,
                      filename = "LIHC_Volcano_Exp2.png",
                      x.cut = 2,
                      y.cut = 10^-3,
                      names = rownames(dataDEGs),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      show.names = "significant",
                      #max.overlaps = 1,
                      title = "Volcano plot (Primary Tumor vs Solid Tissue Normal)",
                      width = 10)
