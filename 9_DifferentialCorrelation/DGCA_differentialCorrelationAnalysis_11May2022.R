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
Patient_exprs_methyl<-read.table("Final_ExpressNormLog2_methyl_data.csv", sep=",", header =  T, check.names=FALSE)
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

# merge two data ("Primary Tumor", and "Solid Tissue Normal") according to common gene name
final_tum_norm_data<-merge(t_Tumor_exp, t_Normal_exp, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names")


View(head(final_tum_norm_data))
rownames(final_tum_norm_data)<-final_tum_norm_data[,1]
final_tum_norm_data<-final_tum_norm_data[,-1]
#design_matrix<-final_tum_norm_data["sample_type_2",]
dim(final_tum_norm_data)

group<-final_tum_norm_data["sample_type_2",]
group<-t(group)
rownames(group)<-NULL
View(head(group))
design_matrix<-model.matrix(~ group + 0)

colnames(design_matrix)<-c("PrimaryTumor","SolidTissueNormal")
View(design_matrix)

# drop the row name "sample_type_2"
final_tum_norm_data<-final_tum_norm_data[!(row.names(final_tum_norm_data) %in% c("sample_type_2")), ]
dim(final_tum_norm_data)
View(final_tum_norm_data)

View(sapply(final_tum_norm_data, class))


#--------------------2: START Differential Correlation Analysis using DDCA -------------
####-- Analysis _1 (selected correlation)
####-- Calculate the Differential Correlation of one gene to all genes present in the matrxi 
####-- below list of genes and probes we are interested-------
#interested_gene<-c("JMJD5", "cg02871891", "cg03101936",  "cg02411582", "cg27118526", "cg10221365", "cg06329197", "cg00636124", "cg16752029")
interested_gene<-c("MFSD2A","JMJD5","cg03101936" ,"cg27118526" ,"cg07361448" ,"cg16752029" ,"cg00636124" ,"cg00632019" ,"cg02871891" ,"cg10221365" ,"cg02411582" ,"cg08572679" ,"cg06329197")
# get all correlation result as lists
output_cor_genes <- list()
output_cor_genes_2 <-list()

for (cor_gene in interested_gene){
  print(cor_gene)
  # differential gene expression 
  output_cor_1 <- paste("ddcor_res_",cor_gene,  sep="")
  #print(output_cor_1)
  # Always put 'compare = c("SolidTissueNormal","PrimaryTumor")' 
  # otherwise gain of correlation and loss of correlation will wrong
  output_cor_1 <- ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                                compare = c("SolidTissueNormal","PrimaryTumor"),
                                adjust = "perm", nPerm = 5, splitSet = cor_gene)
  ## Append all result as lists
  output_cor_genes[[cor_gene]]<-output_cor_1
  
  # Print the all correlation 
  output_file_1 <- paste("correlation-",cor_gene,"-vs-All_gene.csv",  sep="")
  write.csv(output_cor_1, output_file_1)
  
  # filter the correlation having p-value below 0.5
  #output_cor_2 <- output_cor_1  %>% filter(PrimaryTumor_pVal < 0.05) %>% filter(SolidTissueNormal_pVal < 0.05)
  
  # filter the correlation having p-value below 0.5 in either group
  output_cor_2 <- output_cor_1  %>% filter(PrimaryTumor_pVal < 0.05 | SolidTissueNormal_pVal < 0.05) %>% filter (zScoreDiff < -1 | zScoreDiff > 1)
  
  output_file_2 <- paste("correlation_filtere_pval05_eithergroup-",cor_gene,"-vs-All_gene.csv",  sep="")
  write.csv(output_cor_2, output_file_2)
  
  ## Append all result with pvalue 0.5 as lists
  output_cor_genes_2[[cor_gene]]<-output_cor_2
}

# to get the correlation result of JMJD5 with all genes
output_cor_genes[["JMJD5"]]
output_cor_genes_2[["JMJD5"]]
typeof(output_cor_genes)


# All genes
all_gene<-unique(c(output_cor_genes_2[["JMJD5"]]$Gene1, output_cor_genes_2[["cg03101936"]]$Gene1, output_cor_genes_2[["cg02871891"]]$Gene1, output_cor_genes_2[["cg02411582"]]$Gene1))
write.csv(all_gene, "All_gene.csv")
# find the common gene in all dataframe having negative correlation with KDM8
negative_correl<-c("JMJD5","cg03101936", "cg02871891", "cg02411582", "cg10221365" , "cg07361448" )
positive_correl<-c("JMJD5","cg00632019", "cg00636124")
                   
vv<-inner_join(output_cor_genes_2[["JMJD5"]], output_cor_genes_2[["cg03101936"]], by = "Gene1" )
library(plyr)
#yx<-join_all(list(output_cor_genes_2[["JMJD5"]], output_cor_genes_2[["cg03101936"]], output_cor_genes_2[["cg02871891"]],output_cor_genes_2[["cg02411582"]], output_cor_genes_2[["cg10221365"]], output_cor_genes_2[["cg07361448"]]  ), by = "Gene1", type='inner')
yx<-join_all(list( output_cor_genes_2[["cg03101936"]], output_cor_genes_2[["JMJD5"]]  ), by = "Gene1", type='inner')

write.csv(yx, "common.csv")
View(yx)
#---------------------- END DECA analysis--------------------------#

#### Analysis _3 (correlation plot) ##### Important
# To plot the differential correlations between RTN4 and its top target, you can use this function:
# [1]========FOR  KDM8 vs Prbes
# Make plots
plot_list <- list()
interested_gene_correlation<-c("cg03101936" ,"cg27118526" ,"cg07361448" ,"cg16752029" ,"cg00636124" ,"cg00632019" ,"cg02871891" ,"cg10221365" ,"cg02411582" ,"cg08572679" ,"cg06329197")
#interested_gene_correlation<-c( "cg02871891", "cg03101936",  "cg02411582", "cg27118526", "cg10221365", "cg06329197", "cg00636124", "cg16752029")
  for (cor_gene_plot in interested_gene_correlation){
    print(cor_gene_plot)
    ylabel <- paste("DNA Methylation ",cor_gene_plot,  sep="")
    p = plotCors(inputMat = final_tum_norm_data, design = design_matrix,
           compare = c("PrimaryTumor","SolidTissueNormal"),
           #geneA =  "JMJD5", geneB= cor_gene_plot, log = FALSE,
           #xlab = "log2Expression KDM8", ylab = ylabel)
  
           geneA =  "DLK1", geneB= cor_gene_plot, log = FALSE,
           xlab = "log2Expression DLK1", ylab = ylabel)
   
   plot_list[[cor_gene_plot]] = p
  }


# Save plots to jpeg. Makes a separate file for each plot.
for (cor_gene_plot in interested_gene_correlation){
  print(cor_gene_plot)
  myfilename <- paste("KDM8_",cor_gene_plot,"_plot.jpeg",  sep="")
  jpeg(file = myfilename, width=11, height=5, units="in", res=600)
  print(plot_list[[cor_gene_plot]])
  dev.off()
}



# [2]========FOR  KDM8 vs Genes
# Make plots
plot_list <- list()
interested_gene_correlation<-c( "NFKBIL2","GPM6A", "RAD54B","HELLS", "MFSD2A" )
for (cor_gene_plot in interested_gene_correlation){
  print(cor_gene_plot)
  ylabel <- paste("log2 Expression ",cor_gene_plot,  sep="")
  p = plotCors(inputMat = final_tum_norm_data, design = design_matrix,
               compare = c("PrimaryTumor","SolidTissueNormal"),
               geneA =  "JMJD5", geneB= cor_gene_plot, log = FALSE,
               xlab = "log2 Expression KDM8", ylab = ylabel)
  
  plot_list[[cor_gene_plot]] = p
}


# Save plots to jpeg. Makes a separate file for each plot.
for (cor_gene_plot in interested_gene_correlation){
  print(cor_gene_plot)
  myfilename <- paste("KDM8_",cor_gene_plot,"_plot.jpeg",  sep="")
  jpeg(file = myfilename, width=11, height=5, units="in", res=600)
  print(plot_list[[cor_gene_plot]])
  dev.off()
}




# [3]========FOR  DEGs vs Prbes
# Make plots
plot_list <- list()
interested_gene_correlation<-c("cg03101936" ,"cg27118526" ,"cg07361448" ,"cg16752029" ,"cg00636124" ,"cg00632019" ,"cg02871891" ,"cg10221365" ,"cg02411582" ,"cg08572679" ,"cg06329197")

for (cor_gene_plot in interested_gene_correlation){
  print(cor_gene_plot)
  ylabel <- paste("DNA Methylation ",cor_gene_plot,  sep="")
  p = plotCors(inputMat = final_tum_norm_data, design = design_matrix,
               compare = c("PrimaryTumor","SolidTissueNormal"),
               geneA =  "JMJD5", geneB= cor_gene_plot, log = FALSE,
               xlab = "log2Expression KDM8", ylab = ylabel)
  
  plot_list[[cor_gene_plot]] = p
}


# Save plots to jpeg. Makes a separate file for each plot.
for (cor_gene_plot in interested_gene_correlation){
  print(cor_gene_plot)
  myfilename <- paste("KDM8_",cor_gene_plot,"_plot.jpeg",  sep="")
  jpeg(file = myfilename, width=11, height=5, units="in", res=600)
  print(plot_list[[cor_gene_plot]])
  dev.off()
}



#--------------------2: START Gene Ontology analysis of Deferentially correlated gene -------------
####-- Input data: output_cor_genes[["JMJD5"]]-------

all<-names(output_cor_genes_2)
for (geneName in all){
    print(geneName)
  
  # create name of object eg"dcorGO_JMJD5"
  geneddcorGO <- paste("ddcorGO_",geneName,  sep="")
  geneddcorGO = ddcorGO(output_cor_genes_2[[geneName]], universe = rownames(final_tum_norm_data), 
                        pval_GO_cutoff = 0.05, pval_gene_thresh = 0.05, ddcor_find_significant = TRUE,
                        gene_ontology = "all", HGNC_clean = TRUE, HGNC_switch = TRUE, annotation = "org.Hs.eg.db", calculateVariance = TRUE)

#View(ddcorGO_geneName)
# Print significant_gain_of_correlation_genes
gain_cor_gene_output <- paste("gain_of_correlation_genes_",geneName,".csv",  sep="")
gain_cor_gene<-geneddcorGO[["significant_gain_of_correlation_genes"]]
write.csv(gain_cor_gene, gain_cor_gene_output)

#Print significant_loss_of_correlation_genes
loss_cor_gene_output <- paste("loss_of_correlation_genes_",geneName,".csv",  sep="")
loss_cor_gene<-geneddcorGO[["significant_loss_of_correlation_genes"]]
write.csv(loss_cor_gene, loss_cor_gene_output)

## name of the output files for enrichment analysis of "gain of correlation gene"
enrich_BP_gain_cor_gene <- paste("enrich_BP_gain_of_correlation_genes",geneName,".csv",  sep="")
enrich_CC_gain_cor_gene <- paste("enrich_CC_gain_of_correlation_genes",geneName,".csv",  sep="")
enrich_MF_gain_cor_gene <- paste("enrich_MF_gain_of_correlation_genes",geneName,".csv",  sep="")

write.csv(geneddcorGO[["enrichment_significant_gain_of_correlation_genes"]][["BP"]], enrich_BP_gain_cor_gene)
write.csv(geneddcorGO[["enrichment_significant_gain_of_correlation_genes"]][["CC"]], enrich_CC_gain_cor_gene)
write.csv(geneddcorGO[["enrichment_significant_gain_of_correlation_genes"]][["MF"]], enrich_MF_gain_cor_gene)

## name of the output files for enrichment analysis for "loss of correlation gene"
enrich_BP_loss_cor_gene <- paste("enrich_BP_loss_of_correlation_genes",geneName,".csv",  sep="")
enrich_CC_loss_cor_gene <- paste("enrich_CC_loss_of_correlation_genes",geneName,".csv",  sep="")
enrich_MF_loss_cor_gene <- paste("enrich_MF_loss_of_correlation_genes",geneName,".csv",  sep="")

write.csv(geneddcorGO[["enrichment_significant_loss_of_correlation_genes"]][["BP"]], enrich_BP_loss_cor_gene)
write.csv(geneddcorGO[["enrichment_significant_loss_of_correlation_genes"]][["CC"]], enrich_CC_loss_cor_gene)
write.csv(geneddcorGO[["enrichment_significant_loss_of_correlation_genes"]][["MF"]], enrich_MF_loss_cor_gene)



}

help(plotGOTwoGroups)

dfList1 = geneddcorGO$enrichment_significant_gain_of_correlation_genes
typeof(dfList1)
dfList2 = geneddcorGO$enrichment_significant_loss_of_correlation_genes
library (plyr)
dfList1 <- ldply (dfList1, data.frame)
dfList2 <- ldply (dfList2, data.frame)
mergeGO<-merge(dfList1, dfList2, by="Term", all=TRUE)
mergeGO2<-full_join(dfList1, dfList2, by="Term")
common_GO<-inner_join(dfList1, dfList2, by="GOBPID")

dim((mergeGO))
dim((mergeGO2))
View(common_GO)



View((mergeGO2))
dfList1<-read.table("enrich_BP_loss_of_correlation_genesJMJD5.csv", sep = ",", header = TRUE)
dfList2<-read.table("enrich_BP_gain_of_correlation_genesJMJD5.csv", sep =",", header = TRUE)

plotGOTwoGroups(dfList1, dfList2, nTerms = 5, minSize = 40,
                maxSize = 1000, labelsCol = "Ontology", adjustPVals = TRUE,
                plotrix_gap = 20, GOTermTypes = c("BP", "CC", "MF"), pValCutoff = 0.05,
                filterSignificant = FALSE, filterSigThresh = 0.05,
                labels = c("Corr Class 1", "GO Term Name", "Corr Class 2"),
                fill_zero_cats = FALSE)

plotGOTwoGroups(dfList1 = geneddcorGO$enrichment_significant_gain_of_correlation_genes,
                dfList2 = geneddcorGO$enrichment_significant_loss_of_correlation_genes,
                nTerms = 5, minSize = 2,
                maxSize = 1000, labelsCol = "Ontology", adjustPVals =  TRUE,
                plotrix_gap = 20, GOTermTypes = c("BP", "CC", "MF"), pValCutoff = 0.1,
                #plotrix_gap = 20, GOTermTypes = c("MF"), pValCutoff = 0.05,
                filterSignificant = FALSE, filterSigThresh = 0.05,
                labels = c("Gain in Correlation", "GO Term Name", "Loss in Correlation"),
                fill_zero_cats = FALSE)

####----------------- END GO Analysis ------------------------

## https://cran.r-project.org/web/packages/MEGENA/vignettes/MEGENA_pipeline_10062016.html
library(MEGENA)
#names(output_cor_genes)
#ddcor_res<-output_cor_genes$JMJD5 # showing error

ddcor_res = ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                     compare = c("PrimaryTumor","SolidTissueNormal"),
                     adjust = "none", heatmapPlot = FALSE, nPerm = 0, nPairs = "all")
                     

# ----------ERROR NOTE-----
# the "ddcor_res" contains "NA" value which gives error while runing ddMEGENA
# remove na in r - remove rows - na.omit function / option
ddcor_res <- na.omit(ddcor_res) 
str(ddcor_res)

megena = ddMEGENA(ddcor_res, adjusted = FALSE, evalCompactness = TRUE)
names(megena_res)
View(megena_res)

# ------ Module-based differential correlation-----------
modules<-megena_res[["modules"]]
moduleDC_res = moduleDC(inputMat = final_tum_norm_data, design = design_matrix,
                        compare = c("PrimaryTumor","SolidTissueNormal"), genes = modules$Genes,
                        labels = modules$Modules, nPerm = 50, number_DC_genes = 3,
                        dCorAvgMethod = "median")
head(moduleDC_res)

# ------END  Module-based differential correlation-----------


#------- plot module---------
# module hub table
module_hub_table<-megena_res$summary[["module.table"]]
write.csv(module_hub_table, "Megena_modulehub_table_all-cor.csv")

### ---following not working
summary.output<-megena_res[["summary"]]
g <- graph.data.frame(el,directed = FALSE)
?graph.data.frame
pnet.obj <- plot_module(output.summary = summary.output, PFN = gg,subset.module = "c1_2",
                        layout = "kamada.kawai",label.hubs.only = TRUE,
                        gene.set = NULL,color.code =  "grey",
                        output.plot = FALSE,out.dir = "modulePlot",col.names = c("magenta","green","cyan"),label.scaleFactor = 20,
                        hubLabel.col = "black",hubLabel.sizeProp = 1,show.topn.hubs = Inf,show.legend = TRUE)

### ---END following not working

View(module_hub_table)
#------- plot module hierarchy---------
colnames(module_table)[1] <- "id" # first column of module table must be labelled as "id".
hierarchy.obj <- plot_module_hierarchy(module.table = module_table,label.scaleFactor = 0.15,
                                        arrow.size = 0.03,node.label.color = "blue")

print(hierarchy.obj[[1]])

# module gene table
module_genes<-megena_res[["summary"]][["modules"]]

# to export nested list to text feile 
library(erer)
write.list(z = module_genes, file = "Megena_module_gene_table_all-cor.csv")

#-------END plot module hierarchy---------


output_cor_genes[["JMJD5"]]




names(megena_res)
View(megena_res$summary)

# Save an object to a file megena_res.rds
saveRDS(megena_res, file = "megena_res.rds")
# Restore the object
# megena_data<-readRDS(file = "megena_res.rds")

##### -------------- END MEGENA ----------


ddcor_all = ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                     compare = c("PrimaryTumor","SolidTissueNormal"),
                     adjust = "none", nPerm = 0, nPairs = 100)
head(ddcor_all)

#### Analysis _2 (specific gene JMJD5 correlation)
ddcor_res_JMJD5 = ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                     compare = c("PrimaryTumor","SolidTissueNormal"),
                     adjust = "perm", nPerm = 5, splitSet = "JMJD5")
View(ddcor_res_JMJD5)
write.csv(ddcor_res_JMJD5, "JMJD_correlation.csv")
ddcorGO_JMJD5 = ddcorGO(ddcor_res_JMJD5, universe = rownames(final_tum_norm_data), 
                             pval_GO_cutoff = 0.05, pval_gene_thresh = 0.05, ddcor_find_significant = TRUE,
                             gene_ontology = "all", HGNC_clean = TRUE, HGNC_switch = TRUE, annotation = "org.Hs.eg.db", calculateVariance = TRUE)

View(ddcorGO_JMJD5)
gain_gene<-ddcorGO_JMJD5[["significant_gain_of_correlation_genes"]]

gain_gene[gain_gene %in% get("GO:0016570", revmap(org.Hs.egGO))]
ls(org.Hs.egALIAS2EG)
keytypes(org.Hs.eg.db)

#####--------------- Mapping of our gene in Gene Ontology identified by "ddcorGO"----
# get the gene list GAIN of correlation GENE
gain_gene<-ddcorGO_JMJD5[["significant_gain_of_correlation_genes"]]

# to update the gene symbol in "gain_gene"
library(HGNChelper) # to check the gene name is old or updated
gain_gene_1<-checkGeneSymbols(gain_gene) # "geneList1" contain the list of gene, result show the new name
gain_gene_2<-gain_gene_1$Suggested.Symbol
gain_gene_2<-gain_gene_2[!is.na(gain_gene_2)] # remove "NA" gene from list
gain_gene_2 # final updated gene list


# convert the Gene Symbol to Entrez ID
Gene_information<-select(org.Hs.eg.db,
       keys = gain_gene_2,
       columns=c("ENTREZID","SYMBOL","GENENAME", "GO"),
       keytype="SYMBOL")

# get the uniqe entrez id
geneIds<-unique(Gene_information$ENTREZID)

# identify the entrez id in a GO term
interested_GO<-c("GO:0016570")
for (myGO in interested_GO){
  GOgene<-geneIds[geneIds %in% get(myGO, org.Hs.egGO2ALLEGS)]
  #GOgene<-geneIds[geneIds %in% get("GO:0016570", org.Hs.egGO2ALLEGS)]
   
  # convert the entrez id to GeneSymbol
  GOgene2<-Gene_information %>% filter(ENTREZID == GOgene) 
  GOgene2<-unique(GOgene2[,c(1,3)])

# unique gene information
GOgene3<-distinct(GOgene2, SYMBOL, .keep_all = TRUE)

output_file_2 <- paste("GO-",myGO,"-GainGene.csv",  sep="")
write.csv(GOgene3, output_file_2)

print(GOgene3)
}



table(unique(GOgene2)$SYMBOL)
gain_gene[gain_gene %in% org.Hs.egALIAS2EG]
org.Hs.egALIAS
probeSetSummary(GO)
Xn[genes %in% get("GO:0007399", revmap(org.Hs.eg.db))]
probeSetSummary(ddcorGO_JMJD5)

geneIds[geneIds %in% get("GO:0007338", revmap(org.Hs.egGO))]




 #### Analysis _4 (box plot)
plotVals(inputMat = final_tum_norm_data, design = design_matrix,
         compare = c("PrimaryTumor","SolidTissueNormal"), gene = "RAD54B") 

#### Analysis_5 (GO)
#BiocManager::install("MEGENA")
library(GOstats, quietly = TRUE)
library(HGNChelper, quietly = TRUE)
library(org.Hs.eg.db, quietly = TRUE)

ddcorGO_res_JMJD5 = ddcorGO(ddcor_res_JMJD5, universe = rownames(final_tum_norm_data), 
              gene_ontology = "all", HGNC_clean = TRUE, HGNC_switch = TRUE, annotation = "org.Hs.eg.db", calculateVariance = TRUE)

gain_gene<-ddcorGO_res_JMJD5[["significant_gain_of_correlation_genes"]]

geneIds[geneIds %in% get("GO:0007338", revmap(org.Hs.egGO))]
Xn[genes %in% get("GO:0007399", revmap(org.Hs.eg.db))]

str(ddcorGO_res)
View(ddcorGO_res_JMJD5)
ddcorGO_cg03101936 = ddcorGO(ddcor_res_cg03101936, universe = rownames(final_tum_norm_data), 
                             pval_GO_cutoff = 0.05, pval_gene_thresh = 0.05, ddcorGO_res =TRUE, ddcor_find_significant = TRUE,
                             gene_ontology = "all", HGNC_clean = TRUE, HGNC_switch = TRUE, annotation = "org.Hs.eg.db", calculateVariance = TRUE)

probeSetSummary(ddcorGO_cg03101936[["enrichment_significant_gain_of_correlation_genes"]][["BP"]])

ddcorGO_cg03101936
write.csv(ddcorGO_res[["significant_gain_of_correlation_genes"]], "significant_gain_of_correlation_genes.csv")
write.csv(ddcorGO_res[["significant_loss_of_correlation_genes"]], "significant_loss_of_correlation_genes.csv")

genes<-ddcorGO_res[["significant_gain_of_correlation_genes"]]
moduleGO(genes, universe = rownames(final_tum_norm_data), HGNC_clean = TRUE, HGNC_switch = TRUE,
         gene_ontology = "all", pval_GO_cutoff = 1, annotation = "org.Hs.eg.db",
         conditional = FALSE, calculateVariance = FALSE)

write.csv(ddcorGO_res[["enrichment_significant_gain_of_correlation_genes"]][["BP"]], "BP_gain_of_correlation_genes.csv")
write.csv(ddcorGO_res[["enrichment_significant_gain_of_correlation_genes"]][["CC"]], "CC_gain_of_correlation_genes.csv")
write.csv(ddcorGO_res[["enrichment_significant_gain_of_correlation_genes"]][["MF"]], "MF_gain_of_correlation_genes.csv")

write.csv(ddcorGO_res[["enrichment_significant_loss_of_correlation_genes"]][["BP"]], "BP_loss_of_correlation_genes.csv")
write.csv(ddcorGO_res[["enrichment_significant_loss_of_correlation_genes"]][["CC"]], "CC_loss_of_correlation_genes.csv")
write.csv(ddcorGO_res[["enrichment_significant_loss_of_correlation_genes"]][["MF"]], "MF_loss_of_correlation_genes.csv")

### Network integration
library(MEGENA)
?ddMEGENA
megena_res = ddMEGENA(ddcor_all, adjusted = FALSE, evalCompactness = TRUE)
megena_res<-ddMEGENA(ddcor_all, adjusted = TRUE, pval_gene_thresh = 0.05,
         evalCompactness = TRUE, nPerm = 100, hubPVal = 0.05,
         modulePVal = 0.05, minModSize = 10, maxModSize = 100,
         saveOutput = FALSE, parallelize = FALSE, nCores = 4)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


#------
library(gplots, quietly = TRUE)
darmanis_top =  filterGenes(final_tum_norm_data, 
                            filterTypes = c("central", "dispersion"), filterCentralPercentile = 0.50, 
                            filterDispersionPercentile = 0.50)
ddcor_res = ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                     compare = c("PrimaryTumor","SolidTissueNormal"),
                     #adjust = "none", heatmapPlot = TRUE, nPerm = 0, splitSet = "JMJD5")
                     adjust = "none", heatmapPlot = TRUE, nPerm = 0, nPairs = "all")


#-----
#### Analysis _A (specific for methylation probes correlation)
ddcor_res_cg02871891 = ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                                compare = c("PrimaryTumor","SolidTissueNormal"),
                           adjust = "perm", nPerm = 5, splitSet = "cg02871891")
View(ddcor_res_cg02871891)
ddcor_res_cg02871891<-ddcor_res_cg02871891  %>% filter(PrimaryTumor_pVal < 0.05) %>%
                                                filter(SolidTissueNormal_pVal < 0.05)
write.csv(ddcor_res_cg02871891, "cg02871891_correlation.csv")

plotCors(inputMat = final_tum_norm_data, design = design_matrix,
         compare = c("PrimaryTumor","SolidTissueNormal"),
         geneA = "cg02871891", geneB = "CSMD2")

plotVals(inputMat = final_tum_norm_data, design = design_matrix,
         compare = c("PrimaryTumor","SolidTissueNormal"), gene = "cg02871891")

###-------

ddcor_res_cg03101936 = ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                                compare = c("PrimaryTumor","SolidTissueNormal"),
                                adjust = "perm", nPerm = 5, splitSet = "cg03101936")

ddcor_res_ccg03101936<-ddcor_res_cg03101936  %>% filter(PrimaryTumor_pVal < 0.05) %>%
                                                filter(SolidTissueNormal_pVal < 0.05)
View(ddcor_res_cg03101936)
write.csv(ddcor_res_cg03101936, "cg03101936_correlation.csv")



dim(ddcor_res_JMJD5)

View(ddcor_res_cg02411582)
write.csv(ddcor_res_cg02411582, "cg02411582_correlation.csv")

###-------
ddcor_res_cg16752029 = ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                                compare = c("PrimaryTumor","SolidTissueNormal"),
                                adjust = "perm", nPerm = 5, splitSet = "cg16752029")
View(ddcor_res_cg16752029)
write.csv(ddcor_res_cg16752029, "cg16752029_correlation.csv")


###-------
ddcor_res_cg27118526 = ddcorAll(inputMat = final_tum_norm_data, design = design_matrix,
                                compare = c("PrimaryTumor","SolidTissueNormal"),
                                adjust = "perm", nPerm = 5, splitSet = "cg27118526")
View(ddcor_res_cg27118526)
write.csv(ddcor_res_cg27118526, "cg27118526_correlation.csv")

#### Analysis _3 (correlation plot)
# To plot the differential correlations between RTN4 and its top target, you can use this function:
plotCors(inputMat = final_tum_norm_data, design = design_matrix,
         compare = c("PrimaryTumor","SolidTissueNormal"),
         geneA = "BARX1", geneB = "cg27118526")

#### Analysis _4 (box plot)
plotVals(inputMat = final_tum_norm_data, design = design_matrix,
         compare = c("PrimaryTumor","SolidTissueNormal"), gene = "BARX1") 

ddcor_res_cg27118526  %>% filter(PrimaryTumor_pVal < 0.05) %>%
                          filter(SolidTissueNormal_pVal < 0.05)

dim(final_tum_norm_data)
