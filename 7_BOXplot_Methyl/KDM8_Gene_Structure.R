library(ggplot2)

# Gene structure data
gene_structure <- data.frame(
  seqnames = c("chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16", "chr16"),
  start = c(27214422, 27214687, 27214750, 27214781, 27214832, 27215034, 27215225, 27215246, 27215281, 27215306, 27215372, 27217587, 27219548, 27227277, 27228679, 27230286),
  end = c(27214422, 27214687, 27214750, 27214781, 27214832, 27215034, 27215225, 27215246, 27215281, 27215306, 27215372, 27217587, 27219548, 27227277, 27228679, 27230286),
  width = rep(1, 16),
  probeID = c("cg10221365", "cg24206694", "cg16752029", "cg08572679", "cg06329197", "cg10214171", "cg09140090", "cg01755539", "cg00636124", "cg24705286", "cg27118526", "cg00632019", "cg09632858", "cg03101936", "cg02871891", "cg07361448"),
  Gene_Symbol = rep("KDM8", 16)
)

# Methylation probe data
methylation_probes <- data.frame(
  probeID = c("cg10221365", "cg24206694", "cg16752029", "cg08572679", "cg06329197", "cg10214171", "cg09140090", "cg01755539", "cg00636124", "cg24705286", "cg27118526", "cg00632019", "cg09632858", "cg03101936", "cg02871891", "cg07361448"),
  methylation = runif(16)
)

# Create a ggplot object for gene structure visualization
gene_plot <- ggplot(gene_structure) +
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1), fill = "lightblue") +
  scale_x_continuous(limits = c(27214000, 27233000)) +
  labs(x = "Genomic Position", y = "Methylation Probe") +
  theme_bw()

# Add the methylation probe positions as points on the gene structure plot
gene_plot_with_methylation <- gene_plot +
  geom_point(data = methylation_probes, aes(x = start, y = 1, color = methylation), size = 2) +
  scale_color_gradient(low = "blue", high = "red")

# Display the plot
gene_plot_with_methylation








############ END #########
# Load the required packages
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)


# Load the gene annotation data
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Get gene structure information for JMJD5 gene
# Gene id for KDM8 is "79831"
GeneID <- "79831"

gene <- genes(txdb, filter=list(gene_id=GeneID))
transcripts(txdb, filter=list(gene_id=GeneID))
exons <- exonsBy(txdb, filter=list(gene_id=GeneID))
introns <- intronsByTranscript(txdb, filter=list(gene_id=GeneID))
promoters <- promoters(txdb, filter=list(gene_id=GeneID))

# gene id for KDM8 is "79831"
GeneID<-"79831"
gene<-genes(txdb, columns="gene_id", filter=list(gene_id=GeneID), single.strand.genes.only=TRUE)

gene <- genes(txdb, filter = GenomicFeatures::geneSymbol == gene_symbol)
exons <- exonsBy(txdb, filter = GenomicFeatures::geneSymbol == gene_symbol)
introns <- intronsByTranscript(txdb, filter = GenomicFeatures::geneSymbol == gene_symbol)
promoters <- promoters(txdb, upstream = 2000, downstream = 200, filter = GenomicFeatures::geneSymbol == gene_symbol)

# Access the desired information from the gene structure objects
gene_info <- as.data.frame(gene)
exon_ranges <- ranges(exons)
intron_ranges <- ranges(introns)
promoter_ranges <- ranges(promoters)

# Print the extracted information
print(gene_info)
print(exon_ranges)
print(intron_ranges)
print(promoter_ranges)




====


# Load the required package
library(GenomicFeatures)

# Load the gene annotation data
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene





# Get gene structure information
genes <- genes(txdb)
exons <- exons(txdb)
introns <- intronsByTranscript(txdb)
promoters <- promoters(txdb, upstream = 2000, downstream = 200)

# Access the desired information from the gene structure objects
# For example, to get the ranges of exons:
exon_ranges <- ranges(exons)

keytypes(txdb)
keys(txdb)
columns(txdb)

#genes(txdb, columns="gene_id", filter=NULL, single.strand.genes.only=TRUE)

# gene id for KDM8 is "79831"
GeneID<-"79831"
gene<-genes(txdb, columns="gene_id", filter=list(gene_id=GeneID), single.strand.genes.only=TRUE)

# Access the desired information from the gene structure objects
gene_info <- as.data.frame(gene)





gene_symbol <- "JMJD5"

gene_info <- select(txdb, keys = GeneID, keytype = "GENEID", 
                    columns = c("TXSTART", "TXEND", "TXCHROM"))

exons <- exonsBy(txdb, gene = gene_info)
introns <- intronsByTranscript(txdb, gene = gene_info)
promoters <- promoters(txdb, gene = gene_info)

library(ggbio)
library(ggbio)

# Create a GRanges object with the gene structure information
gene_structure <- GRanges(
  seqnames = c("chr16"),
  ranges = IRanges(start = c(27214422, 27214687, 27214750, 27214781, 27214832, 27215034, 27215225, 27215246, 27215281, 27215306, 27215372, 27217587, 27219548, 27227277, 27228679, 27230286, 27232615),
                   end = c(27214422, 27214687, 27214750, 27214781, 27214832, 27215034, 27215225, 27215246, 27215281, 27215306, 27215372, 27217587, 27219548, 27227277, 27228679, 27230286, 27232615))
)

# Create a GRanges object with the methylation probe information
methylation_probes <- GRanges(
  seqnames = c("chr16"),
  ranges = IRanges(start = c(27214422, 27214687, 27214750, 27214781, 27214832, 27215034, 27215225, 27215246, 27215281, 27215306, 27215372, 27217587, 27219548, 27227277, 27228679, 27230286, 27232615),
                   end = c(27214422, 27214687, 27214750, 27214781, 27214832, 27215034, 27215225, 27215246, 27215281, 27215306, 27215372, 27217587, 27219548, 27227277, 27228679, 27230286, 27232615))
)

# Create a ggplot object for gene structure visualization
gene_plot <- autoplot(gene_structure, geom = "rect", fill = "lightblue") +
  scale_x_continuous(limits = c(27214000, 27233000)) +
  labs(x = "Genomic Position", y = "Methylation Probe") +
  theme_bw()

# Add the methylation probe positions as points on the gene structure plot
gene_plot_with_methylation <- gene_plot +
  geom_point(data = methylation_probes, aes(x = start, y = 1), size = 2, color = "red")

# Display the plot
gene_plot_with_methylation
