# Trajectory analysis with slingshot

library (slingshot)
library (SingleCellExperiment)
library (zinbwave)
library (batchelor)
library (uwot)
library (mclust)
library (RColorBrewer)
library (magrittr)
library (Cairo)
#BiocManager::install("")
#devtools.install_version("ggplot2", version = "0.9.1"
#packageurl <- "https://bioconductor.org/packages/3.12/bioc/src/contrib/Archive/BiocParallel/BiocParallel_1.24.0.tar.gz"
#install.packages(packageurl, repos=NULL, type="source")

# import data
gene_data = t (read.table ('/data/liver_csc/data/merged_matrix_7_trim.tsv', header = TRUE, row.names = 1, sep = "\t", as.is = TRUE))
sample_data = read.table ('/data/liver_csc/data/sample_set_7_trim.tsv', header = TRUE, row.names = 1, sep = "\t", as.is = TRUE)
sample_factors = factor (sample_data[1])

# Create SingleCellExperiment object
sce_data = SingleCellExperiment(assays = List(counts = gene_data), colData=sample_data)

# ZINB-WaVE
# filter out low coverage genes, at least 5 reads in at least 5 samples
filter <- rowSums(assay(sce_data)>5)>5
sce_data <- sce_data[filter,]

# keep the top 1000 most variable genes
assay(sce_data) %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sce_data)
vars <- sort(vars, decreasing = TRUE)
sce_data <- sce_data[names(vars)[1:1000],]

# compute low-dimensional representation with 2 factors
sce_data_zinb <- zinbwave(sce_data, K = 2, X="~study", epsilon=1000)

#reducedDims(sce_data) <- SimpleList(zinb = sce_data_zinb)

# use Mclust to cluster data
cl1 = Mclust(reducedDims(sce_data_zinb)$zinbwave)$classification
colData(sce_data_zinb)$GMM = cl1

# use slingshot
sce_data_zinb <- slingshot(sce_data_zinb, clusterLabels = 'GMM', reducedDim = 'zinbwave')

summary(sce_data_zinb$slingPseudotime_1)

# Save/load Rdata
save(sce_data_zinb, file="/data/liver_csc/data/sce_data_zinb_csc_adult_prot_1.Rdata")
#load("sce_data_zinb.Rdata")

# Plot pseudotime
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce_data_zinb$slingPseudotime_1, breaks=100)]
#mycolors <- colorRampPalette(brewer.pal(8, 'Paired'))(15)
mycolors <- colorRampPalette(brewer.pal(8, 'Set1'))(3)

pdf(file = "/home/borzik01/liver_csc/data/trajectory_csc_adult_prot_1.pdf")
plot(reducedDims(sce_data_zinb)$zinbwave, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce_data_zinb), lwd=2, col='black')
dev.off()

pdf(file = "/home/borzik01/liver_csc/data/trajectory_csc_adult_prot_2.pdf")
plot(reducedDims(sce_data_zinb)$zinbwave, col = brewer.pal(9,'Set1')[sce_data_zinb$class], pch=16, asp = 1)
lines(SlingshotDataSet(sce_data_zinb), lwd=2, type = 'lineages', col = 'black')
dev.off()

#pdf(file = "/home/borzik01/liver_csc/data/trajectory_class.pdf")
Cairo::Cairo(
  30, #length
  15, #width
  file = paste("/home/borzik01/liver_csc/data/trajectory csc adult prot class 1", ".png", sep = ""),
  type = "png", #tiff
  bg = "white", #white or transparent depending on your requirement 
  dpi = 600,
  units = "cm" #you can change to pixels etc 
)
par(xpd = T, mar = par()$mar + c(0,0,0,20))
plot(reducedDims(sce_data_zinb)$zinbwave, col = mycolors[factor(sce_data_zinb$class)], pch=16, asp = 1, cex=.5)
lines(SlingshotDataSet(sce_data_zinb), lwd=1, type = 'lineages', col = 'black', cex=1)
legend(5, 2, legend = levels(factor(sce_data_zinb$class)), col = mycolors[seq_along(levels(factor(sce_data_zinb$class)))], pch=16)
#par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

pdf(file = "/home/borzik01/liver_csc/data/trajectory_csc_adult_prot_study.pdf")
plot(reducedDims(sce_data_zinb)$zinbwave, col = factor(sce_data_zinb$study), pch=16, asp = 1)
lines(SlingshotDataSet(sce_data_zinb), lwd=2, type = 'lineages', col = 'black')
dev.off()
