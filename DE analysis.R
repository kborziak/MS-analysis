#EdgeR for prostate data analysis

library (edgeR)
library (foreach)

# import data
gene_data = read.table ('/data/prostate/data_files/prostate_counts_matrix_3', header = TRUE, row.names = 1, sep = "\t", as.is = TRUE)
sample_data = read.table ('/home/borzik01/prostate_cancer/data/prostate_sample_set_3', header = TRUE, row.names = 1, sep = "\t", as.is = TRUE)

#EdgeR

# cell types
group <- factor(sample_data$comp1)

# group design for contrasts
#cdr <- scale(colMeans(gene_data > 0))
cdr <- scale(colMeans(y$counts > 0))

#design <- model.matrix(~0+group)
design <- model.matrix(~cdr+group)

colnames (design) <- levels (group)
#design <- model.matrix(~0+group + sample_data$study)
#colnames (design) <- c(levels (group), "study")
#design <- model.matrix(~0+group, data=y$samples)

# create DGEList object
y <- DGEList(counts=gene_data,group=group)

# free space
rm (gene_data)

# calculate TMM normalization
y <- calcNormFactors(y)

# estimate the NB dispersion 
y <- estimateDisp(y, design, robust=TRUE)

# estimate values of the GLM coefficients for each gene
fit <- glmQLFit(y, design, robust=TRUE)
head (fit$coefficients)

# Output fit coefficients
#write.table (fit$coefficients, file = sprintf ("/home/borzik01/prostate_cancer/data/fit_coefficients_1", qlf$comparison), row.names = TRUE, sep = "\t")

# carry out the QL F-test
con <- makeContrasts((cancer_ep) - (bph_ep), levels=design)

qlf <- glmQLFTest(fit, contrast=con)

# adjust p-values
qlf$table$bonferroni = p.adjust(qlf$table$PValue, method="bonferroni")

# Output QL F-test results
write.table (qlf$table, file = sprintf ("/home/borzik01/prostate_cancer/data/DE_cancer_bph_ep_3", qlf$comparison), row.names = TRUE, sep = "\t")

# Calculate fitted log CPM values
logCPM.fit = cpm (y, normalized.lib.sizes = TRUE, log = TRUE)
logCPM.fit_round = round(logCPM.fit, 3)

# Output logCPM
write.table (logCPM.fit_round, file = sprintf ("/home/borzik01/prostate_cancer/data/logCPM_round_2"), row.names = TRUE, sep = "\t")
