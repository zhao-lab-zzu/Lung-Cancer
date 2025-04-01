# install.packages("readxl")
# install.packages("dplyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"))
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages(c("WGCNA", "stringr", "reshape2"),repos=site)

library(readxl)
library(dplyr)
library(WGCNA)

setwd('~/project/WGCNA/')
trait_data <- read.csv("group_test.csv")
id_set <- trait_data[trait_data$Group != "Control", ]$Sample_ID

clean_trait_data_M <- trait_data[!is.na(trait_data$M), ]
clean_trait_data_N <- trait_data[!is.na(trait_data$N), ]
clean_trait_data_group <- trait_data[!is.na(trait_data$Group), ]
clean_trait_data_DMM <- trait_data[!is.na(trait_data$DMM), ]

data <- read.csv("NSCLCVSSCLCVSControl_diff_MS2_metabolites.csv")
datExpr <- data %>% select(starts_with("LC_"))
datExpr <- as.data.frame(datExpr)
datExpr[] <- lapply(datExpr, as.numeric)
gene_matrix <- t(datExpr)
colnames(gene_matrix) <- data[["MS2.name"]]

CLUSTER_ALL = TRUE
# gene_matrix_cleaned 
# clean_trait_data_M 
if (CLUSTER_ALL) {
  gene_matrix_cleaned <- gene_matrix[id_set, ]
  clean_trait_data_M <- clean_trait_data_M[clean_trait_data_M$Sample_ID %in% rownames(gene_matrix_cleaned), ]
  clean_trait_data_N <- clean_trait_data_N[clean_trait_data_N$Sample_ID %in% rownames(gene_matrix_cleaned), ]
  clean_trait_data_group <- clean_trait_data_group[clean_trait_data_group$Sample_ID %in% rownames(gene_matrix_cleaned), ]
  clean_trait_data_DMM <- clean_trait_data_DMM[clean_trait_data_DMM$Sample_ID %in% rownames(gene_matrix_cleaned), ]
} else {
  ##
  stop("Not implemented")
  common_samples <- intersect(rownames(gene_matrix), clean_trait_data$Sample_ID)
  gene_matrix_cleaned <- gene_matrix[common_samples, ]
  clean_trait_data <- clean_trait_data[clean_trait_data$Sample_ID %in% common_samples, ]
  clean_trait_data <- clean_trait_data[match(common_samples, clean_trait_data$Sample_ID), ]  
}

compute_MEs <- function(gene_matrix_cleaned) {
  gsg <- goodSamplesGenes(gene_matrix_cleaned, verbose = 3)
  if (!gsg$allOK) {
    # Remove genes or samples that don't pass the filter
    gene_matrix_cleaned <- gene_matrix_cleaned[gsg$goodSamples, gsg$goodGenes]
  }
  
  # Choose a set of soft-thresholding powers to test
  powers <- c(1:20)
  
  # Pick the best soft-thresholding power
  sft <- pickSoftThreshold(gene_matrix_cleaned, powerVector = powers, verbose = 5)
  
  # Plot the results
  sizeGrWindow(9, 5)
  par(mfrow = c(1, 2))
  
  # Scale-free fit index (R^2) as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, cex = 0.9, col = "red")
  abline(h = 0.90, col = "red")  # Line indicating an R^2 cut-off
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", 
       main = "Mean connectivity")
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")
  
  sft$powerEstimate
  
  
  # Build the network using blockwiseModules function
  net <- blockwiseModules(gene_matrix_cleaned, power = sft$powerEstimate,
                          TOMType = "unsigned", minModuleSize = 30,
                          reassignThreshold = 0, mergeCutHeight = 0.25,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE, saveTOMFileBase = "TOM",
                          verbose = 3)
  
  # Plot the dendrogram of genes and module colors
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors", dendroLabels = FALSE, hang = 0.03)
  
  # Check the dimensions to ensure alignment
  # dim(gene_matrix_cleaned)
  # dim(clean_trait_data)
  
  # MEs <- net$MEs
  # head(MEs)[1:5, 1:5]
  # geneTree <- net$dendrograms[[1]]
  
  
  # Calculate module eigengenes
  MEs0 <- moduleEigengenes(gene_matrix_cleaned, colors = moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  # MM <- cor(gene_matrix_cleaned, MEs, use = "p")
  # MMPvalue <- corPvalueStudent(MM, nSamples = nrow(gene_matrix_cleaned))
  # # Plot histogram of Module Membership (MM) values
  # hist(abs(MM), breaks = 30, main = "Distribution of Module Membership", xlab = "Module Membership")
  # 
  # hub_genes <- which(abs(MM) > 0.7)
  # hub_genes
  
  # Assume 'moduleColor' contains the color assignments for each gene, e.g., moduleColor = c("blue", "brown", "yellow", ...)
  # 'datExpr' is the gene expression data matrix.
  
  # Loop through each unique module color
  for (module in names(table(moduleColors))) {
    
    # Select genes that belong to the current module
    moduleGenes = which(moduleColors == module)
    
    # Extract expression data for the genes in the current module
    moduleExpr = gene_matrix_cleaned[, moduleGenes]
    
    # Calculate the module eigengene
    moduleEigengene = MEs[, paste0("ME", module)]
    
    # Calculate module membership (kME) for the current module
    kME = cor(moduleExpr, moduleEigengene, use = "p")
    
    # Calculate intramodular connectivity (kIN)
    adjMat = adjacency(moduleExpr, power = sft$powerEstimate)  # Adjust power as needed
    kIN = rowSums(adjMat)
    
    # Print or save results for the current module
    cat("Module:", module, "\n")
    cat("Number of genes:", length(moduleGenes), "\n")
    
    # Optional: Filter hub genes and print/save them
    # hubGenes = data.frame(Gene = colnames(gene_matrix_cleaned)[moduleGenes],
    # ModuleMembership = kME,
    # IntramodularConnectivity = kIN)
    hubGenes = data.frame(ModuleMembership = kME,
                          IntramodularConnectivity = kIN)
    hubGenes = hubGenes[hubGenes$ModuleMembership > 0.7 & hubGenes$IntramodularConnectivity > quantile(kIN, 0.90), ]
    
    print(head(hubGenes))  # Show hub genes for the current module
  }
  
  return(MEs)
}

# table(gene_matrix_cleaned)
MEs <- compute_MEs(gene_matrix_cleaned)
write.csv(MEs, "MEs.csv")




# Merge all 'M1*' values into a single category 'M1'
clean_trait_data_M$M <- ifelse(grepl("^M0", clean_trait_data_M$M), "0", clean_trait_data_M$M)
clean_trait_data_M$M <- ifelse(grepl("^M1", clean_trait_data_M$M), "1", clean_trait_data_M$M)

clean_trait_data_N$N <- ifelse(grepl("^N2", clean_trait_data_N$N), "1_3", clean_trait_data_N$N)
clean_trait_data_N$N <- ifelse(grepl("^N3", clean_trait_data_N$N), "1_3", clean_trait_data_N$N)
clean_trait_data_N$N <- ifelse(grepl("^N0", clean_trait_data_N$N), "0", clean_trait_data_N$N)
clean_trait_data_N$N <- ifelse(grepl("^N1", clean_trait_data_N$N), "1_3", clean_trait_data_N$N)


# Check if the merging worked correctly
table(clean_trait_data_M$M)
table(clean_trait_data_N$N)
table(clean_trait_data_DMM$DMM)
table(clean_trait_data_group$Group)

# Assume that the 'M' column is a factor (categorical variable)
clean_trait_data_M$M <- as.factor(clean_trait_data_M$M)
clean_trait_data_N$N <- as.factor(clean_trait_data_N$N)
clean_trait_data_DMM$DMM <- as.factor(clean_trait_data_DMM$DMM)
clean_trait_data_group$Group <- as.factor(clean_trait_data_group$Group)


# Plot the module-trait relationships
sizeGrWindow(18, 6)
par(mar = c(5, 8, 4, 2)) 

MEs1 <- MEs[clean_trait_data_M$Sample_ID, ]
# all(clean_trait_data$Sample_ID == rownames(MEs))
# Use model.matrix to perform one-hot encoding on the M column
one_hot <- model.matrix(~M-1, data = clean_trait_data_M)
# View the one-hot encoded matrix
head(one_hot)
# Correlate module eigengenes with the trait data
moduleTraitCor <- cor(MEs1, one_hot, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(one_hot))
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
combined_matrix <- moduleTraitCor
combined_x <- colnames(one_hot)
combined_text <- textMatrix

# labeledHeatmap(Matrix = moduleTraitCor, 
#                xLabels = colnames(one_hot), 
#                yLabels = names(MEs1), 
#                ySymbols = names(MEs1),
#                colorLabels = FALSE, 
#                colors = blueWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE, 
#                cex.text = 0.5,
#                zlim = c(-1, 1))

MEs1 <- MEs[clean_trait_data_N$Sample_ID, ]
# all(clean_trait_data$Sample_ID == rownames(MEs))
# Use model.matrix to perform one-hot encoding on the M column
one_hot <- model.matrix(~N-1, data = clean_trait_data_N)
# View the one-hot encoded matrix
head(one_hot)
# Correlate module eigengenes with the trait data
moduleTraitCor <- cor(MEs1, one_hot, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(one_hot))
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

combined_matrix <- cbind(combined_matrix, moduleTraitCor)
combined_x <- c(combined_x, colnames(one_hot))
combined_text <- cbind(combined_text, textMatrix)


MEs1 <- MEs[clean_trait_data_DMM$Sample_ID, ]
# all(clean_trait_data$Sample_ID == rownames(MEs))
# Use model.matrix to perform one-hot encoding on the M column
one_hot <- model.matrix(~DMM-1, data = clean_trait_data_DMM)
# View the one-hot encoded matrix
head(one_hot)
# Correlate module eigengenes with the trait data
moduleTraitCor <- cor(MEs1, one_hot, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(one_hot))
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
combined_matrix <- cbind(combined_matrix, moduleTraitCor)
combined_x <- c(combined_x, colnames(one_hot))
combined_text <- cbind(combined_text, textMatrix)

labeledHeatmap(Matrix = combined_matrix, 
               xLabels = combined_x, 
               yLabels = names(MEs1), 
               ySymbols = names(MEs1),
               colorLabels = FALSE, 
               colors = blueWhiteRed(50),
               textMatrix = combined_text,
               setStdMargins = FALSE, 
               cex.text = 0.5,
               zlim = c(-1, 1))


################ 菌属丰度 ##############
abundance_data <- read.csv('relative_abundance_with_genus.csv',row.names = 1)
abundance_data <- abundance_data[!grepl("unclassified", rownames(abundance_data), ignore.case = TRUE), ]
abundance_data$Taxonomy <- rownames(abundance_data)
taxonomy_data <- abundance_data %>%
  select(Taxonomy) %>%
  separate(Taxonomy,into = c('domain','men','gang','mu','ke','shu'),sep = ';',extra = "drop") %>%
  select(shu,men)
taxonomy_data$species <- rownames(taxonomy_data)

trans_abundance <- abundance_data %>%
  select(-Taxonomy) %>%
  t() %>%
  as.data.frame()

mesamples <- rownames(MEs)
tasamples <- rownames(trans_abundance)


new_abu <- trans_abundance[match(mesamples, rownames(trans_abundance)),]


modCor = cor(MEs,new_abu,use = 'p',method = 'pearson')
modP = corPvalueStudent(modCor,length(mesamples))


threshold_cor <- 0.3
threshold_p <- 0.05

significant_pairs <- (abs(modCor) > threshold_cor) & (modP < threshold_p)


filtered_modCor <- modCor
filtered_modCor[!significant_pairs] <- NA

index <- which(!is.na(filtered_modCor), arr.ind = TRUE)
col_index <- as.data.frame(index)$col

my_modCor <- modCor[,col_index]
my_modP <- modP[,col_index]

col_names <- colnames(my_modCor)


genus_names <- gsub(".*g__([^;]*);.*", "\\1", col_names)
genus_names <- gsub(".*g__([^;]*)$", "\\1", genus_names)

print(genus_names)
colnames(my_modCor) <- genus_names

rdata <- read.csv('Figure4/Figure4c_r.csv',row.names = 1)
pdata <- read.csv('Figure4/Figure4c_p.csv',row.names = 1)

dmm_r <- rdata %>%
  select(DMM.Cluster1,DMM.Cluster2)
dmm_p <- pdata %>%
  select(DMM.Cluster1,DMM.Cluster2)
my_modCor <- merge(my_modCor,dmm_r,by='row.names',all = TRUE)
rownames(my_modCor)<- my_modCor$Row.names
my_modCor <- my_modCor %>%
  select(-Row.names)

my_modP <- merge(my_modP,dmm_p,by='row.names',all= TRUE)
rownames(my_modP) <- my_modP$Row.names
my_modP <- my_modP %>%
  select(-Row.names)

display_numbers <- matrix(
  paste0(sprintf("%.2f", as.matrix(my_modCor)), "\n(", sprintf("%.2g", as.matrix(my_modP)), ")"),
  nrow = nrow(my_modCor)
)


annotation_colors <- data.frame(
  Module = factor(rownames(my_modCor)),
  row.names = rownames(my_modCor)
)


module_colors <- c(
  MEgreen = "green",
  MEred = "red",
  MEblack = "black",
  MEyellow = "yellow",
  MEmagenta = "magenta",
  MEblue = "blue",
  MEturquoise = "turquoise",
  MEbrown = "brown",
  MEpink = "pink",
  MEgrey = "grey"
)


ann_colors <- list(
  Module = module_colors
)
setwd('~/project/sun/new_fig/')
p1 <- ComplexHeatmap::pheatmap(
  as.matrix(my_modCor),

  display_numbers = display_numbers,
  number_format = "%.2f",
  #color = colorRampPalette(c("#7b68ee", "white", "#ff6347"))(100),
  color = colorRampPalette(c("#0d8cff", "white", "#ff3200"))(100), # "#ff3200", "#0d8cff", 
  breaks = seq(-1, 1, length.out = 101), 
  annotation_row = annotation_colors,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize = 10,
  fontsize_number = 10,
  row_names_side = 'left',
  border_color = NA,  
  annotation_names_row = FALSE,  
  annotation_legend = FALSE,    
  #filename = 'taxo_wgcna.pdf',
  fontfamily = "Arial", fontface = 'bold',
  name = NULL,
  #width = 10,height = 8
)
{
pdf('taxo_wgcna.pdf',width = 10,height = 6)
p1
dev.off()
}



########################## 表征 ###################
sample_groups <- read.csv('group_testnew.csv',row.names = 1)
biao_data <- sample_groups %>%
  select(Gender,Age,Ab,BMI,Smoking_history)
biao_data <- biao_data[match(rownames(MEs), rownames(biao_data)),] 



data1 <- biao_data %>%
  select(Age,BMI)
age_bmi <- cor(MEs,data1,use='p',method = "spearman")
age_bmi_p <- corPvalueStudent(age_bmi,length(rownames(MEs)))

data2 <- biao_data %>%
  select(Gender,Ab,Smoking_history)


data2$Gender <- as.numeric(data2$Gender == "M")
data2$Smoking_history <- as.numeric(data2$Smoking_history == "yes")
data2$Ab <- as.numeric(data2$Ab == "yes") 


data2_cor <- bicor(MEs, data2, use = "pairwise.complete.obs")
data2_p <- corPvalueStudent(data2_cor,length(rownames(MEs)))




t_data <- sample_groups %>%
  select(T) %>%
  filter(T != "")
t_data$T <- as.numeric(t_data$T == "T1"|t_data$T == 'T2')


merged_t_data <- merge(MEs, t_data, by = "row.names", all = FALSE)
rownames(merged_t_data) <- merged_t_data$Row.names
merged_t_data$Row.names <- NULL


t_cor <- bicor(merged_t_data%>%select(-T), merged_t_data%>%select(T), use = "pairwise.complete.obs")
t_p <- corPvalueStudent(t_cor,length(rownames(merged_t_data)))


m_data <- sample_groups %>%
  select(M) %>%
  filter(M != "")
m_data$M <- as.numeric(m_data$M == "M1")


merged_m_data <- merge(MEs, m_data, by = "row.names", all = FALSE)
rownames(merged_m_data) <- merged_m_data$Row.names
merged_m_data$Row.names <- NULL


m_cor <- bicor(merged_m_data%>%select(-M), merged_m_data%>%select(M), use = "pairwise.complete.obs")
m_p <- corPvalueStudent(m_cor,length(rownames(merged_m_data)))


n_data <- sample_groups %>%
  select(N) %>%
  filter(N != "")
n_data$N <- as.numeric(n_data$N != "N0")


merged_n_data <- merge(MEs, n_data, by = "row.names", all = FALSE)
rownames(merged_n_data) <- merged_n_data$Row.names
merged_n_data$Row.names <- NULL


n_cor <- bicor(merged_n_data%>%select(-N), merged_n_data%>%select(N), use = "pairwise.complete.obs")
n_p <- corPvalueStudent(n_cor,length(rownames(merged_n_data)))


cor_data <- cbind(age_bmi,data2_cor,t_cor,n_cor,m_cor)
p_data <- cbind(age_bmi_p,data2_p,t_p,n_p,m_p)



display_numbers <- matrix(
  paste0(sprintf("%.2f", as.matrix(cor_data)), "\n(", sprintf("%.2g", as.matrix(p_data)), ")"),
  nrow = nrow(cor_data)
)


annotation_colors <- data.frame(
  Module = factor(rownames(cor_data)),
  row.names = rownames(cor_data)
)


module_colors <- c(
  MEgreen = "green",
  MEred = "red",
  MEblack = "black",
  MEyellow = "yellow",
  MEmagenta = "magenta",
  MEblue = "blue",
  MEturquoise = "turquoise",
  MEbrown = "brown",
  MEpink = "pink",
  MEgrey = "grey"
)


ann_colors <- list(
  Module = module_colors
)

p1 <- ComplexHeatmap::pheatmap(
  as.matrix(cor_data),

  display_numbers = display_numbers,
  #color = colorRampPalette(c("#7b68ee", "white", "#ff6347"))(100),
  color = colorRampPalette(c("#0d8cff", "white", "#ff3200"))(100), # "#ff3200", "#0d8cff", 
  breaks = seq(-1, 1, length.out = 101), 
  annotation_row = annotation_colors,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  ##
  show_rownames = TRUE,
  cellwidth = 40,  
  cellheight = 22, 
  fontsize = 10,
  fontsize_number = 10,
  row_names_side = 'left',
  border_color = NA,  
  annotation_names_row = FALSE,  
  annotation_legend = FALSE,    
  #filename = 'taxo_wgcna.pdf',
  
  fontfamily = "Arial", fontface = 'bold',
  name = NULL,
  #width = 10,height = 8
)
{
  pdf('wgcna_new.pdf',width = 8,height = 6)
  p1
  dev.off()
}
