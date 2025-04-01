setwd('~/project/sun/new_fig/')
group <- read.csv('Figure2a-e/Figure2h_group_DMM.csv')
rawdata <- read.csv('Figure2a-e/Figure2h_ASV_DMM.csv',row.names = 1)
data <- t(rawdata[,-ncol(rawdata)])
data = as.data.frame(data)
data = data %>% 
  mutate(across(everything(), as.numeric))
data_cleaned <- data %>%
  mutate(across(everything(), ~ifelse(is.infinite(.), NA, .))) %>%  # 将无限值转为NA
 mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))  # 用均值替换NA
pca_result <- prcomp(data, scale. = TRUE)
summary(pca_result)
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Species = group$DMM
)

var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)


ggplot(pca_data, aes(x = PC1, y = PC2, color = Species)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "PCA Scatter Plot",
    x = paste0("PC1 (", pc1_var, "% variance explained)"),
    y = paste0("PC2 (", pc2_var, "% variance explained)")
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  stat_ellipse(type = "norm", level = 0.95) 

print("Principal Component Loadings:")
print(pca_result$rotation[,1:2])
print(paste("PC1 variance explained:", pc1_var, "%"))
print(paste("PC2 variance explained:", pc2_var, "%"))

########################### PLSDA ######################

library(mixOmics)  
library(ggplot2)   
library(dplyr)    

# calculate PLS-DA scores
calculate_plsda_scores <- function(X, Y, ncomp = 2, scale = TRUE) {
  
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X should be a matrix or data frame")
  }
  
  X <- as.matrix(X)
  Y <- as.factor(Y)
  plsda_model <- plsda(X, Y, 
                       ncomp = ncomp,
                       scale = scale)
  
  scores <- plsda_model$variates$X
  
  explained_variance <- plsda_model$explained_variance$X
  
  loadings <- plsda_model$loadings$X

  results <- list(
    scores = scores,
    loadings = loadings,
    explained_variance = explained_variance,
    model = plsda_model
  )
  
  return(results)
}




########################## F4a ####################
setwd('~/project/sun/new_fig/')
rawdata <- read.csv('Figure4/Figure4a_PLSDA.csv',row.names = 1)
data <- t(rawdata)
group_data <- read.csv('Figure4/Figure4a_group.csv')

# PERMANOVA
library(vegan)
dist_matrix <- vegdist(data)
perm_test <- adonis2(dist_matrix ~ group_data$Group, permutations = 999)
p_value <- perm_test$`Pr(>F)`[1]
results <- calculate_plsda_scores(data, group_data$Group)
scores <- as.data.frame(results$scores)
scores$Group <- group_data$Group
model <- results$model
ggplot(scores, aes(x = `comp 1`, y = `comp 2`, color = Group)) +
  geom_point(size = 2, alpha = 0.8,shape = 16) +
  stat_ellipse(level = 0.95, type = "t") +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  labs(
    title = "",
    x = paste0("PC1 (", round(model$explained_variance$X[1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(model$explained_variance$X[2] * 100, 2), "%)"),
    color = "Group"
  ) +
  scale_color_manual(values = c("#4169E1", "#ff7070",'#00C9A7')) +
annotate("text", 
         x = 40, 
         y = -10,
         label = paste("PERMANOVA\np =", formatC(p_value, digits = 3, format = "e")),
         hjust = 0,
         vjust = 1)

ggsave('F4a.pdf',width = 6,height = 4)

###################### F5a ######################
setwd('~/project/sun/new_fig/')
rawdata <- read.csv('Figure5/Figure5a_PCA.csv',row.names = 1)
data <- t(rawdata)
group_data <- read.csv('Figure5/Figure5a_Stage.csv')
group_data <- group_data[match(rownames(data),group_data$Sample_ID),]
  
pca_result <- prcomp(data, scale. = TRUE)
summary(pca_result)
pca_data <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  Stage = group_data$Stage
) %>%
  filter(!is.na(Stage))

var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)

mycolor <- c("#ffc3c3", "#a3a7de", "#4dbbd5","#cd5c5b", "#f0d491")
ggplot(pca_data, aes(x = PC1, y = PC2, color = Stage, fill = Stage)) +
  geom_point(size = 3, alpha = 0.9,shape = 16) +
  theme_classic() +
  labs(
    title = "",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  scale_fill_manual(values = mycolor) +
  scale_color_manual(values = mycolor) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    text = element_text(family = 'Arial',face = 'bold')
  ) +
  stat_ellipse(type = "norm", level = 0.95) 
ggsave('F5a.pdf',width = 4,height = 3)


##################### Fs1c ################################
setwd('~/project/sun/new_fig/补充图/')
library(ropls)
rawdata <- read.csv('S1/S1c_PCoA.csv',row.names = 1)
data <- as.data.frame(t(rawdata))
group_data <- read.csv('S1/S1c_group.csv',row.names = 1)


X <- as.matrix(data)
Y <- factor(group_data$Group)  

plsda_model <- opls(X, Y, scaleC = 'pareto', predI = 3, permI = 200)



########################## Fs2c ######################
setwd('~/project/sun/new_fig/补充图/')
rawdata <- read.csv('S2/S2c_ASV.csv',row.names = 1)
data <- as.data.frame(t(rawdata))
group_data <- read.csv('S2/S2c_group.csv')

df_group <- group_data %>%
  filter(Sample_ID %in% rownames(data))
df_group <- df_group[match(rownames(data),df_group$Sample_ID),]
pcoa_solution <- function(file1,group_file,sample_col,group_col,colors0){
  library(vegan)
  rawdata <- read.csv(file1,row.names = 1)
  data <- as.data.frame(t(rawdata))
  group_data <- read.csv(group_file)
  # comom_samples <- intersect(rownames(data),group_data[[sample_col]])
  # data <- data %>%
  #   filter(row.names(.) %in% comom_samples)
  # df_group <- group_data %>%
  #   filter(sample_col %in% rownames(data))
  
  df_group <- group_data[match(rownames(data),group_data[[sample_col]]),]
  ## calulate distance
  distance <- vegdist(data,method = 'bray')
  pcoa <- cmdscale(distance,k=2,eig=TRUE)
  plot_data <- data.frame({pcoa$points})[1:2]
  plot_data$group <- df_group[[group_col]]
  names(plot_data)[1:2] <- c('PCoA1','PCoA2')
  eig <- pcoa$eig
  
  # PERMANOVA
  perm <- adonis2(distance ~ df_group[[group_col]], permutations = 999)

  # R2,p
  r2 <- perm$R2[1]
  p_value <- perm$`Pr(>F)`[1]

  p_text <- paste("P =", format(p_value, digits = 3))
  r2_text <- paste("R² =", format(r2, digits = 3))
  stat_text <- paste(p_text, r2_text, sep = "\n")
  
  return(list(plot_data=plot_data,eig=eig,stat_text=stat_text))
}
plotdot <-function(plot_data,eig,stat_text){
  p <- ggplot(plot_data,aes(x=PCoA1,y=PCoA2,
                            #shape=group,
                            color=group)) +
    geom_point(size=2) +
    scale_color_manual(values = colors0) +
    stat_ellipse(type = "norm", level = 0.95, alpha = 0.9) +  
    theme_minimal() +
    labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digits=4),"%)",sep=""),
         y=paste("PCoA2(",format(100*eig[2]/sum(eig),digits=4),"%)",sep=""),
         color = 'Response',shape='Response') +
    annotate("text", 
             x = max(plot_data$PCoA1), 
             y = max(plot_data$PCoA2),
             label = stat_text,
             hjust = 0, 
             vjust = 1,
             size = 3,
             fontface = "bold") +
    theme(
      axis.line = element_line('black'),
      panel.grid = element_blank(),
      #plot.title = element_text(hjust = 0,vjust = 1.5,size = 8),
      text = element_text(family = 'Arial',face = 'bold'))
  return(p)
}
run_pcoa <- function(data,group_data,sample_col,group_col,colors0){
  result <- pcoa_solution(data,group_data,sample_col,group_col)
  p <- plotdot(result$plot_data,result$eig,result$stat_text)
  return(p)
}

colors0 <- c("#ff7a7f", "#46b2e5")
run_pcoa('S2/S2c_ASV.csv','S2/S2c_group.csv','Sample_ID','DMM',colors0)
ggsave('Fs2c.pdf',width = 6,height = 4)

########################## Fs1c ######################
library(ggplot2)
setwd('~/project/sun/new_fig/补充图/')
run_pcoa('S1/S1c_PCoA.csv','S1/S1c_group.csv','sample','class',colors0)
result <- pcoa_solution('S1/S1c_PCoA.csv','S1/S1c_group.csv','Sample_ID','Group',colors0)
colors0 <- c("#7b68ee", "#ffc025",
              "#efa7ca","#930054")
p <- plotdot(result$plot_data,result$eig,result$stat_text)
p
ggsave('Fs1c1.pdf',width = 5,height = 4)


######################### Fs1d ######################
rawdata <- read.csv('S1/S1d_PCoA.csv')
group_data <- read.csv('S1/S1d_group.csv')
run_pcoa('S1/S1d_PCoA.csv','S1/S1d_group.csv','Sample_ID','Group',colors0)
ggsave('Fs1d1.pdf',width = 5,height = 3)

######################## Fs5a ########################
setwd('~/project/sun/new_fig/补充图/')
rawdata <- read.csv('S5/S5a1_PCA_POS.csv')
group_data <- read.csv('S5/S5a_group.csv')
colors0 <- c("#ffc3c3","#a3a7de","#4dbbd5",'#eb8029')
run_pcoa('S5/S5a1_PCA_POS.csv','S5/S5a_group.csv','Sample_ID','Group',colors0)
ggsave('Fs5a1.pdf',width = 5,height = 3)

run_pcoa('S5/S5a2_PCA_NEG.csv','S5/S5a_group.csv','Sample_ID','Group',colors0)
ggsave('Fs5a2.pdf',width = 5,height = 3)

##################### Fs5b ########################
setwd('~/project/sun/new_fig/补充图/')
rawdata <- read.csv('S5/S5b_PCA.csv')
group_data <- read.csv('S5/S5b_group.csv')

run_pcoa('S5/S5b_PCA.csv','S5/S5b_group.csv','Sample_ID','Group',colors0)
ggsave('Fs5b.pdf',width = 5,height = 3)


###################### F3b ########################
setwd('~/project/sun/new_fig/')
colors0 <-  c("#7b68ee", "#ffc025",
              "#87cefa", "#efa7ca")
run_pcoa('Figure3/Figure3b.csv','Figure3/Figure3b_group.csv','Sample_ID','Stage',colors0)
ggsave('F3b.pdf',width = 5,height = 3)
