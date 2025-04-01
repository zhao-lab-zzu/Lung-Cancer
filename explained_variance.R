setwd('~/project/sun/new_fig/')
library(tidyverse)
library(vegan)
library(ggplot2)
library(patchwork)

# Load example serum metabolomics data
serum_data <- read.csv('ASV_test.csv', row.names = 1) %>%
  select(-Taxonomy)
serum_data <- as.data.frame(t(serum_data))
serum_data[] <- lapply(serum_data, as.numeric)

pheno_data <- read.csv('group_test.csv', row.names = 1)
# Calculate the Bray-Curtis distance matrix for serum metabolomics data
serum_dist <- vegdist(serum_data, method = "bray")

# Analyze each phenotypic variable
results <- data.frame(Variable = character(), R2 = numeric(), Correlation = numeric())

for (variable in colnames(pheno_data)) {
  complete_rows <- rownames(pheno_data)[!is.na(pheno_data[[variable]])]
  if(length(complete_rows) < 2) {
    # Skip the variable if there are too few samples after removing NAs
    next
  }
  # Create a subset containing only complete data
  complete_data <- pheno_data[complete_rows, , drop=FALSE]
  
  # Use rownames to select the corresponding subset of the distance matrix
  sub_serum_dist <- as.dist(as.matrix(serum_dist)[complete_rows, complete_rows])
  
  # Calculate R² for PERMANOVA
  formula <- as.formula(paste("sub_serum_dist ~", variable))
  adonis_result <- adonis2(formula, data = complete_data, permutations = 999)
  R2 <- adonis_result$R2[1]
  
  # Calculate correlation
  if (is.numeric(complete_data[[variable]])) {
    dist_vec <- as.vector(as.matrix(sub_serum_dist))
    var_dist <- as.vector(as.matrix(dist(complete_data[[variable]])))/max(dist(complete_data[[variable]]))
    correlation <- cor(dist_vec, var_dist, method = "spearman")
  } else {
    means <- tapply(as.vector(as.matrix(sub_serum_dist)), 
                    rep(complete_data[[variable]], each = nrow(complete_data)), 
                    mean)
    correlation <- diff(range(means))/diff(range(as.vector(as.matrix(sub_serum_dist))))
  }
  
  results <- rbind(results, 
                   data.frame(Variable = variable,
                              R2 = R2,
                              Correlation = correlation))
}

# Sort by R² value and add an order column
results <- results[order(-results$R2), ]
results$order <- 1:nrow(results)

# Create the first plot
p1 <- ggplot(results, aes(x = 1, y = reorder(Variable, -order))) +
  geom_tile(aes(fill = Correlation), width = 0.8, height = 0.8, color = "white") +  
  scale_fill_gradient(
    low = "#d02633",    
    high = "#4B9CD3",   
    limits = c(min(results$Correlation), max(results$Correlation)),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 15,      
      barheight = 0.8,    
      ticks.length = 0.2, 
      label.theme = element_text(size = 7), 
      title.theme = element_text(size = 7)   
    )
  ) +
  coord_equal() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    text = element_text(family = 'Arial', face = 'bold')
  )
p1

# Create the second plot
p2 <- ggplot(results, aes(y = reorder(Variable, -order))) +
  geom_bar(aes(x = R2), stat = "identity", fill = "black", width = 0.8) +
  scale_x_continuous(limits = c(0, 0.1), 
                     breaks = seq(0, 0.1, 0.02),
                     expand = c(0,0),
                     position = 'top') +
  labs(x = expression("Explained Variance in BC distance (R"^2*")")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(),
    panel.grid.major.y = element_blank(),
    # Move the x-axis to the top
    axis.x.position = "top",
    text = element_text(family = 'Arial', face = 'bold')
  )
p2

# Combine the plots using patchwork
combined_plot <- p1 + p2 + 
  plot_layout(widths = c(0.03, 0.3)) &
  theme(plot.margin = margin(0, 0, 0, ))

# Display the plot
print(combined_plot)
ggsave('varince_explained.pdf', width = 7, height = 7, plot = combined_plot,
       scale = 1,       
       device = cairo_pdf  
       )
