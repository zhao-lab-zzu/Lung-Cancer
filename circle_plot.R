# code for Circular barplot showing differentially abundant fungal genera across clinical stages
setwd('~/project/sun/new_fig/')
###################### Function to calculate fold change #####################
# Function to calculate logFC, handling zero values
calculate_logFC <- function(exp_values, ctrl_values, pseudocount = NULL) {
  # Automatically determine pseudocount value (if not provided)
  if(is.null(pseudocount)) {
    all_nonzero <- c(exp_values[exp_values > 0], ctrl_values[ctrl_values > 0])
    if(length(all_nonzero) > 0) {
      pseudocount <- min(all_nonzero) / 2  # Use half of the smallest non-zero value as pseudocount
    } else {
      pseudocount <- 0.01  # Default value if all are zero
    }
  }
  
  # Calculate mean values
  exp_mean <- mean(exp_values, na.rm = TRUE)
  ctrl_mean <- mean(ctrl_values, na.rm = TRUE)
  
  # Calculate logFC (log2 fold change)
  if(exp_mean == 0 && ctrl_mean == 0) {
    return(0)  # No change if both groups are zero
  } else if(ctrl_mean == 0) {
    return(log2((exp_mean + pseudocount) / pseudocount))
  } else if(exp_mean == 0) {
    return(log2(pseudocount / (ctrl_mean + pseudocount)))
  } else {
    return(log2(exp_mean / ctrl_mean))
  }
}

####################### Data processing #######################
abundance_data <- read.csv('relative_abundance_with_genus.csv',row.names = 1)
abundance_data <- abundance_data[!grepl("unclassified", rownames(abundance_data), ignore.case = TRUE), ]
abundance_data$Taxonomy <- rownames(abundance_data)
taxonomy_data <- abundance_data %>%
  select(Taxonomy) %>%
  separate(Taxonomy,into = c('domain','men','gang','mu','ke','shu'),sep = ';',extra = "drop") %>%
  select(shu,men)
taxonomy_data$species <- rownames(taxonomy_data)

sample_groups <- read.csv('group_testnew.csv',row.names = 1) %>%
  mutate(Stage = ifelse(Stage=="",S_Stage,Stage))
sample_groups <- sample_groups %>%
  mutate(Stage = ifelse(Stage=="","Control",Stage))
group_column <- "Stage"

abundance_data <- abundance_data %>%
  select(-Taxonomy)
## Convert to samples X species
transposed_data <- abundance_data %>%
  t() %>%
  as.data.frame()

###################### Differential significance analysis #######################
grouped_data <- transposed_data %>%
  cbind(group = sample_groups[['Group']])

results <- data.frame(species = colnames(transposed_data))
results[['p_value']] <- sapply(colnames(transposed_data), function(species) {
  # Get species abundance for the first group and the current group
  first_group_values <- grouped_data[grouped_data$group == 'NSCLC', species]
  current_group_values <- grouped_data[grouped_data$group == 'Control', species]
  
  
  # Perform Mann-Whitney test
  tryCatch({
    test_result <- wilcox.test(current_group_values, first_group_values)
    p_value <- test_result$p.value
  },error = function(e){
    p_value <<- 1
  })
  # Handle NaN p-values
  if(is.na(p_value) || is.nan(p_value)) {
    p_value <- 1  # Set NaN to 1, indicating no significant difference
  }
  # print(group_name)
  # print(p_value)
  # Calculate mean values to determine change direction
  first_mean <- mean(first_group_values)
  current_mean <- mean(current_group_values)
  return(p_value)
})

data_nsclc <- results %>%
  left_join(taxonomy_data,by = 'species') %>%
  filter(p_value < 0.5)
results <- data.frame(species = colnames(transposed_data))

results[['p_value']] <- sapply(colnames(transposed_data), function(species) {
  # Get species abundance for the first group and the current group
  first_group_values <- grouped_data[grouped_data$group == 'SCLC', species]
  current_group_values <- grouped_data[grouped_data$group == 'Control', species]
  
  
  # Perform Mann-Whitney test
  tryCatch({
    test_result <- wilcox.test(current_group_values, first_group_values)
    p_value <- test_result$p.value
  },error = function(e){
    p_value <<- 1
  })
  # Handle NaN p-values
  if(is.na(p_value) || is.nan(p_value)) {
    p_value <- 1  # Set NaN to 1, indicating no significant difference
  }
  # print(group_name)
  # print(p_value)
  # Calculate mean values to determine change direction
  first_mean <- mean(first_group_values)
  current_mean <- mean(current_group_values)
  return(p_value)
})

data_sclc <- results %>%
  left_join(taxonomy_data,by = 'species') %>%
  filter(p_value < 0.5)

data <- rbind(data_nsclc,data_sclc) %>%
  distinct(shu,.keep_all = TRUE)

filtered_data <- data %>%
  arrange(men) %>%
  filter(p_value<0.5)

################## Used for circular plot fold change analysis ###################
grouped_data <- transposed_data %>%
  cbind(group = sample_groups[[group_column]])



filted_trans_data <- transposed_data %>%
  select(filtered_data$species)



## logFC
new_results <- data.frame(species = colnames(filted_trans_data))
## p-value labels
new_results_pvalues <- data.frame(species = colnames(filted_trans_data))
for (group_name in unique(grouped_data$group)[-5]){
  new_results[[group_name]] <- sapply(colnames(filted_trans_data), function(species){
    first_group_values <- grouped_data[grouped_data$group == 'Control', species]
    current_group_values <- grouped_data[grouped_data$group == group_name, species]
    # Calculate logFC
    logFC <- calculate_logFC(current_group_values, first_group_values)
    # Wilcoxon test logic
    if(length(unique(first_group_values[!is.na(first_group_values)])) <= 1 || 
       length(unique(current_group_values[!is.na(current_group_values)])) <= 1 ||
       sum(!is.na(first_group_values)) < 2 ||
       sum(!is.na(current_group_values)) < 2) {
      return(0)
    }
    return(logFC)
  })
  
  
  new_results_pvalues[[group_name]] <- sapply(colnames(filted_trans_data), function(species){
    first_group_values <- grouped_data[grouped_data$group == 'Control', species]
    current_group_values <- grouped_data[grouped_data$group == group_name, species]
    
    # Check if we have enough data for statistical testing
    if(length(unique(first_group_values[!is.na(first_group_values)])) <= 1 || 
       length(unique(current_group_values[!is.na(current_group_values)])) <= 1 ||
       sum(!is.na(first_group_values)) < 2 ||
       sum(!is.na(current_group_values)) < 2) {
      return(NA)
    }
    
    # Perform Wilcoxon test
    test_result <- wilcox.test(current_group_values, first_group_values, exact = FALSE)
    p_value <- test_result$p.value
    
    # Add significance stars based on p-value
    if(is.na(p_value)) {
      return(NA)
    } else if(p_value < 0.001) {
      return("***")
    } else if(p_value < 0.01) {
      return("**")
    } else if(p_value < 0.05) {
      return("*")
    } else {
      return("")
    }
  })
}


plot_data <- filtered_data %>%
  left_join(new_results,by = 'species')
plot_data$shu <- gsub("g__","",plot_data$shu)
  
p_data <- filtered_data %>%
  left_join(new_results_pvalues,by='species')
p_data$shu <- gsub("g__","",p_data$shu)

####################### Plotting ######################
# Prepare circular heatmap data
# 1. Extract columns I-IV as matrix
mat <- as.matrix(plot_data[, c('Extensive stage','Limited stage',"IV", "III", "II", "I")])
rownames(mat) <- plot_data$shu

sig_mat <- as.matrix(p_data[,c('Extensive stage','Limited stage',"IV", "III", "II", "I")])
rownames(sig_mat)<- p_data$shu

sig_symbols <- matrix(
  sig_mat,
  nrow = nrow(sig_mat),
  ncol = ncol(sig_mat),
  dimnames = dimnames(sig_mat)
)
# 2. Extract men column data
men_data <- plot_data$men
names(men_data) <- plot_data$shu

# 3. Prepare species names
species_names <- plot_data$shu
names(species_names) <- plot_data$shu

# Set colors
# 1. Heatmap colors - Use blue-white-red scheme similar to the provided image
col_fun <- colorRamp2(c(min(mat), min(mat)/2, 0, max(mat)/2, max(mat)), c("#0364cb", "#02bbee", "#f0f0f8", "#e58c8e", "#e32924"))

unique_men <- unique(men_data)
men_palette <-  c("#7b68ee", "#ffc025","#3eab49",
                  "#87cefa", "#efa7ca","#a0725e")
men_colors <- setNames(men_palette,unique_men)
# Start creating circular heatmap
library(circlize)
circos.clear()

# Ensure sector names match completely
species <- plot_data$shu
gap_size <- c(rep(1, length(species)-1), 20)
names(gap_size) <- species

# Set parameters - place gap.after first
circos.par(gap.after = gap_size,
           #gap.degree = 1, 
           cell.padding = c(0, 0, 0, 0),
           track.margin = c(0.001, 0.001),
           start.degree = 90)
##################### save #############
{
pdf('F3-cicle_all.pdf',width = 6,height = 6)
# Initialize sectors - must be done after setting parameters
circos.initialize(factors = species, xlim = c(0, 1))


## Add labels
circos.track(ylim = c(0, 2), panel.fun = function(x, y) {
  sector.index = CELL_META$sector.index
  # Adjust angle to make text perpendicular to the circle
  theta = circlize(CELL_META$xcenter, 0.5, sector.index)[1, 1] %% 340
  text.facing = ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
  adj.x = ifelse(text.facing == 'clockwise',0,1)
  circos.text(CELL_META$xcenter, 2, sector.index, 
              facing = text.facing, 
              adj = c(adj.x, 0.5), 
              niceFacing = TRUE,
              cex = 0.4,
              family='Arial',
              font = 2)
}, track.height = 0.02, bg.border = NA)




# Add Stage track
track_labels <- c('Extensive stage','Limited stage',"IV", "III", "II", "I")
for(i in 1:ncol(mat)) {
  column_data <- mat[, i]
  circos.track(
    track.index = i + 1,  # Add 1 to account for men track
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector.index = CELL_META$sector.index
      value = column_data[sector.index]
      # Get significance symbol for this cell
      p_symbol = sig_symbols[sector.index, i]
      
      circos.rect(CELL_META$xlim[1], 0,
                  CELL_META$xlim[2], 1,
                  col = col_fun(value),
                  border = NA)
      
      # Add p-value symbol in the cell
      if(!is.na(p_symbol) && p_symbol != "") {
        # Determine text color - decide based on cell color intensity
        # Use white text for dark background; black text for light background
        bg_value <- value
        text_col <- ifelse(bg_value < 0, "white", "black")
        
        circos.text(
          mean(CELL_META$xlim),  # Cell center X
          0.5,                   # Cell center Y
          p_symbol,              # Significance symbol
          cex = 0.3,             # Text size
          col = 'black',        # Text color
          family = 'Arial',
          font = 2               # Bold
        )
      }
      
      # Add track label only at the first sector (starting position)
      if(sector.index == species[1]) {
        # Add label on the left side of the track
        circos.text(
          CELL_META$xlim[1] - 0.5,  # Slightly offset to the left
          0.5,                       # Vertically centered
          track_labels[i],           # Use corresponding column name as label
          facing = "bending.inside",      # Label orientation
          adj = c(1, 0.5),           # Alignment: right-aligned, vertically centered
          cex = 0.4,                 # Text size
          family = 'Arial',
          font = 2                   # Bold
        )
      }
    },
    track.height = 0.05,
    bg.border = NA
  )
  
}

# Add men track
circos.track(
  ylim = c(0, 2),
  panel.fun = function(x, y) {
    sector.index = CELL_META$sector.index
    men_value = men_data[sector.index]
    circos.rect(CELL_META$xlim[1], 0,
                CELL_META$xlim[2], 1,
                col = men_colors[men_value], 
                border = NA)
  }, 
  track.height = 0.05, 
  bg.border = NA
)
# Note: Here the order is reversed, the innermost ring is column IV, the outermost ring is column I
#circos.heatmap(mat, col = col_fun, track.height = 0.1)


## Use circlize package's built-in circular legend (gradient color ring)
breaks = seq(min(mat), max(mat), length.out = 5)
lgd = Legend(title = "logFC",
             col_fun = col_fun,
                          at = round(breaks, 2),
                          legend_width = unit(4, "cm"),
                          legend_height = unit(4, "cm"),
                          type = "colorbar",
                          direction = "vertical",
             # Font settings
             title_gp = gpar( fontface = "bold", fontfamily = "Arial"),  # Title font
             labels_gp = gpar(fontface = 'bold', fontfamily = "Arial"),  # Tick label font
             border = FALSE)
legend("bottomright", 
       inset = c(0.5, 0.4),  # Slightly offset upwards
       legend = names(men_colors), 
       fill = men_colors, 
       title = "Phylum", 
       text.font = 2,
       title.font = 2,
       bty = "n")


# Use grid package's viewport for positioning
library(grid)
# Add legend at the bottom of the plot
pushViewport(viewport(x = 0.5, y = 0.6, width = 0.8, height = 0.1))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
popViewport()
dev.off()
}


######################### Added bar plot ####################
## NSCLC, SCLC vs Control significant differences
nsclc <- data_nsclc %>%
  filter(p_value < 0.05)
sclc <- data_sclc %>%
  filter(p_value < 0.05)
## Union genus names
shus_all <- unique(nsclc$shu,sclc$shu)
## Intersection genus names
shus <- intersect(nsclc$shu,sclc$shu)

grouped_data <- transposed_data %>%
  cbind(group = sample_groups[['Stage']])

sample_group <- grouped_data %>%
  select(group)
all_data <- cbind(filted_trans_data,sample_group)

long_data <- all_data %>%
  pivot_longer(
    cols = -group,
    names_to = 'species',
    values_to = 'value'
  )


plot_data2 <- long_data %>%
  left_join(taxonomy_data %>% select(species,shu),by='species')
## Select significantly different genus
plot_data2 <- plot_data2 %>%
  filter(shu %in% shus_all)

groups <- unique(plot_data2$group)
comparisons <- combn(groups, 2, simplify = FALSE)

# Retain only combinations containing "Control"
control_comparisons <- Filter(function(pair) "Control" %in% pair, comparisons)

colors <- c("#a0725e", "#b9002a", "#70b5a3", "#006737", "#ee92a9", "#9dc4f0",
  "#ffa300", "#004e9c", "#634f93", "#e8c300", "#852d17")
plot_data2$group <- factor(plot_data2$group,levels = c('Control',"I", "II", "III", "IV",'Limited stage','Extensive stage'))
#plot_data2$group <- factor(plot_data2$group,levels = c('Control',"NSCLC", 'SCLC'))

plot_data3 <- plot_data2 %>%
  filter(shu %in% c('g__Asterotremella','g__Verticillium','g__Candida',
                       'g__Neosartorya','g__Thermoascus','g__Tilletiopsis')
         )

plot_data3 <- plot_data2 %>%
  filter(shu %in% c('g__Resupinatus','g__Schizophyllum')
  )

ggplot(plot_data3,aes(x=group,y=value,
                fill = group,
                color = group)
) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8),
             #size = 0.6) +  # Add points outside the box only
  geom_boxplot(width = 0.5,
               position = position_dodge(width = 0.8),
               outliers = FALSE,
               medcol = 'black',
  ) +
  
  stat_boxplot(geom = "errorbar",width = 0,show.legend = TRUE,
               
               position = position_dodge(width = 0.8),
  ) +
  stat_summary(fun = median, geom = "crossbar", 
               width = 0.5, color = "black",size=0.1) +
  ylim(0,0.01) +
  facet_wrap2(~ shu, 
              scales = "free",
              ncol = 3
  ) +
  
  theme_classic() +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  # scale_y_continuous(expand = c(0,0)) +
  #scale_color_manual(values = colors) +  # Manually specify colors
  labs(x = "", y = "Observed_ASVs",color = NULL,fill = NULL) +
  # stat_compare_means(
  #   label = 'p.format',  #"p.format", #"p.signif",
  #   method = "wilcox.test",
  #   comparisons = control_comparisons,
  #   bracket.size = 0.5,
  #   # label.y = max(data$value) + 0.2,
  #   show.legend = FALSE,
  #   hide.ns = TRUE,
  #   family = 'Arial',fontface = 'bold') +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = 'none',
    text = element_text(family = 'Arial', face = 'bold')
  )


ggsave('all_species5.pdf',width = 4,height = 2)
