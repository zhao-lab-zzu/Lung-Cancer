# code for NCM analysis
library(minpack.lm)
library(Hmisc)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Main function for NCM analysis
ncm_beta_analysis <- function(otu, Title) {
  # 1. Data preprocessing
  r_mean <- mean(apply(otu, 1, sum))  
  n_mean <- apply(otu, 2, mean)       
  n_mean <- n_mean[n_mean != 0]
  p <- n_mean / r_mean                
  
  # Calculate occurrence frequency
  spp.bi <- 1 * (otu > 0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  # Merge and filter data
  C <- merge(p, freq, by = 0)
  C <- C[order(C[, 2]), ]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]
  
  p <- C.0[, 2]
  freq <- C.0[, 3]
  names(p) <- C.0[, 1]
  names(freq) <- C.0[, 1]
  
  # 2. Model fitting
  d <- 1 / r_mean
  N <- ncol(otu)
  
  m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE),
                 start = list(m = 0.1))
  
  # Get confidence intervals
  m.ci <- confint(m.fit, 'm', level = 0.95)
  
  # Calculate predicted values
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), 
                     lower.tail = FALSE)
  
  # Predict confidence intervals
  pred.ci <- binconf(freq.pred * nrow(otu), nrow(otu), 
                     alpha = 0.05, method = "wilson", return.df = TRUE)
  
  # Calculate R-squared
  Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
  
  # 3. Integrate data
  results_df <- data.frame(
    p = p,
    freq = freq,
    freq.pred = freq.pred,
    Lower = pred.ci[, 2],
    Upper = pred.ci[, 3]
  )
  
  # Classify species
  results_df$group <- case_when(
    results_df$freq <= results_df$Lower ~ "Below",
    results_df$freq >= results_df$Upper ~ "Above",
    TRUE ~ "Neutral"
  )
  
  # 4. Create visualizations
  # Main plot
  main_plot <- ggplot(results_df, aes(x = log(p), y = freq)) +
    # Scatter plot
    geom_point(aes(color = group), size = 2) +
    # Prediction line
    geom_line(aes(y = freq.pred), color = "#4B4453") +
    # Confidence intervals
    geom_line(aes(y = Lower), color = "#4B4453", linetype = "dashed") +
    geom_line(aes(y = Upper), color = "#4B4453", linetype = "dashed") +
    scale_color_manual(values = colors) +
    ylim(0, 1) +
    ggtitle(Title) +
    # Labels
    labs(x = "log (Mean Relative Abundance)",
         y = "Occurrence Frequency",
         color = '') +
    annotate("text", x = min(log(results_df$p)), y = 0.9,
             label = sprintf("R² = %.3f\nNm = %.3f",
                             Rsqr,
                             N * coef(m.fit)
                             ),
             hjust = 0,
             family = 'Arial', fontface = 'bold'
             ) +
    theme_bw() +
    theme(
      legend.position = c(0.9, 0.15),
      legend.key = element_blank(),
      legend.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(size = 1.5),
      plot.title = element_text(hjust = 0.5),
      text = element_text(family = 'Arial', face = 'bold')
    )
  
  # Calculate group proportions
  group_counts <- table(results_df$group)
  group_props <- prop.table(group_counts) * 100
  
  # Create pie chart data
  pie_data <- data.frame(
    group = names(group_props),
    percentage = as.numeric(group_props)
  )
  library(ggrepel)
  # Pie chart
  pie_plot <- ggplot(pie_data, aes(x = "", y = percentage, fill = group)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_label_repel(aes(label = sprintf("%.1f%%", percentage),
                         y = cumsum(percentage) - percentage / 2),
                     size = 4,
                     family = 'Arial',
                     fontface = "bold",
                     show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    theme_void() +
    labs(title = '') +
    theme(legend.position = "none")
  
  # Combine plots
  combined_plot <- main_plot +
    annotation_custom(
      ggplotGrob(pie_plot),
      xmin = min(log(results_df$p)), xmax = min(log10(results_df$p)) + 0.5,
      ymin = 0.25, ymax = 1
    )
  
  # 5. Return results
  return(list(
    plot = combined_plot,
    results = results_df,
    m_value = coef(m.fit),
    m_ci = m.ci,
    Rsqr = Rsqr,
    Nm = N * coef(m.fit)
  ))
}


setwd('~/project/sun/new_fig/')
rawdata <- read.csv('抽平ASV_test.csv', row.names = 1)
data <- as.data.frame(t(rawdata))
data[] <- lapply(data, as.numeric)  

group <- read.csv('group_testnew.csv')
colors <- c("#FF6F91", "#C34A36", "#4B4453")
colors <- c("#2171B5", "#ffc025", "#87cefa")
colors <- c("#00375a", "#904e00", "#87cefa")

# Run function for analysis
run <- function(data, group_d, group_col, type_filter, Title) {
  # Check the format of type_filter to determine if it's a category list or range condition
  if (is.list(type_filter) && all(c("min", "max") %in% names(type_filter))) {
    # Range filtering
    samples <- group_d %>%
      filter(!!sym(group_col) >= type_filter$min & 
               !!sym(group_col) < type_filter$max) %>%
      pull(Sample_ID)
  } else {
    # Category filtering
    samples <- group_d %>%
      filter(!!sym(group_col) %in% type_filter) %>%
      pull(Sample_ID)
  }
  
  # Filter data and perform analysis
  otu <- data %>%
    filter(row.names(.) %in% samples)
  results <- ncm_beta_analysis(otu, Title)
  return(results$plot)
}

# Generate plots
plot_list <- list()
plot_list[[1]] <- run(data, group, 'Group', c('Control'), 'Control')
plot_list[[2]] <- run(data, group, 'Group', c('NSCLC'), 'NSCLC')
plot_list[[3]] <- run(data, group, 'Group', c('SCLC'), 'SCLC')
plot_list[[4]] <- run(data, group, 'Group', c('NSCLC', 'SCLC'), 'All_tumor')
plot_list[[5]] <- run(data, group, 'His', c('Adenocarcinoma'), 'Adenocarcinoma')
plot_list[[6]] <- run(data, group, 'His', c('Squamous carcinoma'), 'Squamous carcinoma')

plot_list[[7]] <- run(data, group, 'T', c('T1'), 'T1')
plot_list[[8]] <- run(data, group, 'T', c('T2'), 'T2')
plot_list[[9]] <- run(data, group, 'T', c('T3'), 'T3')
plot_list[[10]] <- run(data, group, 'T', c('T4'), 'T4')

plot_list[[11]] <- run(data, group, 'N', c('N0'), 'N0')
plot_list[[12]] <- run(data, group, 'N', c('N1'), 'N1')
plot_list[[13]] <- run(data, group, 'N', c('N2'), 'N2')
plot_list[[14]] <- run(data, group, 'N', c('N3'), 'N3')

plot_list[[15]] <- run(data, group, 'Stage', c('I'), 'Stage_I')
plot_list[[16]] <- run(data, group, 'Stage', c('II'), 'Stage_II')
plot_list[[17]] <- run(data, group, 'Stage', c('III'), 'Stage_III')
plot_list[[18]] <- run(data, group, 'Stage', c('IV'), 'Stage_IV')

plot_list[[19]] <- run(data, group, 'S_Stage', c('Extensive stage'), 'Extensive stage')
plot_list[[20]] <- run(data, group, 'S_Stage', c('Limited stage'), 'Limited stage')

plot_list[[21]] <- run(data, group, 'Gender', c('F'), 'Female')
plot_list[[22]] <- run(data, group, 'Gender', c('M'), 'Male')

plot_list[[23]] <- run(data, group, 'DMM', c(1), 'DMM1')
plot_list[[24]] <- run(data, group, 'DMM', c(2), 'DMM2')

library(patchwork)
p <- wrap_plots(plot_list, ncol = 4)
p
ggsave('NCM.pdf', width = 17, height = 23, plot = p)
ggsave('NCM1.pdf', width = 5, height = 4)

plot_list1 <- list()
plot_list1[[1]] <- run(data, group, 'BMI', list(min = 0, max = 18.5), 'Low BMI')
plot_list1[[2]] <- run(data, group, 'BMI', list(min = 18.5, max = 23), 'Medium BMI')
plot_list1[[3]] <- run(data, group, 'BMI', list(min = 23, max = 100), 'High BMI')

plot_list1[[4]] <- run(data, group, 'Age', list(min = 0, max = 60), 'Age young')
plot_list1[[5]] <- run(data, group, 'Age', list(min = 60, max = 1000), 'Age old')

plot_list1[[6]] <- run(data, group, 'M', c('M0'), 'M0')
plot_list1[[7]] <- run(data, group, 'M', c('M1'), 'M1')

plot_list1[[8]] <- run(data, group, 'Ab', c('yes'), 'Ab yes')
plot_list1[[9]] <- run(data, group, 'Ab', c('no'), 'Ab no')

plot_list1[[10]] <- run(data, group, 'Smoking', c('Currently smoking'), 'Currently smoking')
plot_list1[[11]] <- run(data, group, 'Smoking', c('Former smoker'), 'Former smoker')
plot_list1[[12]] <- run(data, group, 'Smoking', c('Never smoked'), 'Never smoked')

plot_list1[[13]] <- run(data, group, 'Smoking_history', c('yes'), 'Smoked')
plot_list1[[14]] <- run(data, group, 'Smoking_history', c('no'), 'No Smoked')

p <- wrap_plots(plot_list1, ncol = 4)

ggsave('NCM1.pdf', width = 17, height = 23 * 4 / 6, plot = p)