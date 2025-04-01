# Load required packages
library(pls)
library(caret)
library(ggplot2)

#' OSCoresPLS-PLSDA implementation
#' @param X Feature matrix (numeric matrix)
#' @param y Class labels (factor vector)
#' @param ncomp Number of PLS components
#' @param north Number of orthogonal signal correction iterations
#' @return Returns a list containing model information
oscorespls_plsda <- function(X, y, ncomp = 2, north = 1) {
  # Check input
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.factor(y)) y <- as.factor(y)
  
  # Center data
  # X_scaled <- scale(X, center = TRUE, scale = TRUE)
  X_scaled <- X
  y_dummy <- model.matrix(~y - 1)
  
  # Initialization
  X_corrected <- X_scaled
  
  for (i in 1:north) {
    pls_orth <- plsr(y_dummy ~ X_corrected, ncomp = 1, method = "oscorespls")
    T <- scores(pls_orth)
    P <- loadings(pls_orth)
    Xorth <- T %*% t(P)
    X_corrected <- X_corrected - Xorth
  }
  
  # PLS-DA
  final_model <- plsr(y_dummy ~ X_corrected, ncomp = ncomp, method = "oscorespls")
  
  scores_data <- scores(final_model)
  
  return(list(
    model = final_model,
    X_corrected = X_corrected,
    scores = scores_data,
    center = attr(X_scaled, "scaled:center"),
    scale = attr(X_scaled, "scaled:scale"),
    classes = levels(y)
  ))
}

#' Predict function
#' @param X_new New data matrix
#' @return Predicted class labels
predict_oscorespls <- function(model, X_new) {
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  # Scale
  X_scaled <- scale(X_new, 
                    center = model$center, 
                    scale = model$scale)
  pred <- predict(model$model, X_scaled, ncomp = model$model$ncomp)
  pred_class <- factor(model$classes[max.col(pred)],
                       levels = model$classes)
  
  return(pred_class)
}

#' Model evaluation
#' @param true_labels True class labels
#' @param pred_labels Predicted class labels
#' @return Evaluation results
evaluate_model <- function(true_labels, pred_labels) {
  # Calculate confusion matrix
  conf_mat <- confusionMatrix(pred_labels, true_labels)
  
  # Return evaluation results
  return(list(
    accuracy = conf_mat$overall["Accuracy"],
    kappa = conf_mat$overall["Kappa"],
    conf_matrix = conf_mat$table
  ))
}

demo_oscorespls_plsda <- function() {
  set.seed(42)
  n_samples <- 100
  n_features <- 50
  # Generate random feature matrix
  X <- matrix(rnorm(n_samples * n_features), n_samples, n_features)
  
  # Generate binary labels
  beta <- rnorm(n_features)
  y <- factor(ifelse(X %*% beta + rnorm(n_samples) * 0.1 > 0, "Class1", "Class2"))
  train_idx <- sample(1:n_samples, 0.7 * n_samples)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Train the model
  model <- oscorespls_plsda(X_train, y_train, ncomp = 2, north = 1)
  
  # Predict
  pred_train <- predict_oscorespls(model, X_train)
  pred_test <- predict_oscorespls(model, X_test)
  
  # Model evaluation
  train_eval <- evaluate_model(y_train, pred_train)
  test_eval <- evaluate_model(y_test, pred_test)
  
  cat("\nTraining set evaluation results:\n")
  print(train_eval)
  
  cat("\nTest set evaluation results:\n")
  print(test_eval)
}

#' Draw scores plot
#' @param model oscorespls_plsda model
#' @param y Label vector
#' @param components Indices of components to plot
#' @return ggplot object
plot_scores <- function(model, y, components = c(1,2)) {
  # Scores
  scores_mat <- model$scores[, components]
  
  # Create data frame for ggplot
  plot_data <- data.frame(
    Score1 = scores_mat[, 1],
    Score2 = scores_mat[, 2],
    Class = y
  )
  
  # Create ggplot
  p <- ggplot(plot_data, aes(x = Score1, y = Score2, color = Class)) +
    geom_point(size = 3, alpha = 0.6) +
    xlab(paste0("t[", components[1], "]")) +
    ylab(paste0("t[", components[2], "]")) +
    theme_bw() +
    ggtitle("PLSDA Scores Plot") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  return(p)
}

demo_oscorespls_plsda()

set.seed(42)
n_samples <- 100
n_features <- 50
X <- matrix(rnorm(n_samples * n_features), n_samples, n_features)
beta <- rnorm(n_features)
y <- factor(ifelse(X %*% beta + rnorm(n_samples) * 0.1 > 0, "Class1", "Class2"))

# Train-test split
model <- oscorespls_plsda(X, group$Group, ncomp = 2, north = 1)
scores_plot <- plot_scores(model, group$Group)
print(scores_plot)
