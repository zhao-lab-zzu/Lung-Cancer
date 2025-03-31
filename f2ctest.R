# 加载需要的包
library(pls)
library(caret)
library(ggplot2)

#' OSCoresPLS-PLSDA实现
#' @param X 特征矩阵(数值型矩阵)
#' @param y 类别标签(因子型向量)
#' @param ncomp PLS组分数量
#' @param north 正交信号修正迭代次数
#' @return 返回一个包含模型信息的列表
oscorespls_plsda <- function(X, y, ncomp = 2, north = 1) {
  # 检查输入
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.factor(y)) y <- as.factor(y)
  
  # 中心化数据
  # X_scaled <- scale(X, center = TRUE, scale = TRUE)
  X_scaled <- X
  y_dummy <- model.matrix(~y - 1)
  
  # 初始化
  X_corrected <- X_scaled
  
  # 正交信号修正
  for (i in 1:north) {
    # 计算PLS模型
    pls_orth <- plsr(y_dummy ~ X_corrected, ncomp = 1, method = "oscorespls")
    
    # 提取分数和载荷
    T <- scores(pls_orth)
    P <- loadings(pls_orth)
    
    # 计算正交分量
    Xorth <- T %*% t(P)
    
    # 从原始数据中移除正交分量
    X_corrected <- X_corrected - Xorth
  }
  
  # 训练最终的PLS-DA模型
  final_model <- plsr(y_dummy ~ X_corrected, ncomp = ncomp, method = "oscorespls")
  
  # 返回模型对象
  # 获取scores
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

#' 预测新样本
#' @param model oscorespls_plsda函数返回的模型对象
#' @param X_new 新的特征矩阵
#' @return 预测的类别标签
predict_oscorespls <- function(model, X_new) {
  # 检查输入
  if (!is.matrix(X_new)) X_new <- as.matrix(X_new)
  
  # 使用训练集的参数进行标准化
  X_scaled <- scale(X_new, 
                    center = model$center, 
                    scale = model$scale)
  
  # 预测
  pred <- predict(model$model, X_scaled, ncomp = model$model$ncomp)
  
  # 将预测结果转换为类别标签
  pred_class <- factor(model$classes[max.col(pred)],
                       levels = model$classes)
  
  return(pred_class)
}

#' 模型评估函数
#' @param true_labels 真实标签
#' @param pred_labels 预测标签
#' @return 评估指标列表
evaluate_model <- function(true_labels, pred_labels) {
  # 计算混淆矩阵
  conf_mat <- confusionMatrix(pred_labels, true_labels)
  
  # 返回评估结果
  return(list(
    accuracy = conf_mat$overall["Accuracy"],
    kappa = conf_mat$overall["Kappa"],
    conf_matrix = conf_mat$table
  ))
}

# 示例使用
demo_oscorespls_plsda <- function() {
  # 生成模拟数据
  set.seed(42)
  n_samples <- 100
  n_features <- 50
  
  # 生成特征矩阵
  X <- matrix(rnorm(n_samples * n_features), n_samples, n_features)
  
  # 生成响应变量（二分类）
  beta <- rnorm(n_features)
  y <- factor(ifelse(X %*% beta + rnorm(n_samples) * 0.1 > 0, "Class1", "Class2"))
  
  # 划分训练集和测试集
  train_idx <- sample(1:n_samples, 0.7 * n_samples)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # 训练模型
  model <- oscorespls_plsda(X_train, y_train, ncomp = 2, north = 1)
  
  # 预测
  pred_train <- predict_oscorespls(model, X_train)
  pred_test <- predict_oscorespls(model, X_test)
  
  # 评估模型
  train_eval <- evaluate_model(y_train, pred_train)
  test_eval <- evaluate_model(y_test, pred_test)
  
  # 打印结果
  cat("\n训练集评估结果:\n")
  print(train_eval)
  
  cat("\n测试集评估结果:\n")
  print(test_eval)
}

#' 绘制PLSDA得分散点图
#' @param model oscorespls_plsda函数返回的模型对象
#' @param y 类别标签
#' @param components 要绘制的组分，默认为c(1,2)
#' @return ggplot对象
plot_scores <- function(model, y, components = c(1,2)) {
  # 提取scores
  scores_mat <- model$scores[, components]
  
  # 创建数据框
  plot_data <- data.frame(
    Score1 = scores_mat[, 1],
    Score2 = scores_mat[, 2],
    Class = y
  )
  
  # 创建散点图
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

# 运行示例
demo_oscorespls_plsda()

# 创建示例数据的散点图
set.seed(42)
n_samples <- 100
n_features <- 50
X <- matrix(rnorm(n_samples * n_features), n_samples, n_features)
beta <- rnorm(n_features)
y <- factor(ifelse(X %*% beta + rnorm(n_samples) * 0.1 > 0, "Class1", "Class2"))

# 训练模型并绘图
model <- oscorespls_plsda(X, group$Group, ncomp = 2, north = 1)
scores_plot <- plot_scores(model, group$Group)
print(scores_plot)
