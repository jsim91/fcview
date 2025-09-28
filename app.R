library(xgboost) # xgboost regression
library(caret) # test/train splitting, mae, rmse
library(glmnet) # mlr regression
library(pls) # pc regression
library(ggplot2)
library(ggpubr)

load("D:/fluorescence/obare/uganda_bcell/review_outs_40000/ug_fcs_obj_label_small.RData")

data_meta <- obj_lab_app$metadata
data_freq <- obj_lab_app$leiden$abundance
colnames(data_freq) <- make.names(names = colnames(data_freq), unique = TRUE)
colnames(data_meta) <- make.names(names = colnames(data_meta), unique = TRUE)
data_freq <- as.data.frame(data_freq)
data_freq$patient_ID <- gsub(pattern = '_[0-9]+\\-[A-Za-z]+\\-[0-9]+\\.fcs$', replacement = '', x = row.names(data_freq))

data <- merge(x = data_meta, y = data_freq, by = 'patient_ID')
data <- data[!is.na(data$age),]
data <- data[,which(!colnames(data) %in% c('patient_ID','run_date','group','sex','hiv'))]

# data is similar to meta_sample


# general diagnostics helper; MLR and ENet assume homoskedasticity and normality of residuals
assemble_diagnostics_plot <- function(model_name,
                                      feat_imp_scaled,
                                      y_true = NULL,
                                      y_pred = NULL,
                                      model_obj = NULL,
                                      include_perm = FALSE,
                                      perm_rsq = NULL,
                                      obs_rsq = NULL,
                                      obs_rmse = NULL,
                                      obs_mae = NULL,
                                      p_val = NULL) {
  p_feat <- plot_feature_importance(feat_imp_scaled, model_label = model_name, order_label = "RawImportance")
  
  if (!is.null(model_obj) && inherits(model_obj, c("lm","glm"))) {
    p_diag <- plot_residual_diagnostics(model = model_obj, model_name = model_name)
  } else if (!is.null(y_true) && !is.null(y_pred)) {
    p_diag <- plot_residual_diagnostics(y_true = y_true, y_pred = y_pred, model_name = model_name)
  } else {
    p_diag <- NULL
  }
  
  p_perm <- if (include_perm && !is.null(perm_rsq)) {
    plot_permutation_histogram(perm_rsq, obs_rsq, obs_rmse, obs_mae, p_val, model_name)
  } else NULL
  
  plots <- list(p_perm, p_feat, p_diag)
  plots <- plots[!sapply(plots, is.null)]
  ggpubr::ggarrange(plotlist = plots, ncol = 1, nrow = length(plots))
}

# General-purpose splitter
make_train_test_split <- function(data, outcome_col, p = 0.7, seed = 123) {
  # Partition indices
  set.seed(seed); parts <- createDataPartition(data[[outcome_col]], p = p, list = FALSE)
  
  train <- data[parts[,1], ]
  test  <- data[-parts[,1], ]
  
  # Predictor/response matrices
  outcome_colnum <- which(colnames(data) == outcome_col)
  
  train_x <- data.matrix(train[, -outcome_colnum, drop = FALSE])
  train_y <- train[[outcome_col]]
  
  test_x  <- data.matrix(test[, -outcome_colnum, drop = FALSE])
  test_y  <- test[[outcome_col]]
  
  list(
    train = train,
    test = test,
    train_x = train_x,
    train_y = train_y,
    test_x = test_x,
    test_y = test_y
  )
}

# permutation hist helper
plot_permutation_histogram <- function(perm_rsq, obs_rsq, obs_rmse, obs_mae, p_val, model_name) {
  metrics_text <- paste0(
    "Observed R² = ", round(obs_rsq, 3), "\n",
    "RMSE = ", round(obs_rmse, 3), "\n",
    "MAE = ", round(obs_mae, 3), "\n",
    "p-value = ", signif(p_val, 3)
  )
  
  ggplot(data.frame(perm_rsq), aes(x = perm_rsq)) +
    geom_histogram(fill = "lightblue", color = "black", bins = 30, linewidth = 0.3) +
    geom_vline(xintercept = obs_rsq, color = "red", size = 1.2) +
    annotate("text", x = min(perm_rsq), 
             y = max(hist(perm_rsq, plot = FALSE)$counts),
             label = metrics_text, hjust = 0, vjust = 1, size = 3.5) +
    labs(title = paste0("Permutation Test: R² Distribution (", model_name, ")"),
         x = "R² under null (permuted labels)", y = "Count") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
}

# feature importance helper
# Replace your plotting helper with this programmatic-safe version
plot_feature_importance <- function(feat_imp_scaled,
                                    model_label = "Model",
                                    order_label = "RawImportance") {
  df <- as.data.frame(feat_imp_scaled)
  required <- c("Feature", "RawImportance", "ScaledImportance")
  if (!all(required %in% names(df)) || nrow(df) == 0) {
    return(ggplot() +
             annotate("text", x = 0, y = 0, label = "No feature importance available") +
             theme_void())
  }
  ord <- order(-df[["RawImportance"]])
  df <- df[ord, , drop = FALSE]
  df$Feature <- factor(df$Feature, levels = rev(df$Feature))
  
  long <- tidyr::pivot_longer(
    df,
    cols = c("RawImportance", "ScaledImportance"),
    names_to = "Type",
    values_to = "Importance"
  )
  
  ggplot(long, aes(x = Feature, y = Importance, fill = Type)) +
    geom_col(position = position_dodge(width = 0.7)) +
    coord_flip() +
    scale_fill_manual(values = c("RawImportance" = "steelblue",
                                 "ScaledImportance" = "red3"),
                      name = "Importance type") +
    labs(title = paste0(model_label, " feature importance: raw vs scaled"),
         subtitle = paste("Ordered by:", order_label),
         x = "Feature", y = "Importance") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
}

## general purpose plot test predictions
plot_model_performance <- function(train_y, train_pred, test_y, test_pred,
                                   model_name = "Model") {
  # ---- Metrics ----
  # Train R²
  sst_train <- sum((train_y - mean(train_y))^2)
  sse_train <- sum((train_y - train_pred)^2)
  train_r2 <- 1 - sse_train / sst_train
  
  # Test R²
  sst_test <- sum((test_y - mean(test_y))^2)
  sse_test <- sum((test_y - test_pred)^2)
  test_r2 <- 1 - sse_test / sst_test
  
  rmse_test <- sqrt(mean((test_y - test_pred)^2))
  mae_test  <- mean(abs(test_y - test_pred))
  
  # ---- Scatter plot ----
  df <- data.frame(
    Actual = test_y,
    Predicted = test_pred
  )
  
  metrics_text <- paste0(
    "Train R² = ", round(train_r2, 3), "\n",
    "Test R² = ", round(test_r2, 3), "\n",
    "RMSE = ", round(rmse_test, 3), "\n",
    "MAE = ", round(mae_test, 3)
  )
  
  p <- ggplot(df, aes(x = Actual, y = Predicted)) +
    geom_point(shape = 21, fill = "grey40", color = "black", stroke = 0.2, size = 2, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
    geom_hline(yintercept = mean(test_y), color = "blue", linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf, label = metrics_text,
             hjust = -0.05, vjust = 1.05, size = 3.5) +
    labs(title = paste0(model_name, ": Actual vs Predicted (Test Set)"),
         x = "Actual values", y = "Predicted values") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  list(
    train_r2 = train_r2,
    test_r2 = test_r2,
    rmse_test = rmse_test,
    mae_test = mae_test,
    plot = p
  )
}

# residual diagnostic helper
plot_residual_diagnostics <- function(y_true = NULL, y_pred = NULL, 
                                      model = NULL, model_name = "Model") {
  library(ggplot2)
  library(ggpubr)
  library(car)
  
  # If model is provided and supports lm methods
  if (!is.null(model) && inherits(model, c("lm", "glm"))) {
    resid_vec   <- resid(model)
    fitted_vec  <- fitted(model)
    std_resid   <- rstandard(model)
    leverage    <- hatvalues(model)
    cooksd      <- cooks.distance(model)
    
    # Tests
    ks_res  <- ks.test(resid_vec, "pnorm", mean = mean(resid_vec), sd = sd(resid_vec))
    ncv_res <- car::ncvTest(model)
    
    # Threshold for Cook’s distance
    n <- length(resid_vec)
    cook_threshold <- 4 / n
    
    # Residuals vs Fitted
    p_resid_fit <- ggplot(data.frame(fitted = fitted_vec, resid = resid_vec),
                          aes(x = fitted, y = resid)) +
      geom_point(alpha = 0.7, fill = 'grey40', pch = 21, color = 'black', stroke = 0.05, size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      annotate("text", x = -Inf, y = Inf,
               label = paste0("ncvTest p = ", signif(ncv_res$p, 3)),
               hjust = -0.05, vjust = 1.05, size = 3.5) +
      labs(title = "Residuals vs Fitted – Skedasticity",
           x = "Fitted values", y = "Residuals") +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    # Normal Q-Q
    p_qq <- ggplot(data.frame(resid = resid_vec), aes(sample = resid)) +
      stat_qq(fill = "steelblue", pch = 21, color = 'black', stroke = 0.05, size = 2) +
      stat_qq_line(color = "red", linetype = "dashed") +
      annotate("text", x = -Inf, y = Inf,
               label = paste0("KS test p = ", signif(ks_res$p.value, 3)),
               hjust = -0.05, vjust = 1.05, size = 3.5) +
      labs(title = "Normal Q–Q plot – Normality",
           x = "Theoretical quantiles", y = "Sample quantiles") +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    # Scale-Location
    p_scale <- ggplot(data.frame(fitted = fitted_vec,
                                 sqrt_std_resid = sqrt(abs(std_resid))),
                      aes(x = fitted, y = sqrt_std_resid)) +
      geom_point(alpha = 0.7, pch = 21, fill = 'grey40', color = 'black', stroke = 0.05, size = 2) +
      geom_smooth(se = FALSE, color = "red") +
      labs(title = "Scale–Location",
           x = "Fitted values", y = "√|Standardized residuals|") +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    # Residuals vs Leverage
    p_leverage <- ggplot(data.frame(leverage = leverage,
                                    resid = resid_vec,
                                    cook = cooksd),
                         aes(x = leverage, y = resid)) +
      geom_point(aes(size = cook), alpha = 0.7, pch = 21, fill = 'grey40', color = 'black', stroke = 0.1) +
      geom_smooth(se = FALSE, color = "red") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey50") +
      geom_hline(yintercept = cook_threshold, linetype = "dashed", color = "blue") +
      scale_size_continuous(name = "Cook's D") +
      labs(title = "Residuals vs Leverage",
           x = "Leverage", y = "Residuals") +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    return(ggarrange(p_resid_fit, p_qq, p_scale, p_leverage, ncol = 2, nrow = 2))
  }
  
  # ---- Generic fallback for models without lm diagnostics ----
  if (is.null(y_true) || is.null(y_pred)) {
    stop("For non-lm models, please provide y_true and y_pred.")
  }
  
  resid_vec  <- y_true - y_pred
  fitted_vec <- y_pred
  
  # Residuals vs Fitted
  p_resid_fit <- ggplot(data.frame(fitted = fitted_vec, resid = resid_vec),
                        aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.7, fill = 'grey40', pch = 21, color = 'black', stroke = 0.05, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste("Residuals vs Fitted –", model_name),
         x = "Fitted values", y = "Residuals") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  # Normal Q-Q
  p_qq <- ggplot(data.frame(resid = resid_vec), aes(sample = resid)) +
    stat_qq(fill = "steelblue", pch = 21, color = 'black', stroke = 0.05, size = 2) +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(title = paste("Normal Q–Q plot –", model_name),
         x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  # Histogram of residuals
  p_hist <- ggplot(data.frame(resid = resid_vec), aes(x = resid)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    labs(title = paste("Residual Distribution –", model_name),
         x = "Residuals", y = "Count") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  
  ggarrange(p_resid_fit, p_qq, p_hist, ncol = 2, nrow = 2)
}

# importance scale helper
scale_importance <- function(feat_imp, data, outcome_col, order_by = NULL) {
  data <- as.data.frame(data)
  predictors <- data[ , setdiff(names(data), outcome_col), drop = FALSE]
  cor_mat <- cor(predictors, use = "pairwise.complete.obs")
  
  # Identify the raw importance column
  if (!is.null(order_by) && order_by %in% names(feat_imp)) {
    raw_col <- order_by
  } else if ("%IncMSE" %in% names(feat_imp)) {
    raw_col <- "%IncMSE"
  } else if ("Coefficient" %in% names(feat_imp)) {
    raw_col <- "Coefficient"
  } else if ("Gain" %in% names(feat_imp)) {
    raw_col <- "Gain"
  } else {
    stop("No suitable importance column found in feat_imp")
  }
  
  # Extract raw values
  raw_vals <- feat_imp[[raw_col]]
  
  # Ensure Feature column exists
  if (!"Feature" %in% names(feat_imp)) {
    stop("feat_imp must contain a 'Feature' column")
  }
  
  # Order by descending raw importance
  feat_imp <- feat_imp[order(-raw_vals), ]
  raw_vals <- feat_imp[[raw_col]]
  
  # Apply correlation-adjusted scaling
  scaled <- numeric(length(raw_vals))
  for (i in seq_along(raw_vals)) {
    f <- feat_imp$Feature[i]
    if (i == 1) {
      scaled[i] <- raw_vals[i]
    } else {
      higher_feats <- feat_imp$Feature[1:(i-1)]
      if (f %in% colnames(cor_mat)) {
        max_corr <- max(abs(cor_mat[f, higher_feats]), na.rm = TRUE)
      } else {
        max_corr <- 0
      }
      scaled[i] <- raw_vals[i] * (1 - max_corr)
    }
  }
  
  # Return unified table
  out <- data.frame(
    Feature = feat_imp$Feature,
    RawImportance = raw_vals,
    ScaledImportance = scaled,
    stringsAsFactors = FALSE
  )
  out
}

# make.names helper
sanitize_predictors <- function(x, outcome_name = NULL, return_map = FALSE) {
  df <- as.data.frame(x)
  
  # Ensure column names exist and are unique
  if (is.null(colnames(df))) {
    orig_names <- paste0("V", seq_len(ncol(df)))
  } else {
    orig_names <- colnames(df)
  }
  new_names <- make.names(orig_names, unique = TRUE)
  colnames(df) <- new_names
  
  # Drop outcome if requested
  if (!is.null(outcome_name) && outcome_name %in% colnames(df)) {
    df <- df[, setdiff(colnames(df), outcome_name), drop = FALSE]
  }
  
  if (return_map) {
    name_map <- setNames(orig_names, new_names)  # sanitized → original
    return(list(data = df, name_map = name_map))
  } else {
    return(df)
  }
}
make_mm <- function(train_df, test_df, outcome_col) {
  train_s <- sanitize_predictors(train_df, outcome_name = outcome_col)
  test_s  <- sanitize_predictors(test_df, outcome_name = outcome_col)
  
  all_df <- rbind(train_s, test_s)
  mm_all <- model.matrix(~ . - 1, data = all_df)
  n_tr <- nrow(train_s)
  x_tr <- mm_all[seq_len(n_tr), , drop = FALSE]
  x_te <- mm_all[(n_tr + 1):nrow(mm_all), , drop = FALSE]
  list(train_x = x_tr, test_x = x_te)
}






## xgboost modeling
xgb_permutation_eval <- function(data, outcome_col = "age",
                                 validation = c("train_test", "cv"),
                                 p = 0.7, k = 5,
                                 n_perm = 500, seed = 123,
                                 nrounds = 500,
                                 max_depth = 6,
                                 early_stopping_rounds = 10,
                                 order_by = "Gain",
                                 tolerance = 1e-4) {
  validation <- match.arg(validation)
  set.seed(seed)
  
  if (validation == "train_test") {
    split_obj <- make_train_test_split(data, outcome_col, p, seed)
    train_x <- sanitize_predictors(split_obj$train_x)
    test_x  <- sanitize_predictors(split_obj$test_x)
    colnames(test_x) <- colnames(train_x)
    train_y <- split_obj$train_y
    test_y  <- split_obj$test_y
    
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_x), label = train_y)
    dtest  <- xgboost::xgb.DMatrix(data = as.matrix(test_x),  label = test_y)
    
    watchlist <- list(train = dtrain, eval = dtest)
    final <- xgboost::xgb.train(
      data = dtrain,
      objective = "reg:squarederror",
      nrounds = nrounds,
      max_depth = max_depth,
      watchlist = watchlist,
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    train_pred <- predict(final, dtrain)
    test_pred  <- predict(final, dtest)
    
    sst <- sum((test_y - mean(test_y))^2)
    sse <- sum((test_y - test_pred)^2)
    obs_rsq <- 1 - sse/sst
    obs_rmse <- sqrt(mean((test_y - test_pred)^2))
    obs_mae  <- mean(abs(test_y - test_pred))
    
    imp <- xgboost::xgb.importance(model = final)
    feat_imp <- data.frame(Feature = imp$Feature)
    feat_imp[[order_by]] <- imp[[order_by]]
    
    train_df_for_scale <- data.frame(train_x, outcome = train_y)
    names(train_df_for_scale)[ncol(train_df_for_scale)] <- outcome_col
    
    feat_imp_scaled <- scale_importance(
      feat_imp,
      data = train_df_for_scale,
      outcome_col = outcome_col
    )
    
    perm_rsq <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      y_perm <- sample(train_y)
      dtrain_perm <- xgboost::xgb.DMatrix(data = as.matrix(train_x), label = y_perm)
      perm_model <- xgboost::xgb.train(
        data = dtrain_perm,
        objective = "reg:squarederror",
        nrounds = final$best_iteration,
        max_depth = max_depth,
        verbose = 0
      )
      y_pred_perm <- predict(perm_model, dtest)
      sst_perm <- sum((test_y - mean(test_y))^2)
      sse_perm <- sum((test_y - y_pred_perm)^2)
      perm_rsq[i] <- 1 - sse_perm/sst_perm
    }
    p_val <- mean(perm_rsq >= obs_rsq)
    
    p1 <- plot_permutation_histogram(perm_rsq, obs_rsq, obs_rmse, obs_mae, p_val, "XGBoost")
    p2 <- plot_feature_importance(feat_imp_scaled, model_label = "XGBoost", order_label = order_by)
    diagnostics_plot <- plot_residual_diagnostics(y_true = test_y, y_pred = test_pred, model_name = "XGBoost")
    combined_plot <- ggpubr::ggarrange(p1, p2, diagnostics_plot, ncol = 1, nrow = 3)
    
    perf <- plot_model_performance(train_y, train_pred, test_y, test_pred, "XGBoost")
    
    return(list(
      observed_metrics = list(R2 = obs_rsq, RMSE = obs_rmse, MAE = obs_mae),
      feature_importance = feat_imp_scaled,
      permuted_rsq = perm_rsq,
      p_value = p_val,
      model = final,
      permutation_importance_diagnostics_plot = combined_plot,
      performance_plot = perf$plot,
      train_r2 = perf$train_r2,
      test_r2 = perf$test_r2,
      rmse_test = perf$rmse_test,
      mae_test = perf$mae_test
    ))
    
  } else if (validation == "cv") {
    folds <- createFolds(data[[outcome_col]], k = k, list = TRUE)
    all_preds <- numeric(nrow(data))
    fold_metrics <- lapply(seq_along(folds), function(i) {
      test_idx <- folds[[i]]
      train_idx <- setdiff(seq_len(nrow(data)), test_idx)
      train_fold <- data[train_idx, ]
      test_fold  <- data[test_idx, ]
      
      train_x <- sanitize_predictors(train_fold[, setdiff(names(train_fold), outcome_col), drop = FALSE])
      test_x  <- sanitize_predictors(test_fold[, setdiff(names(test_fold), outcome_col), drop = FALSE])
      colnames(test_x) <- colnames(train_x)
      train_y <- train_fold[[outcome_col]]
      test_y  <- test_fold[[outcome_col]]
      
      dtrain <- xgboost::xgb.DMatrix(data = as.matrix(train_x), label = train_y)
      dtest  <- xgboost::xgb.DMatrix(data = as.matrix(test_x),  label = test_y)
      
      watchlist <- list(train = dtrain, eval = dtest)
      model <- xgboost::xgb.train(
        data = dtrain,
        objective = "reg:squarederror",
        nrounds = nrounds,
        max_depth = max_depth,
        watchlist = watchlist,
        early_stopping_rounds = early_stopping_rounds,
        verbose = 0
      )
      y_pred <- predict(model, dtest)
      all_preds[test_idx] <<- y_pred
      
      sst <- sum((test_y - mean(test_y))^2)
      sse <- sum((test_y - y_pred)^2)
      r2  <- 1 - sse/sst
      rmse <- sqrt(mean((test_y - y_pred)^2))
      mae  <- mean(abs(test_y - y_pred))
      
      list(R2 = r2, RMSE = rmse, MAE = mae)
    })
    
    mean_r2   <- mean(sapply(fold_metrics, `[[`, "R2"))
    mean_rmse <- mean(sapply(fold_metrics, `[[`, "RMSE"))
    mean_mae  <- mean(sapply(fold_metrics, `[[`, "MAE"))
    
    perf <- plot_model_performance(
      train_y = data[[outcome_col]],
      train_pred = rep(mean(data[[outcome_col]]), nrow(data)),
      test_y = data[[outcome_col]],
      test_pred = all_preds,
      model_name = paste0("XGBoost (", k, "-fold CV)")
    )
    
    train_x <- sanitize_predictors(data[, setdiff(names(data), outcome_col)])
    train_y <- data[[outcome_col]]
    dtrain_full <- xgboost::xgb.DMatrix(data = as.matrix(train_x), label = train_y)
    final_model <- xgboost::xgb.train(
      data = dtrain_full,
      objective = "reg:squarederror",
      nrounds = nrounds,
      max_depth = max_depth,
      verbose = 0
    )
    imp <- xgboost::xgb.importance(model = final_model)
    feat_imp <- data.frame(Feature = imp$Feature)
    feat_imp[[order_by]] <- imp[[order_by]]
    
    final_df_for_scale <- data.frame(train_x, outcome = train_y)
    names(final_df_for_scale)[ncol(final_df_for_scale)] <- outcome_col
    
    feat_imp_scaled <- scale_importance(
      feat_imp,
      data = final_df_for_scale,
      outcome_col = outcome_col
    )
    
    p2 <- plot_feature_importance(feat_imp_scaled, model_label = "XGBoost", order_label = order_by)
    diagnostics_plot <- plot_residual_diagnostics(y_true = data[[outcome_col]], y_pred = all_preds, model_name = paste0("XGBoost (", k, "-fold CV)"))
    combined_plot <- ggpubr::ggarrange(p2, diagnostics_plot, ncol = 1, nrow = 2)
    
    return(list(
      cv_metrics = list(mean_R2 = mean_r2, mean_RMSE = mean_rmse, mean_MAE = mean_mae),
      fold_metrics = fold_metrics,
      feature_importance = feat_imp_scaled,
      permutation_importance_diagnostics_plot = combined_plot,
      performance_plot = perf$plot,
      model = final_model
    ))
  }
}

# ---- Train/test mode ----
xgb_results_tt <- xgb_permutation_eval(
  data = data,
  outcome_col = "age",
  validation = "train_test",
  p = 0.7,            # proportion of data for training
  max_depth = 3,
  nrounds = 200,
  n_perm = 200,
  order_by = "Gain"   # can also use "Cover" or "Frequency"
)

# Inspect metrics
xgb_results_tt$observed_metrics
xgb_results_tt$train_r2
xgb_results_tt$test_r2
xgb_results_tt$p_value
xgb_results_tt$best_iter

# Feature importance
head(xgb_results_tt$raw_feature_importance, 5)
head(xgb_results_tt$scaled_feature_importance, 5)

# Plots
xgb_results_tt$performance_plot
xgb_results_tt$permutation_importance_plot


# ---- Cross-validation mode ----
xgb_results_cv <- xgb_permutation_eval(
  data = data,
  outcome_col = "age",
  validation = "cv",
  k = 5,              # number of folds
  max_depth = 3,
  nrounds = 200,
  order_by = "Gain"
)

# Inspect CV metrics
xgb_results_cv$cv_metrics        # mean R², RMSE, MAE across folds
xgb_results_cv$fold_metrics      # per-fold metrics

# Feature importance from final model
xgb_results_cv$raw_feature_importance
xgb_results_cv$scaled_feature_importance

# Plots
xgb_results_cv$performance_plot
xgb_results_cv$permutation_importance_plot


















## multiple linear regression
mlr_permutation_eval <- function(data, outcome_col = "age",
                                 validation = c("train_test", "cv"),
                                 p = 0.7, k = 5,
                                 n_perm = 500, seed = 123,
                                 tolerance = 1e-4) {
  validation <- match.arg(validation)
  set.seed(seed)
  
  if (validation == "train_test") {
    split_obj <- make_train_test_split(data, outcome_col, p, seed)
    train <- split_obj$train
    test  <- split_obj$test
    
    train_x <- sanitize_predictors(split_obj$train_x)
    test_x  <- sanitize_predictors(split_obj$test_x)
    colnames(test_x) <- colnames(train_x)
    
    train_y <- split_obj$train_y
    test_y  <- split_obj$test_y
    
    formula <- as.formula(paste(outcome_col, "~ ."))
    mlr_model <- lm(formula, data = train)
    
    train_pred <- predict(mlr_model, train)
    test_pred  <- predict(mlr_model, as.data.frame(test_x))
    
    sst <- sum((test_y - mean(test_y))^2)
    sse <- sum((test_y - test_pred)^2)
    obs_rsq <- 1 - sse/sst
    obs_rmse <- sqrt(mean((test_y - test_pred)^2))
    obs_mae  <- mean(abs(test_y - test_pred))
    
    coefs <- coef(mlr_model)[-1]
    feat_imp <- data.frame(
      Feature = names(coefs),
      Coefficient = as.numeric(coefs),
      stringsAsFactors = FALSE
    )
    feat_imp_scaled <- scale_importance(
      feat_imp,
      data = train,
      outcome_col = outcome_col
    )
    
    perm_rsq <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      y_perm <- sample(train[[outcome_col]])
      temp <- train
      temp[[outcome_col]] <- y_perm
      perm_model <- lm(formula, data = temp)
      y_pred_perm <- predict(perm_model, as.data.frame(test_x))
      sst_perm <- sum((test_y - mean(test_y))^2)
      sse_perm <- sum((test_y - y_pred_perm)^2)
      perm_rsq[i] <- 1 - sse_perm/sst_perm
    }
    p_val <- mean(perm_rsq >= obs_rsq)
    
    p1 <- plot_permutation_histogram(perm_rsq, obs_rsq, obs_rmse, obs_mae, p_val, "MLR")
    p2 <- plot_feature_importance(feat_imp_scaled, model_label = "MLR", order_label = "RawImportance")
    diagnostics_plot <- plot_residual_diagnostics(model = mlr_model, model_name = "MLR")
    combined_plot <- ggarrange(p1, p2, diagnostics_plot, ncol = 1, nrow = 3)
    
    perf <- plot_model_performance(train_y, train_pred, test_y, test_pred, "MLR")
    
    return(list(
      observed_metrics = list(R2 = obs_rsq, RMSE = obs_rmse, MAE = obs_mae),
      feature_importance = feat_imp_scaled,
      permuted_rsq = perm_rsq,
      p_value = p_val,
      model = mlr_model,
      permutation_importance_diagnostics_plot = combined_plot,
      performance_plot = perf$plot,
      train_r2 = perf$train_r2,
      test_r2 = perf$test_r2,
      rmse_test = perf$rmse_test,
      mae_test = perf$mae_test
    ))
    
  } else if (validation == "cv") {
    folds <- createFolds(data[[outcome_col]], k = k, list = TRUE)
    
    all_preds <- numeric(nrow(data))
    fold_metrics <- lapply(seq_along(folds), function(i) {
      test_idx <- folds[[i]]
      train_idx <- setdiff(seq_len(nrow(data)), test_idx)
      
      train_fold <- data[train_idx, ]
      test_fold  <- data[test_idx, ]
      
      formula <- as.formula(paste(outcome_col, "~ ."))
      mlr_model <- lm(formula, data = train_fold)
      y_pred <- predict(mlr_model, test_fold)
      all_preds[test_idx] <<- y_pred
      
      y_test <- test_fold[[outcome_col]]
      sst <- sum((y_test - mean(y_test))^2)
      sse <- sum((y_test - y_pred)^2)
      r2  <- 1 - sse/sst
      rmse <- sqrt(mean((y_test - y_pred)^2))
      mae  <- mean(abs(y_test - y_pred))
      
      list(R2 = r2, RMSE = rmse, MAE = mae)
    })
    
    mean_r2   <- mean(sapply(fold_metrics, `[[`, "R2"))
    mean_rmse <- mean(sapply(fold_metrics, `[[`, "RMSE"))
    mean_mae  <- mean(sapply(fold_metrics, `[[`, "MAE"))
    
    perf <- plot_model_performance(
      train_y = data[[outcome_col]],
      train_pred = rep(mean(data[[outcome_col]]), nrow(data)),
      test_y = data[[outcome_col]],
      test_pred = all_preds,
      model_name = paste0("MLR (", k, "-fold CV)")
    )
    
    final_model <- lm(as.formula(paste(outcome_col, "~ .")), data = data)
    coefs <- coef(final_model)[-1]
    feat_imp <- data.frame(
      Feature = names(coefs),
      Coefficient = as.numeric(coefs),
      stringsAsFactors = FALSE
    )
    feat_imp_scaled <- scale_importance(
      feat_imp,
      data = data,
      outcome_col = outcome_col
    )
    
    p2 <- plot_feature_importance(feat_imp_scaled, model_label = "MLR", order_label = "RawImportance")
    diagnostics_plot <- plot_residual_diagnostics(model = final_model, model_name = "MLR")
    combined_plot <- ggarrange(p2, diagnostics_plot, ncol = 1, nrow = 2)
    
    return(list(
      cv_metrics = list(mean_R2 = mean_r2, mean_RMSE = mean_rmse, mean_MAE = mean_mae),
      fold_metrics = fold_metrics,
      feature_importance = feat_imp_scaled,
      permutation_importance_diagnostics_plot = combined_plot,
      performance_plot = perf$plot,
      model = final_model
    ))
  }
}

# ---- Train/test mode ----
mlr_results_tt <- mlr_permutation_eval(
  data = data,
  outcome_col = "age",
  validation = "train_test",
  p = 0.7,
  n_perm = 200
)

# Inspect metrics
mlr_results_tt$observed_metrics
mlr_results_tt$train_r2
mlr_results_tt$test_r2
mlr_results_tt$p_value

# Feature importance
mlr_results_tt$feature_importance

# Plots
mlr_results_tt$performance_plot
mlr_results_tt$permutation_importance_diagnostics_plot


# ---- Cross-validation mode ----
mlr_results_cv <- mlr_permutation_eval(
  data = data,
  outcome_col = "age",
  validation = "cv",
  k = 5
)

# Inspect CV metrics
mlr_results_cv$cv_metrics
mlr_results_cv$fold_metrics

# Feature importance from final model
mlr_results_cv$feature_importance

# Plots
mlr_results_cv$performance_plot
mlr_results_cv$permutation_importance_diagnostics_plot

























# elastic net
elasticnet_permutation_eval <- function(data, outcome_col = "age",
                                        validation = c("train_test", "cv"),
                                        p = 0.7, k = 5,
                                        alpha = 0.5, n_perm = 500,
                                        seed = 123, tolerance = 1e-4) {
  validation <- match.arg(validation)
  set.seed(seed)
  
  if (validation == "train_test") {
    split_obj <- make_train_test_split(data, outcome_col, p, seed)
    train_y <- split_obj$train_y
    test_y  <- split_obj$test_y
    
    mm <- make_mm(split_obj$train_x, split_obj$test_x, outcome_col)
    train_x <- mm$train_x
    test_x  <- mm$test_x
    
    cv_fit <- cv.glmnet(x = train_x, y = train_y, alpha = alpha,
                        family = "gaussian", standardize = TRUE)
    best_lambda <- cv_fit$lambda.min
    best_model <- glmnet(x = train_x, y = train_y, alpha = alpha,
                         lambda = best_lambda, family = "gaussian",
                         standardize = TRUE)
    
    train_pred <- as.numeric(predict(best_model, newx = train_x, s = best_lambda))
    test_pred  <- as.numeric(predict(best_model, newx = test_x,  s = best_lambda))
    
    sst <- sum((test_y - mean(test_y))^2)
    sse <- sum((test_y - test_pred)^2)
    obs_rsq <- 1 - sse/sst
    obs_rmse <- sqrt(mean((test_y - test_pred)^2))
    obs_mae  <- mean(abs(test_y - test_pred))
    
    perf <- plot_model_performance(
      train_y = train_y,
      train_pred = train_pred,
      test_y = test_y,
      test_pred = test_pred,
      model_name = paste0("Elastic Net (alpha=", alpha, ")")
    )
    
    coefs <- as.matrix(coef(best_model, s = best_lambda))[-1, , drop = FALSE]
    feat_imp <- data.frame(
      Feature = rownames(coefs),
      Coefficient = as.numeric(coefs),
      stringsAsFactors = FALSE
    )
    feat_imp <- feat_imp[abs(feat_imp$Coefficient) >= tolerance, , drop = FALSE]
    
    perm_rsq <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      y_perm <- sample(train_y)
      cv_fit_perm <- cv.glmnet(x = train_x, y = y_perm, alpha = alpha,
                               family = "gaussian", standardize = TRUE)
      best_lambda_perm <- cv_fit_perm$lambda.min
      y_pred_perm <- as.numeric(predict(cv_fit_perm$glmnet.fit,
                                        newx = test_x, s = best_lambda_perm))
      sst_perm <- sum((test_y - mean(test_y))^2)
      sse_perm <- sum((test_y - y_pred_perm)^2)
      perm_rsq[i] <- 1 - sse_perm/sst_perm
    }
    p_val <- mean(perm_rsq >= obs_rsq)
    
    p1 <- plot_permutation_histogram(
      perm_rsq, obs_rsq, obs_rmse, obs_mae, p_val,
      paste0("Elastic Net (alpha=", alpha, ")")
    )
    if (nrow(feat_imp) > 0) {
      p2 <- ggplot(feat_imp, aes(x = Feature, y = Coefficient)) +
        geom_col(fill = "steelblue", color = 'black', linewidth = 0.05) +
        coord_flip() +
        labs(title = paste0("Elastic Net coefficients (alpha=", alpha, ")"),
             x = "Feature", y = "Coefficient") +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    } else {
      p2 <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No predictors above tolerance") +
        theme_void()
    }
    diagnostics_plot <- plot_residual_diagnostics(
      y_true = test_y,
      y_pred = test_pred,
      model_name = paste0("Elastic Net (alpha=", alpha, ")")
    )
    combined_plot <- ggarrange(p1, p2, diagnostics_plot, ncol = 1, nrow = 3)
    
    return(list(
      observed_metrics = list(R2 = obs_rsq, RMSE = obs_rmse, MAE = obs_mae),
      feature_importance = feat_imp,
      permuted_rsq = perm_rsq,
      p_value = p_val,
      best_lambda = best_lambda,
      alpha = alpha,
      tolerance = tolerance,
      model = best_model,
      permutation_importance_plot = combined_plot,
      performance_plot = perf$plot,
      train_r2 = perf$train_r2,
      test_r2 = perf$test_r2,
      rmse_test = perf$rmse_test,
      mae_test = perf$mae_test
    ))
    
  } else if (validation == "cv") {
    folds <- createFolds(data[[outcome_col]], k = k, list = TRUE)
    all_preds <- numeric(nrow(data))
    
    fold_metrics <- lapply(seq_along(folds), function(i) {
      test_idx <- folds[[i]]
      train_idx <- setdiff(seq_len(nrow(data)), test_idx)
      
      train_fold <- data[train_idx, ]
      test_fold  <- data[test_idx, ]
      
      mm <- make_mm(train_fold[, setdiff(names(train_fold), outcome_col), drop = FALSE],
                    test_fold[, setdiff(names(test_fold), outcome_col), drop = FALSE],
                    outcome_col)
      x_tr <- mm$train_x
      x_te <- mm$test_x
      y_tr <- train_fold[[outcome_col]]
      y_te <- test_fold[[outcome_col]]
      
      cv_fit <- cv.glmnet(x = x_tr, y = y_tr, alpha = alpha,
                          family = "gaussian", standardize = TRUE)
      best_lambda <- cv_fit$lambda.min
      model <- glmnet(x = x_tr, y = y_tr, alpha = alpha,
                      lambda = best_lambda, family = "gaussian",
                      standardize = TRUE)
      y_pred <- as.numeric(predict(model, newx = x_te, s = best_lambda))
      all_preds[test_idx] <<- y_pred
      
      sst <- sum((y_te - mean(y_te))^2)
      sse <- sum((y_te - y_pred)^2)
      r2  <- 1 - sse/sst
      rmse <- sqrt(mean((y_te - y_pred)^2))
      mae  <- mean(abs(y_te - y_pred))
      
      list(R2 = r2, RMSE = rmse, MAE = mae)
    })
    
    mean_r2   <- mean(sapply(fold_metrics, `[[`, "R2"))
    mean_rmse <- mean(sapply(fold_metrics, `[[`, "RMSE"))
    mean_mae  <- mean(sapply(fold_metrics, `[[`, "MAE"))
    
    perf <- plot_model_performance(
      train_y = data[[outcome_col]],
      train_pred = rep(mean(data[[outcome_col]]), nrow(data)),
      test_y = data[[outcome_col]],
      test_pred = all_preds,
      model_name = paste0("Elastic Net (", k, "-fold CV, alpha=", alpha, ")")
    )
    
    full_pred <- sanitize_predictors(data[, setdiff(names(data), outcome_col), drop = FALSE])
    mm_full <- model.matrix(~ . - 1, data = full_pred)
    y_full <- data[[outcome_col]]
    cv_fit_final <- cv.glmnet(x = mm_full, y = y_full, alpha = alpha,
                              family = "gaussian", standardize = TRUE)
    best_lambda_final <- cv_fit_final$lambda.min
    final_model <- glmnet(x = mm_full, y = y_full, alpha = alpha,
                          lambda = best_lambda_final, family = "gaussian",
                          standardize = TRUE)
    
    coefs <- as.matrix(coef(final_model, s = best_lambda_final))[-1, , drop = FALSE]
    feat_imp <- data.frame(
      Feature = rownames(coefs),
      Coefficient = as.numeric(coefs),
      stringsAsFactors = FALSE
    )
    feat_imp <- feat_imp[abs(feat_imp$Coefficient) >= tolerance, , drop = FALSE]
    
    if (nrow(feat_imp) > 0) {
      p2 <- ggplot(feat_imp, aes(x = Feature, y = Coefficient)) +
        geom_col(fill = "steelblue", color = 'black', linewidth = 0.05) +
        coord_flip() +
        labs(title = paste0("Elastic Net coefficients (alpha=", alpha, ")"),
             x = "Feature", y = "Coefficient") +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    } else {
      p2 <- ggplot() +
        annotate("text", x = 0, y = 0, label = "No predictors above tolerance") +
        theme_void()
    }
    
    return(list(
      cv_metrics = list(mean_R2 = mean_r2, mean_RMSE = mean_rmse, mean_MAE = mean_mae),
      fold_metrics = fold_metrics,
      feature_importance = feat_imp,
      diagnostics_plot = p2,
      performance_plot = perf$plot,
      model = final_model,
      best_lambda = best_lambda_final, 
      alpha = alpha, 
      tolerance = tolerance
    ))
  }
}

# ---- Train/test mode ----
enet_results_tt <- elasticnet_permutation_eval(
  data = data,
  outcome_col = "age",
  validation = "train_test",
  p = 0.7,
  alpha = 0.5,
  n_perm = 200,
  tolerance = 1e-4
)

# Inspect metrics
enet_results_tt$observed_metrics
enet_results_tt$train_r2
enet_results_tt$test_r2
enet_results_tt$p_value

# Feature importance (signed coefficients)
enet_results_tt$feature_importance

# Plots
enet_results_tt$permutation_importance_plot   # permutation histogram + coefficients
enet_results_tt$performance_plot              # actual vs predicted scatter


# ---- Cross-validation mode ----
enet_results_cv <- elasticnet_permutation_eval(
  data = data,
  outcome_col = "age",
  validation = "cv",
  k = 5,
  alpha = 0.5,
  n_perm = 200,
  tolerance = 1e-4
)

# Inspect CV metrics
enet_results_cv$cv_metrics        # mean R², RMSE, MAE across folds
enet_results_cv$fold_metrics      # per-fold metrics

# Feature importance from final model
enet_results_cv$feature_importance

# Plots
enet_results_cv$permutation_importance_plot   # coefficients plot
enet_results_cv$performance_plot              # CV predictions vs actual scatter


















rf_permutation_eval <- function(data, outcome_col = "age",
                                validation = c("train_test", "cv"),
                                p = 0.7, k = 5,
                                ntree = 500, n_perm = 500,
                                seed = 123) {
  validation <- match.arg(validation)
  set.seed(seed)
  
  if (validation == "train_test") {
    split_obj <- make_train_test_split(data, outcome_col, p, seed)
    train_x <- sanitize_predictors(split_obj$train_x)
    test_x  <- sanitize_predictors(split_obj$test_x)
    colnames(test_x) <- colnames(train_x)
    train_y <- split_obj$train_y
    test_y  <- split_obj$test_y
    
    tuned <- suppressMessages(
      capture.output(
        tuneRF(
          x = train_x, y = train_y,
          ntreeTry = ntree,
          mtryStart = floor(sqrt(ncol(train_x))),
          stepFactor = 1.5,
          improve = 0.01,
          trace = FALSE,
          doBest = FALSE, 
          plot = FALSE
        )
      )
    )
    tuned <- tryCatch(as.data.frame(tuned), error = function(e) NULL)
    if (!is.null(tuned) && nrow(tuned) > 0 && !all(is.na(tuned$OOBError))) {
      best_mtry <- tuned$mtry[which.min(tuned$OOBError)]
    } else {
      best_mtry <- floor(sqrt(ncol(train_x)))
    }
    
    rf_fit <- randomForest(x = train_x, y = train_y,
                           mtry = best_mtry, ntree = ntree, importance = TRUE)
    train_pred <- predict(rf_fit, train_x)
    test_pred  <- predict(rf_fit, test_x)
    
    sst <- sum((test_y - mean(test_y))^2)
    sse <- sum((test_y - test_pred)^2)
    obs_rsq <- 1 - sse/sst
    obs_rmse <- sqrt(mean((test_y - test_pred)^2))
    obs_mae  <- mean(abs(test_y - test_pred))
    
    feat_imp <- as.data.frame(importance(rf_fit))
    feat_imp$Feature <- rownames(feat_imp)
    rownames(feat_imp) <- NULL
    
    train_df_for_scale <- data.frame(train_x, outcome = train_y)
    names(train_df_for_scale)[ncol(train_df_for_scale)] <- outcome_col
    
    feat_imp_scaled <- scale_importance(
      feat_imp,
      data = train_df_for_scale,
      outcome_col = outcome_col
    )
    
    perm_rsq <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      y_perm <- sample(train_y)
      rf_perm <- randomForest(x = train_x, y = y_perm,
                              mtry = best_mtry, ntree = ntree)
      y_pred_perm <- predict(rf_perm, test_x)
      sst_perm <- sum((test_y - mean(test_y))^2)
      sse_perm <- sum((test_y - y_pred_perm)^2)
      perm_rsq[i] <- 1 - sse_perm/sst_perm
    }
    p_val <- mean(perm_rsq >= obs_rsq)
    
    combined_plot <- assemble_diagnostics_plot(
      model_name = "Random Forest",
      feat_imp_scaled = feat_imp_scaled,
      y_true = test_y,
      y_pred = test_pred,
      include_perm = TRUE,
      perm_rsq = perm_rsq,
      obs_rsq = obs_rsq,
      obs_rmse = obs_rmse,
      obs_mae = obs_mae,
      p_val = p_val
    )
    
    perf <- plot_model_performance(train_y, train_pred, test_y, test_pred, "Random Forest")
    
    return(list(
      observed_metrics = list(R2 = obs_rsq, RMSE = obs_rmse, MAE = obs_mae),
      raw_feature_importance = feat_imp,
      scaled_feature_importance = feat_imp_scaled,
      permuted_rsq = perm_rsq,
      p_value = p_val,
      best_mtry = best_mtry,
      permutation_importance_diagnostics_plot = combined_plot,
      performance_plot = perf$plot,
      train_r2 = perf$train_r2,
      test_r2 = perf$test_r2,
      rmse_test = perf$rmse_test,
      mae_test = perf$mae_test,
      model = rf_fit
    ))
    
  } else if (validation == "cv") {
    folds <- createFolds(data[[outcome_col]], k = k, list = TRUE)
    train_x <- sanitize_predictors(data[ , setdiff(names(data), outcome_col)])
    train_y <- data[[outcome_col]]
    
    all_preds <- numeric(nrow(data))
    fold_metrics <- lapply(seq_along(folds), function(i) {
      test_idx <- folds[[i]]
      train_idx <- setdiff(seq_len(nrow(data)), test_idx)
      
      x_train <- train_x[train_idx, , drop = FALSE]
      y_train <- train_y[train_idx]
      x_test  <- train_x[test_idx, , drop = FALSE]
      y_test  <- train_y[test_idx]
      
      rf_fit <- randomForest(x = x_train, y = y_train,
                             mtry = floor(sqrt(ncol(x_train))),
                             ntree = ntree)
      y_pred <- predict(rf_fit, x_test)
      all_preds[test_idx] <<- y_pred
      
      sst <- sum((y_test - mean(y_test))^2)
      sse <- sum((y_test - y_pred)^2)
      r2  <- 1 - sse/sst
      rmse <- sqrt(mean((y_test - y_pred)^2))
      mae  <- mean(abs(y_test - y_pred))
      
      list(R2 = r2, RMSE = rmse, MAE = mae)
    })
    
    mean_r2   <- mean(sapply(fold_metrics, `[[`, "R2"))
    mean_rmse <- mean(sapply(fold_metrics, `[[`, "RMSE"))
    mean_mae  <- mean(sapply(fold_metrics, `[[`, "MAE"))
    
    perf <- plot_model_performance(
      train_y = train_y,
      train_pred = rep(mean(train_y), length(train_y)),
      test_y = train_y,
      test_pred = all_preds,
      model_name = paste0("Random Forest (", k, "-fold CV)")
    )
    
    rf_final <- randomForest(x = train_x, y = train_y,
                             mtry = floor(sqrt(ncol(train_x))),
                             ntree = ntree, importance = TRUE)
    feat_imp <- as.data.frame(importance(rf_final))
    feat_imp$Feature <- rownames(feat_imp)
    rownames(feat_imp) <- NULL
    
    final_df_for_scale <- data.frame(train_x, outcome = train_y)
    names(final_df_for_scale)[ncol(final_df_for_scale)] <- outcome_col
    
    feat_imp_scaled <- scale_importance(
      feat_imp,
      data = final_df_for_scale,
      outcome_col = outcome_col
    )
    
    combined_plot <- assemble_diagnostics_plot(
      model_name = paste0("Random Forest (", k, "-fold CV)"),
      feat_imp_scaled = feat_imp_scaled,
      y_true = train_y,
      y_pred = all_preds
    )
    
    return(list(
      cv_metrics = list(mean_R2 = mean_r2, mean_RMSE = mean_rmse, mean_MAE = mean_mae),
      fold_metrics = fold_metrics,
      raw_feature_importance = feat_imp,
      scaled_feature_importance = feat_imp_scaled,
      permutation_importance_diagnostics_plot = combined_plot,
      performance_plot = perf$plot,
      model = rf_final
    ))
  }
}

# ---- Train/test mode ----
rf_results_tt <- rf_permutation_eval(
  data = data,
  outcome_col = "age",
  validation = "train_test",
  p = 0.7,        # proportion of data for training
  ntree = 500,
  n_perm = 200,
  seed = 123
)

# Inspect metrics
rf_results_tt$observed_metrics
rf_results_tt$train_r2
rf_results_tt$test_r2
rf_results_tt$p_value
rf_results_tt$best_mtry

# Feature importance
head(rf_results_tt$raw_feature_importance, 5)
head(rf_results_tt$scaled_feature_importance, 5)

# Plots
rf_results_tt$performance_plot
rf_results_tt$permutation_importance_plot


# ---- Cross-validation mode ----
rf_results_cv <- rf_permutation_eval(
  data = data,
  outcome_col = "age",
  validation = "cv",
  k = 5,          # number of folds
  ntree = 500,
  seed = 123
)

# Inspect CV metrics
rf_results_cv$cv_metrics        # mean R², RMSE, MAE across folds
rf_results_cv$fold_metrics      # per-fold metrics

# Feature importance from final model
rf_results_cv$raw_feature_importance
rf_results_cv$scaled_feature_importance

# Plots
rf_results_cv$performance_plot
rf_results_cv$permutation_importance_plot
