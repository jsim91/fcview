# App used to explore FCSimple results
# run by: shiny::runApp(appDir = 'app.R_directory')
# where 'app.R_directory' is the location of this script, app.R

# ---- Packages ----
suppressPackageStartupMessages({
  suppressWarnings({
    library(shiny)
    library(shinyWidgets)
    library(plotly)
    library(ggplot2)
    library(ggrastr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ComplexHeatmap)
    library(circlize)
    library(viridis)
    library(scales)
    library(glmnet)
    library(Boruta)
    library(colourpicker)
    library(caret) # train(), trainControl, resampling
    library(nnet) # multinom() for multi-class logistic regression
    library(pROC)
    library(yardstick)
    library(rsample) # vfold_cv(), initial_split (if you want tidy resampling)
    library(boot) # bootstrapping CIs for AUC
    library(rlang) # tidy evaluation for !!!syms
    library(randomForest)
    library(gbm)
    library(xgboost)
    library(broom)
    library(pdp)
    library(iml)
  })
})

options(shiny.maxRequestSize = 100000 * 1024^2)

# ---- Helpers & validation ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

pct_clip <- function(x, p = c(0.01, 0.99)) {
  q <- quantile(x, probs = p, na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

make_safe_names <- function(df, predictors) {
  safe <- make.names(colnames(df))
  name_map <- setNames(colnames(df), safe)  # safe → original
  colnames(df) <- safe
  predictors_safe <- make.names(predictors)
  list(df = df, predictors = predictors_safe, map = name_map)
}

map_safe_to_original <- function(df, map) {
  if (is.null(df) || is.null(map)) return(df)
  if ("Feature" %in% names(df)) {
    df$Feature <- dplyr::recode(df$Feature, !!!map)
  }
  if ("Predictor" %in% names(df)) {
    df$Predictor <- dplyr::recode(df$Predictor, !!!map)
  }
  df
}

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
    theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5))
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
  
  pl <- ggplot(long, aes(x = Feature, y = Importance, fill = Type)) +
    geom_col(position = position_dodge(width = 0.7)) +
    coord_flip() +
    scale_fill_manual(values = c("RawImportance" = "steelblue",
                                 "ScaledImportance" = "red3"),
                      name = "Importance type") +
    labs(title = paste0(model_label, " feature importance: raw vs scaled"),
         subtitle = paste("Ordered by:", order_label),
         x = "Feature", y = "Importance") +
    theme_bw(base_size = 15) + theme(plot.title = element_text(hjust = 0.5))
  nfeats <- unique(long$Feature)
  return(list(pl, nfeats))
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
    geom_point(shape = 21, fill = "grey40", color = "black", stroke = 0.2, size = 4, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 0.9) +
    geom_hline(yintercept = mean(test_y), color = "blue", linetype = "dotted", linewidth = 0.9) +
    annotate("text", x = -Inf, y = Inf, label = metrics_text,
             hjust = -0.05, vjust = 1.05, size = 6) +
    labs(title = paste0(model_name, ": Actual vs Predicted (Test Set)"),
         x = "Actual values", y = "Predicted values") +
    theme_bw(base_size = 16) +
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
    relheight <- p2[[2]]
    combined_plot <- ggpubr::ggarrange(p1, p2[[1]], diagnostics_plot, ncol = 1, nrow = 3, heights = c(1,3,2))
    
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

mlr_permutation_eval <- function(data, outcome_col = "age",
                                 validation = c("train_test", "cv"),
                                 p = 0.7, k = 5,
                                 n_perm = 500, seed = 123,
                                 tolerance = 1e-4) {
  validation <- match.arg(validation)
  set.seed(seed)
  
  # Build design matrix
  X <- model.matrix(~ . - 1, data = data[, setdiff(names(data), outcome_col), drop = FALSE])
  storage.mode(X) <- "double"
  y <- data[[outcome_col]]
  
  if (validation == "train_test") {
    idx <- caret::createDataPartition(y, p = p, list = FALSE)
    train_x <- X[idx, , drop = FALSE]
    test_x  <- X[-idx, , drop = FALSE]
    train_y <- y[idx]
    test_y  <- y[-idx]
    
    ctrl <- caret::trainControl(verboseIter = FALSE)
    mlr_model <- caret::train(x = train_x, y = train_y, method = "lm", trControl = ctrl)
    
    train_pred <- predict(mlr_model, newdata = train_x)
    test_pred  <- predict(mlr_model, newdata = test_x)
    
    sst <- sum((test_y - mean(test_y))^2)
    sse <- sum((test_y - test_pred)^2)
    obs_rsq <- 1 - sse/sst
    obs_rmse <- sqrt(mean((test_y - test_pred)^2))
    obs_mae  <- mean(abs(test_y - test_pred))
    
    coefs <- coef(mlr_model$finalModel)[-1]
    feat_imp <- data.frame(Feature = names(coefs), Coefficient = as.numeric(coefs))
    feat_imp_scaled <- scale_importance(feat_imp,
                                        data = data.frame(train_x, outcome = train_y),
                                        outcome_col = "outcome")
    
    # Permutation test
    perm_rsq <- replicate(n_perm, {
      y_perm <- sample(train_y)
      perm_model <- caret::train(x = train_x, y = y_perm, method = "lm", trControl = ctrl)
      y_pred_perm <- predict(perm_model, newdata = test_x)
      sst_perm <- sum((test_y - mean(test_y))^2)
      sse_perm <- sum((test_y - y_pred_perm)^2)
      1 - sse_perm/sst_perm
    })
    p_val <- mean(perm_rsq >= obs_rsq)
    
    perf <- plot_model_performance(train_y, train_pred, test_y, test_pred, "MLR")
    
    return(list(
      observed_metrics = list(R2 = obs_rsq, RMSE = obs_rmse, MAE = obs_mae),
      feature_importance = feat_imp_scaled,
      permuted_rsq = perm_rsq,
      p_value = p_val,
      model = mlr_model$finalModel,
      permutation_importance_diagnostics_plot = NULL,
      performance_plot = perf$plot,
      train_r2 = perf$train_r2,
      test_r2 = perf$test_r2,
      rmse_test = perf$rmse_test,
      mae_test = perf$mae_test
    ))
  }
  
  # CV branch (similar, but loop over folds)
  folds <- caret::createFolds(y, k = k, list = TRUE)
  all_preds <- numeric(length(y))
  fold_metrics <- lapply(folds, function(test_idx) {
    train_idx <- setdiff(seq_along(y), test_idx)
    x_tr <- X[train_idx, , drop = FALSE]
    y_tr <- y[train_idx]
    x_te <- X[test_idx, , drop = FALSE]
    y_te <- y[test_idx]
    
    model <- caret::train(x = x_tr, y = y_tr, method = "lm")
    y_pred <- predict(model, newdata = x_te)
    all_preds[test_idx] <<- y_pred
    
    sst <- sum((y_te - mean(y_te))^2)
    sse <- sum((y_te - y_pred)^2)
    list(R2 = 1 - sse/sst,
         RMSE = sqrt(mean((y_te - y_pred)^2)),
         MAE = mean(abs(y_te - y_pred)))
  })
  
  mean_r2   <- mean(sapply(fold_metrics, `[[`, "R2"))
  mean_rmse <- mean(sapply(fold_metrics, `[[`, "RMSE"))
  mean_mae  <- mean(sapply(fold_metrics, `[[`, "MAE"))
  
  final_model <- caret::train(x = X, y = y, method = "lm")
  coefs <- coef(final_model$finalModel)[-1]
  feat_imp <- data.frame(Feature = names(coefs), Coefficient = as.numeric(coefs))
  feat_imp_scaled <- scale_importance(feat_imp,
                                      data = data.frame(X, outcome = y),
                                      outcome_col = "outcome")
  
  return(list(
    cv_metrics = list(mean_R2 = mean_r2, mean_RMSE = mean_rmse, mean_MAE = mean_mae),
    fold_metrics = fold_metrics,
    feature_importance = feat_imp_scaled,
    performance_plot = NULL,
    model = final_model$finalModel
  ))
}

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

normalize_regression_result <- function(res, validation) {
  # Default placeholders
  out <- list(
    metrics = NULL,
    feature_importance = res$feature_importance %||% res$scaled_feature_importance %||% NULL,
    performance_plot = res$performance_plot %||% NULL,
    diagnostics_plot = res$permutation_importance_diagnostics_plot %||% res$diagnostics_plot %||% NULL,
    model = res$model %||% NULL
  )
  
  if (validation == "train_test") {
    out$metrics <- list(
      Train_R2 = res$train_r2 %||% NA,
      Test_R2  = res$test_r2 %||% NA,
      RMSE     = res$rmse_test %||% NA,
      MAE      = res$mae_test %||% NA,
      p_value  = res$p_value %||% NA
    )
  } else if (validation == "cv") {
    cm <- res$cv_metrics %||% list()
    out$metrics <- list(
      Mean_R2   = cm$mean_R2 %||% NA,
      Mean_RMSE = cm$mean_RMSE %||% NA,
      Mean_MAE  = cm$mean_MAE %||% NA
    )
  }
  
  return(out)
}

run_regression_model <- function(model_type, df, outcome, predictors,
                                 validation = "train_test", p = 0.7, k = 5,
                                 alpha = 0.5, n_perm = 100, seed = 123) {
  df <- df[, c(outcome, predictors), drop = FALSE]
  df <- na.omit(df)
  
  # Sanitize predictors
  X <- model.matrix(~ . - 1, data = df[, predictors, drop = FALSE])
  storage.mode(X) <- "double"
  y <- df[[outcome]]
  
  safe_names <- colnames(X)
  name_map <- setNames(predictors, safe_names)
  df_safe <- data.frame(outcome = y, X, check.names = FALSE)
  names(df_safe)[1] <- outcome
  
  raw <- switch(model_type,
                "xgbTree" = xgb_permutation_eval(df_safe, outcome_col = outcome,
                                                 validation = validation, p = p, k = k,
                                                 n_perm = n_perm, seed = seed),
                "lm" = mlr_permutation_eval(df_safe, outcome_col = outcome,
                                            validation = validation, p = p, k = k,
                                            n_perm = n_perm, seed = seed),
                "glmnet" = elasticnet_permutation_eval(df_safe, outcome_col = outcome,
                                                       validation = validation, p = p, k = k,
                                                       alpha = alpha, n_perm = n_perm, seed = seed),
                "rf" = rf_permutation_eval(df_safe, outcome_col = outcome,
                                           validation = validation, p = p, k = k,
                                           n_perm = n_perm, seed = seed),
                stop("Unknown model type")
  )
  
  raw$name_map <- name_map
  normalize_regression_result(raw, validation)
}

# Assert that a given metadata frame is sample-level
assert_sample_level <- function(meta_df, tab_name = "Unknown") {
  if (is.null(meta_df)) {
    showNotification(paste0("[", tab_name, "] metadata is NULL."), type = "error")
    return(FALSE)
  }
  # Heuristic: sample-level metadata should have one row per patient_ID
  if ("patient_ID" %in% colnames(meta_df)) {
    dup_ids <- meta_df$patient_ID[duplicated(meta_df$patient_ID)]
    if (length(dup_ids) > 0) {
      showNotification(
        paste0("[", tab_name, "] Warning: duplicate patient_IDs detected: ",
               paste(unique(dup_ids), collapse = ", ")),
        type = "warning", duration = 10
      )
    }
  }
  TRUE
}

resetResults <- function(resultReactive, cacheReactive = NULL) {
  resultReactive(NULL)
  if (!is.null(cacheReactive)) cacheReactive(NULL)
}

align_metadata_abundance <- function(metadata, abundance, notify = NULL) {
  # Extract patient_ID from abundance rownames
  patient_ids <- gsub(
    pattern = "_[0-9]+\\-[A-Za-z]+\\-[0-9]+.*$",
    replacement = "",
    x = rownames(abundance)
  )
  
  abund_df <- as.data.frame(abundance)
  abund_df$patient_ID <- patient_ids
  
  print('dim metadata: ')
  print(dim(metadata))
  print('dim abund_df: ')
  print(dim(abund_df))
  print('head metadata: ')
  print(head(metadata,2))
  print('head abund_df: ')
  print(head(abund_df,2))
  
  # Duplicate checks
  if (anyDuplicated(metadata$patient_ID)) {
    msg <- "Duplicate patient_IDs found in metadata. Consider deduplicating with distinct()."
    warning(msg)
    if (!is.null(notify)) notify(msg, type = "warning")
  }
  if (anyDuplicated(abund_df$patient_ID)) {
    msg <- "Duplicate patient_IDs found in abundance rownames. Check your input files."
    warning(msg)
    if (!is.null(notify)) notify(msg, type = "warning")
  }

  # Mismatch checks
  meta_ids <- unique(metadata$patient_ID)
  abund_ids <- unique(abund_df$patient_ID)

  missing_in_abund <- setdiff(meta_ids, abund_ids)
  missing_in_meta  <- setdiff(abund_ids, meta_ids)

  if (length(missing_in_abund) > 0) {
    msg <- paste("These patient_IDs are in metadata but not in abundance:",
                 paste(missing_in_abund, collapse = ", "))
    warning(msg)
    if (!is.null(notify)) notify(msg, type = "warning")
  }
  if (length(missing_in_meta) > 0) {
    msg <- paste("These patient_IDs are in abundance but not in metadata:",
                 paste(missing_in_meta, collapse = ", "))
    warning(msg)
    if (!is.null(notify)) notify(msg, type = "warning")
  }

  # Merge with metadata (metadata is the anchor)
  merged <- dplyr::left_join(metadata, abund_df, by = "patient_ID")
  print('head merged: ')
  print(head(merged,5))
  
  return(merged)
}

expand_predictors <- function(meta_sample, abundance_sample,
                              outcome, predictors,
                              cluster_subset = NULL,
                              mode = c("regression","classification","fs"),
                              encode_factors = FALSE,
                              notify = NULL) {
  mode <- match.arg(mode)
  
  # Split meta vs pseudo-predictor
  pred_meta <- setdiff(intersect(predictors, colnames(meta_sample)), "leiden_cluster")
  
  # Handle clusters
  cluster_predictors <- character(0)
  if ("leiden_cluster" %in% predictors) {
    if (is.null(abundance_sample)) {
      stop("Abundance matrix not available, but leiden_cluster was selected.")
    }
    all_clusters <- colnames(abundance_sample)
    cluster_predictors <- if (!is.null(cluster_subset) && length(cluster_subset) > 0) {
      intersect(cluster_subset, all_clusters)
    } else {
      all_clusters
    }
  }
  
  predictors_final <- c(pred_meta, cluster_predictors)
  if (length(predictors_final) == 0) {
    stop("Select at least one predictor.")
  }
  
  # Subset metadata
  keep_cols <- unique(c("patient_ID", outcome, pred_meta))
  meta_sub <- meta_sample[, keep_cols, drop = FALSE]
  
  # If no cluster predictors, skip abundance merge
  if (length(cluster_predictors) == 0) {
    merged <- meta_sub
  } else {
    abund_sub <- abundance_sample[, cluster_predictors, drop = FALSE]
    merged <- align_metadata_abundance(metadata = meta_sub,
                                       abundance = abund_sub,
                                       notify = notify)
  }
  
  # Filter NAs
  merged <- merged[!is.na(merged[[outcome]]), , drop = FALSE]
  
  # Defensive check
  missing <- setdiff(c(outcome, predictors_final), colnames(merged))
  if (length(missing) > 0) {
    stop("These selected columns are not in the merged data: ",
         paste(missing, collapse = ", "))
  }
  
  # Mode-specific encoding (unchanged from your current version)
  if (mode == "regression") {
    X <- merged[, predictors_final, drop = FALSE]
    
    if (isTRUE(encode_factors)) {
      # Build design matrix (safe names)
      mm <- model.matrix(~ . - 1, data = X)
      storage.mode(mm) <- "double"
      colnames(mm) <- gsub("leiden_cluster", "leiden_cluster:", colnames(mm))
      
      # Build mapping safe → original
      # predictors_final is the actual expanded set (meta + cluster columns)
      if (length(predictors_final) != ncol(mm)) {
        # Defensive: if model.matrix expanded factors into multiple dummies,
        # replicate the original predictor name across those new columns
        expanded_map <- rep(predictors_final, times = sapply(X, function(col) {
          if (is.factor(col) || is.character(col)) {
            nlevels(as.factor(col))
          } else {
            1
          }
        }))
        name_map <- setNames(expanded_map, colnames(mm))
      } else {
        name_map <- setNames(predictors_final, colnames(mm))
      }
      
      merged_out <- data.frame(
        patient_ID = merged$patient_ID,
        merged[[outcome]],
        mm,
        check.names = FALSE
      )
      names(merged_out)[2] <- outcome  # restore outcome column name
      
      return(list(
        df = merged_out,
        predictors = colnames(mm),
        map = name_map
      ))
    }
    
    # No encoding: assume numeric only
    return(list(
      df = merged[, c("patient_ID", outcome, predictors_final), drop = FALSE],
      predictors = predictors_final,
      map = setNames(predictors_final, predictors_final)
    ))
  }
  
  # Classification / FS
  return(list(df = merged[, c("patient_ID", outcome, predictors_final), drop = FALSE],
              predictors = predictors_final,
              map = setNames(predictors_final, predictors_final)))
}

clean_dummy_names <- function(nms) {
  gsub(pattern = 'leiden_cluster', replacement = 'leiden_cluster:', x = nms)
}

spearman_test <- function(df, freq_col = "freq", cont_var) {
  ok <- complete.cases(df[[freq_col]], df[[cont_var]])
  if (!any(ok)) return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
  ct <- suppressWarnings(cor.test(df[[freq_col]][ok], df[[cont_var]][ok], method = "spearman"))
  data.frame(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

compute_multiROC <- function(preds) {
  if (!is.factor(preds$obs)) preds$obs <- factor(preds$obs)
  class_levels <- levels(preds$obs)
  if (length(class_levels) < 2) stop("Need ≥2 classes for ROC")
  
  # Build truth and prob data.frames
  truth_df <- as.data.frame(
    sapply(class_levels, function(cls) as.integer(preds$obs == cls), simplify = "matrix")
  )
  names(truth_df) <- paste0(class_levels, "_true")
  
  prob_df <- preds[, paste0(".pred_", class_levels), drop = FALSE]
  names(prob_df) <- paste0(class_levels, "_pred")
  
  # Combine
  multi_df <- cbind(truth_df, prob_df)
  
  # Force to plain numeric matrix
  multi_mat <- as.matrix(multi_df)
  storage.mode(multi_mat) <- "numeric"
  
  message("compute_multiROC: multi_mat head")
  print(utils::head(multi_mat))
  
  # Run multiROC on the matrix
  roc_res <- multiROC::multi_roc(multi_mat, force_diag = TRUE)
  attr(roc_res, "multi_df") <- multi_df
  roc_res
}

# Validate minimal structure
validateInput <- function(obj, id_col = NULL) {
  required <- c("data", "source", "metadata", "leiden")
  miss <- setdiff(required, names(obj))
  if (length(miss)) stop("Missing required elements: ", paste(miss, collapse = ", "))
  
  if (!is.matrix(obj$data)) stop("data must be a matrix (cells × features).")
  n <- nrow(obj$data)
  
  if (length(obj$source) != n) stop("source length must equal nrow(data).")
  if (!is.data.frame(obj$metadata)) stop("metadata must be a data.frame.")
  if (!is.list(obj$leiden) || is.null(obj$leiden$clusters)) stop("leiden must be a list with 'clusters'.")
  if (length(obj$leiden$clusters) != n) stop("leiden$clusters length must equal nrow(data).")
  
  # Optional elements: if present, check shape
  if ("umap" %in% names(obj)) {
    if (!is.data.frame(obj$umap$coordinates) || nrow(obj$umap$coordinates) != n)
      stop("umap$coordinates must align with cells.")
  }
  if ("tsne" %in% names(obj)) {
    if (!is.data.frame(obj$tsne$coordinates) || nrow(obj$tsne$coordinates) != n)
      stop("tsne$coordinates must align with cells.")
  }
  if ("run_date" %in% names(obj)) {
    if (length(obj$run_date) != n) stop("run_date must be length nrow(data).")
  }
  
  # id_col presence check (if provided)
  if (!is.null(id_col) && !(id_col %in% colnames(obj$metadata))) {
    stop("Selected ID column '", id_col, "' not found in metadata.")
  }
  invisible(TRUE)
}

# Presence checks
hasUMAP <- function(obj) "umap" %in% names(obj) && !is.null(obj$umap$coordinates)
hasTSNE <- function(obj) "tsne" %in% names(obj) && !is.null(obj$tsne$coordinates)
hasHeatmap <- function(obj) "leiden_heatmap" %in% names(obj) && !is.null(obj$leiden_heatmap$heatmap_tile_data)
hasRunDate <- function(obj) "run_date" %in% names(obj) && length(obj$run_date) == nrow(obj$data)
hasClusterMapping <- function(obj) "cluster_mapping" %in% names(obj)

# Auto-detect metadata ID column that matches source
guess_id_col <- function(metadata, source_vec) {
  # Prioritized candidates by name
  candidates <- intersect(tolower(colnames(metadata)),
                          c("patientid","patient_id","patient","source","subject","sample","id"))
  if (length(candidates)) {
    hit <- candidates[1]
    return(colnames(metadata)[tolower(colnames(metadata)) == hit][1])
  }
  # Fallback: choose column with max overlap
  overlaps <- sapply(metadata, function(col) {
    if (!is.atomic(col)) return(0)
    length(intersect(as.character(col), as.character(source_vec)))
  })
  if (all(overlaps == 0)) return(NULL)
  colnames(metadata)[which.max(overlaps)]
}

build_design_matrix <- function(df) {
  # model.matrix will one-hot encode factors/characters; drop intercept
  mm <- model.matrix(~ . - 1, data = df)
  # Ensure a plain numeric matrix
  storage.mode(mm) <- "double"
  mm
}

# ---- Embedding module UI ----
EmbeddingUI <- function(id, title = "UMAP") {
  ns <- NS(id)
  tagList(
    h3(title),
    fluidRow(
      column(
        3,
        pickerInput(ns("color_by"), "Color by", choices = NULL),
        checkboxInput(ns("show_labels"), "Show cluster labels", value = FALSE),
        pickerInput(ns("split_by"), "Facet by", choices = NULL,
                    options = list(`none-selected-text` = "None")),
        uiOutput(ns("split_levels_ui")),
        selectInput(
          ns("max_facets"), "Facet columns",
          choices = c(1, 2, 3, 4), selected = 2
        ),
        actionButton(ns("plot_facets"), "Plot facets"),
        hr(),
        # Moved export button here:
        downloadButton(ns("export_embed_pdf"), "Export embedding as PDF"),
        hr()
      ),
      column(
        9,
        plotlyOutput(ns("embed_plot"), height = "650px")
      )
    )
  )
}

# ---- Embedding module server ----
EmbeddingServer <- function(id, embedding_name, coords, expr, meta_cell, clusters, cluster_map,
                            active_tab) {
  moduleServer(id, function(input, output, session) {
    message(sprintf("EmbeddingServer %s started", embedding_name))
    message(sprintf("coords NULL? %s | expr NULL? %s | meta_cell NULL? %s",
                    is.null(coords), is.null(expr), is.null(meta_cell)))
    
    # One-time picker initialization in a reactive context
    initialized <- FALSE
    plot_cache_gg <- reactiveVal(NULL)
    facet_cols_saved <- reactiveVal(2)
    
    observeEvent(list(expr(), meta_cell(), clusters()), ignoreInit = FALSE, {
      if (initialized) return()
      expr_val     <- expr()
      meta_val     <- meta_cell()
      clusters_val <- clusters()
      
      req(!is.null(expr_val), !is.null(meta_val))
      
      # Add leiden_cluster column as factor if missing
      if (!"leiden_cluster" %in% colnames(meta_val) && !is.null(clusters_val$assignments)) {
        meta_val$leiden_cluster <- factor(clusters_val$assignments)
      }
      
      numeric_markers <- colnames(expr_val)
      meta_cols <- setdiff(colnames(meta_val), c(".cell"))
      
      # Ensure leiden_cluster is present and comes right after markers
      if (!"leiden_cluster" %in% meta_cols && "leiden_cluster" %in% colnames(meta_val)) {
        meta_cols <- c("leiden_cluster", setdiff(meta_cols, "leiden_cluster"))
      }
      
      # Sort metadata portion alphabetically (excluding leiden_cluster)
      meta_cols_sorted <- sort(setdiff(meta_cols, "leiden_cluster"))
      
      # Final order: markers → leiden_cluster → sorted metadata
      ordered_choices <- c(numeric_markers, "leiden_cluster", meta_cols_sorted)
      
      # Continuous and categorical choices from metadata
      cont_choices <- sort(meta_cols[sapply(meta_val[meta_cols], is.numeric)])
      factor_cols  <- meta_cols[sapply(meta_val[meta_cols], is.factor)]
      char_cols    <- meta_cols[sapply(meta_val[meta_cols], is.character)]
      categorical_choices <- sort(c(factor_cols, char_cols))
      
      # Default unit_var
      unit_default <- if ("PatientID" %in% meta_cols) "PatientID" else meta_cols[1]
      
      updatePickerInput(session, "color_by",
                        choices = ordered_choices,
                        selected = if (length(numeric_markers)) numeric_markers[1] else ordered_choices[1])
      updatePickerInput(session, "split_by", choices = c("", categorical_choices), selected = "")
      updatePickerInput(session, "group_var", choices = c("", sort(meta_cols)), selected = "")
      updatePickerInput(session, "cont_var", choices = c("", cont_choices), selected = "")
      updatePickerInput(session, "unit_var", choices = sort(meta_cols), selected = unit_default)
      
      message(sprintf("Picker inputs initialized for %s (reactive observer)", embedding_name))
      initialized <<- TRUE
    })
    
    observeEvent(input$split_by, {
      split_var <- input$split_by
      meta_val <- meta_cell()
      if (!is.null(split_var) && nzchar(split_var) && split_var %in% colnames(meta_val)) {
        levs <- sort(unique(as.character(meta_val[[split_var]])))
        updatePickerInput(session, "split_levels", choices = levs, selected = levs)
      }
    })
    
    output$split_levels_ui <- renderUI({
      req(input$split_by, nzchar(input$split_by))
      pickerInput(ns("split_levels"), "Facet categories to show",
                  choices = NULL, multiple = TRUE)
    })
    
    ns <- session$ns
    
    # Build plotting dataframe defensively, with downsampling for rendering
    df <- reactive({
      coords_val       <- coords()
      expr_val         <- expr()
      meta_val         <- meta_cell()
      clusters_val     <- clusters()
      cluster_map_val  <- cluster_map()
      
      req(!is.null(coords_val))
      dd <- as.data.frame(coords_val)
      
      coords_full <- as.data.frame(coords())
      names(coords_full)[1:2] <- c("UMAP1", "UMAP2")
      coords_full$cluster <- factor(clusters()$assignments)  # factor for discrete colors
      
      color_text_add <- coords_full %>%
        group_by(cluster) %>%
        summarise(
          UMAP1 = mean(UMAP1, na.rm = TRUE),
          UMAP2 = mean(UMAP2, na.rm = TRUE),
          .groups = "drop"
        )
      
      if (ncol(dd) < 2 || nrow(dd) == 0) {
        message(sprintf("[%s] df(): coords invalid — ncol=%s nrow=%s",
                        embedding_name, ncol(dd), nrow(dd)))
        return(tibble::tibble(x = numeric(0), y = numeric(0), .cell = integer(0)))
      }
      
      names(dd)[1:2] <- c("x", "y")
      
      # IMPORTANT: .cell = original row number BEFORE downsampling
      dd$.cell <- seq_len(nrow(dd))
      
      # Add cluster and celltype columns
      if (!is.null(clusters_val) && !is.null(clusters_val$assignments)) {
        if (length(clusters_val$assignments) == nrow(dd)) {
          dd$cluster <- clusters_val$assignments
        } else {
          message(sprintf("[%s] df(): cluster length=%s != nrow(coords)=%s",
                          embedding_name, length(clusters_val$assignments), nrow(dd)))
          dd$cluster <- NA_integer_
        }
      } else {
        dd$cluster <- NA_integer_
      }
      
      if (!is.null(cluster_map_val) &&
          all(c("cluster", "celltype") %in% names(cluster_map_val))) {
        dd$celltype <- as.character(cluster_map_val$celltype[
          match(dd$cluster, cluster_map_val$cluster)
        ])
      } else {
        dd$celltype <- as.character(dd$cluster)
      }
      
      # ---- Downsample for plotting only ----
      if (nrow(dd) > 100000) {
        set.seed(123)  # reproducible sampling
        keep_idx <- sample.int(nrow(dd), 100000)
        dd <- dd[keep_idx, , drop = FALSE]
        message(sprintf("[%s] df(): downsampled to %d rows for plotting",
                        embedding_name, nrow(dd)))
      }
      
      tibble::as_tibble(dd)
    })
    
    current_sel <- reactiveVal(integer(0))
    
    plot_cache_base <- reactiveVal(NULL)  # base plot
    plot_cache      <- reactiveVal(NULL)  # final plot with overlays
    
    # Helper: clip values to reference distribution percentiles (1%–99%)
    clip_to_ref <- function(values, ref, probs = c(0.01, 0.99)) {
      qs <- stats::quantile(ref, probs = probs, na.rm = TRUE)
      pmin(pmax(values, qs[1]), qs[2])
    }
    
    # Add split_by to triggers for this observer
    observeEvent(
      list(df(), expr(), meta_cell(), input$color_by, input$plot_facets),
      {
        expr_val <- expr()
        meta_val <- meta_cell()
        req(expr_val, meta_val)
        
        numeric_markers <- colnames(expr_val)
        meta_cols <- colnames(meta_val)
        valid_cols <- c(numeric_markers, meta_cols)
        
        color_by <- input$color_by
        if (is.null(color_by) || !(color_by %in% valid_cols)) {
          color_by <- if (length(numeric_markers)) numeric_markers[1] else meta_cols[1]
        }
        
        dd <- df()
        if (nrow(dd) == 0 || all(is.na(dd$x)) || all(is.na(dd$y))) {
          plot_cache_base(
            plotly_empty(type = "scatter", mode = "markers", source = ns("embed")) %>%
              layout(
                xaxis = list(title = paste0(embedding_name, " 1")),
                yaxis = list(title = paste0(embedding_name, " 2"))
              )
          )
          return()
        }
        
        # --- Always add the colour column to dd ---
        if (color_by %in% numeric_markers) {
          dd$.color_val <- as.numeric(expr_val[dd$.cell, color_by])
        } else if (color_by %in% meta_cols) {
          dd$.color_val <- meta_val[[color_by]][dd$.cell]
        } else {
          dd$.color_val <- NA
        }
        
        # Handle faceting
        split_var <- input$split_by
        if (!is.null(split_var) && nzchar(split_var)) {
          dd[[split_var]] <- meta_val[[split_var]][dd$.cell]
          shown_levels <- input$split_levels
          if (!is.null(shown_levels) && length(shown_levels) > 0) {
            dd <- dd[dd[[split_var]] %in% shown_levels, , drop = FALSE]
          }
          dd <- dd[!is.na(dd[[split_var]]), , drop = FALSE]
          dd[[split_var]] <- droplevels(as.factor(dd[[split_var]]))
          if (nrow(dd) == 0) {
            plot_cache_base(NULL); plot_cache(NULL)
            showNotification("No data points available for the selected facet categories.", type = "error")
            return()
          }
          
          # Balanced downsampling
          max_total <- input$max_cells_upload %||% 100000
          facet_counts <- table(dd[[split_var]])
          n_facets <- length(facet_counts)
          target_per_facet <- floor(max_total / n_facets)
          target_per_facet <- min(target_per_facet, min(facet_counts))
          if (target_per_facet < 1) {
            showNotification("Selected facets have too few points to plot.", type = "error")
            plot_cache_base(NULL); plot_cache(NULL)
            return()
          }
          set.seed(123)
          dd <- dd %>%
            group_by(.data[[split_var]]) %>%
            slice_sample(n = target_per_facet) %>%
            ungroup()
          
          # Build faceted ggplot
          gg <- ggplot(dd, aes(x = x, y = y, color = .data[[".color_val"]])) +
            ggrastr::geom_point_rast(size = 0.25, alpha = 0.5) +
            guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
            facet_wrap(as.formula(paste("~", split_var)), ncol = as.numeric(input$max_facets)) +
            (if (is.numeric(dd$.color_val)) scale_color_viridis_c() else scale_color_viridis_d()) +
            theme_minimal() +
            guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
            theme(legend.position = "right") +
            labs(
              x = paste0(embedding_name, " 1"),
              y = paste0(embedding_name, " 2"),
              color = color_by
            )
          
          plot_cache_gg(gg)
          p_base <- ggplotly(gg, tooltip = "text") %>% layout(dragmode = "pan")
          
        } else {
          # Single-panel ggplot for export
          gg <- ggplot(dd, aes(x = x, y = y, color = .data[[".color_val"]])) +
            geom_point(size = 0.25, alpha = 0.25) +
            (if (is.numeric(dd$.color_val)) scale_color_viridis_c() else scale_color_viridis_d()) +
            theme_minimal() + 
            guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) + 
            labs(
              x = paste0(embedding_name, " 1"),
              y = paste0(embedding_name, " 2"),
              color = color_by
            )
          plot_cache_gg(gg)
          
          # Compute colours for plotly
          if (is.numeric(dd$.color_val)) {
            rng <- range(dd$.color_val, na.rm = TRUE)
            if (!is.finite(rng[1]) || !is.finite(rng[2])) rng <- c(0, 1)
            cols <- scales::col_numeric(
              viridis::viridis(256),
              domain = rng
            )(dd$.color_val)
          } else {
            vals <- as.character(dd$.color_val)
            levs <- sort(unique(vals))
            pal <- viridis::viridis(length(levs))
            cols <- setNames(pal, levs)[vals]
          }
          
          p_base <- plot_ly(
            data = dd,
            x = ~x,
            y = ~y,
            type = "scatter",
            mode = "markers",
            marker = list(color = cols, size = 3),
            source = ns("embed"),
            customdata = ~.cell
          ) %>% layout(
            xaxis = list(title = paste0(embedding_name, " 1")),
            yaxis = list(title = paste0(embedding_name, " 2")),
            dragmode = "zoom"
          )
        }
        
        plot_cache_base(p_base)
        plot_cache(p_base)
      }
    )
    
    # Update overlays without rebuilding points
    observeEvent(input$show_labels, {
      p <- plot_cache_base()
      req(p)
      
      # Add labels if toggled on
      if (isTRUE(input$show_labels)) {
        coords_full <- as.data.frame(coords())
        names(coords_full)[1:2] <- c("x", "y")
        coords_full$cluster <- factor(clusters()$assignments)
        
        label_df <- coords_full %>%
          group_by(cluster) %>%
          summarise(
            x = mean(x, na.rm = TRUE),
            y = mean(y, na.rm = TRUE),
            .groups = "drop"
          )
        
        annots <- lapply(seq_len(nrow(label_df)), function(i) {
          list(
            x = label_df$x[i],
            y = label_df$y[i],
            xref = "x",
            yref = "y",
            text = as.character(label_df$cluster[i]),
            showarrow = FALSE,
            xanchor = "center",
            yanchor = "middle",
            align = "center",
            font = list(color = "black", size = 18),
            bgcolor = "rgba(255,255,255,0.85)",
            bordercolor = "rgba(0,0,0,0)",
            borderpad = 2,
            opacity = 1
          )
        })
        
        p <- p %>% layout(annotations = annots)
      }
      
      plot_cache(p)
    })
    
    output$embed_plot <- renderPlotly({
      req(plot_cache())
      plot_cache()
    })
    
    output$export_embed_pdf <- downloadHandler(
      filename = function() {
        embed_name <- tolower(embedding_name)  # "umap" or "tsne"
        color_var  <- tolower(input$color_by %||% "color")
        
        if (!is.null(input$split_by) && nzchar(input$split_by)) {
          facet_var <- tolower(input$split_by)
          paste0(embed_name, "_", facet_var, "_", color_var, ".pdf")
        } else {
          paste0(embed_name, "_", color_var, ".pdf")
        }
      },
      content = function(file) {
        gg <- plot_cache_gg()
        if (is.null(gg) || !inherits(gg, "ggplot")) {
          showNotification("No ggplot available to export. Try re‑plotting before exporting.", type = "error")
          return()
        }
        
        if (!is.null(input$split_by) && nzchar(input$split_by)) {
          # Faceted → dynamic sizing
          n_facets <- length(unique(df()[[input$split_by]]))
          if (n_facets < 1) n_facets <- 1
          
          ncol_facets <- suppressWarnings(as.numeric(input$max_facets))
          if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
          
          nrow_facets <- ceiling(n_facets / ncol_facets)
          if (!is.finite(nrow_facets) || nrow_facets < 1) nrow_facets <- 1
          
          pdf_width  <- 6 * ncol_facets
          pdf_height <- 6 * nrow_facets
        } else {
          # Single panel → fixed size
          pdf_width  <- 8
          pdf_height <- 8
        }
        
        ggsave(file, plot = gg, device = cairo_pdf,
               width = pdf_width, height = pdf_height, units = "in")
      },
      contentType = "application/pdf"
    )
    
    # Ensure the plot renders even when the tab is hidden (prevents some init races)
    outputOptions(output, "embed_plot", suspendWhenHidden = FALSE)
    
    observeEvent(input$clear_selection, current_sel(integer(0)))
  })
}

# --- Main UI ---
ui <- navbarPage(
  "FCView",
  id = "main_tab",
  header = tags$head(
    tags$style(HTML("
      /* Grey out disabled tabs */
      #main_tab li.disabled > a,
      #main_tab li.disabled > a:hover {
        color: #aaa !important;
        cursor: not-allowed;
        text-decoration: none;
      }
    ")),
    tags$script(HTML("
      (function() {
        var tabsLocked = true;
        function disableTabs() {
          tabsLocked = true;
          $('#main_tab li').addClass('disabled');
        }
        function enableTabs() {
          tabsLocked = false;
          $('#main_tab li').removeClass('disabled');
        }
        $(document).on('show.bs.tab click', '#main_tab a[data-toggle=\"tab\"]', function(e) {
          if (tabsLocked) {
            e.preventDefault();
            e.stopImmediatePropagation();
            return false;
          }
        });
        Shiny.addCustomMessageHandler('enableTabs', function(enable) {
          if (enable) {
            enableTabs();
          } else {
            disableTabs();
            var $home = $('#main_tab a[data-toggle=\"tab\"]').first();
            if ($home.length) $home.tab('show');
          }
        });
        $(document).ready(disableTabs);
        $(document).on('shiny:connected', disableTabs);
      })();
    "))
  ),
  tabPanel("Home",
           sidebarLayout(
             sidebarPanel(
               fileInput("rdata_upload", "Upload .RData (FCSimple analysis object)", accept = ".RData"),
               numericInput("max_cells_upload", "Max cells to read in", value = 300000, min = 1000, step = 1000),
               helpText("If the uploaded dataset has more cells than this number, it will be randomly downsampled after upload. This is done to speed up UMAP and tSNE facet plotting."), 
               width = 3
             ),
             mainPanel(
               h3("Dataset overview"),
               verbatimTextOutput("ds_summary"),
               h3("Available metadata"),
               tableOutput("meta_overview")
             )
           )
  ),
  
  tabPanel("UMAP", EmbeddingUI("umap", title = "UMAP")),
  tabPanel("tSNE", EmbeddingUI("tsne", title = "tSNE")),
  
  tabPanel(
    "Heatmap",
    fluidRow(
      column(3,
             checkboxInput("cluster_rows", "Cluster rows", value = TRUE),
             checkboxInput("cluster_columns", "Cluster columns", value = TRUE),
             selectInput("heatmap_theme", "Heatmap color theme",
                         choices = c("viridis", "heat", "greyscale"),
                         selected = "viridis"),
             br(),
             downloadButton("export_heatmap_pdf", "Export heatmap as PDF")
      ),
      column(9, plotOutput("cluster_heatmap", height = "700px"))
    )
  ),
  tabPanel("Testing",
           h4("Test Settings"),
           fluidRow(
             column(
               3,
               conditionalPanel(
                 condition = "output.hasClusterMap",
                 pickerInput("test_entity", "Entity", choices = c("Clusters", "Celltypes"), selected = "Clusters")
               ),
               pickerInput("group_var", "Categorical metadata", choices = NULL, options = list(`none-selected-text` = "None")),
               pickerInput("cont_var", "Continuous metadata", choices = NULL, options = list(`none-selected-text` = "None")),
               radioButtons("test_type", "Test",
                            choices = c("Wilcoxon (2-group)", "Kruskal–Wallis (multi-group)", "Spearman (continuous)")),
               selectInput("p_adj_method", "P‑value adjustment method",
                           choices = c("BH", "bonferroni", "BY", "fdr"), selected = "BH"),
               actionButton("run_test", "Run tests"),
               br(), br(),
               conditionalPanel(
                 condition = "output.hasResults",
                 downloadButton("export_results", "Export results as CSV"),
                 br(), br(),
                 actionButton("reset_test", "Clear Results")
               )
             ),
             column(
               9,
               conditionalPanel(
                 condition = "output.hasResults",
                 h4("Test Results"),
                 tableOutput("test_table")
               ), 
               textOutput("test_cleared_msg")
             )
           )
  ), 
  tabPanel("Categorical",
           h4("Plot Settings"),
           fluidRow(
             column(
               3,
               pickerInput("cat_entity", "Entity", choices = c("Clusters", "Celltypes"), selected = "Clusters"),
               pickerInput("cat_group_var", "Categorical metadata", choices = NULL, options = list(`none-selected-text` = "None")),
               radioButtons("cat_test_type", "Test", choices = c("Wilcoxon (2-group)", "Kruskal–Wallis (multi-group)")),
               selectInput("cat_p_adj_method", "P‑value adjustment method",
                           choices = c("BH", "bonferroni", "BY", "fdr"), selected = "BH"),
               checkboxInput("cat_use_adj_p", "Plot adjusted pvalues", value = TRUE),
               selectInput("cat_max_facets", "Facet columns", choices = 2:6, selected = 4),
               selectInput("cat_plot_type", "Plot type", choices = c("Boxplot" = "box", "Violin" = "violin"), selected = "box"),
               radioButtons("cat_points", "Show data points",
                            choices = c("Draw" = "draw", "Draw with jitter" = "jitter", "Do not draw" = "none"),
                            selected = "draw"),
               br(), 
               actionButton("cat_populate_colors", "Populate colors for selected group variable"),
               br(), 
               uiOutput("cat_color_pickers_ui"),
               br(), 
               actionButton("generate_cat_plots", "Generate plots"),
               br(), br(),
               conditionalPanel(
                 condition = "output.hasCatResults",
                 downloadButton("export_cat_pdf", "Export boxplots as PDF"),
                 br(), br(),
                 actionButton("reset_cat", "Clear Results")
               )
             ),
             column(
               9,
               conditionalPanel(
                 condition = "output.hasCatResults",
                 h4("Categorical Plots"),
                 plotOutput("categorical_plot")
               ), 
               textOutput("cat_cleared_msg")
             )
           )
  ), 
  tabPanel("Continuous",
           h4("Plot Settings"),
           fluidRow(
             column(
               3,
               pickerInput("cont_entity", "Entity", choices = c("Clusters", "Celltypes"), selected = "Clusters"),
               pickerInput("cont_group_var", "Continuous metadata", choices = NULL, options = list(`none-selected-text` = "None")),
               selectInput("cont_p_adj_method", "P‑value adjustment method",
                           choices = c("BH", "bonferroni", "BY", "fdr"), selected = "BH"),
               checkboxInput("cont_use_adj_p", "Plot adjusted pvalues", value = TRUE),
               checkboxInput("cont_transpose", "Transpose axes", value = FALSE),
               selectInput("cont_max_facets", "Facet columns", choices = 2:6, selected = 4),
               actionButton("generate_cont_plots", "Generate plots"),
               br(), br(),
               conditionalPanel(
                 condition = "output.hasContResults",
                 downloadButton("export_cont_pdf", "Export scatter plots as PDF"),
                 br(), br(),
                 actionButton("reset_cont", "Clear Results")
               )
             ),
             column(
               9,
               conditionalPanel(
                 condition = "output.hasContResults",
                 h4("Continuous Plots"),
                 plotOutput("continuous_plot")
               ), 
               textOutput("cont_cleared_msg")
             )
           )
  ), 
  tabPanel("Feature Selection",
           h4("Model Settings"),
           fluidRow(
             column(
               3,
               pickerInput("fs_method", "Method", 
                           choices = c("Ridge Regression", "Elastic Net", "Random Forest (Boruta)")),
               pickerInput("fs_outcome", "Outcome variable", choices = NULL),
               pickerInput("fs_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
               conditionalPanel(
                 condition = "input.fs_predictors.includes('leiden_cluster')",
                 pickerInput("fs_leiden_subset", "Select clusters to include",
                             choices = NULL, multiple = TRUE)
               ), 
               conditionalPanel(
                 condition = "input.fs_method == 'Elastic Net'",
                 sliderInput("fs_alpha", "Alpha (0 = Ridge, 1 = Lasso)", 
                             min = 0, max = 1, value = 0.5, step = 0.05)
               ),
               actionButton("run_fs", "Run Feature Selection"),
               br(), br(),
               conditionalPanel(
                 condition = "output.hasFSResults",
                 downloadButton("export_fs_results", "Export results as CSV"),
                 br(), br(),
                 actionButton("reset_fs", "Clear Results")
               )
             ),
             column(
               9,
               conditionalPanel(
                 condition = "output.hasFSResults",
                 h4("Summary Plot"),
                 plotOutput("fs_plot", height = "550px"),
                 h4("Selected Features"),
                 tableOutput("fs_results"), 
                 h4("Details"),
                 verbatimTextOutput("fs_summary")
               ),
               textOutput("fs_cleared_msg")
             )
           )
  ), 
  tabPanel("Classification",
           h4("Model Settings"),
           fluidRow(
             column(
               3,
               pickerInput("lm_outcome", "Outcome variable", choices = NULL),
               pickerInput("lm_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
               conditionalPanel(
                 condition = "input.lm_predictors.includes('leiden_cluster')",
                 pickerInput("lm_leiden_subset", "Select clusters to include", choices = NULL, multiple = TRUE)
               ),
               radioButtons("lm_model_type", "Model type",
                            choices = c("Logistic Regression", "Elastic Net", "Random Forest")),
               conditionalPanel(
                 condition = "input.lm_model_type == 'Elastic Net'",
                 sliderInput("lm_alpha", "Elastic Net alpha (0 = Ridge, 1 = Lasso)",
                             min = 0, max = 1, value = 0.5, step = 0.05)
               ),
               radioButtons("lm_validation", "Validation strategy",
                            choices = c("Train/Test split", "k-fold CV")),
               conditionalPanel(
                 condition = "input.lm_validation == 'Train/Test split'",
                 sliderInput("lm_train_frac", "Train fraction", min = 0.5, max = 0.95,
                             value = 0.7, step = 0.05)
               ),
               conditionalPanel(
                 condition = "input.lm_validation == 'k-fold CV'",
                 numericInput("lm_k", "Number of folds", value = 10, min = 2, max = 30)
               ),
               actionButton("run_lm", "Run Model"),
               br(), br(),
               conditionalPanel(
                 condition = "output.hasLMResults",
                 downloadButton("export_lm_zip", "Download All Results (ZIP)"),
                 br(), br(),
                 actionButton("reset_lm", "Clear Results")
               )
             ),
             column(
               9,
               conditionalPanel(
                 condition = "output.hasLMResults",
                 fluidRow(
                   column(
                     width = 5,
                     h4("Model Summary"),
                     verbatimTextOutput("lm_summary"),
                     h4("Performance Metrics"),
                     tableOutput("lm_perf_table")
                   ),
                   column(
                     width = 7,
                     h4("ROC Curve"),
                     plotOutput("lm_roc_plot", height = "500px"),
                     h4("Model Features"),
                     tableOutput("lm_features")
                   )
                 )
               ),
               textOutput("lm_cleared_msg")
             )
           )
  ), 
  tabPanel("Regression",
           h4("Model Settings"),
           fluidRow(
             column(
               3,
               pickerInput("reg_outcome", "Continuous outcome variable", choices = NULL),
               checkboxInput("reg_allow_categorical", 
                             "Allow categorical encoding", 
                             value = TRUE),
               pickerInput("reg_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
               conditionalPanel(
                 condition = "input.reg_predictors.includes('leiden_cluster')",
                 pickerInput("reg_leiden_subset", "Select clusters to include", choices = NULL, multiple = TRUE)
               ),
               radioButtons("reg_model_type", "Model type",
                            choices = c("Linear Regression" = "lm",
                                        "Elastic Net" = "glmnet",
                                        "Random Forest Regression" = "rf",
                                        "Extreme Gradient Boosting (XGBoost)" = "xgbTree")),
               conditionalPanel(
                 condition = "input.reg_model_type == 'glmnet'",
                 sliderInput("reg_alpha", "Elastic Net alpha (0 = Ridge, 1 = Lasso)",
                             min = 0, max = 1, value = 0.5, step = 0.05)
               ),
               radioButtons("reg_validation", "Validation strategy",
                            choices = c("Train/Test split" = "split",
                                        "k-fold CV" = "cv")),
               conditionalPanel(
                 condition = "input.reg_validation == 'split'",
                 sliderInput("reg_train_frac", "Train fraction", min = 0.5, max = 0.95,
                             value = 0.7, step = 0.05)
               ),
               conditionalPanel(
                 condition = "input.reg_validation == 'cv'",
                 numericInput("reg_k", "Number of folds", value = 5, min = 2, max = 20)
               ),
               actionButton("run_reg", "Run Regression Model"),
               br(), br(),
               conditionalPanel(
                 condition = "output.hasRegResults",
                 downloadButton("export_reg_zip", "Download All Results (ZIP)"),
                 br(), br(),
                 actionButton("reset_reg", "Clear Results")
               )
             ),
             column(
               9,
               conditionalPanel(
                 condition = "output.hasRegResults",
                 fluidRow(
                   column(
                     width = 6,
                     # h4("Model Summary"),
                     # verbatimTextOutput("reg_summary"),
                     h4("Performance Metrics"),
                     tableOutput("reg_perf_table"), 
                     h4("Residual Diagnostics"),
                     plotOutput("reg_resid_plot", height = "700px")
                   ),
                   column(
                     width = 6,
                     h4("Predicted vs Observed"),
                     plotOutput("reg_pred_plot", height = "375px"),
                     h4("Model Features"),
                     tableOutput("reg_features")
                     # Removed Tree-based Interpretability (PDP/SHAP) since not implemented
                   )
                 )
               ),
               textOutput("reg_cleared_msg")
             )
           )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  # App-wide stores
  rv <- reactiveValues(
    obj = NULL,
    expr = NULL,
    meta_cell = NULL,
    clusters = NULL,
    cluster_map = NULL,
    UMAP = NULL,
    tSNE = NULL,
    cluster_heat = NULL,
    pop_size = NULL,
    rep_used = NULL,
    data_ready = FALSE
  )
  cat_plot_cache <- reactiveVal(NULL)
  cont_plot_cache <- reactiveVal(NULL)
  test_results_rv <- reactiveVal(NULL)
  cat_state <- reactiveVal(NULL)
  cont_state <- reactiveVal(NULL)
  fs_state <- reactiveVal(NULL)
  lm_state <- reactiveVal(NULL)
  reg_state <- reactiveVal(NULL)
  # rv$log <- reactiveVal(character())
  
  # Disable tabs at startup
  # observe({
  #   if (!isTRUE(rv$data_ready)) {
  #     session$sendCustomMessage("enableTabs", FALSE)
  #   }
  # })
  
  observeEvent(input$main_tab, {
    message("Tab changed to: ", input$main_tab)
  })
  
  # Upload RData and initialize datasets immediately (no mapping button)
  observeEvent(input$rdata_upload, {
    session$sendCustomMessage("enableTabs", FALSE)  # lock tabs during load
    req(input$rdata_upload)
    
    # Load all objects from the uploaded .RData
    e <- new.env(parent = emptyenv())
    load(input$rdata_upload$datapath, envir = e)
    
    # Find all objects with the required structure
    candidates <- ls(e)
    matches <- vapply(candidates, function(nm) {
      cand <- e[[nm]]
      is.list(cand) && all(c("data", "source", "metadata", "leiden") %in% names(cand))
    }, logical(1))
    
    if (sum(matches) == 1) {
      obj <- e[[candidates[matches][1]]]
    } else if (sum(matches) > 1) {
      showNotification(
        paste("Multiple valid objects found in .RData:",
              paste(candidates[matches], collapse = ", ")),
        type = "error"
      )
      return()
    } else {
      showNotification("No valid input object found in .RData", type = "error")
      return()
    }
    
    # Ensure embeddings are data.frames
    if ("tsne" %in% names(obj) && is.matrix(obj$tsne$coordinates)) {
      obj$tsne$coordinates <- as.data.frame(obj$tsne$coordinates)
    }
    if ("umap" %in% names(obj) && is.matrix(obj$umap$coordinates)) {
      obj$umap$coordinates <- as.data.frame(obj$umap$coordinates)
    }
    
    # --- Downsample cells if needed (NEVER downsample sample-level metadata) ---
    n_cells <- nrow(obj$data)
    max_cells <- input$max_cells_upload %||% 300000
    
    if (n_cells > max_cells) {
      set.seed(123)
      keep_idx <- sort(sample.int(n_cells, max_cells))
      
      obj$data   <- obj$data[keep_idx, , drop = FALSE]
      obj$source <- obj$source[keep_idx]
      
      if (!is.null(obj$leiden$clusters) && length(obj$leiden$clusters) == n_cells) {
        obj$leiden$clusters <- obj$leiden$clusters[keep_idx]
      }
      if (!is.null(obj$umap$coordinates) && nrow(obj$umap$coordinates) == n_cells) {
        obj$umap$coordinates <- obj$umap$coordinates[keep_idx, , drop = FALSE]
      }
      if (!is.null(obj$tsne$coordinates) && nrow(obj$tsne$coordinates) == n_cells) {
        obj$tsne$coordinates <- obj$tsne$coordinates[keep_idx, , drop = FALSE]
      }
      if (!is.null(obj$run_date) && length(obj$run_date) == n_cells) {
        obj$run_date <- obj$run_date[keep_idx]
      }
      
      message(sprintf("Downsampled from %d to %d cells at upload stage", n_cells, max_cells))
      showNotification(sprintf("Downsampled from %d to %d cells", n_cells, max_cells), type = "message")
    }
    
    # Validate the (possibly downsampled) object
    validateInput(obj, id_col = NULL)
    
    expr <- obj$data
    run_date <- obj$run_date %||% NULL
    
    # --- Build meta_cell (per-cell) with robust patient_ID mapping for UMAP/tSNE tabs ---
    meta_cell <- data.frame(
      source  = as.character(obj$source),
      RunDate = if (!is.null(run_date)) run_date else NA
    )
    
    # Prepare metadata IDs (character, unique, deduplicated rows)
    if (!("patient_ID" %in% colnames(obj$metadata))) {
      showNotification("metadata does not contain 'patient_ID' column.", type = "error")
      return()
    }
    obj$metadata$patient_ID <- as.character(obj$metadata$patient_ID)
    
    metadata_unique <- obj$metadata %>%
      dplyr::distinct(patient_ID, .keep_all = TRUE)
    
    ids <- unique(metadata_unique$patient_ID)
    if (length(ids) == 0) {
      showNotification("No patient_ID values found in metadata.", type = "error")
      return()
    }
    
    # Escape regex metacharacters in IDs, then sort by length (longest first)
    ids_escaped <- stringr::str_replace_all(ids, "([\\^$.|?*+()\\[\\]{}\\\\])", "\\\\\\1")
    ids_escaped <- ids_escaped[order(nchar(ids_escaped), decreasing = TRUE)]
    pattern <- paste0("(", paste0(ids_escaped, collapse = "|"), ")")
    
    # Extract patient_ID from source via regex built from metadata IDs
    meta_cell$patient_ID <- stringr::str_extract(meta_cell$source, pattern)
    
    # Post-extraction sanity check (IDs extracted but not present in metadata)
    not_in_meta <- setdiff(unique(na.omit(meta_cell$patient_ID)), metadata_unique$patient_ID)
    if (length(not_in_meta) > 0) {
      showNotification(
        paste0("Extracted patient_IDs not found in metadata: ", paste(not_in_meta, collapse = ", ")),
        type = "warning",
        duration = NULL
      )
    }
    
    # Join metadata by patient_ID (safe: metadata_unique has 1 row per ID)
    meta_cell <- meta_cell %>%
      dplyr::left_join(metadata_unique, by = "patient_ID")
    meta_cell[] <- lapply(meta_cell, function(col) {
      if (is.character(col) && all(grepl("^\\s*-?\\d*(\\.\\d+)?\\s*$", col[!is.na(col)]))) {
        as.numeric(col)
      } else {
        col
      }
    })
    
    unmatched <- sum(is.na(meta_cell$patient_ID))
    if (unmatched > 0) {
      showNotification(
        paste0("Warning: ", unmatched, " cells could not be matched to metadata (patient_ID missing)."),
        type = "warning",
        duration = NULL
      )
      message("First few unmatched source entries:")
      print(utils::head(meta_cell$source[is.na(meta_cell$patient_ID)], 10))
    }
    
    if (!("PatientID" %in% names(meta_cell))) {
      meta_cell$PatientID <- meta_cell$patient_ID
    }
    
    if ("RunDate" %in% names(meta_cell)) {
      meta_cell$RunDate <- as.factor(meta_cell$RunDate)
    }
    
    # Add leiden_cluster factor if available
    if (!"leiden_cluster" %in% names(meta_cell) &&
        !is.null(obj$leiden$clusters) &&
        length(obj$leiden$clusters) == nrow(obj$data)) {
      meta_cell$leiden_cluster <- factor(obj$leiden$clusters,
                                         levels = sort(unique(obj$leiden$clusters)))
    }
    
    clusters <- list(
      assignments = obj$leiden$clusters,
      settings    = obj$leiden$settings %||% list()
    )
    cluster_map <- if (!is.null(obj$cluster_mapping)) obj$cluster_mapping else NULL
    
    UMAP <- if (!is.null(obj$umap)) list(coords = obj$umap$coordinates,
                                         settings = obj$umap$settings) else NULL
    tSNE <- if (!is.null(obj$tsne)) list(coords = obj$tsne$coordinates,
                                         settings = obj$tsne$settings) else NULL
    
    cluster_heat <- if (!is.null(obj$leiden_heatmap)) obj$leiden_heatmap$heatmap_tile_data else NULL
    pop_size     <- if (!is.null(obj$leiden_heatmap)) obj$leiden_heatmap$population_size else NULL
    rep_used     <- if (!is.null(obj$leiden_heatmap)) obj$leiden_heatmap$rep_used else NA
    
    # --- Add per-sample abundance matrix and canonical per-sample metadata for downstream tabs ---
    if (!is.null(obj$leiden$abundance) && is.matrix(obj$leiden$abundance)) {
      clusters$abundance <- obj$leiden$abundance
      rv$abundance_sample <- clusters$abundance
      
      if (is.null(rownames(rv$abundance_sample)) || any(!nzchar(rownames(rv$abundance_sample)))) {
        showNotification(
          "leiden$abundance has missing rownames; cannot map to metadata. Please set rownames to source strings.",
          type = "error",
          duration = NULL
        )
        message("Abundance matrix rownames are missing or empty; mapping to metadata will fail.")
      } else {
        message(sprintf("Abundance matrix loaded: %d sources × %d entities",
                        nrow(rv$abundance_sample), ncol(rv$abundance_sample)))
      }
    } else {
      clusters$abundance <- NULL
      rv$abundance_sample <- NULL
      showNotification("No leiden$abundance matrix found in upload; abundance-based tabs will be disabled.", type = "warning")
    }
    
    # Canonical per-sample metadata (prefer direct sample-level object if provided)
    # If obj already includes a per-sample metadata frame, use it; otherwise derive from metadata_unique.
    if (!is.null(obj$metadata_sample) && is.data.frame(obj$metadata_sample) &&
        "patient_ID" %in% colnames(obj$metadata_sample)) {
      rv$meta_sample <- obj$metadata_sample %>%
        dplyr::distinct(patient_ID, .keep_all = TRUE)
    } else {
      rv$meta_sample <- metadata_unique  # one row per patient_ID
    }
    
    # Store in rv (cell-level objects preserved for UMAP/tSNE/Heatmap)
    rv$expr         <- expr
    rv$meta_cell    <- meta_cell
    rv$clusters     <- clusters
    rv$cluster_map  <- cluster_map
    rv$UMAP         <- UMAP
    rv$tSNE         <- tSNE
    rv$cluster_heat <- cluster_heat
    rv$pop_size     <- pop_size
    rv$rep_used     <- rep_used
    
    message("Upload complete: expr rows=", nrow(rv$expr),
            " meta_cell rows=", nrow(rv$meta_cell),
            " meta_sample rows=", nrow(rv$meta_sample),
            " abundance_sample rows=", if (!is.null(rv$abundance_sample)) nrow(rv$abundance_sample) else 0,
            " UMAP coords=", if (!is.null(rv$UMAP)) nrow(rv$UMAP$coords) else 0,
            " tSNE coords=", if (!is.null(rv$tSNE)) nrow(rv$tSNE$coords) else 0)
    
    showNotification("Data loaded and initialized.", type = "message")
    
    # Initialize FS/LM pickers from per-sample objects
    if (!is.null(rv$meta_sample)) {
      meta_cols <- colnames(rv$meta_sample)
      categorical_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.factor(x) || is.character(x))])
      predictor_choices <- sort(meta_cols)
      
      # Feature Selection
      updatePickerInput(session, "fs_outcome", choices = categorical_choices, selected = NULL)
      updatePickerInput(session, "fs_predictors", choices = unique(c(predictor_choices, "leiden_cluster")), selected = NULL)
      updatePickerInput(session, "fs_leiden_subset",
                        choices = if (!is.null(rv$abundance_sample)) colnames(rv$abundance_sample) else character(0),
                        selected = character(0))
      
      # Classification
      updatePickerInput(session, "lm_outcome", choices = categorical_choices, selected = NULL)
      updatePickerInput(session, "lm_predictors", choices = unique(c(predictor_choices, "leiden_cluster")), selected = NULL)
      updatePickerInput(session, "lm_leiden_subset",
                        choices = if (!is.null(rv$abundance_sample)) colnames(rv$abundance_sample) else character(0),
                        selected = character(0))
      
      # Regression
      continuous_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.numeric(x) || is.integer(x))])
      updatePickerInput(session, "reg_outcome",
                        choices = continuous_choices,
                        selected = if (length(continuous_choices) > 0) continuous_choices[1])
      updatePickerInput(session, "reg_predictors", choices = unique(c(predictor_choices, "leiden_cluster")), selected = NULL)
      updatePickerInput(session, "reg_leiden_subset",
                        choices = if (!is.null(rv$abundance_sample)) colnames(rv$abundance_sample) else character(0),
                        selected = character(0))
    }
    
    rv$data_ready <- TRUE
    session$sendCustomMessage("enableTabs", TRUE)
  })
  
  observeEvent(list(rv$meta_sample, input$reg_allow_categorical), {
    req(rv$meta_sample)
    meta_cols <- colnames(rv$meta_sample)
    
    numeric_cols <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.numeric(x) || is.integer(x))])
    categorical_cols <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.factor(x) || is.character(x))])
    
    if (isTRUE(input$reg_allow_categorical)) {
      predictor_choices <- sort(c(numeric_cols, categorical_cols))
    } else {
      predictor_choices <- numeric_cols
    }
    
    updatePickerInput(session, "reg_predictors",
                      choices = unique(c(predictor_choices, "leiden_cluster")),
                      selected = NULL
    )
  })
  
  observeEvent(input$main_tab, {
    if (!isTRUE(rv$data_ready) && !identical(input$main_tab, "Home")) {
      updateNavbarPage(session, "main_tab", selected = "Home")
    }
  })
  
  # UI-facing flag for conditionalPanel (no nested reactive)
  output$hasClusterMap <- reactive({
    !is.null(rv$cluster_map) && all(c("cluster", "celltype") %in% names(rv$cluster_map))
  })
  outputOptions(output, "hasClusterMap", suspendWhenHidden = FALSE)
  
  observeEvent(rv$cluster_map, {
    cm <- rv$cluster_map
    if (!is.null(cm) && all(c("cluster","celltype") %in% names(cm))) {
      # nothing to do
    } else {
      updatePickerInput(session, "test_entity", selected = "Clusters")
    }
  }, ignoreInit = TRUE)
  
  observeEvent(rv$cluster_map, {
    cm <- rv$cluster_map
    if (!is.null(cm) && all(c("cluster","celltype") %in% names(cm))) {
      # nothing to do
    } else {
      updatePickerInput(session, "cat_entity", selected = "Clusters")
    }
  }, ignoreInit = TRUE)
  
  # Auto-detect categorical vs continuous metadata
  observeEvent(rv$meta_sample, {
    meta_cols <- colnames(rv$meta_sample)
    categorical_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.character(x) || is.factor(x))])
    continuous_choices  <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.numeric(x) || is.integer(x))])
    
    updatePickerInput(session, "group_var", choices = c("", categorical_choices), selected = "")
    updatePickerInput(session, "cont_var",  choices = c("", continuous_choices),  selected = "")
  }, ignoreInit = TRUE)
  
  observeEvent(rv$meta_sample, {
    meta_cols <- colnames(rv$meta_sample)
    continuous_choices <- sort(meta_cols[sapply(rv$meta_sample, is.numeric)])
    updatePickerInput(session, "cont_group_var", choices = c("", continuous_choices), selected = "")
  }, ignoreInit = TRUE)
  
  observeEvent(rv$meta_sample, {
    meta_cols <- sort(colnames(rv$meta_sample))
    updatePickerInput(session, "model_outcome", choices = meta_cols)
    updatePickerInput(session, "model_predictors", choices = meta_cols)
    updatePickerInput(session, "model_covariates", choices = meta_cols)
    updatePickerInput(session, "model_random", choices = meta_cols)
  }, ignoreInit = TRUE)
  
  # Launch embedding modules as soon as data is ready
  observeEvent(list(rv$UMAP, rv$data_ready), {
    req(rv$data_ready, rv$UMAP, rv$expr, rv$meta_cell, rv$clusters)
    message("Launching UMAP module")
    EmbeddingServer(
      "umap", "UMAP",
      reactive(rv$UMAP$coords),
      reactive(rv$expr),
      reactive(rv$meta_cell),
      reactive(rv$clusters),
      reactive(rv$cluster_map),
      reactive(input$main_tab)
    )
  }, ignoreInit = TRUE)
  
  observeEvent(list(rv$tSNE, rv$data_ready), {
    req(rv$data_ready, rv$tSNE, rv$expr, rv$meta_cell, rv$clusters)
    message("Launching tSNE module")
    EmbeddingServer(
      "tsne", "tSNE",
      reactive(rv$tSNE$coords),
      reactive(rv$expr),
      reactive(rv$meta_cell),
      reactive(rv$clusters),
      reactive(rv$cluster_map),
      reactive(input$main_tab)
    )
  }, ignoreInit = TRUE)
  
  # Home summaries
  output$ds_summary <- renderPrint({
    req(rv$expr, rv$clusters)
    list(
      n_cells = nrow(rv$expr),
      n_markers = ncol(rv$expr),
      n_clusters = length(unique(rv$clusters$assignments)),
      embeddings_available = c(UMAP = !is.null(rv$UMAP), tSNE = !is.null(rv$tSNE)),
      has_cluster_heatmap = !is.null(rv$cluster_heat),
      rep_used = rv$rep_used
    )
  })
  
  # Metadata overview
  output$meta_overview <- renderTable({
    req(rv$meta_sample)
    data.frame(
      name = colnames(rv$meta_sample),
      type = sapply(rv$meta_sample, function(x) class(x)[1]),
      example = sapply(rv$meta_sample, function(x) paste(utils::head(unique(x), 3), collapse = ", "))
    )
  }, sanitize.text.function = function(x) x)
  
  heatmap_obj <- reactive({
    req(rv$cluster_heat)
    M <- rv$cluster_heat
    
    # Clean row/col names
    if (!is.null(rownames(M))) rownames(M) <- gsub("\\n", " ", rownames(M))
    if (!is.null(colnames(M))) colnames(M) <- gsub("\\n", " ", colnames(M))
    
    # Annotation
    ranno <- NULL
    if (!is.null(rv$pop_size)) {
      size_vals <- rv$pop_size[, 1]
      size_col_fun <- circlize::colorRamp2(
        c(min(size_vals, na.rm = TRUE), max(size_vals, na.rm = TRUE)),
        c("white", "red")
      )
      ranno <- rowAnnotation(
        Size = size_vals,
        col = list(Size = size_col_fun),
        gp = grid::gpar(col = "black", lwd = 0.5)
      )
    }
    
    # Palette
    palette_choice <- switch(
      input$heatmap_theme,
      "viridis" = {
        col_seq <- seq(min(M), max(M), length.out = n <- 256)
        circlize::colorRamp2(col_seq, viridis(n))
      },
      "heat" = {
        col_seq <- seq(min(M), max(M), length.out = n <- 256)
        circlize::colorRamp2(
          col_seq,
          colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(n)
        )
      },
      "greyscale" = {
        col_seq <- seq(min(M), max(M), length.out = n <- 256)
        circlize::colorRamp2(
          col_seq,
          colorRampPalette(c("black", "grey50", "white"))(n)
        )
      }
    )
    
    # Return the Heatmap object
    Heatmap(
      M,
      name = "expr",
      cluster_rows = isTRUE(input$cluster_rows),
      cluster_columns = isTRUE(input$cluster_columns),
      right_annotation = ranno,
      row_names_side = "left",
      rect_gp = gpar(lwd = 0.33, col = "black"),
      col = palette_choice
    )
  })
  
  output$cluster_heatmap <- renderPlot({
    draw(heatmap_obj())
  })
  
  output$export_heatmap_pdf <- downloadHandler(
    filename = function() {
      paste0("cluster_heatmap_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      # Get the underlying matrix from your reactive
      M <- rv$cluster_heat  # or heatmap_obj()$matrix if you store it there
      
      # Calculate dimensions using the same factors as fcs_plot_heatmap
      pdf_width  <- (ncol(M) * 0.33) + 3.25
      pdf_height <- (nrow(M) * 0.3)  + 2.25
      
      pdf(file, width = pdf_width, height = pdf_height)
      draw(heatmap_obj())
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  # Store results + adj_col name from the run
  run_tests <- eventReactive(input$run_test, {
    req(rv$meta_sample, rv$abundance_sample)
    assert_sample_level(rv$meta_sample, "Testing")
    
    test_type_run   <- input$test_type
    p_adj_method_run <- input$p_adj_method
    group_var_run   <- input$group_var
    cont_var_run    <- input$cont_var
    test_entity_run <- input$test_entity
    
    abund0 <- rv$abundance_sample
    if (is.null(abund0)) {
      showNotification("No abundance matrix available.", type = "error")
      return(list(df = NULL, adj_col = NULL))
    }
    
    # Aggregate to celltypes if requested
    abund <- abund0
    if (test_entity_run == "Celltypes" && !is.null(rv$cluster_map)) {
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) {
        showNotification("No overlapping clusters to aggregate into celltypes.", type = "error")
        return(list(df = NULL, adj_col = NULL))
      }
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }
    
    # Merge with per-sample metadata
    meta_sub <- rv$meta_sample %>% dplyr::select(patient_ID, dplyr::everything())
    abund_df <- as.data.frame(abund)
    abund_df$patient_ID <- stringr::str_extract(string = rownames(abund_df), pattern = paste0('(',paste0(meta_sub$patient_ID,collapse='|'),')'))
    merged <- merge(x = meta_sub, y = abund_df, by = 'patient_ID')
    # Long format
    abund_long <- merged %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq")
    
    # Run tests per entity
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        if (test_type_run == "Wilcoxon (2-group)") {
          if (!nzchar(group_var_run)) return(data.frame(test = "wilcox", p = NA, n = nrow(.x)))
          
          g_raw <- .x[[group_var_run]]
          ok <- !is.na(g_raw)
          g <- droplevels(factor(g_raw[ok]))
          freq_ok <- .x$freq[ok]
          
          if (length(levels(g)) != 2) {
            # Not exactly 2 groups after NA removal
            return(data.frame(test = "wilcox", p = NA, n = sum(ok)))
          }
          
          wt <- suppressWarnings(wilcox.test(freq_ok ~ g))
          data.frame(test = "wilcox", n = sum(ok), p = wt$p.value)
        } else if (test_type_run == "Kruskal–Wallis (multi-group)") {
          if (!nzchar(group_var_run)) {
            return(data.frame(test = "kruskal", p = NA, n = nrow(.x)))
          }
          
          g_raw <- .x[[group_var_run]]
          ok <- !is.na(g_raw)
          g <- droplevels(factor(g_raw[ok]))
          freq_ok <- .x$freq[ok]
          
          if (length(levels(g)) < 2) {
            # Not enough groups with data
            return(data.frame(test = "kruskal", p = NA, n = sum(ok)))
          }
          
          kw <- suppressWarnings(kruskal.test(freq_ok ~ g))
          data.frame(test = "kruskal", n = sum(ok), p = kw$p.value)
        } else if (test_type_run == "Spearman (continuous)") {
          if (!nzchar(cont_var_run)) {
            return(data.frame(test = "spearman", rho = NA, p = NA, n = nrow(.x)))
          }
          
          x <- .x$freq
          y <- suppressWarnings(as.numeric(.x[[cont_var_run]]))
          
          ok <- complete.cases(x, y)
          if (sum(ok) < 3) {
            return(data.frame(test = "spearman", rho = NA, p = NA, n = sum(ok)))
          }
          
          ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
          data.frame(test = "spearman", rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
        }
      }) %>%
      dplyr::ungroup()
    
    # Adjust p-values
    adj_col <- NULL
    if (nrow(res) && "p" %in% names(res) && nzchar(p_adj_method_run)) {
      adj_col <- paste0(tolower(p_adj_method_run), "_padj")
      res[[adj_col]] <- p.adjust(res$p, method = p_adj_method_run)
    }
    
    list(df = res, adj_col = adj_col)
  })
  
  observeEvent(input$run_test, {
    res <- run_tests()
    test_results_rv(res)
    output$test_cleared_msg <- renderText(NULL)
  })
  
  observeEvent(input$reset_test, {
    test_results_rv(NULL)        # clear
    showNotification("Testing results cleared.", type = "message", duration = 5)
    output$test_cleared_msg <- renderText("Results cleared. Run a new test to see results here.")
  })
  
  output$hasResults <- reactive({
    run <- test_results_rv()
    !is.null(run) && !is.null(run$df) && nrow(run$df) > 0
  })
  outputOptions(output, "hasResults", suspendWhenHidden = FALSE)
  
  output$test_table <- renderTable({
    run <- test_results_rv()
    df <- req(run$df)
    adj_col <- run$adj_col
    
    # sort while numeric
    if (!is.null(adj_col) && adj_col %in% names(df) && is.numeric(df[[adj_col]])) {
      df <- df[order(df[[adj_col]], na.last = TRUE), ]
    }
    
    # format for display
    df_display <- df
    num_cols <- intersect(c("p", adj_col), names(df_display))
    for (col in num_cols) {
      if (is.numeric(df_display[[col]])) {
        df_display[[col]] <- formatC(df_display[[col]], format = "f", digits = 3)
      }
    }
    if ("rho" %in% names(df_display) && is.numeric(df_display$rho)) {
      df_display$rho <- formatC(df_display$rho, format = "f", digits = 2)
    }
    
    df_display
  }, sanitize.text.function = function(x) x)
  
  output$export_results <- downloadHandler(
    filename = function() {
      info <- rv$last_test_info
      if (is.null(info)) return("results.csv")
      fname <- paste(info$entity, info$test, info$metadata, sep = "_")
      fname <- gsub(" ", "_", fname)
      fname <- tolower(fname)
      paste0(fname, ".csv")
    },
    content = function(file) {
      run <- run_tests()
      df <- run$df
      req(!is.null(df), nrow(df) > 0)
      write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  get_categorical_meta <- function(meta_df) {
    cols <- colnames(meta_df)
    sort(cols[sapply(meta_df, function(x) is.character(x) || is.factor(x))])
  }
  
  get_continuous_meta <- function(meta_df) {
    cols <- colnames(meta_df)
    sort(cols[sapply(meta_df, function(x) is.numeric(x) || is.integer(x))])
  }
  
  observeEvent(rv$meta_sample, {
    updatePickerInput(session, "group_var", choices = c("", get_categorical_meta(rv$meta_sample)), selected = "")
    updatePickerInput(session, "cont_var",  choices = c("", get_continuous_meta(rv$meta_sample)),  selected = "")
  })
  
  observeEvent(rv$meta_sample, {
    updatePickerInput(session, "cat_group_var", choices = c("", get_categorical_meta(rv$meta_sample)), selected = "")
  })
  
  observeEvent(rv$meta_sample, {
    updatePickerInput(session, "cont_group_var", choices = c("", get_continuous_meta(rv$meta_sample)), selected = "")
  })
  
  # Helper to make safe input IDs from arbitrary group values
  sanitize_id <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
    x
  }
  
  # Observe any manual color changes and persist into rv$cat_colors
  observe({
    # If rv$cat_colors exists, watch its keys and update from inputs
    req(!is.null(rv$cat_colors))
    new_colors <- rv$cat_colors
    for (grp in names(rv$cat_colors)) {
      input_name <- paste0("cat_color_", sanitize_id(grp))
      if (!is.null(input[[input_name]])) {
        new_colors[[grp]] <- input[[input_name]]
      }
    }
    rv$cat_colors <- new_colors
  })
  
  cat_plot_data <- eventReactive(input$generate_cat_plots, {
    req(rv$meta_sample, rv$abundance_sample)
    assert_sample_level(rv$meta_sample, "Categorical")
    
    abund0 <- rv$abundance_sample
    if (is.null(abund0)) {
      showNotification("No abundance matrix available.", type = "error")
      return(NULL)
    }
    
    # Aggregate to celltypes if needed
    abund <- abund0
    if (input$cat_entity == "Celltypes" && !is.null(rv$cluster_map)) {
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) return(NULL)
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }
    
    # Merge with per-sample metadata
    meta_sub <- rv$meta_sample %>% dplyr::select(patient_ID, dplyr::everything())
    abund_df <- as.data.frame(abund, check.names = FALSE, stringsAsFactors = FALSE)
    abund_df$patient_ID <- stringr::str_extract(string = rownames(abund_df), pattern = paste0('(',paste0(meta_sub$patient_ID,collapse='|'),')'))
    # merged <- dplyr::left_join(meta_sub, abund_df, by = "patient_ID")
    merged <- merge(x = meta_sub, y = abund_df, by = 'patient_ID')
    
    # Long format and clean
    abund_long <- merged %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq") %>%
      dplyr::mutate(entity = gsub(pattern = "\\n", replacement = " ", x = entity)) %>%
      dplyr::filter(!is.na(freq))
    
    # Capture plot-type and point-mode at Generate time
    plot_type_selected <- input$cat_plot_type %||% "box"
    point_mode_selected <- input$cat_points %||% "draw"
    
    # Run test per entity
    test_type <- input$cat_test_type
    group_var <- input$cat_group_var
    
    if (is.null(group_var) || !nzchar(group_var)) {
      res <- abund_long %>%
        dplyr::group_by(entity) %>%
        dplyr::summarise(p = NA_real_, .groups = "drop")
    } else {
      res <- abund_long %>%
        dplyr::group_by(entity) %>%
        dplyr::group_modify(~ {
          g_raw <- .x[[group_var]]
          ok_rows <- !is.na(g_raw)
          if (!any(ok_rows)) return(data.frame(p = NA_real_))
          g <- droplevels(factor(g_raw[ok_rows]))
          freq_ok <- .x$freq[ok_rows]
          if (length(unique(g)) < 2) return(data.frame(p = NA_real_))
          if (test_type == "Wilcoxon (2-group)" && length(unique(g)) == 2) {
            wt <- suppressWarnings(wilcox.test(freq_ok ~ g))
            data.frame(p = wt$p.value)
          } else if (test_type == "Kruskal–Wallis (multi-group)") {
            kw <- kruskal.test(freq_ok ~ as.factor(g))
            data.frame(p = kw$p.value)
          } else {
            data.frame(p = NA_real_)
          }
        }) %>%
        dplyr::ungroup()
    }
    
    # Adjust p-values if requested
    if (nrow(res) && nzchar(input$cat_p_adj_method) && "p" %in% names(res)) {
      res$padj <- p.adjust(res$p, method = input$cat_p_adj_method)
    }
    
    # Save info for export
    rv$last_cat_info <- list(
      entity = tolower(input$cat_entity %||% "clusters"),
      group = tolower(input$cat_group_var %||% "group"),
      test_raw = input$cat_test_type %||% "test"
    )
    
    # Return everything needed for plotting
    list(
      data = abund_long,
      results = res,
      group_var = input$cat_group_var,
      use_adj_p = input$cat_use_adj_p,
      facet_cols = as.numeric(input$cat_max_facets),
      plot_type = plot_type_selected,
      point_mode = point_mode_selected
    )
  })
  
  observeEvent(input$generate_cat_plots, {
    cp <- cat_plot_data()     # your existing eventReactive
    cat_state(cp)
    output$cat_cleared_msg <- renderText(NULL)
  })
  
  observeEvent(input$reset_cat, {
    cat_state(NULL)           # clear the state
    cat_plot_cache(NULL)      # also clear cached ggplot
    showNotification("Categorical plots cleared.", type = "message", duration = 5)
    output$cat_cleared_msg <- renderText("Results cleared. Generate new plots to see them here.")
  })
  
  output$hasCatResults <- reactive({
    cp <- cat_state()
    !is.null(cp) && !is.null(cp$data) && nrow(cp$data) > 0
  })
  outputOptions(output, "hasCatResults", suspendWhenHidden = FALSE)
  
  # Helper to sanitize IDs for input names
  sanitize_id <- function(x) {
    x <- as.character(x)
    gsub("[^A-Za-z0-9_\\-]", "_", x)
  }
  
  build_color_pickers_sample <- function(meta_df, group_var, plotted_df = NULL, input_prefix = "cat_color") {
    # Prefer group levels from plotted_df to match exactly what’s on the plot
    levels_vec <- NULL
    if (!is.null(plotted_df) && !is.null(group_var) && nzchar(group_var) && group_var %in% colnames(plotted_df)) {
      levels_vec <- unique(as.character(plotted_df[[group_var]]))
    }
    if (is.null(levels_vec) || length(levels_vec) == 0) {
      if (!is.null(group_var) && nzchar(group_var) && (group_var %in% colnames(meta_df))) {
        levels_vec <- sort(unique(as.character(meta_df[[group_var]])))
      } else {
        return(list(ui = NULL, colors = NULL, error = "No valid grouping variable selected for colors."))
      }
    }
    levels_vec <- levels_vec[!is.na(levels_vec)]
    if (length(levels_vec) == 0) {
      return(list(ui = NULL, colors = NULL, error = "No non-missing group levels found to populate colors for."))
    }
    default_pal <- viridis::viridis(length(levels_vec))
    names(default_pal) <- levels_vec
    
    ui_list <- lapply(seq_along(levels_vec), function(i) {
      lv <- levels_vec[i]
      input_id <- paste0(input_prefix, "_", sanitize_id(lv))
      colourpicker::colourInput(
        inputId = input_id,
        label = paste0("Color for ", lv),
        value = default_pal[i],
        showColour = "both"
      )
    })
    
    list(
      ui = tagList(tags$div(style = "max-height: 300px; overflow-y: auto; padding-right: 6px;", ui_list)),
      colors = setNames(as.character(default_pal), levels_vec),
      error = NULL
    )
  }
  
  # Observer using sample-level metadata and plotted groups
  observeEvent(input$cat_populate_colors, {
    req(rv$meta_sample)
    cp <- NULL
    # Try to use last plotted data (cat_state) for exact levels
    try({ cp <- cat_state() }, silent = TRUE)
    plotted_df <- if (!is.null(cp) && is.list(cp) && !is.null(cp$data)) cp$data else NULL
    group_var <- if (!is.null(cp)) cp$group_var else input$cat_group_var
    
    res <- build_color_pickers_sample(rv$meta_sample, group_var, plotted_df, "cat_color")
    if (!is.null(res$error)) {
      showNotification(res$error, type = "error")
      output$cat_color_pickers_ui <- renderUI(NULL)
      rv$cat_colors <- NULL
      return()
    }
    output$cat_color_pickers_ui <- renderUI(res$ui)
    rv$cat_colors <- res$colors
  })
  
  # Observe manual color changes
  observe({
    req(!is.null(rv$cat_colors))
    new_colors <- rv$cat_colors
    for (grp in names(rv$cat_colors)) {
      input_name <- paste0("cat_color_", sanitize_id(grp))
      if (!is.null(input[[input_name]])) {
        new_colors[[grp]] <- input[[input_name]]
      }
    }
    rv$cat_colors <- new_colors
  })
  
  # Reset colors when group_var changes
  observeEvent(input$cat_group_var, {
    output$cat_color_pickers_ui <- renderUI(NULL)
    rv$cat_colors <- NULL
  })
  
  
  output$categorical_plot <- renderPlot({
    # cp <- cat_plot_data(); req(cp)
    cp <- cat_state(); req(cp)
    abund_long <- cp$data
    res <- cp$results
    group_var <- cp$group_var
    use_adj_p <- cp$use_adj_p
    facet_cols <- cp$facet_cols
    plot_type <- cp$plot_type %||% input$cat_plot_type %||% "box"
    point_mode <- cp$point_mode %||% input$cat_points %||% "draw"
    
    # Validate grouping variable
    if (is.null(group_var) || !nzchar(group_var) || !(group_var %in% colnames(abund_long))) {
      showNotification("No valid grouping variable selected for categorical plotting.", type = "error")
      return(invisible(NULL))
    }
    
    # Remove rows with NA grouping values for plotting
    abund_long_plot <- abund_long %>% dplyr::filter(!is.na(.data[[group_var]]))
    
    if (nrow(abund_long_plot) == 0) {
      showNotification("No datapoints after removing NA frequencies or NA group values; nothing to plot.", type = "warning")
      return(invisible(NULL))
    }
    
    # Prepare p-value annotation dataframe
    p_df <- res %>%
      dplyr::mutate(
        p_to_show = if (isTRUE(use_adj_p) && "padj" %in% names(res)) padj else p,
        label = paste0("p = ", signif(p_to_show, 3)),
        x = length(unique(abund_long_plot[[group_var]])) / 2 + 0.5,
        y = tapply(abund_long_plot$freq, abund_long_plot$entity, max, na.rm = TRUE)[entity] * 1.15
      )
    
    # Determine group levels present in plotting data
    grp_levels <- unique(as.character(abund_long_plot[[group_var]]))
    if (length(grp_levels) == 0) {
      showNotification("No valid group levels found for plotting.", type = "error")
      return(invisible(NULL))
    }
    
    # Resolve manual colors non-reactively (isolate => changing pickers won't auto-trigger replot)
    colors_named <- isolate({
      if (!is.null(rv$cat_colors) && length(rv$cat_colors) > 0) {
        matched <- rv$cat_colors[names(rv$cat_colors) %in% grp_levels]
        missing_lvls <- setdiff(grp_levels, names(matched))
        if (length(missing_lvls) > 0) {
          filler <- viridis::viridis(length(missing_lvls))
          names(filler) <- missing_lvls
          matched <- c(matched, as.character(filler))
        }
        cols_ord <- as.character(matched[grp_levels])
        names(cols_ord) <- grp_levels
        cols_ord
      } else {
        pal <- viridis::viridis(length(grp_levels))
        setNames(as.character(pal), grp_levels)
      }
    })
    
    # Build ggplot
    gg <- ggplot2::ggplot(abund_long_plot, ggplot2::aes(x = .data[[group_var]], y = freq))
    
    if (identical(plot_type, "violin")) {
      gg <- gg + ggplot2::geom_violin(ggplot2::aes(fill = .data[[group_var]]), alpha = 0.7, scale = "width", trim = TRUE)
    } else {
      gg <- gg + ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[group_var]]), alpha = 0.7, outlier.shape = NA)
    }
    
    # Add points per captured point_mode (points use pch=21 with black border and fill mapped to group)
    if (identical(point_mode, "draw")) {
      gg <- gg + ggplot2::geom_point(
        ggplot2::aes(fill = .data[[group_var]]),
        position = ggplot2::position_jitter(width = 0.04, height = 0),
        pch = 21, size = 2.6, stroke = 0.3, color = "black", alpha = 0.9, inherit.aes = TRUE
      )
    } else if (identical(point_mode, "jitter")) {
      gg <- gg + ggplot2::geom_jitter(
        ggplot2::aes(fill = .data[[group_var]]),
        width = 0.15, height = 0, pch = 21, size = 2.6, stroke = 0.3, color = "black", alpha = 0.9
      )
    } else {
      # none -> no point layer
    }
    
    gg <- gg +
      ggplot2::facet_wrap(~entity, ncol = facet_cols, scales = "free_y") +
      ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::labs(x = group_var, y = "Frequency") +
      ggplot2::theme(
        strip.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 1.1, b = 1.1)),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1) # rotate facet x-axis labels 45 degrees
      )
    
    # Apply manual fill scale; points use black border so no color scale required for border
    gg <- gg + ggplot2::scale_fill_manual(values = colors_named)
    
    # Add p-value text annotations (filter to plotted entities)
    if (!is.null(p_df) && nrow(p_df) > 0) {
      p_df_plot <- p_df %>% dplyr::filter(entity %in% unique(abund_long_plot$entity))
      if (nrow(p_df_plot) > 0) {
        gg <- gg + ggplot2::geom_text(
          data = p_df_plot,
          ggplot2::aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          size = 5
        )
      }
    }
    
    # Cache for export
    cat_plot_cache(gg)
    
    gg
  },
  height = function() {
    gg <- cat_plot_cache()
    cp <- cat_plot_data()
    if (is.null(gg) || is.null(cp)) return(400)
    n_facets <- length(unique(gg$data$entity))
    ncol_facets <- cp$facet_cols
    if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
    nrow_facets <- ceiling(n_facets / ncol_facets)
    200 * nrow_facets
  },
  width = function() {
    gg <- cat_plot_cache()
    cp <- cat_plot_data()
    if (is.null(gg) || is.null(cp)) return(400)
    ncol_facets <- cp$facet_cols
    if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
    225 * ncol_facets
  })
  
  observeEvent(input$cat_group_var, {
    output$cat_color_pickers_ui <- renderUI(NULL)
    rv$cat_colors <- NULL
  })
  
  output$export_cat_pdf <- downloadHandler(
    filename = function() {
      info <- rv$last_cat_info
      if (is.null(info)) return("categorical_plots.pdf")
      test <- tolower(info$test_raw)
      test <- gsub("\\s+", "_", test)
      test <- gsub("_\\(.*\\)", "", test)
      test <- gsub("kruskal_wallis", "kruskal-wallis", test)
      paste0("categorical_", info$entity, "_", info$group, "_", test, ".pdf")
    },
    content = function(file) {
      gg <- cat_plot_cache()
      cp <- cat_plot_data()
      if (is.null(gg) || is.null(cp)) {
        showNotification("No plot available to export. Please generate the plots first.", type = "error")
        return()
      }
      n_facets <- length(unique(gg$data$entity))
      ncol_facets <- cp$facet_cols
      if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
      nrow_facets <- ceiling(n_facets / ncol_facets)
      if (!is.finite(nrow_facets) || nrow_facets < 1) nrow_facets <- 1
      pdf_width <- 4 * ncol_facets
      pdf_height <- 3 * nrow_facets
      ggsave(file, plot = gg, device = cairo_pdf, width = pdf_width, height = pdf_height, units = "in")
    },
    contentType = "application/pdf"
  )
  
  cont_plot_data <- eventReactive(input$generate_cont_plots, {
    req(rv$meta_sample, rv$abundance_sample)
    assert_sample_level(rv$meta_sample, "Continuous")
    
    abund0 <- rv$abundance_sample
    if (is.null(abund0)) {
      showNotification("No abundance matrix available.", type = "error")
      return(NULL)
    }
    
    # Aggregate to celltypes if requested
    abund <- abund0
    if (input$cont_entity == "Celltypes" && !is.null(rv$cluster_map)) {
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) return(NULL)
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }
    
    # Merge with per-sample metadata
    meta_sub <- rv$meta_sample %>% dplyr::select(patient_ID, dplyr::everything())
    abund_df <- as.data.frame(abund, check.names = FALSE, stringsAsFactors = FALSE)
    abund_df$patient_ID <- stringr::str_extract(string = rownames(abund_df), pattern = paste0('(',paste0(meta_sub$patient_ID,collapse='|'),')'))
    # merged <- dplyr::left_join(meta_sub, abund_df, by = "patient_ID")
    merged <- merge(x = meta_sub, y = abund_df, by = 'patient_ID')
    
    # Long format; remove NA freqs immediately
    abund_long <- merged %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq") %>%
      dplyr::mutate(entity = gsub("\\n", " ", entity)) %>%
      dplyr::filter(!is.na(freq))
    
    # Capture inputs at Generate time
    cont_var <- input$cont_group_var
    transpose_flag <- isTRUE(input$cont_transpose)
    test_padj_method <- input$cont_p_adj_method
    
    # Run Spearman per entity
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        if (is.null(cont_var) || !nzchar(cont_var) || !(cont_var %in% colnames(.x))) {
          return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
        }
        ok <- complete.cases(.x$freq, .x[[cont_var]])
        if (!any(ok)) return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
        ct <- suppressWarnings(cor.test(.x$freq[ok], .x[[cont_var]][ok], method = "spearman"))
        data.frame(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
      }) %>%
      dplyr::ungroup()
    
    # Adjust p-values if requested
    if (nrow(res) && nzchar(test_padj_method) && "p" %in% names(res)) {
      res$padj <- p.adjust(res$p, method = test_padj_method)
    }
    
    # Save info for export
    rv$last_cont_info <- list(
      entity = tolower(input$cont_entity %||% "clusters"),
      group = tolower(cont_var %||% "group"),
      test_raw = "spearman",
      transpose = transpose_flag
    )
    
    # Return everything needed for plotting
    list(
      data = abund_long,
      results = res,
      cont_var = cont_var,
      use_adj_p = input$cont_use_adj_p,
      facet_cols = as.numeric(input$cont_max_facets),
      transpose = transpose_flag
    )
  })
  
  observeEvent(input$generate_cont_plots, {
    cp <- cont_plot_data()    # your existing eventReactive
    cont_state(cp)
    output$cont_cleared_msg <- renderText(NULL)
  })
  
  observeEvent(input$reset_cont, {
    cont_state(NULL)          # clear the state
    cont_plot_cache(NULL)     # clear cached ggplot
    showNotification("Continuous plots cleared.", type = "message", duration = 5)
    output$cont_cleared_msg <- renderText("Results cleared. Generate new plots to see them here.")
  })
  
  output$hasContResults <- reactive({
    cp <- cont_state()
    !is.null(cp) && !is.null(cp$data) && nrow(cp$data) > 0
  })
  outputOptions(output, "hasContResults", suspendWhenHidden = FALSE)
  
  output$continuous_plot <- renderPlot({
    # cp <- cont_plot_data(); req(cp)
    cp <- cont_state(); req(cp)
    abund_long <- cp$data
    res <- cp$results
    cont_var <- cp$cont_var
    use_adj_p <- cp$use_adj_p
    facet_cols <- cp$facet_cols
    transpose_flag <- cp$transpose %||% FALSE
    
    # Validate continuous variable
    if (is.null(cont_var) || !nzchar(cont_var) || !(cont_var %in% colnames(abund_long))) {
      showNotification("No valid continuous metadata variable selected for continuous plotting.", type = "error")
      return(invisible(NULL))
    }
    
    # Remove rows with NA in either axis
    plot_df <- abund_long %>% dplyr::filter(!is.na(freq) & !is.na(.data[[cont_var]]))
    if (nrow(plot_df) == 0) {
      showNotification("No datapoints available for plotting after removing missing values.", type = "warning")
      return(invisible(NULL))
    }
    
    # Prepare p-value annotation dataframe
    p_df <- res %>%
      dplyr::mutate(
        p_to_show = if (isTRUE(use_adj_p) && "padj" %in% names(res)) padj else p,
        label = paste0("p = ", signif(p_to_show, 3), "\n", "rho = ", signif(rho, 3))
      )
    
    # Choose aesthetics and label positions
    if (!transpose_flag) {
      aes_pt <- ggplot2::aes(x = freq, y = .data[[cont_var]])
      smooth_aes <- ggplot2::aes(x = freq, y = .data[[cont_var]])
      x_lab <- "Abundance"; y_lab <- cont_var
      p_df <- p_df %>%
        dplyr::mutate(
          x = tapply(plot_df$freq, plot_df$entity, function(v) mean(range(v, na.rm = TRUE)))[entity],
          y = tapply(plot_df[[cont_var]], plot_df$entity, max, na.rm = TRUE)[entity] * 1.05
        )
    } else {
      aes_pt <- ggplot2::aes(x = .data[[cont_var]], y = freq)
      smooth_aes <- ggplot2::aes(x = .data[[cont_var]], y = freq)
      x_lab <- cont_var; y_lab <- "Abundance"
      p_df <- p_df %>%
        dplyr::mutate(
          x = tapply(plot_df[[cont_var]], plot_df$entity, function(v) mean(range(v, na.rm = TRUE)))[entity],
          y = tapply(plot_df$freq, plot_df$entity, max, na.rm = TRUE)[entity] * 1.05
        )
    }
    
    # Build ggplot with your preferred aesthetics
    gg <- ggplot2::ggplot(plot_df, mapping = aes_pt) +
      ggplot2::geom_point(alpha = 0.75, pch = 21, color = 'black', fill = 'grey40',
                          stroke = 0.1, size = 3) +
      ggplot2::geom_smooth(mapping = smooth_aes, method = "lm", se = FALSE, color = "red2") +
      ggplot2::facet_wrap(~entity, ncol = facet_cols, scales = "free") +
      ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      ggplot2::scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = x_lab, y = y_lab) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(
        margin = ggplot2::margin(t = 1.1, b = 1.1))
      )
    
    # Add p/rho annotations
    if (!is.null(p_df) && nrow(p_df) > 0) {
      p_df_plot <- p_df %>% dplyr::filter(entity %in% unique(plot_df$entity))
      if (nrow(p_df_plot) > 0) {
        gg <- gg + ggplot2::geom_text(
          data = p_df_plot,
          ggplot2::aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          size = 5,
          lineheight = 0.85
        )
      }
    }
    
    cont_plot_cache(gg)
    gg
  },
  height = function() {
    gg <- cont_plot_cache(); cp <- cont_plot_data()
    if (is.null(gg) || is.null(cp)) return(400)
    n_facets <- length(unique(gg$data$entity))
    ncol_facets <- cp$facet_cols
    if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
    nrow_facets <- ceiling(n_facets / ncol_facets)
    200 * nrow_facets
  },
  width = function() {
    gg <- cont_plot_cache(); cp <- cont_plot_data()
    if (is.null(gg) || is.null(cp)) return(400)
    ncol_facets <- cp$facet_cols
    if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
    225 * ncol_facets
  })
  
  output$export_cont_pdf <- downloadHandler(
    filename = function() {
      info <- rv$last_cont_info
      if (is.null(info)) return("continuous_plots.pdf")
      paste0("continuous_", info$entity, "_", info$group, "_spearman.pdf")
    },
    content = function(file) {
      gg <- cont_plot_cache()
      cp <- cont_plot_data()
      if (is.null(gg) || is.null(cp)) {
        showNotification("No plot available to export. Please generate the plots first.", type = "error")
        return()
      }
      n_facets <- length(unique(gg$data$entity))
      ncol_facets <- cp$facet_cols
      if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
      nrow_facets <- ceiling(n_facets / ncol_facets)
      if (!is.finite(nrow_facets) || nrow_facets < 1) nrow_facets <- 1
      pdf_width <- 4 * ncol_facets
      pdf_height <- 3 * nrow_facets
      ggsave(file, plot = gg, device = cairo_pdf, width = pdf_width, height = pdf_height, units = "in")
    },
    contentType = "application/pdf"
  )
  
  # Populate Feature Selection dropdowns from same metadata source as other tabs
  observeEvent(rv$meta_sample, {
    meta_cols <- sort(colnames(rv$meta_sample))
    updatePickerInput(session, "fs_outcome", choices = meta_cols)
    updatePickerInput(session, "fs_predictors", choices = c(meta_cols, "leiden_cluster"))
  }, ignoreInit = TRUE)
  
  observeEvent(rv$clusters$abundance, {
    cluster_names <- colnames(rv$clusters$abundance)
    updatePickerInput(session, "fs_leiden_subset", choices = cluster_names)
  })
  
  run_fs <- function() {
    req(rv$meta_sample, rv$abundance_sample, input$fs_outcome, input$fs_predictors)
    assert_sample_level(rv$meta_sample, "Feature Selection")
    
    # --- Unified predictor expansion ---
    expanded <- expand_predictors(
      meta_sample = rv$meta_sample,
      abundance_sample = rv$abundance_sample,
      outcome = input$fs_outcome,
      predictors = input$fs_predictors,
      cluster_subset = input$fs_leiden_subset,
      mode = "fs",
      notify = showNotification
    )
    
    merged <- expanded$df
    predictors_final <- expanded$predictors
    rv$fs_name_map <- expanded$map   # safe → original
    
    if (length(predictors_final) == 0) {
      showNotification("Select at least one predictor.", type = "error")
      return(NULL)
    }
    
    # Filter missingness
    merged <- merged[!is.na(merged[[input$fs_outcome]]), ]
    y_raw <- merged[[input$fs_outcome]]
    X_raw <- merged[, predictors_final, drop = FALSE]
    
    n_before <- nrow(merged)
    complete_rows <- complete.cases(data.frame(X_raw, .y = y_raw))
    dropped_ids <- merged$patient_ID[!complete_rows]
    
    X_raw <- X_raw[complete_rows, , drop = FALSE]
    y_raw <- y_raw[complete_rows]
    merged <- merged[complete_rows, , drop = FALSE]
    
    n_after <- sum(complete_rows)
    n_dropped <- n_before - n_after
    if (n_after < 3) {
      showNotification("Too few samples after filtering for feature selection.", type = "error")
      return(NULL)
    }
    
    y <- if (is.character(y_raw)) factor(y_raw) else y_raw
    
    # User controls
    method <- input$fs_method %||% "Elastic Net"
    tolerance <- 1e-4
    seed_val <- input$fs_seed %||% 123
    reps <- input$fs_reps %||% 1
    boruta_maxruns <- input$fs_maxruns %||% 500
    set.seed(seed_val)
    
    clean_dummy_names <- function(nms) {
      gsub(pattern = "leiden_cluster", replacement = "leiden_cluster:", x = nms)
    }
    
    # --- Boruta branch ---
    if (method == "Random Forest (Boruta)") {
      all_imp <- list()
      for (r in seq_len(reps)) {
        set.seed(seed_val + r - 1)
        X_boruta <- model.matrix(~ . - 1, data = X_raw)
        colnames(X_boruta) <- clean_dummy_names(colnames(X_boruta))
        X_boruta <- as.data.frame(X_boruta)
        
        if (is.factor(y)) {
          bor <- Boruta::Boruta(x = X_boruta, y = y, doTrace = 0, maxRuns = boruta_maxruns)
        } else {
          bor <- Boruta::Boruta(x = X_boruta, y = as.numeric(y), doTrace = 0, maxRuns = boruta_maxruns)
        }
        imp <- Boruta::attStats(bor)
        imp$meanImp[!is.finite(imp$meanImp)] <- 0
        all_imp[[r]] <- imp
      }
      merged_imp <- Reduce(function(a, b) {
        common <- intersect(rownames(a), rownames(b))
        a[common, "meanImp"] <- (a[common, "meanImp"] + b[common, "meanImp"]) / 2
        a
      }, all_imp)
      merged_imp <- merged_imp[abs(merged_imp$meanImp) >= tolerance, , drop = FALSE]
      
      sel <- rownames(merged_imp)[merged_imp$decision %in% c("Confirmed", "Tentative")]
      res_df <- data.frame(
        Feature = rownames(merged_imp),
        ImportanceMean = merged_imp$meanImp,
        Decision = merged_imp$decision,
        stringsAsFactors = FALSE
      )
      
      return(list(
        method = "Boruta",
        results = res_df,
        selected = sel,
        outcome = y,
        predictors = colnames(X_raw),
        tolerance = tolerance,
        details = list(samples_before = n_before,
                       samples_after = n_after,
                       samples_dropped = n_dropped,
                       dropped_ids = dropped_ids),
        merged = merged
      ))
    }
    
    # --- glmnet branch ---
    family <- if (is.factor(y)) {
      if (nlevels(y) == 2) "binomial" else "multinomial"
    } else {
      "gaussian"
    }
    
    Xmat <- model.matrix(~ . - 1, data = X_raw)
    colnames(Xmat) <- clean_dummy_names(colnames(Xmat))
    storage.mode(Xmat) <- "double"
    
    added_dummy <- FALSE
    if (ncol(Xmat) == 1) {
      Xmat <- cbind(Xmat, `__DUMMY__` = 0)
      added_dummy <- TRUE
    }
    
    alpha_val <- if (method == "Ridge Regression") 0 else (input$fs_alpha %||% 0.5)
    nfolds_val <- input$fs_nfolds %||% 5
    if (nfolds_val > n_after) nfolds_val <- max(3, floor(n_after / 2))
    
    coef_list <- list()
    for (r in seq_len(reps)) {
      set.seed(seed_val + r - 1)
      cvfit <- glmnet::cv.glmnet(
        x = Xmat, y = y, family = family,
        alpha = alpha_val, nfolds = nfolds_val
      )
      coef_obj <- coef(cvfit, s = "lambda.min")
      if (family == "multinomial") {
        coef_df <- do.call(rbind, lapply(names(coef_obj), function(cls) {
          cm <- as.matrix(coef_obj[[cls]])
          data.frame(Feature = rownames(cm), Class = cls, Coef = as.numeric(cm[, 1]), stringsAsFactors = FALSE)
        }))
        coef_df <- coef_df[coef_df$Feature != "(Intercept)", , drop = FALSE]
        if (added_dummy) coef_df <- coef_df[coef_df$Feature != "__DUMMY__", , drop = FALSE]
        coef_list[[r]] <- coef_df
      } else {
        cm <- as.matrix(coef_obj)
        coef_df <- data.frame(Feature = rownames(cm), Coef = as.numeric(cm[, 1]), stringsAsFactors = FALSE)
        coef_df <- coef_df[coef_df$Feature != "(Intercept)", , drop = FALSE]
        if (added_dummy) coef_df <- coef_df[coef_df$Feature != "__DUMMY__", , drop = FALSE]
        coef_list[[r]] <- coef_df
      }
    }
    
    if (family == "multinomial") {
      all_df <- do.call(rbind, coef_list)
      agg <- all_df %>%
        dplyr::group_by(Feature, Class) %>%
        dplyr::summarise(Coef = median(Coef, na.rm = TRUE), .groups = "drop")
      agg <- agg[abs(agg$Coef) >= tolerance, , drop = FALSE]
      res_df <- agg[order(agg$Coef, decreasing = TRUE), ]
      selected <- unique(head(res_df$Feature, 50))
    } else {
      all_df <- do.call(rbind, coef_list)
      agg <- all_df %>%
        dplyr::group_by(Feature) %>%
        dplyr::summarise(Coef = median(Coef, na.rm = TRUE), .groups = "drop")
      agg <- agg[abs(agg$Coef) >= tolerance, , drop = FALSE]
      res_df <- agg[order(agg$Coef, decreasing = TRUE), ]
      selected <- head(res_df$Feature, 50)
    }
    
    return(list(
      method = if (alpha_val == 0) "Ridge" else "Elastic Net",
      family = family,
      results = res_df,
      selected = selected,
      outcome = y,
      predictors = colnames(Xmat),
      tolerance = tolerance,
      details = list(samples_before = n_before,
                     samples_after = n_after,
                     samples_dropped = n_dropped,
                     dropped_ids = dropped_ids),
      merged = merged
    ))
  }
  
  output$fs_data_head <- renderTable({
    # res <- run_fs(); req(res)
    res <- fs_state(); req(res)
    head(res$merged, 5)
  }, sanitize.text.function = function(x) x)
  
  output$fs_plot <- renderPlot({
    # res <- run_fs(); req(res)
    res <- fs_state(); req(res)
    tol <- res$tolerance
    
    if (identical(res$method, "Boruta")) {
      df <- res$results
      df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No Boruta features pass the tolerance threshold.\nNothing to plot.", cex = 1.1)
        return()
      }
      df <- df[order(df$ImportanceMean, decreasing = TRUE), ]
      top_n <- head(df, 30)
      ggplot2::ggplot(top_n,
                      ggplot2::aes(x = reorder(Feature, ImportanceMean),
                                   y = ImportanceMean,
                                   fill = Decision)) +
        ggplot2::geom_col(color = 'black', lwd = 0.4) +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Boruta feature importance",
                      x = "Feature", y = "Mean importance") +
        ggplot2::scale_fill_manual(values = c('Confirmed' = 'green4',
                                              'Rejected' = 'red4',
                                              'Tentative' = 'grey40')) +
        ggplot2::theme_bw(base_size = 14)
    } else {
      df <- res$results
      df <- df[abs(df$Coef) >= tol, , drop = FALSE]
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "All coefficients are zero after regularization.\nTry lowering alpha or adding predictors.", cex = 1.1)
        return()
      }
      if ("Class" %in% names(df)) {
        top_n <- df %>%
          dplyr::group_by(Class) %>%
          dplyr::arrange(dplyr::desc(Coef), .by_group = TRUE) %>%
          dplyr::slice_head(n = 20)
        ggplot2::ggplot(top_n,
                        ggplot2::aes(x = reorder(Feature, Coef),
                                     y = Coef,
                                     fill = Coef)) +
          ggplot2::geom_col(color = 'black', lwd = 0.4) +
          ggplot2::facet_wrap(~Class, scales = "free_y") +
          ggplot2::coord_flip() +
          ggplot2::labs(title = paste(res$method, "coefficients (lambda.min)"),
                        x = "Feature", y = "Coefficient") +
          ggplot2::scale_fill_gradient(low = 'blue3', high = 'red3') +
          ggplot2::theme_bw(base_size = 14)
      } else {
        df <- df[order(df$Coef, decreasing = TRUE), ]
        top_n <- head(df, 30)
        ggplot2::ggplot(top_n,
                        ggplot2::aes(x = reorder(Feature, Coef),
                                     y = Coef,
                                     fill = Coef)) +
          ggplot2::geom_col(color = 'black', lwd = 0.4) +
          ggplot2::coord_flip() +
          ggplot2::labs(title = paste(res$method, "coefficients (lambda.min)"),
                        x = "Feature", y = "Coefficient") +
          ggplot2::scale_fill_gradient(low = 'blue3', high = 'red3') +
          ggplot2::theme_bw(base_size = 14)
      }
    }
  })
  
  output$fs_results <- renderTable({
    res <- fs_state(); req(res)
    df <- res$results
    tol <- res$tolerance
    
    if (identical(res$method, "Boruta")) {
      df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
      df <- df[order(df$ImportanceMean, decreasing = TRUE), , drop = FALSE]
      
      if (nrow(df) == 0) {
        return(data.frame(Message = "No Boruta features passed the tolerance threshold."))
      }
      
    } else {
      if ("Coef" %in% names(df)) {
        df <- df[abs(df$Coef) >= tol, , drop = FALSE]
        df <- df[order(df$Coef, decreasing = TRUE), , drop = FALSE]
        
        if (nrow(df) == 0) {
          return(data.frame(Message = "All coefficients shrank to zero after regularization."))
        }
      }
    }
    
    # Map safe names back to originals if available
    df <- map_safe_to_original(df, rv$fs_name_map)
    
    df
  }, sanitize.text.function = function(x) x)
  
  output$fs_summary <- renderPrint({
    # res <- run_fs(); req(res)
    res <- fs_state(); req(res)
    det <- res$details %||% list(samples_before = NA,
                                 samples_after = NA,
                                 samples_dropped = NA,
                                 dropped_ids = character(0))
    
    cat("Samples before filtering:", det$samples_before, "\n")
    cat("Samples after filtering:", det$samples_after, "\n")
    cat("Samples dropped:", det$samples_dropped, "\n\n")
    
    if (length(det$dropped_ids) > 0) {
      cat("Dropped patient_IDs:\n")
      print(det$dropped_ids)
      cat("\n")
    }
    
    cat("Selected features (top):\n")
    if (identical(res$method, "Boruta")) {
      df <- res$results
      df <- df[order(df$ImportanceMean, decreasing = TRUE), , drop = FALSE]
      if (nrow(df) == 0) {
        cat("No Boruta features passed the tolerance threshold.\n")
      } else {
        print(utils::head(df$Feature, 20))
      }
    } else {
      df <- res$results
      if ("Coef" %in% names(df)) {
        df <- df[order(df$Coef, decreasing = TRUE), , drop = FALSE]
        if (nrow(df) == 0) {
          cat("All coefficients shrank to zero after regularization.\n")
        } else {
          print(utils::head(df$Feature, 20))
        }
      }
    }
  })
  
  output$hasFSResults <- reactive({
    res <- fs_state()
    !is.null(res) && !is.null(res$results) && nrow(res$results) > 0
  })
  outputOptions(output, "hasFSResults", suspendWhenHidden = FALSE)
  
  output$export_fs_results <- downloadHandler(
    filename = function() {
      # res <- run_fs(); req(res)
      res <- fs_state(); req(res)
      method <- res$method %||% "FeatureSelection"
      outcome <- input$fs_outcome %||% "outcome"
      paste0(method, "_feature_selection_with_outcome_", outcome, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      # res <- run_fs(); req(res)
      res <- fs_state(); req(res)
      df <- res$results
      tol <- res$tolerance
      
      if (identical(res$method, "Boruta")) {
        df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
        df <- df[order(df$ImportanceMean, decreasing = TRUE), , drop = FALSE]
        if (nrow(df) == 0) {
          utils::write.csv(data.frame(Message = "No Boruta features passed the tolerance threshold."),
                           file, row.names = FALSE)
          return()
        }
      } else {
        if ("Coef" %in% names(df)) {
          df <- df[abs(df$Coef) >= tol, , drop = FALSE]
          df <- df[order(df$Coef, decreasing = TRUE), , drop = FALSE]
          if (nrow(df) == 0) {
            utils::write.csv(data.frame(Message = "All coefficients shrank to zero after regularization."),
                             file, row.names = FALSE)
            return()
          }
        }
      }
      utils::write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  observeEvent(list(rv$meta_sample, rv$abundance_sample), {
    req(rv$meta_sample, rv$abundance_sample)
    
    meta_cols <- colnames(rv$meta_sample)
    categorical_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.factor(x) || is.character(x))])
    predictor_choices <- sort(meta_cols)
    
    updatePickerInput(session, "lm_outcome", choices = categorical_choices, selected = NULL)
    updatePickerInput(session, "lm_predictors", choices = c(predictor_choices, "leiden_cluster"), selected = NULL)
    updatePickerInput(session, "lm_leiden_subset",
                      choices = colnames(rv$abundance_sample),
                      selected = character(0))
  }, ignoreInit = TRUE)
  
  observeEvent(input$reset_fs, {
    fs_state(NULL)
    showNotification("Feature Selection results cleared.", type = "message", duration = 5)
    output$fs_cleared_msg <- renderText("Results cleared. Run Feature Selection again to see results here.")
  })
  
  observeEvent(input$run_fs, {
    res <- run_fs()     # your existing eventReactive or function
    fs_state(res)
    output$fs_cleared_msg <- renderText(NULL)
  })
  
  run_lm <- function() {
    req(rv$meta_sample, rv$abundance_sample, input$lm_outcome, input$lm_predictors)
    assert_sample_level(rv$meta_sample, "Classification")
    
    meta_patient <- rv$meta_sample
    
    # Outcome
    if (!(input$lm_outcome %in% colnames(meta_patient))) {
      showNotification(paste0("Outcome '", input$lm_outcome, "' not found in metadata."), type = "error")
      return(NULL)
    }
    outcome_raw <- meta_patient[[input$lm_outcome]]
    outcome <- factor(outcome_raw)
    
    # Sanitize levels
    orig_levels <- levels(outcome)
    safe_levels <- make.names(orig_levels)
    levels(outcome) <- safe_levels
    rv$lm_label_map <- setNames(orig_levels, safe_levels)
    
    if (nlevels(outcome) < 2) {
      showNotification("Outcome must have ≥2 classes.", type = "error")
      return(NULL)
    }
    
    # --- Unified predictor expansion ---
    expanded <- expand_predictors(
      meta_sample = rv$meta_sample,
      abundance_sample = rv$abundance_sample,
      outcome = input$lm_outcome,
      predictors = input$lm_predictors,
      cluster_subset = input$lm_leiden_subset,
      mode = "classification",
      notify = showNotification
    )
    
    df <- expanded$df
    predictors_final <- expanded$predictors
    rv$lm_name_map <- expanded$map   # safe → original
    
    if (length(predictors_final) == 0) {
      showNotification("Select at least one predictor.", type = "error")
      return(NULL)
    }
    
    # Filter missingness
    df <- df[!is.na(df[[input$lm_outcome]]), ]
    X <- df[, predictors_final, drop = FALSE]
    y <- factor(df[[input$lm_outcome]], levels = levels(outcome))
    
    n_before <- nrow(df)
    complete_rows <- complete.cases(data.frame(X, .y = y))
    dropped_ids <- df$patient_ID[!complete_rows]
    
    X <- X[complete_rows, , drop = FALSE]
    y <- y[complete_rows]
    df <- df[complete_rows, , drop = FALSE]
    
    n_after <- sum(complete_rows)
    n_dropped <- n_before - n_after
    if (n_after < 10) {
      showNotification("Too few samples after filtering. Model not run.", type = "error")
      return(NULL)
    }
    
    # Null model baseline
    majority_class <- names(which.max(table(y)))
    null_acc <- mean(y == majority_class)
    
    # Choose model method
    method <- switch(input$lm_model_type,
                     "Logistic Regression" = "multinom",
                     "Elastic Net"         = "glmnet",
                     "Random Forest"       = "rf")
    cat("Running model type:", method, "\n")
    
    validation <- input$lm_validation
    split_info <- NULL
    
    # --- Train/Test split ---
    if (validation == "Train/Test split") {
      set.seed(123)
      train_frac <- input$lm_train_frac %||% 0.7
      idx <- caret::createDataPartition(y, p = train_frac, list = FALSE)
      trainX <- X[idx, , drop = FALSE]
      testX  <- X[-idx, , drop = FALSE]
      trainY <- y[idx]
      testY  <- y[-idx]
      
      split_info <- list(
        n_train = length(trainY),
        n_test  = length(testY),
        train_counts = table(trainY),
        test_counts  = table(testY)
      )
      
      if (method == "glmnet") {
        trainMat <- model.matrix(~ . - 1, data = trainX)
        testMat  <- model.matrix(~ . - 1, data = testX)
        storage.mode(trainMat) <- "double"
        storage.mode(testMat)  <- "double"
        alpha_val <- input$lm_alpha %||% 0.5
        model <- caret::train(
          x = trainMat, y = trainY,
          method = "glmnet",
          trControl = caret::trainControl(classProbs = TRUE, verboseIter = TRUE),
          tuneGrid = expand.grid(alpha = alpha_val,
                                 lambda = 10^seq(-3, 1, length = 20))
        )
        probs <- predict(model, newdata = testMat, type = "prob")
        pred_class <- predict(model, newdata = testMat, type = "raw")
      } else {
        suppressWarnings(
          model <- caret::train(
            x = trainX, y = trainY,
            method = method,
            trControl = caret::trainControl(classProbs = TRUE, verboseIter = TRUE)
          )
        )
        probs <- predict(model, newdata = testX, type = "prob")
        pred_class <- predict(model, newdata = testX, type = "raw")
      }
      
      # Build preds_tbl
      if (is.vector(probs)) {
        pos <- levels(trainY)[2]
        probs <- data.frame(setNames(list(as.numeric(probs)), pos), check.names = FALSE)
      }
      preds_tbl <- dplyr::bind_cols(
        truth = testY,
        dplyr::rename(probs, !!!setNames(names(probs), paste0(".pred_", names(probs)))),
        pred = pred_class
      )
      names(preds_tbl)[1] <- "obs"
      preds_tbl$obs <- factor(preds_tbl$obs, levels = levels(y))
      
      if (!"pred" %in% names(preds_tbl)) {
        prob_cols <- grep("^\\.pred_", names(preds_tbl), value = TRUE)
        classes <- gsub("^\\.pred_", "", prob_cols)
        max_idx <- apply(as.matrix(preds_tbl[, prob_cols, drop = FALSE]), 1, which.max)
        preds_tbl$pred <- factor(classes[max_idx], levels = classes)
      }
      
      return(list(model = model, preds = preds_tbl, null_acc = null_acc,
                  split_info = split_info,
                  details = list(samples_before = n_before,
                                 samples_after = n_after,
                                 samples_dropped = n_dropped,
                                 dropped_ids = dropped_ids)))
    }
    
    # --- CV ---
    ctrl <- caret::trainControl(
      method = "cv",
      number = input$lm_k,
      classProbs = TRUE,
      savePredictions = "final",
      verboseIter = TRUE
    )
    
    if (method == "glmnet") {
      Xmat <- model.matrix(~ . - 1, data = X)
      storage.mode(Xmat) <- "double"
      alpha_val <- input$lm_alpha %||% 0.5
      model <- caret::train(
        x = Xmat, y = y,
        method = "glmnet",
        trControl = ctrl,
        tuneGrid = expand.grid(alpha = alpha_val,
                               lambda = 10^seq(-3, 1, length = 20))
      )
    } else {
      suppressWarnings(
        model <- caret::train(
          x = X, y = y,
          method = method,
          trControl = ctrl
        )
      )
    }
    
    preds <- model$pred
    if (is.null(preds) || !nrow(preds)) {
      showNotification("No predictions available from resampling. Try Train/Test split.", type = "error")
      return(NULL)
    }
    
    level_names <- levels(y)
    prob_cols <- intersect(level_names, names(preds))
    rename_map <- setNames(prob_cols, paste0(".pred_", prob_cols))
    preds <- dplyr::rename(preds, !!!rename_map)
    preds$obs <- factor(preds$obs, levels = levels(y))
    
    if (!"pred" %in% names(preds)) {
      prob_cols <- grep("^\\.pred_", names(preds), value = TRUE)
      classes <- gsub("^\\.pred_", "", prob_cols)
      max_idx <- apply(as.matrix(preds[, prob_cols, drop = FALSE]), 1, which.max)
      preds$pred <- factor(classes[max_idx], levels = classes)
    }
    
    return(list(model = model, preds = preds, null_acc = null_acc,
                details = list(samples_before = n_before,
                               samples_after = n_after,
                               samples_dropped = n_dropped,
                               dropped_ids = dropped_ids)))
  }
  
  observeEvent(input$run_lm, {
    res <- run_lm()
    req(res)
    
    preds <- res$preds
    req(preds)
    
    n_classes <- nlevels(preds$obs)
    roc_plot <- NULL
    
    if (n_classes == 2) {
      # Binary ROC with pROC
      classes <- levels(preds$obs)
      positive <- classes[2]
      roc_obj <- pROC::roc(
        response = preds$obs,
        predictor = preds[[paste0(".pred_", positive)]],
        levels = rev(classes)
      )
      roc_plot <- ggplot2::ggplot(data.frame(
        fpr = 1 - roc_obj$specificities,
        tpr = roc_obj$sensitivities
      ), ggplot2::aes(x = fpr, y = tpr)) +
        ggplot2::geom_line(color = "blue", linewidth = 1) +
        ggplot2::geom_abline(linetype = "dashed", color = "grey50") +
        ggplot2::labs(title = sprintf("Binary ROC Curve (AUC = %.3f)", pROC::auc(roc_obj)),
                      x = "False Positive Rate", y = "True Positive Rate") +
        ggplot2::theme_minimal(base_size = 14)
    } else {
      # Multiclass ROC with yardstick
      long_preds <- preds %>%
        dplyr::select(obs, starts_with(".pred_"))
      roc_curves <- yardstick::roc_curve(long_preds, truth = obs, dplyr::starts_with(".pred_"))
      roc_plot <- ggplot2::ggplot(roc_curves,
                                  ggplot2::aes(x = 1 - specificity, y = sensitivity, color = .level)) +
        ggplot2::geom_path(linewidth = 1) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
        ggplot2::labs(title = "Multiclass ROC Curves (yardstick)", color = "Class") +
        ggplot2::theme_bw(base_size = 14)
    }
    
    # Store everything in lm_state, including cached ROC plot
    lm_state(c(res, list(roc_plot = roc_plot)))
    
    output$lm_cleared_msg <- renderText(NULL)
  })
  
  observeEvent(input$reset_lm, {
    lm_state(NULL)
    showNotification("Logistic Modeling results cleared.", type = "message", duration = 5)
    output$lm_cleared_msg <- renderText("Results cleared. Run Logistic Modeling again to see results here.")
  })
  
  output$lm_roc_plot <- renderPlot({
    s <- lm_state()
    req(s, s$roc_plot)
    s$roc_plot
  })
  
  build_perf_table <- function(res, rv) {
    req(res, res$preds)
    preds <- res$preds
    req(preds)
    
    n_classes <- nlevels(preds$obs)
    
    obs <- preds$obs
    pred <- preds$pred
    if (is.null(pred)) {
      pred_cols <- grep("^\\.pred_", names(preds), value = TRUE)
      classes <- gsub("^\\.pred_", "", pred_cols)
      max_idx <- apply(as.matrix(preds[, pred_cols, drop = FALSE]), 1, which.max)
      pred <- factor(classes[max_idx], levels = classes)
    }
    
    model_acc <- mean(pred == obs, na.rm = TRUE)
    successes <- sum(pred == obs, na.rm = TRUE)
    n <- length(obs)
    
    binom_p <- tryCatch({
      binom.test(successes, n, p = res$null_acc, alternative = "greater")$p.value
    }, error = function(e) NA)
    
    majority_safe <- names(which.max(table(obs)))
    majority_orig <- rv$lm_label_map[[majority_safe]]
    
    perf_table <- data.frame(
      Metric = c("Null Accuracy (most frequent class)", "Model Accuracy", "Binomial p-value vs Null"),
      Value  = c(paste0(round(res$null_acc, 3), " (", majority_orig, ")"),
                 round(model_acc, 3),
                 signif(binom_p, 3)),
      stringsAsFactors = FALSE
    )
    
    if (n_classes == 2) {
      classes <- levels(preds$obs)
      positive <- classes[2]
      roc_obj <- pROC::roc(response = preds$obs,
                           predictor = preds[[paste0(".pred_", positive)]],
                           levels = rev(classes))
      auc_val <- as.numeric(pROC::auc(roc_obj))
      perf_table <- rbind(perf_table,
                          data.frame(Metric = "AUC (pROC)", Value = round(auc_val, 3), stringsAsFactors = FALSE))
    } else {
      long_preds <- preds %>%
        dplyr::select(obs, starts_with(".pred_"))
      
      auc_macro <- yardstick::roc_auc(long_preds, truth = obs,
                                      dplyr::starts_with(".pred_"), estimator = "macro")
      auc_macro_wt <- yardstick::roc_auc(long_preds, truth = obs,
                                         dplyr::starts_with(".pred_"), estimator = "macro_weighted")
      auc_ht <- yardstick::roc_auc(long_preds, truth = obs,
                                   dplyr::starts_with(".pred_"), estimator = "hand_till")
      
      perf_table <- rbind(
        perf_table,
        data.frame(Metric = "Macro AUC (yardstick)", Value = round(auc_macro$.estimate, 3), stringsAsFactors = FALSE),
        data.frame(Metric = "Macro-weighted AUC (yardstick)", Value = round(auc_macro_wt$.estimate, 3), stringsAsFactors = FALSE),
        data.frame(Metric = "Hand-Till AUC (yardstick)", Value = round(auc_ht$.estimate, 3), stringsAsFactors = FALSE)
      )
    }
    
    if (!is.null(res$split_info)) {
      train_summary <- paste(names(res$split_info$train_counts),
                             res$split_info$train_counts, collapse = "; ")
      test_summary  <- paste(names(res$split_info$test_counts),
                             res$split_info$test_counts,  collapse = "; ")
      extra_rows <- data.frame(
        Metric = c("Train size", "Test size", "Train breakdown", "Test breakdown"),
        Value  = c(res$split_info$n_train, res$split_info$n_test, train_summary, test_summary),
        stringsAsFactors = FALSE
      )
      perf_table <- rbind(perf_table, extra_rows)
    }
    
    det <- res$details
    extra_rows2 <- data.frame(
      Metric = c("Samples before filtering", "Samples after filtering", "Samples dropped"),
      Value  = c(det$samples_before, det$samples_after, det$samples_dropped),
      stringsAsFactors = FALSE
    )
    perf_table <- rbind(perf_table, extra_rows2)
    if (length(det$dropped_ids) > 0) {
      perf_table <- rbind(perf_table,
                          data.frame(Metric = "Dropped patient_IDs",
                                     Value = paste(unique(det$dropped_ids), collapse = ", "),
                                     stringsAsFactors = FALSE))
    }
    
    perf_table
  }
  
  output$lm_perf_table <- renderTable({
    res <- lm_state(); req(res)
    build_perf_table(res, rv)
  })
  
  output$lm_summary <- renderPrint({
    # res <- lm_results(); req(res)
    res <- lm_state(); req(res)
    det <- res$details %||% list(samples_before = NA,
                                 samples_after = NA,
                                 samples_dropped = NA,
                                 dropped_ids = character(0))
    
    cat("Samples before filtering:", det$samples_before, "\n")
    cat("Samples after filtering:", det$samples_after, "\n")
    cat("Samples dropped:", det$samples_dropped, "\n\n")
    
    if (length(det$dropped_ids) > 0) {
      cat("Dropped patient_IDs:\n")
      print(det$dropped_ids)
      cat("\n")
    }
    
    # Outcome distribution
    preds <- res$preds
    if (!is.null(preds) && "obs" %in% names(preds)) {
      cat("Outcome distribution (observed):\n")
      print(table(preds$obs))
      cat("\n")
    }
    
    # Train/Test breakdown if available
    if (!is.null(res$split_info)) {
      cat("Train breakdown:\n")
      print(res$split_info$train_counts)
      cat("\nTest breakdown:\n")
      print(res$split_info$test_counts)
      cat("\n")
    }
    
    cat("Null accuracy (baseline):", round(res$null_acc, 3), "\n")
  })
  
  order_features <- function(df, raw_col) {
    if (is.null(df) || !raw_col %in% colnames(df)) return(df)
    
    # Always work with a plain data.frame
    df <- as.data.frame(df)
    
    # Flag intercept rows
    intercept_row <- grepl("\\(Intercept\\)", df$feature, fixed = TRUE)
    df_intercept <- df[intercept_row, , drop = FALSE]
    df_notintercept <- df[!intercept_row, , drop = FALSE]
    
    # Sort non‑intercepts by absolute raw value
    if (nrow(df_notintercept) > 0) {
      raw_vals_notint <- as.numeric(df_notintercept[[raw_col]])
      df_notintercept <- df_notintercept[order(abs(raw_vals_notint), decreasing = TRUE), , drop = FALSE]
    }
    
    # Sort intercepts (in case there are multiple)
    if (nrow(df_intercept) > 0) {
      raw_vals_int <- as.numeric(df_intercept[[raw_col]])
      df_intercept <- df_intercept[order(abs(raw_vals_int), decreasing = TRUE), , drop = FALSE]
      df <- rbind(df_intercept, df_notintercept)
    } else {
      df <- df_notintercept
    }
    
    rownames(df) <- NULL
    df
  }
  
  # --- Extract coefficients for glm / glmnet ---
  # Coefficients for glmnet, glm, multinom with consistent schema and column order
  coef_table <- reactive({
    s <- lm_state()
    req(s, s$model)
    
    if (s$model$method == "glmnet") {
      coefs <- tryCatch({
        coef(s$model$finalModel, s = s$model$bestTune$lambda)
      }, error = function(e) NULL)
      if (is.null(coefs)) return(NULL)
      
      if (is.list(coefs)) {
        # Multiclass glmnet: list of matrices, each matrix rows = predictors
        df_list <- lapply(names(coefs), function(cls) {
          mat <- as.matrix(coefs[[cls]])
          data.frame(
            class = cls,                     # outcome category
            feature = rownames(mat),         # predictor
            ScaledCoefficient = as.numeric(mat),
            stringsAsFactors = FALSE
          )
        })
        df <- do.call(rbind, df_list)
      } else {
        # Binary glmnet: single matrix, rows = predictors
        mat <- as.matrix(coefs)
        df <- data.frame(
          feature = rownames(mat),
          ScaledCoefficient = as.numeric(mat),
          stringsAsFactors = FALSE
        )
      }
      
      # Back-transform to unscaled coefficients if caret preProcess was used
      if (!is.null(s$model$preProcess)) {
        pp <- s$model$preProcess
        sds <- pp$std
        df$RawCoefficient <- df$ScaledCoefficient
        for (feat in intersect(names(sds), df$feature)) {
          idx <- df$feature == feat
          df$RawCoefficient[idx] <- df$ScaledCoefficient[idx] / sds[feat]
        }
      } else {
        df$RawCoefficient <- df$ScaledCoefficient
      }
      
      # Enforce column order
      if ("class" %in% names(df)) {
        df <- df[, c("class", "feature", "RawCoefficient", "ScaledCoefficient")]
      } else {
        df <- df[, c("feature", "RawCoefficient", "ScaledCoefficient")]
      }
      df <- order_features(df, "RawCoefficient")
      rownames(df) <- NULL
      df
    } else if (s$model$method %in% c("glm", "multinom")) {
      coefs <- coef(s$model$finalModel)
      
      if (is.matrix(coefs)) {
        # Multiclass multinom: rows = outcome classes, columns = predictors
        df <- as.data.frame(coefs)
        df$class <- rownames(coefs)  # outcome categories from row names
        
        df_long <- tidyr::pivot_longer(
          df,
          cols = setdiff(names(df), "class"),
          names_to = "feature",          # predictors from column names
          values_to = "RawCoefficient"
        )
        
        df_long$ScaledCoefficient <- df_long$RawCoefficient
        df <- df_long[, c("class", "feature", "RawCoefficient", "ScaledCoefficient")]
      } else {
        # Binary logistic regression: named vector of coefficients (predictors)
        df <- data.frame(
          feature = names(coefs),
          RawCoefficient = unname(coefs),
          stringsAsFactors = FALSE
        )
        df$ScaledCoefficient <- df$RawCoefficient
        df <- df[, c("feature", "RawCoefficient", "ScaledCoefficient")]
      }
      df <- order_features(df, "RawCoefficient")
      rownames(df) <- NULL
      df
    } else {
      NULL
    }
  })
  
  # --- Extract variable importance for random forest ---
  rf_importance <- reactive({
    s <- lm_state()
    req(s, s$model)
    
    if (s$model$method == "rf") {
      imp_scaled <- caret::varImp(s$model, scale = TRUE)$importance
      imp_scaled$feature <- rownames(imp_scaled)
      colnames(imp_scaled)[1] <- "ScaledImportance"
      
      rf_model <- s$model$finalModel
      imp_raw <- randomForest::importance(rf_model)
      if ("%IncMSE" %in% colnames(imp_raw)) {
        raw_vals <- imp_raw[, "%IncMSE"]
      } else {
        raw_vals <- imp_raw[, 1]
      }
      imp_raw_df <- data.frame(feature = rownames(imp_raw),
                               RawImportance = raw_vals,
                               stringsAsFactors = FALSE)
      
      df <- dplyr::left_join(imp_scaled, imp_raw_df, by = "feature")
      df <- df[, c("feature", "RawImportance", "ScaledImportance")]
      df <- order_features(df, "RawImportance")
      rownames(df) <- NULL
      df
    } else {
      NULL
    }
  })
  
  output$lm_features <- renderTable({
    s <- lm_state()
    req(s, s$model)
    
    if (s$model$method %in% c("glm", "glmnet", "multinom")) {
      df <- coef_table()
    } else if (s$model$method == "rf") {
      df <- rf_importance()
    } else {
      df <- data.frame(Message = "Feature importance not available for this model type.")
    }
    
    # Map safe names back to originals if available
    df <- map_safe_to_original(df, rv$lm_name_map)
    
    df
  }, digits = 3)
  
  output$hasLMResults <- reactive({
    res <- lm_state()
    !is.null(res) && !is.null(res$preds) && nrow(res$preds) > 0
  })
  outputOptions(output, "hasLMResults", suspendWhenHidden = FALSE)
  
  output$export_lm_zip <- downloadHandler(
    filename = function() {
      model <- gsub("\\s+", "_", tolower(input$lm_model_type %||% "model"))
      validation <- gsub("\\s+", "_", tolower(input$lm_validation %||% "validation"))
      outcome <- gsub("\\s+", "_", tolower(input$lm_outcome %||% "outcome"))
      paste0(model, "_", validation, "_", outcome, ".zip")
    },
    content = function(file) {
      s <- lm_state()
      req(s, s$model)
      
      tmpdir <- tempdir()
      files <- c()
      
      # 1. Model Summary (key-value CSV)
      summary_file <- file.path(tmpdir, "model_summary.csv")
      summary_kv <- data.frame(
        Key = c("Model type", "Outcome variable", "Predictors",
                "Validation strategy", "Null accuracy",
                "Samples before filtering", "Samples after filtering", "Samples dropped"),
        Value = c(input$lm_model_type,
                  input$lm_outcome,
                  paste(input$lm_predictors, collapse = "; "),
                  input$lm_validation,
                  sprintf("%.3f", s$null_acc %||% NA),
                  s$details$samples_before %||% NA,
                  s$details$samples_after %||% NA,
                  s$details$samples_dropped %||% NA),
        stringsAsFactors = FALSE
      )
      write.csv(summary_kv, summary_file, row.names = FALSE)
      files <- c(files, summary_file)
      
      # 2. Performance Metrics (exactly as in UI)
      perf_file <- file.path(tmpdir, "performance_metrics.csv")
      perf_df <- tryCatch({
        build_perf_table(s, rv)
      }, error = function(e) data.frame(Message = "Error extracting performance metrics"))
      write.csv(perf_df, perf_file, row.names = FALSE)
      files <- c(files, perf_file)
      
      # 3. Model Features (coefficients / importance with enforced column order)
      feat_file <- file.path(tmpdir, "model_features.csv")
      feat_df <- NULL
      if (s$model$method %in% c("glm", "glmnet", "multinom")) {
        feat_df <- coef_table()
      } else if (s$model$method == "rf") {
        feat_df <- rf_importance()
      }
      
      if (!is.null(feat_df) && nrow(feat_df) > 0) {
        # Enforce column order
        if (all(c("class","feature","RawCoefficient","ScaledCoefficient") %in% names(feat_df))) {
          feat_df <- feat_df[, c("class","feature","RawCoefficient","ScaledCoefficient")]
        } else if (all(c("feature","RawCoefficient","ScaledCoefficient") %in% names(feat_df))) {
          feat_df <- feat_df[, c("feature","RawCoefficient","ScaledCoefficient")]
        } else if (all(c("feature","RawImportance","ScaledImportance") %in% names(feat_df))) {
          feat_df <- feat_df[, c("feature","RawImportance","ScaledImportance")]
        }
        write.csv(feat_df, feat_file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No feature importance available"),
                  feat_file, row.names = FALSE)
      }
      files <- c(files, feat_file)
      
      # 4. ROC Curve (PDF)
      roc_file <- file.path(tmpdir, "roc_curve.pdf")
      if (!is.null(s$roc_plot)) {
        ggsave(roc_file, plot = s$roc_plot,
               device = cairo_pdf, width = 6, height = 6, units = "in")
      } else {
        pdf(roc_file); plot.new(); text(0.5, 0.5, "ROC plot not available"); dev.off()
      }
      files <- c(files, roc_file)
      
      # Bundle into zip
      zip::zip(zipfile = file, files = files, mode = "cherry-pick")
    },
    contentType = "application/zip"
  )
  
  # Populate Regression tab dropdowns from sample-level metadata
  observe({
    req(rv$meta_sample)
    meta_cols <- colnames(rv$meta_sample)
    
    # Continuous outcomes: numeric or integer metadata
    continuous_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.numeric(x) || is.integer(x))])
    
    # Predictors: all sample-level metadata columns (same as classification tab)
    predictor_choices <- sort(meta_cols)
    
    updatePickerInput(session, "reg_outcome",
                      choices = continuous_choices,
                      selected = if (length(continuous_choices) > 0) continuous_choices[1])
    
    updatePickerInput(session, "reg_predictors",
                      choices = c(predictor_choices, "leiden_cluster"),
                      selected = NULL)
  })
  
  reg_results <- eventReactive(input$run_reg, {
    req(rv$meta_sample, input$reg_outcome, input$reg_predictors)
    
    expanded <- expand_predictors(
      meta_sample = rv$meta_sample,
      abundance_sample = rv$abundance_sample,
      outcome = input$reg_outcome,
      predictors = input$reg_predictors,
      cluster_subset = input$reg_leiden_subset,
      mode = "regression",
      encode_factors = isTRUE(input$reg_allow_categorical),
      notify = showNotification
    )
    
    df <- expanded$df
    preds <- expanded$predictors
    rv$reg_name_map <- expanded$map   # safe → original
    
    # Outcome must be continuous with variance > 0
    y <- suppressWarnings(as.numeric(df[[input$reg_outcome]]))
    if (!is.finite(sd(y, na.rm = TRUE)) || sd(y, na.rm = TRUE) == 0 || length(y) < 2) {
      stop("Outcome has insufficient variability or too few samples after merging/filters.")
    }
    
    run_regression_model(
      model_type = input$reg_model_type,
      df = df,
      outcome = input$reg_outcome,
      predictors = preds,
      validation = switch(input$reg_validation,
                          "split" = "train_test",
                          "cv"    = "cv"),
      p = input$reg_train_frac,
      k = input$reg_k,
      alpha = input$reg_alpha,
      n_perm = input$reg_n_perm %||% 100
    )
  })
  
  observeEvent(reg_results(), {
    reg_state(reg_results())
  })
  
  # output$reg_summary <- renderPrint({
  #   res <- reg_state()
  #   if (is.null(res)) return("No results yet.")
  #   res$metrics
  # })
  
  output$reg_perf_table <- renderTable({
    res <- reg_state()
    if (is.null(res)) return(NULL)
    as.data.frame(res$metrics)
  })
  
  output$reg_pred_plot <- renderPlot({
    res <- reg_state()
    req(res$performance_plot)
    res$performance_plot
  })
  
  output$reg_resid_plot <- renderPlot({
    res <- reg_state()
    req(res$diagnostics_plot)
    res$diagnostics_plot
  })
  
  output$reg_features <- renderTable({
    res <- reg_state()
    req(res)
    
    df <- res$feature_importance %||% data.frame()
    df <- map_safe_to_original(df, res$name_map)
    
    df
  }, sanitize.text.function = function(x) x)
  
  output$hasRegResults <- reactive({ !is.null(reg_state()) })
  outputOptions(output, "hasRegResults", suspendWhenHidden = FALSE)
  
  # Clear Regression results
  observeEvent(input$reset_reg, {
    reg_state(NULL)   # clear the stored results
    
    showNotification("Regression results cleared.", type = "message", duration = 5)
    
    output$reg_cleared_msg <- renderText(
      "Results cleared. Run a new regression to see results here."
    )
  })
  
  observeEvent(reg_state(), {
    if (!is.null(reg_state())) {
      output$reg_cleared_msg <- renderText(NULL)
    }
  })
  
  output$export_reg_zip <- downloadHandler(
    filename = function() {
      model <- gsub("\\s+", "_", tolower(input$reg_model_type %||% "model"))
      validation <- gsub("\\s+", "_", tolower(input$reg_validation %||% "validation"))
      outcome <- gsub("\\s+", "_", tolower(input$reg_outcome %||% "outcome"))
      paste0("regression_", model, "_", validation, "_", outcome, ".zip")
    },
    content = function(file) {
      s <- reg_state(); req(s)
      
      tmpdir <- tempdir()
      files <- c()
      
      # 1. Model Summary
      # summary_file <- file.path(tmpdir, "reg_summary.csv")
      # write.csv(broom::glance(s$model), summary_file, row.names = FALSE)
      # files <- c(files, summary_file)
      
      # 2. Performance Metrics
      perf_file <- file.path(tmpdir, "reg_performance.csv")
      write.csv(s$perf, perf_file, row.names = FALSE)
      files <- c(files, perf_file)
      
      # 3. Model Features
      feat_file <- file.path(tmpdir, "reg_features.csv")
      feat_df <- NULL
      method <- input$reg_model_type
      if (method %in% c("lm", "glmnet")) {
        df <- broom::tidy(s$model$finalModel)
        df <- df[, c("term", "estimate")]
        colnames(df) <- c("feature_safe", "RawCoefficient")
        df$feature <- s$safe_to_orig[df$feature_safe] %||% df$feature_safe
        df$ScaledCoefficient <- df$RawCoefficient
        feat_df <- order_features(df[, c("feature", "RawCoefficient", "ScaledCoefficient")],
                                  "RawCoefficient")
      } else if (method == "rf") {
        imp <- caret::varImp(s$model)$importance
        imp$feature_safe <- rownames(imp)
        imp$feature <- s$safe_to_orig[imp$feature_safe] %||% imp$feature_safe
        colnames(imp)[1] <- "ScaledImportance"
        rf_imp <- randomForest::importance(s$model$finalModel)
        raw_vals <- if ("%IncMSE" %in% colnames(rf_imp)) rf_imp[, "%IncMSE"] else rf_imp[, 1]
        imp$RawImportance <- raw_vals[match(imp$feature_safe, rownames(rf_imp))]
        feat_df <- imp[, c("feature", "RawImportance", "ScaledImportance")]
        feat_df <- order_features(feat_df, "RawImportance")
      } else if (method %in% c("gbm", "xgbTree")) {
        imp <- caret::varImp(s$model)$importance
        imp$feature_safe <- rownames(imp)
        imp$feature <- s$safe_to_orig[imp$feature_safe] %||% imp$feature_safe
        colnames(imp)[1] <- "ScaledImportance"
        imp$RawImportance <- imp$ScaledImportance
        feat_df <- imp[, c("feature", "RawImportance", "ScaledImportance")]
        feat_df <- order_features(feat_df, "RawImportance")
      }
      
      if (!is.null(feat_df) && nrow(feat_df) > 0) {
        write.csv(feat_df, feat_file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No feature importance available"),
                  feat_file, row.names = FALSE)
      }
      files <- c(files, feat_file)
      
      # 4. Predicted vs Observed plot
      pred_plot_file <- file.path(tmpdir, "pred_vs_obs.pdf")
      pdf(pred_plot_file)
      df <- data.frame(obs = s$obs, preds = s$preds)
      print(
        ggplot(df, aes(x = obs, y = preds)) +
          geom_point(alpha = 0.6) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
          labs(x = "Observed", y = "Predicted") +
          theme_minimal()
      )
      dev.off()
      files <- c(files, pred_plot_file)
      
      # 5. Residual plot
      resid_plot_file <- file.path(tmpdir, "residuals.pdf")
      pdf(resid_plot_file)
      resid <- s$obs - s$preds
      df <- data.frame(preds = s$preds, resid = resid)
      print(
        ggplot(df, aes(x = preds, y = resid)) +
          geom_point(alpha = 0.6) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          labs(x = "Predicted", y = "Residuals") +
          theme_minimal()
      )
      dev.off()
      files <- c(files, resid_plot_file)
      
      # 6. Partial Dependence Plot (tree-based only)
      if (method %in% c("rf", "gbm", "xgbTree")) {
        pdp_file <- file.path(tmpdir, "partial_dependence.pdf")
        pdf(pdp_file)
        feat_orig <- input$reg_feature_pdp
        if (!is.null(feat_orig) && feat_orig %in% names(s$safe_to_orig)) {
          feat_safe <- names(s$safe_to_orig)[match(feat_orig, s$safe_to_orig)]
          pd <- pdp::partial(s$model, pred.var = feat_safe, train = s$data)
          print(autoplot(pd) + theme_minimal() +
                  labs(title = paste("Partial Dependence:", feat_orig)))
        }
        dev.off()
        files <- c(files, pdp_file)
        
        # 7. SHAP Plot (example for first observation)
        shap_file <- file.path(tmpdir, "shap_plot.pdf")
        pdf(shap_file)
        X <- s$data[, s$safe_preds, drop = FALSE]
        predictor <- iml::Predictor$new(s$model, data = X, y = s$data[[s$safe_outcome]])
        shap <- iml::Shapley$new(predictor, x.interest = X[1, , drop = FALSE])
        plot(shap)
        dev.off()
        files <- c(files, shap_file)
      }
      
      # Bundle into zip
      zip::zip(zipfile = file, files = files, mode = "cherry-pick")
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)

