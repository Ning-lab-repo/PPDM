# ==============================================================================
# Multi-Omics Mediation Analysis for Physical Activity and Disease
# 
# Description: This script performs structural equation modeling (SEM) to 
# investigate the mediation effects of metabolomics and inflammatory markers 
# in the relationship between physical activity and disease risk. It includes
# Common Antecedent Factor (CAF) analysis and bidirectional mediation modeling.
# 
# Dependencies: lavaan, tidyverse, ggplot2, gridExtra, knitr, broom, 
#               purrr, furrr, progress, plyr
# ==============================================================================

# Load required libraries
library(lavaan)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(knitr)
library(broom)
library(purrr)
library(furrr)
library(progress)
library(plyr)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set virtual paths for GitHub
CONFIG <- list(
  data_dir = "./data",
  results_dir = "./results/multiomics_mediation",
  
  # Data file names
  static_features = "all_static_features.csv",
  icd_diagnoses = "all_icd_f_group_1000_10y_wear.csv",
  movement_features = "all_movement_features.csv",
  metabolomics_data = "merged_metabolomics_data.csv",
  blood_data = "cleaned_blood_data.csv",
  
  # Analysis parameters
  age_threshold = 60,  # Age threshold for age stratification
  min_disease_cases = 50,  # Minimum cases for analysis
  min_sample_size = 50,  # Minimum sample size per model
  significance_threshold = 0.05,  # p-value threshold for significance
  max_caf_combinations = 100,  # Maximum CAF combinations to test
  outlier_sd_threshold = 3,  # Standard deviations for outlier detection
  
  # Activity states
  activity_states = c("SB", "LPA", "MVPA", "Sleep"),
  days = c("workday", "weekend"),
  
  # Movement variables (processed)
  movement_vars = c(
    "SB_time_workday", "LPA_time_workday", "MVPA_time_workday", "Sleep_time_workday",
    "SB_time_weekend", "LPA_time_weekend", "MVPA_time_weekend", "Sleep_time_weekend"
  ),
  
  # Covariates
  covariates = c("Sex", "age", "BMI"),
  
  # Model fit thresholds
  cfi_threshold = 0.90,
  rmsea_threshold = 0.08,
  srmr_threshold = 0.08
)

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

load_data <- function(config) {
  # Load all required datasets with virtual paths
  
  cat("Loading data from:", config$data_dir, "\n")
  
  # Construct file paths
  files <- list(
    static = file.path(config$data_dir, config$static_features),
    icd = file.path(config$data_dir, config$icd_diagnoses),
    movement = file.path(config$data_dir, config$movement_features),
    metabolomics = file.path(config$data_dir, config$metabolomics_data),
    blood = file.path(config$data_dir, config$blood_data)
  )
  
  # Check if files exist
  missing_files <- sapply(files, function(f) !file.exists(f))
  if (any(missing_files)) {
    warning("Missing data files: ", 
            paste(names(files)[missing_files], collapse = ", "))
  }
  
  # Read data
  data_list <- list()
  for (name in names(files)) {
    if (file.exists(files[[name]])) {
      data_list[[name]] <- read.csv(files[[name]], stringsAsFactors = FALSE)
      cat(sprintf("  Loaded %-15s: %d rows × %d columns\n", 
                  name, nrow(data_list[[name]]), ncol(data_list[[name]])))
    } else {
      data_list[[name]] <- NULL
      cat(sprintf("  Missing: %s\n", name))
    }
  }
  
  return(data_list)
}

# ==============================================================================
# DATA PROCESSING FUNCTIONS
# ==============================================================================

process_movement_data <- function(movement_df) {
  # Process raw movement data into time proportions
  
  cat("Processing movement data...\n")
  pb <- progress_bar$new(total = nrow(movement_df), format = "[:bar] :percent :eta")
  
  results <- lapply(1:nrow(movement_df), function(i) {
    pb$tick()
    row <- movement_df[i, ]
    participant_id <- row$Participant.ID
    
    # Extract workday data (first 96 columns: 24 hours × 4 states)
    workday_data <- matrix(as.numeric(row[1:96]), nrow = 24, ncol = 4, byrow = TRUE)
    
    # Extract weekend data (next 96 columns)
    weekend_data <- matrix(as.numeric(row[97:192]), nrow = 24, ncol = 4, byrow = TRUE)
    
    # Calculate time proportions (column sums)
    workday_props <- colSums(workday_data) / 24
    weekend_props <- colSums(weekend_data) / 24
    
    # Calculate weekly weighted average (5 workdays + 2 weekend days)
    weekly_avg <- (workday_props * 5 + weekend_props * 2) / 7
    
    # Return as data frame row
    data.frame(
      Participant.ID = participant_id,
      
      # Workday proportions
      SB_time_workday = workday_props[1],
      LPA_time_workday = workday_props[2],
      MVPA_time_workday = workday_props[3],
      Sleep_time_workday = workday_props[4],
      
      # Weekend proportions
      SB_time_weekend = weekend_props[1],
      LPA_time_weekend = weekend_props[2],
      MVPA_time_weekend = weekend_props[3],
      Sleep_time_weekend = weekend_props[4],
      
      # Weekly averages
      SB_time_avg = weekly_avg[1],
      LPA_time_avg = weekly_avg[2],
      MVPA_time_avg = weekly_avg[3],
      Sleep_time_avg = weekly_avg[4],
      
      stringsAsFactors = FALSE
    )
  })
  
  movement_processed <- bind_rows(results)
  cat(sprintf("Processed movement data: %d rows × %d columns\n", 
              nrow(movement_processed), ncol(movement_processed)))
  
  return(movement_processed)
}

prepare_analysis_data <- function(movement_processed, static_df, metabolomics_df, 
                                  blood_df, icd_df, disease, config) {
  # Prepare combined dataset for analysis
  
  # Check if disease exists in ICD data
  if (!disease %in% colnames(icd_df)) {
    cat(sprintf("Warning: Disease %s not found in ICD data\n", disease))
    return(NULL)
  }
  
  # Select covariates
  covariate_cols <- c("Participant.ID", config$covariates)
  covariate_data <- static_df %>%
    select(any_of(covariate_cols))
  
  # Standardize BMI column name if it exists
  if ("Body.mass.index..BMI.._x" %in% colnames(covariate_data)) {
    covariate_data <- covariate_data %>%
      rename(BMI = Body.mass.index..BMI.._x)
  }
  
  # Process disease data - create binary variable
  disease_binary <- paste0(disease, "_binary")
  icd_processed <- icd_df %>%
    mutate(!!disease_binary := ifelse(!!sym(disease) == 1, 1, 0))
  
  # Merge all datasets
  merged_data <- movement_processed %>%
    inner_join(covariate_data, by = "Participant.ID")
  
  # Add metabolomics data if available
  if (!is.null(metabolomics_df)) {
    merged_data <- merged_data %>%
      inner_join(metabolomics_df, by = "Participant.ID")
  }
  
  # Add blood data if available
  if (!is.null(blood_df)) {
    merged_data <- merged_data %>%
      inner_join(blood_df, by = "Participant.ID")
  }
  
  # Add disease data
  merged_data <- merged_data %>%
    inner_join(icd_processed %>% 
                 select(Participant.ID, !!disease_binary), 
               by = "Participant.ID")
  
  # Process sex variable
  if ("Sex" %in% colnames(merged_data)) {
    merged_data$Sex <- ifelse(merged_data$Sex %in% c("Male", 0, "0"), 1, 0)
    # Impute missing sex with mode
    if (any(is.na(merged_data$Sex))) {
      sex_mode <- as.numeric(names(which.max(table(merged_data$Sex))))
      merged_data$Sex[is.na(merged_data$Sex)] <- sex_mode
    }
  }
  
  # Clean column names (replace dots with underscores)
  colnames(merged_data) <- gsub("\\.", "_", colnames(merged_data))
  
  # Check disease case count
  disease_count <- sum(merged_data[[disease_binary]], na.rm = TRUE)
  if (disease_count < config$min_disease_cases) {
    cat(sprintf("Warning: Insufficient cases for disease %s (%d cases)\n", 
                disease, disease_count))
    return(NULL)
  }
  
  cat(sprintf("Prepared data for %s: %d rows, %d columns, %d cases\n", 
              disease, nrow(merged_data), ncol(merged_data), disease_count))
  
  return(merged_data)
}

# ==============================================================================
# PREPROCESSING AND MODEL FITTING FUNCTIONS
# ==============================================================================

robust_preprocess <- function(data, vars, outlier_sd = NULL) {
  # Preprocess variables with outlier handling and standardization
  
  if (is.null(outlier_sd)) {
    outlier_sd <- CONFIG$outlier_sd_threshold
  }
  
  processed_data <- data
  
  for (var in vars) {
    if (var %in% colnames(processed_data)) {
      x <- processed_data[[var]]
      
      if (is.numeric(x) && length(unique(x)) > 5) {
        # Winsorize outliers
        mean_val <- mean(x, na.rm = TRUE)
        sd_val <- sd(x, na.rm = TRUE)
        lower <- mean_val - outlier_sd * sd_val
        upper <- mean_val + outlier_sd * sd_val
        
        x[x < lower] <- lower
        x[x > upper] <- upper
        
        # Standardize
        processed_data[[var]] <- scale(x)
      }
    }
  }
  
  return(processed_data)
}

safe_sem_fit <- function(model_spec, data, outcome, estimator = "WLSMV", 
                         max_iter = 1000) {
  # Safely fit SEM model with error handling
  
  tryCatch({
    # Fit model with specified estimator
    if (estimator == "WLSMV") {
      fit <- sem(model_spec, data = data, estimator = "WLSMV", 
                 ordered = outcome, meanstructure = TRUE,
                 optim.method = "NLMINB", control = list(iter.max = max_iter))
    } else {
      fit <- sem(model_spec, data = data, estimator = estimator, 
                 meanstructure = TRUE,
                 control = list(iter.max = max_iter))
    }
    
    # Check convergence
    if (lavInspect(fit, "converged")) {
      # Check for extreme coefficients
      param_est <- parameterEstimates(fit)
      extreme_coefs <- any(abs(param_est$est) > 10, na.rm = TRUE)
      
      if (extreme_coefs && estimator == "WLSMV") {
        cat("  Extreme coefficients detected, trying ML estimator\n")
        # Try ML estimator as alternative
        fit <- sem(model_spec, data = data, estimator = "ML", 
                   meanstructure = TRUE,
                   control = list(iter.max = max_iter))
      }
      
      return(fit)
    } else {
      cat("  Model did not converge\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("  Model fitting error:", e$message, "\n")
    return(NULL)
  })
}

# ==============================================================================
# COMMON ANTECEDENT FACTOR (CAF) ANALYSIS
# ==============================================================================

build_caf_model_spec <- function(indicators, movement_var, outcome, covariates) {
  # Build CAF model specification
  
  if (length(indicators) < 2) {
    stop("At least 2 indicators required for CAF model")
  }
  
  model_spec <- ""
  
  # 1. Set latent variable scale (fix first loading to 1)
  first_indicator <- indicators[1]
  other_indicators <- indicators[-1]
  
  model_spec <- paste0(model_spec, "Latent =~ 1*", first_indicator, "\n")
  if (length(other_indicators) > 0) {
    model_spec <- paste0(model_spec, 
                         "Latent ~ ", paste(other_indicators, collapse = " + "), "\n")
  }
  
  # 2. Structural model
  model_spec <- paste0(model_spec, "Latent ~ a * PA_time\n")
  model_spec <- paste0(model_spec, outcome, " ~ b * Latent\n")
  model_spec <- paste0(model_spec, outcome, " ~ c * PA_time\n")
  
  # 3. Add covariates
  if (length(covariates) > 0) {
    cov_str <- paste(covariates, collapse = " + ")
    model_spec <- paste0(model_spec, "Latent ~ ", cov_str, "\n")
    model_spec <- paste0(model_spec, outcome, " ~ ", cov_str, "\n")
  }
  
  # 4. Effect definitions
  model_spec <- paste0(model_spec, "indirect := a * b\n")
  model_spec <- paste0(model_spec, "total := c + indirect\n")
  
  return(model_spec)
}

caf_analysis_single_set <- function(data, movement_var, biomarkers, outcome, 
                                    covariates, config) {
  # Perform CAF analysis for a single set of biomarkers
  
  if (length(biomarkers) < 4) {
    cat(sprintf("  Insufficient biomarkers (%d), skipping CAF analysis\n", 
                length(biomarkers)))
    return(NULL)
  }
  
  best_model <- NULL
  best_fit_score <- -Inf
  
  # Preprocess variables
  all_continuous_vars <- c(movement_var, biomarkers, 
                          setdiff(covariates, "Sex"))
  processed_data <- robust_preprocess(data, all_continuous_vars, 
                                     config$outlier_sd_threshold)
  
  # Rename movement variable
  colnames(processed_data)[colnames(processed_data) == movement_var] <- "PA_time"
  
  # Generate combinations (limited to avoid combinatorial explosion)
  n_combinations <- min(config$max_caf_combinations, 
                        choose(length(biomarkers), 4))
  set.seed(123)  # For reproducibility
  
  cat(sprintf("  Testing %d indicator combinations\n", n_combinations))
  
  for (i in 1:n_combinations) {
    indicators <- sample(biomarkers, 4)
    
    # Build and fit model
    model_spec <- build_caf_model_spec(indicators, "PA_time", outcome, covariates)
    
    # Prepare analysis data
    analysis_vars <- c("PA_time", outcome, indicators, covariates)
    analysis_df <- processed_data %>% 
      select(all_of(analysis_vars)) %>% 
      drop_na()
    
    if (nrow(analysis_df) < config$min_sample_size) {
      next
    }
    
    # Fit model
    fit <- safe_sem_fit(model_spec, analysis_df, outcome)
    
    if (!is.null(fit) && lavInspect(fit, "converged")) {
      fit_measures <- fitMeasures(fit)
      
      # Check fit indices
      cfi <- fit_measures["cfi"]
      rmsea <- fit_measures["rmsea"]
      srmr <- fit_measures["srmr"]
      
      if (!is.na(cfi) && !is.na(rmsea) && !is.na(srmr) &&
          cfi > config$cfi_threshold && 
          rmsea < config$rmsea_threshold && 
          srmr < config$srmr_threshold) {
        
        # Check parameter estimates
        param_est <- parameterEstimates(fit)
        coefs_reasonable <- all(abs(param_est$est) < 5, na.rm = TRUE)
        
        if (coefs_reasonable) {
          # Calculate fit score
          fit_score <- cfi - (rmsea / config$rmsea_threshold) - 
            (srmr / config$srmr_threshold)
          
          if (fit_score > best_fit_score) {
            best_fit_score <- fit_score
            best_model <- list(
              movement_var = movement_var,
              indicators = indicators,
              model_spec = model_spec,
              fit = fit,
              fit_measures = fit_measures,
              n_obs = nrow(analysis_df),
              fit_score = fit_score
            )
          }
        }
      }
    }
  }
  
  if (!is.null(best_model)) {
    cat(sprintf("  Best CAF model: %d indicators, CFI=%.3f, RMSEA=%.3f\n", 
                length(best_model$indicators), 
                best_model$fit_measures["cfi"],
                best_model$fit_measures["rmsea"]))
    return(best_model)
  } else {
    cat("  No acceptable CAF model found\n")
    return(NULL)
  }
}

# ==============================================================================
# BIOMARKER SELECTION FUNCTIONS
# ==============================================================================

load_precomputed_biomarkers <- function(movement_var, disease, input_dir) {
  # Load precomputed biomarkers from Python output
  
  movement_mapping <- c(
    "SB_time_workday" = "SB_workday",
    "LPA_time_workday" = "LPA_workday",
    "MVPA_time_workday" = "MVPA_workday",
    "Sleep_time_workday" = "Sleep_workday",
    "SB_time_weekend" = "SB_weekend",
    "LPA_time_weekend" = "LPA_weekend",
    "MVPA_time_weekend" = "MVPA_weekend",
    "Sleep_time_weekend" = "Sleep_weekend"
  )
  
  r_movement_var <- movement_mapping[movement_var]
  
  # Construct file paths
  blood_file <- file.path(input_dir, paste0(r_movement_var, "_blood_top20.csv"))
  meta_file <- file.path(input_dir, paste0(r_movement_var, "_meta_top20.csv"))
  
  biomarkers <- list(blood = NULL, meta = NULL)
  
  # Load blood biomarkers
  if (file.exists(blood_file)) {
    blood_df <- read.csv(blood_file, stringsAsFactors = FALSE)
    biomarkers$blood <- blood_df$biomarker
    cat(sprintf("  Loaded %d blood biomarkers\n", length(biomarkers$blood)))
  }
  
  # Load metabolomics biomarkers
  if (file.exists(meta_file)) {
    meta_df <- read.csv(meta_file, stringsAsFactors = FALSE)
    biomarkers$meta <- meta_df$biomarker
    cat(sprintf("  Loaded %d metabolomics biomarkers\n", length(biomarkers$meta)))
  }
  
  return(biomarkers)
}

# ==============================================================================
# COMPREHENSIVE CAF ANALYSIS
# ==============================================================================

comprehensive_caf_analysis <- function(data, movement_vars, outcome, covariates, 
                                       disease, biomarkers_dir, config) {
  # Comprehensive CAF analysis for all movement variables
  
  cat(sprintf("Comprehensive CAF analysis for disease: %s\n", disease))
  cat(rep("=", 70), "\n", sep = "")
  
  best_models <- list()
  
  # Set up parallel processing
  plan(multisession, workers = min(4, availableCores() - 1))
  
  for (movement_var in movement_vars) {
    cat(sprintf("Analyzing movement variable: %s\n", movement_var))
    
    # Load precomputed biomarkers
    biomarkers <- load_precomputed_biomarkers(movement_var, disease, biomarkers_dir)
    
    if (length(biomarkers$meta) >= 4 && length(biomarkers$blood) >= 4) {
      cat(sprintf("  Metabolomics: %d, Blood: %d biomarkers\n", 
                  length(biomarkers$meta), length(biomarkers$blood)))
      
      # Perform CAF analysis in parallel
      models <- future_map(
        list(biomarkers$meta, biomarkers$blood),
        ~ caf_analysis_single_set(data, movement_var, .x, outcome, covariates, config),
        .options = furrr_options(seed = TRUE)
      )
      
      best_models[[movement_var]] <- list(
        metabolomics = models[[1]],
        blood = models[[2]]
      )
      
      # Report results
      if (!is.null(models[[1]])) {
        cat(sprintf("  Best metabolomics CAF: %d indicators, CFI=%.3f\n",
                    length(models[[1]]$indicators), 
                    models[[1]]$fit_measures["cfi"]))
      } else {
        cat("  No acceptable metabolomics CAF model\n")
      }
      
      if (!is.null(models[[2]])) {
        cat(sprintf("  Best blood CAF: %d indicators, CFI=%.3f\n",
                    length(models[[2]]$indicators), 
                    models[[2]]$fit_measures["cfi"]))
      } else {
        cat("  No acceptable blood CAF model\n")
      }
    } else {
      cat(sprintf("  Insufficient biomarkers, skipping %s\n", movement_var))
      best_models[[movement_var]] <- NULL
    }
  }
  
  # Close parallel processing
  plan(sequential)
  
  return(list(
    disease = disease,
    best_models = best_models
  ))
}

# ==============================================================================
# BIDIRECTIONAL MEDIATION ANALYSIS
# ==============================================================================

build_bidirectional_model_spec <- function(direction, meta_indicators, 
                                           blood_indicators, outcome, covariates) {
  # Build bidirectional mediation model specification
  
  if (length(meta_indicators) < 2 || length(blood_indicators) < 2) {
    stop("At least 2 indicators required for each latent variable")
  }
  
  model_spec <- ""
  
  # 1. Define latent variables
  model_spec <- paste0(model_spec, "Metabolism =~ 1*", meta_indicators[1], "\n")
  if (length(meta_indicators) > 1) {
    model_spec <- paste0(model_spec, "Metabolism ~ ", 
                        paste(meta_indicators[-1], collapse = " + "), "\n")
  }
  
  model_spec <- paste0(model_spec, "Inflammation =~ 1*", blood_indicators[1], "\n")
  if (length(blood_indicators) > 1) {
    model_spec <- paste0(model_spec, "Inflammation ~ ", 
                        paste(blood_indicators[-1], collapse = " + "), "\n")
  }
  
  # 2. Structural model based on direction
  if (direction == "meta_to_blood") {
    # Direction 1: PA → Metabolism → Inflammation → Disease
    model_spec <- paste0(model_spec, "Metabolism ~ a1 * PA_time\n")
    model_spec <- paste0(model_spec, "Inflammation ~ a2 * PA_time + d1 * Metabolism\n")
    model_spec <- paste0(model_spec, outcome, " ~ c_prime * PA_time + b1 * Metabolism + b2 * Inflammation\n")
    
    # Indirect effects
    model_spec <- paste0(model_spec, "indirect_metabolic_only := a1 * b1\n")
    model_spec <- paste0(model_spec, "indirect_inflammatory_only := a2 * b2\n")
    model_spec <- paste0(model_spec, "indirect_meta_to_blood := a1 * d1 * b2\n")
    model_spec <- paste0(model_spec, "total_indirect := indirect_metabolic_only + indirect_inflammatory_only + indirect_meta_to_blood\n")
  } else {
    # Direction 2: PA → Inflammation → Metabolism → Disease
    model_spec <- paste0(model_spec, "Inflammation ~ a1 * PA_time\n")
    model_spec <- paste0(model_spec, "Metabolism ~ a2 * PA_time + d2 * Inflammation\n")
    model_spec <- paste0(model_spec, outcome, " ~ c_prime * PA_time + b1 * Metabolism + b2 * Inflammation\n")
    
    # Indirect effects
    model_spec <- paste0(model_spec, "indirect_metabolic_only := a2 * b1\n")
    model_spec <- paste0(model_spec, "indirect_inflammatory_only := a1 * b2\n")
    model_spec <- paste0(model_spec, "indirect_blood_to_meta := a1 * d2 * b1\n")
    model_spec <- paste0(model_spec, "total_indirect := indirect_metabolic_only + indirect_inflammatory_only + indirect_blood_to_meta\n")
  }
  
  # 3. Add covariates
  if (length(covariates) > 0) {
    cov_str <- paste(covariates, collapse = " + ")
    model_spec <- paste0(model_spec, sprintf("Metabolism ~ %s\n", cov_str))
    model_spec <- paste0(model_spec, sprintf("Inflammation ~ %s\n", cov_str))
    model_spec <- paste0(model_spec, sprintf("%s ~ %s\n", outcome, cov_str))
  }
  
  model_spec <- paste0(model_spec, "total_effect := c_prime + total_indirect\n")
  
  return(model_spec)
}

bidirectional_mediation_single <- function(data, movement_var, meta_indicators, 
                                           blood_indicators, outcome, covariates, config) {
  # Perform bidirectional mediation analysis for single movement variable
  
  cat(sprintf("  Movement variable: %s\n", movement_var))
  
  # Preprocess all variables
  all_continuous_vars <- c(movement_var, meta_indicators, blood_indicators, 
                          setdiff(covariates, "Sex"))
  processed_data <- robust_preprocess(data, all_continuous_vars, 
                                     config$outlier_sd_threshold)
  
  # Rename movement variable
  colnames(processed_data)[colnames(processed_data) == movement_var] <- "PA_time"
  
  direction_results <- list()
  directions <- c("meta_to_blood", "blood_to_meta")
  
  for (direction in directions) {
    cat(sprintf("    Direction: %s\n", 
                ifelse(direction == "meta_to_blood", 
                       "Metabolism → Inflammation", 
                       "Inflammation → Metabolism")))
    
    tryCatch({
      # Build model specification
      model_spec <- build_bidirectional_model_spec(direction, meta_indicators, 
                                                   blood_indicators, outcome, covariates)
      
      # Prepare analysis data
      analysis_vars <- c("PA_time", outcome, meta_indicators, blood_indicators, covariates)
      analysis_df <- processed_data %>% 
        select(all_of(analysis_vars)) %>% 
        drop_na()
      
      if (nrow(analysis_df) < config$min_sample_size) {
        cat("    Insufficient sample size\n")
        next
      }
      
      # Fit model
      fit <- safe_sem_fit(model_spec, analysis_df, outcome)
      
      if (!is.null(fit) && lavInspect(fit, "converged")) {
        # Get model results
        fit_measures <- fitMeasures(fit)
        param_est <- parameterEstimates(fit, standardized = TRUE)
        
        # Check parameter validity
        coefs_reasonable <- all(abs(param_est$est) < 5, na.rm = TRUE)
        
        if (coefs_reasonable) {
          # Extract key paths
          extract_effect <- function(label, param_df) {
            row <- param_df[param_df$label == label, ]
            if (nrow(row) > 0) {
              list(
                estimate = row$est[1],
                std_estimate = row$std.all[1],
                p_value = row$pvalue[1]
              )
            } else {
              list(estimate = NA, std_estimate = NA, p_value = NA)
            }
          }
          
          # Extract path coefficients
          path_labels <- if (direction == "meta_to_blood") {
            c("a1", "a2", "d1", "b1", "b2", "c_prime")
          } else {
            c("a1", "a2", "d2", "b1", "b2", "c_prime")
          }
          
          path_coeffs <- lapply(path_labels, function(label) {
            extract_effect(label, param_est)
          })
          names(path_coeffs) <- path_labels
          
          # Extract indirect effects
          indirect_labels <- if (direction == "meta_to_blood") {
            c("indirect_metabolic_only", "indirect_inflammatory_only", 
              "indirect_meta_to_blood", "total_indirect")
          } else {
            c("indirect_metabolic_only", "indirect_inflammatory_only", 
              "indirect_blood_to_meta", "total_indirect")
          }
          
          indirect_effects <- lapply(indirect_labels, function(label) {
            extract_effect(label, param_est)
          })
          names(indirect_effects) <- indirect_labels
          
          # Store results
          direction_results[[direction]] <- list(
            direction = direction,
            fit = fit,
            fit_measures = fit_measures,
            param_est = param_est,
            path_coeffs = path_coeffs,
            indirect_effects = indirect_effects,
            n_obs = nrow(analysis_df)
          )
          
          cat(sprintf("    Success: CFI=%.3f, RMSEA=%.3f\n", 
                      fit_measures["cfi"], fit_measures["rmsea"]))
          
          # Report significant indirect effects
          for (effect_name in names(indirect_effects)[1:3]) {
            effect <- indirect_effects[[effect_name]]
            if (!is.na(effect$p_value) && effect$p_value < config$significance_threshold) {
              cat(sprintf("      Significant: %s = %.4f (p=%.4f)\n", 
                          effect_name, effect$estimate, effect$p_value))
            }
          }
        } else {
          cat("    Unreasonable coefficients\n")
        }
      } else {
        cat("    Model fitting failed\n")
      }
    }, error = function(e) {
      cat(sprintf("    Error: %s\n", e$message))
    })
  }
  
  return(direction_results)
}

comprehensive_bidirectional_mediation <- function(data, caf_results, outcome, covariates, config) {
  # Comprehensive bidirectional mediation analysis
  
  disease <- caf_results$disease
  cat(sprintf("\nBidirectional mediation analysis for disease: %s\n", disease))
  cat(rep("=", 70), "\n", sep = "")
  
  mediation_results <- list()
  
  for (movement_var in names(caf_results$best_models)) {
    best_models <- caf_results$best_models[[movement_var]]
    
    if (!is.null(best_models$metabolomics) && !is.null(best_models$blood)) {
      # Get best indicator combinations
      meta_indicators <- best_models$metabolomics$indicators
      blood_indicators <- best_models$blood$indicators
      
      cat(sprintf("Movement variable: %s\n", movement_var))
      cat(sprintf("  Metabolomics indicators: %d\n", length(meta_indicators)))
      cat(sprintf("  Blood indicators: %d\n", length(blood_indicators)))
      
      # Perform bidirectional mediation
      direction_results <- bidirectional_mediation_single(
        data, movement_var, meta_indicators, blood_indicators, 
        outcome, covariates, config
      )
      
      # Select best direction based on model fit
      if (length(direction_results) == 2) {
        cfi1 <- direction_results[["meta_to_blood"]]$fit_measures["cfi"]
        cfi2 <- direction_results[["blood_to_meta"]]$fit_measures["cfi"]
        
        if (cfi1 >= cfi2) {
          best_direction <- "meta_to_blood"
          best_direction_name <- "Metabolism → Inflammation"
        } else {
          best_direction <- "blood_to_meta"
          best_direction_name <- "Inflammation → Metabolism"
        }
        
        mediation_results[[movement_var]] <- list(
          movement_var = movement_var,
          meta_indicators = meta_indicators,
          blood_indicators = blood_indicators,
          best_direction = best_direction,
          best_direction_name = best_direction_name,
          all_directions = direction_results,
          best_result = direction_results[[best_direction]]
        )
        
        cat(sprintf("  Best direction: %s (CFI=%.3f)\n", 
                    best_direction_name, 
                    mediation_results[[movement_var]]$best_result$fit_measures["cfi"]))
      } else if (length(direction_results) == 1) {
        # Only one direction succeeded
        best_direction <- names(direction_results)[1]
        best_direction_name <- ifelse(best_direction == "meta_to_blood", 
                                      "Metabolism → Inflammation", 
                                      "Inflammation → Metabolism")
        
        mediation_results[[movement_var]] <- list(
          movement_var = movement_var,
          meta_indicators = meta_indicators,
          blood_indicators = blood_indicators,
          best_direction = best_direction,
          best_direction_name = best_direction_name,
          all_directions = direction_results,
          best_result = direction_results[[best_direction]]
        )
        
        cat(sprintf("  Available direction: %s (CFI=%.3f)\n", 
                    best_direction_name, 
                    mediation_results[[movement_var]]$best_result$fit_measures["cfi"]))
      } else {
        cat("  No successful bidirectional models\n")
      }
    } else {
      cat(sprintf("Missing CAF models for %s, skipping mediation\n", movement_var))
    }
  }
  
  return(list(
    disease = disease,
    mediation_results = mediation_results
  ))
}

# ==============================================================================
# RESULT EXPORT FUNCTIONS
# ==============================================================================

export_simplified_statistics <- function(caf_results, mediation_results, output_dir, config) {
  # Export simplified model statistics
  
  disease <- caf_results$disease
  cat(sprintf("\nExporting simplified statistics for disease: %s\n", disease))
  
  all_results <- list()
  
  # Extract CAF model statistics
  for (movement_var in names(caf_results$best_models)) {
    best_models <- caf_results$best_models[[movement_var]]
    
    # Metabolomics CAF
    if (!is.null(best_models$metabolomics)) {
      model <- best_models$metabolomics
      fit_measures <- model$fit_measures
      param_est <- parameterEstimates(model$fit, standardized = TRUE)
      
      # Extract effects
      direct_effect <- param_est[param_est$label == "c", ]
      indirect_effect <- param_est[param_est$label == "indirect", ]
      total_effect <- param_est[param_est$label == "total", ]
      
      all_results[[paste0(movement_var, "_metabolomics_CAF")]] <- data.frame(
        disease = disease,
        movement_var = movement_var,
        model_type = "CAF_Metabolomics",
        n_indicators = length(model$indicators),
        CFI = round(fit_measures["cfi"], 3),
        RMSEA = round(fit_measures["rmsea"], 3),
        SRMR = round(fit_measures["srmr"], 3),
        direct_effect = ifelse(nrow(direct_effect) > 0, 
                               round(direct_effect$est[1], 4), NA),
        direct_effect_p = ifelse(nrow(direct_effect) > 0, 
                                 round(direct_effect$pvalue[1], 6), NA),
        indirect_effect = ifelse(nrow(indirect_effect) > 0, 
                                 round(indirect_effect$est[1], 4), NA),
        indirect_effect_p = ifelse(nrow(indirect_effect) > 0, 
                                   round(indirect_effect$pvalue[1], 6), NA),
        total_effect = ifelse(nrow(total_effect) > 0, 
                              round(total_effect$est[1], 4), NA),
        total_effect_p = ifelse(nrow(total_effect) > 0, 
                                round(total_effect$pvalue[1], 6), NA),
        stringsAsFactors = FALSE
      )
    }
    
    # Blood CAF
    if (!is.null(best_models$blood)) {
      model <- best_models$blood
      fit_measures <- model$fit_measures
      param_est <- parameterEstimates(model$fit, standardized = TRUE)
      
      # Extract effects
      direct_effect <- param_est[param_est$label == "c", ]
      indirect_effect <- param_est[param_est$label == "indirect", ]
      total_effect <- param_est[param_est$label == "total", ]
      
      all_results[[paste0(movement_var, "_blood_CAF")]] <- data.frame(
        disease = disease,
        movement_var = movement_var,
        model_type = "CAF_Blood",
        n_indicators = length(model$indicators),
        CFI = round(fit_measures["cfi"], 3),
        RMSEA = round(fit_measures["rmsea"], 3),
        SRMR = round(fit_measures["srmr"], 3),
        direct_effect = ifelse(nrow(direct_effect) > 0, 
                               round(direct_effect$est[1], 4), NA),
        direct_effect_p = ifelse(nrow(direct_effect) > 0, 
                                 round(direct_effect$pvalue[1], 6), NA),
        indirect_effect = ifelse(nrow(indirect_effect) > 0, 
                                 round(indirect_effect$est[1], 4), NA),
        indirect_effect_p = ifelse(nrow(indirect_effect) > 0, 
                                   round(indirect_effect$pvalue[1], 6), NA),
        total_effect = ifelse(nrow(total_effect) > 0, 
                              round(total_effect$est[1], 4), NA),
        total_effect_p = ifelse(nrow(total_effect) > 0, 
                                round(total_effect$pvalue[1], 6), NA),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Extract mediation model statistics
  for (movement_var in names(mediation_results$mediation_results)) {
    mediation_result <- mediation_results$mediation_results[[movement_var]]
    best_result <- mediation_result$best_result
    
    if (!is.null(best_result)) {
      fit_measures <- best_result$fit_measures
      param_est <- best_result$param_est
      
      # Extract key effects
      direct_effect <- param_est[param_est$label == "c_prime", ]
      total_indirect <- param_est[param_est$label == "total_indirect", ]
      total_effect <- param_est[param_est$label == "total_effect", ]
      
      # Extract specific indirect effects
      indirect_metabolic <- param_est[param_est$label == "indirect_metabolic_only", ]
      indirect_inflammatory <- param_est[param_est$label == "indirect_inflammatory_only", ]
      
      if (mediation_result$best_direction == "meta_to_blood") {
        indirect_chain <- param_est[param_est$label == "indirect_meta_to_blood", ]
      } else {
        indirect_chain <- param_est[param_est$label == "indirect_blood_to_meta", ]
      }
      
      all_results[[paste0(movement_var, "_mediation")]] <- data.frame(
        disease = disease,
        movement_var = movement_var,
        model_type = paste0("Mediation_", mediation_result$best_direction),
        CFI = round(fit_measures["cfi"], 3),
        RMSEA = round(fit_measures["rmsea"], 3),
        SRMR = round(fit_measures["srmr"], 3),
        best_direction = mediation_result$best_direction,
        direct_effect = ifelse(nrow(direct_effect) > 0, 
                               round(direct_effect$est[1], 4), NA),
        direct_effect_p = ifelse(nrow(direct_effect) > 0, 
                                 round(direct_effect$pvalue[1], 6), NA),
        total_indirect_effect = ifelse(nrow(total_indirect) > 0, 
                                       round(total_indirect$est[1], 4), NA),
        total_indirect_p = ifelse(nrow(total_indirect) > 0, 
                                  round(total_indirect$pvalue[1], 6), NA),
        total_effect = ifelse(nrow(total_effect) > 0, 
                              round(total_effect$est[1], 4), NA),
        total_effect_p = ifelse(nrow(total_effect) > 0, 
                                round(total_effect$pvalue[1], 6), NA),
        indirect_metabolic_effect = ifelse(nrow(indirect_metabolic) > 0, 
                                           round(indirect_metabolic$est[1], 4), NA),
        indirect_metabolic_p = ifelse(nrow(indirect_metabolic) > 0, 
                                      round(indirect_metabolic$pvalue[1], 6), NA),
        indirect_inflammatory_effect = ifelse(nrow(indirect_inflammatory) > 0, 
                                              round(indirect_inflammatory$est[1], 4), NA),
        indirect_inflammatory_p = ifelse(nrow(indirect_inflammatory) > 0, 
                                         round(indirect_inflammatory$pvalue[1], 6), NA),
        indirect_chain_effect = ifelse(nrow(indirect_chain) > 0, 
                                       round(indirect_chain$est[1], 4), NA),
        indirect_chain_p = ifelse(nrow(indirect_chain) > 0, 
                                  round(indirect_chain$pvalue[1], 6), NA),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all results
  if (length(all_results) > 0) {
    stats_df <- do.call(plyr::rbind.fill, all_results)
    
    # Save to CSV
    csv_file <- file.path(output_dir, "simplified_model_statistics.csv")
    write.csv(stats_df, csv_file, row.names = FALSE)
    cat(sprintf("Simplified statistics saved to: %s\n", csv_file))
    cat(sprintf("Contains %d model statistics\n", nrow(stats_df)))
    
    return(stats_df)
  } else {
    cat("No model statistics to export\n")
    return(NULL)
  }
}

export_complete_path_coefficients <- function(caf_results, mediation_results, output_dir) {
  # Export complete path coefficients
  
  disease <- caf_results$disease
  cat(sprintf("\nExporting complete path coefficients for disease: %s\n", disease))
  
  all_paths <- list()
  
  # Extract CAF model path coefficients
  for (movement_var in names(caf_results$best_models)) {
    best_models <- caf_results$best_models[[movement_var]]
    
    # Metabolomics CAF
    if (!is.null(best_models$metabolomics)) {
      param_est <- parameterEstimates(best_models$metabolomics$fit, standardized = TRUE)
      param_est$disease <- disease
      param_est$movement_var <- movement_var
      param_est$model_type <- "CAF_Metabolomics"
      all_paths[[paste0(movement_var, "_metabolomics_CAF")]] <- param_est
    }
    
    # Blood CAF
    if (!is.null(best_models$blood)) {
      param_est <- parameterEstimates(best_models$blood$fit, standardized = TRUE)
      param_est$disease <- disease
      param_est$movement_var <- movement_var
      param_est$model_type <- "CAF_Blood"
      all_paths[[paste0(movement_var, "_blood_CAF")]] <- param_est
    }
  }
  
  # Extract mediation model path coefficients
  for (movement_var in names(mediation_results$mediation_results)) {
    mediation_result <- mediation_results$mediation_results[[movement_var]]
    
    for (direction_name in names(mediation_result$all_directions)) {
      direction_result <- mediation_result$all_directions[[direction_name]]
      
      if (!is.null(direction_result)) {
        param_est <- parameterEstimates(direction_result$fit, standardized = TRUE)
        param_est$disease <- disease
        param_est$movement_var <- movement_var
        param_est$model_type <- paste0("Mediation_", direction_name)
        all_paths[[paste0(movement_var, "_mediation_", direction_name)]] <- param_est
      }
    }
  }
  
  # Combine all path coefficients
  if (length(all_paths) > 0) {
    paths_df <- do.call(plyr::rbind.fill, all_paths)
    
    # Save complete path coefficients
    complete_file <- file.path(output_dir, "complete_path_coefficients.csv")
    write.csv(paths_df, complete_file, row.names = FALSE)
    cat(sprintf("Complete path coefficients saved to: %s\n", complete_file))
    
    # Save significant path coefficients
    significant_paths <- paths_df %>%
      filter(pvalue < CONFIG$significance_threshold) %>%
      arrange(disease, movement_var, model_type, pvalue)
    
    if (nrow(significant_paths) > 0) {
      sig_file <- file.path(output_dir, "significant_path_coefficients.csv")
      write.csv(significant_paths, sig_file, row.names = FALSE)
      cat(sprintf("Significant path coefficients saved to: %s\n", sig_file))
      cat(sprintf("Contains %d significant paths (p < %.3f)\n", 
                  nrow(significant_paths), CONFIG$significance_threshold))
    } else {
      cat("No significant path coefficients found\n")
    }
    
    return(list(
      complete = paths_df,
      significant = significant_paths
    ))
  } else {
    cat("No path coefficients to export\n")
    return(NULL)
  }
}

# ==============================================================================
# MAIN ANALYSIS PIPELINE
# ==============================================================================

analyze_single_disease <- function(disease, config, biomarkers_dir = "./biomarkers") {
  # Main pipeline for single disease analysis
  
  cat(sprintf("\n" + rep("=", 70) + "\n"))
  cat(sprintf("ANALYZING DISEASE: %s\n", disease))
  cat(rep("=", 70), "\n", sep = "")
  
  # Create output directory
  disease_dir <- file.path(config$results_dir, disease)
  if (!dir.exists(disease_dir)) {
    dir.create(disease_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", disease_dir))
  }
  
  # Load data
  cat("\n1. Loading data...\n")
  data_list <- load_data(config)
  
  # Check required data
  required <- c("movement", "static", "icd")
  missing <- setdiff(required, names(data_list))
  if (length(missing) > 0) {
    cat(sprintf("Missing required data: %s\n", paste(missing, collapse = ", ")))
    return(NULL)
  }
  
  # Process movement data
  cat("\n2. Processing movement data...\n")
  movement_processed <- process_movement_data(data_list$movement)
  
  # Prepare analysis data
  cat("\n3. Preparing analysis data...\n")
  outcome <- paste0(disease, "_binary")
  
  analysis_data <- prepare_analysis_data(
    movement_processed = movement_processed,
    static_df = data_list$static,
    metabolomics_df = data_list$metabolomics,
    blood_df = data_list$blood,
    icd_df = data_list$icd,
    disease = disease,
    config = config
  )
  
  if (is.null(analysis_data)) {
    return(NULL)
  }
  
  # Perform CAF analysis
  cat("\n4. Performing CAF analysis...\n")
  caf_results <- comprehensive_caf_analysis(
    data = analysis_data,
    movement_vars = config$movement_vars,
    outcome = outcome,
    covariates = config$covariates,
    disease = disease,
    biomarkers_dir = biomarkers_dir,
    config = config
  )
  
  # Perform bidirectional mediation analysis
  cat("\n5. Performing bidirectional mediation analysis...\n")
  mediation_results <- comprehensive_bidirectional_mediation(
    data = analysis_data,
    caf_results = caf_results,
    outcome = outcome,
    covariates = config$covariates,
    config = config
  )
  
  # Export results
  cat("\n6. Exporting results...\n")
  simplified_stats <- export_simplified_statistics(
    caf_results = caf_results,
    mediation_results = mediation_results,
    output_dir = disease_dir,
    config = config
  )
  
  path_coeffs <- export_complete_path_coefficients(
    caf_results = caf_results,
    mediation_results = mediation_results,
    output_dir = disease_dir
  )
  
  # Return summary
  list(
    disease = disease,
    output_dir = disease_dir,
    n_CAF_models = length(caf_results$best_models),
    n_mediation_models = length(mediation_results$mediation_results),
    simplified_stats = simplified_stats,
    path_coeffs = path_coeffs
  )
}

main_multiomics_analysis <- function(diseases, config, biomarkers_dir = "./biomarkers") {
  # Main pipeline for multiple diseases
  
  cat(rep("=", 70), "\n", sep = "")
  cat("MULTI-OMICS MEDIATION ANALYSIS PIPELINE\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # Set memory limit
  options(future.globals.maxSize = 2 * 1024 * 1024 * 1024)  # 2GB
  
  # Create results directory
  if (!dir.exists(config$results_dir)) {
    dir.create(config$results_dir, recursive = TRUE)
  }
  
  # Analyze each disease
  all_results <- list()
  
  for (disease in diseases) {
    cat(sprintf("\nProcessing disease: %s\n", disease))
    
    # Clean memory
    gc()
    
    # Analyze single disease
    result <- analyze_single_disease(
      disease = disease,
      config = config,
      biomarkers_dir = biomarkers_dir
    )
    
    if (!is.null(result)) {
      all_results[[disease]] <- result
    }
    
    # Clean memory
    gc()
  }
  
  # Print summary
  cat(rep("=", 70), "\n", sep = "")
  cat("ANALYSIS SUMMARY\n")
  cat(rep("=", 70), "\n", sep = "")
  
  for (disease in names(all_results)) {
    result <- all_results[[disease]]
    cat(sprintf("  %-10s: %d CAF models, %d mediation models\n",
                disease, result$n_CAF_models, result$n_mediation_models))
    cat(sprintf("            Output: %s\n", result$output_dir))
  }
  
  cat(sprintf("\nTotal diseases analyzed: %d\n", length(all_results)))
  cat(sprintf("Results saved to: %s\n", config$results_dir))
  
  return(all_results)
}

# ==============================================================================
# EXAMPLE EXECUTION
# ==============================================================================

if (FALSE) {  # Set to TRUE to run example
  # Example diseases (ICD codes)
  example_diseases <- c("K74")  # Liver fibrosis and cirrhosis
  
  # Optional: add more diseases
  # example_diseases <- c("N18", "F32", "E11", "J44", "K74")
  
  # Update configuration if needed
  CONFIG$data_dir <- "./example_data"  # Update with your data directory
  CONFIG$results_dir <- "./example_results"
  
  # Run analysis
  results <- main_multiomics_analysis(
    diseases = example_diseases,
    config = CONFIG,
    biomarkers_dir = "./example_biomarkers"  # Directory with precomputed biomarkers
  )
  
  # Print session info
  sessionInfo()
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

explore_data_structure <- function(config) {
  # Explore data structure and dimensions
  
  cat("Exploring data structure...\n")
  
  data_list <- load_data(config)
  
  cat("\nData summary:\n")
  for (name in names(data_list)) {
    if (!is.null(data_list[[name]])) {
      cat(sprintf("  %-15s: %d rows × %d columns\n", 
                  name, nrow(data_list[[name]]), ncol(data_list[[name]])))
      
      # Show first few column names for large datasets
      if (ncol(data_list[[name]]) > 10) {
        cat(sprintf("    Columns: %s, ...\n", 
                    paste(colnames(data_list[[name]])[1:5], collapse = ", ")))
      }
    }
  }
  
  # Check participant overlap
  cat("\nParticipant overlap:\n")
  ids <- list()
  for (name in names(data_list)) {
    if (!is.null(data_list[[name]]) && "Participant.ID" %in% colnames(data_list[[name]])) {
      ids[[name]] <- data_list[[name]]$Participant.ID
    }
  }
  
  if (length(ids) > 0) {
    common_ids <- Reduce(intersect, ids)
    cat(sprintf("  Common participants: %d\n", length(common_ids)))
  }
  
  return(data_list)
}

print_analysis_parameters <- function(config) {
  # Print analysis parameters
  
  cat("Analysis Parameters:\n")
  cat(rep("-", 50), "\n", sep = "")
  
  params <- c(
    "min_disease_cases", "min_sample_size", "significance_threshold",
    "cfi_threshold", "rmsea_threshold", "srmr_threshold"
  )
  
  for (param in params) {
    cat(sprintf("  %-25s: %s\n", param, config[[param]]))
  }
  
  cat("\nMovement variables to analyze:\n")
  for (var in config$movement_vars) {
    cat(sprintf("  - %s\n", var))
  }
}

# ==============================================================================
# SESSION INFO
# ==============================================================================

print_session_info <- function() {
  # Print session information
  
  cat("\n" + rep("=", 70) + "\n")
  cat("SESSION INFORMATION\n")
  cat(rep("=", 70) + "\n")
  
  info <- sessionInfo()
  cat("R version:", R.version.string, "\n")
  cat("Platform:", R.version$platform, "\n\n")
  
  cat("Loaded packages:\n")
  packages <- c("lavaan", "tidyverse", "ggplot2", "gridExtra", "knitr", 
                "broom", "purrr", "furrr", "progress", "plyr")
  
  for (pkg in packages) {
    if (pkg %in% names(info$otherPkgs)) {
      version <- info$otherPkgs[[pkg]]$Version
      cat(sprintf("  %-15s v%s\n", pkg, version))
    }
  }
}

# ==============================================================================
# MAIN EXECUTION (COMMENTED OUT FOR SAFETY)
# ==============================================================================

# Uncomment and modify to run the analysis
# main_multiomics_analysis(
#   diseases = c("K74"),
#   config = CONFIG,
#   biomarkers_dir = "./biomarkers"
# )