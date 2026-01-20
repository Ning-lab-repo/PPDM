# ==============================================================================
# Protein-Mediation Analysis for Physical Activity and Disease Risk
# 
# Description: This script performs structural equation modeling (SEM) to 
# investigate the mediation effects of proteins in the relationship between 
# physical activity and disease risk. It analyzes multiple diseases, activity 
# types, and proteins in parallel.
# 
# Dependencies: lavaan, tidyverse, progress, furrr, purrr
# ==============================================================================

# Load required libraries
library(lavaan)
library(tidyverse)
library(progress)
library(furrr)
library(purrr)

# Set memory limits for large datasets
options(future.globals.maxSize = 5 * 1024 * 1024 * 1024)  # 5GB

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set virtual paths for GitHub
CONFIG <- list(
  data_dir = "./data",
  results_dir = "./results/protein_mediation",
  
  # Data file names
  static_features = "all_static_features.csv",
  icd_diagnoses = "all_icd_f_group_1000_10y_wear.csv",
  movement_features = "all_movement_features.csv",
  protein_data = "merged_data.csv",
  
  # Analysis parameters
  age_threshold = 60,  # Age threshold for age stratification
  min_disease_cases = 50,  # Minimum cases for analysis
  min_sample_size = 50,  # Minimum sample size per model
  indirect_effect_threshold = 0.05,  # p-value threshold for significance
  max_proteins_per_disease = NULL,  # NULL = analyze all proteins
  
  # Activity states
  activity_states = c("SB", "LPA", "MVPA", "Sleep"),
  days = c("workday", "weekend"),
  
  # Covariates
  covariates = c("Sex", "age", "BMI")
)

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

load_data <- function(config) {
  # Load all required datasets with virtual paths
  
  cat("Loading data from:", config$data_dir, "\n")
  
  # Construct file paths
  files <- list(
    imm = file.path(config$data_dir, config$static_features),
    icd = file.path(config$data_dir, config$icd_diagnoses),
    movement = file.path(config$data_dir, config$movement_features),
    prot = file.path(config$data_dir, config$protein_data)
  )
  
  # Check if files exist
  missing_files <- sapply(files, function(f) !file.exists(f))
  if (any(missing_files)) {
    stop("Missing data files: ", 
         paste(names(files)[missing_files], collapse = ", "))
  }
  
  # Read data
  data_list <- list(
    imm = read.csv(files$imm, stringsAsFactors = FALSE),
    icd = read.csv(files$icd, stringsAsFactors = FALSE),
    movement = read.csv(files$movement, stringsAsFactors = FALSE),
    prot = read.csv(files$prot, stringsAsFactors = FALSE)
  )
  
  # Print data dimensions
  cat("\nData dimensions:\n")
  for (name in names(data_list)) {
    cat(sprintf("  %-15s: %d rows × %d columns\n", 
                name, nrow(data_list[[name]]), ncol(data_list[[name]])))
  }
  
  return(data_list)
}

# ==============================================================================
# DATA PROCESSING FUNCTIONS
# ==============================================================================

process_movement_data <- function(movement_df) {
  # Process raw movement data into time proportions
  
  results <- list()
  total_rows <- nrow(movement_df)
  
  cat("Processing movement data for", total_rows, "participants...\n")
  pb <- progress_bar$new(total = total_rows, format = "[:bar] :percent :eta")
  
  for (i in 1:total_rows) {
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
    
    # Store results
    result <- list(
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
      Sleep_time_avg = weekly_avg[4]
    )
    
    results[[i]] <- result
  }
  
  movement_processed <- bind_rows(results)
  cat(sprintf("\nProcessed movement data: %d rows × %d columns\n", 
              nrow(movement_processed), ncol(movement_processed)))
  
  return(movement_processed)
}

prepare_analysis_data <- function(movement_processed, imm_df, prot_df, 
                                  icd_df, disease, config) {
  # Prepare combined dataset for analysis
  
  # Select covariates
  covariate_cols <- c("Participant.ID", config$covariates)
  covariate_data <- imm_df %>%
    select(any_of(covariate_cols))
  
  # Standardize BMI column name if it exists
  if ("Body.mass.index..BMI.._x" %in% colnames(covariate_data)) {
    covariate_data <- covariate_data %>%
      rename(BMI = Body.mass.index..BMI.._x)
  }
  
  # Process disease data - create binary variable
  if (!disease %in% colnames(icd_df)) {
    cat("Warning: Disease", disease, "not found in ICD data\n")
    return(NULL)
  }
  
  disease_binary <- paste0(disease, "_binary")
  icd_processed <- icd_df %>%
    mutate(!!disease_binary := ifelse(!!sym(disease) == 1, 1, 0))
  
  # Merge all datasets
  merged_data <- movement_processed %>%
    inner_join(covariate_data, by = "Participant.ID") %>%
    inner_join(prot_df, by = "Participant.ID") %>%
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
    cat("Warning: Insufficient cases for disease", disease, 
        "(", disease_count, "cases)\n")
    return(NULL)
  }
  
  cat(sprintf("Prepared data for %s: %d rows, %d columns, %d cases\n", 
              disease, nrow(merged_data), ncol(merged_data), disease_count))
  
  return(merged_data)
}

# ==============================================================================
# MEDIATION ANALYSIS FUNCTIONS
# ==============================================================================

single_protein_mediation <- function(data, movement_var, protein_var, 
                                     outcome, covariates, config) {
  # Perform single protein mediation analysis
  
  tryCatch({
    # Build mediation model
    model_spec <- ""
    
    # Structural paths
    model_spec <- paste0(model_spec, sprintf("%s ~ a * PA_time\n", protein_var))
    model_spec <- paste0(model_spec, 
                         sprintf("%s ~ b * %s + c * PA_time\n", 
                                 outcome, protein_var))
    
    # Add covariates
    if (length(covariates) > 0) {
      cov_str <- paste(covariates, collapse = " + ")
      model_spec <- paste0(model_spec, 
                           sprintf("%s ~ %s\n", protein_var, cov_str))
      model_spec <- paste0(model_spec, 
                           sprintf("%s ~ %s\n", outcome, cov_str))
    }
    
    # Define effects
    model_spec <- paste0(model_spec, "indirect := a * b\n")
    model_spec <- paste0(model_spec, "total := c + indirect\n")
    
    # Prepare analysis dataset
    analysis_vars <- c(movement_var, outcome, protein_var, covariates)
    analysis_df <- data %>%
      select(all_of(analysis_vars)) %>%
      drop_na()
    
    # Rename movement variable for model
    colnames(analysis_df)[colnames(analysis_df) == movement_var] <- "PA_time"
    
    # Standardize continuous variables
    continuous_vars <- c("PA_time", protein_var)
    for (var in continuous_vars) {
      if (var %in% colnames(analysis_df) && 
          length(unique(analysis_df[[var]])) > 5) {
        analysis_df[[var]] <- scale(analysis_df[[var]])
      }
    }
    
    # Check sample size
    if (nrow(analysis_df) < config$min_sample_size) {
      return(NULL)
    }
    
    # Fit SEM model
    fit <- sem(model_spec, data = analysis_df, 
               estimator = "WLSMV", 
               ordered = outcome, 
               meanstructure = TRUE)
    
    # Extract model results
    fit_measures <- fitMeasures(fit)
    param_est <- parameterEstimates(fit, standardized = TRUE)
    
    # Extract path coefficients
    extract_path <- function(label, param_df) {
      path_row <- param_df[param_df$label == label, ]
      if (nrow(path_row) > 0) {
        list(
          estimate = path_row$est[1],
          std_estimate = path_row$std.all[1],
          pvalue = path_row$pvalue[1]
        )
      } else {
        list(estimate = NA, std_estimate = NA, pvalue = NA)
      }
    }
    
    # Extract all paths
    a_path <- extract_path("a", param_est)
    b_path <- extract_path("b", param_est)
    c_path <- extract_path("c", param_est)
    indirect_path <- extract_path("indirect", param_est)
    total_path <- extract_path("total", param_est)
    
    # Return results
    list(
      movement_var = movement_var,
      protein = protein_var,
      outcome = outcome,
      n_obs = nrow(analysis_df),
      
      # Path estimates
      a_estimate = a_path$estimate,
      a_std_estimate = a_path$std_estimate,
      a_pvalue = a_path$pvalue,
      
      b_estimate = b_path$estimate,
      b_std_estimate = b_path$std_estimate,
      b_pvalue = b_path$pvalue,
      
      c_estimate = c_path$estimate,
      c_std_estimate = c_path$std_estimate,
      c_pvalue = c_path$pvalue,
      
      indirect_estimate = indirect_path$estimate,
      indirect_std_estimate = indirect_path$std_estimate,
      indirect_pvalue = indirect_path$pvalue,
      
      total_estimate = total_path$estimate,
      total_std_estimate = total_path$std_estimate,
      total_pvalue = total_path$pvalue,
      
      # Model fit indices
      CFI = fit_measures["cfi"],
      RMSEA = fit_measures["rmsea"],
      SRMR = fit_measures["srmr"],
      
      # Full model objects (optional)
      # fit = fit,
      # fit_measures = fit_measures,
      # param_est = param_est
    )
    
  }, error = function(e) {
    # Return NULL for failed models
    return(NULL)
  })
}

analyze_movement_variable <- function(movement_var, data, outcome, 
                                      covariates, protein_vars, config) {
  # Analyze all proteins for a single movement variable
  
  cat("  Analyzing movement variable:", movement_var, "\n")
  
  # Process proteins in parallel
  results <- future_map(protein_vars, function(protein_var) {
    single_protein_mediation(
      data = data,
      movement_var = movement_var,
      protein_var = protein_var,
      outcome = outcome,
      covariates = covariates,
      config = config
    )
  }, .options = furrr_options(seed = TRUE))
  
  # Name results by protein
  names(results) <- protein_vars
  
  # Remove NULL results
  results <- results[!sapply(results, is.null)]
  
  cat(sprintf("    Completed: %d/%d protein models\n", 
              length(results), length(protein_vars)))
  
  return(results)
}

comprehensive_mediation_analysis <- function(data, movement_vars, outcome, 
                                             covariates, protein_vars, config) {
  # Comprehensive mediation analysis for all movement variables and proteins
  
  cat("Starting comprehensive protein mediation analysis\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # Limit number of proteins if specified
  if (!is.null(config$max_proteins_per_disease) && 
      length(protein_vars) > config$max_proteins_per_disease) {
    set.seed(123)
    protein_vars <- sample(protein_vars, config$max_proteins_per_disease)
    cat(sprintf("Randomly selected %d proteins for analysis\n", 
                config$max_proteins_per_disease))
  }
  
  total_models <- length(movement_vars) * length(protein_vars)
  cat(sprintf("Total models to run: %d (%d movement vars × %d proteins)\n", 
              total_models, length(movement_vars), length(protein_vars)))
  
  # Analyze all movement variables in parallel
  all_results <- future_map(movement_vars, function(movement_var) {
    analyze_movement_variable(
      movement_var = movement_var,
      data = data,
      outcome = outcome,
      covariates = covariates,
      protein_vars = protein_vars,
      config = config
    )
  }, .options = furrr_options(seed = TRUE))
  
  # Name results by movement variable
  names(all_results) <- movement_vars
  
  cat(sprintf("\nAnalysis complete! Analyzed %d movement variables × %d proteins = %d models\n",
              length(movement_vars), length(protein_vars), total_models))
  
  return(all_results)
}

# ==============================================================================
# RESULT PROCESSING FUNCTIONS
# ==============================================================================

save_mediation_results <- function(mediation_results, output_dir, disease, config) {
  # Save mediation results to CSV files
  
  cat(sprintf("\nSaving results for disease: %s\n", disease))
  
  all_results <- list()
  
  # Extract results from nested structure
  for (movement_var in names(mediation_results)) {
    movement_results <- mediation_results[[movement_var]]
    
    for (protein_var in names(movement_results)) {
      result <- movement_results[[protein_var]]
      
      if (!is.null(result)) {
        result_row <- data.frame(
          disease = disease,
          movement_variable = movement_var,
          protein = protein_var,
          n_obs = result$n_obs,
          
          # Path a: activity -> protein
          a_estimate = round(result$a_estimate, 6),
          a_std_estimate = round(result$a_std_estimate, 6),
          a_pvalue = round(result$a_pvalue, 6),
          
          # Path b: protein -> disease
          b_estimate = round(result$b_estimate, 6),
          b_std_estimate = round(result$b_std_estimate, 6),
          b_pvalue = round(result$b_pvalue, 6),
          
          # Path c: activity -> disease (direct)
          c_estimate = round(result$c_estimate, 6),
          c_std_estimate = round(result$c_std_estimate, 6),
          c_pvalue = round(result$c_pvalue, 6),
          
          # Indirect effect
          indirect_estimate = round(result$indirect_estimate, 6),
          indirect_std_estimate = round(result$indirect_std_estimate, 6),
          indirect_pvalue = round(result$indirect_pvalue, 6),
          
          # Total effect
          total_estimate = round(result$total_estimate, 6),
          total_std_estimate = round(result$total_std_estimate, 6),
          total_pvalue = round(result$total_pvalue, 6),
          
          # Model fit
          CFI = round(result$CFI, 4),
          RMSEA = round(result$RMSEA, 4),
          SRMR = round(result$SRMR, 4),
          
          stringsAsFactors = FALSE
        )
        
        all_results[[paste0(movement_var, "_", protein_var)]] <- result_row
      }
    }
  }
  
  # Combine all results
  if (length(all_results) > 0) {
    results_df <- bind_rows(all_results)
    
    # Save complete results
    complete_file <- file.path(output_dir, "complete_mediation_results.csv")
    write.csv(results_df, complete_file, row.names = FALSE)
    cat(sprintf("Complete results saved to: %s\n", complete_file))
    cat(sprintf("Contains %d protein mediation models\n", nrow(results_df)))
    
    # Filter significant results
    sig_threshold <- config$indirect_effect_threshold
    significant_results <- results_df %>%
      filter(indirect_pvalue < sig_threshold) %>%
      arrange(indirect_pvalue)
    
    if (nrow(significant_results) > 0) {
      sig_file <- file.path(output_dir, "significant_mediation_results.csv")
      write.csv(significant_results, sig_file, row.names = FALSE)
      cat(sprintf("Significant results saved to: %s\n", sig_file))
      cat(sprintf("Found %d significant mediation effects (p < %.3f)\n", 
                  nrow(significant_results), sig_threshold))
      
      # Display top results
      cat("\nTop 10 most significant mediation effects:\n")
      top10 <- head(significant_results, 10)
      for (i in 1:nrow(top10)) {
        row <- top10[i, ]
        cat(sprintf("  %2d. %-15s | %-20s: indirect=%.4f, p=%.6f\n", 
                    i, row$movement_variable, row$protein, 
                    row$indirect_estimate, row$indirect_pvalue))
      }
    } else {
      cat(sprintf("No significant mediation effects found (p < %.3f)\n", sig_threshold))
    }
    
    return(list(
      complete_results = results_df,
      significant_results = significant_results
    ))
  } else {
    cat("Warning: No valid mediation results found\n")
    return(NULL)
  }
}

analyze_single_disease <- function(disease, movement_vars, covariates, 
                                   data_list, config) {
  # Main function for single disease analysis
  
  cat(sprintf("\n>>> Analyzing disease: %s\n", disease))
  
  # Create output directory
  disease_dir <- file.path(config$results_dir, disease)
  if (!dir.exists(disease_dir)) {
    dir.create(disease_dir, recursive = TRUE)
  }
  
  # Prepare outcome variable name
  outcome <- paste0(disease, "_binary")
  
  # Prepare analysis data
  analysis_data <- prepare_analysis_data(
    movement_processed = data_list$movement_processed,
    imm_df = data_list$imm,
    prot_df = data_list$prot,
    icd_df = data_list$icd,
    disease = disease,
    config = config
  )
  
  if (is.null(analysis_data)) {
    return(NULL)
  }
  
  # Identify protein variables
  non_protein_vars <- c("Participant_ID", covariates, outcome, movement_vars)
  all_vars <- colnames(analysis_data)
  protein_vars <- setdiff(all_vars, non_protein_vars)
  
  # Filter for numeric variables with sufficient variation
  numeric_protein_vars <- c()
  for (var in protein_vars) {
    var_data <- analysis_data[[var]]
    if (is.numeric(var_data) && 
        length(unique(na.omit(var_data))) > 10 &&
        sum(!is.na(var_data)) > config$min_sample_size) {
      numeric_protein_vars <- c(numeric_protein_vars, var)
    }
  }
  
  cat(sprintf("Identified %d protein variables for analysis\n", 
              length(numeric_protein_vars)))
  
  if (length(numeric_protein_vars) == 0) {
    cat("Warning: No valid protein variables found\n")
    return(NULL)
  }
  
  # Run comprehensive mediation analysis
  mediation_results <- comprehensive_mediation_analysis(
    data = analysis_data,
    movement_vars = movement_vars,
    outcome = outcome,
    covariates = covariates,
    protein_vars = numeric_protein_vars,
    config = config
  )
  
  # Save results
  csv_results <- save_mediation_results(
    mediation_results = mediation_results,
    output_dir = disease_dir,
    disease = disease,
    config = config
  )
  
  # Return summary
  list(
    disease = disease,
    output_dir = disease_dir,
    n_models = if (!is.null(csv_results$complete_results)) 
      nrow(csv_results$complete_results) else 0,
    n_significant = if (!is.null(csv_results$significant_results)) 
      nrow(csv_results$significant_results) else 0
  )
}

# ==============================================================================
# MAIN ANALYSIS PIPELINE
# ==============================================================================

main_protein_mediation_analysis <- function(diseases, config, n_workers = 4) {
  # Main pipeline for protein mediation analysis
  
  cat("PROTEIN MEDIATION ANALYSIS PIPELINE\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # Set up parallel processing
  plan(multisession, workers = n_workers)
  cat(sprintf("Using %d parallel workers\n", n_workers))
  
  # Load and process data
  cat("\n1. Loading and processing data...\n")
  raw_data <- load_data(config)
  
  # Process movement data
  raw_data$movement_processed <- process_movement_data(raw_data$movement)
  
  # Define movement variables
  movement_vars <- c(
    "SB_time_workday", "LPA_time_workday", "MVPA_time_workday", "Sleep_time_workday",
    "SB_time_weekend", "LPA_time_weekend", "MVPA_time_weekend", "Sleep_time_weekend"
  )
  
  # Define covariates
  covariates <- c("Sex", "age", "BMI")
  
  # Create results directory
  if (!dir.exists(config$results_dir)) {
    dir.create(config$results_dir, recursive = TRUE)
  }
  
  # Run analysis for each disease in parallel
  cat(sprintf("\n2. Analyzing %d diseases in parallel...\n", length(diseases)))
  
  results <- future_map(diseases, function(disease) {
    # Clean memory
    gc()
    
    # Analyze single disease
    disease_result <- analyze_single_disease(
      disease = disease,
      movement_vars = movement_vars,
      covariates = covariates,
      data_list = raw_data,
      config = config
    )
    
    # Clean memory
    gc()
    
    return(disease_result)
  }, .options = furrr_options(seed = TRUE))
  
  # Name results by disease
  names(results) <- diseases
  
  # Remove NULL results
  results <- results[!sapply(results, is.null)]
  
  # Stop parallel processing
  plan(sequential)
  
  # Print summary
  cat("\n" + rep("=", 70) + "\n")
  cat("ANALYSIS SUMMARY\n")
  cat(rep("=", 70) + "\n")
  
  for (disease in names(results)) {
    result <- results[[disease]]
    cat(sprintf("  %-10s: %4d models, %3d significant, saved to: %s\n",
                disease, result$n_models, result$n_significant, result$output_dir))
  }
  
  cat(sprintf("\nTotal diseases analyzed: %d\n", length(results)))
  cat(sprintf("Results saved to: %s\n", config$results_dir))
  
  return(results)
}

# ==============================================================================
# EXECUTION SCRIPT
# ==============================================================================

# Example execution
if (TRUE) {  # Set to FALSE to skip execution when sourcing
  # Define diseases to analyze (ICD codes)
  diseases_to_analyze <- c("N18", "F32", "E11", "J44", "K74")
  
  # Optional: smaller subset for testing
  # diseases_to_analyze <- c("N18", "F32")
  
  # Update configuration if needed
  CONFIG$max_proteins_per_disease <- 100  # Limit for testing
  
  # Set number of parallel workers (adjust based on your system)
  n_workers <- min(6, parallel::detectCores() - 2)
  
  # Run the analysis
  analysis_results <- main_protein_mediation_analysis(
    diseases = diseases_to_analyze,
    config = CONFIG,
    n_workers = n_workers
  )
  
  cat("\n" + rep("*", 70) + "\n")
  cat("ANALYSIS COMPLETED SUCCESSFULLY\n")
  cat(rep("*", 70) + "\n")
}

# ==============================================================================
# HELPER FUNCTIONS FOR DATA EXPLORATION
# ==============================================================================

explore_data <- function(config) {
  # Quick data exploration function
  
  cat("Data Exploration\n")
  cat(rep("-", 50), "\n", sep = "")
  
  # Load data
  data_list <- load_data(config)
  
  # Basic statistics
  cat("\n1. Basic Statistics:\n")
  
  # Movement data
  movement_cols <- ncol(data_list$movement)
  cat(sprintf("  Movement features: %d columns\n", movement_cols))
  if (movement_cols >= 192) {
    cat("    Expected format: 192 columns (48 hours × 4 states)\n")
  }
  
  # Protein data
  cat(sprintf("  Protein variables: %d\n", ncol(data_list$prot) - 1))
  
  # Disease data
  disease_cols <- colnames(data_list$icd)[-1]  # Exclude Participant.ID
  cat(sprintf("  Disease variables: %d\n", length(disease_cols)))
  
  # Sample top diseases by frequency
  if (length(disease_cols) > 0) {
    disease_counts <- colSums(data_list$icd[-1], na.rm = TRUE)
    top_diseases <- head(sort(disease_counts, decreasing = TRUE), 10)
    
    cat("\n2. Top 10 diseases by frequency:\n")
    for (i in 1:length(top_diseases)) {
      cat(sprintf("  %2d. %-10s: %5d cases\n", 
                  i, names(top_diseases)[i], top_diseases[i]))
    }
  }
  
  # Participant overlap
  common_ids <- Reduce(intersect, list(
    data_list$imm$Participant.ID,
    data_list$icd$Participant.ID,
    data_list$movement$Participant.ID,
    data_list$prot$Participant.ID
  ))
  
  cat(sprintf("\n3. Participant overlap: %d common participants\n", 
              length(common_ids)))
  
  return(data_list)
}

# ==============================================================================
# SESSION INFO
# ==============================================================================

print_session_info <- function() {
  # Print session information for reproducibility
  
  cat("\n" + rep("=", 70) + "\n")
  cat("SESSION INFORMATION\n")
  cat(rep("=", 70) + "\n")
  
  cat("R version:", R.version.string, "\n")
  cat("Platform:", R.version$platform, "\n")
  cat("\nLoaded packages:\n")
  
  # Get package versions
  packages <- c("lavaan", "tidyverse", "progress", "furrr", "purrr")
  for (pkg in packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("  %-15s v%s\n", pkg, packageVersion(pkg)))
    }
  }
  
  # Memory usage
  mem_usage <- pryr::mem_used()
  cat(sprintf("\nMemory usage: %.1f MB\n", mem_usage / 1024^2))
}

# Uncomment to explore data
# explore_data(CONFIG)

# Uncomment to print session info
# print_session_info()