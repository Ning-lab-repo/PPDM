# ==============================================================================
# Mendelian Randomization (MR) Analysis Pipeline
# 
# Description: Batch MR analysis for multiple exposures and outcomes using
# TwoSampleMR package. Processes GWAS summary statistics, performs LD clumping,
# harmonizes data, and runs multiple MR methods.
# 
# Dependencies: dplyr, TwoSampleMR, purrr, readr
# ==============================================================================

# Load required libraries
library(dplyr)
library(TwoSampleMR)
library(purrr)
library(readr)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set virtual paths for GitHub
CONFIG <- list(
  # Input directories
  gwas_dir = "./data/gwas/exposures",           # Exposure GWAS files
  ld_dir = "./data/ld_clumped",                 # LD clumped results
  outcome_dir = "./data/gwas/outcomes",         # Outcome GWAS files (.gz format)
  
  # Output directory
  output_dir = "./results/mr_analysis",
  
  # File patterns
  exposure_pattern = "^assoc\\.regenie\\.merged_.+\\.txt$",
  outcome_pattern = "^finngen_R12_.+\\.gz$",
  ld_pattern = "^plink_all_ld_clumped_.+_ivs\\.txt$",
  
  # Analysis parameters
  mr_methods = c("mr_ivw", "mr_ivw_mre", "mr_weighted_median", "mr_egger_regression"),
  pval_threshold = 0.05,
  
  # Column names in GWAS files (adjust based on your data format)
  exposure_cols = list(
    snp = "ID",
    effect_allele = "ALLELE1",
    other_allele = "ALLELE0",
    beta = "BETA",
    se = "SE",
    pval = "LOG10P",           # Note: This is -log10(p) in original data
    eaf = "A1FREQ"
  ),
  
  outcome_cols = list(
    snp = "rsids",
    effect_allele = "alt",
    other_allele = "ref",
    beta = "beta",
    se = "sebeta",
    pval = "pval",
    eaf = "af_alt"
  )
)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

create_output_directory <- function(path) {
  # Create output directory if it doesn't exist
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    cat("Created directory:", path, "\n")
  }
  return(path)
}

read_exposure_gwas <- function(file_path, ld_snps, config) {
  # Read and process exposure GWAS data
  
  tryCatch({
    # Read GWAS file
    gwas_data <- read.table(
      file_path,
      header = TRUE,
      stringsAsFactors = FALSE
    )
    
    # Filter by LD-clumped SNPs
    gwas_filtered <- gwas_data %>%
      filter(.data[[config$exposure_cols$snp]] %in% ld_snps)
    
    if (nrow(gwas_filtered) == 0) {
      cat("  No SNPs overlap with LD file\n")
      return(NULL)
    }
    
    # Convert p-value from -log10(p) to p-value
    gwas_filtered <- gwas_filtered %>%
      mutate(
        pval = 10^(-.data[[config$exposure_cols$pval]]),
        effect_allele = as.character(.data[[config$exposure_cols$effect_allele]]),
        other_allele = as.character(.data[[config$exposure_cols$other_allele]])
      ) %>%
      select(
        SNP = .data[[config$exposure_cols$snp]],
        effect_allele,
        other_allele,
        beta = .data[[config$exposure_cols$beta]],
        se = .data[[config$exposure_cols$se]],
        pval,
        eaf = .data[[config$exposure_cols$eaf]]
      )
    
    return(gwas_filtered)
    
  }, error = function(e) {
    cat("  Error reading exposure GWAS:", e$message, "\n")
    return(NULL)
  })
}

read_outcome_gwas <- function(file_path, exposure_snps, config) {
  # Read and process outcome GWAS data (compressed format)
  
  tryCatch({
    # Read compressed GWAS file
    gwas_data <- read.table(
      gzfile(file_path),
      header = TRUE,
      sep = "\t",
      fill = TRUE,
      quote = "",
      stringsAsFactors = FALSE,
      comment.char = ""
    )
    
    # Filter by SNPs present in exposure data
    gwas_filtered <- gwas_data %>%
      filter(.data[[config$outcome_cols$snp]] %in% exposure_snps)
    
    if (nrow(gwas_filtered) == 0) {
      cat("  No SNPs overlap with exposure data\n")
      return(NULL)
    }
    
    # Process columns
    gwas_filtered <- gwas_filtered %>%
      mutate(
        pval = as.numeric(.data[[config$outcome_cols$pval]]),
        beta = as.numeric(.data[[config$outcome_cols$beta]]),
        se = as.numeric(.data[[config$outcome_cols$se]]),
        eaf = as.numeric(.data[[config$outcome_cols$eaf]]),
        effect_allele = as.character(.data[[config$outcome_cols$effect_allele]]),
        other_allele = as.character(.data[[config$outcome_cols$other_allele]])
      ) %>%
      select(
        SNP = .data[[config$outcome_cols$snp]],
        effect_allele,
        other_allele,
        beta,
        se,
        pval,
        eaf
      )
    
    return(gwas_filtered)
    
  }, error = function(e) {
    cat("  Error reading outcome GWAS:", e$message, "\n")
    return(NULL)
  })
}

perform_mr_analysis <- function(exposure_data, outcome_data, result_dir, config) {
  # Perform full MR analysis pipeline
  
  # Save formatted data
  exposure_file <- file.path(result_dir, "exposure_mr_format.csv")
  outcome_file <- file.path(result_dir, "outcome_mr_format.csv")
  
  write.csv(exposure_data, exposure_file, row.names = FALSE, quote = FALSE)
  write.csv(outcome_data, outcome_file, row.names = FALSE, quote = FALSE)
  
  # Load data for MR analysis
  tryCatch({
    exposure_dat <- read_exposure_data(exposure_file, sep = ",")
    outcome_dat <- read_outcome_data(outcome_file, sep = ",")
    
    # Harmonize data
    dat <- harmonise_data(exposure_dat, outcome_dat)
    
    if (nrow(dat) == 0) {
      cat("  No SNPs after harmonization\n")
      return(NULL)
    }
    
    # Save harmonized data
    write.csv(dat, file.path(result_dir, "harmonised_data.csv"), row.names = FALSE)
    
    # Main MR analysis
    mr_results <- mr(dat, method_list = config$mr_methods)
    
    # Heterogeneity tests
    heterogeneity <- mr_heterogeneity(dat)
    
    # Pleiotropy test (Egger intercept)
    pleiotropy <- mr_pleiotropy_test(dat)
    
    # Save results
    write.csv(mr_results, file.path(result_dir, "mr_results.csv"), row.names = FALSE)
    write.csv(heterogeneity, file.path(result_dir, "heterogeneity_results.csv"), row.names = FALSE)
    write.csv(pleiotropy, file.path(result_dir, "pleiotropy_results.csv"), row.names = FALSE)
    
    return(list(
      mr_results = mr_results,
      heterogeneity = heterogeneity,
      pleiotropy = pleiotropy,
      harmonised_data = dat
    ))
    
  }, error = function(e) {
    cat("  MR analysis failed:", e$message, "\n")
    return(NULL)
  })
}

# ==============================================================================
# BATCH PROCESSING FUNCTIONS
# ==============================================================================

process_exposure_outcome_pair <- function(gwas_file, outcome_file, config) {
  # Process single exposure-outcome pair
  
  # Extract file tags for naming
  exposure_tag <- gsub("^assoc\\.regenie\\.merged_|\\.txt$", "", basename(gwas_file))
  outcome_tag <- gsub("^finngen_R12_|\\.gz$", "", basename(outcome_file))
  
  cat("\nProcessing:", exposure_tag, "vs", outcome_tag, "\n")
  cat(rep("-", 50), "\n", sep = "")
  
  # Find corresponding LD file
  ld_pattern <- gsub("_", ".*", exposure_tag)  # Flexible matching
  ld_files <- list.files(
    config$ld_dir, 
    pattern = paste0("plink_all_ld_clumped_.*", exposure_tag, ".*_ivs\\.txt$"),
    full.names = TRUE
  )
  
  if (length(ld_files) == 0) {
    cat("  No LD file found for", exposure_tag, "\n")
    return(NULL)
  }
  
  # Use first matching LD file
  ld_file <- ld_files[1]
  cat("  Using LD file:", basename(ld_file), "\n")
  
  # Read LD-clumped SNPs
  ld_snps <- tryCatch({
    readLines(ld_file)
  }, error = function(e) {
    cat("  Error reading LD file:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(ld_snps) || length(ld_snps) == 0) {
    cat("  No SNPs in LD file\n")
    return(NULL)
  }
  
  # Read exposure data
  exposure_data <- read_exposure_gwas(gwas_file, ld_snps, config)
  if (is.null(exposure_data)) {
    return(NULL)
  }
  
  cat("  Exposure SNPs:", nrow(exposure_data), "\n")
  
  # Read outcome data
  outcome_data <- read_outcome_gwas(outcome_file, exposure_data$SNP, config)
  if (is.null(outcome_data)) {
    return(NULL)
  }
  
  cat("  Outcome SNPs:", nrow(outcome_data), "\n")
  
  # Create result directory
  result_dir <- file.path(config$output_dir, paste0(exposure_tag, "_vs_", outcome_tag))
  create_output_directory(result_dir)
  
  # Perform MR analysis
  mr_results <- perform_mr_analysis(exposure_data, outcome_data, result_dir, config)
  
  if (is.null(mr_results)) {
    return(NULL)
  }
  
  # Return summary
  list(
    exposure = exposure_tag,
    outcome = outcome_tag,
    result_dir = result_dir,
    n_snps = nrow(mr_results$harmonised_data),
    mr_results = mr_results$mr_results,
    heterogeneity = mr_results$heterogeneity,
    pleiotropy = mr_results$pleiotropy
  )
}

# ==============================================================================
# MAIN ANALYSIS PIPELINE
# ==============================================================================

run_batch_mr_analysis <- function(config) {
  # Main batch MR analysis pipeline
  
  cat("\n" + rep("=", 70) + "\n")
  cat("MENDELIAN RANDOMIZATION BATCH ANALYSIS\n")
  cat(rep("=", 70) + "\n")
  
  # Create output directory
  create_output_directory(config$output_dir)
  
  # Get file lists
  gwas_files <- list.files(
    config$gwas_dir, 
    pattern = config$exposure_pattern, 
    full.names = TRUE
  )
  
  outcome_files <- list.files(
    config$outcome_dir, 
    pattern = config$outcome_pattern, 
    full.names = TRUE
  )
  
  cat("\nFound files:\n")
  cat("  Exposures:", length(gwas_files), "\n")
  cat("  Outcomes:", length(outcome_files), "\n")
  
  if (length(gwas_files) == 0 || length(outcome_files) == 0) {
    cat("Insufficient files for analysis\n")
    return(NULL)
  }
  
  # Run all combinations
  all_results <- list()
  total_combinations <- length(gwas_files) * length(outcome_files)
  current_combo <- 0
  
  for (gwas_file in gwas_files) {
    for (outcome_file in outcome_files) {
      current_combo <- current_combo + 1
      cat(sprintf("\n[%d/%d] ", current_combo, total_combinations))
      
      # Process pair
      result <- process_exposure_outcome_pair(gwas_file, outcome_file, config)
      
      if (!is.null(result)) {
        all_results[[paste0(result$exposure, "_vs_", result$outcome)]] <- result
        cat("  ✓ Completed\n")
      } else {
        cat("  ✗ Failed\n")
      }
    }
  }
  
  # Generate summary
  if (length(all_results) > 0) {
    summary_df <- generate_summary_table(all_results, config)
    save_summary_results(summary_df, config$output_dir)
    print_summary_statistics(all_results)
    
    return(list(
      results = all_results,
      summary = summary_df,
      output_dir = config$output_dir
    ))
  } else {
    cat("\nNo successful MR analyses completed\n")
    return(NULL)
  }
}

generate_summary_table <- function(results_list, config) {
  # Generate summary table from all results
  
  cat("\nGenerating summary table...\n")
  
  summary_data <- lapply(names(results_list), function(result_name) {
    result <- results_list[[result_name]]
    
    # Extract key statistics
    mr_res <- result$mr_results
    het_res <- result$heterogeneity
    pleio_res <- result$pleiotropy
    
    # Find rows for different methods
    ivw_row <- which(mr_res$method == "Inverse variance weighted")
    ivw_mre_row <- which(mr_res$method == "Inverse variance weighted (multiplicative random effects)")
    wm_row <- which(mr_res$method == "Weighted median")
    egger_row <- which(mr_res$method == "MR Egger")
    
    # Find heterogeneity row
    het_row <- which(het_res$method == "Inverse variance weighted")
    
    # Create summary row
    data.frame(
      Exposure = result$exposure,
      Outcome = result$outcome,
      N_SNPs = result$n_snps,
      
      # IVW results
      IVW_beta = if (length(ivw_row) > 0) mr_res$b[ivw_row] else NA,
      IVW_se = if (length(ivw_row) > 0) mr_res$se[ivw_row] else NA,
      IVW_pval = if (length(ivw_row) > 0) mr_res$pval[ivw_row] else NA,
      
      # IVW-MRE results
      IVW_MRE_beta = if (length(ivw_mre_row) > 0) mr_res$b[ivw_mre_row] else NA,
      IVW_MRE_se = if (length(ivw_mre_row) > 0) mr_res$se[ivw_mre_row] else NA,
      IVW_MRE_pval = if (length(ivw_mre_row) > 0) mr_res$pval[ivw_mre_row] else NA,
      
      # Weighted median results
      WM_beta = if (length(wm_row) > 0) mr_res$b[wm_row] else NA,
      WM_se = if (length(wm_row) > 0) mr_res$se[wm_row] else NA,
      WM_pval = if (length(wm_row) > 0) mr_res$pval[wm_row] else NA,
      
      # MR Egger results
      Egger_beta = if (length(egger_row) > 0) mr_res$b[egger_row] else NA,
      Egger_se = if (length(egger_row) > 0) mr_res$se[egger_row] else NA,
      Egger_pval = if (length(egger_row) > 0) mr_res$pval[egger_row] else NA,
      Egger_intercept = if (!is.null(pleio_res)) pleio_res$egger_intercept else NA,
      Egger_intercept_pval = if (!is.null(pleio_res)) pleio_res$pval else NA,
      
      # Heterogeneity
      Cochran_Q = if (length(het_row) > 0) het_res$Q[het_row] else NA,
      Q_pval = if (length(het_row) > 0) het_res$Q_pval[het_row] else NA,
      Q_df = if (length(het_row) > 0) het_res$Q_df[het_row] else NA,
      
      stringsAsFactors = FALSE
    )
  })
  
  # Combine all summaries
  summary_df <- do.call(rbind, summary_data)
  
  # Add significance flags
  summary_df <- summary_df %>%
    mutate(
      IVW_significant = ifelse(IVW_pval < config$pval_threshold, "Yes", "No"),
      WM_significant = ifelse(WM_pval < config$pval_threshold, "Yes", "No"),
      Egger_significant = ifelse(Egger_pval < config$pval_threshold, "Yes", "No"),
      Heterogeneity = ifelse(Q_pval < config$pval_threshold, "Present", "Absent"),
      Pleiotropy = ifelse(Egger_intercept_pval < config$pval_threshold, "Present", "Absent")
    )
  
  return(summary_df)
}

save_summary_results <- function(summary_df, output_dir) {
  # Save summary results to CSV
  
  if (nrow(summary_df) == 0) {
    cat("No results to save\n")
    return(NULL)
  }
  
  summary_file <- file.path(output_dir, "summary_results.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  
  cat("Summary results saved to:", summary_file, "\n")
  
  # Also save simplified version
  simplified <- summary_df %>%
    select(Exposure, Outcome, N_SNPs, 
           IVW_beta, IVW_pval, IVW_significant,
           WM_beta, WM_pval, WM_significant,
           Heterogeneity, Pleiotropy)
  
  simplified_file <- file.path(output_dir, "simplified_summary.csv")
  write.csv(simplified, simplified_file, row.names = FALSE)
  
  return(list(
    full_summary = summary_file,
    simplified_summary = simplified_file
  ))
}

print_summary_statistics <- function(results_list) {
  # Print summary statistics
  
  cat("\n" + rep("=", 70) + "\n")
  cat("ANALYSIS SUMMARY STATISTICS\n")
  cat(rep("=", 70) + "\n")
  
  n_analyses <- length(results_list)
  
  if (n_analyses == 0) {
    cat("No analyses completed\n")
    return(NULL)
  }
  
  # Count significant results
  significant_counts <- list(
    IVW = 0,
    WM = 0,
    Egger = 0
  )
  
  total_snps <- 0
  
  for (result in results_list) {
    mr_res <- result$mr_results
    
    ivw_row <- which(mr_res$method == "Inverse variance weighted")
    wm_row <- which(mr_res$method == "Weighted median")
    egger_row <- which(mr_res$method == "MR Egger")
    
    if (length(ivw_row) > 0 && mr_res$pval[ivw_row] < 0.05) {
      significant_counts$IVW <- significant_counts$IVW + 1
    }
    
    if (length(wm_row) > 0 && mr_res$pval[wm_row] < 0.05) {
      significant_counts$WM <- significant_counts$WM + 1
    }
    
    if (length(egger_row) > 0 && mr_res$pval[egger_row] < 0.05) {
      significant_counts$Egger <- significant_counts$Egger + 1
    }
    
    total_snps <- total_snps + result$n_snps
  }
  
  # Print statistics
  cat("\nCompleted analyses:", n_analyses, "\n")
  cat("Average SNPs per analysis:", round(total_snps / n_analyses, 1), "\n")
  
  cat("\nSignificant results (p < 0.05):\n")
  cat("  IVW:", significant_counts$IVW, "/", n_analyses, 
      sprintf("(%.1f%%)\n", significant_counts$IVW / n_analyses * 100))
  cat("  Weighted Median:", significant_counts$WM, "/", n_analyses,
      sprintf("(%.1f%%)\n", significant_counts$WM / n_analyses * 100))
  cat("  MR Egger:", significant_counts$Egger, "/", n_analyses,
      sprintf("(%.1f%%)\n", significant_counts$Egger / n_analyses * 100))
}

# ==============================================================================
# EXAMPLE EXECUTION
# ==============================================================================

if (FALSE) {  # Set to TRUE to run example
  # Run batch MR analysis
  results <- run_batch_mr_analysis(CONFIG)
  
  # Print session info
  sessionInfo()
}

# ==============================================================================
# HELPER FUNCTIONS FOR DATA EXPLORATION
# ==============================================================================

explore_data_files <- function(config) {
  # Explore available data files
  
  cat("Exploring data files...\n")
  
  # List exposure files
  exposure_files <- list.files(config$gwas_dir, pattern = config$exposure_pattern)
  cat("\nExposure files (", length(exposure_files), "):\n", sep = "")
  for (file in head(exposure_files, 10)) {
    cat("  ", file, "\n")
  }
  if (length(exposure_files) > 10) {
    cat("  ... and", length(exposure_files) - 10, "more\n")
  }
  
  # List outcome files
  outcome_files <- list.files(config$outcome_dir, pattern = config$outcome_pattern)
  cat("\nOutcome files (", length(outcome_files), "):\n", sep = "")
  for (file in head(outcome_files, 10)) {
    cat("  ", file, "\n")
  }
  if (length(outcome_files) > 10) {
    cat("  ... and", length(outcome_files) - 10, "more\n")
  }
  
  # List LD files
  ld_files <- list.files(config$ld_dir, pattern = config$ld_pattern)
  cat("\nLD clumped files (", length(ld_files), "):\n", sep = "")
  for (file in head(ld_files, 5)) {
    cat("  ", file, "\n")
  }
  if (length(ld_files) > 5) {
    cat("  ... and", length(ld_files) - 5, "more\n")
  }
  
  # Calculate total combinations
  total_combos <- length(exposure_files) * length(outcome_files)
  cat("\nTotal exposure-outcome combinations:", total_combos, "\n")
  
  return(list(
    exposures = exposure_files,
    outcomes = outcome_files,
    ld_files = ld_files
  ))
}

# ==============================================================================
# SESSION INFO FUNCTION
# ==============================================================================

print_mr_session_info <- function() {
  # Print MR-specific session information
  
  cat("\n" + rep("=", 70) + "\n")
  cat("MR ANALYSIS SESSION INFORMATION\n")
  cat(rep("=", 70) + "\n")
  
  cat("R version:", R.version.string, "\n")
  cat("Platform:", R.version$platform, "\n\n")
  
  cat("Key package versions:\n")
  packages <- c("TwoSampleMR", "dplyr", "purrr", "readr")
  
  for (pkg in packages) {
    try({
      version <- packageVersion(pkg)
      cat(sprintf("  %-15s v%s\n", pkg, version))
    }, silent = TRUE)
  }
  
  # TwoSampleMR specific info
  cat("\nTwoSampleMR available methods:\n")
  methods <- mr_method_list()
  cat("  ", length(methods$obj), "MR methods available\n")
  cat("  ", length(methods$heterogeneity), "heterogeneity tests available\n")
}

# ==============================================================================
# MAIN EXECUTION (COMMENTED OUT FOR SAFETY)
# ==============================================================================

# Uncomment to run the analysis
# explore_data_files(CONFIG)
# print_mr_session_info()
# results <- run_batch_mr_analysis(CONFIG)