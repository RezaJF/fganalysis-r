#' @import lme4
#' @import dplyr
#' @importFrom dplyr %>% group_by filter mutate select left_join distinct summarise
NULL

#' @title Calculate BLUP slopes for lab measurements over age
#' @description Implements a linear mixed model (LMM) to calculate Best Linear Unbiased Predictors (BLUPs)
#' for individual-specific slopes of lab value changes over age, following the methodology from
#' Wiegrebe et al. (2024) Nature Communications.
#' @param data Either a `drug.response` object or a data frame with lab measurements.
#' If a data frame, it must contain columns: FINNGENID, OMOP_CONCEPT_ID, EVENT_AGE, and MEASUREMENT_VALUE_HARMONIZED
#' @param output_dir Directory where output files will be saved. Defaults to current directory
#' @param min_measurements Minimum number of measurements per individual to include in analysis (default: 2)
#' @param include_sex Logical indicating whether to include sex as a fixed effect in the model (default: TRUE).
#' If TRUE and SEX column is not found, an error will be raised with instructions.
#' @param debug_dir Optional path to a directory where debugging information (problematic models and data) will be saved. If `NULL` (default), no debug files are saved.
#' @param drug_exposed_only Logical indicating whether to restrict analysis to individuals with recorded drug purchases (default: FALSE).
#' When TRUE, separate analyses are performed for each unique ATC code.
#' @param calculate_post_variance Logical indicating whether to calculate variance of lab values in the post-drug period (default: FALSE).
#' @param calculate_qc Logical indicating whether to calculate quality control metrics including correlation between fixed-effect and BLUP slopes (default: FALSE).
#' @param normalize_variance Logical indicating whether to add quantile normalized variance to output files (default: FALSE).
#' @param save_model Logical indicating whether to save the fitted lmer model object as an RDS file for each OMOP concept (default: FALSE).
#' @param plot_blup_correlation Logical indicating whether to create a scatter plot comparing BLUP slopes with fixed-effect slopes, including correlation coefficient and p-value (default: FALSE). Requires ggpubr package.
#' @return A list containing BLUP results for each OMOP_CONCEPT_ID
#' @details The model fitted is:
#' lab_value_i,t = β0 + β1*sex_i + β2*age_i,t + γ0i + γ1i*age_i,t + ε_i,t
#' where γ0i and γ1i are random intercept and slope for individual i.
#' Sex is coded following PLINK/REGENIE standard: 1=Male, 2=Female, 0=Missing.
#'
#' The function handles convergence issues by:
#' 1. Scaling both age and lab values by their respective standard deviations
#' 2. Using the bobyqa optimizer with increased iterations
#' 3. Trying two model specifications:
#'    - Full model with correlated random intercepts and slopes
#'    - Model with uncorrelated random intercepts and slopes
#' 4. If neither converges, the analysis is skipped for that OMOP concept
#'
#' The output slopes are in original units (lab value change per year).
#' @export
calculate_blup_slopes <- function(data, output_dir = ".",
                                  min_measurements = 2, include_sex = TRUE,
                                  debug_dir = NULL,
                                  drug_exposed_only = FALSE,
                                  calculate_post_variance = FALSE,
                                  calculate_qc = FALSE,
                                  normalize_variance = FALSE,
                                  save_model = FALSE,
                                  plot_blup_correlation = FALSE) {

  # Check input type and extract lab data accordingly
  is_drug_response <- inherits(data, "drug.reponse")

  if (!is_drug_response && !is.data.frame(data)) {
    stop("Input must be either a drug.reponse object or a data frame with lab measurements.")
  }

  # Validate data frame input
  if (!is_drug_response) {
    required_cols <- c("FINNGENID", "OMOP_CONCEPT_ID", "EVENT_AGE", "MEASUREMENT_VALUE_HARMONIZED")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) {
      stop("Lab measurement data frame is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }

  # Check if lme4 is available
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for this function. Please install it with: install.packages('lme4')")
  }

  # Check if ggpubr is available when plotting is requested
  if (plot_blup_correlation && !requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required for plot_blup_correlation = TRUE. Please install it with: install.packages('ggpubr')")
  }

  # Extract lab measurements based on input type
  if (is_drug_response) {
    lab_data <- data$all_measurements %>%
      filter(!is.na(.data$MEASUREMENT_VALUE_HARMONIZED))

    # Extract period definitions for variance calculation
    after_period <- data$lab_response_period$after_period

    # Extract drug purchases if needed
    all_drug_purchases <- data$all_drug_purchases
  } else {
    # Direct lab measurement input
    lab_data <- data %>%
      filter(!is.na(.data$MEASUREMENT_VALUE_HARMONIZED))

    # No period definitions available for direct lab input
    after_period <- NULL
    all_drug_purchases <- NULL

    # Disable drug-specific features for direct lab input
    if (drug_exposed_only || calculate_post_variance) {
      warning("drug_exposed_only and calculate_post_variance are not supported with direct lab measurement input. These options will be ignored.")
      drug_exposed_only <- FALSE
      calculate_post_variance <- FALSE
    }
  }

  # Check for SEX column if include_sex is TRUE
  if (include_sex) {
    if (!"SEX" %in% colnames(lab_data)) {
      if (is_drug_response) {
        stop("SEX column not found in drug_response object. ",
             "Please run create_drug_response() with covariates = conn$cov_pheno ",
             "and covariate_cols = c('SEX') to include sex data.")
      } else {
        stop("SEX column not found in lab measurement data. ",
             "Please add a SEX column to your data frame or set include_sex = FALSE.")
      }
    }

    # Standardise sex coding to PLINK/REGENIE format (1=Male, 2=Female, 0=Missing)
    lab_data <- lab_data %>%
      mutate(SEX_CODED = case_when(
        toupper(.data$SEX) %in% c("M", "male", "MALE", "1") ~ 1L,
        toupper(.data$SEX) %in% c("F", "female", "FEMALE", "2") ~ 2L,
        TRUE ~ 0L # All other values (including NA) are coded as missing
      ))
  } else {
    # If not including sex, create dummy variable where everyone is male (1)
    lab_data$SEX_CODED <- 1L
  }

  # Get unique OMOP concept IDs
  concept_ids <- unique(lab_data$OMOP_CONCEPT_ID)
  analysis_groups <- data.frame(
    concept_id = concept_ids,
    atc_code = NA_character_,
    stringsAsFactors = FALSE
  )

  # Initialize results list
  blup_results <- list()

  # Process each OMOP concept ID separately
  for (concept_id in concept_ids) {

    cat(paste0("Processing OMOP_CONCEPT_ID: ", concept_id, "\n"))

    # Filter data for current concept
    concept_data <- lab_data %>%
      filter(.data$OMOP_CONCEPT_ID == concept_id)

    # Count measurements per individual
    measurement_counts <- concept_data %>%
      group_by(.data$FINNGENID) %>%
      summarise(n_measurements = n()) %>%
      filter(.data$n_measurements >= min_measurements)

    # Filter to individuals with sufficient measurements
    analysis_data <- concept_data %>%
      filter(.data$FINNGENID %in% measurement_counts$FINNGENID)

    # Check if there's enough data
    n_individuals <- length(unique(analysis_data$FINNGENID))
    if (n_individuals < 10) {
      warning(paste0("Only ", n_individuals, " individuals for OMOP_CONCEPT_ID ",
                     concept_id, ". Skipping analysis."))
      next
    }

    # Standardize variables to improve convergence
    # Scale age by its SD
    age_mean <- mean(analysis_data$EVENT_AGE, na.rm = TRUE)
    age_sd <- sd(analysis_data$EVENT_AGE, na.rm = TRUE)
    if (age_sd == 0) {
      warning(paste0("No variation in age for OMOP_CONCEPT_ID ", concept_id, ". Skipping analysis."))
      next
    }
    analysis_data$EVENT_AGE_SCALED <- (analysis_data$EVENT_AGE - age_mean) / age_sd

    # Scale lab values by their SD
    lab_mean <- mean(analysis_data$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)
    lab_sd <- sd(analysis_data$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)
    if (lab_sd == 0) {
      warning(paste0("No variation in lab values for OMOP_CONCEPT_ID ", concept_id, ". Skipping analysis."))
      next
    }
    analysis_data$LAB_VALUE_SCALED <- (analysis_data$MEASUREMENT_VALUE_HARMONIZED - lab_mean) / lab_sd

    # Fit the linear mixed model
    # Model: lab_value ~ sex + age + (age | FINNGENID)
    # This includes random intercepts and random slopes for age
    model_converged <- FALSE
    tryCatch({
      # First try the full model with correlated random effects
      if (include_sex) {
        # Model with sex as fixed effect
        lmm_model <- lme4::lmer(
          LAB_VALUE_SCALED ~ SEX_CODED + EVENT_AGE_SCALED + (EVENT_AGE_SCALED | FINNGENID),
          data = analysis_data,
          REML = TRUE,
          control = lme4::lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000),
            calc.derivs = FALSE,
            check.nobs.vs.nRE = "ignore"
          )
        )
      } else {
        # Model without sex effect
        lmm_model <- lme4::lmer(
          LAB_VALUE_SCALED ~ EVENT_AGE_SCALED + (EVENT_AGE_SCALED | FINNGENID),
          data = analysis_data,
          REML = TRUE,
          control = lme4::lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000),
            calc.derivs = FALSE,
            check.nobs.vs.nRE = "ignore"
          )
        )
      }

      # Check convergence
      if (!lme4::isSingular(lmm_model) && length(lmm_model@optinfo$conv$lme4$messages) == 0) {
        cat(paste0("  Successfully fitted full model for OMOP_CONCEPT_ID: ", concept_id, "\n"))
        model_converged <- TRUE
      } else {
        warning(paste0("Convergence issues for full model, trying uncorrelated random effects for OMOP_CONCEPT_ID: ", concept_id))

        # Try model with uncorrelated random effects
        if (include_sex) {
          lmm_model <- lme4::lmer(
            LAB_VALUE_SCALED ~ SEX_CODED + EVENT_AGE_SCALED + (1 | FINNGENID) + (0 + EVENT_AGE_SCALED | FINNGENID),
            data = analysis_data,
            REML = TRUE,
            control = lme4::lmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 200000),
              check.nobs.vs.nRE = "ignore"
            )
          )
        } else {
          lmm_model <- lme4::lmer(
            LAB_VALUE_SCALED ~ EVENT_AGE_SCALED + (1 | FINNGENID) + (0 + EVENT_AGE_SCALED | FINNGENID),
            data = analysis_data,
            REML = TRUE,
            control = lme4::lmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 200000),
              check.nobs.vs.nRE = "ignore"
            )
          )
        }

        # Check if this model converged
        if (!lme4::isSingular(lmm_model) && length(lmm_model@optinfo$conv$lme4$messages) == 0) {
          cat(paste0("  Successfully fitted uncorrelated random effects model for OMOP_CONCEPT_ID: ", concept_id, "\n"))
          model_converged <- TRUE
        }
      }

      if (!model_converged) {
        stop("Model failed to converge")
      }

      # Extract BLUPs (random effects)
      random_effects <- tryCatch({
        ranef_obj <- lme4::ranef(lmm_model)
        if (!"FINNGENID" %in% names(ranef_obj)) {
          stop("FINNGENID not found in random effects structure")
        }
        ranef_obj[["FINNGENID"]]
      }, error = function(e) {
        stop(paste0("Failed to extract random effects: ", e$message))
      })

      # Get the slopes from the scaled model
      # The slopes are in units of: SD(lab_value) / SD(age)
      if ("EVENT_AGE_SCALED" %in% colnames(random_effects)) {
        scaled_slopes <- random_effects[, "EVENT_AGE_SCALED"]
      } else {
        # This shouldn't happen if model converged properly
        stop("No random slopes found in converged model")
      }

      # Convert slopes back to original units: lab_value per year
      # scaled_slope = (lab_change/lab_sd) / (age_change/age_sd)
      # original_slope = scaled_slope * (lab_sd / age_sd)
      original_slopes <- scaled_slopes * (lab_sd / age_sd)

      blup_slopes <- data.frame(
        FINNGENID = rownames(random_effects),
        slope = original_slopes,
        stringsAsFactors = FALSE
      )

      # Get all unique FINNGENIDs from the original data
      all_finngenids <- unique(lab_data$FINNGENID)

      # Create output dataframe with all individuals
      # Those not in the analysis get NA slopes
      output_df <- data.frame(
        FID = all_finngenids,
        IID = all_finngenids,
        stringsAsFactors = FALSE
      )

      # Add the slope column with appropriate name
      slope_col_name <- paste0(concept_id, "_slope")
      output_df[[slope_col_name]] <- NA_real_

      # Fill in the calculated slopes
      match_idx <- match(blup_slopes$FINNGENID, output_df$FID)
      output_df[[slope_col_name]][match_idx] <- blup_slopes$slope

      # Save to file
      output_file <- file.path(output_dir, paste0(concept_id, "_DF13.tsv"))
      write.table(output_df,
                  file = output_file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)

      cat(paste0("Saved results to: ", output_file, "\n"))

      # Save model if requested
      model_file <- NULL
      if (save_model) {
        model_file <- file.path(output_dir, paste0(concept_id, "_model.rds"))
        saveRDS(lmm_model, file = model_file)
        cat(paste0("Saved model to: ", model_file, "\n"))
      }

      # Create BLUP vs fixed-effect correlation plot if requested
      plot_file <- NULL
      blup_fixed_correlation <- NULL
      if (plot_blup_correlation || calculate_qc) {
        # Calculate fixed-effect slopes for comparison
        fixed_effect_slopes <- analysis_data %>%
          group_by(.data$FINNGENID) %>%
          do({
            df <- .
            tryCatch({
              if (include_sex && nrow(df) >= 3) {
                lm_fit <- lm(MEASUREMENT_VALUE_HARMONIZED ~ EVENT_AGE + SEX_CODED, data = df)
              } else if (nrow(df) >= 3) {
                lm_fit <- lm(MEASUREMENT_VALUE_HARMONIZED ~ EVENT_AGE, data = df)
              } else {
                return(data.frame(fixed_slope = NA_real_))
              }

              if (length(coef(lm_fit)) >= 2 && !is.na(coef(lm_fit)[2])) {
                data.frame(fixed_slope = coef(lm_fit)[2])  # Age coefficient
              } else {
                data.frame(fixed_slope = NA_real_)
              }
            }, error = function(e) {
              data.frame(fixed_slope = NA_real_)
            })
          }) %>%
          ungroup() %>%
          filter(!is.na(.data$fixed_slope))

        # Merge BLUP and fixed-effect slopes
        comparison_data <- blup_slopes %>%
          left_join(fixed_effect_slopes, by = "FINNGENID") %>%
          filter(!is.na(.data$slope) & !is.na(.data$fixed_slope))

        if (nrow(comparison_data) > 3) {
          # Calculate correlation
          cor_test <- cor.test(comparison_data$slope, comparison_data$fixed_slope)
          blup_fixed_correlation <- list(
            correlation = cor_test$estimate,
            p_value = cor_test$p.value,
            n_pairs = nrow(comparison_data)
          )

          # Create plot if requested
          if (plot_blup_correlation) {
            library(ggpubr)
            p <- ggpubr::ggscatter(comparison_data,
                                   x = "fixed_slope", y = "slope",
                                   add = "reg.line", conf.int = TRUE,
                                   cor.coef = TRUE, cor.method = "pearson",
                                   xlab = "Fixed-Effect Slope (OLS)",
                                   ylab = "BLUP Slope",
                                   title = paste0("BLUP vs Fixed-Effect Slopes: ", concept_id),
                                   subtitle = paste0("n = ", nrow(comparison_data),
                                                   ", r = ", round(cor_test$estimate, 3),
                                                   ", p = ", format.pval(cor_test$p.value, digits = 3))) +
              ggplot2::theme_bw()

            plot_file <- file.path(output_dir, paste0(concept_id, "_blup_correlation.pdf"))
            ggplot2::ggsave(plot_file, plot = p, width = 8, height = 6)
            cat(paste0("Saved correlation plot to: ", plot_file, "\n"))
          }
        } else {
          warning(paste0("Insufficient data for correlation analysis for OMOP_CONCEPT_ID: ", concept_id))
        }
      }

      # Store results
      blup_results[[as.character(concept_id)]] <- list(
        model = lmm_model,
        blup_slopes = blup_slopes,
        n_individuals = n_individuals,
        output_file = output_file,
        model_file = model_file,
        plot_file = plot_file,
        blup_fixed_correlation = blup_fixed_correlation
      )

    }, error = function(e) {
      warning(paste0("Error fitting model for OMOP_CONCEPT_ID ", concept_id,
                     ": ", e$message))
    })
  }

  cat("BLUP calculation completed.\n")

  # Return results invisibly
  invisible(blup_results)
}

#' @title Summarize BLUP slope results
#' @description Provides summary statistics for BLUP slope calculations
#' @param blup_results Results from calculate_blup_slopes function
#' @return A data frame with summary statistics for each OMOP concept
#' @export
summarize_blup_results <- function(blup_results) {

  if (length(blup_results) == 0) {
    return(data.frame(
      OMOP_CONCEPT_ID = character(),
      n_individuals = integer(),
      mean_slope = numeric(),
      sd_slope = numeric(),
      min_slope = numeric(),
      max_slope = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  summary_df <- do.call(rbind, lapply(names(blup_results), function(concept_id) {
    result <- blup_results[[concept_id]]
    slopes <- result$blup_slopes$slope

    data.frame(
      OMOP_CONCEPT_ID = concept_id,
      n_individuals = result$n_individuals,
      mean_slope = mean(slopes, na.rm = TRUE),
      sd_slope = sd(slopes, na.rm = TRUE),
      min_slope = min(slopes, na.rm = TRUE),
      max_slope = max(slopes, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  return(summary_df)
}