#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#' @import stringr
#' @importFrom stats median lm quantile sd
#' @import utils
#' @import grDevices
#' @import UpSetR
#NULL

#' @title data object returned from drug response analyse (create_drug_response)
#' @param responses data frame with response data
#' @param lab_measurements data frame with all lab measurements
#' @param drug_purchases data frame with all drug purchases
#' @param before_period vector with two elements, start and end of the before period
#' @param after_period vector with two elements, start and end of the after period
#' @return object of class drug.response
#' @export
drug.response <- function(responses, lab_measurements, drug_purchases, before_period, after_period) {
  obj <- list(responses=responses, all_measurements=lab_measurements, all_drug_purchases=drug_purchases,
              lab_response_period=list(before_period=before_period, after_period=after_period))
  class(obj) <- "drug.reponse"
  return(obj)
}

#' @title Create drug response object
#' @param conn fg_data_connection object
#' @param lablist vector of lab measurement concept IDs
#' @param druglist vector of drug ATC codes. The ATC codes are matched with the first part of the code (e.g. A01*).
#' @param before_period vector with two elements, start and end of the before period
#' @param after_period vector with two elements, start and end of the after period
#' @param finngen_ids vector of FINNGENIDs to filter the data
#' @param remove_outliers_sd integer, if defined, remove outliers from the lab measurements. The value is the number of standard deviations from the mean to define an outlier (i.e. 1, 2, 3, 4, 5, 6)
#' @param covariates Optional data frame with covariate data (e.g., sex, age at death)
#' @param covariate_cols Optional character vector of column names to add from the covariates data frame
#' @return drug.response object
#' @export
create_drug_response <- function(conn, lablist, druglist,
before_period, after_period, finngen_ids=NULL, remove_outliers_sd=NULL,
covariates = NULL,
covariate_cols = NULL) {

  if(!inherits(conn, "fg_data_connection")) {
    stop("conn must be a fg_data_connection object")
  }

  print("Querying lab measurements...")
  lab_measurements <- get_lab_measurements(all_labs=conn$labs, lablist=lablist, finngen_ids=finngen_ids, require_values = TRUE)

  if (!is.null(remove_outliers_sd)) {
    if (!is.numeric(remove_outliers_sd) || remove_outliers_sd < 1 || remove_outliers_sd > 6) {
      stop("remove_outliers_sd must be an integer between 1 and 6")
    }

    mean_val <- mean(lab_measurements$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)
    sd_val <- sd(lab_measurements$MEASUREMENT_VALUE_HARMONIZED, na.rm = TRUE)

    lower_bound <- mean_val - (remove_outliers_sd * sd_val)
    upper_bound <- mean_val + (remove_outliers_sd * sd_val)

    original_rows <- nrow(lab_measurements)
    lab_measurements <- lab_measurements %>%
      filter(.data$MEASUREMENT_VALUE_HARMONIZED >= lower_bound & .data$MEASUREMENT_VALUE_HARMONIZED <= upper_bound)

    removed_rows <- original_rows - nrow(lab_measurements)
    print(paste("Removed", removed_rows, "outliers based on", remove_outliers_sd, "standard deviations from the mean."))
  }

  all_fg_ids <- unique(c(lab_measurements$FINNGENID, finngen_ids))
  print("Querying purchases...")
  drug_purchases <- get_drug_purchases(conn, druglist, all_fg_ids)

  dr_first_purchase <- drug_purchases %>% group_by(.data$FINNGENID) %>% arrange(.data$EVENT_AGE) %>%
    summarize(n = n(), first_drug_age = first(.data$EVENT_AGE), first_drug=first(.data$ATC))

  lab_measurements <- left_join(lab_measurements, dr_first_purchase, by = "FINNGENID")
  lab_measurements <- lab_measurements %>% mutate(time_to_drug = .data$first_drug_age - .data$EVENT_AGE)

  print("generating response summary...")

  lab_response <- generate_response_summary(lab_measurements, before_period, after_period)
  cat("Number of individuals with response data: ", nrow(lab_response), "\n")

  # Add covariates if provided
  if (!is.null(covariates) && !is.null(covariate_cols)) {
    print(paste("Adding covariates:", paste(covariate_cols, collapse = ", ")))

    # Ensure FINNGENID is in the columns for the join
    cols_to_select <- unique(c("FINNGENID", covariate_cols))

    # Check that all requested columns exist in the covariates dataframe
    missing_cols <- setdiff(cols_to_select, colnames(covariates))
    if (length(missing_cols) > 0) {
      stop(paste("The following `covariate_cols` are not in the `covariates` dataframe:",
                 paste(missing_cols, collapse = ", ")))
    }

    # Select the requested columns and ensure one row per FINNGENID
    # Collect the data if it's a lazy table
    cov_data_to_join <- covariates %>%
      select(all_of(cols_to_select)) %>%
      distinct(.data$FINNGENID, .keep_all = TRUE)

    # If covariates is a lazy table (database connection), collect it
    if (inherits(covariates, "tbl_lazy") || inherits(covariates, "tbl_sql")) {
      cov_data_to_join <- dplyr::collect(cov_data_to_join)
    }

    # Join with the final response summary data frame
    lab_response <- lab_response %>%
      left_join(cov_data_to_join, by = "FINNGENID")

    # Also join with the full measurements data frame
    lab_measurements <- lab_measurements %>%
      left_join(cov_data_to_join, by = "FINNGENID")
  }

  return(drug.response(responses = lab_response, lab_measurements=lab_measurements,
                              drug_purchases=drug_purchases, before_period = before_period, after_period = after_period))
}

quant_text <- function(vector) {
  paste0( paste(names(quantile(vector)), sep="\t"), ":", paste(quantile(vector),sep="\t"), collapse=" ")
}


#' @title Summarize drug response
#' @description Summarize drug response data created with create_drug_response. writes plots and tables to disk
#' @param drug_response drug.response object
#' @param out_file_prefix prefix for output files
#' @return NULL
#' @export
summarize_drug_response <- function(drug_response, out_file_prefix) {

  labs <- drug_response$all_measurements %>% filter(!is.na(.data$MEASUREMENT_VALUE_HARMONIZED))
  responses <- drug_response$responses
  drugs <- drug_response$all_drug_purchases

  pdf(paste0(out_file_prefix, ".pdf"), width=10, height=6)

  inds_with_lab <- length(unique(labs$FINNGENID))
  n_lab_meas <- nrow(labs)

  inds_with_drugs <- length(unique(drugs$FINNGENID))
  n_drugs_meas <- nrow(drugs)

  inds_in_analysis <- nrow(responses)
  n_range_before <- quant_text(responses$n_before)
  n_range_after <- quant_text(responses$n_before)

  range_baseline_age <- quant_text(responses$n_before)

  plot(ggtexttable(data.frame("Group"=c("labs","drugs","in analysis"),
                         "N"= c(inds_with_lab, inds_with_drugs, inds_in_analysis),
                         "N events" = c(n_lab_meas, n_drugs_meas, inds_in_analysis)), rows=NULL))


  per_drug <- responses %>% group_by(.data$first_drug) %>%
    summarise(N=n(),p=summary(lm("response ~ 1", data=pick(.data$FINNGENID,.data$response)))$coefficients[1,4],sd=sd(.data$response),
    response=mean(.data$response),
    purch_age_dist=quant_text(.data$baseline_age))

  all_resp <- rbind(per_drug, data.frame(first_drug="All drugs", N=inds_in_analysis,
            response=mean(responses$response), p=summary(lm("response ~ 1", data=responses))$coefficients[1,4],
            purch_age_dist=quant_text(responses$baseline_age), sd=sd(responses$response)))

  write.table(all_resp %>% arrange(desc(.data$N)) %>%
    select(.data$first_drug, .data$N, .data$response, .data$p, .data$purch_age_dist),
    paste0(out_file_prefix, "_responses_by_drug.txt"), sep="\t", row.names=FALSE, quote=FALSE)


  plot(ggtexttable(responses %>% group_by(.data$first_drug) %>%
    summarise(n_purch=n(), n_indiv=length(unique(.data$FINNGENID)),
     p=summary(lm("response ~ 1", data=pick(.data$FINNGENID,.data$response)))$coefficients[1,4],
     response=mean(.data$response),
     purch_age_dist=quant_text(.data$baseline_age)) %>%
     select(.data$first_drug, .data$n_purch, .data$response, .data$p, .data$purch_age_dist) %>%
    arrange(desc(.data$n_purch))))


  begin <- ceiling(max(min(-labs$time_to_drug, na.rm = TRUE), drug_response$lab_response_period$before_period[1]))
  end <- ceiling(min(max(-labs$time_to_drug, na.rm = TRUE), drug_response$lab_response_period$after_period[2]))
  labs$bin <- cut(-labs$time_to_drug,
                  breaks=seq(begin, end, by=.25),
                  include.lowest = TRUE)

  plot(ggplot(labs %>% filter(!is.na(.data$bin))) + geom_boxplot(aes(x=.data$bin, y=.data$MEASUREMENT_VALUE_HARMONIZED)) +
      labs(x="Time to drug purchase (years)", y="Lab measurement") +
      ggtitle("Lab measurements before and after drug purchase")) +  theme_bw() +
      theme(axis.text.x = element_text(angle = 45, size=20)) +


  write.table(labs %>% group_by(.data$bin) %>%
    summarise(n=n(), mean=mean(.data$MEASUREMENT_VALUE_HARMONIZED), sd=sd(.data$MEASUREMENT_VALUE_HARMONIZED)),
    paste0(out_file_prefix, "_labs_by_time_to_drug.txt"), sep="\t", row.names=FALSE, quote=FALSE)


  plot(ggplot(responses) + geom_histogram(aes(x=.data$response)) + theme_bw() +
      labs(x="Response (after - before)", y="Count") +
      ggtitle("Distribution of drug response"))

  fit <- summary(lm(after ~ before, data=responses))

  slope <- format(fit$coefficients["before","Estimate"], digits=2)
  r2 <- format(fit$r.squared, digits=2)
  p <- format(fit$coefficients["before","Pr(>|t|)"], digits=2, scientific=TRUE)

  suppressWarnings(print(ggplot(responses, aes(x=.data$before, y=.data$after)) + geom_point() + geom_smooth(method="lm") +
      geom_abline(slope=1) +
      ggtitle(paste0("Before vs. after values. Slope: ", slope, " R2: ", r2, " p: ", p))))

  dev.off()

  print(paste0("Created summary plots and tables with prefix: ", out_file_prefix))

}


#' @title Generate response summary
#' @param lab_measurements data frame with lab measurements
#' @param before_period vector with two elements, start and end of the before period
#' @param after_period vector with two elements, start and end of the after period
#' @param summary_function function to summarize the lab measurements (default is median)
#' @return data frame with response summary
#' @export
generate_response_summary <- function(lab_measurements, before_period, after_period, summary_function=median) {

  lab_measurements <- lab_measurements %>% mutate(lab_period = case_when(
    dplyr::between(.data$time_to_drug, after_period[1], after_period[2]) ~ 'Before',
    dplyr::between(.data$time_to_drug, before_period[1], before_period[2]) ~ 'After',
    TRUE ~ NA_character_
  ))

  lab_response <- lab_measurements %>% dplyr::filter(!is.na(.data$lab_period) & !is.na(.data$MEASUREMENT_VALUE_HARMONIZED)) %>%
    dplyr::group_by(.data$FINNGENID) %>%
    dplyr::summarize(before = summary_function(.data$MEASUREMENT_VALUE_HARMONIZED[.data$lab_period=='Before'], na.rm=TRUE),
              after = summary_function(.data$MEASUREMENT_VALUE_HARMONIZED[.data$lab_period=='After'], na.rm=TRUE),
              n_before = length(.data$MEASUREMENT_VALUE_HARMONIZED[.data$lab_period=='Before']),
              n_after = length(.data$MEASUREMENT_VALUE_HARMONIZED[.data$lab_period=='After']),
              baseline_age = first(.data$first_drug_age),
              first_drug = first(.data$first_drug),
              response = .data$after - .data$before) %>% dplyr::filter(!is.na(.data$response))

  return(lab_response)
}

#' @title Get lab measurements from FinnGen data
#' @param all_labs data frame with lab measurements
#' @param lablist vector of lab measurement concept IDs
#' @param require_values logical, if TRUE, only return rows with non-missing MEASUREMENT_VALUE_HARMONIZED
#' @param return_cols vector of column names to return
#' @param finngen_ids vector of FINNGENIDs to filter the data
#' @param lazy logical, if TRUE, return a lazy tbl object
#' @return data frame with lab measurements
#' @export
get_lab_measurements <- function(all_labs, lablist, require_values=TRUE, return_cols=c("FINNGENID","OMOP_CONCEPT_ID", "EVENT_AGE", "MEASUREMENT_VALUE_HARMONIZED"),
                                 finngen_ids=NULL, lazy=FALSE) {

  return_cols <- unique(c("OMOP_CONCEPT_ID", return_cols))

  labs <- all_labs %>% select(all_of(return_cols)) %>% dplyr::filter(.data$OMOP_CONCEPT_ID %in% lablist)

  if (!is.null(finngen_ids)) {
    labs <- labs %>% dplyr::filter(.data$FINNGENID %in% finngen_ids)
  }
  if (require_values) {
    labs <- labs %>% dplyr::filter(!is.na(.data$MEASUREMENT_VALUE_HARMONIZED))
  }

  ifelse(lazy, return(labs), return(dplyr::collect(labs)))
}


#' @title Get drug purchases from FinnGen data
#' @param conn finngen data connection object
#' @param druglist vector of drug ATC codes. The ATC codes are matched with the first part of the code (e.g. A01*)
#' @param finngen_ids vector of FINNGENIDs to filter the data. leave empty to get all
#' @param return_cols vector of column names to return
#' @param lazy logical, if TRUE, return a lazy tbl object
#' @return data frame with drug purchases
#' @export
get_drug_purchases <- function(conn, druglist, finngen_ids=NULL,
                               return_cols=c("FINNGENID","EVENT_AGE","APPROX_EVENT_DAY", ATC="CODE1", REIMB_CODE="CODE2", VNR="CODE3", N_PACKS="CODE4"),
                               lazy=FALSE) {

  ## check that conn is a fg_data_connection object and has pheno data name
  if (!inherits(conn, "fg_data_connection")) {
    stop("conn must be a fg_data_connection object")
  }

  if (!"pheno" %in% names(conn)) {
    stop("conn must contain 'pheno' data")
  }

  drugs_regex <- paste0("^(",
                        paste0(druglist, collapse = '|'),
                      ")")

  all_phenos <- conn$pheno

  drugs <- all_phenos %>% dplyr::filter(.data$SOURCE=="PURCH" & str_detect(.data$CODE1, drugs_regex)) %>% select(all_of(return_cols))

  if (!is.null(finngen_ids)) {
    drugs <- drugs %>% dplyr::filter(.data$FINNGENID %in% finngen_ids)
  }

  if ("vnr" %in% names(conn)) {
    columns <- c("VNR","Substance","MedicineName","PackageSize","DDDPerPack","Dosage","DosageUnit")
    vnr <- conn$vnr %>% select(all_of(columns))
    drugs <- left_join(drugs, vnr, by="VNR", copy=TRUE)
  }

  ifelse(lazy, return(drugs), return(dplyr::collect(drugs)))
}


#' @title Get first drug purchase from FinnGen data
#' @param all_phenos data frame with drug purchases
#' @param druglist vector of drug ATC codes. The ATC codes are matched with the first part of the code (e.g. A01*)
#' @param finngen_ids vector of FINNGENIDs to filter the data. leave empty to get all
#' @param return_cols vector of column names to return
#' @param lazy logical, if TRUE, return a lazy tbl object
#' @return data frame with first drug purchases for each FINNGENID
#' @export
get_first_purchase <- function(all_phenos, druglist, finngen_ids=NULL, return_cols=c("FINNGENID","EVENT_AGE","CODE1"),
                               lazy=FALSE) {

  first_purch <- get_drug_purchases(all_phenos, druglist, finngen_ids, return_cols, lazy=TRUE) %>%
    group_by(.data$FINNGENID) %>%
    filter(.data$EVENT_AGE == min(.data$EVENT_AGE)) %>% distinct(.data$EVENT_AGE, .keep_all = TRUE) %>%
    ungroup() %>%
    select(all_of(return_cols))
}

#' @title Generate an UpSet plot of drug purchase combinations
#' @description Creates and saves an UpSet plot to visualize the overlap of purchased drug ATC codes.
#' @param drug_response A `drug.response` object created by `create_drug_response`.
#' @param out_file_prefix A string to use as the prefix for the output PDF file.
#' @return NULL
#' @export
summarize_drug_purchases_upset <- function(drug_response, out_file_prefix) {
  if (!inherits(drug_response, "drug.reponse")) {
    stop("Input must be a drug.reponse object.")
  }

  upset_data <- drug_response$all_drug_purchases %>%
    select("FINNGENID", "ATC") %>%
    distinct() %>%
    mutate(value = 1) %>%
    tidyr::pivot_wider(names_from = "ATC", values_from = "value", values_fill = 0) %>%
    select(-"FINNGENID") %>%
    as.data.frame()

  pdf(paste0(out_file_prefix, "_upset_plot.pdf"), width = 10, height = 7)

  UpSetR::upset(upset_data,
        nsets = ncol(upset_data),
        nintersects = 20,
        mb.ratio = c(0.6, 0.4),
        order.by = "freq",
        decreasing = TRUE,
        text.scale = 1.2,
        mainbar.y.label = "Number of FINNGEN IDs",
        sets.x.label = "Patients per Drug Subtype",
        set_size.show = TRUE)

  dev.off()

  print(paste0("UpSet plot saved to ", out_file_prefix, "_upset_plot.pdf"))
}

#' @title Plot Distribution of Lab Values Before and After Drug Use
#' @description Creates a boxplot comparing the distribution of lab values
#' before and after the first drug purchase, faceted by drug type.
#' @param drug_response A `drug.response` object.
#' @param remove_outliers A logical indicating whether to remove outliers
#' using the 1.5 * IQR rule. Defaults to `FALSE`.
#' @return A `ggplot` object.
#' @export
plot_lab_value_distribution <- function(drug_response, remove_outliers = FALSE) {
  if (!inherits(drug_response, "drug.reponse")) {
    stop("Input must be a drug.reponse object.")
  }

  # Define periods from the response object for consistency
  before_period_def <- drug_response$lab_response_period$before_period
  after_period_def <- drug_response$lab_response_period$after_period

  lab_data_periods <- drug_response$all_measurements %>%
    filter(!is.na(.data$first_drug_age) & !is.na(.data$MEASUREMENT_VALUE_HARMONIZED)) %>%
    mutate(period = case_when(
      between(.data$time_to_drug, before_period_def[1], before_period_def[2]) ~ "Before",
      between(.data$time_to_drug, after_period_def[1], after_period_def[2]) ~ "After",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(.data$period))

  plot_data <- lab_data_periods
  if (remove_outliers) {
    plot_data <- plot_data %>%
      group_by(.data$first_drug, .data$period) %>%
      mutate(
        Q1 = quantile(.data$MEASUREMENT_VALUE_HARMONIZED, 0.25, na.rm = TRUE),
        Q3 = quantile(.data$MEASUREMENT_VALUE_HARMONIZED, 0.75, na.rm = TRUE),
        IQR = .data$Q3 - .data$Q1
      ) %>%
      filter(
        .data$MEASUREMENT_VALUE_HARMONIZED >= (.data$Q1 - 1.5 * .data$IQR) &
        .data$MEASUREMENT_VALUE_HARMONIZED <= (.data$Q3 + 1.5 * .data$IQR)
      ) %>%
      ungroup()
  }

  # Generate plot
  p <- ggplot(plot_data, aes(x = .data$period, y = .data$MEASUREMENT_VALUE_HARMONIZED, fill = .data$period)) +
    ggplot2::geom_boxplot(outlier.shape = if (remove_outliers) NA else 19) +
    ggpubr::stat_compare_means(method = "t.test", label = "p.format", paired = FALSE,
                       label.x = 1.5, label.y.npc = 0.9) +
    labs(
      title = "Distribution of Lab Values Before and After First Drug Purchase",
      x = "Period Relative to Drug Purchase",
      y = "Harmonised Measurement Value"
    ) +
    theme_minimal() +
    facet_wrap(~.data$first_drug, scales = "free_y")

  return(p)
}

#' @title Calculate BLUP slopes for lab measurements over age
#' @description Implements a linear mixed model (LMM) to calculate Best Linear Unbiased Predictors (BLUPs)
#' for individual-specific slopes of lab value changes over age, following the methodology from
#' Wiegrebe et al. (2024) Nature Communications.
#' @param drug_response A `drug.response` object containing lab measurements
#' @param sex_data Optional data frame with columns FINNGENID and SEX (coded as 0/1 or M/F)
#' @param output_dir Directory where output files will be saved. Defaults to current directory
#' @param min_measurements Minimum number of measurements per individual to include in analysis (default: 2)
#' @return A list containing BLUP results for each OMOP_CONCEPT_ID
#' @details The model fitted is:
#' lab_value_i,t = β0 + β1*sex_i + β2*age_i,t + γ0i + γ1i*age_i,t + ε_i,t
#' where γ0i and γ1i are random intercept and slope for individual i
#' @export
#' @import lme4
#' @importFrom dplyr %>% group_by filter mutate select left_join distinct
calculate_blup_slopes <- function(drug_response, sex_data = NULL, output_dir = ".",
                                  min_measurements = 2) {

  if (!inherits(drug_response, "drug.reponse")) {
    stop("Input must be a drug.reponse object.")
  }

  # Check if lme4 is available
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for this function. Please install it with: install.packages('lme4')")
  }

  # Extract lab measurements
  lab_data <- drug_response$all_measurements %>%
    filter(!is.na(.data$MEASUREMENT_VALUE_HARMONIZED))

  # Get unique OMOP concept IDs
  concept_ids <- unique(lab_data$OMOP_CONCEPT_ID)

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

    # Add sex data if provided
    if (!is.null(sex_data)) {
      # Standardize sex coding to 0/1
      sex_data <- sex_data %>%
        mutate(SEX_CODED = case_when(
          toupper(.data$SEX) %in% c("M", "male", "MALE", "1") ~ 1,
          toupper(.data$SEX) %in% c("F", "female", "FEMALE", "0") ~ 0,
          TRUE ~ NA_real_
        ))

      analysis_data <- analysis_data %>%
        left_join(sex_data %>% select(.data$FINNGENID, .data$SEX_CODED),
                  by = "FINNGENID")
    } else {
      # If no sex data provided, create dummy variable
      analysis_data$SEX_CODED <- 0
    }

    # Remove individuals with missing sex if sex data was provided
    if (!is.null(sex_data)) {
      analysis_data <- analysis_data %>%
        filter(!is.na(.data$SEX_CODED))
    }

    # Check if there's enough data
    n_individuals <- length(unique(analysis_data$FINNGENID))
    if (n_individuals < 10) {
      warning(paste0("Only ", n_individuals, " individuals for OMOP_CONCEPT_ID ",
                     concept_id, ". Skipping analysis."))
      next
    }

    # Fit the linear mixed model
    # Model: lab_value ~ sex + age + (age | FINNGENID)
    # This includes random intercepts and random slopes for age
    tryCatch({
      if (!is.null(sex_data)) {
        # Model with sex as fixed effect
        lmm_model <- lme4::lmer(
          MEASUREMENT_VALUE_HARMONIZED ~ SEX_CODED + EVENT_AGE + (EVENT_AGE | FINNGENID),
          data = analysis_data,
          REML = TRUE
        )
      } else {
        # Model without sex effect
        lmm_model <- lme4::lmer(
          MEASUREMENT_VALUE_HARMONIZED ~ EVENT_AGE + (EVENT_AGE | FINNGENID),
          data = analysis_data,
          REML = TRUE
        )
      }

      # Extract BLUPs (random effects)
      random_effects <- lme4::ranef(lmm_model)$FINNGENID

      # The second column contains the random slopes (γ1i)
      blup_slopes <- data.frame(
        FINNGENID = rownames(random_effects),
        slope = random_effects[, "EVENT_AGE"],
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

      # Store results
      blup_results[[as.character(concept_id)]] <- list(
        model = lmm_model,
        blup_slopes = blup_slopes,
        n_individuals = n_individuals,
        output_file = output_file
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