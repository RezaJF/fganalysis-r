#' @title Get lab measurements from FinnGen data
#' @param all_labs data frame with lab measurements
#' @param lablist vector of lab measurement concept IDs
#' @param require_values logical, if TRUE, only return rows with non-missing MEASUREMENT_VALUE_HARMONIZED
#' @param return_cols vector of column names to return
#' @param finngen_ids vector of FINNGENIDs to filter the data
#' @param lazy logical, if TRUE, return a lazy tbl object
#' @param covariates Optional data frame with covariate data (e.g., sex, age at death)
#' @param covariate_cols Optional character vector of column names to add from the covariates data frame
#' @return data frame with lab measurements
#' @export
#' @importFrom dplyr %>% filter select collect all_of left_join distinct mutate
#' @import stringr
get_lab_measurements <- function(all_labs, lablist, require_values=TRUE,
                                 return_cols=c("FINNGENID","OMOP_CONCEPT_ID", "EVENT_AGE", "MEASUREMENT_VALUE_HARMONIZED"),
                                 finngen_ids=NULL, lazy=FALSE,
                                 covariates=NULL, covariate_cols=NULL) {

  # Store the base columns from all_labs (without covariate columns)
  # Covariate columns will be added later via the join
  base_return_cols <- unique(c("OMOP_CONCEPT_ID", return_cols))

  # Only select columns that exist in all_labs at this point
  # Filter out any covariate_cols from base_return_cols if they were accidentally included
  if (!is.null(covariate_cols)) {
    base_return_cols <- setdiff(base_return_cols, covariate_cols)
    # But ensure we keep the default return_cols
    base_return_cols <- unique(c("OMOP_CONCEPT_ID", "FINNGENID", "EVENT_AGE", "MEASUREMENT_VALUE_HARMONIZED",
                                  setdiff(return_cols, covariate_cols)))
  }

  # Ensure robust matching regardless of column/vector types (character vs numeric)
  # The OMOP_CONCEPT_ID column may be stored as DECIMAL in the parquet file
  # We need to cast it to character BEFORE filtering to avoid type conversion errors
  # Convert lablist to character for consistent comparison
  lablist_chr <- as.character(lablist)

  # Cast OMOP_CONCEPT_ID to character in the database query before filtering
  # This ensures the comparison works regardless of the storage type
  labs <- all_labs %>%
    mutate(OMOP_CONCEPT_ID = as.character(.data$OMOP_CONCEPT_ID)) %>%
    select(all_of(base_return_cols)) %>%
    dplyr::filter(.data$OMOP_CONCEPT_ID %in% lablist_chr)

  if (!is.null(finngen_ids)) {
    labs <- labs %>% dplyr::filter(.data$FINNGENID %in% finngen_ids)
  }
  if (require_values) {
    labs <- labs %>% dplyr::filter(!is.na(.data$MEASUREMENT_VALUE_HARMONIZED))
  }

  # Add covariates if provided
  if (!is.null(covariates) && !is.null(covariate_cols)) {
    # Ensure FINNGENID is in the columns for the join
    cols_to_select <- unique(c("FINNGENID", covariate_cols))

    # Check that all requested columns exist in the covariates dataframe
    missing_cols <- setdiff(cols_to_select, colnames(covariates))
    if (length(missing_cols) > 0) {
      stop(paste("The following `covariate_cols` are not in the `covariates` dataframe:",
                 paste(missing_cols, collapse = ", ")))
    }

    # Select the requested columns and ensure one row per FINNGENID
    cov_data_to_join <- covariates %>%
      select(all_of(cols_to_select)) %>%
      distinct(.data$FINNGENID, .keep_all = TRUE)

    # Join covariates with lab measurements
    # If labs is a lazy table (e.g., DuckDB), allow copying the RHS to the same src
    copy_needed <- inherits(labs, "tbl_lazy") || inherits(labs, "tbl_sql")
    labs <- labs %>%
      left_join(cov_data_to_join, by = "FINNGENID", copy = copy_needed)
  }

  if (lazy) {
    labs
  } else {
    dplyr::collect(labs)
  }
}


#' @title Get drug purchases from FinnGen data
#' @param all_phenos data frame with drug purchases
#' @param druglist vector of drug ATC codes. The ATC codes are matched with the first part of the code (e.g. A01*)
#' @param finngen_ids vector of FINNGENIDs to filter the data. leave empty to get all
#' @param return_cols vector of column names to return
#' @param lazy logical, if TRUE, return a lazy tbl object
#' @return data frame with drug purchases
#' @export
#' @importFrom dplyr %>% filter select collect rename
#' @importFrom rlang sym :=
#' @import stringr
get_drug_purchases <- function(all_phenos, druglist, finngen_ids=NULL,
                               return_cols=c("FINNGENID","EVENT_AGE", ATC="CODE1", REIMB_CODE="CODE2", VNR="CODE3", N_PACKS="CODE4"),
                               lazy=FALSE) {

  drugs_regex <- paste0("^(",
                        paste0(druglist, collapse = '|'),
                      ")")

    drugs <- all_phenos %>% dplyr::filter(.data$SOURCE=="PURCH" & str_detect(.data$CODE1, drugs_regex))

  # Handle column selection and renaming
  if (is.null(names(return_cols))) {
    # No renaming, just use columns as is
    available_cols <- intersect(return_cols, colnames(drugs))
    drugs <- drugs %>% select(all_of(available_cols))
  } else {
    # Handle renaming - extract actual column names and rename mapping
    actual_cols <- c()
    rename_map <- list()

    for (i in seq_along(return_cols)) {
      if (names(return_cols)[i] == "") {
        # No rename, just select
        actual_cols <- c(actual_cols, return_cols[i])
      } else {
        # Rename from value to name
        actual_cols <- c(actual_cols, return_cols[i])
        rename_map[[names(return_cols)[i]]] <- return_cols[i]
      }
    }

    # Select only available columns
    available_cols <- intersect(actual_cols, colnames(drugs))
    drugs <- drugs %>% select(all_of(available_cols))

    # Apply renaming for columns that exist
    for (new_name in names(rename_map)) {
      old_name <- rename_map[[new_name]]
      if (old_name %in% colnames(drugs)) {
        drugs <- drugs %>% rename(!!sym(new_name) := !!sym(old_name))
      }
    }
  }

  if (!is.null(finngen_ids)) {
    drugs <- drugs %>% dplyr::filter(.data$FINNGENID %in% finngen_ids)
  }

  if (lazy) {
    drugs
  } else {
    dplyr::collect(drugs)
  }
}


#' @title Get first drug purchase from FinnGen data
#' @param all_phenos data frame with drug purchases
#' @param druglist vector of drug ATC codes. The ATC codes are matched with the first part of the code (e.g. A01*)
#' @param finngen_ids vector of FINNGENIDs to filter the data. leave empty to get all
#' @param return_cols vector of column names to return
#' @param lazy logical, if TRUE, return a lazy tbl object
#' @return data frame with first drug purchases for each FINNGENID
#' @export
#' @importFrom dplyr %>% group_by filter distinct ungroup select collect
get_first_purchase <- function(all_phenos, druglist, finngen_ids=NULL,
                               return_cols=c("FINNGENID","EVENT_AGE","CODE1"),
                               lazy=FALSE) {

  first_purch <- get_drug_purchases(all_phenos, druglist, finngen_ids, return_cols, lazy=TRUE) %>%
    group_by(.data$FINNGENID) %>%
    filter(.data$EVENT_AGE == min(.data$EVENT_AGE)) %>% distinct(.data$EVENT_AGE, .keep_all = TRUE) %>%
    ungroup() %>%
    select(all_of(return_cols))

  if (lazy) {
    first_purch
  } else {
    dplyr::collect(first_purch)
  }
}