#' @title ATC Code Mapping Functions
#' @description Functions to handle ATC code alterations and mappings over time
#' @import dplyr
#' @importFrom rjson fromJSON
NULL

# Package-level cache for ATC mappings
.atc_cache <- new.env(parent = emptyenv())


#' @title Load ATC Mappings
#' @description Loads ATC mappings from packaged JSON file
#' @param custom_file Optional path to a custom mapping file (for testing or updates)
#' @return List of ATC mappings
#' @export
load_atc_mappings <- function(custom_file = NULL) {

  # Check if already in memory cache
  if (exists("mappings", envir = .atc_cache)) {
    return(get("mappings", envir = .atc_cache))
  }

  # Determine which file to load
  if (!is.null(custom_file) && file.exists(custom_file)) {
    mapping_file <- custom_file
    message(sprintf("Loading ATC mappings from custom file: %s", custom_file))
  } else {
    # Try multiple locations in order of preference
    # Optionally check for environment variable
    env_mapping_file <- Sys.getenv("FGANALYSIS_ATC_MAPPING_FILE", unset = NA)
    env_file <- if (!is.na(env_mapping_file) && nzchar(env_mapping_file)) env_mapping_file else NULL
    possible_files <- c(
      "config/atc_mappings.json",  # Local development path
      env_file,  # Environment variable
      system.file("config", "atc_mappings.json", package = "fganalysis")  # Installed package
    )

    mapping_file <- NULL
    for (file in possible_files) {
      if (!is.null(file) && file.exists(file)) {
        mapping_file <- file
        break
      }
    }

    if (is.null(mapping_file)) {
      warning("ATC mappings file not found. Using empty mappings.")
      mappings <- list(
        mappings = list(),
        version = "0.0.0",
        generated_date = Sys.Date(),
        source_url = "none",
        total_mappings = 0
      )
      assign("mappings", mappings, envir = .atc_cache)
      return(mappings)
    }
  }

  # Load JSON file
  tryCatch({
    json_content <- readLines(mapping_file, warn = FALSE)
    mappings <- rjson::fromJSON(paste(json_content, collapse = ""))

    # Store in memory cache
    assign("mappings", mappings, envir = .atc_cache)

    if (!is.null(mappings$total_mappings) && mappings$total_mappings > 0) {
      message(sprintf("Loaded %d ATC mappings (version %s, generated %s)",
                     mappings$total_mappings,
                     mappings$version,
                     mappings$generated_date))
    }

    return(mappings)

  }, error = function(e) {
    warning(sprintf("Error loading ATC mappings: %s", e$message))
    # Return empty mappings as fallback
    mappings <- list(
      mappings = list(),
      version = "0.0.0",
      generated_date = Sys.Date(),
      source_url = "none",
      total_mappings = 0
    )
    assign("mappings", mappings, envir = .atc_cache)
    return(mappings)
  })
}

#' @title Expand ATC Codes with Historical and Related Codes
#' @description Takes a vector of ATC codes and expands them to include all historical and related codes
#' @param atc_codes Character vector of ATC codes to expand
#' @param include_hierarchical If TRUE, also includes hierarchical relationships (e.g., A10A includes A10AA, A10AB)
#' @param verbose If TRUE, prints detailed expansion information
#' @param custom_mapping_file Optional path to custom mapping file
#' @return Character vector of expanded ATC codes
#' @export
expand_atc_codes <- function(atc_codes, include_hierarchical = TRUE, verbose = TRUE, custom_mapping_file = NULL) {

  # Load mappings
  mappings_data <- load_atc_mappings(custom_file = custom_mapping_file)
  mappings <- mappings_data$mappings

  # Initialize expanded set with original codes
  expanded_codes <- character(0)

  if (verbose) {
    message("=== ATC Code Expansion ===")
    message(sprintf("Expanding %d input ATC code(s) using historical mappings", length(atc_codes)))
    message(sprintf("Mappings database contains %d entries (last updated: %s)",
                   length(mappings), mappings_data$generated_date))
  }

  # Process each input code
  for (code in atc_codes) {
    code_expanded <- c(code)  # Always include the original code

    # Check if this code has mappings
    if (code %in% names(mappings)) {
      mapping <- mappings[[code]]

      # Add historical codes (handle both list and vector formats)
      historical <- if (is.list(mapping$historical)) unlist(mapping$historical) else mapping$historical
      if (length(historical) > 0) {
        code_expanded <- c(code_expanded, historical)
      }

      # Add related/new codes (handle both list and vector formats)
      related <- if (is.list(mapping$related)) unlist(mapping$related) else mapping$related
      if (length(related) > 0) {
        code_expanded <- c(code_expanded, related)
      }

      # Add current code if different
      if (!is.null(mapping$current) && mapping$current != code) {
        code_expanded <- c(code_expanded, mapping$current)
      }

      if (verbose && length(code_expanded) > 1) {
        message(sprintf("  %s -> expanded to %d codes: [%s]",
                       code, length(code_expanded),
                       paste(code_expanded, collapse = ", ")))
      }
    } else if (verbose) {
      message(sprintf("  %s -> no alterations found (using original code only)", code))
    }

    # Handle hierarchical expansion if requested
    if (include_hierarchical) {
      # For each expanded code, also include partial matches
      # e.g., if we have A10AA01, also search for A10AA, A10A, A10
      hierarchical_codes <- character(0)
      for (exp_code in code_expanded) {
        # Generate hierarchical codes by progressively removing characters
        for (i in seq(nchar(exp_code) - 1, 3, -1)) {
          hierarchical_codes <- c(hierarchical_codes, substr(exp_code, 1, i))
        }
      }
      code_expanded <- unique(c(code_expanded, hierarchical_codes))
    }

    expanded_codes <- c(expanded_codes, code_expanded)
  }

  # Remove duplicates
  expanded_codes <- unique(expanded_codes)

  if (verbose) {
    message(sprintf("=== Expansion Complete ==="))
    message(sprintf("Total: %d input code(s) expanded to %d unique code(s)",
                   length(atc_codes), length(expanded_codes)))
    if (length(expanded_codes) > length(atc_codes)) {
      new_codes <- setdiff(expanded_codes, atc_codes)
      message(sprintf("Added codes: [%s]", paste(head(new_codes, 10), collapse = ", ")))
      if (length(new_codes) > 10) {
        message(sprintf("  ... and %d more", length(new_codes) - 10))
      }
    }
  }

  return(expanded_codes)
}

#' @title Get All Related ATC Codes
#' @description Returns all codes related to a given ATC code (historical, current, and related)
#' @param atc_code Single ATC code to query
#' @param mappings Optional pre-loaded mappings (if NULL, will load from cache)
#' @return List with detailed mapping information
#' @export
get_atc_relationships <- function(atc_code, mappings = NULL) {

  if (is.null(mappings)) {
    mappings_data <- load_atc_mappings()
    mappings <- mappings_data$mappings
  }

  if (atc_code %in% names(mappings)) {
    return(mappings[[atc_code]])
  } else {
    return(list(
      original = atc_code,
      current = atc_code,
      historical = character(0),
      related = character(0),
      alteration_years = character(0),
      message = "No alterations found for this code"
    ))
  }
}

#' @title Clear ATC Mapping Cache
#' @description Clears the in-memory cache of ATC mappings
#' @export
clear_atc_cache <- function() {
  if (exists("mappings", envir = .atc_cache)) {
    rm("mappings", envir = .atc_cache)
    message("ATC mapping cache cleared")
  }
}
