#' @import DBI
#' @import duckdb
#' @import duckdbfs
#' @import bigrquery
#' @import rjson
#' @importFrom data.table fread
#NULL
get_bigquery_dbplyr <- function(projectid,dataset, table) {
   con <- dbConnect(
    bigrquery::bigquery(),
    project = projectid,
    dataset=dataset
  )
  return(tbl(con, table))
}

get_parquet_dbplyr <- function(parquet_file) {
  conn = dbConnect(duckdb::duckdb())
  return(duckdb::tbl_file(conn, parquet_file))
}

get_duckdb_dbplyr <- function(duckdb_file, table) {
  conn = dbConnect(duckdb::duckdb(duckdb_file), read_only=TRUE)
  return(tbl(conn, table))
}

fg_data_connection <- function(connections) {
    if (!is.list(connections)) {
        stop("Connections must be a list")
    }
    if (!all(c("pheno", "labs") %in% names(connections))) {
        stop("Connections list must contain 'pheno' and 'labs'")
    }
    ## do further checks on the connections as needed.
    class(connections) <- "fg_data_connection"
    return(connections)
}

#' @export
print.fg_data_connection <- function(x,...) {
    cat("FinnGen data connection object\n")
    cat("Contains lazy loaded data frames of FinnGen phenotype data e.g.:\n")
    cat("pheno : longitudinal service sector data\n")
    cat("labs: kanta lab values\n")
    cat("minimum: minimum phenotypes\n")
    cat("cov_pheno: covariate phenotypes\n")
    cat("Connections are lazy loaded, use dplyr verbs to load data and use function collect() to localize results \n")
    cat("All atrributes of the class:\n")
    print(attributes(x))
}

call_connect <- function(conf) {
    req_tags <- c("path", "type")
    if (!all(req_tags %in% names(conf))) {
        stop(paste("Configuration must contain these tag [", paste(req_tags,collapse=","),"]"))
    }

    path <- conf$path
    typestring <- conf$type

    if (typestring == "parquet") {
        # Use duckdb directly to avoid extension loading issues
        dat <- tryCatch({
            # Try duckdbfs first with disabled extensions
            options(duckdb.allow_unsigned_extensions = FALSE)
            duckdbfs::open_dataset(path, format="parquet")
        }, error = function(e) {
            # Fallback to direct duckdb connection if duckdbfs fails
            message("duckdbfs failed, using direct duckdb connection")
            conn <- dbConnect(duckdb::duckdb())
            tbl(conn, paste0("read_parquet('", path, "')"))
        })
    } else if (typestring == "parquet-hive") {
        # Use duckdb directly to avoid extension loading issues
        dat <- tryCatch({
            # Try duckdbfs first with disabled extensions
            options(duckdb.allow_unsigned_extensions = FALSE)
            duckdbfs::open_dataset(path, format="parquet", hive_style=TRUE)
        }, error = function(e) {
            # Fallback to direct duckdb connection if duckdbfs fails
            message("duckdbfs failed, using direct duckdb connection for hive-style parquet")
            conn <- dbConnect(duckdb::duckdb())
            # Hive-style parquet with one level of partitioning (SOURCE=*)
            tbl(conn, paste0("read_parquet('", path, "/*/*.parquet', hive_partitioning=true)"))
        })
    } else if (typestring == "tsv")
    {
        dat <- fread(path, sep="\t")
    } else {
        stop(paste("Unsupported connection type given in configuration file:", typestring,
                   ". Supported types are: parquet, parquet-hive, tsv"))
    }

    if ("recodings" %in% names(conf)) {
        for (recoding in conf$recodings) {
            if (!("column" %in% names(recoding)) || !("function" %in% names(recoding))) {
                stop("Recoding must contain 'column' and 'function' keys")
            }

            column <- recoding$column
            function_name <- recoding[["function"]]
            print(function_name)

            dat[[column]] <- do.call(function_name, list(dat[[column]]))
        }
    }
    return(dat)
}

#' Connect to FinnGen data
#' @param path_to_conf Path to the configuration file (JSON format)
#' @return A fg.data.connection object containing the connections to data
#' @export
connect_fgdata <-  function(path_to_conf) {
    json_data <- rjson::fromJSON(file = path_to_conf)
    req_confs <- c("pheno", "labs")
    connections <- list()
    # Extract the values from the JSON object
    if (! all(req_confs %in% names(json_data))) {
        stop(paste("Elements not found in configuration data. Configuration file must contain the following keys:", paste(req_confs, collapse = ", ")))
    }
    print("loading connections")
    for (conf in names(json_data)) {
        connections[[conf]] <- call_connect(json_data[[conf]])
    }

    return(fg_data_connection( connections))
}

#' Create a mock fg_data_connection object for testing
#' @description Creates a mock connection object using data frames instead of database connections.
#' This is useful for testing, development, and when working with small datasets that can fit in memory.
#' @param pheno_data Data frame containing phenotype data (required)
#' @param labs_data Data frame containing lab measurement data (required)
#' @param minimum_data Optional data frame containing minimum phenotype data
#' @param cov_pheno_data Optional data frame containing covariate phenotype data
#' @param endpoint_data Optional data frame containing endpoint data
#' @param vnr_data Optional data frame containing VNR data
#' @param long_anthropometric_data Optional data frame containing longitudinal anthropometric data
#' @return A fg_data_connection object that can be used with package functions
#' @export
#' @examples
#' # Create minimal mock connection
#' mock_conn <- create_mock_connection(
#'   pheno_data = data.frame(FINNGENID = "FG1", EVENT_AGE = 50, CODE1 = "A01", SOURCE = "PURCH"),
#'   labs_data = data.frame(FINNGENID = "FG1", EVENT_AGE = 49, MEASUREMENT_VALUE_HARMONIZED = 100, OMOP_CONCEPT_ID = "L1")
#' )
#' 
#' # Create comprehensive mock connection
#' mock_conn_full <- create_mock_connection(
#'   pheno_data = pheno_df,
#'   labs_data = labs_df,
#'   minimum_data = minimum_df,
#'   cov_pheno_data = cov_df
#' )
create_mock_connection <- function(pheno_data, labs_data, 
                                  minimum_data = NULL, cov_pheno_data = NULL,
                                  endpoint_data = NULL, vnr_data = NULL,
                                  long_anthropometric_data = NULL) {
  
  # Validate required data
  if (is.null(pheno_data) || is.null(labs_data)) {
    stop("pheno_data and labs_data are required")
  }
  
  # Create connections list with required data
  connections <- list(
    pheno = pheno_data,
    labs = labs_data
  )
  
  # Add optional data if provided
  if (!is.null(minimum_data)) {
    connections$minimum <- minimum_data
  }
  if (!is.null(cov_pheno_data)) {
    connections$cov_pheno <- cov_pheno_data
  }
  if (!is.null(endpoint_data)) {
    connections$endpoint <- endpoint_data
  }
  if (!is.null(vnr_data)) {
    connections$vnr <- vnr_data
  }
  if (!is.null(long_anthropometric_data)) {
    connections$long_anthropometric <- long_anthropometric_data
  }
  
  return(fg_data_connection(connections))
}


