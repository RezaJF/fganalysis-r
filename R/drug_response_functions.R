#' @import dplyr
#' @import ggplot2
#' @import ggpubr
#' @import stringr
#' @importFrom stats median lm quantile sd
#' @import utils
#' @import grDevices
#' @importFrom checkmate assertClass assertNumeric assertCharacter assertLogical
# NULL

#' @title data object returned from drug response analyse (create_drug_response)
#' @param responses data frame with response data
#' @param lab_measurements data frame with all lab measurements
#' @param drug_purchases data frame with all drug purchases
#' @param before_period vector with two elements, start and end of the before period
#' @param after_period vector with two elements, start and end of the after period
#' @return object of class drug.response
#' @export
drug.response <- function(responses, lab_measurements, drug_purchases, before_period, after_period) {
    obj <- list(
        responses = responses, all_measurements = lab_measurements, all_drug_purchases = drug_purchases,
        lab_response_period = list(before_period = before_period, after_period = after_period)
    )
    class(obj) <- "drug.reponse"
    return(obj)
}

#' @title Create drug response object
#' @param conn fg_data_connection object
#' @param lablist vector of lab measurement concept IDs
#' @param druglist vector of drug ATC codes. The ATC codes are matched with the first part of the code (e.g. A01*).
#' @param before_period vector with two elements, start and end of the before period
#' @param after_period vector with two elements, start and end of the after period
#' @param summary_function function to summarize the lab measurements (default is median, use mean for mean etc.)
#' @param filter_min_max vector of length 2, min and max values to filter lab measurements (default c(-Inf, Inf))
#' @param finngen_ids vector of FINNGENIDs to filter the data
#' @param use_lab_free_text_values logical, if TRUE, use lab measurements combining reported values and those from free text values (default TRUE)
#' @param use_only_reimbursement_drugs logical, if TRUE, use only reimbursement data (i.e. longitudinal data from purchase registry) (default FALSE)
#' @return drug.response object
#' @export
create_drug_response <- function(
    conn, lablist, druglist,
    before_period, after_period, summary_function = median, filter_min_max = c(-Inf, Inf), finngen_ids = NULL,
    use_lab_free_text_values = TRUE, use_only_reimbursement_drugs = FALSE) {
    if (!inherits(conn, "fg_data_connection")) {
        stop("conn must be a fg_data_connection object")
    }

    assertCharacter(lablist, any.missing = FALSE)
    assertCharacter(druglist, any.missing = FALSE)
    assertNumeric(filter_min_max, len = 2, any.missing = FALSE)
    assertNumeric(before_period, len = 2, any.missing = FALSE)
    assertNumeric(after_period, len = 2, any.missing = FALSE)
    if (!inherits(conn, "fg_data_connection")) {
        stop("conn must be a fg_data_connection object")
    }
    assertLogical(use_lab_free_text_values)
    assertLogical(use_only_reimbursement_drugs)

    # assert_class( filter_function, "dplyr::filter" )

    print("Querying lab measurements...")
    lab_measurements <- get_lab_measurements(
        all_labs = conn$labs, lablist = lablist, finngen_ids = finngen_ids,
        require_values = TRUE, use_freetext_values = use_lab_free_text_values)

    orig <- nrow(lab_measurements)
    lab_measurements <- lab_measurements %>% filter(.data$VALUE >= filter_min_max[1] &
        .data$VALUE <= filter_min_max[2])

    print(paste0("Filtered ", orig - nrow(lab_measurements), " measurements due to min and max filters"))

    all_fg_ids <- unique(c(lab_measurements$FINNGENID), finngen_ids)
    print("Querying purchases...")
    drug_purchases <- get_drug_purchases(conn, druglist, all_fg_ids, use_only_reimbursement = use_only_reimbursement_drugs)
    print(paste0("Number of drug purchases: ", nrow(drug_purchases)))
    dr_first_purchase <- drug_purchases %>%
        group_by(.data$FINNGENID) %>%
        arrange(.data$EVENT_AGE) %>%
        summarize(n = n(), first_drug_age = first(.data$EVENT_AGE), first_drug = first(.data$ATC))

    lab_measurements <- left_join(lab_measurements, dr_first_purchase, by = "FINNGENID")
    lab_measurements <- lab_measurements %>% mutate(time_to_drug = .data$first_drug_age - .data$EVENT_AGE)

    print("generating response summary...")

    lab_response <- generate_response_summary(lab_measurements, before_period, after_period,
        summary_function = summary_function
    )
    print(paste0("Number of individuals with response data: ", nrow(lab_response)))

    return(drug.response(
        responses = lab_response, lab_measurements = lab_measurements,
        drug_purchases = drug_purchases, before_period, after_period
    ))
}

quant_text <- function(vector) {
    paste0(paste(names(quantile(vector)), sep = "\t"), ":", paste(quantile(vector), sep = "\t"), collapse = " ")
}


#' @title Summarize drug response
#' @description Summarize drug response data created with create_drug_response. writes plots and tables to disk
#' @param drug_response drug.response object
#' @param out_file_prefix prefix for output files
#' @return NULL
#' @export
summarize_drug_response <- function(drug_response, out_file_prefix) {
    labs <- drug_response$all_measurements %>% filter(!is.na(.data$VALUE))
    responses <- drug_response$responses
    drugs <- drug_response$all_drug_purchases

    pdf(paste0(out_file_prefix, ".pdf"), width = 10, height = 6)

    inds_with_lab <- length(unique(labs$FINNGENID))
    n_lab_meas <- nrow(labs)

    inds_with_drugs <- length(unique(drugs$FINNGENID))
    n_drugs_meas <- nrow(drugs)

    inds_in_analysis <- nrow(responses)
    n_range_before <- quant_text(responses$n_before)
    n_range_after <- quant_text(responses$n_before)

    range_baseline_age <- quant_text(responses$n_before)

    n_no_pre <- drug_response$responses %>% filter(is.na(response) & n_before == 0)
    n_no_pos <- drug_response$responses %>% filter(is.na(response) & n_after == 0)

    plot(ggtexttable(data.frame(
        "Group" = c("labs", "drugs", "in analysis", "no_pre_value", "no_post_value"),
        "N" = c(inds_with_lab, inds_with_drugs, inds_in_analysis, nrow(n_no_pre), nrow(n_no_pos)),
        "N events" = c(n_lab_meas, n_drugs_meas, inds_in_analysis, sum(n_no_pre$n_after), sum(n_no_pos$n_before))
    ), rows = NULL))


    responses <- responses %>% filter(!is.na(response))
    per_drug <- responses %>%
        group_by(.data$first_drug) %>%
        summarise(
            N = n(), p = summary(lm("response ~ 1", data = pick(.data$FINNGENID, .data$response)))$coefficients[1, 4], sd = sd(.data$response),
            response = mean(.data$response),
            purch_age_dist = quant_text(.data$baseline_age)
        )

    all_resp <- rbind(per_drug, data.frame(
        first_drug = "All drugs", N = inds_in_analysis,
        response = mean(responses$response), p = summary(lm("response ~ 1", data = responses))$coefficients[1, 4],
        purch_age_dist = quant_text(responses$baseline_age), sd = sd(responses$response)
    ))

    write.table(
        all_resp %>% arrange(desc(.data$N)) %>%
            select(.data$first_drug, .data$N, .data$response, .data$p, .data$purch_age_dist),
        paste0(out_file_prefix, "_responses_by_drug.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )


    plot(ggtexttable(responses %>% group_by(.data$first_drug) %>%
        summarise(
            n_purch = n(), n_indiv = length(unique(.data$FINNGENID)),
            p = summary(lm("response ~ 1", data = pick(.data$FINNGENID, .data$response)))$coefficients[1, 4],
            response = mean(.data$response),
            purch_age_dist = quant_text(.data$baseline_age)
        ) %>%
        select(.data$first_drug, .data$n_purch, .data$response, .data$p, .data$purch_age_dist) %>%
        arrange(desc(.data$n_purch))))


    begin <- ceiling(max(min(-labs$time_to_drug, na.rm = TRUE), drug_response$lab_response_period$before_period[1]))
    end <- ceiling(min(max(-labs$time_to_drug, na.rm = TRUE), drug_response$lab_response_period$after_period[2]))
    labs$bin <- cut(-labs$time_to_drug,
        breaks = seq(begin, end, by = .25),
        include.lowest = TRUE
    )

    plot(ggplot(labs %>% filter(!is.na(.data$bin))) +
        geom_boxplot(aes(x = .data$bin, y = .data$VALUE)) +
        labs(x = "Time to drug purchase (years)", y = "Lab measurement") +
        ggtitle("Lab measurements before and after drug purchase")) + theme_bw() +
        theme(axis.text.x = element_text(angle = 45, size = 20)) 


    uniq_drugs <- unique(responses$first_drug)

    for( drug in uniq_drugs){
        labs_sub <- labs %>% filter(!is.na(.data$bin))  
        p <- ggplot(labs_sub)  + geom_boxplot(aes(x = .data$bin, y = .data$VALUE)) +
        labs(x = "Time to drug purchase (years)", y = "Lab measurement") +
        ggtitle(paste("Lab measurements before and after drug purchase for drug ", drug)) + theme_bw() +
        theme(axis.text.x = element_text(size = 10))
        plot(p)
    }

    write.table(
        labs %>% group_by(.data$bin) %>%
            summarise(n = n(), mean = mean(.data$VALUE), sd = sd(.data$VALUE)),
        paste0(out_file_prefix, "_labs_by_time_to_drug.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )


    plot(ggplot(responses) +
        geom_histogram(aes(x = .data$response)) +
        theme_bw() +
        labs(x = "Response (after - before)", y = "Count") +
        ggtitle("Distribution of drug response"))

    
    plot(ggplot(responses) +
        geom_histogram(aes(x = .data$response_percent)) +
        theme_bw() +
        labs(x = "Response%", y = "Count") +
        ggtitle("Distribution of % drug response"))

    
    for( drug in uniq_drugs){
        labs_sub <- responses %>% filter( .data$first_drug == drug)
        p <- ggplot(labs_sub) +
        geom_histogram(aes(x = .data$response)) +
        theme_bw() +
        labs(x = "Response (after - before)", y = "Count") +
        ggtitle( paste0("Distribution of drug response for drug ", drug))   
        plot(p)

        p <- ggplot(labs_sub) +
        geom_histogram(aes(x = .data$response_percent)) +
        theme_bw() +
        labs(x = "Response%", y = "Count") +
        ggtitle(paste0("Distribution of % drug response for drug ", drug))
        plot(p)

    }



    fit <- summary(lm(after ~ before, data = responses))

    slope <- format(fit$coefficients["before", "Estimate"], digits = 2)
    r2 <- format(fit$r.squared, digits = 2)
    p <- format(fit$coefficients["before", "Pr(>|t|)"], digits = 2, scientific = TRUE)

    suppressWarnings(print(ggplot(responses, aes(x = .data$before, y = .data$after)) +
        geom_point() +
        geom_smooth(method = "lm") +
        geom_abline(slope = 1) +
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
generate_response_summary <- function(lab_measurements, before_period, after_period, summary_function = median) {
    
    
    lab_measurements <- lab_measurements %>% mutate(lab_period = case_when(
        dplyr::between(.data$time_to_drug,  -before_period[2],-before_period[1]) ~ "Before",
        dplyr::between(.data$time_to_drug, -after_period[2], -after_period[1]) ~ "After",
        TRUE ~ NA_character_
    ))

    lab_response <- lab_measurements %>%
        dplyr::filter(!is.na(.data$lab_period) & !is.na(.data$VALUE)) %>%
        dplyr::group_by(.data$FINNGENID) %>%
        dplyr::summarize(
            n_before = length(.data$VALUE[.data$lab_period == "Before"]),
            n_after = length(.data$VALUE[.data$lab_period == "After"]),
            before = ifelse(n_before>0,summary_function(.data$VALUE[.data$lab_period == "Before"], na.rm = TRUE),NA),
            after = ifelse(n_after>0,summary_function(.data$VALUE[.data$lab_period == "After"], na.rm = TRUE),NA),
            baseline_age = first(.data$first_drug_age),
            first_drug = first(.data$first_drug),
            response = ifelse(!is.na(.data$after) & !is.na(.data$before), .data$after - .data$before, NA),
            response_percent = response/ .data$before * 100
        )
    # %>% dplyr::filter(!is.na(.data$response))

    return(lab_response)
}

#' @title Get lab measurements from FinnGen data
#' @param all_labs data frame with lab measurements
#' @param lablist vector of lab measurement concept IDs
#' @param require_values logical, if TRUE, only return rows with non-missing VALUE
#' @param use_freetext_values logical, if False, use only MEASUREMENT_VALUE_HARMONIZED column for lab values and not combining from reported values and extracted from free text
#' (default True)
#' @param finngen_ids vector of FINNGENIDs to filter the data
#' @param lazy logical, if TRUE, return a lazy tbl object
#' @return data frame with lab measurements, lab value will be in column named VALUE, regardless of original column name 
#' (e.g. MEASUREMENT_VALUE_HARMONIZED or MEASUREMENT_VALUE_MERGED)
#' @export
get_lab_measurements <- function(all_labs, lablist, require_values = TRUE, use_freetext_values = TRUE,
                                 finngen_ids = NULL, lazy = FALSE) {

    if(! use_freetext_values){
        lab_value_column <- "MEASUREMENT_VALUE_HARMONIZED"
    }
    else {
        lab_value_column <- "MEASUREMENT_VALUE_MERGED"
    }

    return_cols <- c("OMOP_CONCEPT_ID","FINNGENID", "OMOP_CONCEPT_ID", "EVENT_AGE", "VALUE"=lab_value_column)

    labs <- all_labs %>%
        select(all_of(return_cols)) %>%
        dplyr::filter(.data$OMOP_CONCEPT_ID %in% lablist)

    if (!is.null(finngen_ids)) {
        labs <- labs %>% dplyr::filter(.data$FINNGENID %in% finngen_ids)
    }
    if (require_values) {
        labs <- labs %>% dplyr::filter(!is.na(.data$VALUE))
    }

    ifelse(lazy, return(labs), return(dplyr::collect(labs)))
}


#' @title Get drug purchases from FinnGen data
#' @param conn finngen data connection object
#' @param druglist vector of drug ATC codes. The ATC codes are matched with the first part of the code (e.g. A01*)
#' @param finngen_ids vector of FINNGENIDs to filter the data. leave empty to get all
#' @param use_only_reimbursement logical, if TRUE, use only reimbursement data (default FALSE) and combine reimbursement and delivery data if available in conn
#' @param lazy logical, if TRUE, return a lazy tbl object
#' @return data frame with drug purchases
#' @export
get_drug_purchases <- function(conn, druglist, finngen_ids = NULL, use_only_reimbursement = FALSE,
                               lazy = FALSE) {
    ## check that conn is a fg_data_connection object and has pheno data name
    if (!inherits(conn, "fg_data_connection")) {
        stop("conn must be a fg_data_connection object")
    }

    if (!"pheno" %in% names(conn)) {
        stop("conn must contain 'pheno' data")
    }

    drugs_regex <- paste0("^", paste(druglist, collapse = "|^"), collapse = "")

    all_phenos <- conn$pheno

    if ("drug_events" %in% names(conn) & !use_only_reimbursement) {
        print("Using drug data combining reimbursement and delivery data.")
        drugs <- conn$drug_events %>% 
                dplyr::filter( str_detect( .data$ATC, drugs_regex ) ) 
    } else {
        print("Using drug data from reimbursement data only! i.e. longitudinal data from purchase registry.")
        return_cols <- c("FINNGENID", "EVENT_AGE", "APPROX_EVENT_DAY", ATC = "CODE1", REIMB_CODE = "CODE2", VNR = "CODE3", N_PACKS = "CODE4")
        drugs <- all_phenos %>%
            dplyr::filter(.data$SOURCE == "PURCH" & str_detect( .data$CODE1, drugs_regex )) %>%
            select(all_of(return_cols))
    }

    if (!is.null(finngen_ids)) {
        drugs <- drugs %>% dplyr::filter(.data$FINNGENID %in% finngen_ids)
    }

    if ("vnr" %in% names(conn)) {
        columns <- c("VNR", "Substance", "MedicineName", "PackageSize", "DDDPerPack", "Dosage", "DosageUnit")
        vnr <- conn$vnr %>% select(all_of(columns))
        drugs <- left_join(drugs, vnr, by = "VNR", copy = TRUE)
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
get_first_purchase <- function(all_phenos, druglist, finngen_ids = NULL, return_cols = c("FINNGENID", "EVENT_AGE", "CODE1"),
                               lazy = FALSE) {
    first_purch <- get_drug_purchases(all_phenos, druglist, finngen_ids, return_cols, lazy = TRUE) %>%
        group_by(.data$FINNGENID) %>%
        filter(.data$EVENT_AGE == min(.data$EVENT_AGE)) %>%
        distinct(.data$EVENT_AGE, .keep_all = TRUE) %>%
        ungroup() %>%
        select(all_of(return_cols))
}
