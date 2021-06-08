library(dplyr)
library(tidyr)

#' Load metadata from a sample info csv file
#'
#' @param path path of the sample info file
#' @param condition_list a string vector describing the different conditions
#'
#' @return tibble containing the csv file with an extra sample_ID column
#'
#' @examples
#' \dontrun{
#' load_metadata("sample_info/eggpreserve.csv")
#' }
#'
#' @importFrom rlang .data
#' @export
load_metadata <- function(path, condition_list = "condition") {
  readr::read_csv(path) %>%
    dplyr::mutate(sample_ID = stringr::str_remove(base::basename(path), "_(cutadapt|final_trimming).fastq.gz*")) %>%
    tidyr::separate(.data$conditions, condition_list, sep = "-", remove = F)
}

#' Load data from the prost compressed_by_annotation output file
#'
#' When read counts for multiple annotations, the counts are split equaly between
#' each annotation
#'
#' @param path path of the data path
#'
#' @return tibble containing the raw and normalized count by annotated miRNA of each sample
#'
#' @examples
#' \dontrun{
#' load_prost_compressed_by_annotation("prost_11/eggpreserve_compressed_by_annotation.tsv")
#' }
#'
#' @importFrom rlang .data
#' @export
load_prost_compressed_by_annotation <- function(path) {
  readr::read_tsv(path) %>%
    # Remove irrelevant lines
    dplyr::filter(.data$MainSeqMatchesAnnotationFile != "Reverse") %>%
    # Split read counts for multiple annotations
    dplyr::mutate(MainSeqMatchesAnnotationFile =
             dplyr::if_else(.data$MainSeqMatchesAnnotationFile == "Multiple", "Splitted", .data$MainSeqMatchesAnnotationFile)) %>%
    tidyr::separate_rows(.data$omy_miRNA, sep = ",") %>%
    dplyr::ungroup() %>%
    # Compute a normalization factor for multiple annotations
    dplyr::group_by(.data$omy_miRNA) %>%
    dplyr::mutate(MainSeqMatchesAnnotationFile = factor(.data$MainSeqMatchesAnnotationFile,
                                                 levels = c("Yes", "No", "Splitted"))) %>%
    dplyr::arrange(.data$MainSeqMatchesAnnotationFile, .by_group = TRUE) %>%
    dplyr::mutate_if(is.numeric, ~ tibble(count = .x, norm_factor = rep(.x[1], dplyr::n()))) %>%
    dplyr::group_by(.data$Anno_idx) %>%
    dplyr::mutate_if(tibble::is_tibble, ~ tibble(count = .x$count, norm_factor = .x$norm_factor / sum(.x$norm_factor))) %>%
    dplyr::mutate_if(tibble::is_tibble, ~ dplyr::if_else(is.na(.x$norm_factor), 0, .x$count*.x$norm_factor)) %>%
    # Summarize multiple annotations
    dplyr::group_by(.data$omy_miRNA) %>%
    dplyr::mutate_if(is.numeric, sum) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ .x[1]))
}

#' Sum the read counts for every sample by condition
#'
#' The condition is from the conditions column in sample_info
#'
#' @param compressed_by_annotation compressed_by_annotation tibble
#' @param sample_info sample_info tibble
#' @param condition_list a string vector describing the different conditions
#' @param normalized whether to use raw or normalized from prost, default raw
#'
#' @return a tibble with a line by condition and a column by annotated miRNA
#'
#' @importFrom rlang .data
read_count_per_condition <- function(compressed_by_annotation, sample_info, condition_list = "conditions", normalized = F) {
  compressed_by_annotation %>%
    dplyr::select(dplyr::matches("^omy_miRNA$") |
             #XOR normalized
             (((where(~!normalized) | dplyr::contains("_norm")) & !(where(~!normalized) & dplyr::contains("_norm"))) &
                where(is.numeric))) %>%
    tidyr::pivot_longer(-.data$omy_miRNA) %>%
    dplyr::mutate(name = stringr::str_split(.data$name, "_", simplify = T)[,1]) %>%
    tidyr::pivot_wider(names_from = .data$omy_miRNA, values_from = .data$value) %>%
    dplyr::full_join(sample_info %>% dplyr::select(.data$sample_ID, .data$condition_list),
              by = c("name" = "sample_ID")) %>%
    dplyr::select(-.data$name) %>%
    dplyr::group_by_at(dplyr::vars(dplyr::one_of(condition_list)))
    dplyr::summarise(dplyr::across(dplyr::everything(), sum))
}

#' Sum the total read counts for each sample
#'
#' @param compressed_by_annotation compressed_by_annotation tibble
#'
#' @return a tibble with a column sample and a column total_count
#'
#' @importFrom rlang .data
#' @export
raw_read_count_per_sample <- function(compressed_by_annotation) {
  compressed_by_annotation %>%
    dplyr::select(!dplyr::contains("_norm") & where(is.numeric)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "sample_ID", values_to = "total_count") %>%
    dplyr::group_by(.data$sample_ID) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), sum)) %>%
    dplyr::arrange(.data$total_count) %>%
    dplyr::mutate(sample_ID = stringr::str_split(.data$sample_ID, "_(final|cut)", simplify = T)[,1])
}

#' List the miRNAs counted in every condition
#'
#' @param count_per_condition count_per_condition tibble
#' @param threshold min threshold to be considered as present
#'
#' @return a single column tibble with the list of miRNAs
#'
#' @importFrom rlang .data
mir_in_every_condition <- function(count_per_condition, threshold = 10) {
  count_per_condition %>%
    tidyr::pivot_longer(-.data$conditions, names_to = "omy_miRNA", values_to = "count") %>%
    dplyr::filter(.data$count > threshold) %>%
    tidyr::pivot_wider(names_from = .data$conditions, values_from = .data$count) %>%
    tidyr::drop_na() %>%
    dplyr::select(.data$omy_miRNA)
}

#' List the miRNAs counted in at least one condition
#'
#' @param count_per_condition count_per_condition tibble
#' @param threshold min threshold to be considered as present
#'
#' @return a single column tibble with the list of miRNAs
#'
#' @importFrom rlang .data
mir_in_any_condition <- function(count_per_condition, threshold = 10) {
  count_per_condition %>%
    tidyr::pivot_longer(-.data$condition, names_to = "omy_miRNA", values_to = "count") %>%
    dplyr::filter(.data$count > threshold) %>%
    dplyr::select(.data$omy_miRNA) %>%
    dplyr::distinct()
}

#' Select samples to keep only samples with a minimum number of reads relative to
#' every samples.
#'
#' The samples kept have at least a log total count superior to the mean log count
#' minus a user defined numerb of sigma (default 3)
#'
#' @param compressed_by_annotation compressed_by_annotation tibble
#' @param sample_info sample_info tibble
#' @param number_of_sigma number of standard error inferior to the mean a log total count can be before being filtered out
#'
#' @importFrom rlang .data
#' @export
select_samples <- function(compressed_by_annotation, sample_info, number_of_sigma = 3) {
compressed_by_annotation %>%
  raw_read_count_per_sample() %>%
  dplyr::right_join(sample_info, by = "sample_ID") %>%
  dplyr::select(.data$sample_ID, .data$total_count, .data$conditions) %>%
  dplyr::mutate(log_count = log2(.data$total_count)) %>%
  dplyr::mutate(keep = .data$log_count > mean(.data$log_count) - number_of_sigma*stats::sd(.data$log_count)) %>%
  dplyr::filter(.data$keep) %>%
  dplyr::pull(.data$sample_ID)
}

#' Subset the sample_info tibble with a sample list
#'
#' @param sample_info sample_info tibble
#' @param sample_list list of sample to keep
#'
#' @importFrom rlang .data
#' @export
subset_samples_sample_info <- function(sample_info, sample_list) {
  sample_info %>% dplyr::filter(.data$sample_ID %in% sample_list)
}

#' Subset the compressed_by_annotation tibble with a sample list
#'
#' @param compressed_by_annotation compressed_by_annotation tibble
#' @param sample_list list of sample to keep
#' @export
subset_samples_compressed_by_annotation <- function(compressed_by_annotation, sample_list) {
  compressed_by_annotation %>%
    dplyr::select(!where(is.numeric) | 
                    dplyr::matches(paste0(sample_list, "_cutadapt")) | 
                    dplyr::matches(paste0(sample_list, "_final_trimming")))
}

#' Subset the compressed_by_annotation tibble to keep only mir with a minimum rpm
#'
#' @param compressed_by_annotation compressed_by_annotation tibble
#' @param sample_info sample_info tibble
#' @param condition_col string, name of the condition column in sample_info
#' @param min_rpm numeric, minimum rpm required for at least some samples
#' @param min_total_count numeric, minimum total count required
#' @param min_prop numeric, minimum proportion of samples in the smallest group
#'
#' @importFrom rlang .data
#' @export
select_mir_compressed_by_annotation <- function(compressed_by_annotation, sample_info, condition_col = "condition",
                                                min_rpm = 10, min_total_count = 0, min_prop = 0.7) {
  n_samples <- sample_info %>%
    dplyr::select(condition_col) %>%
    dplyr::group_by_all() %>%
    dplyr::summarise(tot = dplyr::n()) %>%
    dplyr::pull(.data$tot) %>% min

  mir_to_keep <- compressed_by_annotation %>%
    dplyr::select("omy_miRNA", 
                  dplyr::matches(paste(sample_info$sample_ID, "cutadapt", sep = "_")), 
                  dplyr::matches(paste(sample_info$sample_ID, "final_trimming", sep = "_"))
                  ) %>%
    tidyr::pivot_longer(-.data$omy_miRNA) %>%
    dplyr::group_by(.data$name) %>%
    dplyr::mutate(rpm = .data$value/sum(.data$value)*1e6) %>%
    dplyr::group_by(.data$omy_miRNA) %>%
    dplyr::filter(sum(.data$value) > min_total_count) %>%
    dplyr::filter(sort(.data$rpm, decreasing = T)[round(min_prop*n_samples)] > min_rpm) %>%
    dplyr::pull(.data$omy_miRNA) %>% unique

  compressed_by_annotation %>% dplyr::filter(.data$omy_miRNA %in% mir_to_keep)
}

#' Load data from the prost compressed_by_location output file
#'
#' @param path path to the file
#'
#' @importFrom rlang .data
load_prost_compressed_by_genomic_location_miRNA_only <- function(path) {
  columns_to_load <- readr::cols(
    Loc_idx = readr::col_character(),
    BinStarter = readr::col_character(),
    Locations = readr::col_character(),
    CIGARs_5pto3p = readr::col_character(),
    omy_miRNA = readr::col_character(),
    `Gapped?` = readr::col_logical()
  )
  readr::read_tsv(path, col_type = columns_to_load) %>%
    dplyr::filter(!is.na(.data$omy_miRNA)) %>%
    tidyr::separate_rows(.data$omy_miRNA, sep = ",")
}


#' Get location information for specific mirs
#'
#' @param compressed_by_genomic_location compressed_by_genomic_location tibble
#' @param sample_info sample_info tibble
#' @param mir_list vector of micro RNA to be searched
#'
#' @importFrom rlang .data
#' @export
mir_locations <- function(compressed_by_genomic_location, sample_info, mir_list) {
  compressed_by_genomic_location %>%
    dplyr::select("omy_miRNA", "Loc_idx","Locations", "CIGARs_5pto3p",
                  dplyr::matches(paste(sample_info$sample_ID, "cutadapt", sep = "_")), 
                  dplyr::matches(paste(sample_info$sample_ID, "final_trimming", sep = "_"))) %>%
    dplyr::filter(.data$omy_miRNA %in% mir_list) %>%
    tidyr::pivot_longer(where(is.numeric)) %>%
    dplyr::group_by(.data$Loc_idx, .data$Locations, .data$CIGARs_5pto3p, .data$omy_miRNA) %>%
    dplyr::summarise(total = sum(.data$value)) %>% dplyr::arrange(.data$omy_miRNA, .data$total)
}
