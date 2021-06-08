#' Convert the reads tibble into a DESeq dataset
#'
#' @param reads reads tibble
#' @param sample_info sample info tibble with a conditions column
#' @param condition_list a string vector describing the different conditions
#' @param ds_formula a formula describing the design of the analysis
#'
#' @return DESeq dataset
#'
#' @importFrom rlang .data
#' @export
dds_from_compressed_by_annotation <- function(reads,
                                              sample_info,
                                              condition_list = "conditions",
                                              ds_formula = ~ conditions) {
  coldata <- sample_info %>% dplyr::arrange(!!(rlang::sym(condition_list[1]))) %>%
    as.matrix() %>%
    magrittr::set_rownames(dplyr::pull(dplyr::arrange(sample_info, !!(rlang::sym(condition_list[1]))), .data$sample_ID))
  cts <- reads %>%
    dplyr::select(dplyr::matches(paste(coldata[,"sample_ID"], "cutadapt$", sep ="_")),
                  dplyr::matches(paste(coldata[,"sample_ID"], "final_trimming$", sep="_"))) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix() %>%
    magrittr::set_rownames(reads$omy_miRNA)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                        colData = coldata,
                                        design = ds_formula)
  if (is.factor(sample_info$conditions)) {
    dds[["conditions"]] <- factor(dds[["conditions"]], levels = intersect(x = levels(sample_info$conditions),
                                                                          y = dds[["conditions"]]))
  }
  for (condition in condition_list) {
    if (is.factor(dplyr::pull(sample_info, condition))) {
      dds[[condition]] <- factor(dds[[condition]], levels = intersect(x = levels(dplyr::pull(sample_info, condition)),
                                                                      y = dds[[condition]]))
    }
  }
  return(dds)
}

#' Compute all the resampling combinaisons to subset biosamples
#'
#' @param sample_info sample info tibble with a conditions column
#' @param condition_name a string describing the name of the condition which is a column of sample_info
#' @param subset_size the size of the subsets to resample, if subset_proportion is set to NA
#' @param subset_proportion the proportion of the condition to resample, if set to NA the size will be constant at subset_size
#'
#' @return list of all subset_size combinaison for each condition
#'
#' @importFrom rlang .data
#' @export
compute_all_resampling_combinaisons <- function(sample_info, condition_name = "conditions", 
                                                subset_size = 4, subset_proportion = 0.66) {
  sample_info %>% 
    dplyr::select("sample_ID", condition_name) %>% 
    dplyr::group_by(dplyr::across(condition_name)) %>%
    dplyr::summarise(.data$sample_ID %>% utils::combn(dplyr::if_else(is.na(subset_proportion), 
                                                                     subset_size,
                                                                     round(subset_proportion*dplyr::n()), 
                                                                     subset_size)) %>% 
                                                        t %>% tidyr::as_tibble(.name_repair = "universal")) %>%
    tidyr::nest() %>% 
    dplyr::mutate(data = purrr::map(data, function(x) {select(x, where(~ !all(is.na(.x))))})) %>%
    tibble::deframe()
}

#' Run the analysis function for n_resample random subset of biosamples
#'
#' @param reads reads tibble
#' @param sample_info sample info tibble with a conditions column
#' @param condition_list a string vector describing the different conditions
#' @param ds_formula a formula describing the design of the analysis
#' @param analysis_function custom analysis function using the reads, sample_info, condition_list and ds_formula
#' @param control_combinaisons a tibble with control biosample subsets to analyse 
#'     (a line is a subset, a column is a biosample, column names doesnt matter)
#' @param treatment_combinaisons a tibble with treatment biosample subsets to analyse 
#'     (a line is a subset, a column is a biosample, column names doesnt matter)
#' @param n_resample the number of resampling to analyse
#'
#' @importFrom rlang .data
#' @export
compute_resampling <- function(reads, sample_info, condition_list, ds_formula, analysis_function, 
                               control_combinaisons, treatment_combinaisons, n_resample = 50) {
  tidyr::crossing(
      dplyr::slice_sample(control_combinaisons, n = 10*n_resample), 
      dplyr::slice_sample(treatment_combinaisons, n = 10*n_resample), .name_repair = "universal") %>% 
    dplyr::sample_n(n_resample) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    tidyr::pivot_longer(-.data$id) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::summarise(sample_list = list(.data$value)) %>%
    dplyr::mutate(analysis = purrr::map(.data$sample_list,
                                        ~ analysis_function(subset_samples_compressed_by_annotation(reads, .x),
                                                            subset_samples_sample_info(sample_info, .x),
                                                            condition_list, ds_formula)))
}
