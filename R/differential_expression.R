#' Convert the compressed_by_annotation tibble into a DESeq dataset
#'
#' @param compressed_by_annotation compressed_by_annotation tibble
#' @param sample_info sample info tibble with a condition column
#' @param condition_list a string vector describing the different conditions
#' @param ds_formula a formula describing the design of the analysis
#'
#' @return DESeq dataset
#'
#' @importFrom rlang .data
#' @export
dds_from_compressed_by_annotation <- function(compressed_by_annotation,
                                              sample_info,
                                              condition_list = "conditions",
                                              ds_formula = ~ conditions) {
  coldata <- sample_info %>% dplyr::arrange(!!(rlang::sym(condition_list[1]))) %>%
    as.matrix() %>%
    magrittr::set_rownames(dplyr::pull(dplyr::arrange(sample_info, !!(rlang::sym(condition_list[1]))), .data$sample_ID))
  cts <- compressed_by_annotation %>%
    dplyr::select(paste(coldata[,"sample_ID"], "cutadapt", sep ="_")) %>%
    dplyr::mutate_all(as.integer) %>%
    as.matrix() %>%
    magrittr::set_rownames(compressed_by_annotation$omy_miRNA)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                        colData = coldata,
                                        design = ds_formula)
  if (is.factor(sample_info$conditions)) {
    dds[["conditions"]] <- factor(dds[["conditions"]], levels = levels(sample_info$conditions))
  }
  for (condition in condition_list) {
    if (is.factor(dplyr::pull(sample_info, condition))) {
      dds[[condition]] <- factor(dds[[condition]], levels = levels(dplyr::pull(sample_info, condition)))
    }
  }
  return(dds)
}
