library(dplyr)
library(tidyr)
library(ggplot2)

#' Histogram of the count distribution
#'
#' The range of value is from 0 to 100,
#'
#' @param compressed_by_annotation compressed_by_annotation tibble
#' @param normalized whether to use raw or normalized from prost, default normalized
#' @param xlim A numeric vector of length two providing limits of the scale
#' @param title title to add to the plot, default : "Project"
#'
#' @return ggplot histogram
#'
#' @importFrom rlang .data
cpm_distribution_histogram <- function(compressed_by_annotation, normalized = T, xlim = c(0,100), title = "Project") {
  compressed_by_annotation %>%
    dplyr::select(where(is.numeric) &
             #XOR normalized
             ((where(~!normalized) | dplyr::contains("_norm")) & !(where(~!normalized) & dplyr::contains("_norm")))) %>%
    tidyr::pivot_longer(dplyr::everything()) %>%
    ggplot2::ggplot(ggplot2::aes(x=.data$value)) +
      ggplot2::geom_histogram(bins = 50) +
      ggplot2::scale_x_continuous(limits = xlim, oob = scales::squish) +
      ggplot2::ggtitle(paste("miRNA counts distribution", title, sep = " - "))
}

#' Vector of ggplot colors
#'
#' @param n number of colors
#' @export
  gg_color <- function(n) {
    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
  }

#' Heatmap sample to miRMA
#'
#' Warning : requires metadata to have "condition" and "project_name" entries
#'
#' @param dds Normalized DESeq dataset
#' @param deseq_transformation DESeq transformation (e.g. normTransform(dds))
#' @param cluster_cols boolean values determining if columns should be clustered or hclust object
#' @param condition_list a string vector describing the different conditions
#' @param scale character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param title title to add to the plot, default : "Project"
#'
#' @return pheatmap heatmap
#' @export
sample_to_mir_heatmap <- function(dds, deseq_transformation, cluster_cols = T, condition_list = "conditions", scale="row", title = "Project"){
  select <- order(rowMeans(DESeq2::counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:nrow(dds)]
  mat <- SummarizedExperiment::assay(deseq_transformation)[select,]
  ann <- SummarizedExperiment::colData(deseq_transformation)[,condition_list, drop = F] %>%
    as.data.frame
  colnames(mat) <- colnames(mat)
  rownames(ann) <- rownames(ann)
  ann_col <- stats::setNames(list(stats::setNames(gg_color(length(unique(ann[[1]]))), unique(ann[[1]]))),
                      condition_list[1])
  pheatmap::pheatmap(mat,
                     cluster_rows=FALSE,
                     #show_rownames=FALSE,
                     cluster_cols=cluster_cols,
                     annotation_colors = ann_col,
                     annotation_col=ann,
                     scale = scale,
                     main = paste("Sample to miRNA", title, sep = " - "))
}

#' Heatmap sample to sample
#'
#' @param dds Normalised DESeq dataset
#' @param deseq_transformation DESeq transformation (e.g. normTransform(dds))
#' @param title title to add to the plot, default : "Project"
#'
#' @return pheatmap heatmap
#' @export
sample_to_sample_heatmap <- function(dds, deseq_transformation, title = "Project"){
  sampleDists <- stats::dist(t(SummarizedExperiment::assay(deseq_transformation)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colnames(sampleDistMatrix) <- NULL
  colors <- grDevices::colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,
                     col=colors,
                     main = paste("Sample to miRNA", title, sep = " - "))
}

#' PCA plot of the samples
#'
#' @param deseq_transformation DESeq transformation (e.g. normTransform(dds))
#' @param condition_list a string vector describing the different conditions
#' @param title title to add to the plot, default : "Project"
#' @param add_theme if TRUE then a theme is added to the figure, labels, title, colors...
#' @param draw_arrows if TRUE then miR contribution arrows are added to the figure
#' @param n_arrows number of arrow to add to the figure, up to 10
#' @param draw_ellipses if TRUE then ellipses around conditions are added to the figure
#' @param draw_sample_label if TRUE then sample names are added to the figure
#'
#' @return ggplot point
#'
#' @importFrom rlang .data
#' @export
pca_plot <- function(deseq_transformation, condition_list = "conditions", title = "Project",
                     add_theme = T,
                     draw_arrows = T,
                     n_arrows = 10,
                     draw_ellipses = T,
                     draw_sample_label = F){

  res.pca = FactoMineR::PCA(graph=F,
                            X = t(SummarizedExperiment::assay(deseq_transformation)),
                            scale.unit = F)

  pcatab <- data.frame(factoextra::get_pca_ind(res.pca)$coord[,1:2]) %>%
    tidyr::as_tibble(rownames = "sample_ID") %>%
    dplyr::right_join(tidyr::as_tibble(SummarizedExperiment::colData(deseq_transformation)[,condition_list, drop = F],
                         rownames = "sample_ID"), by = "sample_ID")

  pcacontrib <- data.frame(factoextra::get_pca_var(res.pca)$coord[,1:2]) %>%
    tidyr::as_tibble(rownames = "omy_miRNA") %>%
    dplyr::arrange(dplyr::desc(.data$Dim.1**2 + .data$Dim.2**2)) %>% dplyr::slice(1:10)

  n_color <- length(unique(dplyr::pull(pcatab, condition_list[1])))

  pcacontrib_arrows <- list(
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[1], yend = 5*pcacontrib$Dim.2[1]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[1], y = 6*pcacontrib$Dim.2[1], label = pcacontrib$omy_miRNA[1]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[2], yend = 5*pcacontrib$Dim.2[2]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[2], y = 6*pcacontrib$Dim.2[2], label = pcacontrib$omy_miRNA[2]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[3], yend = 5*pcacontrib$Dim.2[3]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[3], y = 6*pcacontrib$Dim.2[3], label = pcacontrib$omy_miRNA[3]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[4], yend = 5*pcacontrib$Dim.2[4]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[4], y = 6*pcacontrib$Dim.2[4], label = pcacontrib$omy_miRNA[4]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[5], yend = 5*pcacontrib$Dim.2[5]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[5], y = 6*pcacontrib$Dim.2[5], label = pcacontrib$omy_miRNA[5]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[6], yend = 5*pcacontrib$Dim.2[6]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[6], y = 6*pcacontrib$Dim.2[6], label = pcacontrib$omy_miRNA[6]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[7], yend = 5*pcacontrib$Dim.2[7]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[7], y = 6*pcacontrib$Dim.2[7], label = pcacontrib$omy_miRNA[7]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[8], yend = 5*pcacontrib$Dim.2[8]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[8], y = 6*pcacontrib$Dim.2[8], label = pcacontrib$omy_miRNA[8]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[9], yend = 5*pcacontrib$Dim.2[9]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[9], y = 6*pcacontrib$Dim.2[9], label = pcacontrib$omy_miRNA[9]), color = "lightgrey", size = 3),
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, y = 0, xend = 5*pcacontrib$Dim.1[10], yend = 5*pcacontrib$Dim.2[10]), color = "lightgrey", arrow = ggplot2::arrow(length = ggplot2::unit(0.25, "cm"))),
    ggplot2::geom_text(mapping = ggplot2::aes(x = 6*pcacontrib$Dim.1[10], y = 6*pcacontrib$Dim.2[10], label = pcacontrib$omy_miRNA[10]), color = "lightgrey", size = 3))

  ellipse_points <- pcatab %>%
    dplyr::select(condition_list[1], "Dim.1", "Dim.2") %>%
    dplyr::mutate(dplyr::across(where(is.character), as.factor)) %>%
    as.data.frame() %>%
    FactoMineR::coord.ellipse(bary = T) %>%
    "$"("res")

  pca_theme <- list(
    ggplot2::xlab(paste0("Dimension 1 (",round(res.pca$eig[1,2],digits = 2),"%)")),
    ggplot2::ylab(paste0("Dimension 2 (",round(res.pca$eig[2,2],digits = 2),"%)")),
    ggplot2::theme_minimal(),
    ggplot2::geom_hline(yintercept = 0,linetype="dashed"),
    ggplot2::geom_vline(xintercept = 0,linetype="dashed"),
    ggplot2::theme(legend.position = "bottom"),
    ggplot2::scale_color_manual(values = gg_color(n_color)),
    ggplot2::ggtitle(paste("PCA of miRNA normalized counts", title, sep = " - "))
  )

  color_condition <- condition_list[1]
  shape_condition <- condition_list[2]

  if(is.na(shape_condition)) {
    base_plot <- ggplot2::ggplot(data = pcatab, ggplot2::aes_string(label="sample_ID", x="Dim.1", y="Dim.2",
                                                           colour=color_condition)) }
  else {
    base_plot <- ggplot2::ggplot(data = pcatab, ggplot2::aes_string(label="sample_ID", x="Dim.1", y="Dim.2",
                                                           colour=color_condition, shape = shape_condition))
  }
  if(add_theme) {
    base_plot <- base_plot + pca_theme
  }
  if(draw_arrows) {
    base_plot <- base_plot + pcacontrib_arrows[1:(2*n_arrows)]
  }
  if(draw_ellipses) {
    base_plot <- base_plot +
      ggplot2::geom_polygon(data=ellipse_points, ggplot2::aes_string(x="Dim.1",y="Dim.2",fill=condition_list[1]),
                            alpha=0.15,inherit.aes = F)
  }
  if(draw_sample_label) {
    base_plot <- base_plot + ggrepel::geom_text_repel(size=3)
  }

  base_plot +
    ggplot2::geom_point(size=2)
}

#' Boxplot of individual miRNA count compared by the first condition
#'
#' @param dds DESeq dataset
#' @param mir_of_interest List of mir to plot, ideally 25 or less
#' @param condition_list a string vector describing the different conditions
#' @param title title to add to the plot, default : "Project"
#' @param y_label string to add as a label on the y axis
#' @export
boxplot_mir <- function(dds, mir_of_interest,
                        condition_list = "conditions",
                        title = "Project",
                        y_label = "normalised counts - log scale"){
  dds %>%
    SummarizedExperiment::assay() %>%
    tidyr::as_tibble(rownames = "omy_miRNA") %>%
    dplyr::filter(.data$omy_miRNA %in% mir_of_interest) %>%
    tidyr::pivot_longer(-.data$omy_miRNA, names_to = "sample_ID", values_to = "count") %>%
    dplyr::left_join(
      dds %>%
        SummarizedExperiment::colData() %>%
        tidyr::as_tibble() %>%
        dplyr::select(dplyr::all_of(c("sample_ID", condition_list))),
      by = "sample_ID"
    ) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = condition_list[1], y = "count", color = condition_list[1])) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter() +
    ggplot2::facet_wrap("omy_miRNA", scales = "free") +
    ggplot2::labs(y = y_label) +
    ggplot2::ggtitle(paste("Individual miRNA count by condition - ", title, sep =""))
}

#' Boxplot of individual miRNA count compared by the first condition
#'
#' @param dds DESeq dataset
#' @param mir name of a single mir to plot (eg : "omy-miR-375-3p")
#' @param condition_list a string vector describing the different conditions
#' @param title title to add to the plot, default is the mir name
#' @param y_label string to add as a label on the y axis
#' @export
boxplot_one_mir <- function(dds, mir,
                            title = NA,
                            condition_list = "conditions",
                            y_label = "normalised counts - log scale"){
  if (is.na(title)) {
    title <- mir
  }

  dds %>%
    SummarizedExperiment::assay() %>%
    tidyr::as_tibble(rownames = "omy_miRNA") %>%
    dplyr::filter(.data$omy_miRNA == mir) %>%
    tidyr::pivot_longer(-.data$omy_miRNA, names_to = "sample_ID", values_to = "count") %>%
    dplyr::left_join(
      dds %>%
        SummarizedExperiment::colData() %>%
        tidyr::as_tibble() %>%
        dplyr::select(dplyr::all_of(c("sample_ID", condition_list))),
      by = "sample_ID"
    ) %>%
    ggplot2::ggplot(ggplot2::aes_string(x = condition_list[1], y = "count", color = condition_list[1])) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter() +
    ggplot2::labs(y = y_label) +
    ggplot2::ggtitle(title)
}

#' Barplot describing resampling analysis, for each miR found differential at least one time the number of time it was found differential
#'
#' @param reference_analysis tibble from compute_resampling function
#' @param resampling_results analysis whit every biosamples
#' @param n_resample the number of resampling to analyse
#' @param threshold_padj maximum adjusted pvalue to be considered differentially expressed (default 0.05)
#' @param threshold_lfc minimum log2 fold change to be considered differentially expressed (default 0.58, corresponds to 1.5 fold change)
#' @param title
#'
#' @importFrom rlang .data
#' @export
resampling_barplot <- function(reference_analysis, resampling_results, n_resample = 50, 
                               threshold_padj = 0.05, threshold_lfc = 0.58,
                               title = "") {
  resampling_results %>% 
    dplyr::mutate(condition = purrr::map(.data$analysis, 
                                         function(x) {x %>% 
                                           dplyr::filter(padj < threshold_padj) %>% 
                                           dplyr::filter(abs(log2FoldChange) > threshold_lfc) %>%
                                           dplyr::select(omy_miRNA, padj)})) %>%
    dplyr::select(.data$condition) %>% 
    dplyr::mutate(id = dplyr::row_number()) %>% 
    tidyr::unnest() %>% 
    tidyr::pivot_wider(names_from = .data$omy_miRNA, values_from = .data$padj) %>%
    dplyr::mutate(id = 0) %>%
    dplyr::group_by(.data$id) %>%
    dplyr::summarise_all(list(total = ~ sum(!is.na(.x))/n_resample*100, pval_mean = ~ mean(.x, na.rm = T))) %>%
    dplyr::select(-.data$id) %>%
    tidyr::pivot_longer(dplyr::everything()) %>%
    tidyr::extract(.data$name, into = c("omy_miRNA", "metric"), regex = "(.*p)_(.*)") %>%
    tidyr::pivot_wider(names_from = .data$metric, values_from = .data$value) %>%
    dplyr::left_join(reference_analysis %>% dplyr::filter(.data$padj < threshold_padj) %>% 
                       dplyr::select(.data$omy_miRNA) %>% dplyr::mutate(in_ref = T)) %>% 
    dplyr::mutate(in_ref = dplyr::if_else(is.na(.data$in_ref), "False", "True")) %>%
    dplyr::arrange(.data$total) %>% 
    dplyr::mutate(omy_miRNA = factor(.data$omy_miRNA, levels = .data$omy_miRNA)) %>%
    dplyr::filter(.data$total >= 5) %>%
    ggplot2::ggplot() + 
    ggplot2::geom_bar(ggplot2::aes_string(x = "omy_miRNA", y = "total", fill = "in_ref"), stat = "identity") + 
    ggplot2::labs(fill = "Differential when\nconsidering\nevery samples", y = "Differential in resampling (%)") + 
    ggplot2::coord_flip() + 
    ggplot2::scale_fill_manual(values = c("True" = "#F8766D", "False" = "grey")) + 
    ggplot2::ggtitle(title) + 
    ggplot2::scale_y_continuous(limits = c(0,100))
}

#' Heatmap describing resampling analysis
#'
#' @param sample_info sample_info tibble
#' @param resampling_results analysis whit every biosamples
#' @param subset_size number of sample for each condition in a resampling
#' @param threshold_padj maximum adjusted pvalue to be considered differentially expressed (default 0.05)
#' @param threshold_lfc minimum log2 fold change to be considered differentially expressed (default 0.58, corresponds to 1.5 fold change)
#' @param title
#' @param condition_name condition used to facet the heatmap
#'
#' @importFrom rlang .data
#' @export
resampling_heatmap <- function(sample_info, resampling_results, subset_size = 4, 
                               threshold_padj = 0.05, threshold_lfc = 0.58,
                               title = "", condition_name = "conditions") {
  resampling_results %>% 
    dplyr::mutate(condition = purrr::map(.data$analysis, 
                                         function(x) {x %>% 
                                             dplyr::filter(padj < threshold_padj) %>% 
                                             dplyr::filter(log2FoldChange > threshold_lfc) %>%
                                             dplyr::select(omy_miRNA, padj)})) %>%
    dplyr::select(.data$sample_list, .data$condition) %>% 
    tidyr::unnest(.data$sample_list) %>%
    dplyr::group_by(.data$sample_list) %>% dplyr::mutate(max_sample = dplyr::n()) %>% 
    tidyr::unnest(.data$condition) %>%
    dplyr::group_by(.data$sample_list, .data$omy_miRNA) %>%
    dplyr::summarise(tot = dplyr::n(), max_sample = .data$max_sample[1]) %>%
    dplyr::group_by(.data$omy_miRNA) %>% 
    dplyr::mutate(tot_miRNA = sum(.data$tot)/(2*subset_size)) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(.data$tot_miRNA) %>%
    dplyr::mutate(omy_miRNA = factor(.data$omy_miRNA, levels = unique(.data$omy_miRNA))) %>%
    dplyr::left_join(sample_info %>% dplyr::select("sample_ID", condition_name), by = c("sample_list" = "sample_ID")) %>%
    ggplot2::ggplot() + 
    ggplot2::geom_raster(ggplot2::aes(x = sample_list, y = omy_miRNA, fill = tot/max_sample*100)) + 
    ggplot2::labs(x = "Biosample", fill = "Portion of analysis\nincluding a biosample\nwhere a miRNA is found\ndifferentially expressed (%)") +
    ggplot2::scale_fill_gradient2(limits = c(0, 100)) + 
    ggplot2::theme(panel.grid = ggplot2::element_blank()) + 
    ggplot2::ggtitle(title) + 
    ggplot2::facet_wrap(condition_name, scales = "free_x")
}


#' Plot level of expression in every organ of omy miRNA 
#'
#' @param omy_miRNA list of miRNA patterns that will be interpreted as a regex
#' @param plot_title title of the plot
#' @param path_MiRNAOrigin data path of a MiRNAOrigin askomics table
#' @param path_MiRNA data path of a MiRNA askomics table
#' @param path_Organ data path of a Organ askomics table
#'
#' @importFrom rlang .data
#' @export
plot_mir_origin_from_juanchich_data <- function(omy_miRNA, plot_title, 
                                                path_MiRNAOrigin = "~/Documents/phenomir/askomics/datatables/MiRNAOrigin.csv",
                                                path_MiRNA = "~/Documents/phenomir/askomics/datatables/MiRNA.csv",
                                                path_Organ = "~/Documents/phenomir/askomics/datatables/Organ.csv") {
  data <- readr::read_csv(path_MiRNAOrigin, 
                          col_types = readr::cols(
                            MiRNAOrigin = readr::col_double(),
                            `referenceLevelOf@MiRNA` = readr::col_character(),
                            referenceLevel = readr::col_double(),
                            `referenceLevelFrom@Organ` = readr::col_character())) %>% 
    dplyr::left_join(readr::read_csv(path_MiRNA, 
                                     col_types = readr::cols(
                                       MiRNA = readr::col_character(),
                                       sequence = readr::col_character(),
                                       miRNAName = readr::col_character())), 
                     by = c("referenceLevelOf@MiRNA" = "MiRNA")) %>% 
    dplyr::left_join(readr::read_csv(path_Organ, 
                                     col_types = readr::cols(
                                       Organ = readr::col_character(),
                                       `rdfs:label` = readr::col_character())), 
                     by = c("referenceLevelFrom@Organ" = "Organ")) %>% 
    dplyr::mutate(organ = .data$`rdfs:label`) %>% 
    dplyr::filter(stringr::str_detect(.data$miRNAName, omy_miRNA))
  if(nrow(data) == 0) {
    ggplot2::ggplot() + 
      ggplot2::theme_void() + 
      ggplot2::geom_text(ggplot2::aes(x = 0, y = 0, label= paste0(plot_title,"\nNo data")))
  } else {
    ggplot2::ggplot(data) + 
      ggplot2::geom_bar(ggplot2::aes_string(x = "organ", y = "referenceLevel"), stat = "identity") +
      ggplot2::facet_wrap("miRNAName", scales = "free", ncol = 2) + ggplot2::coord_flip() + 
      ggplot2::ggtitle(plot_title)
  }
}