% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_prost_data.R
\name{select_samples}
\alias{select_samples}
\title{Select samples to keep only samples with a minimum number of reads relative to
every samples.}
\usage{
select_samples(compressed_by_annotation, sample_info, number_of_sigma = 3)
}
\arguments{
\item{compressed_by_annotation}{compressed_by_annotation tibble}

\item{sample_info}{sample_info tibble}

\item{number_of_sigma}{number of standard error inferior to the mean a log total count can be before being filtered out}
}
\description{
The samples kept have at least a log total count superior to the mean log count
minus a user defined numerb of sigma (default 3)
}
