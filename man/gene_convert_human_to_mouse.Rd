% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_utils.R
\name{gene_convert_human_to_mouse}
\alias{gene_convert_human_to_mouse}
\title{Use biomaRT to convert human gene to mouse gene}
\usage{
gene_convert_human_to_mouse(genelist)
}
\arguments{
\item{genelist}{Gene list you want to convert}
}
\description{
For the human gene symbols that we can not find corresponding mouse symbols,
we will simply convert the first letter into Capital and return it
}
\examples{
gene_convert_human_to_mouse("TP53")
}
\keyword{standardize}
