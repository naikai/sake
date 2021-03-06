% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsne_utils.R
\name{plot_tsne}
\alias{plot_tsne}
\title{Plot t-SNE results}
\usage{
plot_tsne(tsne_out, color = NULL, alpha = 1, title = "tSNE",
  brewer = "naikai", add.centroids = F, add.line = F,
  add.errorbar = F, add.label = F, label.size = 3, conf = F,
  add.legend = F, save.plot = T, real.plot = T, point.size = 3,
  centroid.size = 6)
}
\arguments{
\item{tsne_out}{t-SNE results from run_tSNE}

\item{color}{color scheme for each dot}

\item{alpha}{alpha value for each dot}

\item{title}{figure title}

\item{brewer}{brewer color}
}
\description{
This function allows you to plot t-SNE results
}
\examples{
tsne_out <- run_tsne(mtcar, perplexity=3)
tsne_out <- parse_tsne_res(tsne_out)
plot_tsne(tsne_out)
}
\keyword{tsne}
