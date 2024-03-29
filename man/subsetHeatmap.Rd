% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_funclib_movAPA.r
\name{subsetHeatmap}
\alias{subsetHeatmap}
\title{Subset heatmapResults object}
\usage{
subsetHeatmap(
  heat,
  padjThd = 0.1,
  valueThd = NULL,
  upThd = NULL,
  dnThd = NULL,
  returnTF = FALSE
)
}
\arguments{
\item{heat}{a heatmapResults object}

\item{padjThd}{a cutoff for padj}

\item{valueThd}{a cutoff for value, to filter by abs(..)>=valueThd}

\item{upThd}{a cutoff to filter value by ...>=upThd (up-regulated)}

\item{dnThd}{a cutoff to filter by value by ...<=dnThd (down-regulated)}

\item{returnTF}{If TRUE, then return a matrix with TRUE or FALSE, otherwise return a heatmapResults object.}
}
\description{
heatmapResults is a generic class to store differential analyses results, which is part of a movRes object.
subsetHeatmap is to filter a heatmapResults by padj and/or value thresholds.
}
