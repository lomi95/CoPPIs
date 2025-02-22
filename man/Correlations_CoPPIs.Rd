% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Correlations_CoPPIs.R
\name{Correlations_CoPPIs}
\alias{Correlations_CoPPIs}
\title{Title}
\usage{
Correlations_CoPPIs(
  list.groups,
  names_of_groups,
  categories.interactome,
  gCOR.groups.l,
  min_n.corr,
  signifCorr,
  correctionCorr,
  corr.test,
  gInter,
  compute_weights
)
}
\arguments{
\item{list.groups}{list of dataset}

\item{names_of_groups}{names of the groups}

\item{categories.interactome}{interactome for each cateogry}

\item{gCOR.groups.l}{if a precedent analysis was run, an object gCOR.groups.l may be used to speed up the computation}

\item{min_n.corr}{minimum number of observation pairs for computing the correlation}

\item{signifCorr}{significance threshold for correlation pvalue}

\item{correctionCorr}{pvalue adjusting method}

\item{corr.test}{type of correlation test (Spearman or Pearson)ù}

\item{gInter}{graphs of interactions}

\item{compute_weights}{boolean, the weight transformation to the non significant
correlation need to be computed?}
}
\value{
a list containing information about correlation of conditions
}
\description{
Title
}
