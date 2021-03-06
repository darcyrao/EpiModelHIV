% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod.nwfeatures.WHAMP.R
\name{nwfeatures_msm_whamp}
\alias{nwfeatures_msm_whamp}
\title{Network diagnostics.}
\usage{
nwfeatures_msm_whamp(dat, at)
}
\arguments{
\item{dat}{Master data list object of class \code{dat} containing networks,
individual-level attributes, and summary statistics.}

\item{at}{Current time step.}
}
\value{
This function returns the \code{dat} object with the updated partnership list for diagnostics
on \code{dat$plist}.
}
\description{
Module function for tracking features of the sexual newtork over time.
}
\details{
At each time step, the cross-network degree distribution is calculated and stored, along with
statistics about race, region, and age mixing, and the duration of active ties.
}
\keyword{WHAMP}
\keyword{diagnostics}
\keyword{module}
\keyword{network}
