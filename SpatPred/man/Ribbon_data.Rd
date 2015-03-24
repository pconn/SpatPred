\docType{data}
\name{Ribbon_data}
\alias{Ribbon_data}
\title{Ribbon seal count data and covariate data from the eastern Bering Sea}
\description{
  This data set includes objects summarizing effort, aerial survey counts, and covariates available for explaining ribbon seal apparent abundance in the eastern Bering Sea. Counts were obtained April 20-27, and covariates represent average values for this time period 2012 .  The Ribbon_data workspace includes the following objects:
  \itemize{
    \item{\code{Adj}} An adjacency matrix describing the spatial proximity of sampling units (sample units that are queen's move neighbors receive a 1, those that are not neighbors receive a 0)
    \item{\code{Adj2}} An adjacency matrix for a larger spatial extent (not used)
    \item{\code{Grid}} A list vector, each element of which provides a spatial polygons dataframe encoding the spatial position and covariates associated with sample units on each day.  For example, \code{Data$Grid[[1]]} holds covariate data for day 1 (4 April 2012).  Covariates provided include proportion land cover (\code{land_cover}), standardized distance from sample unit centroid to mainland (\code{dist_mainland}), standardized distance from land (including islands; \code{dist_land}), standardized distance from 1000m shelf break (\code{dist_shelf}), standardized depth (\code{depth}), proportion sea ice coverage (obtained from NSIDC; \code{ice_conc}), standardized distance from 90 percent sea ice contour (\code{dist_contour}), standardized distance from the southern sea ice edge (\code{dist_edge}), and Piatt and Springer ecoregion \code{Ecoregion}.
    \item{\code{Meta}} Provides some basic metadata
  }
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}
\keyword{data}

