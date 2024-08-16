#' Survival of mice with glioma under different treatments
#'
#' This data set stems from an animal experiment described in Stankovic (2018).
#' In particular, the data in question is shown in Figure 6J and 6K.
#'
#' @details
#' The data were read directly from the survival plots in the publication with the help of Plot Digitizer, version 2.6.9.
#'
#' @format
#' \describe{
#'   \item{Figure}{The figure in the publication where the data is shown}
#'   \item{Time}{Survival in days}
#'   \item{Status}{Right-censor status: 1 means observed event}
#'   \item{Group}{Experimental group identifier}
#'   \item{Colour}{Colour used in the Stankovic publication to mark this group}
#' }
#' @source Dudvarski Stankovic N, Bicker F, Keller S, et al. EGFL7 enhances surface expression of integrin a5b1 to promote angiogenesis in malignant brain tumors. EMBO Mol Med. 2018;10(9):e8420. doi:10.15252/emmm.201708420 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6127886/
"stankovic"


#' Small data sets from different publications
#'
#' @description
#' Most data sets come from publications about parameter estimation in Weibull models.
#' See the references in the section "Source" below.
#'
#' @aliases rockette
#' @details
#' These small data sets are provided as numeric vectors.
#' \describe{
#'   \item{`rockette`:}{Artificial sample of length 4 given by Rockette. The maximum likelihood function has two stationary points, none of them is the global maximum.}
#'   \item{`fatigue`:}{Fatigue times of ten bearings of a specific type in hours.}
#'   \item{`susquehanna`:}{Maximum flood levels (in millions of cubic feet per second) for the Susquehanna River of Harrisburg (Pennsylvania, USA) over 20 4-year periods.}
#'   \item{`pollution`:}{Beach pollution levels in South Wales (measured in number of coliform per 100 ml) on 20 days over a 5-week period.}
#' }
#'
#' @source McCool, J.I., 1974. Inferential techniques for Weibull populations. Technical Report TR 74-0180, Wright Patterson Air Force Base, Ohio.
#' @source Rockette, H., 1974. Maximum Likelihood Estimation with the Weibull Model.
#' @source Dumonceaux, R. and Antle, C. E., 1973. Discrimination between the lognormal and the Weibull distributions. Technometrics, 15, 923-926.
#' @source Steen, P. J. and Stickler, D. J., 1976. A Sewage Pollution Study of Beaches from Cardiff to Ogmore. Report January 1976, Cardiff: Department of Applied Biology, UWIST.
"publication_examples"


#' @rdname publication_examples
"fatigue"

#' @rdname publication_examples
"susquehanna"

#' @rdname publication_examples
"pollution"
