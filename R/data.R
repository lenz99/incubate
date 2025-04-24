#' Survival of mice with glioma under different treatments
#'
#' This data set stems from an animal experiment described in Stankovic (2018).
#' In particular, the data in question is shown in Figure 6J and 6K.
#'
#' @details
#' The data were read directly from the survival plots in the publication with
#' the help of Plot Digitizer, version 2.6.9.
#'
#' @format A data frame with 45 rows and 5 variables:
#' \describe{
#'   \item{Figure}{The figure in the publication where the data is shown}
#'   \item{Time}{Survival in days}
#'   \item{Status}{Right-censor status: 1 means observed event}
#'   \item{Group}{Experimental group identifier}
#'   \item{Colour}{Colour used in the Stankovic publication to mark this group}
#' }
#' @source Dudvarski Stankovic N, Bicker F, Keller S, et al. EGFL7 enhances surface expression of integrin a5b1 to promote angiogenesis in malignant brain tumors. EMBO Mol Med. 2018;10(9):e8420. doi:10.15252/emmm.201708420 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6127886/
"stankovic"


#' Relapse-free survival of melanoma patients under adjuvant treatment
#'
#' Data stem from a double-blind, placebo-controlled, phase 3 trial
#' where 870 patients with completely resected, stage III melanoma with BRAF V600E or V600K mutations
#' were randomly assigned to receive oral dabrafenib plus trametinib (combination therapy, 438 patients)
#' or two matched placebo tablets (432 patients) for 12 months.
#'
#' @details
#' The data were digitized by Sean Devlin from the survival plot Fig 1A in the
#' publication. Therefore, the data given here are **not** a 100% faithful
#' representation of the original published data and some deviations are to be
#' expected.
#'
#' @format A data frame with 870 rows and 4 variables:
#' \describe{
#'   \item{ID}{artificially generated patient ID}
#'   \item{time}{Time to relapse-free survival in months}
#'   \item{status}{Status of observation, encoded as 0 for right-censoring vs 1 for RFS-event}
#'   \item{trtmt}{Treatment group: Dabrafenib+Trametinib vs Placebo}
#' }
#' @source Long GV, Hauschild A, Santinami M, et al. Adjuvant Dabrafenib plus Trametinib in Stage III BRAF-Mutated Melanoma. N Engl J Med. 2017;377(19):1813-1823. doi:10.1056/NEJMoa1708539
#' @source Devlin SM and O'Quigley J, The nph2ph-transform: applications to the statistical analysis of completed clinical trials, arXiv:2407.18905, 2024. doi:10.48550/arXiv.2407.18905.
"long2017"


#' Small data sets from miscellaneous publications
#'
#' @description
#' Most data sets come from publications about parameter estimation in Weibull models.
#' See the references in the section "Source" below.
#'
#' @aliases rockette
#' @details
#' The following small data sets are provided as numeric vectors.
#' \describe{
#'   \item{`rockette`:}{Artificial sample of length 4 given by Rockette. The maximum likelihood function has two stationary points, none of them is the global maximum.}
#'   \item{`fatigue`:}{Fatigue times of ten bearings of a specific type in hours.}
#'   \item{`susquehanna`:}{Maximum flood levels (in millions of cubic feet per second) for the Susquehanna River of Harrisburg (Pennsylvania, USA) over 20 4-year periods.}
#'   \item{`pollution`:}{Beach pollution levels in South Wales (measured in number of coliform per 100 ml) on 20 days over a 5-week period.}
#'   \item{`graphite`:}{Breaking stress (in MPa x 10^6) of 41 beam specimens cut from a single graphite H590 block, from a reliability study reported by Margetson & Cooper (1984), cited by Cheng & Stephen (1989)}
#' }
#'
#' @references McCool, J.I., 1974. Inferential techniques for Weibull populations. Technical Report TR 74-0180, Wright Patterson Air Force Base, Ohio.
#' @references Rockette, H., 1974. Maximum Likelihood Estimation with the Weibull Model.
#' @references Dumonceaux, R. and Antle, C. E., 1973. Discrimination between the lognormal and the Weibull distributions. Technometrics, 15, 923-926.
#' @references Steen, P. J. and Stickler, D. J., 1976. A Sewage Pollution Study of Beaches from Cardiff to Ogmore. Report January 1976, Cardiff: Department of Applied Biology, UWIST.
#' @references Cheng, R.C.H. and Stephen, M.A., 1989. A Goodness of Fit Test Using Moran’s Statistic with Estimated Parameters. Biometrika, 76, 386-392.
"publication_examples"


#' @rdname publication_examples
"fatigue"

#' @rdname publication_examples
"susquehanna"

#' @rdname publication_examples
"pollution"

#' @rdname publication_examples
"graphite"
