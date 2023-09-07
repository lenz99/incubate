#!/usr/bin/env Rscript
# mkuhn, 2022-07-07
# gather individual simulation results from temporary files


# setup --------------------------------------------------------------------

library("incubate")
cat("incubate package version: ", toString(packageVersion("incubate")), "\n")

library('rlang')
library('dplyr', warn.conflicts = FALSE)
library('tidyr', warn.conflicts = FALSE)
library("purrr", warn.conflicts = FALSE)
library('glue')
suppressPackageStartupMessages(library('R.utils'))
stopifnot(packageVersion("purrr") > "1.0.0") #for list_flatten()


cmdArgs <- R.utils::commandArgs(trailingOnly=TRUE,
                                asValues = TRUE,
                                excludeReserved = FALSE, excludeEnvVars = TRUE,
                                defaults = list(resultsDir=file.path(getwd(), "results"),
                                                resultsTag="MS",
                                                type = "test"))

if (any(c('help', 'h') %in% names(cmdArgs))) {
  cat('Gather Monte-Carlo simulation results from temporary RDS result files and save it as a common list.\n')
  cat('Temporary data files have a date tag in their name.\n')
  cat('Parameter options are:\n')
  cat('  --help\t print this help\n')
  cat('  --resultsDir=\t specify the directory where to find and also where to put the result files. Defaults to sub-directory "results" of working directory.\n')
  cat('  --resultsTag=\t specify a name suffix for results file. Default is "MS".\n')
  cat('  --type=\t what type of results to gather? "test" (default) or "confint"\n')
  cat('  --removeTemp\t flag to clean temporary results file after they have been saved.\n')

  quit(save = 'no')
}

myRemoveTemp <- isTRUE(any(c('r', 'removeTemp') %in% names(cmdArgs)))

myResultsDir <- cmdArgs[["resultsDir"]]
stopifnot( is.character(myResultsDir), dir.exists(myResultsDir) )

myResultsTag <- cmdArgs[["resultsTag"]]
stopifnot( is.character(myResultsTag), length(myResultsTag) == 1L, nzchar(myResultsTag) )

myType <- cmdArgs[["type"]]
stopifnot( is.character(myType), length(myType) == 1L, nzchar(myType) )
myType <- match.arg(arg = tolower(myType), choices = c("test", "confint"))

# temporary results files have a date tag in their file name
simResFileNames <- list.files(myResultsDir,
                              pattern = paste0('simRes_',myType,'_[234]\\d+.+[.]rds$'),
                              full.names=TRUE)

if (! length(simResFileNames)) {
	cat("No matching temporary result files for ", myType,
	    " were found in sub-directory ", myResultsDir, "!\n")
	q(save="no")
}


# create namespace for the different runs as vector
RUN_NS <- tidyr::crossing(L1=LETTERS, L2=LETTERS) %>%
  dplyr::mutate(L=paste0(L1, L2), .keep = "none") %>%
  dplyr::pull(L)


# previous results --------------------------------------------------------

# look at previous results (already saved)
resData <- NULL
indOffset <- 0L
RES_FILEN <- file.path(myResultsDir, paste0('simRes_', myType, '_', myResultsTag,'.rds'))

# prepare index offset (if there are previous results)
if (file.exists(RES_FILEN)) {
	resData <- readRDS(RES_FILEN)

	if ( ! is.list(resData) || is.null(names(resData))) {
	  # no results data!
		cat('\nResults file is not a list! We start over from scratch!\n')
		resData <- NULL
	} else {
	  # results data found!
	  cat(glue("There have been {length(resData)} entries saved already!"), "\n")
	  # run-numbers only for tests (currently not for confint)
	  if (myType == "test") {
	    indOffset <- max(0L,
	                     which(RUN_NS %in% purrr::map_chr(resData,
	                                                      .f = ~ { substr(.x[['run']][[1L]], start = 1L, stop = nchar(RUN_NS[[1L]])) })),
	                     na.rm = TRUE)
	    cat(glue("Offset index from saved results is {indOffset}."), "\n")
	  }#fi
	}#fi resData
}


# read in temporary results data --------------------------------------------------

# read in temporary result file, process the results
# @param rdsFN file name of RDS-file containing results
# @param ind index number of given rds filename
# @return named list containing the unnested data
readResultFile <- function(rdsFN, ind) {
	rdsF <- readRDS(rdsFN)
	rdsFC <- comment(rdsF)
	mdList <- eval(parse(text = rdsFC))
	stopifnot(is.list(mdList), all(c('host', 'time') %in% names(mdList)))

	resName <- paste(mdList[['host']], mdList[['time']], sep='||')

	unnestVar <- if (myType == "confint") "ci_res" else  "results"
	stopifnot(all(unnestVar %in% names(rdsF)))
	resData <- rdsF %>%
		tidyr::unnest(cols = all_of(unnestVar))

	if (myType == "test") {
	  stopifnot(indOffset + ind <= length(RUN_NS))
	  resData <- resData %>%
	    # prepend run namespace, like 'AS'
	    dplyr::mutate(run = paste0(RUN_NS[[indOffset + ind]], run))
	}

	# pass on comment to unnested dataframe
	comment(resData) <- rdsFC

	rlang::list2(!!resName := resData)
}


resCandidates <- purrr::imap(.x = simResFileNames, .f = readResultFile) %>%
  # drop imap's additional list level
  purrr::list_flatten()



# save all results data ---------------------------------------------------

if (is.null(resData)) {
	cat('\nStart with fresh results from scratch!\n')
	saveRDS(resCandidates, file = RES_FILEN)
} else {
  # check for duplicates
	resDuplicates <- intersect(names(resData), names(resCandidates))
	if (length(resDuplicates)) {
		cat(glue("These temporary result dataframes are already stored in the {myResultsTag}-results file:\n  * ",
		         "{paste(resDuplicates, collapse = '\n  * ')}", .trim = FALSE), "\n\n")
		cat('Please clean up temporary results files that are already saved in result list, first!\n')
		q(save='no')
	}
	cat("\nAdd result candidates to existing result list!\n")
	saveRDS(c(resData, resCandidates), file = RES_FILEN)
}

if (myRemoveTemp) {
  cat("\n")
  cat("About to remove temporary result files!\n")
  fnRmv <- file.remove(simResFileNames)
  cat("Successfully removed ", sum(fnRmv), " temporary RDS-files.\n")
  if (sum(fnRmv)) {
    cat("\n  * ")
    cat(paste(simResFileNames[fnRmv], collapse = "\n  * "))
  }
  cat("\n")
}

cat("\n~fine~\n")
