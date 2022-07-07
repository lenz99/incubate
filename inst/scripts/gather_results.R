#!/usr/bin/env Rscript
# mkuhn, 2022-07-07
# gather individual simulation results

library('rlang')
suppressPackageStartupMessages(library('purrr'))
suppressPackageStartupMessages(library('dplyr'))
library('tidyr')
library('glue')
suppressPackageStartupMessages(library('R.utils'))



cmdArgs <- R.utils::commandArgs(trailingOnly=TRUE,
                                asValues = TRUE,
                                excludeReserved = FALSE, excludeEnvVars = TRUE,
                                defaults = list(resultsDir=file.path(getwd(), 'results'),
                                                resultsName='MS',
                                                type = 'test'))

if (any(c('help', 'h') %in% names(cmdArgs))){
  cat('Gather Monte-Carlo simulation results and save it as a common list.\n')
  cat('Parameter options are:\n')
  cat('  --help\t print this help\n')
  cat('  --resultsDir=\t specify the directory where to find and to put the result files. Defaults to directory "results".\n')
  cat('  --resultsName=\t specify a name suffix for results file. Default is "MS".\n')
  cat('  --type=\t specify if we gather results for "test" (default) or "confint".\n')
  cat('  --removeTemp\t flag to clean temporary results file after they have been saved.\n')
  quit(save = 'no')
}

myRemoveTemp <- isTRUE(any(c('r', 'removeTemp') %in% names(cmdArgs)))

myResultsDir <- cmdArgs[["resultsDir"]]
stopifnot( is.character(myResultsDir), dir.exists(myResultsDir) )

myResultsName <- cmdArgs[["resultsName"]]
stopifnot( is.character(myResultsName), length(myResultsName) == 1L, nzchar(myResultsName) )

myType <- cmdArgs[["type"]]
stopifnot( is.character(myType), length(myType) == 1L, nzchar(myType) )
myType <- match.arg(arg = tolower(myType), choices = c("test", "confint"))

simResFileNames <- list.files(myResultsDir,
                              pattern = paste0('simRes_',myType,'_[234]\\d+.+[.]rds$'),
                              full.names=TRUE)

if (! length(simResFileNames)) {
	cat('No temporary result files found!\n')
	q(save='no')
}


# run-namespace
RUN_NS <- tidyr::crossing(L1=LETTERS, L2=LETTERS) %>%
  dplyr::transmute(L=paste0(L1, L2)) %>%
  dplyr::pull(L)

# look at previous results (already saved)
resData <- NULL
indOffset <- 0L
RES_FILEN <- file.path(myResultsDir, paste0('simRes_', myType, '_', myResultsName,'.rds'))
if (file.exists(RES_FILEN)){
	resData <- readRDS(RES_FILEN)
	if ( ! is.list(resData) || is.null(names(resData))){
		cat('\nResults file is not a list! Start over from scratch!\n')
		resData <- NULL
	} else {
		cat(glue('There are already {length(resData)} entries saved!'), '\n')
		indOffset <- max(0L, which(RUN_NS %in% purrr::map_chr(resData, ~ { substr(.x[['run']][[1L]], start = 1L, stop = 2L) })), na.rm = TRUE)
		cat(glue('Offset index is {indOffset}.'), '\n')
	}
}

# read in temporary result files, process the results
# @return list name (from metadata) and data (unnested)
readResultFile <- function(rdsFN, ind) {
	rdsF <- readRDS(rdsFN)
	rdsFC <- comment(rdsF)
	mdList <- eval(parse(text = rdsFC))
	stopifnot( is.list(mdList), all( c('host', 'time') %in% names(mdList)) )

	resName <- paste(mdList[['host']], mdList[['time']], sep='||')
	stopifnot( indOffset + ind <= length(RUN_NS) )
	resData <- rdsF %>%
		tidyr::unnest(P) %>%
		dplyr::mutate(run = paste0(RUN_NS[[indOffset + ind]], run))

	# pass on comment to unnested dataframe
	comment(resData) <- rdsFC

	rlang::list2(!!resName := resData)
}

resCandidates <- purrr::imap(.x = simResFileNames, .f = readResultFile) %>%
  # drop list level from map
  purrr::flatten()

if (is.null(resData)){
	cat('\nStart with fresh results from scratch!\n')
	saveRDS(resCandidates, file = RES_FILEN)
} else {
	resDuplicates <- intersect(names(resData), names(resCandidates))
	if ( length(resDuplicates) ){
		cat(glue("Temporary results are already stored in results list:",
		         "{paste(resDuplicates, collapse = '\n  * ')}"), "\n\n")
		cat('Please clean up temporary results files that are already saved in result list, first!\n')
		q(save='no')
	}
	cat('\nAdd result candidates to existing result list!\n')
	saveRDS(c(resData, resCandidates), file = RES_FILEN)
}

if (myRemoveTemp){
  cat('\nAbout to remove temporary result files!\n')
  file.remove(simResFileNames)
}

cat('\n~fine~\n')
