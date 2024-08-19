# S3-integration functions

#' @export
print.incubate_fit <- function(x, ...) {
  coe <- coef(x)
  rangeTime <- if (x[["twoGroup"]]) {
    ns <- x[["nobs"]]
    paste(
      round(sort(c(x[["data"]]$x[[1L]], x[["data"]]$y[[1L]]))[[1L]], 4),
      round(sort(c(x[["data"]]$x[[ns[["x"]]]], x[["data"]]$y[[ns[["y"]]]]))[[2L]], 4),
      sep = " to ")
  } else {
    paste(round(x$data[[1L]],4), round(x[["data"]][[length(x$data)]],4), sep = " to ")
  }
  cat(glue::glue_data(x, .sep = "\n",
                      "Fit a {distO$dist_name}{c('', ' with two delay phases')[[1L+twoPhase]]} through{c('', ' profiled')[[1L+optimizer$profiled]]} {switch(method,
                      MPSE = 'Maximum Product of Spacings Estimation (MPSE)', MLEn = 'naive Maximum Likelihood Estimation (MLEn)',
                      MLEw = 'weighted Maximum Likelihood Estimation (MLEw)',
                      MLEc = 'corrected Maximum Likelihood Estimation (MLEc)', '???')} for {c('a single group', 'two independent groups')[[1L+twoGroup]]}.",
                      "Data: {if (twoGroup) paste(nobs, collapse = ' and ') else nobs[[1L]]} observations, ranging from {rangeTime}",
                      "Criterion: {signif(criterion,3)}",
                      "Fitted coefficients: {if (is.null(coe)) '-' else paste(paste('\n  ', names(coe)), signif(coe,5L), sep = ': ', collapse = ' ')}"),
      "\n")
}

#' Coefficients of a delay-model fit.
#' @param object object that is a `incubate_fit`
#' @param transformed flag. Do we request the transformed parameters as used within the optimization?
#' @param group character string to request the canonical parameter for one group
#' @param ... further arguments, currently not used.
#' @return named coefficient vector
#' @export
coef.incubate_fit <- function(object, transformed = FALSE, group = NULL, ...) {
  stopifnot( inherits(object, "incubate_fit") )
  transformed <- isTRUE(transformed)

  rlang::env_get(rlang::fn_env(object$objFun), nm = "extractPars")(purrr::chuck(object, !!! if (transformed) list("optimizer", "parOpt") else "par"),
                                                                   group = group, isOpt = transformed, transform = FALSE, named = TRUE)
}

#' @export
summary.incubate_fit <- function(object, ...) {
  print(object)
}



#' Plot a fitted delay-model object of class `incubate_fit`
#'
#' The fitted delay-model is plotted: a Kaplan-Meier survival curve is shown together with the parametric model fit.
#'
#' @details
#' This function requires the `ggplot2`-package to be installed.
#'
#' @param x a fitted delay-model
#' @param y not used
#' @param title character. Optionally, provide a title to the plot.
#' @param subtitle character. Optionally, provide a subtitle to the plot. By default the coefficients are shown.
#' @param xlim numeric. Optionally, limits for the x-axis (time). If unspecified starts from 0 to last observation.
#' @export
plot.incubate_fit <- function(x, y, title, subtitle, xlim, ...) {
  # parameter y comes from the plot-generic but y is not used here.
  stopifnot(inherits(x, "incubate_fit"))

  rlang::check_installed(pkg = 'ggplot2', reason = 'to draw plots', version = '3.3')

  distO <- x$distO
  cumFun <- distO$cdf

  # add time = 0 per group
  kmFit0 <- survival::survfit0(x[["kmFit"]], start.time = 0)
  kmFit0 <- tibble(group = if (is.null(kmFit0$strata)) "x" else rep.int(c("x", "y"), times = kmFit0$strata),
                   time = kmFit0$time,
                   n.risk = kmFit0$n.risk,
                   n.event = kmFit0$n.event,
                   n.censor = kmFit0$n.censor,
                   surv = kmFit0$surv,
                   evrate = 1 - surv)


  # add estimated delay model
  p <- if (x[["twoGroup"]]) {
    ggplot2::ggplot(data = kmFit0,
                    mapping = ggplot2::aes(x = .data$time, y = .data$evrate, col = .data$group)) +
      ggplot2::geom_function(mapping = ggplot2::aes(col = rep.int("x", NROW(kmFit0))),
                             fun = cumFun, args = coef(x, group = "x"), linetype = "dashed") +
      ggplot2::geom_function(mapping = ggplot2::aes(col = rep.int("y", NROW(kmFit0))),
                             fun = cumFun, args = coef(x, group = "y"), linetype = "dashed")
  } else {
    ggplot2::ggplot(data = kmFit0,
                    mapping = ggplot2::aes(x = .data$time, y = .data$evrate)) +
      ggplot2::geom_function(inherit.aes = FALSE, fun = cumFun, args = coef(x, group = "x"), linetype = "dashed")
  }

  p <- p +
    # kaplan meier step function
    ggplot2::geom_step() +
    # mark (right-)censored observations
    ggplot2::geom_point(data = function(.x) .x[.x$n.censor > 0,], shape = 3L)


  if (missing(title)) title <- glue::glue_data(x,
                                               "Fitted {distO$dist_name} {c('model ', 'models ')[[1L+twoGroup]]}",
                                               "{c('', 'with two delay phases')[[1L+twoPhase]]}")
  coefPrint <- function(gr) {
    co <- coef.incubate_fit(x, group = gr)
    paste(names(co), signif(co, 4), sep = ": ", collapse = ", ")
  }
  if (missing(subtitle)) subtitle <- if (x[["twoGroup"]]) paste(coefPrint("x"), coefPrint("y"), sep = " - ") else coefPrint("x")

  if (missing(xlim) || is.null(xlim)) xlim <- c(0L, NA)

  p +
    #ggplot2::xlim(0L, NA) +
    # transforms "after_stat" which matters for stat_ecdf
    ggplot2::coord_trans(y = "reverse", xlim = xlim) +
    ggplot2::labs(x = 'Time', y = 'Cumulative prop. of events',
                  col = if (x[["twoGroup"]]) 'Group' else NULL,
                  title = title, subtitle = subtitle)

}




#' Transform observed data to unit interval
#'
#' The transformation used is the probability integral transform:
#' the cumulative distribution function with the estimated parameters of the model fit takes the data into the 0-1 interval.
#' All available data in the model fit is transformed. Censored observations lead to censored back-transformed observations as well.
#'
#' @note
#' This S3-method implementation is quite different from its default method that allows for non-standard evaluation on data frames, primarily intended for interactive use.
#' But the name `transform` fits so nicely to the intended purpose that it is re-used for the probability integral transform, here.
#'
#' @param _data a fitted model object of class `incubate_fit`
#' @param ... currently ignored
#' @return The transformed data, either a vector (for single group) or a list with entries x and y (in two group scenario)
#' @export
transform.incubate_fit <- function(`_data`, ...) {
  stopifnot(inherits(`_data`, "incubate_fit"))

  cdfFun <- `_data`$distO$cdf

  twoGroup <- isTRUE(`_data`$twoGroup)
  isSurv <- isTRUE(`_data`$cens$isSurv)

  x <- if (twoGroup) `_data`$data$x else `_data`$data

  tr <- NULL

  if (isSurv) {
    # currently, handle right-censored case only
    stopifnot( attr(x, which = "type", exact = TRUE) == "right")
    tr <- Surv(time = rlang::exec(cdfFun, !!! c(list(q=x[,1L]), coef(`_data`, group = "x"))),
               event = x[, "status"], type = "right")
    if (twoGroup) tr <- list(x = tr,
                             y = Surv(time = rlang::exec(cdfFun, !!! c(list(q=`_data`$data$y[,1L]), coef(`_data`, group = "y"))),
                                      event = `_data`$data$y[, "status"], type = "right"))
  } else {
    tr <- rlang::exec(cdfFun, !!! c(list(q=x), coef(`_data`, group = "x")))
    if (twoGroup) tr <- list(x = tr,
                             y = rlang::exec(cdfFun, !!! c(list(q=`_data`$data$y), coef(`_data`, group = "y"))))
  }

  tr
}
