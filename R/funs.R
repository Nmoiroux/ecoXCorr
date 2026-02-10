#' Plot a cross-correlation map (CCM) from lagged regression results
#'
#' This function visualises the strength and direction of associations between
#' a response variable and a lagged predictor across multiple lag windows \insertCite{currieroCrossCorrelationMaps2005}{ecoXCorr}, using
#' the output of \code{\link{fit_models_by_lag}}. The resulting plot is a
#' two-dimensional "cross-correlation map", where each tile represents a
#' lag window defined by \code{lag_start} and \code{lag_end}.
#'
#' The colour of each tile corresponds to \code{model_outcome}. Positive associations are shown in red, negative
#' associations in blue, and non-significant or filtered values (using \code{threshold_p}) are shown in grey.
#'
#' The lag window yielding the maximum absolute value of \code{model_outcome} is highlighted
#' with a coloured border.
#'
#' @param data A data.frame produced by \code{\link{fit_models_by_lag}}, containing
#'   at least the columns \code{lag_start}, \code{lag_end}, \code{r2},
#'   \code{p_value}, and \code{sign}.
#' @param model_outcome Character string specifying the model's outcomes to plot. Either:
#' \itemize{
#'   \item \code{"d_aic"} for the delta AIC (compared to the null model) (default),
#'   \item \code{"R2sign"} for the signed coefficient of determination,
#'   computed as the marginal or classical \eqn{R^2} multiplied by the sign of the estimated effect,
#'   \item \code{"R2"} for the coefficient of determination (\eqn{R^2}),
#'   \item \code{"betas"} for the estimated beta parameter (slope) of the linear predictor or,
#'   \item \code{"weight"} for the Akaike weight.
#'   }
#' @param threshold_p Numeric value giving the p-value threshold above which
#'   associations are masked (set to \code{NA}) in the plot. Filtering is performed on
#'   the adjusted (for multiple testing) p-values \code{p_adj}. (Default is \code{1},
#'   meaning that no filtering is applied.
#'
#' @return A \code{ggplot2} object representing the cross-correlation map.
#'
#' @details
#' The x-axis corresponds to \code{lag_start} (displayed in reverse order), and
#' the y-axis corresponds to \code{lag_end}. Tiles are coloured using a diverging
#' colour scale centred on zero. Lag windows with \code{p_value >= threshold_p}
#' are not displayed and appear in grey.
#'
#' This function does not perform any modelling itself; it is intended solely
#' for visualising results obtained from \code{\link{fit_models_by_lag}}.
#'
#' @references
#'  \insertAllCited{}
#'
#' @seealso \code{\link{fit_models_by_lag}}, \code{\link[ggplot2]{ggplot}}
#'
#' @examples
#' res_glm <- ecoXCorr(
#' meteo_data    = meteoMPL2023,
#' response_data = albopictusMPL2023,
#' date_col_meteo   = "date",
#' date_col_resp    = "date",
#' value_cols    = "rain_sum",
#' agg_fun       = "sum",
#' response      = "individualCount",
#' lag_unit      = 7,
#' max_lag       = 8,
#' family        = "poisson"
#' )
#'
#' plotCCM(res_glm, model_outcome = "R2sign", threshold_p = 0.05)
#' plotCCM(res_glm, model_outcome = "R2")
#' plotCCM(res_glm, model_outcome = "d_aic")
#' plotCCM(res_glm, model_outcome = "betas", threshold_p = 0.05)
#'
#' @export
plotCCM <- function(data,
                    model_outcome = c("d_aic","R2sign","R2","betas","weight"),
										threshold_p = 1){

  model_outcome <- match.arg(model_outcome)
  indicator <- model_outcome

  # get R2 name from data
  R2_type <- names(data)[which(grepl('^R2(?!s)',names(data), perl=TRUE) ==TRUE)]
  if (indicator == "R2"){
    indicator <- R2_type
  }


  # filter data according to p-value
	data[data$p_adj>threshold_p, indicator] <- NA

	# find max absolute indicator
	abs_ind <- abs(data[indicator]) # find max absolute of value(s) of indicator
	index_max <- which(abs_ind == max(abs_ind, na.rm = T)) # index of max


	# legend name
	if (indicator == "R2sign"){
	  name_legend <- paste("signed", R2_type, sep="\n")
	} else if (indicator == "d_aic"){
	  name_legend <- "delta AIC"
	} else {
	  name_legend <- indicator
	}

	# legend limits
	limits <- c(min(data[indicator],na.rm=T),max(data[indicator],na.rm=T))

	# value to plot
	data$value <- data[,indicator]
	maxdata <- data[index_max,]

  # plot
	plot <- ggplot(data = data, aes(lag_start, lag_end, fill = value)) +
	  geom_tile() +
	  geom_tile(data = maxdata , color = "deeppink3", linewidth = 1, show.legend = FALSE)+
	  scale_x_reverse() +
	  scale_fill_gradient2(
			low = "blue",
			high = "red",
			mid = "white",
			midpoint = 0,
			limits = limits,
			name = name_legend,
			na.value = "grey"
		) +
	  theme_bw()

	# message
	message(paste0("Cross correlation maps showing ",model_outcome," of models fitted on ", unique(data$response)," with lagged effects of ", unique(data$predictor) ," as predictor. Pink-bordered square highlight the time lag with the highest absolute value of ",model_outcome,". Grey squares show intervals whith adjusted (for multiple testing) p-values > ",threshold_p,"."))

	# return plot
	return(plot)
}


#' Fit regression models by lag window on aggregated meteorological predictors
#'
#' This function fits a regression model separately for each lag window
#' defined by the \code{lag_start} and \code{lag_end} columns of the input
#' data frame. For each lag window, the model is fitted using
#' observations corresponding to different reference dates (\code{date}),
#' and summary statistics (betas, sign of effect, \eqn{R^2}, delta AIC, Akaike weight, sample size, p-value, p-value adjusted for multiple testing) are returned
#' for the specified predictor.
#'
#' Both fixed-effect and mixed-effect models are supported.
#'
#' The modelling function used depends on the \code{random} and \code{family} arguments:
#' \itemize{
#'   \item \code{random} is not specified (default):  \code{\link[stats]{glm}}
#'   \item \code{random} is note an empty string OR \code{family} is a valid glmmTMB family: \code{\link[glmmTMB]{glmmTMB}}
#' }
#'
#' For mixed-effects models, marginal \eqn{R^2} (Nakagawa) is returned. For fixed-effects
#' models, appropriate \eqn{R^2} is used (see \code{?performance::r2()}).
#' Depending on model specification (depending on \code{random} and/or \code{family}), \eqn{R^2} for \code{\link[glmmTMB]{glmmTMB}} models may not be computed: the returned error or warning is printed in the console.
#'
#' @param data A data frame containing, at minimum, the columns
#'   \code{lag_start}, \code{lag_end}, \code{date}, the response variable,
#'   the predictor variable(s) and optional random-effect variables.
#'
#' @param response Character string giving the name of the response variable.
#'
#' @param predictors Character vector of predictor names. Currently, only
#'   a single predictor is supported; providing more than one predictor
#'   will result in an error.
#'
#' @param random Optional character string specifying random-effects terms or
#'   covariance structure to be added to the model formula (without a leading \code{+}), e.g.
#'   \code{"(1 | site/year)"}, \code{"(1 | site) + (1 | year)"} or \code{"ar1(times + 0 | group)"} (\code{?glmmTMB::glmmTMB}).
#'   If empty (default), a fixed-effect model is fitted (using \code{glm()}.
#'
#' @param family Character string. The name of a family function
#'   to be used in GLM or GLMM models. Default to "gaussian" (Linear model).
#'   see \code{?stats::family} and \code{?glmmTMB::family_glmmTMB}
#'
#' @param min_n Minimum number of observations required to fit a model.
#'   (Currently not enforced; retained for future extensions.)
#'
#' @param track If TRUE, lag window is printed in the console before model fitting.
#'
#' @param ... Additional arguments passed to the underlying modelling
#'   function (\code{glm}, or \code{glmmTMB::glmmTMB}).
#'
#' @details
#' For each unique combination of \code{lag_start} and \code{lag_end}, the
#' function:
#' \enumerate{
#'   \item Subsets the data to the corresponding lag window,
#'   \item Removes rows with missing values in the response or predictor,
#'   \item Fits the specified model,
#'   \item Extracts beta parameter of the linear predictor,
#'   \item Extracts the p-value of the predictor effect,
#'   \item Computes AIC for the specified and null models,
#'   \item Computes the appropriate model\eqn{R^2} (marginal Nakagawa \eqn{R^2} for mixed models),
#'   \item Records the sign of the estimated effect and the sample size,
#'   \item Computes delta-AIC and Akaike weight of each model \insertCite{vandepolIdentifyingBestClimatic2016}{ecoXCorr},
#'   \item Computes adjusted p-value using the False Discovery Rate method \insertCite{benjaminiControlFalseDiscovery2001}{ecoXCorr}
#' }
#'
#' The returned table is suitable for lag-window screening, heatmap
#' visualisation, or sensitivity analyses in epidemiological or ecological
#' studies.
#'
#' @return A data frame with one row per lag window, containing:
#'   \describe{
#'     \item{lag_start}{Start lag index of the aggregation window.}
#'     \item{lag_end}{End lag index of the aggregation window.}
#'     \item{response}{Name of the response variable.}
#'     \item{predictor}{Name of the predictor variable.}
#'     \item{R2...}{Coefficient of determination (marginal \eqn{R^2} for mixed models).}
#'     \item{betas}{Estimated beta parameter of the linear predictor.}
#'     \item{sign}{Sign of the estimated predictor effect (-1 or +1).}
#'     \item{d_aic}{AIC reduction compared to the null model.}
#'     \item{n}{Number of observations used to fit the model.}
#'     \item{p_value}{P-value associated with the predictor effect.}
#'     \item{weight}{Akaike weight.}
#'     \item{p_adj}{P-value adjusted for multiple testing.}
#'   }
#'
#' @references
#'  \insertAllCited{}
#'
#' @seealso
#' \code{\link[glmmTMB]{glmmTMB}},
#' \code{\link[performance]{r2}},
#' \code{\link[performance]{r2_nakagawa}}
#'
#' @examples
#' sampling_dates <- unique(albopictusMPL2023$date)
#'
#' met_agg <- aggregate_lagged_intervals(
#' data       = meteoMPL2023,
#' date_col   = "date",
#' value_cols = c("rain_sum", "temp_mean"),
#' d          = sampling_dates,
#' i          = 7,
#' m          = 8
#' )
#'
#' albo_lag <- merge(met_agg, albopictusMPL2023, by = "date", all = TRUE)
#'
#' res_glm <- fit_models_by_lag(
#' data       = albo_lag,
#' response   = "individualCount",
#' predictors = "rain_sum_sum",
#' random     = "",
#' family     = "poisson"
#' )
#'
#' head(res_glm)
#'
#' @export
fit_models_by_lag <- function(data,
															response,
															predictors,
															random = "", # Random-effects terms to be added to the formulae, wihtout initial "+", e.g. "(a|b/c)+(a|d)"
															family = "gaussian",
															min_n = 10,
															track = F,
															...) {

	out <- data


	# Ensure response exists
	stopifnot(response %in% names(out))
	stopifnot(all(predictors %in% names(out)))
	stopifnot(is.character(family))
	stopifnot(family %in% c(.glmmtmb_family,.glm_family))


	# Bi- or multi-variate model ?
	if (length(predictors) >1 ){
		multiv <- TRUE
		stop("The function does not support multiple predictors")
	}

	# Unique lag windows
	lag_windows <- unique(out[, c("lag_start", "lag_end")])
	n_mod_to_fit <- nrow(lag_windows)

	# mixed-effect model ?
	if (nchar(random) > 0 ){
		mixed <- TRUE
		if (n_mod_to_fit > 30){
			message(
				"There are ",	paste0(n_mod_to_fit)," models to be fitted, which may take some time..."
			)
		}
	} else {
		mixed <- FALSE
		if (n_mod_to_fit > 100){
		  message(
		    "There are ",	paste0(n_mod_to_fit)," models to be fitted, which may take some time..."
		  )
		}
	}

	# Build formula
	if (mixed == TRUE) {
	  fml <- as.formula(
	    paste(response, "~", paste(c(predictors,random), collapse = " + ")))
	  fml_null <- as.formula(
	    paste(response, "~", paste(c("1",random), collapse = " + ")))
	} else {
	  fml <- as.formula(
	    paste(response, "~", paste(predictors, collapse = " + ")))
	  fml_null <- as.formula(
	    paste(response, "~ 1"))
	}

	mod_fun <- ifelse(mixed == T | (family %in% .glmmtmb_family), "glmmTMB", "glm")

	message(paste0("The fitted models will be of the form: ",mod_fun,"(",deparse1(fml),", family=",family,")"))


	results <- list()

	warns <- c()
	k <- 1

	for (i in seq_len(nrow(lag_windows))) {

		ls <- lag_windows$lag_start[i]
		le <- lag_windows$lag_end[i]

		# Subset data for this lag window
		dat <- out[out$lag_start == ls &
							 	out$lag_end   == le, ]


		n <- nrow(dat)
		#if (n < min_n) next
		if (var(dat[,predictors[1]])==0) next # will lead to singularities so next


		if (track == T){
			print(lag_windows[i,])
		}

		# Fit null model (only required one time), allways the same
		if (i==1 & mod_fun == "glmmTMB"){
		  fit_null <- glmmTMB(fml_null, data = dat, family = family, ...)
		} else if(i==1 & mod_fun == "glm"){
		  fit_null <- glm(fml_null, data = dat, family = family, ...)
		}

		# Fit model and extract betas and p-value
		if (mod_fun == "glmmTMB"){
			fit <- glmmTMB(fml, data = dat, family = family, ...)
			sm <- summary(fit)
			pval <- sm$coefficients$cond[predictors[1],4]
			betas <- sm$coefficients$cond[predictors[1],1]
			sign <- sign(betas)
		} else if(mod_fun == "glm"){
			fit <- glm(fml, data = dat, family = family, ...)
			sm <- summary(fit)
			pval <- sm$coefficients[predictors[1],4]
			betas <- sm$coefficients[predictors[1],1]
			sign <- sign(betas)
		}

    # Compute R2
		r2_res <- tryCatch({
		  r2(fit, null_model = fit_null)
		},error = function(e) {
		  warns <<- c(warns, paste0("Error in R2 computing (performance::r2()):",
		                            conditionMessage(e)))
		  return(NA)
		},warning = function(w){
		  warns <<- c(warns, paste0("Warning in R2 computing (performance::r2()):",
		                            conditionMessage(w)))
		  return(NA)
		} )

		if (mixed == TRUE & !all(is.na(r2_res))){
		  r2 <- r2_res[2]
		}else if(mixed == FALSE & !all(is.na(r2_res))){
		  r2 <- r2_res[1]
		}else if(all(is.na(r2_res))){
		  r2 <- NA
		}

		# compute R2sign
		R2sign <- unlist(r2) * sign


		# compute delta AIC
		aic <- performance_aic(fit)
		aic_null <- performance_aic(fit_null)
		d_aic <- aic - aic_null


		# fill results dataframe
		results[[k]] <- data.frame(
			lag_start = ls,
			lag_end   = le,
			response  = response,
			predictor = predictors[1],
			R2        = r2,
			R2sign    = R2sign,
			betas     = betas,
			sign      = sign,
			d_aic     = d_aic,
			n         = n,
			p_value   = pval
		)

		k <- k + 1
	}


	results <- do.call(rbind, results)

	# Akaike weight
	results$weight <- exp(-0.5*results$d_aic) / (sum(exp(-0.5*results$d_aic)))

	# Adjust p-value for multiple testing
	results$p_adj <- p.adjust(results$p_value, method = "BY") # adjust p_value for multiple testing


	rownames(results) <- 1:nrow(results)

	message("Done !")

	message(unique(warns))

	return(results)

}




#' Aggregate meteorological time series over multiple lagged time intervals
#'
#' This function computes aggregated values of one or several meteorological
#' time series over all possible lagged time intervals defined relative to one
#' or more reference dates. For each reference date \code{d}, all intervals
#' \eqn{[d - k \times i,\; d - (l-1) \times i)} are generated,
#' where \code{i} is the interval length (in days) and \code{k, l} range from 1 to \code{m} with
#' \code{k >= l}. Each interval is then used to aggregate the specified
#' meteorological variables using one or more summary functions.
#'
#' The function supports multiple reference dates, multiple variables, and
#' multiple aggregation functions, and returns all combinations as additional
#' columns in the output data frame.
#'
#' Reference dates for which at least one interval contains no observations are
#' reported in the console as having missing data. Reference dates for which at
#' least one interval partially lies outside the temporal bounds of the input
#' time series are reported as having truncated intervals.
#'
#' @param data A data.frame containing the meteorological time series.
#' @param date_col Character string giving the name of the date column in
#'   \code{data}. The column must be of class \code{Date} or \code{POSIXct}.
#' @param value_cols Character vector giving the names of numeric variables
#'   to be aggregated (e.g. rainfall, temperature).
#' @param d Vector of reference *dates*. Can be of class \code{Date} or coercible
#'   to \code{Date}. Aggregations are computed independently for each date.
#' @param i Integer giving the length of the base time *interval*, expressed in days
#'   (e.g. \code{1} for daily data, \code{7} for weekly intervals,
#'   \code{14} for fortnightly intervals).
#' @param m Integer giving the *maximum* lag (number of intervals) to consider.
#'   All combinations of lag windows with \code{1 <= lag_end <= lag_start <= m}
#'   are evaluated.
#' @param funs Named list of aggregation functions to apply to each variable.
#'   Each function must accept a numeric vector as first argument. The names
#'   of the list are used to construct output column names
#'   (e.g. \code{rain_mean}, \code{temp_max}). Defaults to mean, min, max and sum.
#' @param na.rm Logical indicating whether missing values should be removed
#'   before aggregation. Passed to the aggregation functions (default: \code{TRUE}).
#'
#' @return A data.frame with one row per reference date and lag window, containing:
#'   \itemize{
#'     \item \code{date}: reference date
#'     \item \code{start}, \code{end}: start (inclusive) and end (exclusive) of
#'       the aggregation interval
#'     \item \code{lag_start}, \code{lag_end}: lag indices defining the interval
#'     \item One column per combination of variable and aggregation function
#'       (e.g. \code{rain_mean}, \code{temperature_sum})
#'   }
#'
#' @details
#' Intervals are defined as left-closed and right-open
#' (\code{[start, end)}). An interval is considered truncated if it extends
#' beyond the temporal bounds of the input time series. An interval is considered
#' missing if no observations fall within it.
#'
#' Console messages are printed to inform the user of reference dates for which
#' missing data or truncated intervals occurred.
#'
#' @examples
#' sampling_dates <- unique(albopictusMPL2023$date)
#'
#' met_agg <- aggregate_lagged_intervals(
#' data       = meteoMPL2023,
#' date_col   = "date",
#' value_cols = c("rain_sum", "temp_mean"),
#' d          = sampling_dates,
#' i          = 7,
#' m          = 8
#' )
#'
#' head(met_agg)
#'
#' @export
aggregate_lagged_intervals <- function(data,date_col,value_cols,d,
																			 i = 1, # integer, in days (7 for weekly intervals, 14 for fortnight, 30 for months...)
																			 m,
																			 funs = list(mean = mean,
																			 						min = min,
																			 						max  = max,
																			 						sum  = sum),
																			 na.rm = TRUE) {


	# Sanity checks

	stopifnot(inherits(data[[date_col]], c("Date", "POSIXct")))
	stopifnot(all(value_cols %in% names(data)))
	stopifnot(all(sapply(data[value_cols], is.numeric)))
	stopifnot(is.list(funs), !is.null(names(funs)))
	stopifnot(m >= 1, i >= 1)
	stopifnot(i == as.integer(i)) # check if its an integer

	d <- as.Date(d)


	# Time unit handling

	step <- i

	# Time series limits
	ts_min <- min(as.Date(data[[date_col]]), na.rm = TRUE)
	ts_max <- max(as.Date(data[[date_col]]), na.rm = TRUE)

	# Containers
	results <- list()
	missing_dates   <- as.Date(character(0))
	truncated_dates <- as.Date(character(0))


	# Loop over reference dates

	for (dd in d) {

		intervals <- list()
		k <- 1

		for (end_lag in 1:m) {
			end_date <- dd - (end_lag - 1) * step

			for (start_lag in end_lag:m) {
				start_date <- dd - start_lag * step

				intervals[[k]] <- data.frame(
					date     = dd,
					start     = start_date,
					end       = end_date,
					lag_start = start_lag,
					lag_end   = end_lag,
					stringsAsFactors = FALSE
				)
				k <- k + 1
			}
		}

		intervals <- do.call(rbind, intervals)

		has_missing   <- FALSE
		has_truncated <- FALSE


		# Aggregation

		agg_list <- lapply(seq_len(nrow(intervals)), function(j) {

			# Check if interval exceeds time series bounds
			if (intervals$start[j] < ts_min || intervals$end[j] > ts_max) {
				has_truncated <<- TRUE
			}

			idx <- data[[date_col]] >= intervals$start[j] &
				data[[date_col]] <  intervals$end[j]

			out <- c()

			for (v in value_cols) {

				x <- data[[v]][idx]

				if (length(x) == 0) {
					has_missing <<- TRUE
					vals <- rep(NA_real_, length(funs))
				} else {
					vals <- sapply(funs, function(f)
						do.call(f, list(x, na.rm = na.rm))
					)
				}

				names(vals) <- paste(v, names(funs), sep = "_")
				out <- c(out, vals)
			}

			out
		})

		agg_df <- as.data.frame(do.call(rbind, agg_list))
		rownames(agg_df) <- NULL

		results[[length(results) + 1]] <- cbind(intervals, agg_df)

		if (has_missing) {
			missing_dates <- c(missing_dates, dd)
		}

		if (has_truncated) {
			truncated_dates <- c(truncated_dates, dd)
		}
	}


	# Final assembly

	out <- do.call(rbind, results)

	out$date <- as.Date(out$date, origin = "1970-01-01")
	out$start <- as.Date(out$start, origin = "1970-01-01")
	out$end   <- as.Date(out$end,   origin = "1970-01-01")
	#out$interval <- out$end - out$start # might be used for plots


	# Console messages

	if (length(missing_dates) > 0) {
		message(
			"Missing data: at least one interval contained no observations for reference date(s): ",
			paste(format(unique(missing_dates)), collapse = ", ")
		)
	}

	if (length(truncated_dates) > 0) {
		message(
			"Truncated intervals: some statistics computed on partial intervals (outside time series bounds) for reference date(s): ",
			paste(format(unique(truncated_dates)), collapse = ", ")
		)
	}

	out
}

#' Run a complete ecoXCorr analysis: aggregation + lagged modelling
#'
#' This wrapper function combines \code{\link{aggregate_lagged_intervals}}
#' and \code{\link{fit_models_by_lag}} into a single workflow. It aggregates
#' environmental predictors over multiple lag windows relative to sampling
#' dates, merges them with a response dataset, and fits regression models
#' separately for each lag window.
#'
#' @param meteo_data Data frame containing meteorological time series.
#' @param response_data Data frame containing the response variable and
#'   sampling dates.
#' @param date_col_meteo Name of the date column in \code{meteo_data}.
#' @param date_col_resp Name of the date column in \code{response_data}.
#' @param value_cols Name of one meteorological variables to aggregate.
#' @param response Name of the response variable.
#' @param agg_fun Name (character string) of the aggregation function. Function must accept a numeric vector as first argument. Default to \code{"mean"}.
#' @param lag_unit Length of the base lag interval (in days).
#' @param max_lag Maximum number of lag intervals.
#' @param random Optional random-effects structure (passed to
#'   \code{fit_models_by_lag}).
#' @param family Model family (GLM or glmmTMB).
#' @param na.rm Logical indicating whether NA values are removed before aggregation.
#' @param track Logical; if TRUE, prints lag windows during model fitting.
#' @param ... Additional arguments passed to the model fitting function.
#'
#' @inherit fit_models_by_lag details
#' @inherit fit_models_by_lag return
#' @inherit fit_models_by_lag references
#'
#'
#' @examples
#' res_glmm <- ecoXCorr(
#' meteo_data    = meteoMPL2023,
#' response_data = albopictusMPL2023,
#' date_col_meteo   = "date",
#' date_col_resp    = "date",
#' value_cols    = "rain_sum",
#' agg_fun       = "sum",
#' response      = "individualCount",
#' lag_unit      = 7,
#' max_lag       = 8,
#' random        = "(1|area/trap)",
#' family        = "nbinom2"
#' )
#'
#' head(res_glmm)
#'
#' @export
ecoXCorr <- function(
    meteo_data,
    response_data,
    date_col_meteo   = "date",
    date_col_resp    = "date",
    value_cols,
    agg_fun = "mean",
    response,
    lag_unit = 1,
    max_lag,
    random = "",
    family = "gaussian",
    na.rm = TRUE,
    track = FALSE,
    ...
) {

  ## --- Checks
  stopifnot(length(value_cols)==1)
  stopifnot(length(agg_fun)==1 & is.character(agg_fun)==T)

  # create argument funs for `aggregate_lagged_intervals`
  funs <- list()
  funs[[1]] <- match.fun(agg_fun)
  names(funs) <- agg_fun


  ## --- Reference dates
  sampling_dates <- unique(response_data[[date_col_resp]])

  ## --- Aggregate lagged predictors
  met_agg <- aggregate_lagged_intervals(
    data       = meteo_data,
    date_col   = date_col_meteo,
    value_cols = value_cols,
    d          = sampling_dates,
    i          = lag_unit,
    m          = max_lag,
    funs       = funs,
    na.rm      = na.rm
  )

  # predictors names
  predictors <- paste0(value_cols,"_", agg_fun)

  stopifnot(predictors %in% names(met_agg))


  ## --- Merge with response data
  merged_data <- merge(
    met_agg,
    response_data,
    by.x = "date",
    by.y = date_col_resp,
    all  = TRUE
  )

  ## --- Fit models by lag window
  res <- fit_models_by_lag(
    data       = merged_data,
    response   = response,
    predictors = predictors,
    random     = random,
    family     = family,
    track      = track,
    ...
  )

  res
}

#' Launch the ecoXCorr Shiny application
#'
#' This function launches an interactive Shiny application allowing users
#' to run a complete ecoXCorr workflow (aggregation, lagged modelling and
#' visualisation) using either example datasets included in the package
#' or user-provided data.
#'
#' @export
ecoXCorrApp <- function() {
  runApp(
    system.file("shiny/app.R", package = "ecoXCorr"),
    launch.browser = TRUE
  )
}
