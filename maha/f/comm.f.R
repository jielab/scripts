if (!interactive() && grDevices::dev.cur() == 1L) {
	grDevices::pdf(NULL)
}

maha_split_steps <- function(x) {
	if (is.null(x) || length(x) == 0 || !nzchar(x[1])) return(character())
	unique(trimws(unlist(strsplit(x[1], "[,;[:space:]]+"))))
}

maha_expand_steps <- function(steps, available_steps, dependencies = list()) {
	if (!length(steps) || "all" %in% steps) return(available_steps)
	out <- character()
	add_step <- function(step) {
		for (dep in dependencies[[step]] %||% character()) add_step(dep)
		if (!step %in% out) out <<- c(out, step)
	}
	for (step in steps) add_step(step)
	available_steps[available_steps %in% out]
}

maha_configure_steps <- function(available_steps, dependencies = list()) {
	args <- commandArgs(trailingOnly = TRUE)
	step_arg <- grep("^--steps=", args, value = TRUE)
	steps_raw <- if (length(step_arg)) sub("^--steps=", "", step_arg[[length(step_arg)]]) else Sys.getenv("MAHA_STEPS", unset = "")
	requested <- maha_split_steps(steps_raw)
	if (!length(requested)) requested <- "all"
	if (any(args %in% c("--list-steps", "list-steps")) || any(requested == "list")) {
		cat("Available steps:\n", paste0("  - ", available_steps, collapse = "\n"), "\n", sep = "")
		if (interactive()) stop("Step list printed.", call. = FALSE)
		q(save = "no", status = 0)
	}
	unknown <- setdiff(requested, c("all", available_steps))
	if (length(unknown)) stop("Unknown MAHA step(s): ", paste(unknown, collapse = ", "), call. = FALSE)
	selected <- maha_expand_steps(requested, available_steps, dependencies)
	assign(".maha_selected_steps", selected, envir = .GlobalEnv)
	assign(".maha_requested_steps", requested, envir = .GlobalEnv)
	message("MAHA steps: ", paste(selected, collapse = ", "))
	invisible(selected)
}

maha_should_run_step <- function(step) {
	selected <- get0(".maha_selected_steps", envir = .GlobalEnv, ifnotfound = character())
	!length(selected) || step %in% selected
}

maha_run_step <- function(step, expr) {
	if (maha_should_run_step(step)) {
		eval.parent(substitute(expr))
	} else {
		message("Skipping MAHA step: ", step)
	}
	invisible(NULL)
}

maha_run_step_cached <- function(step, outputs, expr, min_bytes = 1) {
	if (!maha_should_run_step(step)) {
		message("Skipping MAHA step: ", step)
		return(invisible(FALSE))
	}
	output_paths <- eval.parent(substitute(outputs))
	force_steps <- maha_split_steps(Sys.getenv("MAHA_FORCE_STEPS", unset = ""))
	info <- suppressWarnings(file.info(output_paths))
	complete <- length(output_paths) > 0 && all(!is.na(info$size) & info$size >= min_bytes)
	if (complete && !step %in% force_steps && !"all" %in% force_steps) {
		message("Skipping completed MAHA step: ", step,
			" (verified outputs: ", paste(basename(output_paths), collapse = ", "), ")")
		return(invisible(FALSE))
	}
	eval.parent(substitute(expr))
	info <- suppressWarnings(file.info(output_paths))
	if (!all(!is.na(info$size) & info$size >= min_bytes)) {
		stop("MAHA step ", step, " finished without all required outputs: ",
			paste(output_paths[is.na(info$size) | info$size < min_bytes], collapse = ", "), call. = FALSE)
	}
	message("Completed and cached MAHA step: ", step)
	invisible(TRUE)
}

ap_tost_one <- function(beta, se, margin_hr = 1.05, alpha = 0.05) {
	margin <- log(margin_hr)
	if (!is.finite(beta) || !is.finite(se) || se <= 0) {
		return(tibble::tibble(
			margin_hr = margin_hr,
			CI90_low_HR = NA_real_, CI90_high_HR = NA_real_,
			TOST_p = NA_real_, equivalent = NA, Equivalence_5pct = NA
		))
	}
	z <- stats::qnorm(1 - alpha)
	ci_low <- beta - z * se
	ci_high <- beta + z * se
	p1 <- stats::pnorm((beta - (-margin)) / se, lower.tail = FALSE)
	p2 <- stats::pnorm((margin - beta) / se, lower.tail = FALSE)
	equivalent <- ci_low > -margin & ci_high < margin
	tibble::tibble(
		margin_hr = margin_hr,
		CI90_low_HR = exp(ci_low),
		CI90_high_HR = exp(ci_high),
		TOST_p = max(p1, p2),
		equivalent = equivalent,
		Equivalence_5pct = equivalent
	)
}

ap_bf01_one <- function(beta, se, prior_95_hr = 1.10) {
	if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NA_real_)
	tau <- log(prior_95_hr) / 1.96
	bf10 <- sqrt(se^2 / (se^2 + tau^2)) * exp((beta^2 * tau^2) / (2 * se^2 * (se^2 + tau^2)))
	1 / bf10
}

ap_power_one <- function(se, margin_hr = 1.05, alpha = 0.05, power = 0.80) {
	if (!is.finite(se) || se <= 0) {
		return(tibble::tibble(
			margin_hr = margin_hr,
			MDE_logHR = NA_real_,
			MDE_HR_lower = NA_real_,
			MDE_HR_upper = NA_real_,
			powered_for_margin = NA,
			Powered_5pct = NA
		))
	}
	z_alpha <- stats::qnorm(1 - alpha / 2)
	z_power <- stats::qnorm(power)
	mde <- (z_alpha + z_power) * se
	powered_for_margin <- exp(mde) <= margin_hr
	tibble::tibble(
		margin_hr = margin_hr,
		MDE_logHR = mde,
		MDE_HR_lower = exp(-mde),
		MDE_HR_upper = exp(mde),
		powered_for_margin = powered_for_margin,
		Powered_5pct = powered_for_margin
	)
}

ap_ci_lab <- function(hr, lo, hi) sprintf("%.2f (%.2f, %.2f)", hr, lo, hi)
ap_p_lab <- function(p) ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
ap_bf_lab <- function(x) ifelse(is.na(x), NA_character_, ifelse(x >= 100, ">100", sprintf("%.2f", x)))

ap_theme <- function(base_size = 14) {
	ggplot2::theme_classic(base_size = base_size) +
		ggplot2::theme(
			plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = base_size * 1.05),
			axis.text = ggplot2::element_text(color = "black"),
			axis.title = ggplot2::element_text(face = "bold"),
			panel.grid.major = ggplot2::element_blank(),
			panel.grid.minor = ggplot2::element_blank(),
			legend.position = "bottom",
			legend.title = ggplot2::element_blank()
		)
}
