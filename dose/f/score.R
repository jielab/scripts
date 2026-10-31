#!/usr/bin/env Rscript
# score: build, distribution, stroke, sensitivity, plot, outcomes, stroke_comparison
# Launch with dose.sh; source this file to register definitions for reuse.
if (!exists("dose_step_code", mode="function")) {
    source(file.path(Sys.getenv("DOSE_SCRIPT_DIR", if (.Platform$OS.type == "windows")
        "D:/scripts/dose/f" else "/mnt/d/scripts/dose/f"), "0.common.R"), local=TRUE)
}

# Read completed pattern jobs in the parent before starting any workers. New
# jobs persist their results and return only an acknowledgement. Unix workers
# use local process pipes, avoiding the TCP connections that can stall in WSL.
dose_cached_pattern_jobs <- function(jobs, load_result, run_job, workers=3L, replace=FALSE) {
    if (length(workers) != 1L || is.na(workers) || !is.finite(workers) ||
        workers < 1L || workers != as.integer(workers))
        stop("DOSE_ASSOC_WORKERS must be a positive integer.")
    results <- if (replace) vector("list", length(jobs)) else lapply(jobs, load_result)
    pending <- which(vapply(results, is.null, logical(1)))
    message("Pattern checkpoints: reusing ", length(jobs)-length(pending), "/", length(jobs),
            "; pending ", length(pending), ".")
    if (!length(pending)) return(results)
    workers <- min(as.integer(workers), length(pending))
    run_pending <- function(i) {
        run_job(jobs[[i]])
        TRUE
    }
    if (.Platform$OS.type == "unix" && workers > 1L) {
        status <- parallel::mclapply(pending, run_pending, mc.cores=workers,
                                    mc.preschedule=FALSE, mc.set.seed=FALSE)
    } else {
        # A serial fallback also avoids socket transport on non-Unix systems.
        status <- lapply(pending, run_pending)
    }
    if (!all(vapply(status, identical, logical(1), TRUE)))
        stop("A pattern worker failed; completed child checkpoints are preserved for resumption.")
    results[pending] <- lapply(jobs[pending], load_result)
    if (any(vapply(results, is.null, logical(1))))
        stop("A pattern worker did not save a matching result checkpoint.")
    results
}

# ===========================================================================
# score / build
# ===========================================================================
dose_register_step("score", "build", expression({
# ======================================================================
# score.R
# Apply the user-specified DOSE rules to RAW intake in servings/day.
#
# DOSE11: eleven separate binary food components, total 0-11.
#         Apples/Pears > 0 qualifies, with no upper limit.
# DOSE9 : nine binary components, total 0-9.
#         Fruit = Dried fruit > 0 OR Apples/Pears > 0.
#         Vegetables = Tomatoes > 0 OR (0 < Green leafy vegetables <= 1.56).
#         Test each food separately, then OR; each merged component <= 1.
#         Never add intakes before applying a merged component's rules.
# DOSE9half: same food criteria as DOSE9, but each paired domain is the
#           mean of its two binary criteria: neither=0, one=0.5, both=1.
# All versions: the same eleven food criteria, including wine and whole grains.
# No dairy, margarine, water or refined-grain component.
# V3 maximum is 9 (seven independent foods plus two paired domains).
# The thresholds below are fixed research specifications, not re-estimated.
#
# Input : ukb/02foodwas/dat.group.rds (foodwas.R)
# Output: ukb/04score/dose11_score.rds
#         ukb/04score/dose9_score.rds
#         ukb/04score/dose9half_score.rds (V3, half-point steps)
#         ukb/04score/dose_score.rds (PRIMARY: dose = dose9half)
#         ukb/04score/DOSE_scoring_rules.csv (component rules)
#         ukb/04score/DOSE_food_scoring_rules.csv (food criteria)
#         ukb/04score/DOSE_score_distribution.csv (all three versions)
#         ukb/04score/figures/Fig_DOSE{11,9,9HALF}.pdf
#
# Base R and the standard grid package only; no disease models are run here.
# ======================================================================

DIR_PROJECT <- dose_output_dir()
DIR_FOOD <- file.path(DIR_PROJECT, "ukb/02foodwas")
DIR_OUT <- file.path(DIR_PROJECT, "ukb/04score")
DIR_FIG <- file.path(DIR_OUT, "figures")
FILE_INPUT <- file.path(DIR_FOOD, "dat.group.rds")

# Display definitions for exported rules and figures. Keep table order.
# The directly readable ifelse() scoring calculation is in calculate_dose().
dose_food_rules <- function(version=c("dose9half", "dose11", "dose9")) {
    version <- match.arg(version)
    rules <- data.frame(
        Food=c("Apples_Pears", "Dried_fruit", "Legumes", "Olive_oil", "Tomatoes",
               "Eggs", "Green_leafy_vegetables", "Whole_grains", "Wine",
               "Processed_meats", "Sweetened_beverages"),
        Label=c("Apples and pears", "Dried fruits", "Legumes", "Olive oil", "Tomatoes",
                "Eggs", "Green leafy vegetables", "Whole grains", "Wine",
                "Processed meats", "Sweetened beverages"),
        Type=c(rep("adequacy", 5), rep("moderation", 5), "restriction"),
        Lower=c(rep(0, 7), 1, 0, 0, 0),
        Upper=c(rep(Inf, 5), 0.5, 1.56, 3, 2, 1, 0),
        Lower_inclusive=c(rep(FALSE, 10), TRUE),
        Upper_inclusive=TRUE,
        stringsAsFactors=FALSE
    )
    rules$Component <- paste0(rules$Food, ".", version)
    rules$Component_label <- rules$Label
    if (version != "dose11") {
        fruit <- rules$Food %in% c("Dried_fruit", "Apples_Pears")
        vegetables <- rules$Food %in% c("Tomatoes", "Green_leafy_vegetables")
        rules$Component[fruit] <- paste0("Fruit.", version)
        rules$Component_label[fruit] <- "Dried fruits, apples, and pears"
        rules$Component[vegetables] <- paste0("Vegetables.", version)
        rules$Component_label[vegetables] <- "Tomatoes and green leafy vegetables"
    }
    rules$Aggregation <- if (version == "dose9half") "mean" else if (version == "dose9") "OR" else "separate"
    rules$Food_weight <- if (version == "dose9half")
        1 / as.numeric(table(rules$Component)[rules$Component]) else 1
    rules$Version <- version
    rules$Score_variable <- paste0("diet.", version, ".sum")
    rules$Unit <- "servings/day"
    rules$Total_maximum <- if (version == "dose11") 11L else 9L
    rules
}

# Short interval labels follow the scoring table; >0-2 means >0 and <=2.
dose_interval_labels <- function(rule) {
    fmt <- function(x) format(x, trim=TRUE, scientific=FALSE)
    if (rule$Type == "restriction") {
        return(list(Low="0", Middle=">0", High=NA_character_,
                    Low_score=1L, Middle_score=0L, Zero=">0", One="0"))
    }
    if (rule$Type == "moderation") {
        lower <- fmt(rule$Lower)
        upper <- fmt(rule$Upper)
        low <- if (rule$Lower == 0) "0" else paste0("0-", lower)
        middle <- paste0(">", lower, "-", upper)
        return(list(Low=low, Middle=middle, High=paste0(">", upper),
                    Low_score=0L, Middle_score=1L,
                    Zero=paste0(low, " OR >", upper), One=middle))
    }
    lower <- fmt(rule$Lower)
    zero <- if (rule$Lower == 0) "0" else paste0("<", lower)
    one <- paste0(if (rule$Lower_inclusive) ">=" else ">", lower)
    list(Low=zero, Middle=one, High=NA_character_,
         Low_score=0L, Middle_score=1L, Zero=zero, One=one)
}

# Simple scoring: ifelse() for each food, then ordinary addition.
# Base R only. Missing criteria stay missing; totals never use na.rm=TRUE.
# OR is known to be 1 when one member qualifies, even if the other is missing.
calculate_dose <- function(dat) {
    foods <- dose_food_rules("dose11")$Food
    required <- c("eid", foods)
    missing_cols <- setdiff(required, names(dat))
    if (length(missing_cols)) {
        stop("Missing required columns: ", paste(missing_cols, collapse=", "))
    }
    if (!nrow(dat)) stop("The scoring input contains no participants.")
    ids <- as.character(dat$eid)
    if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
        stop("eid must be non-missing and unique: one row per participant.")
    }
    unit <- attr(dat, "intake_unit", exact=TRUE)
    if (!is.null(unit) && !identical(unit, "servings/day")) {
        stop("DOSE requires raw servings/day; input intake_unit is ", unit)
    }
    for (food in foods) {
        x <- dat[[food]]
        if (!is.numeric(x)) stop("Expected numeric raw intake for: ", food)
        if (any(is.infinite(x) | x < 0, na.rm=TRUE)) {
            stop("Negative or infinite intake found in: ", food)
        }
    }
    with(dat, {
        apples <- ifelse(Apples_Pears > 0, 1L, 0L)
        dried <- ifelse(Dried_fruit > 0, 1L, 0L)
        legumes <- ifelse(Legumes > 0, 1L, 0L)
        olive <- ifelse(Olive_oil > 0, 1L, 0L)
        tomatoes <- ifelse(Tomatoes > 0, 1L, 0L)
        eggs <- ifelse(Eggs > 0 & Eggs <= 0.5, 1L, 0L)
        leafy <- ifelse(Green_leafy_vegetables > 0 & Green_leafy_vegetables <= 1.56, 1L, 0L)
        whole <- ifelse(Whole_grains > 1 & Whole_grains <= 3, 1L, 0L)
        wine <- ifelse(Wine > 0 & Wine <= 2, 1L, 0L)
        processed <- ifelse(Processed_meats > 0 & Processed_meats <= 1, 1L, 0L)
        beverages <- ifelse(Sweetened_beverages == 0, 1L, 0L)

        # OR is applied to the two binary criteria, never to combined intakes.
        fruit_or <- ifelse(apples == 1 | dried == 1, 1L, 0L)
        vegetables_or <- ifelse(tomatoes == 1 | leafy == 1, 1L, 0L)
        fruit_half <- (apples + dried) / 2
        vegetables_half <- (tomatoes + leafy) / 2
        other7 <- legumes + olive + eggs + whole + wine + processed + beverages

        dose11 <- other7 + apples + dried + tomatoes + leafy
        dose9 <- other7 + fruit_or + vegetables_or
        dose9half <- other7 + fruit_half + vegetables_half
        data.frame(eid=as.character(eid), dose11, dose9, dose9half, dose=dose9half,
                   apples, dried, legumes, olive, tomatoes, eggs, leafy, whole,
                   wine, processed, beverages,
                   fruit_or, vegetables_or, fruit_half, vegetables_half)
    })
}

# Export each alternative with its own components; the arithmetic is above.
score_dose <- function(dat, version=c("dose9half", "dose11", "dose9"), calculated=NULL) {
    version <- match.arg(version)
    x <- if (is.null(calculated)) calculate_dose(dat) else calculated
    rules <- dose_food_rules(version)
    columns <- switch(version,
        dose11=c("apples","dried","legumes","olive","tomatoes","eggs","leafy","whole","wine","processed","beverages"),
        dose9=c("fruit_or","legumes","olive","vegetables_or","eggs","whole","wine","processed","beverages"),
        dose9half=c("fruit_half","legumes","olive","vegetables_half","eggs","whole","wine","processed","beverages"))
    component_names <- unique(rules$Component)
    components <- setNames(x[columns], component_names)
    result <- data.frame(eid=x$eid, total=as.numeric(x[[version]]),
                         components, check.names=FALSE)
    names(result)[2] <- rules$Score_variable[1]
    attr(result, "score_version") <- version
    attr(result, "score_intake_unit") <- "servings/day"
    attr(result, "score_contains_alcohol") <- any(rules$Food %in% c("Wine", "Beer", "Alcohol"))
    attr(result, "score_theoretical_range") <- c(0, rules$Total_maximum[1])
    attr(result, "component_columns") <- component_names
    attr(result, "food_scoring_rules") <- rules
    step <- if (version == "dose9half") 0.5 else 1
    attr(result, "score_step") <- step
    attr(result, "component_steps") <- setNames(vapply(component_names,
        function(nm) min(rules$Food_weight[rules$Component == nm]), numeric(1)), component_names)
    stopifnot(ncol(components) == rules$Total_maximum[1],
              isTRUE(all.equal(result[[2]], rowSums(components), check.attributes=FALSE)),
              all(is.na(result[[2]]) | result[[2]] %in% seq(0, rules$Total_maximum[1], by=step)))
    result
}

# Stable downstream interface. Only this scoring script selects the primary
# version; disease scripts read dose and its metadata, never a versioned name.
make_primary_dose <- function(calculated, scores) {
    primary <- scores$dose9half
    old_names <- attr(primary, "component_columns")
    main_components <- primary[old_names]
    names(main_components) <- sub("[.]dose9half$", ".dose", old_names)
    dat.score <- with(calculated, data.frame(
        eid, dose11, dose9, dose9half, dose=dose9half,
        diet.dose.sum=dose9half, diet.dose11.sum=dose11,
        diet.dose9.sum=dose9, diet.dose9half.sum=dose9half,
        main_components, check.names=FALSE))
    attr(dat.score, "score_version") <- "dose"
    attr(dat.score, "primary_scoring_method") <- "dose9half"
    for (nm in c("score_intake_unit", "score_contains_alcohol", "score_theoretical_range", "score_step"))
        attr(dat.score, nm) <- attr(primary, nm)
    attr(dat.score, "component_columns") <- names(main_components)
    attr(dat.score, "component_steps") <- setNames(attr(primary, "component_steps"), names(main_components))
    rules <- attr(primary, "food_scoring_rules")
    rules$Component <- sub("[.]dose9half$", ".dose", rules$Component)
    rules$Score_variable <- "dose"
    attr(dat.score, "food_scoring_rules") <- rules
    stopifnot(identical(dat.score$dose, dat.score$dose9half),
              identical(dat.score$dose, dat.score$diet.dose.sum),
              isTRUE(all.equal(dat.score$dose, rowSums(main_components), check.attributes=FALSE)))
    dat.score
}

dose_component_rules <- function(rules) {
    result <- lapply(unique(rules$Component), function(component) {
        members <- rules[rules$Component == component, , drop=FALSE]
        labels <- lapply(seq_len(nrow(members)), function(i) dose_interval_labels(members[i, ]))
        combine <- function(field, operator) {
            text <- vapply(labels, function(x) x[[field]], character(1))
            if (nrow(members) == 1) return(text)
            paste(paste0(members$Label, " (", text, ")"), collapse=operator)
        }
        data.frame(Version=members$Version[1], Score_variable=members$Score_variable[1],
                   Component=component, Food=members$Component_label[1],
                   Type=if (nrow(members) > 1) "combined domain" else members$Type,
                   Aggregation=members$Aggregation[1],
                   Zero_points=combine("Zero", " AND "),
                   Half_point=if (members$Aggregation[1] == "mean" && nrow(members) == 2)
                       "Exactly one of the two food criteria is met" else NA_character_,
                   One_point=combine("One", if (members$Aggregation[1] == "mean") " AND " else " OR "),
                   Unit="servings/day", Maximum_points=1L,
                   Total_maximum=members$Total_maximum[1], stringsAsFactors=FALSE)
    })
    do.call(rbind, result)
}

dose_distribution <- function(scores) {
    version <- attr(scores, "score_version")
    maximum <- if (version == "dose11") 11L else 9L
    variable <- paste0("diet.", version, ".sum")
    x <- scores[[variable]]
    levels <- seq(0, maximum, by=if (version == "dose9half") 0.5 else 1)
    result <- data.frame(Version=version, Score_variable=variable,
                         Score=c(as.character(levels), "Missing"),
                         N=c(as.integer(table(factor(x, levels=levels))), sum(is.na(x))))
    result$Percent <- 100 * result$N / nrow(scores)
    result
}

# All figures share the scoring definitions above and the original blue palette.
# Merged members stay on separate food rows, with matching markers in DOSE9.
plot_dose_rules <- function(rules, filename) {
    C <- list(hi="#2F5D8A", lo="#D9E1EA")
    version <- rules$Version[1]
    labels <- rules$Label
    if (version != "dose11") {
        fruit <- rules$Food %in% c("Dried_fruit", "Apples_Pears")
        vegetables <- rules$Food %in% c("Tomatoes", "Green_leafy_vegetables")
        labels[fruit] <- paste0("*", labels[fruit])
        labels[vegetables] <- paste0("**", labels[vegetables])
    }
    fig_width <- 7.95
    if (isTRUE(capabilities("cairo"))) {
        grDevices::cairo_pdf(filename, width=fig_width, height=6.6, family="sans")
    } else {
        grDevices::pdf(filename, width=fig_width, height=6.6, family="Helvetica")
    }
    on.exit(bounded_dev_off(), add=TRUE)
    grid::grid.newpage()
    # Preserve the narrow intake blocks; reserve room for the longer food labels.
    x_left <- 2.85 / fig_width
    x_right <- 1 - 0.17 / fig_width
    x_width <- x_right - x_left
    grid::grid.text(paste0(toupper(version), " (0-", rules$Total_maximum[1], " points)"),
                    x=(x_left + x_right)/2, y=0.974,
                    gp=grid::gpar(fontsize=14, fontface="bold"))
    grid::grid.text("Consumption (servings per day)",
                    x=(x_left + x_right)/2, y=0.925,
                    gp=grid::gpar(fontsize=14))
    draw_block <- function(left, right, y, label, score) {
        grid::grid.rect(x=x_left + x_width*(left + right)/2, y=y,
                        width=x_width*(right - left), height=0.054,
                        gp=grid::gpar(fill=if (score == 1) C$hi else C$lo, col=NA))
        grid::grid.text(label, x=x_left + x_width*(left + right)/2, y=y,
                        gp=grid::gpar(fontsize=12,
                                      col=if (score == 1) "white" else "grey25"))
    }
    for (i in seq_len(nrow(rules))) {
        row <- dose_interval_labels(rules[i, ])
        y <- 0.865 - (i - 1)*0.060
        grid::grid.text(labels[i], x=x_left - 0.012, y=y, just="right",
                        gp=grid::gpar(fontsize=13))
        draw_block(0, 1/3, y, row$Low, row$Low_score)
        if (is.na(row$High)) {
            draw_block(1/3, 1, y, row$Middle, row$Middle_score)
            dividers <- 1/3
        } else {
            draw_block(1/3, 0.74, y, row$Middle, row$Middle_score)
            draw_block(0.74, 1, y, row$High, 0L)
            dividers <- c(1/3, 0.74)
        }
        grid::grid.segments(x0=x_left + x_width*dividers, x1=x_left + x_width*dividers,
                            y0=y - 0.027, y1=y + 0.027,
                            gp=grid::gpar(col="white", lwd=3))
    }
    legend_y <- 0.202
    legend_origin <- (x_left + x_right)/2 - 1.6/fig_width
    legend_x <- function(inches) legend_origin + inches/fig_width
    grid::grid.text(if (version == "dose9half") "Rule:" else "Score:", x=legend_x(0), y=legend_y, just="left", gp=grid::gpar(fontsize=12))
    grid::grid.rect(x=legend_x(0.9), y=legend_y, width=0.49/fig_width, height=0.040,
                    gp=grid::gpar(fill=C$lo, col=NA))
    grid::grid.text(if (version == "dose9half") "Not met" else "0 point", x=legend_x(1.25), y=legend_y, just="left", gp=grid::gpar(fontsize=12))
    grid::grid.rect(x=legend_x(2.35), y=legend_y, width=0.49/fig_width, height=0.040,
                    gp=grid::gpar(fill=C$hi, col=NA))
    grid::grid.text(if (version == "dose9half") "Met" else "1 point", x=legend_x(2.7), y=legend_y, just="left", gp=grid::gpar(fontsize=12))
    footnote <- if (version == "dose9half") {
        paste0("* Fruit and ** vegetables: assess each food criterion separately.\n",
               "Each paired domain: neither met = 0; one met = 0.5; both met = 1 point.\n",
               "Other seven foods: criterion met = 1 point. Block widths are schematic.")
    } else if (version == "dose9") {
        paste0("* Fruit: dried fruits >0 OR apples/pears >0.\n",
               "** Vegetables: tomatoes >0 OR green leafy vegetables >0-1.56.\n",
               "Each merged component: either or both criteria met = 1 point. Block widths are schematic.")
    } else {
        "Eleven separate components, each worth 0 or 1 point. Block widths are schematic."
    }
    grid::grid.text(footnote, x=0.025, y=0.095, just="left",
                    gp=grid::gpar(fontsize=9, col="grey25", lineheight=1.4))
}

run_dose_scoring <- function(input=FILE_INPUT, output_dir=DIR_OUT) {
    if (!file.exists(input)) {
        stop("Input not found: ", input, ". Run f/foodwas.R first.")
    }
    dat <- readRDS(input)
    calculated <- calculate_dose(dat)
    # Finish and validate all versions before writing any scoring outputs.
    versions <- c("dose11", "dose9", "dose9half")
    scores <- setNames(lapply(versions, function(v) score_dose(dat, v, calculated)), versions)
    stopifnot(identical(scores$dose11$eid, scores$dose9$eid),
              identical(scores$dose11$eid, scores$dose9half$eid))
    food_rules <- do.call(rbind, lapply(names(scores), dose_food_rules))
    component_rules <- do.call(rbind, lapply(names(scores), function(v) dose_component_rules(dose_food_rules(v))))
    distribution <- do.call(rbind, lapply(scores, dose_distribution))
    rownames(food_rules) <- rownames(component_rules) <- rownames(distribution) <- NULL
    food_rules$Zero_points <- vapply(seq_len(nrow(food_rules)), function(i) dose_interval_labels(food_rules[i, ])$Zero, character(1))
    food_rules$One_point <- vapply(seq_len(nrow(food_rules)), function(i) dose_interval_labels(food_rules[i, ])$One, character(1))
    food_rules$Criterion_met <- food_rules$One_point
    food_rules$Half_point <- NA_character_
    half_member <- food_rules$Version == "dose9half" & food_rules$Food_weight == .5
    food_rules$Half_point[half_member] <- food_rules$Criterion_met[half_member]
    food_rules$One_point[half_member] <- NA_character_
    food_rules$Scoring_note <- ifelse(food_rules$Version == "dose9half",
        "Add member contributions within each domain; maximum 1 per domain",
        ifelse(food_rules$Version == "dose9", "OR within paired domain; maximum 1 per domain", "Separate binary food component"))
    figure_dir <- file.path(output_dir, "figures")
    dir.create(figure_dir, recursive=TRUE, showWarnings=FALSE)
    for (version in names(scores)) {
        score_file <- file.path(output_dir, paste0(version, "_score.rds"))
        figure_file <- file.path(figure_dir, paste0("Fig_", toupper(version), ".pdf"))
        dose_save_rds(scores[[version]], score_file)
        plot_dose_rules(dose_food_rules(version), figure_file)
        cat(toupper(version), ": ", nrow(dat), " participants; theoretical total 0-",
            if (version == "dose11") 11L else 9L, ".\n", sep="")
        print(distribution[distribution$Version == version, c("Score", "N", "Percent")], row.names=FALSE)
        cat("Score variable: diet.", version, ".sum\nScores: ", score_file,
            "\nScoring figure: ", figure_file, "\n", sep="")
    }
    # Retain all alternatives for sensitivity work, but expose only dose as primary.
    dat.score <- make_primary_dose(calculated, scores)
    dose_save_rds(dat.score, file.path(output_dir, "dose_score.rds"))
    write.csv(component_rules, file.path(output_dir, "DOSE_scoring_rules.csv"), row.names=FALSE)
    write.csv(food_rules, file.path(output_dir, "DOSE_food_scoring_rules.csv"), row.names=FALSE)
    write.csv(distribution, file.path(output_dir, "DOSE_score_distribution.csv"), row.names=FALSE)
    definition <- data.frame(Version=versions, Label=c("V1: 11 separate foods", "V2: 9 components, OR", "V3: 9 components, half points"),
        Maximum=c(11,9,9), Step=c(1,1,0.5), Apple_upper=rep(Inf,3),
        Fruit_vegetable_aggregation=c("separate", "OR", "mean of two binary criteria"))
    definition$Primary <- versions == "dose9half"
    write.csv(definition, file.path(output_dir, "DOSE_version_definitions.csv"), row.names=FALSE)
    valid <- complete.cases(scores$dose9[[2]], scores$dose9half[[2]])
    difference <- scores$dose9[[2]][valid] - scores$dose9half[[2]][valid]
    stopifnot(all(difference %in% c(0,0.5,1)))
    audit <- data.frame(Metric=c("Participants", "Apples_Pears_above_2", "V2_V3_common_nonmissing", "V2_equals_V3", "V2_minus_V3_0.5", "V2_minus_V3_1"),
        N=c(nrow(dat), sum(dat$Apples_Pears > 2, na.rm=TRUE), sum(valid), sum(difference == 0), sum(difference == .5), sum(difference == 1)))
    write.csv(audit, file.path(output_dir, "DOSE_version_audit.csv"), row.names=FALSE)
    cat("PRIMARY: dose = dose9half; downstream exposure = dose; diet.dose.sum is a compatibility alias. All three totals are retained.\n")
    invisible(list(scores=scores, primary=dat.score, score_rules=component_rules,
                   food_rules=food_rules, score_distribution=distribution))
}

if (isTRUE(getOption("dose.score.run", TRUE))) {
    dose_result <- run_dose_scoring()
    dat.score11 <- dose_result$scores$dose11
    dat.score9 <- dose_result$scores$dose9
    dat.score9half <- dose_result$scores$dose9half
    dat.score <- dose_result$primary
    score_rules <- dose_result$score_rules
    score_distribution <- dose_result$score_distribution
}
}))

# ===========================================================================
# score / distribution
# ===========================================================================
dose_register_step("score", "distribution", expression({
# DOSE adherence: score distribution and mean food intake (Fig. 1e style).
# Reference: https://www.nature.com/articles/s41591-025-03891-5/figures/1
# https://github.com/DW-Group/Gene_Metabolites_MedDiet_Dementia/blob/main/Figures/Fig1e.R
# Run after score.R. Reads primary score metadata; does not rerun scoring.
# Both panels use participants with an observed primary DOSE score.
# Food means include zero intake and exclude missing values (not consumers only).
# Eleven original food intakes are displayed separately, even when two food
# criteria contribute to one paired scoring domain. Units: servings/day.
# Optional food-group colors: options(dose.score.plot.colour="food_groups").

run_dose_score_plot <- function(
    project_dir=dose_output_dir(),
    score_file=file.path(project_dir, "ukb/04score/dose_score.rds"),
    food_file=file.path(project_dir, "ukb/02foodwas/dat.group.rds"),
    output_dir=file.path(project_dir, "ukb/04score/figures"),
    colour=getOption("dose.score.plot.colour", "blue")) {
    colour <- match.arg(colour, c("blue", "food_groups"))
    for (pkg in c("ggplot2", "patchwork")) {
        if (!requireNamespace(pkg, quietly=TRUE)) stop("Please install package: ", pkg)
    }
    for (path in c(score_file, food_file)) {
        if (!file.exists(path)) stop("Input not found: ", path, "; run score.R first.")
    }
    scores <- readRDS(score_file)
    foods <- readRDS(food_file)
    if (!all(c("eid", "dose") %in% names(scores)) || !"eid" %in% names(foods))
        stop("Expected eid and primary dose in score input; eid in food input.")
    check_ids <- function(x) {
        ids <- as.character(x$eid)
        if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids))
            stop("Inputs must have unique, non-missing participant IDs.")
    }
    check_ids(scores)
    check_ids(foods)
    rules <- attr(scores, "food_scoring_rules", exact=TRUE)
    limits <- attr(scores, "score_theoretical_range", exact=TRUE)
    step <- attr(scores, "score_step", exact=TRUE)
    if (is.null(rules) || !all(c("Food", "Label") %in% names(rules)) ||
        length(limits) != 2L || any(!is.finite(limits)) || diff(limits) <= 0 ||
        length(step) != 1L || !is.finite(step) || step <= 0)
        stop("Missing or invalid score metadata; rerun score.R.")
    if (!identical(attr(scores, "score_intake_unit"), "servings/day"))
        stop("Expected servings/day in score metadata.")
    food_unit <- attr(foods, "intake_unit", exact=TRUE)
    if (!is.null(food_unit) && !identical(food_unit, "servings/day"))
        stop("Food input must contain raw servings/day.")
    if (anyDuplicated(rules$Food) || !all(rules$Food %in% names(foods)))
        stop("Missing or duplicated foods in score metadata/input.")
    index <- match(as.character(scores$eid), as.character(foods$eid))
    if (anyNA(index)) stop("Some scored participants have no matching food record.")
    if (!is.numeric(scores$dose)) stop("Primary DOSE must be numeric.")
    keep <- !is.na(scores$dose)
    x <- scores$dose[keep]
    if (!length(x) || any(!is.finite(x)) || any(x < limits[1] | x > limits[2]) ||
        any(abs((x - limits[1])/step - round((x - limits[1])/step)) > 1e-8))
        stop("Primary scores are empty or inconsistent with range/step metadata.")
    dat <- foods[index[keep], rules$Food, drop=FALSE]
    if (any(!vapply(dat, is.numeric, logical(1))) ||
        any(vapply(dat, function(z) any(is.infinite(z) | z < 0, na.rm=TRUE), logical(1))))
        stop("Food intakes must be numeric, nonnegative and finite when observed.")

    intake <- do.call(rbind, lapply(seq_len(nrow(rules)), function(i) {
        z <- dat[[rules$Food[i]]]
        n <- sum(!is.na(z))
        if (!n) stop("No observed intake: ", rules$Food[i])
        data.frame(Food=rules$Food[i], Label=rules$Label[i], N=n,
            Missing=sum(is.na(z)), Mean=mean(z, na.rm=TRUE),
            Median=median(z, na.rm=TRUE),
            Q1=unname(quantile(z, .25, na.rm=TRUE)),
            Q3=unname(quantile(z, .75, na.rm=TRUE)), Unit="servings/day")
    }))
    # Same colors identify foods in the same dietary group.
    groups <- c(Apples_Pears="Fruit", Dried_fruit="Fruit", Legumes="Legumes",
        Olive_oil="Oil", Tomatoes="Vegetables", Eggs="Eggs",
        Green_leafy_vegetables="Vegetables", Whole_grains="Whole grains",
        Wine="Wine", Processed_meats="Meat", Sweetened_beverages="Beverages")
    intake$Group <- unname(groups[intake$Food])
    if (anyNA(intake$Group)) stop("Add new food groups to the plotting palette.")
    palette <- c(Fruit="#7D8C3C", Vegetables="#7D8C3C", Legumes="#D89E5F",
        Oil="#FAD02E", Eggs="#D8A8B8", "Whole grains"="#2C7BB6",
        Wine="#A8A7A4", Meat="#9D5C63", Beverages="#69B3A2")
    if (colour == "blue") palette[] <- "#3D6388"
    # Stable order within ties; highest mean appears at the top.
    ord <- order(intake$Mean, seq_len(nrow(intake)))
    intake$Label <- factor(intake$Label, levels=intake$Label[ord])
    score_grid <- seq(limits[1], limits[2], by=step)
    distribution <- data.frame(Score=score_grid,
        N=as.integer(table(factor(x, levels=score_grid))))
    distribution$Percent <- 100 * distribution$N / length(x)
    distribution$Density <- distribution$N / (length(x) * step)
    fmt_n <- function(n) format(n, big.mark=",", scientific=FALSE, trim=TRUE)
    summary <- data.frame(N_input=nrow(scores), N_plotted=length(x),
        N_missing_score=sum(!keep), Mean=mean(x), SD=sd(x), Median=median(x),
        Q1=unname(quantile(x, .25)), Q3=unname(quantile(x, .75)),
        Observed_min=min(x), Observed_max=max(x),
        Theoretical_min=limits[1], Theoretical_max=limits[2], Step=step)

    p_score <- ggplot2::ggplot(distribution, ggplot2::aes(Score, Percent)) +
        ggplot2::geom_col(width=step, fill="#AFC4D5",
                          colour="white", linewidth=.25)
    # Smooth overlay is descriptive only: the DOSE score remains discrete.
    if (length(x) > 1L && length(unique(x)) > 1L) {
        # Large n makes default bandwidth too narrow for half-point scores.
        # At least one scoring step avoids artificial spikes at each lattice point.
        kd <- stats::density(x, bw=max(stats::bw.nrd0(x), step),
            from=limits[1], to=limits[2], n=512)
        # Density * bin width * 100 matches percentage per score step.
        p_score <- p_score + ggplot2::geom_line(
            data=data.frame(Score=kd$x, Percent=kd$y * step * 100),
            colour="#D46A2E", linewidth=.65)
    }
    p_score <- p_score +
        ggplot2::scale_x_continuous(breaks=seq(limits[1], limits[2], by=1),
            limits=limits + c(-step/2, step/2), expand=ggplot2::expansion(mult=.01)) +
        ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0, .12))) +
        ggplot2::labs(title="DOSE score", subtitle=paste0("n = ", fmt_n(length(x))),
            x="DOSE (points)", y="Participants (%)") +
        ggplot2::theme_classic(base_size=12) +
        ggplot2::theme(plot.title=ggplot2::element_text(face="bold", size=16),
            plot.subtitle=ggplot2::element_text(colour="grey35"),
            plot.margin=ggplot2::margin(12, 18, 12, 12))
    p_intake <- ggplot2::ggplot(intake, ggplot2::aes(Mean, Label, fill=Group)) +
        ggplot2::geom_col(width=.78, show.legend=FALSE) +
        ggplot2::geom_text(ggplot2::aes(label=sprintf("%.2f", Mean)),
            hjust=-.2, size=3.8, fontface="bold") +
        ggplot2::scale_fill_manual(values=palette) +
        ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0, .18))) +
        ggplot2::labs(title="DOSE food intake", subtitle="Mean intake (servings/day)",
            x=NULL, y=NULL) +
        ggplot2::theme_void(base_size=12) +
        ggplot2::theme(axis.text.y=ggplot2::element_text(colour="grey12", hjust=1,
                margin=ggplot2::margin(r=8)),
            plot.title=ggplot2::element_text(face="bold", size=16),
            plot.subtitle=ggplot2::element_text(colour="grey35", margin=ggplot2::margin(b=12)),
            plot.margin=ggplot2::margin(12, 18, 12, 12))
    caption <- paste0("Both panels: participants with an observed primary DOSE score (n = ",
        fmt_n(length(x)), "); missing score: ", fmt_n(sum(!keep)), ".\n",
        "Food means use all available intake within this cohort, including zero intake; missing values are excluded.\n",
        "The 11 food intakes are shown separately. The orange curve is a descriptive kernel smooth\n",
        "(bandwidth >= score step), scaled as density x score step x 100 to match Participants (%).")
    p_inset <- p_score +
        ggplot2::labs(title="DOSE score", x=NULL) +
        ggplot2::theme(plot.title=ggplot2::element_text(size=12, face="bold", hjust=.5),
            plot.subtitle=ggplot2::element_text(size=9, hjust=.5),
            axis.title.y=ggplot2::element_text(size=9),
            axis.text=ggplot2::element_text(size=9, colour="grey20"),
            axis.line=ggplot2::element_line(colour="grey45", linewidth=.3),
            axis.ticks=ggplot2::element_line(colour="grey45", linewidth=.3),
            plot.margin=ggplot2::margin(4, 4, 4, 4),
            plot.background=ggplot2::element_rect(fill="white", colour=NA))
    # Panel-relative coordinates place the inset beside the shorter food bars.
    # The full figure legend is exported separately for manuscript assembly.
    combined <- p_intake + patchwork::inset_element(p_inset,
        left=.46, bottom=.035, right=.995, top=.55, align_to="panel")
    dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
    save_plot <- function(plot, name, width, height) {
        pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else "pdf"
        bounded_ggsave(file.path(output_dir, paste0(name, ".pdf")), plot,
            device=pdf_device, width=width, height=height, bg="white")
        bounded_ggsave(file.path(output_dir, paste0(name, ".png")), plot,
            width=width, height=height, dpi=300, bg="white")
    }
    save_plot(combined, "Fig_DOSE_adherence", 8, 5.8)
    save_plot(p_score, "Fig_DOSE_score_histogram", 5.3, 4.5)
    save_plot(p_intake, "Fig_DOSE_food_mean_intake", 6.5, 5)
    writeLines(caption, file.path(output_dir, "Fig_DOSE_adherence_caption.txt"))
    write.csv(intake, file.path(output_dir, "DOSE_food_intake_summary.csv"), row.names=FALSE)
    write.csv(distribution, file.path(output_dir, "DOSE_plot_score_distribution.csv"), row.names=FALSE)
    write.csv(summary, file.path(output_dir, "DOSE_plot_cohort_summary.csv"), row.names=FALSE)
    message("DOSE plots saved to: ", output_dir, " (n = ", fmt_n(length(x)), ")")
    invisible(list(combined=combined, score=p_score, intake=p_intake,
        food_summary=intake, score_distribution=distribution, cohort_summary=summary))
}

if (isTRUE(getOption("dose.score.plot.run", TRUE))) {
    dose_score_plots <- run_dose_score_plot()
}
}))

# ===========================================================================
# score / stroke
# ===========================================================================
dose_register_step("score", "stroke", expression({
# ======================================================================
# DOSE + established dietary patterns: associations with incident stroke.
# Run after score.R, then run score.R for figures/tables.
# Adapted from foodwas.R and dose_OLD/scripts/03_pattern.R.
# Primary outcome/model: ischemic stroke (cvd_stroke_i), M2.
# MODERN grouping: cut(score, quantile(score), include.lowest=TRUE,
#                      right=TRUE); identical scores stay together.
# Reference: https://github.com/Happychrischen/Diet_ML_Dementia/blob/main/Main%20analysis/Pattern_Dementia.R
# The exposure is always dose; s3_1 selects its primary scoring method.
# Alternative scores remain in dose_score.rds but are not analysed here.
# Compare primary DOSE and 13 established patterns on the SAME cohort.
# Score range, step, components and alcohol handling come from RDS metadata.
# Common C-index comparisons use the same covariates as the DOSE M2 model.
# Exact score ties stay together, including half-point DOSE scores.
# Outputs: ukb/04score/pattern/ (earlier s1 results are preserved).
# Dependencies: base R and survival. No external helper scripts required.
# ======================================================================

library(survival)
DIR_PROJECT <- dose_output_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_SCORE <- file.path(DIR_PROJECT, "ukb/04score")
DIR_OUT <- file.path(DIR_SCORE, "pattern")
dir.create(DIR_OUT, recursive=TRUE, showWarnings=FALSE)
PRIMARY_OUTCOME <- "cvd_stroke_i"
MODELS <- c("M1", "M2", "M3")
DOSE_DIETS <- "dose"

input_files <- c(cohort=file.path(DIR_DATA, "dat.base.imp.rds"),
                 config=file.path(DIR_DATA, "analysis_config.rds"),
                 dose=file.path(DIR_SCORE, "dose_score.rds"))
if (any(!file.exists(input_files))) stop("Missing input: ", paste(input_files[!file.exists(input_files)], collapse=", "))
cfg <- readRDS(input_files[["config"]])
dat <- as.data.frame(readRDS(input_files[["cohort"]]))
dat$eid <- as.character(dat$eid)
if (anyNA(dat$eid) || anyDuplicated(dat$eid)) stop("Invalid or duplicated cohort eid.")
dose <- readRDS(input_files[["dose"]])
component_cols <- attr(dose, "component_columns", exact=TRUE)
dose_meta <- list(method=attr(dose, "primary_scoring_method"),
    range=attr(dose, "score_theoretical_range"), step=attr(dose, "score_step"),
    alcohol=attr(dose, "score_contains_alcohol"), component_steps=attr(dose, "component_steps"))
if (!identical(attr(dose, "score_version"), "dose") ||
    !identical(attr(dose, "score_intake_unit"), "servings/day") ||
    length(dose_meta$method) != 1L || is.na(dose_meta$method) ||
    length(dose_meta$range) != 2L || any(!is.finite(dose_meta$range)) || diff(dose_meta$range) <= 0 ||
    length(dose_meta$step) != 1L || !is.finite(dose_meta$step) || dose_meta$step <= 0 ||
    !is.logical(dose_meta$alcohol) || length(dose_meta$alcohol) != 1L || is.na(dose_meta$alcohol) ||
    !length(component_cols) || anyDuplicated(component_cols) ||
    !all(c("eid", "dose", component_cols) %in% names(dose)) ||
    !identical(names(dose_meta$component_steps), component_cols) ||
    any(!is.finite(dose_meta$component_steps) | dose_meta$component_steps <= 0))
    stop("Invalid primary DOSE metadata. Run the current score.R first.")
dose$eid <- as.character(dose$eid)
if (anyNA(dose$eid) || anyDuplicated(dose$eid) || !setequal(dat$eid,dose$eid)) stop("Score/cohort ID mismatch.")
valid_components <- vapply(component_cols, function(v) is.numeric(dose[[v]]) &&
    all(is.na(dose[[v]]) | dose[[v]] %in% seq(0, 1, by=dose_meta$component_steps[[v]])), logical(1))
if (!all(valid_components) || !is.numeric(dose$dose) ||
    !all(is.na(dose$dose) | dose$dose %in% seq(dose_meta$range[1], dose_meta$range[2], by=dose_meta$step)) ||
    !isTRUE(all.equal(dose$dose, rowSums(dose[component_cols]), check.attributes=FALSE)))
    stop("DOSE must match its component sum, range and step metadata.")
if ("diet.dose.sum" %in% names(dose) && !identical(dose$dose, dose$diet.dose.sum))
    stop("The compatibility alias diet.dose.sum differs from dose.")
dose$diet.dose.sum <- dose$dose
new_cols <- c("dose", "diet.dose.sum", component_cols)
dat[new_cols] <- dose[match(dat$eid,dose$eid), new_cols, drop=FALSE]

# Retain the Oxford WebQ MedDiet (medi24); exclude the baseline touchscreen
# MedDiet-original (medito), whose missing scores must not restrict this analysis.
# MAHA includes alcohol.
spec <- data.frame(
    diet=c("medi24", "mind", "ahei", "redii", "fds", "dash", "hpdi", "phdi",
           "digm", "hlcd", "hlfd", "modern", "dose", "medito", "maha"),
    label=c("MedDiet", "MIND", "AHEI", "rEDII", "FDS", "DASH", "hPDI", "PHDI",
            "DI-GM", "HLCD", "HLFD", "MODERN", "DOSE", "MedDiet-original", "MAHA"),
    alcohol_in_score=c(rep(TRUE, 5), rep(FALSE, 7), dose_meta$alcohol, FALSE, TRUE),
    include=TRUE, stringsAsFactors=FALSE
)
spec$include[spec$diet == "medito"] <- FALSE
spec$variable <- paste0("diet.", spec$diet, ".sum")
spec$variable[spec$diet == "dose"] <- "dose"
spec$scoring_method <- ifelse(spec$diet == "dose", dose_meta$method, spec$diet)
spec$minimum <- ifelse(spec$diet == "dose", dose_meta$range[1], NA_real_)
spec$maximum <- ifelse(spec$diet == "dose", dose_meta$range[2], NA_real_)
spec$step <- ifelse(spec$diet == "dose", dose_meta$step, NA_real_)
spec$available <- spec$variable %in% names(dat)
if (any(spec$include & !spec$available)) stop("Requested pattern missing: ", paste(spec$diet[spec$include & !spec$available], collapse=", "))
use <- spec[spec$include & spec$available, ]
Ys <- PRIMARY_OUTCOME
if (!all(Ys %in% cfg$Ys)) stop("Requested outcome is not in analysis_config.rds.")
out_map <- cfg$out_map
outcome_vars <- unlist(lapply(Ys, function(y) paste0(y, c(".t2e", ".Yt2e"))))
all_covs <- unique(c(cfg$covs.m1, cfg$covs.alcohol, cfg$covs.m2.add, cfg$covs.m3.add))
required <- unique(c("eid", outcome_vars, all_covs, use$variable, "diet.dose.sum", component_cols))
if (!all(required %in% names(dat))) stop("Missing analysis variables: ", paste(setdiff(required, names(dat)), collapse=", "))
dat <- dat[required]
rm(dose)
invisible(gc())

model_covs <- function(diet, model) {
    alc <- use$alcohol_in_score[match(diet, use$diet)]
    m1 <- unique(c(cfg$covs.m1, if (!alc) cfg$covs.alcohol))
    co <- switch(model, M1=m1, M2=unique(c(m1, cfg$covs.m2.add)),
           M3=unique(c(m1, cfg$covs.m2.add, cfg$covs.m3.add)))
    if (alc) setdiff(co, cfg$covs.alcohol) else co
}
if (dose_meta$alcohol) for (m in MODELS)
    stopifnot(!any(cfg$covs.alcohol %in% model_covs("dose", m)))

score_info <- list(); group_info <- list()
for (i in seq_len(nrow(use))) {
    diet <- use$diet[i]; x <- dat[[use$variable[i]]]
    if (!is.numeric(x) || any(is.infinite(x))) stop("Invalid score: ", diet)
    br <- unname(quantile(x, probs=c(0, .25, .5, .75, 1), na.rm=TRUE))
    # Do not use ntile() or silently split participants with equal scores.
    if (anyDuplicated(br)) stop("Fewer than four distinct quantile intervals for ", diet, "; inspect before grouping.")
    qvar <- paste0("diet.", diet, ".q")
    dat[[qvar]] <- cut(x, breaks=br, include.lowest=TRUE, right=TRUE, labels=paste0("Q", 1:4))
    contrasts(dat[[qvar]]) <- contr.treatment(levels(dat[[qvar]]), base=1)
    p1090 <- unname(quantile(x, c(.1, .9), na.rm=TRUE))
    sx <- sd(x, na.rm=TRUE)
    if (!is.finite(sx) || sx == 0 || diff(p1090) == 0) stop("Degenerate score scale: ", diet)
    dat[[paste0("diet.", diet, ".p1090")]] <- (x-p1090[1])/diff(p1090)
    dat[[paste0("diet.", diet, ".z")]] <- (x-mean(x, na.rm=TRUE))/sx
    score_info[[diet]] <- data.frame(diet=diet, label=use$label[i], N=sum(!is.na(x)),
                                     mean=mean(x, na.rm=TRUE), SD=sx, P10=p1090[1], P90=p1090[2],
                                     P10_P90_increment=diff(p1090), min=br[1], Q25=br[2], Q50=br[3], Q75=br[4], max=br[5])
    for (j in 1:4) {
        z <- x[which(dat[[qvar]] == paste0("Q", j))]
        group_info[[paste(diet, j)]] <- data.frame(diet=diet, label=use$label[i], category=paste0("Q", j),
            cut_lower=br[j], cut_upper=br[j+1], lower_inclusive=j == 1, upper_inclusive=TRUE,
            score_min=min(z), score_max=max(z), score_median=median(z), N=length(z), Percent=100*length(z)/sum(!is.na(x)))
    }
}
score_info <- do.call(rbind, score_info); rownames(score_info) <- NULL
group_info <- do.call(rbind, group_info); rownames(group_info) <- NULL
write.csv(spec, file.path(DIR_OUT, "pattern_spec.csv"), row.names=FALSE)
write.csv(score_info, file.path(DIR_OUT, "pattern_score_summary.csv"), row.names=FALSE)
write.csv(group_info, file.path(DIR_OUT, "pattern_quartile_definitions.csv"), row.names=FALSE)
dose_save_rds(dat, file.path(DIR_OUT, "dat.pattern.rds"))

# Common complete-case cohort within each outcome for every diet and model.
# Use only included scores for the common cohort. The excluded touchscreen
# MedDiet-original must not remove participants with complete WebQ scores.
cohort_qc <- list(); cohorts <- list()
for (Y in Ys) {
    t <- paste0(Y, ".t2e"); y <- paste0(Y, ".Yt2e")
    ok <- complete.cases(dat[c(t, y, all_covs, use$variable)]) & is.finite(dat[[t]]) & dat[[t]] > 0
    if (any(!dat[[y]][ok] %in% 0:1)) stop("Outcome is not binary: ", Y)
    cohorts[[Y]] <- which(ok)
    cohort_qc[[Y]] <- data.frame(outcome=Y, label=unname(out_map[Y]), input_N=nrow(dat),
        analysis_N=sum(ok), excluded_N=sum(!ok), Events=sum(dat[[y]][ok]), person_years=sum(dat[[t]][ok]))
}
write.csv(do.call(rbind, cohort_qc), file.path(DIR_OUT, "analysis_cohort_QC.csv"), row.names=FALSE)

fit_warnings <- character(0)
fit_cox <- function(d, terms, context) {
    fit <- withCallingHandlers(
        coxph(reformulate(terms, response="Surv(time, event)"), data=d, ties="efron", singular.ok=FALSE),
        warning=function(w) {
            fit_warnings <<- c(fit_warnings, paste(context, conditionMessage(w), sep=": "))
            invokeRestart("muffleWarning")
        }
    )
    if (any(!is.finite(coef(fit)))) stop("Non-finite Cox coefficient: ", context)
    fit
}
extract_effect <- function(fit, term, increment=1) {
    b <- unname(coef(fit)[term]); se <- sqrt(vcov(fit)[term, term])
    data.frame(HR=exp(b*increment), LCI=exp((b-1.959964*se)*increment),
               HCI=exp((b+1.959964*se)*increment), P=2*pnorm(-abs(b/se)))
}
ph_test <- function(fit, term) {
    tryCatch(unname(cox.zph(fit, terms=TRUE)$table[term, "p"]), error=function(e) {
        fit_warnings <<- c(fit_warnings, paste("PH test", term, conditionMessage(e)))
        NA_real_
    })
}
lr_p <- function(loglik_full, loglik_null, df) pchisq(max(0, 2*(loglik_full-loglik_null)), df=df, lower.tail=FALSE)
null_cache <- new.env(parent=emptyenv())
get_null <- function(d, covs, key) {
    if (!exists(key, envir=null_cache, inherits=FALSE)) {
        f <- fit_cox(d, covs, paste(key, "baseline"))
        assign(key, as.numeric(logLik(f)), envir=null_cache)
    }
    get(key, envir=null_cache, inherits=FALSE)
}

run_outcome <- function(Y, diet_subset=use$diet) {
    cache_suffix <- paste(diet_subset, collapse="_")
    cache_file <- file.path(DIR_OUT, paste0("checkpoint_", Y, "_", cache_suffix, ".rds"))
    if (!dose_replace() && file.exists(cache_file)) {
        cached <- readRDS(cache_file)
        if (identical(cached$fingerprint, analysis_fingerprint)) return(cached$result)
    }
    continuous <- list(); categories <- list(); overall <- list(); model_spec <- list(); top_contrasts <- list()
    rcs_curves <- list(); rcs_tests <- list(); rcs_hist <- list(); k <- 0L
    d0 <- dat[cohorts[[Y]], , drop=FALSE]
    d0$time <- d0[[paste0(Y, ".t2e")]]; d0$event <- d0[[paste0(Y, ".Yt2e")]]
    for (i in match(diet_subset, use$diet)) {
        diet <- use$diet[i]; label <- use$label[i]
        d <- d0; d$pheno <- d[[use$variable[i]]]; d$group <- d[[paste0("diet.", diet, ".q")]]
        group_medians <- tapply(d$pheno, d$group, median)
        d$trend <- unname(group_medians[as.character(d$group)])
        si <- score_info[score_info$diet == diet, ]
        for (m in MODELS) {
            k <- k+1L; co <- model_covs(diet, m); context <- paste(Y, diet, m)
            base_key <- paste(Y, m, use$alcohol_in_score[i], sep="_")
            ll0 <- get_null(d, co, base_key)
            fl <- fit_cox(d, c("pheno", co), context)
            ll_linear <- as.numeric(logLik(fl))
            ph_linear <- ph_test(fl, "pheno")
            # A linear change of exposure scale has identical fitted risks.
            # Derive both increments from one fit, avoiding duplicate Cox fits.
            continuous[[k]] <- do.call(rbind, lapply(c("P10-P90", "Per 1-SD"), function(scale) {
                step <- if (scale == "P10-P90") si$P10_P90_increment else si$SD
                cbind(data.frame(diet=diet, label=label, outcome=Y, model=m, scale=scale,
                                 N=fl$n, Events=fl$nevent, increment=step, PH_P=ph_linear, AIC_linear=AIC(fl)),
                      extract_effect(fl, "pheno", step))
            }))
            fc <- fit_cox(d, c("group", co), paste(context, "quartiles"))
            ft <- fit_cox(d, c("trend", co), paste(context, "trend"))
            p_global <- lr_p(as.numeric(logLik(fc)), ll0, 3)
            p_trend <- extract_effect(ft, "trend")$P
            ph_cat <- ph_test(fc, "group")
            bdiff <- unname(coef(fc)["groupQ4"] - coef(fc)["groupQ3"])
            vv <- vcov(fc)
            sediff <- sqrt(vv["groupQ4","groupQ4"] + vv["groupQ3","groupQ3"] - 2*vv["groupQ4","groupQ3"])
            top_contrasts[[k]] <- data.frame(diet=diet, outcome=Y, model=m,
                contrast="Q4 versus Q3", HR=exp(bdiff), LCI=exp(bdiff-1.959964*sediff),
                HCI=exp(bdiff+1.959964*sediff), P=2*pnorm(-abs(bdiff/sediff)))
            overall[[k]] <- data.frame(diet=diet, label=label, outcome=Y, model=m, N=fc$n, Events=fc$nevent,
                                      P_global=p_global, P_trend=p_trend, PH_P=ph_cat, AIC_categorical=AIC(fc))
            categories[[k]] <- do.call(rbind, lapply(1:4, function(j) {
                q <- paste0("Q", j); z <- d[d$group == q, , drop=FALSE]
                ef <- if (j == 1) data.frame(HR=1, LCI=1, HCI=1, P=NA_real_) else extract_effect(fc, paste0("group", q))
                cbind(data.frame(diet=diet, label=label, outcome=Y, model=m, category=q, reference="Q1",
                    N=fc$n, Events=fc$nevent, Group_N=nrow(z), Group_Events=sum(z$event),
                    Person_years=sum(z$time), Incidence_percent=100*mean(z$event),
                    Rate_per_100000=1e5*sum(z$event)/sum(z$time),
                    score_min=min(z$pheno), score_max=max(z$pheno), score_median=median(z$pheno)), ef)
            }))
            model_spec[[k]] <- data.frame(diet=diet, outcome=Y, model=m,
                alcohol_in_score=use$alcohol_in_score[i], covariates=paste(co, collapse=";"), region_stratified=FALSE)

            if (diet %in% DOSE_DIETS && m == "M2") {
                # Three-knot restricted cubic spline, equivalent spline space
                # to rms::rcs(x, c(P10,P50,P90)); linear outside outer knots.
                knots <- unname(quantile(d$pheno, c(.1, .5, .9), type=2))
                if (length(unique(knots)) != 3) stop("DOSE RCS requires three distinct knots.")
                nonlinear_basis <- function(x) {
                    (pmax(x-knots[1], 0)^3 - pmax(x-knots[2], 0)^3*(knots[3]-knots[1])/(knots[3]-knots[2]) +
                     pmax(x-knots[3], 0)^3*(knots[2]-knots[1])/(knots[3]-knots[2]))/(knots[3]-knots[1])^2
                }
                d$nlin <- nonlinear_basis(d$pheno)
                fr <- fit_cox(d, c("pheno", "nlin", co), paste(context, "RCS"))
                ref <- median(d$pheno)
                score_levels <- seq(use$minimum[i], use$maximum[i], by=use$step[i])
                xx <- sort(unique(c(seq(min(d$pheno), max(d$pheno), length.out=201),
                                    score_levels[score_levels >= min(d$pheno) & score_levels <= max(d$pheno)], ref)))
                design <- cbind(xx-ref, nonlinear_basis(xx)-nonlinear_basis(ref))
                idx <- c("pheno", "nlin")
                lp <- drop(design %*% coef(fr)[idx])
                se <- sqrt(pmax(0, rowSums((design %*% vcov(fr)[idx, idx]) * design)))
                key <- paste(Y, diet)
                rcs_curves[[key]] <- data.frame(diet=diet, label=label, outcome=Y, score=xx, HR=exp(lp), LCI=exp(lp-1.959964*se), HCI=exp(lp+1.959964*se))
                rcs_tests[[key]] <- data.frame(diet=diet, label=label, outcome=Y, model=m, N=fr$n, Events=fr$nevent,
                    knot1=knots[1], knot2=knots[2], knot3=knots[3], reference=ref,
                    P_overall=lr_p(as.numeric(logLik(fr)), ll0, 2),
                    P_nonlinear=lr_p(as.numeric(logLik(fr)), ll_linear, 1),
                    AIC_linear=AIC(fl), AIC_RCS=AIC(fr), test="Likelihood ratio")
                h <- table(factor(d$pheno, levels=score_levels))
                rcs_hist[[key]] <- data.frame(diet=diet, outcome=Y, score=score_levels, N=as.integer(h))
            }
        }
        cat("Completed ", Y, ": ", label, " (M1/M2/M3)\n", sep="")
        flush.console()
    }
    result <- list(continuous=continuous, categories=categories, overall=overall,
                   model_spec=model_spec, rcs_curves=rcs_curves, rcs_tests=rcs_tests,
                   rcs_hist=rcs_hist, top_contrasts=top_contrasts, warnings=fit_warnings)
    dose_save_rds(list(fingerprint=analysis_fingerprint, result=result), cache_file)
    result
}
analysis_fingerprint <- list(inputs=tools::md5sum(input_files), spec=use, groups=group_info,
    models=MODELS, dose_diets=DOSE_DIETS, functions=lapply(list(run_outcome, model_covs, fit_cox, extract_effect, ph_test, lr_p, get_null),
                                   function(f) deparse(body(f))))
run_outcome_logged <- function(Y, diet_subset=use$diet) {
    withCallingHandlers(run_outcome(Y, diet_subset), warning=function(w) {
        # R 4.6 may emit this while cox.zph rebuilds the model frame.
        # Treatment contrasts and the resulting PH statistic were verified.
        # Suppress only this known message, not convergence/model warnings.
        if (identical(conditionMessage(w), "contrasts dropped from factor group")) invokeRestart("muffleWarning")
    })
}
# Independent per-pattern jobs, each with its own resumable checkpoint.
# Limit parallelism to three workers by default within the step memory scope.
jobs <- expand.grid(outcome=Ys, diet=use$diet, stringsAsFactors=FALSE)
# Prioritize the primary DOSE job; estimates do not depend on order.
jobs <- jobs[order(!jobs$diet %in% DOSE_DIETS), , drop=FALSE]
job_list <- lapply(seq_len(nrow(jobs)), function(i) as.list(jobs[i, ]))
run_job <- function(job) run_outcome_logged(job$outcome, job$diet)
load_pattern_result <- function(job) {
    path <- file.path(DIR_OUT, paste0("checkpoint_", job$outcome, "_", job$diet, ".rds"))
    cached <- if (file.exists(path)) tryCatch(readRDS(path), error=function(e) NULL) else NULL
    required <- c("continuous", "categories", "overall", "model_spec", "rcs_curves",
                  "rcs_tests", "rcs_hist", "top_contrasts", "warnings")
    if (is.list(cached) && identical(cached$fingerprint, analysis_fingerprint) &&
        is.list(cached$result) && all(required %in% names(cached$result))) cached$result else NULL
}
invisible(gc())
association_results <- dose_cached_pattern_jobs(job_list, load_pattern_result, run_job,
    workers=as.integer(Sys.getenv("DOSE_ASSOC_WORKERS", "3")), replace=dose_replace())

collect <- function(name) unlist(lapply(association_results, function(x) x[[name]]), recursive=FALSE, use.names=FALSE)
continuous <- collect("continuous"); categories <- collect("categories"); overall <- collect("overall")
model_spec <- collect("model_spec"); rcs_curves <- collect("rcs_curves")
rcs_tests <- collect("rcs_tests"); rcs_hist <- collect("rcs_hist")
fit_warnings <- unlist(lapply(association_results, `[[`, "warnings"), use.names=FALSE)
bind <- function(x) { z <- do.call(rbind, x); rownames(z) <- NULL; z }
continuous <- bind(continuous); categories <- bind(categories); overall <- bind(overall)
adjust_fdr <- function(d, p, by) {
    out <- rep(NA_real_, nrow(d))
    for (ix in split(seq_len(nrow(d)), interaction(d[by], drop=TRUE))) out[ix] <- p.adjust(d[[p]][ix], method="BH")
    out
}
continuous$FDR <- adjust_fdr(continuous, "P", c("outcome", "model", "scale"))
categories$FDR <- adjust_fdr(categories, "P", c("outcome", "model", "category"))
overall$FDR_global <- adjust_fdr(overall, "P_global", c("outcome", "model"))
overall$FDR_trend <- adjust_fdr(overall, "P_trend", c("outcome", "model"))
write.csv(continuous, file.path(DIR_OUT, "pattern_continuous.csv"), row.names=FALSE)
write.csv(categories, file.path(DIR_OUT, "pattern_categories.csv"), row.names=FALSE)
write.csv(overall, file.path(DIR_OUT, "pattern_overall_tests.csv"), row.names=FALSE)
write.csv(bind(model_spec), file.path(DIR_OUT, "model_specification.csv"), row.names=FALSE)
rcs_curves <- bind(rcs_curves); rcs_tests <- bind(rcs_tests); rcs_hist <- bind(rcs_hist)
top_contrasts <- bind(collect("top_contrasts"))
write.csv(top_contrasts, file.path(DIR_OUT, "pattern_Q4_vs_Q3.csv"), row.names=FALSE)
# Stable primary-score filenames for plotting and later analyses.
for (what in c("curves", "tests", "hist")) {
    z <- switch(what, curves=rcs_curves, tests=rcs_tests, hist=rcs_hist)
    suffix <- if (what == "hist") "distribution" else what
    write.csv(z[z$diet == "dose", ], file.path(DIR_OUT, paste0("DOSE_RCS_", suffix, ".csv")), row.names=FALSE)
}

# Descriptive prediction comparison: SAME DOSE M2 covariates,
# SAME participants, one raw score per model. Apparent (in-sample) Harrell C,
# not validation of the learned food selection/cutoffs or of future prediction.
dp <- dat[cohorts[[PRIMARY_OUTCOME]], , drop=FALSE]
dp$time <- dp[[paste0(PRIMARY_OUTCOME, ".t2e")]]; dp$event <- dp[[paste0(PRIMARY_OUTCOME, ".Yt2e")]]
common_covs <- model_covs("dose", "M2")
if (dose_meta$alcohol) stopifnot(!any(cfg$covs.alcohol %in% common_covs))
fp0 <- fit_cox(dp, common_covs, "Common prediction baseline")
pred <- data.frame(time=dp$time, event=dp$event, baseline=fp0$linear.predictors)
pred_meta <- list()
for (i in seq_len(nrow(use))) {
    dp$pheno <- dp[[use$variable[i]]]
    f <- fit_cox(dp, c("pheno", common_covs), paste("Common prediction", use$diet[i]))
    pred[[use$diet[i]]] <- f$linear.predictors
    pred_meta[[i]] <- data.frame(diet=use$diet[i], label=use$label[i], N=f$n, Events=f$nevent,
        AIC=AIC(f), Delta_AIC=AIC(f)-AIC(fp0), P_added_score=lr_p(as.numeric(logLik(f)), as.numeric(logLik(fp0)), 1))
}
# The joint IJ covariance retains dependence between predictions on the same people.
cc <- concordance(reformulate(c("baseline", use$diet), response="Surv(time, event)"),
                  data=pred, reverse=TRUE, timewt="n")
cval <- as.numeric(cc$concordance); cv <- as.matrix(cc$var)
stopifnot(length(cval) == nrow(use)+1L)
prediction <- bind(pred_meta)
prediction$C_index <- cval[-1]; prediction$C_SE <- sqrt(diag(cv)[-1])
prediction$C_LCI <- pmax(0, prediction$C_index-1.959964*prediction$C_SE)
prediction$C_HCI <- pmin(1, prediction$C_index+1.959964*prediction$C_SE)
prediction$Baseline_C <- cval[1]
prediction$Delta_C <- cval[-1]-cval[1]
prediction$Delta_C_SE <- sqrt(pmax(0, diag(cv)[-1]+cv[1, 1]-2*cv[-1, 1]))
prediction$Delta_C_LCI <- prediction$Delta_C-1.959964*prediction$Delta_C_SE
prediction$Delta_C_HCI <- prediction$Delta_C+1.959964*prediction$Delta_C_SE
id_dose <- match("dose", use$diet)+1L
prediction$C_DOSE_minus_pattern <- cval[id_dose]-cval[-1]
pair_se <- sqrt(pmax(0, cv[id_dose, id_dose]+diag(cv)[-1]-2*cv[id_dose, -1]))
prediction$C_DOSE_minus_pattern_LCI <- prediction$C_DOSE_minus_pattern-1.959964*pair_se
prediction$C_DOSE_minus_pattern_HCI <- prediction$C_DOSE_minus_pattern+1.959964*pair_se
dimnames(cv) <- list(c("baseline", use$diet), c("baseline", use$diet))
write.csv(cv, file.path(DIR_OUT, "pattern_IS_C_index_covariance.csv"), row.names=TRUE)
prediction$evaluation <- "Apparent / in-sample; not independent validation"
prediction$covariates <- paste(common_covs, collapse=";")
prediction$rank_C <- rank(-prediction$C_index, ties.method="min")
write.csv(prediction, file.path(DIR_OUT, "pattern_IS_prediction_common_M2.csv"), row.names=FALSE)

# Association rankings use matched scales; they are descriptive, not tests
# that one dietary pattern is superior to another.
ranking <- continuous[continuous$outcome == PRIMARY_OUTCOME & continuous$model == "M2" & continuous$scale == "P10-P90", ]
sd_res <- continuous[continuous$outcome == PRIMARY_OUTCOME & continuous$model == "M2" & continuous$scale == "Per 1-SD", ]
q4 <- categories[categories$outcome == PRIMARY_OUTCOME & categories$model == "M2" & categories$category == "Q4", ]
ranking$HR_perSD <- sd_res$HR[match(ranking$diet, sd_res$diet)]
ranking$HR_Q4vsQ1 <- q4$HR[match(ranking$diet, q4$diet)]
ranking$P_Q4vsQ1 <- q4$P[match(ranking$diet, q4$diet)]
ranking$rank_HR_P10P90 <- rank(ranking$HR, ties.method="min")
ranking$rank_HR_perSD <- rank(ranking$HR_perSD, ties.method="min")
ranking$C_index_common_M2 <- prediction$C_index[match(ranking$diet, prediction$diet)]
ranking$Delta_C_common_M2 <- prediction$Delta_C[match(ranking$diet, prediction$diet)]
ranking$rank_C_common_M2 <- prediction$rank_C[match(ranking$diet, prediction$diet)]
ranking <- ranking[order(ranking$HR), ]
write.csv(ranking, file.path(DIR_OUT, "pattern_IS_M2_comparison.csv"), row.names=FALSE)

# Primary DOSE summary: matched scales, quartiles, curve tests and fit.
dose_summary <- ranking[match(DOSE_DIETS, ranking$diet), ]
dose_summary$scoring_method <- use$scoring_method[match(DOSE_DIETS, use$diet)]
dose_summary$maximum <- use$maximum[match(DOSE_DIETS, use$diet)]
dose_summary$step <- use$step[match(DOSE_DIETS, use$diet)]
dose_summary$association_covariates <- vapply(DOSE_DIETS,
    function(v) paste(model_covs(v, "M2"), collapse=";"), character(1))
dose_summary$prediction_covariates <- paste(common_covs, collapse=";")
for (q in paste0("Q", 1:4)) {
    z <- subset(categories, outcome == PRIMARY_OUTCOME & model == "M2" & category == q)
    z <- z[match(DOSE_DIETS, z$diet), ]
    for (col in c("HR", "LCI", "HCI", "Group_N", "Group_Events", "score_min", "score_max"))
        dose_summary[[paste0(q, "_", col)]] <- z[[col]]
}
z <- subset(overall, outcome == PRIMARY_OUTCOME & model == "M2")
z <- z[match(DOSE_DIETS, z$diet), ]
dose_summary$P_trend <- z$P_trend
dose_summary$P_global <- z$P_global
dose_summary$PH_P_categorical <- z$PH_P
dose_summary$AIC_categorical <- z$AIC_categorical
z <- rcs_tests[match(DOSE_DIETS, rcs_tests$diet), ]
dose_summary$P_RCS_overall <- z$P_overall
dose_summary$P_RCS_nonlinear <- z$P_nonlinear
dose_summary$AIC_RCS <- z$AIC_RCS
z <- subset(top_contrasts, outcome == PRIMARY_OUTCOME & model == "M2")
z <- z[match(DOSE_DIETS, z$diet), ]
dose_summary$HR_Q4vsQ3 <- z$HR
dose_summary$P_Q4vsQ3 <- z$P
dose_summary$Q4_lowest_point_estimate <- with(dose_summary, Q4_HR <= pmin(Q1_HR,Q2_HR,Q3_HR))
dose_summary$decreasing_quartile_point_estimates <- with(dose_summary, Q1_HR >= Q2_HR & Q2_HR >= Q3_HR & Q3_HR >= Q4_HR)
dose_summary$evaluation <- "Exploratory in-sample comparison; same participants; no threshold search"
write.csv(dose_summary, file.path(DIR_OUT, "DOSE_IS_M2_summary.csv"), row.names=FALSE)

# Check that adding DOSE reproduces existing s1 associations when inputs match.
previous_file <- file.path(DIR_PROJECT, "ukb/01pattern/pattern_continuous.csv")
if (file.exists(previous_file)) {
    previous <- read.csv(previous_file)
    matched <- merge(continuous, previous, by=c("diet", "outcome", "model", "scale"), suffixes=c("_new", "_s1"))
    matched$abs_logHR_difference <- abs(log(matched$HR_new)-log(matched$HR_s1))
    matched$same_N_events <- matched$N == matched$total & matched$Events == matched$case
    write.csv(matched[c("diet", "outcome", "model", "scale", "same_N_events", "abs_logHR_difference")],
              file.path(DIR_OUT, "s1_reproduction_QC.csv"), row.names=FALSE)
    if (any(!matched$same_N_events | matched$abs_logHR_difference > 1e-6)) warning("Some previous s1 results differ; inspect s1_reproduction_QC.csv.")
}
writeLines(unique(fit_warnings), file.path(DIR_OUT, "model_warnings.txt"))
write.csv(data.frame(input=names(input_files), path=unname(input_files), md5=unname(tools::md5sum(input_files))),
          file.path(DIR_OUT, "input_manifest.csv"), row.names=FALSE)
source_files <- file.path(dose_script_dir(),c("score.R","score.R"))
write.csv(data.frame(script=basename(source_files), path=source_files, md5=unname(tools::md5sum(source_files))),
          file.path(DIR_OUT,"code_manifest.csv"),row.names=FALSE)
capture.output(sessionInfo(), file=file.path(DIR_OUT, "sessionInfo.txt"))
writeLines(c(
    paste0("Primary DOSE and ", sum(use$diet != "dose"), " established dietary patterns: incident ischemic stroke."),
    "Oxford WebQ MedDiet (medi24) is included; baseline touchscreen MedDiet-original (medito) is excluded from the comparison and common-cohort eligibility.",
    paste0("Exposure: dose (compatibility alias diet.dose.sum); scoring method: ", dose_meta$method,
           "; theoretical range ", paste(dose_meta$range, collapse="-"), "; step ", dose_meta$step, "."),
    "Primary scoring rules are defined only in score.R and stored in dose_score.rds.",
    "Alternative totals remain in the scoring RDS; only dose enters the current association analyses.",
    paste0("All ", nrow(use), " scores use the same complete-case participants and outcome. Missing established scores are not imputed here."),
    "Quartiles: R quantile type 7; cut right=TRUE, include.lowest=TRUE. Identical scores remain together.",
    "Continuous HRs: P10-P90 increment and one SD from a single linear Cox fit. P-trend uses within-group medians.",
    paste0("BH FDR families: all ", nrow(use), " scores within each model/scale; category contrasts separately by Q2/Q3/Q4."),
    "Primary DOSE RCS: M2, knots P10/P50/P90 type 2, median reference, likelihood-ratio overall/nonlinear tests.",
    paste0("DOSE contains alcohol: ", dose_meta$alcohol, ". Alcohol adjustment follows the score metadata; alcohol-containing scores omit separate alcohol adjustment."),
    paste0("Common C-index covariates: ", paste(common_covs, collapse=";"), ". Paired IJ covariance accounts for correlated predictions."),
    "C-index, AIC and score comparisons are apparent/in-sample; this is not independent validation of food selection or cutoffs.",
    "Do not select a scoring method solely to force the lowest Q4 HR or P-nonlinear>0.05. No cutoff search is performed.",
    "Primary summary: DOSE_IS_M2_summary.csv. Historical DOSE_versions_* files, if present, belong to the earlier version comparison and are not updated or read by these scripts.",
    "Primary curves/tests: DOSE_RCS_*.csv. All model and plot labels use DOSE.",
    "Full results: pattern_IS_M2_comparison.csv, pattern_categories.csv, pattern_overall_tests.csv, pattern_IS_prediction_common_M2.csv.",
    "Checkpoints require matching input hashes, grouping, model specifications and fitting functions."
), file.path(DIR_OUT, "README_analysis.txt"))
cat("\nAnalysis complete. Primary M2 comparison (ordered by P10-P90 HR):\n")
print(ranking[c("label", "HR", "LCI", "HCI", "P", "HR_perSD", "HR_Q4vsQ1", "C_index_common_M2")], row.names=FALSE)
cat("\nPrimary DOSE M2 summary:\n")
print(dose_summary[c("label", "HR_perSD", "HR_Q4vsQ1", "P_trend", "P_RCS_nonlinear", "C_index_common_M2", "Q4_lowest_point_estimate")], row.names=FALSE)
cat("\nDOSE quartiles (Q1 reference):\n")
print(subset(categories, diet %in% DOSE_DIETS & outcome == PRIMARY_OUTCOME & model == "M2",
             select=c(diet, category, score_min, score_max, Group_N, Group_Events, HR, LCI, HCI, P)), row.names=FALSE)
cat("\nDOSE RCS tests:\n")
print(rcs_tests, row.names=FALSE)
cat("\nResults: ", DIR_OUT, "\n", sep="")
}))

# ===========================================================================
# score / sensitivity
# ===========================================================================
dose_register_step("score", "sensitivity", expression({
# DOSE / incident ischemic stroke: component, risk, sensitivity and subgroup analyses.
# Run: .tools/R/bin/Rscript.exe f/score.R
# Requires current score.R and score.R outputs. Base R + survival.
# DOSE_PROJECT may override the project root; DOSE_RISK_BOOT controls optional
# participant bootstrap CIs for absolute risk (default 0: point estimates only).
# No thresholds are refitted to the outcome. All analyses are exploratory.
# Add/update dietary sensitivities without rerunning completed risk/subgroup models:
# Add/update PRS strata, continuous interaction and joint models only:

suppressPackageStartupMessages(library(survival))
options(contrasts=c('contr.treatment','contr.poly'))
ROOT <- dose_output_dir()
OUT <- file.path(ROOT, 'ukb/04score/IS_sensitivity')
DIETARY_ONLY <- '--dietary-only' %in% commandArgs(trailingOnly=TRUE)
PRS_ONLY <- '--prs-only' %in% commandArgs(trailingOnly=TRUE)
COMPONENTS_ONLY <- '--components-only' %in% commandArgs(trailingOnly=TRUE)
DOSE9_COMPONENTS_ONLY <- '--dose9-components-only' %in% commandArgs(trailingOnly=TRUE)
stopifnot(sum(c(DIETARY_ONLY,PRS_ONLY,COMPONENTS_ONLY,DOSE9_COMPONENTS_ONLY))<=1)
dir.create(file.path(OUT,'figures'), recursive=TRUE, showWarnings=FALSE)
csv <- function(x,n) write.csv(x,file.path(OUT,paste0(n,'.csv')),row.names=FALSE,na='')
logmsg <- function(...) cat(format(Sys.time(),'%H:%M:%S'), ..., '\n')
issues <- data.frame(context=character(), message=character())
note <- function(context,message) {
  issues <<- rbind(issues,data.frame(context,message)); logmsg(context,message)
  csv(issues,'analysis_warnings')
}
files <- c(cohort='ukb/00data/dat.base.imp.rds',
           raw='ukb/00data/dat.base.rds',
           config='ukb/00data/analysis_config.rds',
           score='ukb/04score/dose_score.rds',
           spec='ukb/04score/pattern/pattern_spec.csv')
files <- file.path(ROOT,files) |> setNames(names(files))
stopifnot(all(file.exists(files)))
if(DIETARY_ONLY || PRS_ONLY || COMPONENTS_ONLY || DOSE9_COMPONENTS_ONLY) {
  previous <- read.csv(file.path(OUT,'input_manifest.csv'))
  stopifnot(identical(previous$input,names(files)),
            identical(previous$md5,unname(tools::md5sum(files))))
}
csv(data.frame(input=names(files),path=files,md5=unname(tools::md5sum(files))), 'input_manifest')
cfg <- readRDS(files['config']); score <- readRDS(files['score'])
comps <- attr(score,'component_columns',exact=TRUE)
meta <- attributes(score)
stopifnot(identical(meta$score_version,'dose'),
          identical(meta$score_intake_unit,'servings/day'),length(comps)>1,
          all(c('eid','dose',comps) %in% names(score)),
          isTRUE(all.equal(score$dose,rowSums(score[comps]),check.attributes=FALSE)))
co <- unique(c(cfg$covs.m1,if (!meta$score_contains_alcohol) cfg$covs.alcohol,cfg$covs.m2.add))
co3 <- unique(c(co,cfg$covs.m3.add))
dat <- as.data.frame(readRDS(files['cohort']))
stopifnot(!anyNA(dat$eid),!anyDuplicated(dat$eid),!anyDuplicated(score$eid),
          setequal(as.character(dat$eid),as.character(score$eid)))
dat$eid <- as.character(dat$eid)
new <- setdiff(names(score),'eid'); dat[new] <- score[match(dat$eid,as.character(score$eid)),new]
dat$time <- dat$cvd_stroke_i.t2e; dat$event <- dat$cvd_stroke_i.Yt2e
stopifnot(all(na.omit(dat$event) %in% 0:1),!any(dat$prior_qi_stroke,na.rm=TRUE))
spec <- read.csv(files['spec'])
commonvars <- unique(c('time','event',cfg$covs.all,spec$variable[spec$include & spec$available]))
ok <- complete.cases(dat[commonvars]) & is.finite(dat$time) & dat$time>0
d <- dat[ok,]; rownames(d) <- NULL
stopifnot(nrow(d)>1000,all(complete.cases(d[c('dose',comps,co3)])))
# Match the primary comparison's cohort and grouping; ties are never split.
br <- unname(quantile(dat$dose,c(0,.25,.5,.75,1),na.rm=TRUE))
stopifnot(!anyDuplicated(br))
qgroup <- function(x,breaks=br) cut(x,breaks,include.lowest=TRUE,right=TRUE,labels=paste0('Q',1:4))
dat$group <- qgroup(dat$dose); d$group <- qgroup(d$dose)
full_sd <- sd(dat$dose,na.rm=TRUE)
full_p <- diff(quantile(dat$dose,c(.1,.9),na.rm=TRUE)) |> unname()
csv(data.frame(group=paste0('Q',1:4),lower=br[1:4],upper=br[2:5],
               lower_inclusive=c(TRUE,FALSE,FALSE,FALSE),upper_inclusive=TRUE), 'quartile_definitions')
csv(data.frame(input_N=nrow(dat),analysis_N=nrow(d),events=sum(d$event),
               person_years=sum(d$time),SD=full_sd,P10_P90=full_p,
               scoring_method=meta$primary_scoring_method), 'cohort_QC')
# Avoid retaining unrelated participant variables in every model/subset copy.
extras <- c('start_date','date_death','date_lost','end0','cvd','enroll','region',
            'prs_istroke',paste0('PC',1:10),'diet.dash.sum','diet.ahei.sum')
keep <- unique(c('eid','time','event','group',new,cfg$covs.all,extras))
keep <- intersect(keep,names(dat)); dat <- dat[keep]; d <- d[keep]
invisible(gc())
logmsg('Primary cohort',nrow(d),'IS events',sum(d$event))

prepare <- function(z,terms) {
  v <- unique(c('time','event',all.vars(reformulate(terms))))
  z <- z[complete.cases(z[v]) & is.finite(z$time) & z$time>0,,drop=FALSE]
  droplevels(z)
}
varying <- function(z,vars) vars[vapply(vars,function(v) length(unique(z[[v]]))>1,logical(1))]
fit <- function(z,terms,context) {
  withCallingHandlers({
    f <- coxph(reformulate(terms,response='Surv(time,event)'),data=z,
               ties='efron',singular.ok=FALSE,x=FALSE,y=TRUE,model=FALSE)
    if(any(!is.finite(coef(f)))) stop('Nonfinite coefficients: ',context)
    f
  },warning=function(w) {note(context,conditionMessage(w));invokeRestart('muffleWarning')})
}
effect <- function(f,term,increment=1) {
  b <- unname(coef(f)[term]); s <- sqrt(vcov(f)[term,term])
  data.frame(HR=exp(increment*b),LCI=exp(increment*(b-1.959964*s)),
             HCI=exp(increment*(b+1.959964*s)),P=2*pnorm(-abs(b/s)))
}
ph <- function(f,term) tryCatch(unname(cox.zph(f)$table[term,'p']),
  error=function(e) {note('PH test',conditionMessage(e));NA_real_})
model_specs <- list(); continuous <- list(); categorical <- list(); eligibility <- list()
run <- function(z,label,covs=co,x='dose',categories=TRUE,fixed_scale=TRUE,extra_terms=character()) {
  logmsg('Cox:',label)
  z$exposure <- z[[x]]
  terms <- c('exposure',covs,extra_terms)
  z <- prepare(z,terms); covs <- varying(z,covs)
  eligibility[[label]] <<- data.frame(analysis=label,N=nrow(z),events=sum(z$event),
                                      person_years=sum(z$time),excluded_from_primary=nrow(d)-nrow(z))
  if(nrow(z)<200 || sum(z$event)<20 || sd(z$exposure)==0) {
    note(label,'Not estimable: fewer than 200 participants / 20 events, or invariant score'); return(NULL)
  }
  f <- fit(z,c('exposure',covs,extra_terms),label)
  incr <- c('Per point'=1,'Per 1-SD'=if(fixed_scale) full_sd else sd(z$exposure),
            'P10-P90'=if(fixed_scale) full_p else unname(diff(quantile(z$exposure,c(.1,.9)))))
  pp <- ph(f,'exposure')
  continuous[[label]] <<- do.call(rbind,lapply(names(incr),function(s)
    cbind(data.frame(analysis=label,scale=s,increment=unname(incr[s]),N=f$n,events=f$nevent,PH_P=pp),
          effect(f,'exposure',incr[s]))))
  model_specs[[label]] <<- data.frame(analysis=label,exposure=x,
    covariates=paste(c(covs,extra_terms),collapse=' + '),
    scale_population=if(fixed_scale)'Full source cohort DOSE' else 'Scenario cohort and score')
  if(categories) {
    if(!fixed_scale || x!='dose') {
      qb <- unname(quantile(z$exposure,c(0,.25,.5,.75,1)))
      if(anyDuplicated(qb)) {note(label,'Quartiles unavailable: duplicated boundaries');return(f)}
      z$group <- qgroup(z$exposure,qb)
    }
    cnt <- table(z$group); ev <- tapply(z$event,z$group,sum)
    if(any(cnt==0) || any(ev<5)) {note(label,'Quartiles unavailable: fewer than 5 events in a group');return(f)}
    med <- tapply(z$exposure,z$group,median); z$trend <- unname(med[as.character(z$group)])
    fc <- fit(z,c('group',covs,extra_terms),paste(label,'quartiles'))
    ft <- fit(z,c('trend',covs,extra_terms),paste(label,'trend'))
    pt <- effect(ft,'trend')$P
    categorical[[label]] <<- do.call(rbind,lapply(1:4,function(k) {
      g <- paste0('Q',k)
      e <- if(k==1) data.frame(HR=1,LCI=NA_real_,HCI=NA_real_,P=NA_real_) else effect(fc,paste0('group',g))
      cbind(data.frame(analysis=label,group=g,N=as.integer(cnt[g]),events=as.integer(ev[g]),
                       observed_min=min(z$exposure[z$group==g]),observed_max=max(z$exposure[z$group==g]),P_trend=pt),e)
    }))
  }
  invisible(f)
}
save_models <- function() {
  csv(do.call(rbind,continuous),'sensitivity_continuous')
  csv(do.call(rbind,categorical),'sensitivity_quartiles')
  csv(do.call(rbind,model_specs),'model_specifications')
  csv(do.call(rbind,eligibility),'scenario_cohort_QC')
}

run_dietary_adjustments <- function(d) {
  # 7. Nutrients: correct the legacy code's mislabeled trans-fat field 26155.
  # UKB: https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=100117
  nf <- Sys.getenv('DOSE_NUTRIENT_FILE',file.path(dose_ukb_dir(), "Rdata/diet.fudan.rds"))
  gf <- file.path(ROOT,'ukb/02foodwas/dat.group.rds')
  mf <- file.path(ROOT,'ukb/02foodwas/food_group_map.csv')
  stopifnot(all(file.exists(c(nf,gf,mf))))
  inp <- c(nutrients=nf,food_groups=gf,food_map=mf)
  csv(data.frame(input=names(inp),path=unname(inp),md5=unname(tools::md5sum(inp))),
      'dietary_input_manifest')
  nutrient_spec <- data.frame(nutrient=c('sodium','potassium','fiber','saturated_fat'),
    field=c('w.p26052','w.p26024','w.p26017','w.p26014'),unit=c('mg/day','mg/day','g/day','g/day'),
    source_url=paste0('https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=',c(26052,26024,26017,26014)))
  csv(nutrient_spec,'nutrient_field_definitions')
  nr <- as.data.frame(readRDS(nf)); nr$eid <- as.character(nr$eid)
  stopifnot(!anyNA(nr$eid),!anyDuplicated(nr$eid),all(nutrient_spec$field %in% names(nr)))
  z <- d; idx <- match(z$eid,nr$eid); nq <- list(); nv <- character()
  for(k in seq_len(nrow(nutrient_spec))) {
    x <- nr[[nutrient_spec$field[k]]][idx]
    stopifnot(is.numeric(x),all(is.na(x) | (is.finite(x) & x>=0)))
    sx <- sd(x,na.rm=TRUE); mx <- mean(x,na.rm=TRUE); stopifnot(is.finite(sx),sx>0)
    v <- paste0('nutrient.',nutrient_spec$nutrient[k]); nv <- c(nv,v)
    z[[v]] <- (x-mx)/sx
    nq[[k]] <- data.frame(nutrient=nutrient_spec$nutrient[k],field=nutrient_spec$field[k],
      N=nrow(z),unmatched_ID=sum(is.na(idx)),missing=sum(is.na(x)),mean=mx,SD=sx,
      model_variable=v,scaling='One SD (nuisance covariate only)')
  }
  rm(nr); invisible(gc()); csv(do.call(rbind,nq),'nutrient_covariate_QC')
  matched_baselines <- list()
  add_model <- function(zz,label,cv) {
    zz <- prepare(zz,c('dose',cv))
    # When extra-covariate missingness changes N, separate selection from adjustment.
    if(nrow(zz)!=nrow(d) && !any(vapply(matched_baselines,function(ids)identical(ids,zz$eid),logical(1)))) {
      matched_baselines[[label]] <<- zz$eid
      run(zz,paste0(label,'_matched_M2'),covs=co)
    }
    run(zz,label,covs=cv); save_models()
  }
  for(k in seq_along(nv)) add_model(z,paste0('Additional_',nutrient_spec$nutrient[k]),c(co,nv[k]))
  add_model(z,'Additional_four_nutrients',c(co,nv))

  # 8. All OTHER mapped food groups, not only outcome-selected candidates.
  # Exclude the 11 underlying foods of the current 9-domain score, not DOSE7.
  gr <- as.data.frame(readRDS(gf)); gr$eid <- as.character(gr$eid)
  map <- read.csv(mf,stringsAsFactors=FALSE)
  stopifnot(!anyNA(gr$eid),!anyDuplicated(gr$eid),all(c('Food.Group','name') %in% names(map)))
  source_rules <- meta$food_scoring_rules
  stopifnot(is.data.frame(source_rules),'Food' %in% names(source_rules))
  score_foods <- unique(source_rules$Food)
  all_foods <- unique(map$Food.Group)
  stopifnot(all(score_foods %in% all_foods),all(all_foods %in% names(gr)))
  score_items <- unique(map$name[map$Food.Group %in% score_foods])
  idx <- match(z$eid,gr$eid); fq <- list(); fv <- character()
  for(g in all_foods) {
    x <- gr[[g]][idx]; stopifnot(is.numeric(x),all(is.na(x) | (is.finite(x) & x>=0)))
    np <- sum(!is.na(x)); nu <- length(unique(x[!is.na(x)]))
    qu <- if(np) unname(quantile(x,c(.1,.9),na.rm=TRUE)) else c(NA_real_,NA_real_)
    sx <- sd(x,na.rm=TRUE); inc <- diff(qu)
    overlap <- any(map$name[map$Food.Group==g] %in% score_items)
    reason <- if(g %in% score_foods) 'Current DOSE constituent food' else if(overlap)
      'Shares food items with DOSE' else if(nu<2 || !is.finite(sx) || sx==0) 'Invariant/unavailable' else 'Included'
    fallback <- !is.finite(inc) || inc==0
    v <- paste0('otherfood.',make.names(g))
    fq[[g]] <- data.frame(food_group=g,model_variable=v,N=nrow(z),unmatched_ID=sum(is.na(idx)),
      missing=sum(is.na(x)),unique_values=nu,P10=qu[1],P90=qu[2],SD=sx,
      covariate_increment=if(fallback) sx else inc,
      scaling=if(fallback) 'SD fallback: identical P10/P90' else 'P10-P90',decision=reason)
    if(reason=='Included') {
      stopifnot(!v %in% fv)
      z[[v]] <- (x-if(fallback) mean(x,na.rm=TRUE) else qu[1]) / if(fallback) sx else inc
      fv <- c(fv,v)
    }
  }
  csv(do.call(rbind,fq),'other_foodgroup_covariate_QC')
  stopifnot(length(fv)>0)
  rm(gr); invisible(gc())
  add_model(z,'Additional_other_foodgroups',c(co,fv))
  # More stringent exploratory check, reported separately from the legacy models.
  add_model(z,'Additional_nutrients_and_other_foodgroups',c(co,nv,fv))
  csv(issues,'analysis_warnings')
  capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo_dietary_adjustments.txt'))
  writeLines(c('Dietary sensitivity update',
    'M2 + sodium; potassium; fiber; saturated fat, individually and jointly.',
    'Saturated fat uses w.p26014. Legacy w.p26155 is TRANS fat and is not substituted for saturated fat.',
    paste('Current constituent foods excluded:',paste(score_foods,collapse=', ')),
    paste('Other food groups included:',length(fv)),
    'All other mapped, variable food groups are adjusted jointly; none are selected using IS P values.',
    'Nutrient/food missingness is handled by complete cases; a same-cohort M2 comparator is added if N changes.',
    'DOSE increment and quartile cutoffs remain fixed; only adjustment covariates are scaled.',
    'LOO retains the CURRENT nine DOSE9half components, with own-score increments and quartile boundaries.',
    'These are overadjustment checks, not mediation estimates or proof of causal independence.',
    'Table_S13_Sensitivity.csv combines continuous and quartile estimates including P values and group definitions.'),
    file.path(OUT,'dietary_sensitivity_methods.txt'))
}

run_prs_analysis <- function(d) {
  # Genetic strata and continuous interaction: current DOSE9half, M2 + PC1-PC10.
  # Keep this module independent of aggregate sensitivity tables and risk models.
  logmsg('PRS: genetic strata, continuous interaction and joint analysis')
  pc <- paste0('PC',1:10)
  stopifnot(all(c('prs_istroke',pc) %in% names(d)),is.numeric(d$prs_istroke))
  cv <- unique(c(co,pc))
  pb <- unname(quantile(d$prs_istroke,c(0,1/3,2/3,1),na.rm=TRUE))
  stopifnot(all(is.finite(pb)),!anyDuplicated(pb))
  z <- prepare(d,c('dose','prs_istroke',cv))
  stopifnot(all(is.finite(z$prs_istroke)),all(vapply(z[pc],function(x)all(is.finite(x)),logical(1))))
  z$prs_group <- cut(z$prs_istroke,pb,include.lowest=TRUE,labels=c('Low','Middle','High'))
  pm <- mean(z$prs_istroke); ps <- sd(z$prs_istroke)
  stopifnot(is.finite(ps),ps>0)
  z$prs_z <- (z$prs_istroke-pm)/ps
  z$dose_z <- (z$dose-mean(z$dose))/full_sd
  cv <- varying(z,cv)
  csv(data.frame(primary_N=nrow(d),N=nrow(z),events=sum(z$event),
    excluded=nrow(d)-nrow(z),missing_PRS=sum(is.na(d$prs_istroke)),
    missing_any_PC=sum(!complete.cases(d[pc])),PRS_mean=pm,PRS_SD=ps,
    DOSE_SD=full_sd,DOSE_P10_P90=full_p),'prs_cohort_QC')
  csv(data.frame(level=c('Low','Middle','High'),lower=pb[1:3],upper=pb[2:4],
    lower_inclusive=c(TRUE,FALSE,FALSE),upper_inclusive=TRUE,
    definition_population='Primary cohort with observed PRS'),'prs_stratum_definitions')
  specs <- list()
  pf <- function(zz,terms,label) {
    specs[[label]] <<- data.frame(model=label,formula=paste('Surv(time,event) ~',paste(terms,collapse=' + ')))
    logmsg('PRS Cox:',label)
    fit(zz,terms,label)
  }
  ref <- data.frame(HR=1,LCI=NA_real_,HCI=NA_real_,P=NA_real_)
  cont <- cats <- list()
  for(g in levels(z$prs_group)) {
    zz <- droplevels(z[z$prs_group==g,]); vc <- varying(zz,cv)
    stopifnot(nrow(zz)>=200,sum(zz$event)>=20,all(table(zz$group)>0),
              all(tapply(zz$event,zz$group,sum)>=5))
    ff <- pf(zz,c('dose',vc),paste0(g,' continuous'))
    inc <- c('Per point'=1,'Per 1-SD'=full_sd,'P10-P90'=full_p)
    cont[[g]] <- do.call(rbind,lapply(names(inc),function(s)
      cbind(data.frame(PRS=g,scale=s,increment=unname(inc[s]),N=ff$n,events=ff$nevent),effect(ff,'dose',inc[s]))))
    fc <- pf(zz,c('group',vc),paste0(g,' quartiles'))
    med <- tapply(zz$dose,zz$group,median); zz$trend <- unname(med[as.character(zz$group)])
    ft <- pf(zz,c('trend',vc),paste0(g,' trend'))
    cats[[g]] <- do.call(rbind,lapply(levels(z$group),function(q) {
      ix <- zz$group==q
      cbind(data.frame(PRS=g,DOSE=q,N=sum(ix),events=sum(zz$event[ix]),
        person_years=sum(zz$time[ix]),P_trend=effect(ft,'trend')$P),
        if(q=='Q1') ref else effect(fc,paste0('group',q)))
    }))
  }
  cc <- do.call(rbind,cont); qq <- do.call(rbind,cats)
  csv(cc,'prs_stratified_continuous'); csv(qq,'prs_stratified_quartiles')
  # Identical observations, main effects and covariates in both nested models.
  f0 <- pf(z,c('dose_z','prs_z',cv),'Continuous interaction null')
  f1 <- pf(z,c('dose_z*prs_z',cv),'Continuous interaction alternative')
  lrt <- function(a,b) {
    stopifnot(a$n==b$n,a$nevent==b$nevent,identical(a$y,b$y))
    df <- length(coef(b))-length(coef(a)); stat <- 2*as.numeric(logLik(b)-logLik(a))
    stopifnot(df>0,stat > -1e-6)
    data.frame(N=b$n,events=b$nevent,df=df,LRT=max(0,stat),P_interaction=pchisq(max(0,stat),df,lower.tail=FALSE))
  }
  it <- cbind(data.frame(interaction='DOSE per SD x PRS per SD'),lrt(f0,f1),effect(f1,'dose_z:prs_z'))
  names(it)[names(it)=='P'] <- 'P_Wald'
  csv(it,'prs_continuous_interaction')
  gi <- list()
  for(ex in c('dose_z','group')) {
    a <- pf(z,c(ex,'prs_group',cv),paste(ex,'strata interaction null'))
    b <- pf(z,c(paste0(ex,'*prs_group'),cv),paste(ex,'strata interaction alternative'))
    gi[[ex]] <- cbind(data.frame(exposure=ex),lrt(a,b))
  }
  csv(do.call(rbind,gi),'prs_categorical_interactions')
  # Conditional DOSE slope at selected continuous PRS values (not strata means).
  slopes <- do.call(rbind,lapply(c(-1,0,1),function(v) {
    w <- setNames(rep(0,length(coef(f1))),names(coef(f1)))
    w['dose_z'] <- 1; w['dose_z:prs_z'] <- v
    b <- sum(w*coef(f1)); se <- sqrt(drop(t(w)%*%vcov(f1)%*%w))
    data.frame(PRS_z=v,HR=exp(b),LCI=exp(b-1.959964*se),HCI=exp(b+1.959964*se),P=2*pnorm(-abs(b/se)))
  }))
  csv(slopes,'prs_conditional_DOSE_perSD')
  z$joint <- interaction(z$prs_group,z$group,sep=' | ',drop=TRUE)
  z$joint <- relevel(z$joint,ref='Low | Q1')
  fj <- pf(z,c('joint',cv),'Joint PRS and DOSE')
  joint <- do.call(rbind,lapply(levels(z$joint),function(g) {
    ix <- z$joint==g
    cbind(data.frame(joint=g,N=sum(ix),events=sum(z$event[ix])),
      if(g=='Low | Q1') ref else effect(fj,paste0('joint',g)))
  }))
  csv(joint,'prs_DOSE_joint'); csv(do.call(rbind,specs),'prs_model_specifications')
  hf <- function(x) ifelse(is.na(x$LCI),'Ref.',sprintf('%.3f (%.3f, %.3f)',x$HR,x$LCI,x$HCI))
  jj <- joint[match(paste(qq$PRS,qq$DOSE,sep=' | '),joint$joint),]
  tab <- data.frame(PRS=qq$PRS,DOSE=qq$DOSE,Events_Total=paste0(qq$events,'/',qq$N),
    Within_PRS_HR=hf(qq),Within_PRS_P=qq$P,P_trend=qq$P_trend,
    Joint_HR=hf(jj),Joint_P=jj$P,check.names=FALSE)
  csv(tab,'Table_S15_PRS_DOSE_joint')
  cs <- cc[cc$scale=='P10-P90',]; qs <- qq[qq$DOSE=='Q4',]
  csv(data.frame(PRS=cs$PRS,N=cs$N,events=cs$events,DOSE_P10_P90_HR=hf(cs),P=cs$P,
    Q4_vs_Q1_HR=hf(qs),P_Q4=qs$P,P_trend=qs$P_trend,
    P_continuous_interaction=it$P_interaction),'Table_PRS_summary')
  writeLines(c('# PRS analysis',
    sprintf('Complete-case genetic cohort: %s participants; %s IS events; %s excluded.',nrow(z),sum(z$event),nrow(d)-nrow(z)),
    'Score: current DOSE9half. All models: current M2 + PC1-PC10, Efron ties.',
    'PRS tertiles use observed PRS in the primary cohort; no PRS or PC imputation.',
    'DOSE SD and quartile boundaries remain those of the original source cohort.',
    'Continuous PRS is centered/scaled within the complete-case genetic cohort.',
    'Uses primary IS follow-up to remain comparable with existing subgroup models.',
    'The legacy example stratified baseline hazards by region; this module retains the current M2 specification.',
    'Continuous interaction: nested-model 1-df LRT; HR is a ratio of DOSE per-SD HRs per 1-SD increase in PRS.',
    sprintf('Continuous interaction ratio of HRs: %s; LRT P=%.6g.',hf(it),it$P_interaction),
    'Categorical interaction LRTs (2 and 6 df) are secondary; reported P values are unadjusted.',
    'Within-stratum reference: Q1 of that PRS stratum; joint reference: Low PRS + Q1.',
    'Observational associations do not establish efficacy; nonsignificant interaction does not establish equivalence.'),
    file.path(OUT,'PRS_results.md'))
  capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo_PRS.txt'))
  stopifnot(nrow(cc)==9,nrow(qq)==12,nrow(joint)==12,sum(qq$N)==nrow(z),
    sum(qq$events)==sum(z$event),it$df==1,all(is.finite(cc$HR)),all(is.finite(joint$HR)))
  logmsg('PRS analyses complete; all internal checks passed')
}
run_component_categories <- function(d) {
  logmsg('Component categories: unordered score levels with 0 as reference')
  res <- overall <- counts <- specs <- km <- risk <- list()
  for(v in comps) {
    z <- d; lev <- sort(unique(as.numeric(z[[v]])))
    stopifnot(identical(lev,c(0,1)) || identical(lev,c(0,.5,1)))
    z$component_group <- factor(z[[v]],levels=lev,ordered=FALSE)
    counts[[v]] <- do.call(rbind,lapply(levels(z$component_group),function(g) {
      ix <- z$component_group==g
      data.frame(component=v,score=g,N=sum(ix),events=sum(z$event[ix]),person_years=sum(z$time[ix]))
    }))
    stopifnot(all(counts[[v]]$events>=5))
    for(mode in c('M2','M2_mutually_adjusted')) {
      cv <- c(co,if(mode=='M2_mutually_adjusted')setdiff(comps,v))
      logmsg('Component Cox:',v,mode)
      f <- fit(z,c('component_group',cv),paste(v,mode,'categorical'))
      f0 <- fit(z,cv,paste(v,mode,'null'))
      df <- length(coef(f))-length(coef(f0)); lr <- 2*as.numeric(logLik(f)-logLik(f0))
      stopifnot(f$n==f0$n,identical(f$y,f0$y),df==length(lev)-1,lr > -1e-6)
      key <- paste(v,mode)
      overall[[key]] <- data.frame(component=v,model=mode,N=f$n,events=f$nevent,
        df=df,LRT=max(0,lr),P_overall=pchisq(max(0,lr),df,lower.tail=FALSE),PH_P=ph(f,'component_group'))
      specs[[key]] <- data.frame(component=v,model=mode,
        formula=paste('Surv(time,event) ~',paste(c('component_group',cv),collapse=' + ')))
      res[[key]] <- do.call(rbind,lapply(levels(z$component_group)[-1],function(g)
        cbind(data.frame(component=v,model=mode,comparison=paste(g,'vs 0'),score=as.numeric(g),
          N=f$n,events=f$nevent),effect(f,paste0('component_group',g)))))
    }
    sf <- survfit(Surv(time,event)~component_group,data=z,conf.type='log-log')
    # Exact event-time steps through 15 years, with confidence limits; no censor marks.
    tt <- sort(unique(c(0,z$time[z$event==1 & z$time<=15],15)))
    sm <- summary(sf,times=tt,extend=FALSE)
    km[[v]] <- data.frame(component=v,score=sub('component_group=','',as.character(sm$strata)),
      time=sm$time,survival=sm$surv,lower=sm$lower,upper=sm$upper)
    sr <- summary(sf,times=c(0,5,10,15),extend=TRUE)
    risk[[v]] <- data.frame(component=v,score=sub('component_group=','',as.character(sr$strata)),
      time=sr$time,n_risk=sr$n.risk)
  }
  rr <- do.call(rbind,res); oo <- do.call(rbind,overall); nn <- do.call(rbind,counts)
  rr$P_FDR <- ave(rr$P,rr$model,FUN=function(p)p.adjust(p,'BH'))
  oo$P_overall_FDR <- ave(oo$P_overall,oo$model,FUN=function(p)p.adjust(p,'BH'))
  csv(rr,'component_categorical_IS_associations'); csv(oo,'component_categorical_overall_tests')
  csv(nn,'component_category_counts'); csv(do.call(rbind,specs),'component_categorical_model_specs')
  csv(do.call(rbind,km),'component_KM_estimates'); csv(do.call(rbind,risk),'component_KM_risk_table')
  tab <- rr
  idx <- match(paste(rr$component,rr$score),paste(nn$component,nn$score))
  zero <- nn[nn$score=='0',]; iz <- match(rr$component,zero$component)
  tab$Events_Total <- paste0(nn$events[idx],'/',nn$N[idx])
  tab$Reference_Events_Total <- paste0(zero$events[iz],'/',zero$N[iz])
  tab$HR_95CI <- sprintf('%.3f (%.3f, %.3f)',rr$HR,rr$LCI,rr$HCI)
  j <- match(paste(rr$component,rr$model),paste(oo$component,oo$model))
  tab$P_overall <- oo$P_overall[j]; tab$P_overall_FDR <- oo$P_overall_FDR[j]
  csv(tab,'Table_component_categorical')
  label <- function(v) gsub('_',' ',sub('[.]dose$','',v))
  cols <- c('0'='#777777','0.5'='#D89220','1'='#176B91')
  export <- function(name,fun,w,h) {
    pdf(file.path(OUT,'figures',paste0(name,'.pdf')),width=w,height=h); fun(); bounded_dev_off()
    bounded_png(file.path(OUT,'figures',paste0(name,'.png')),width=w,height=h,units='in',res=180); fun(); bounded_dev_off()
  }
  ymin <- floor(min(vapply(km,function(a)min(a$lower,na.rm=TRUE),numeric(1)))*100)/100
  kmplot <- function() {
    par(mfrow=c(3,3),mar=c(8,4,2.6,1),oma=c(1,0,3,0),mgp=c(2.5,.7,0))
    for(v in comps) {
      a <- km[[v]]; ng <- counts[[v]]; rt <- risk[[v]]
      plot(NA,xlim=c(0,15),ylim=c(ymin,1),xaxs='i',yaxs='i',xaxt='n',
        xlab='Follow-up (years)',ylab='IS-free survival probability',main=label(v))
      axis(1,at=c(0,5,10,15)); abline(h=axTicks(2),col='#EEEEEE')
      for(g in ng$score) {
        b <- a[a$score==g,]
        step <- function(x,y) list(x=rep(x,each=2)[-1],y=rep(y,each=2)[-length(rep(y,each=2))])
        lo <- step(b$time,b$lower); hi <- step(b$time,b$upper)
        polygon(c(lo$x,rev(hi$x)),c(lo$y,rev(hi$y)),border=NA,col=adjustcolor(cols[g],alpha.f=.12))
        lines(b$time,b$survival,type='s',col=cols[g],lwd=1.6)
      }
      legend('bottomleft',legend=paste0(ng$score,' point: ',ng$events,'/',ng$N),
        col=cols[ng$score],lty=1,lwd=1.6,bty='n',cex=.73,title='Score: events / participants')
      mtext('Number at risk',side=1,line=3.5,adj=0,cex=.68)
      for(k in seq_len(nrow(ng))) {
        g <- ng$score[k]; b <- rt[rt$score==g,]
        mtext(paste0(g,':'),side=1,line=3.5+k,adj=-.09,col=cols[g],cex=.67)
        mtext(format(b$n_risk,trim=TRUE),side=1,line=3.5+k,at=b$time,cex=.67,col=cols[g])
      }
    }
    mtext('a  DOSE components and IS-free survival',outer=TRUE,line=1,font=2,cex=1.2)
    mtext('Unadjusted Kaplan-Meier curves with 95% CI; common truncated y-axis. Death is censored.',outer=TRUE,side=1,line=0,cex=.8)
  }
  forest <- function(mode) {
    a <- tab[tab$model==mode,]; a <- a[order(match(a$component,comps),a$score),]
    n <- nrow(a); yy <- rev(seq_len(n))
    layout(matrix(1:3,nrow=1),widths=c(3.2,2.5,3.7))
    par(mar=c(4,0,4,0)); plot.new(); plot.window(xlim=c(0,1),ylim=c(.3,n+1.4))
    text(0,yy,paste(label(a$component),paste0('  (',a$comparison,')')),adj=0,cex=.86)
    text(0,n+1,'Component / comparison',adj=0,font=2,cex=.9)
    par(mar=c(4,.5,4,.5)); plot(a$HR,yy,log='x',xlim=range(c(.6,1.15,a$LCI,a$HCI)),
      ylim=c(.3,n+1.4),yaxt='n',ylab='',xlab='Adjusted HR (95% CI)',bty='n',pch=19,col='#176B91')
    segments(a$LCI,yy,a$HCI,yy,col='#176B91',lwd=1.5); abline(v=1,lty=2,col='grey50')
    par(mar=c(4,0,4,.5)); plot.new(); plot.window(xlim=c(0,1),ylim=c(.3,n+1.4))
    text(.01,yy,a$HR_95CI,adj=0,cex=.84)
    fmt <- function(p)ifelse(p<.001,'<0.001',sprintf('%.3f',p))
    text(.70,yy,fmt(a$P),cex=.84); text(.92,yy,fmt(a$P_overall),cex=.84)
    text(.01,n+1,'HR (95% CI)',adj=0,font=2,cex=.9)
    text(.70,n+1,'P',font=2,cex=.9); text(.92,n+1,'Overall P',font=2,cex=.85)
    mtext(paste('b  DOSE component associations -',mode),side=3,line=1,adj=1,font=2)
  }
  export('component_categorical_KM',kmplot,15,13)
  export('component_categorical_forest_M2',function()forest('M2'),12,7)
  export('component_categorical_forest_mutually_adjusted',function()forest('M2_mutually_adjusted'),12,7)
  writeLines(c('# Component categorical analysis',
    'Fruit and vegetables: unordered 0 / 0.5 / 1 categories; other components: 0 / 1. Reference = 0.',
    'Primary forest: separate M2-adjusted Cox model for each component. Supplement: M2 plus the other original component scores.',
    'No model adjusts for total DOSE. Overall P: nested likelihood-ratio test, 2 df for three-level components and 1 df otherwise.',
    'Contrast P: two-sided Wald test. BH FDR across 11 contrasts and across 9 overall tests, separately by model family, provided in tables.',
    'Panel a: unadjusted Kaplan-Meier IS-free survival and pointwise 95% CI; death censored; not competing-risk cumulative incidence.',
    'Follow-up and covariates match primary analysis. Curves truncated at 15 years; Cox uses all primary follow-up.',
    'Panel b: adjusted categorical HRs. Existing continuous-score component HRs are retained separately as supplementary analyses.',
    'Original component correlations, contributions and leave-one-out definitions remain based on 0 / 0.5 / 1 scores.'),
    file.path(OUT,'component_categorical_methods.md'))
  stopifnot(nrow(rr)==22,nrow(oo)==18,all(is.finite(rr$HR)),all(vapply(counts,function(x)sum(x$N)==nrow(d),logical(1))))
  logmsg('Component categorical models and figures complete')
}
run_dose9_components <- function(d) {
  path <- file.path(ROOT,'ukb/04score/dose9_score.rds')
  s <- readRDS(path); vc <- attr(s,'component_columns')
  stopifnot(identical(attr(s,'score_version'),'dose9'),length(vc)==9,!anyDuplicated(s$eid))
  ix <- match(d$eid,as.character(s$eid)); stopifnot(!anyNA(ix))
  z <- d; z[vc] <- s[ix,vc,drop=FALSE]
  changed <- c('Fruit.dose','Vegetables.dose'); audit <- list(); rows <- list()
  for(v in vc) {
    orig <- sub('[.]dose9$','.dose',v)
    expected <- if(orig %in% changed)as.numeric(d[[orig]]>0) else as.numeric(d[[orig]])
    stopifnot(all(z[[v]] %in% 0:1),identical(as.numeric(z[[v]]),expected))
    audit[[v]] <- data.frame(component=v,original_component=orig,
      changed_values=sum(z[[v]]!=d[[orig]]),N=nrow(z))
    logmsg('DOSE9 component M2:',v)
    f <- fit(z,c(v,co),paste(v,'M2'))
    rows[[v]] <- cbind(data.frame(component=v,model='M2',comparison='1 vs 0',
      N=f$n,events=f$nevent,N0=sum(z[[v]]==0),events0=sum(z$event[z[[v]]==0]),
      N1=sum(z[[v]]==1),events1=sum(z$event[z[[v]]==1])),effect(f,v))
  }
  a <- do.call(rbind,rows); a$P_FDR <- p.adjust(a$P,'BH')
  a$HR_95CI <- sprintf('%.3f (%.3f, %.3f)',a$HR,a$LCI,a$HCI)
  csv(a,'Table_DOSE9_component_M2'); csv(do.call(rbind,audit),'dose9_component_encoding_QC')
  csv(data.frame(path=path,md5=unname(tools::md5sum(path)),covariates=paste(co,collapse=' + ')),
    'dose9_component_model_specification')
  oldpath <- file.path(OUT,'component_IS_associations.csv')
  if(file.exists(oldpath)) {
    old <- read.csv(oldpath); old <- old[old$model=='M2',]
    j <- match(sub('[.]dose9$','.dose',a$component),old$component)
    same <- !sub('[.]dose9$','.dose',a$component) %in% changed
    stopifnot(all.equal(a$HR[same],old$HR[j[same]],tolerance=1e-8),
      all.equal(a$P[same],old$P[j[same]],tolerance=1e-8))
  }
  plotfun <- function() {
    par(mar=c(5,13,3,11),xpd=NA)
    yy <- 9:1
    plot(a$HR,yy,log='x',xlim=c(.65,1.1),ylim=c(.5,10.3),yaxt='n',ylab='',
      xlab='M2-adjusted HR (95% CI): 1 vs 0 point',main='DOSE9 component associations',pch=19,col='#176B91')
    segments(a$LCI,yy,a$HCI,yy,col='#176B91'); abline(v=1,lty=2,col='grey50')
    axis(2,at=yy,labels=gsub('_',' ',sub('[.]dose9$','',a$component)),las=1,tick=FALSE)
    text(1.12,yy,paste0(a$HR_95CI,'; P=',format.pval(a$P,digits=3,eps=.001)),adj=0,cex=.8)
    text(1.12,10,'HR (95% CI); unadjusted P',adj=0,cex=.8,font=2)
  }
  pdf(file.path(OUT,'figures/dose9_component_forest_M2.pdf'),width=12,height=7); plotfun(); bounded_dev_off()
  bounded_png(file.path(OUT,'figures/dose9_component_forest_M2.png'),width=12,height=7,units='in',res=180); plotfun(); bounded_dev_off()
  writeLines(c('DOSE9 component analysis: exploratory comparison using the same primary cohort and M2 covariates.',
    'Fruit and vegetables use OR aggregation: original 0.5 and 1 both map to 1. Other seven components unchanged.',
    'Separate Cox model per component, without total score or other components; all comparisons are 1 vs 0.',
    'BH FDR calculated across the nine component contrasts; no threshold changed based on outcome significance.',
    sprintf('Nominal P<0.05: %s/9; FDR<0.05: %s/9.',sum(a$P<.05),sum(a$P_FDR<.05)),
    'Encoding checked against the saved DOSE9 score; seven unchanged components reproduce previous M2 estimates.'),
    file.path(OUT,'dose9_component_results.txt'))
  logmsg('DOSE9 component models and verification complete')
}
if(DOSE9_COMPONENTS_ONLY) {
  if(file.exists(file.path(OUT,'analysis_warnings.csv')))issues <- read.csv(file.path(OUT,'analysis_warnings.csv'))
  run_dose9_components(d)
  quit(save='no',status=0)
}
if(COMPONENTS_ONLY) {
  if(file.exists(file.path(OUT,'analysis_warnings.csv'))) issues <- read.csv(file.path(OUT,'analysis_warnings.csv'))
  run_component_categories(d)
  quit(save='no',status=0)
}
if(PRS_ONLY) {
  if(file.exists(file.path(OUT,'analysis_warnings.csv')))
    issues <- read.csv(file.path(OUT,'analysis_warnings.csv'))
  run_prs_analysis(d)
  quit(save='no',status=0)
}

if(!DIETARY_ONLY) {
# 1. Component explanation: rules are imported from the scoring implementation.
env <- new.env(); options(dose.score.run=FALSE)
dose_source_step("score", "build", envir=env)
rules <- env$dose_food_rules(meta$primary_scoring_method)
rules$criterion_for_one <- vapply(seq_len(nrow(rules)),function(i)
  env$dose_interval_labels(rules[i,])$One,character(1))
rules$direction <- ifelse(rules$Type=='adequacy','Higher intake qualifies',
  ifelse(rules$Type=='restriction','Zero intake qualifies','Only the specified intermediate intake qualifies'))
csv(rules,'component_food_rules')
cs <- do.call(rbind,lapply(comps,function(v) {
  x <- d[[v]]; rest <- d$dose-x
  data.frame(component=v,mean=mean(x),SD=sd(x),proportion_0=mean(x==0),
    proportion_half=mean(x==.5),proportion_1=mean(x==1),
    correlation_total=cor(x,d$dose),correlation_rest=cor(x,rest),
    mean_score_share=mean(x)/mean(d$dose),
    variance_share=cov(x,d$dose)/var(d$dose))
}))
stopifnot(abs(sum(cs$variance_share)-1)<1e-10,abs(sum(cs$mean_score_share)-1)<1e-10)
csv(cs,'component_contributions')
for(method in c('pearson','spearman')) {
  cm <- cor(d[c(comps,'dose')],method=method)
  csv(data.frame(variable=rownames(cm),cm,check.names=FALSE),paste0('component_correlation_',method))
}
component_results <- list()
for(v in comps) {
  for(mode in c('M2','M2_mutually_adjusted')) {
    cov <- c(co,if(mode=='M2_mutually_adjusted') setdiff(comps,v))
    f <- fit(d,c(v,cov),paste(v,mode))
    component_results[[paste(v,mode)]] <- cbind(data.frame(component=v,model=mode,
      contrast='1 versus 0 score (linear for half-point components)',N=f$n,events=f$nevent),effect(f,v))
  }
}
cr <- do.call(rbind,component_results)
cr$P_FDR <- ave(cr$P,cr$model,FUN=function(x)p.adjust(x,'BH'))
csv(cr,'component_IS_associations')
run_component_categories(d)
run_dose9_components(d)
main <- run(d,'Primary_M2')
for(v in comps) {
  z <- d; z$loo <- z$dose-z[[v]]
  # All LOO models keep M2 for arithmetic comparability; removal of wine is
  # additionally evaluated with alcohol adjustment below.
  run(z,paste0('LOO_',sub('.dose$','',v)),x='loo',fixed_scale=FALSE)
}
wine <- comps[grepl('Wine',comps)]
if(length(wine)==1) {
  z <- d; z$no_wine <- z$dose-z[[wine]]
  run(z,'No_wine_alcohol_adjusted',covs=unique(c(co,cfg$covs.alcohol)),x='no_wine',fixed_scale=FALSE)
}
save_models()

# 3. Sensitivity analyses. All subset analyses use the original DOSE increment
# and original quartile boundaries; alternative scoring/LOO uses its own scale.
run(d,'M3_medications',covs=co3)
run(d,'M2_alcohol_adjusted',covs=unique(c(co,cfg$covs.alcohol)))
run(d,'M2_CVD_adjusted',covs=unique(c(co,'cvd')))
run(d,'Region_stratified',extra_terms='strata(region)')
for(lag in c(2,5)) {
  z <- d[d$time>lag,]; z$time <- z$time-lag
  run(z,paste0('Landmark_',lag,'y'))
}
for(v in c('cvd','primary_ahtn','t2dm')) run(d[which(d[[v]]=='No'),],paste0('Exclude_',v))
run(d[which(d$cvd=='No' & d$primary_ahtn=='No' & d$t2dm=='No'),],'Disease_free')
run(d[which(d$enroll>=2),],'At_least_2_diet_assessments')
run(d,'Additional_DASH',covs=c(co,'diet.dash.sum'))
run(d,'Additional_AHEI',covs=c(co,'diet.ahei.sum'))
for(v in c('dose11','dose9')) run(d,paste0('Alternative_',v),x=v,fixed_scale=FALSE)
z <- dat[complete.cases(dat[c('dose','time','event',co)]) & is.finite(dat$time) & dat$time>0,]
run(z,'DOSE_only_available_cohort')
rm(z,dat); invisible(gc())
raw <- as.data.frame(readRDS(files['raw']))
stopifnot(!anyDuplicated(raw$eid),all(d$eid %in% as.character(raw$eid)))
z <- d; z[co] <- raw[match(d$eid,as.character(raw$eid)),co,drop=FALSE]
# Unknown categories are also excluded in the complete-case sensitivity.
for(v in co) if(is.factor(z[[v]]) || is.character(z[[v]]))
  z[[v]][tolower(as.character(z[[v]])) %in% c('unknown','missing')] <- NA
run(z,'Complete_case_original_covariates')
rm(raw,z); invisible(gc()); save_models()
run_dietary_adjustments(d)

# 4. Subgroups: same DOSE increment and boundaries in every stratum. LRT tests
# interaction in a pooled model, not differences between subgroup P values.
sg <- list(age60=list(value=factor(ifelse(d$age<60,'<60','>=60')),drop=character()),
           sex=list(value=d$sex.b,drop='sex.b'),
           BMI30=list(value=factor(ifelse(d$bmi<30,'<30','>=30')),drop=character()),
           smoking=list(value=d$smoke.c,drop='smoke.c'),
           hypertension=list(value=d$primary_ahtn,drop='primary_ahtn'),
           diabetes=list(value=d$t2dm,drop='t2dm'),
           CVD=list(value=d$cvd,drop='cvd'),
           ethnicity=list(value=d$ethnic.c3,drop='ethnic.c3'),
           diet_assessments=list(value=factor(ifelse(d$enroll<2,'1','>=2')),drop=character()))
if('prs_istroke' %in% names(d)) {
  pb <- unname(quantile(d$prs_istroke,c(0,1/3,2/3,1),na.rm=TRUE))
  if(!anyDuplicated(pb)) sg$PRS_tertiles <- list(value=cut(d$prs_istroke,pb,include.lowest=TRUE,
      labels=c('Low','Middle','High')),drop=character(),add=paste0('PC',1:10))
}
interactions <- list(); subgroup_qc <- list()
for(n in names(sg)) {
  z <- d; z$subgroup <- factor(sg[[n]]$value)
  cov <- unique(c(setdiff(co,sg[[n]]$drop),sg[[n]]$add))
  z <- prepare(z,c('subgroup','dose',cov)); z$exposure <- z$dose
  cnt <- table(z$subgroup); ev <- tapply(z$event,z$subgroup,sum)
  eligible <- names(cnt)[cnt>=200 & ev>=20]
  subgroup_qc[[n]] <- data.frame(subgroup=n,level=names(cnt),N=as.integer(cnt),events=as.integer(ev),
                                 estimable=names(cnt) %in% eligible)
  for(g in eligible) run(z[z$subgroup==g,],paste('Subgroup',n,g,sep='__'),covs=cov)
  if(length(eligible)!=length(cnt) || length(eligible)<2) {
    note(n,'Interaction omitted: insufficient events in one or more levels');next
  }
  cov <- varying(z,cov)
  for(ex in c('exposure','group')) {
    f0 <- fit(z,c(ex,'subgroup',cov),paste(n,ex,'no interaction'))
    f1 <- fit(z,c(paste0(ex,'*subgroup'),cov),paste(n,ex,'interaction'))
    df <- length(coef(f1))-length(coef(f0))
    interactions[[paste(n,ex)]] <- data.frame(subgroup=n,exposure=ex,N=f1$n,events=f1$nevent,
      df=df,LRT=2*(as.numeric(logLik(f1))-as.numeric(logLik(f0))),
      P_interaction=pchisq(max(0,2*(as.numeric(logLik(f1))-as.numeric(logLik(f0)))),df,lower.tail=FALSE))
  }
}
it <- do.call(rbind,interactions)
it$P_interaction_FDR <- ave(it$P_interaction,it$exposure,FUN=function(x)p.adjust(x,'BH'))
csv(it,'subgroup_interactions'); csv(do.call(rbind,subgroup_qc),'subgroup_eligibility')
save_models()
ct <- do.call(rbind,continuous); qt <- do.call(rbind,categorical)
csv(ct[grepl('^Subgroup__',ct$analysis),],'subgroup_continuous')
csv(qt[grepl('^Subgroup__',qt$analysis),],'subgroup_quartiles')
csv(ct[grepl('^LOO_',ct$analysis),],'leave_one_component_out')

# 2. Marginal standardized risks. Reconstruct death before IS/censoring from
# dates and verify consistency against the original outcome follow-up.
run_prs_analysis(d)
logmsg('Standardized absolute risk with death as a competing event')
stopifnot(all(c('date_death','start_date','end0') %in% names(d)))
death_date <- as.Date(d$date_death); start <- as.Date(d$start_date)
end <- as.Date(d$end0)
death_years <- as.numeric(death_date-start)/365.25
# The legacy t2e() helper gives a nonmissing loss date precedence over death,
# rather than taking the earliest date. Bound every subject by end0, which is
# already min(loss, death, administrative end) in foodwas / cohort.
end_years <- as.numeric(end-start)/365.25
old_time <- d$time; old_event <- d$event
d$event <- as.integer(old_event==1 & old_time<=end_years+1e-8)
d$time <- pmin(old_time,end_years)
stopifnot(all(is.finite(d$time)),all(d$time>0))
csv(data.frame(metric=c('Followup_times_corrected','IS_events_recensored'),
               N=c(sum(abs(old_time-d$time)>1e-8),sum(old_event!=d$event))), 'followup_date_QC')
if(any(abs(old_time-d$time)>1e-8) || any(old_event!=d$event)) {
  note('Followup dates',paste(sum(abs(old_time-d$time)>1e-8),
    'follow-up times corrected to earliest loss/death/administrative end for risk estimation; original association models retained.'))
  run(d,'Earliest_censoring_date_corrected')
  save_models()
  ct <- do.call(rbind,continuous); qt <- do.call(rbind,categorical)
}
# A recorded stroke on the date of death is assigned to IS (same-day convention).
death_event <- d$event==0 & !is.na(death_date) & death_date<=end & death_date>start
stopifnot(!any(death_years[death_event]<=0),
          max(abs(death_years[death_event]-d$time[death_event]))<0.02)
stopifnot(!any(d$event==1 & !is.na(death_years) & death_years<d$time-0.02))
d$death_event <- as.integer(death_event)
csv(data.frame(status=c('Censored','IS','Death_before_IS'),
               N=c(sum(d$event==0 & !death_event),sum(d$event),sum(death_event))), 'competing_event_QC')
csv(do.call(rbind,lapply(0:10,function(t) do.call(rbind,lapply(levels(d$group),function(g)
  data.frame(year=t,group=g,N_at_risk=sum(d$group==g & d$time>=t)))))), 'risk_numbers_at_risk')
fg <- fit(d,c('group',co),'Risk_IS')
dd <- d; dd$event <- dd$death_event
fd_ph <- fit(dd,c('group',co),'Risk_death_PH_comparison')
make_death_stratum <- function(z) {
  band <- cut(z$age,c(-Inf,45,50,55,60,65,70,Inf),right=FALSE,labels=paste0('A',1:7))
  factor(paste(z$group,band,sep='_'))
}
fit_death <- function(z,context) {
  z$death_stratum <- make_death_stratum(z)
  z$event <- z$death_event
  f <- fit(z,c(co,'strata(death_stratum)'),context)
  attr(f,'death_strata_levels') <- levels(z$death_stratum)
  f
}
fd <- fit_death(d,'Risk_death_DOSE_age_stratified')
risk_ph_table <- function(fi,fd) do.call(rbind,lapply(c('IS','death'),function(n) {
  z <- as.data.frame(cox.zph(if(n=='IS') fi else fd)$table)
  data.frame(model=n,term=rownames(z),z,row.names=NULL)
}))
csv(risk_ph_table(fg,fd_ph),'risk_model_PH_tests_unstratified')
csv(risk_ph_table(fg,fd),'risk_model_PH_tests')
note('Competing death model','Death baseline hazards stratified by DOSE quartile and baseline age bands (<45, 45-49, 50-54, 55-59, 60-64, 65-69, >=70); continuous age retained within bands. See unstratified PH diagnostic and risk comparison.')
# Cox objects keep a local data symbol; newdata is supplied explicitly whenever
# a model matrix is needed. Baseline hazards use the stored response/LP.
bh <- function(f,times,stratum=NULL) {
  b <- basehaz(f,centered=TRUE)
  if(!is.null(stratum)) {
    stopifnot('strata' %in% names(b),stratum %in% as.character(b$strata))
    b <- b[as.character(b$strata)==stratum,,drop=FALSE]
  }
  c(0,b$hazard)[findInterval(times,b$time)+1L]
}
# Integrate both cause-specific hazards on a fine grid using constant hazards
# within each interval. Exact population averaging, with a grid-convergence QC.
make_risk_kernel <- function() {
  # Optional speedup only. The original base-R algorithm remains the fallback.
  .libPaths(unique(c(file.path(ROOT,'.r-lib'),.libPaths())))
  if(!requireNamespace('Rcpp',quietly=TRUE)) return(NULL)
  oldpath <- Sys.getenv('PATH'); on.exit(Sys.setenv(PATH=oldpath),add=TRUE)
  if(.Platform$OS.type=='windows') {
    rt <- Sys.getenv('RTOOLS45_HOME','D:/software/win/rtools45')
    if(dir.exists(rt)) Sys.setenv(PATH=paste(file.path(rt,'usr/bin'),
      file.path(rt,'x86_64-w64-mingw32.static.posix/bin'),oldpath,sep=.Platform$path.sep))
  }
  en <- new.env()
  tryCatch({
    Rcpp::cppFunction(env=en,code='
      Rcpp::NumericMatrix dose_risk_kernel(Rcpp::NumericVector ri,
          Rcpp::NumericVector rd, Rcpp::NumericVector hi,
          Rcpp::NumericMatrix hd, Rcpp::IntegerVector si) {
        int n=ri.size(), nt=hi.size();
        Rcpp::NumericMatrix out(nt,2);
        std::vector<double> s(n,1.0), f(n,0.0);
        for(int j=1;j<nt;j++) {
          double sumf=0.0, sumn=0.0, dh=hi[j]-hi[j-1];
          for(int i=0;i<n;i++) {
            double a=dh*ri[i], b=(hd(j,si[i])-hd(j-1,si[i]))*rd[i];
            double h=a+b;
            if(h>0) { double e=std::exp(-h); f[i]+=s[i]*(1-e)*a/h; s[i]*=e; }
            sumf+=f[i]; sumn+=-std::expm1(-hi[j]*ri[i]);
          }
          out(j,0)=sumf/n; out(j,1)=sumn/n;
          if(j%20==0) Rcpp::checkUserInterrupt();
        }
        return out;
      }',plugins='cpp11')
    en$dose_risk_kernel
  },error=function(e) {note('Risk kernel',paste('Using base R fallback:',conditionMessage(e)));NULL})
}
risk_kernel <- make_risk_kernel()
# Analytic no-competing-death identity checks the compiled recursion.
if(!is.null(risk_kernel)) {
  test_hi <- c(0,.01,.03,.1); test_ri <- c(.4,1,2)
  kk <- risk_kernel(test_ri,c(1,1,1),test_hi,matrix(0,4,1),rep(0L,3))
  stopifnot(max(abs(kk[,1]-vapply(test_hi,function(h)mean(-expm1(-h*test_ri)),numeric(1))))<1e-12)
  # Independent base-R competing-hazard recursion, including two strata.
  hd_test <- cbind(c(0,.03,.07,.12),c(0,.01,.04,.09))
  ri_test <- c(.4,1,2); rd_test <- c(.8,1.3,.6); si_test <- c(0L,1L,0L)
  kk <- risk_kernel(ri_test,rd_test,test_hi,hd_test,si_test)
  ss <- rep(1,3); ff <- numeric(3); expected <- numeric(4)
  for(j in 2:4) {
    a <- diff(test_hi)[j-1]*ri_test
    b <- (hd_test[j,]-hd_test[j-1,])[si_test+1L]*rd_test; h <- a+b
    ff <- ff+ss*(-expm1(-h))*a/h; ss <- ss*exp(-h); expected[j] <- mean(ff)
  }
  stopifnot(max(abs(kk[,1]-expected))<1e-12)
}
risks <- function(fi,fd,ref,step=.05) {
  tt <- seq(0,10,by=step)
  hi <- bh(fi,tt)
  strat_levels <- attr(fd,'death_strata_levels')
  if(is.null(strat_levels)) hd <- bh(fd,tt) else {
    hb <- basehaz(fd,centered=TRUE)
    stopifnot(setequal(unique(as.character(hb$strata)),strat_levels))
    hd <- sapply(strat_levels,function(s) {
      b <- hb[as.character(hb$strata)==s,,drop=FALSE]
      c(0,b$hazard)[findInterval(tt,b$time)+1L]
    })
  }
  out <- list()
  for(g in levels(ref$group)) {
    nd <- ref; nd$group <- factor(g,levels=levels(ref$group))
    if(!is.null(strat_levels)) {
      nd$death_stratum <- factor(make_death_stratum(nd),levels=strat_levels)
      stopifnot(!anyNA(nd$death_stratum))
      si <- match(as.character(nd$death_stratum),strat_levels)
    }
    ri <- exp(predict(fi,newdata=nd,type='lp',reference='sample'))
    rd <- exp(predict(fd,newdata=nd,type='lp',reference='sample'))
    if(!is.null(risk_kernel)) {
      ans <- risk_kernel(ri,rd,hi,if(is.null(strat_levels)) matrix(hd,ncol=1) else hd,
                         if(is.null(strat_levels)) rep(0L,nrow(nd)) else as.integer(si-1L))
      cifmean <- ans[,1]; netmean <- ans[,2]
    } else {
      surv <- rep(1,nrow(nd)); cif <- numeric(nrow(nd))
      cifmean <- netmean <- numeric(length(tt))
      for(j in 2:length(tt)) {
      a <- (hi[j]-hi[j-1])*ri
      b <- if(is.null(strat_levels)) (hd[j]-hd[j-1])*rd else (hd[j,]-hd[j-1,])[si]*rd
      total <- a+b
      frac <- numeric(length(total)); nz <- total>0; frac[nz] <- a[nz]/total[nz]
      cif <- cif+surv*(-expm1(-total))*frac; surv <- surv*exp(-total)
      cifmean[j] <- mean(cif); netmean[j] <- mean(-expm1(-hi[j]*ri))
      }
    }
    out[[g]] <- data.frame(group=g,year=tt,CIF=cifmean,net_risk=netmean)
  }
  do.call(rbind,out)
}
risk <- risks(fg,fd,d,.05); fine <- risks(fg,fd,d,.025)
chk <- merge(risk,fine,by=c('group','year'),suffixes=c('_coarse','_fine'))
err <- max(abs(chk$CIF_coarse-chk$CIF_fine))
csv(data.frame(coarse_years=.05,fine_years=.025,max_absolute_risk_difference=err), 'risk_integration_QC')
stopifnot(err<1e-5,all(fine$CIF>=0),all(fine$CIF<=fine$net_risk+1e-10),
          all(vapply(split(fine$CIF,fine$group),function(x) all(diff(x)>=-1e-12),logical(1))))
risk <- fine; risk10 <- risk[risk$year==10,]; rownames(risk10) <- NULL
risk_ph <- risks(fg,fd_ph,d,.05)
csv(risk_ph,'adjusted_cumulative_incidence_unstratified_death')
cmp <- merge(risk10,risk_ph[risk_ph$year==10,],by=c('group','year'),suffixes=c('_stratified','_PH'))
cmp$risk_difference_between_death_models <- cmp$CIF_stratified-cmp$CIF_PH
csv(cmp,'death_model_10year_risk_comparison')
risk10$risk_difference_vs_Q1 <- risk10$CIF-risk10$CIF[1]
risk10$risk_ratio_vs_Q1 <- risk10$CIF/risk10$CIF[1]
risk10$CI_method <- 'Not estimated (point estimates only)'
# Optional nonparametric participant bootstrap, refitting BOTH models and
# restandardizing over each resampled cohort. Do not use coefficient-only draws
# as if they included baseline-hazard uncertainty.
B <- as.integer(Sys.getenv('DOSE_RISK_BOOT','0')); stopifnot(!is.na(B),B==0 || B>=100)
if(B>0) {
  set.seed(20260919); boots <- vector('list',B)
  for(b in seq_len(B)) {
    logmsg('Absolute-risk bootstrap',b,'/',B)
    zz <- d[sample.int(nrow(d),replace=TRUE),]; rownames(zz) <- NULL
    bi <- fit(zz,c('group',co),paste('bootstrap IS',b))
    bd <- fit_death(zz,paste('bootstrap death',b))
    rr <- risks(bi,bd,zz,.05); rr <- rr[rr$year==10,]
    rr$RD <- rr$CIF-rr$CIF[1]; rr$RR <- rr$CIF/rr$CIF[1]; rr$replicate <- b
    boots[[b]] <- rr; csv(do.call(rbind,boots[seq_len(b)]),'risk_bootstrap')
  }
  bt <- do.call(rbind,boots)
  for(v in c('CIF','RD','RR')) for(k in 1:4) {
    ci <- quantile(bt[bt$group==risk10$group[k],v],c(.025,.975))
    risk10[k,paste0(v,c('_LCI','_HCI'))] <- ci
  }
  risk10$CI_method <- paste(B,'participant bootstrap percentile replicates')
}
csv(risk,'adjusted_cumulative_incidence'); csv(risk10,'adjusted_10year_risk')
model_summary <- function(f) list(coefficients=coef(f),variance=vcov(f),
  baseline_hazard=basehaz(f,centered=TRUE),centering=f$means,N=f$n,events=f$nevent,
  death_strata_levels=attr(f,'death_strata_levels'))
dose_save_rds(list(IS=model_summary(fg),death=model_summary(fd),reference_N=nrow(d),covariates=co,group_breaks=br,
             method=meta$primary_scoring_method),file.path(OUT,'risk_models.rds'))

# Compact vector + raster figures; no participant-level files are exported.
draw <- function(name,fun,w=9,h=6) {
  pdf(file.path(OUT,'figures',paste0(name,'.pdf')),width=w,height=h); fun(); bounded_dev_off()
  bounded_png(file.path(OUT,'figures',paste0(name,'.png')),width=w,height=h,units='in',res=180); fun(); bounded_dev_off()
}
labs <- sub('[.]dose$','',comps); labs <- gsub('_',' ',labs)
draw('component_correlation',function() {
  cm <- cor(d[c(comps,'dose')]); n <- ncol(cm)
  par(mar=c(10,10,3,1)); image(1:n,1:n,cm,zlim=c(-1,1),
    col=colorRampPalette(c('#2166AC','white','#B2182B'))(101),axes=FALSE,xlab='',ylab='')
  axis(1,1:n,c(labs,'Total DOSE'),las=2,cex.axis=.8);axis(2,1:n,c(labs,'Total DOSE'),las=2,cex.axis=.8)
  for(i in 1:n) for(j in 1:n) text(i,j,sprintf('%.2f',cm[i,j]),cex=.65)
  title('Component and total-score Pearson correlations')
},9,9)
draw('component_contribution',function() {
  par(mfrow=c(1,2),mar=c(5,10,3,1))
  barplot(100*cs$mean_score_share,names.arg=labs,horiz=TRUE,las=1,col='#4C859D',
    xlab='Share of mean total score (%)',main='Mean-score contribution',cex.names=.8)
  barplot(100*cs$variance_share,names.arg=labs,horiz=TRUE,las=1,col='#B87355',
    xlab='Cov(component, total) / Var(total) (%)',main='Variance decomposition',cex.names=.8)
},13,6)
forest <- function(tab,title) {
  tab <- tab[is.finite(tab$LCI) & is.finite(tab$HCI),]; n <- nrow(tab)
  par(mar=c(5,19,3,2)); yy <- rev(seq_len(n))
  plot(tab$HR,yy,log='x',xlim=range(c(.7,1.1,tab$LCI,tab$HCI)),ylim=c(.5,n+.5),
    pch=19,col='#27647B',yaxt='n',ylab='',xlab='Adjusted hazard ratio (95% CI)',main=title)
  segments(tab$LCI,yy,tab$HCI,yy,col='#27647B');abline(v=1,lty=2,col='grey50')
  axis(2,yy,tab$analysis,las=1,cex.axis=.72)
}
draw('leave_one_component_out',function() forest(ct[ct$scale=='Per 1-SD' &
  (ct$analysis=='Primary_M2' | grepl('^LOO_',ct$analysis)),], 'DOSE: leave-one-component-out (per score-specific SD)'),10,7)
ss <- ct[ct$scale=='P10-P90' & !grepl('^(LOO_|Subgroup__)',ct$analysis),]
draw('sensitivity_forest',function() forest(ss,'DOSE sensitivity analyses (P10-P90 increment)'),11,9)
ssg <- ct[ct$scale=='P10-P90' & grepl('^Subgroup__',ct$analysis),]
draw('subgroup_forest',function() forest(ssg,'DOSE subgroup analyses (common P10-P90 increment)'),12,10)
pal <- c('#4477AA','#66AABB','#CC9955','#AA4466')
draw('adjusted_cumulative_incidence',function() {
  par(mar=c(5,5,3,1)); plot(NA,xlim=c(0,10),ylim=c(0,max(risk$CIF)*105),
    xlab='Years since baseline',ylab='Adjusted cumulative incidence of IS (%)',
    main='Population-standardized incidence; competing death included')
  for(k in 1:4) {r <- risk[risk$group==paste0('Q',k),];lines(r$year,100*r$CIF,col=pal[k],lwd=2)}
  legend('topleft',paste0('Q',1:4),col=pal,lwd=2,bty='n')
},8,6)
draw('adjusted_10year_risk',function() {
  par(mar=c(5,5,3,1)); barplot(100*risk10$CIF,names.arg=risk10$group,col=pal,
    ylim=c(0,max(risk10$CIF)*120),ylab='Adjusted 10-year IS risk (%)',main='Population-standardized 10-year risk')
  mtext('Point estimates; death treated as a competing event',side=1,line=3)
},7,5)
csv(issues,'analysis_warnings')
capture.output(sessionInfo(),file=file.path(OUT,'sessionInfo.txt'))
writeLines(c('# DOSE / IS sensitivity analysis',
  paste('Completed:',Sys.time()),paste('Primary scoring method:',meta$primary_scoring_method),
  paste('Primary cohort:',nrow(d),'; IS events:',sum(d$event)),
  '', '## Methods and interpretation',
  '- Primary cohort and quartiles match score.R (common cohort for all dietary patterns).',
  '- M2 is read from analysis_config.rds. It includes age, sex, ethnicity, education, TDI, energy, smoking, PA, BMI, family history, hypertension and diabetes. Primary DOSE includes wine, so M2 omits alcohol; alcohol adjustment is a sensitivity analysis.',
  '- Higher component score means meeting the specified rule. Intermediate-range foods are nonmonotonic in intake. These research rules are not clinical dietary advice.',
  '- Fruit and vegetables each average two binary criteria in DOSE9half. Rules are exported at food level.',
  '- Component/total correlation includes mathematical part-whole correlation. Correlation with the remaining score is also reported. Variance shares sum to 1 and are descriptive, not causal importance.',
  '- Component effects use 1 vs 0 score; for half-point components this imposes a linear score effect. Mutual adjustment includes the other components, never the total score.',
  '- LOO changes the score and uses its own SD, P10-P90 and tie-preserving quartiles. Per-point HRs allow a common-unit comparison. The no-wine model additionally adjusts for alcohol.',
  '- Subsets retain the full-source-cohort DOSE SD/P10-P90 and quartile boundaries. Landmark analyses remove everyone with follow-up <=2/5 years and reset the time origin, estimating conditional associations among survivors still under observation.',
  '- Complete-case analysis restores original pre-imputation covariates and excludes Unknown/Missing levels. Original preparation uses mean/mode imputation, not multiple imputation.',
  '- Adjusted absolute risks average predictions across every member of the primary cohort after setting everyone to each DOSE quartile. They are standardized associations, not established causal effects.',
  '- Death before IS is a competing event; IS on the date of death takes precedence. Two cause-specific Cox models produce the CIF. Net risk 1-S(t) is exported separately and does not account for competing death.',
  '- The death model uses separate baseline hazards for DOSE quartile x baseline age bands (<45,45-49,50-54,55-59,60-64,65-69,>=70), retaining continuous age and other M2 covariates. This addresses the strong DOSE/age PH violation in the unstratified death model. Both diagnostics and the 10-year risk comparison are exported.',
  '- Absolute-risk follow-up is bounded by end0=min(loss, death, administrative end). The legacy helper can prioritize a later loss date over an earlier death; followup_date_QC records corrections and an additional corrected-censoring Cox sensitivity is exported. Original primary-model results are preserved for comparability.',
  '- CIF integration uses a 0.025-year grid and is checked against 0.05 years (maximum permitted absolute discrepancy 0.00001).',
  '- If local Rcpp and a compiler are available, an equivalent C++ integration kernel accelerates population averaging; otherwise base R is used. The kernel is checked against analytic net risk and an independent competing-risk recursion.',
  paste('- Absolute-risk intervals:',if(B==0) 'Not estimated in the default run. These are point estimates; do not describe them as 95% CIs. Set DOSE_RISK_BOOT=200 or higher to run participant-bootstrap intervals.' else paste(B,'participant-bootstrap replicates.')),
  '- Cox estimates include Wald 95% CIs and PH tests. PH violations must be considered when interpreting a single follow-up-averaged HR.',
  '- Subgroups retain continuous age/BMI within age/BMI strata. Constant covariates are removed. PRS tertiles additionally adjust for PC1-PC10. LRT interaction P values and within-family BH FDR are exported for continuous and quartile exposure.',
  '- Subgroup models require >=200 participants and >=20 IS events; quartile models require >=5 events in each cell. Eligibility and warnings explicitly record unavailable analyses.',
  '- Sensitivities test available prespecified alternatives; they do not validate outcome-selected score development in an independent population.',
  '', 'Method reference: https://stat.ethz.ch/CRAN/web/packages/survivalVignettes/vignettes/tutorial.html',
  'Legacy comparison: D:/projects/dose_v1/scripts/03_pattern.R',
  '', '## Files',
  '- component_food_rules, component_contributions, component_correlation_*, component_IS_associations',
  '- leave_one_component_out, sensitivity_continuous, sensitivity_quartiles',
  '- adjusted_cumulative_incidence, adjusted_10year_risk, risk_numbers_at_risk, competing_event_QC',
  '- subgroup_continuous, subgroup_quartiles, subgroup_interactions, subgroup_eligibility',
  '- model_specifications, scenario_cohort_QC, input_manifest, sessionInfo, analysis_warnings',
  '- figures/: PDF and PNG versions of seven figures'),file.path(OUT,'README.md'))
logmsg('COMPLETE:',OUT)
} else {
  read_previous <- function(n) read.csv(file.path(OUT,paste0(n,'.csv')),check.names=FALSE)
  ct <- read_previous('sensitivity_continuous'); qt <- read_previous('sensitivity_quartiles')
  split_ordered <- function(x) split(x,factor(x$analysis,levels=unique(x$analysis)))
  continuous <- split_ordered(ct); categorical <- split_ordered(qt)
  model_specs <- split_ordered(read_previous('model_specifications'))
  eligibility <- split_ordered(read_previous('scenario_cohort_QC'))
  issues <- read_previous('analysis_warnings')
  # Reuse completed LOO fits only with unchanged inputs and the expected nine domains.
  expected_loo <- paste0('LOO_',sub('[.]dose$','',comps))
  stopifnot(setequal(expected_loo,unique(ct$analysis[grepl('^LOO_',ct$analysis)])),
    all(ct$N[ct$analysis %in% expected_loo]==nrow(d)),
    all(ct$events[ct$analysis %in% expected_loo]==sum(d$event)),
    all(vapply(expected_loo,function(n)sum(qt$analysis==n)==4L,logical(1))))
  run_dietary_adjustments(d)
  save_models()
  logmsg('Dietary adjustment models complete; updating merged tables.')
}

# Review-ready summary tables, regenerated entirely from aggregate CSV outputs.
write_review_tables <- function(out) {
  old_locale <- Sys.getlocale('LC_CTYPE')
  on.exit(suppressWarnings(Sys.setlocale('LC_CTYPE',old_locale)),add=TRUE)
  suppressWarnings(Sys.setlocale('LC_CTYPE',if(.Platform$OS.type=='windows')
    'English_United States.utf8' else 'C.UTF-8'))
  read <- function(n) read.csv(file.path(out,paste0(n,'.csv')),check.names=FALSE)
  cc <- read('sensitivity_continuous'); qq <- read('sensitivity_quartiles')
  ii <- read('subgroup_interactions'); rr <- read('adjusted_10year_risk')
  fmt <- function(x) ifelse(is.na(x),'',sprintf('%.3f',x))
  hrfmt <- function(x) if(!nrow(x)) '' else sprintf('%.3f (%.3f, %.3f)',x$HR,x$LCI,x$HCI)
  rows <- lapply(unique(cc$analysis),function(a) {
    x <- cc[cc$analysis==a,]; q <- qq[qq$analysis==a,]
    s <- x[x$scale=='P10-P90',]; p <- x[x$scale=='Per 1-SD',]
    data.frame(Analysis=a,N=s$N,Events=s$events,P10_P90_increment=s$increment,
      HR_P10_P90=hrfmt(s),P=s$P,HR_perSD=hrfmt(p),P_perSD=p$P,Q1='Ref.',
      Q2_vs_Q1=hrfmt(q[q$group=='Q2',]),Q3_vs_Q1=hrfmt(q[q$group=='Q3',]),
      Q4_vs_Q1=hrfmt(q[q$group=='Q4',]),
      P_Q2=if(nrow(q))q$P[q$group=='Q2'] else NA_real_,
      P_Q3=if(nrow(q))q$P[q$group=='Q3'] else NA_real_,
      P_Q4=if(nrow(q))q$P[q$group=='Q4'] else NA_real_,
      P_trend=if(nrow(q))q$P_trend[1] else NA_real_,
      quartile_score_ranges=if(nrow(q))paste(paste0(q$group,': ',q$observed_min,'-',q$observed_max),collapse='; ') else '')
  })
  tab <- do.call(rbind,rows)
  main_order <- c('Primary_M2','M3_medications','M2_alcohol_adjusted','M2_CVD_adjusted',
    'Region_stratified','Landmark_2y','Landmark_5y','Exclude_cvd','Exclude_primary_ahtn','Exclude_t2dm',
    'Disease_free','At_least_2_diet_assessments','Complete_case_original_covariates',
    'DOSE_only_available_cohort','Earliest_censoring_date_corrected','Additional_DASH','Additional_AHEI',
    'Additional_sodium','Additional_potassium','Additional_fiber','Additional_saturated_fat',
    'Additional_four_nutrients','Additional_other_foodgroups','Additional_nutrients_and_other_foodgroups',
    'Alternative_dose11','Alternative_dose9','No_wine_alcohol_adjusted')
  sens <- tab[!grepl('^Subgroup__',tab$Analysis),]
  ord <- unique(c(main_order,setdiff(sens$Analysis,c(main_order,sens$Analysis[grepl('^LOO_',sens$Analysis)])),
                  sens$Analysis[grepl('^LOO_',sens$Analysis)]))
  sens <- sens[order(match(sens$Analysis,ord)),]
  sens$Section <- ifelse(grepl('^LOO_',sens$Analysis),'Leave-one-component-out','Sensitivity analyses')
  write.csv(sens,file.path(out,'Table_sensitivity.csv'),row.names=FALSE)
  # Legacy-style publication table: a single row per scenario, no Cartesian join.
  s13 <- data.frame(Analysis=sens$Analysis,Section=sens$Section,
    `Events/Total`=paste0(sens$Events,'/',sens$N),
    `P10-P90 increment`=sens$P10_P90_increment,
    `P10-P90 HR (95% CI)`=sens$HR_P10_P90,`P continuous`=sens$P,
    `Per-SD HR (95% CI)`=sens$HR_perSD,`P per-SD`=sens$P_perSD,
    Q1=sens$Q1,Q2=sens$Q2_vs_Q1,`P Q2`=sens$P_Q2,
    Q3=sens$Q3_vs_Q1,`P Q3`=sens$P_Q3,Q4=sens$Q4_vs_Q1,`P Q4`=sens$P_Q4,
    `P trend`=sens$P_trend,`Observed group score ranges`=sens$quartile_score_ranges,check.names=FALSE)
  write.csv(s13,file.path(out,'Table_S13_Sensitivity.csv'),row.names=FALSE,na='')
  write.csv(cc[cc$analysis %in% sens$Analysis,],file.path(out,'res_sensitivity_std.csv'),row.names=FALSE)
  write.csv(qq[qq$analysis %in% sens$Analysis,],file.path(out,'res_sensitivity_4g.csv'),row.names=FALSE)
  # Refresh the sensitivity forest using the same aggregate rows as the table.
  fx <- cc[cc$scale=='P10-P90' & cc$analysis %in% sens$Analysis & !grepl('^LOO_',cc$analysis),]
  fx <- fx[order(match(fx$analysis,sens$Analysis)),]
  plot_sens <- function() {
    n <- nrow(fx); yy <- rev(seq_len(n)); par(mar=c(5,20,3,1))
    plot(fx$HR,yy,log='x',xlim=range(c(.6,1.1,fx$LCI,fx$HCI)),ylim=c(.5,n+.5),
      pch=19,col='#27647B',yaxt='n',ylab='',xlab='Adjusted HR (95% CI), P10-P90 increment',
      main='DOSE sensitivity analyses')
    segments(fx$LCI,yy,fx$HCI,yy,col='#27647B'); abline(v=1,lty=2,col='grey50')
    axis(2,yy,gsub('_',' ',fx$analysis),las=1,cex.axis=.7)
  }
  height <- max(9,nrow(fx)*.31)
  pdf(file.path(out,'figures/sensitivity_forest.pdf'),width=12,height=height); plot_sens(); bounded_dev_off()
  bounded_png(file.path(out,'figures/sensitivity_forest.png'),width=12,height=height,units='in',res=180); plot_sens(); bounded_dev_off()
  st <- tab[grepl('^Subgroup__',tab$Analysis),]
  bits <- strsplit(st$Analysis,'__',fixed=TRUE)
  st$Subgroup <- vapply(bits,'[',character(1),2)
  st$Level <- vapply(bits,'[',character(1),3)
  for(e in c('exposure','group')) {
    z <- ii[ii$exposure==e,]; m <- match(st$Subgroup,z$subgroup)
    st[[paste0('P_interaction_',e)]] <- z$P_interaction[m]
    st[[paste0('P_interaction_FDR_',e)]] <- z$P_interaction_FDR[m]
  }
  write.csv(st,file.path(out,'Table_subgroups.csv'),row.names=FALSE)
  main <- cc[cc$analysis=='Primary_M2' & cc$scale=='P10-P90',]
  loo <- cc[grepl('^LOO_',cc$analysis) & cc$scale=='Per 1-SD',]
  sig <- ii[ii$exposure=='exposure' & ii$P_interaction_FDR<.05,]
  lines <- c('# DOSE 与缺血性卒中：补充分析结果',
    sprintf('主队列：%s 人，%s 例 IS；评分采用 DOSE9half，范围 0–9 分。',main$N,main$events),
    '',sprintf('- M2 中 DOSE 从 P10 增至 P90（增加 %.1f 分）的 HR 为 %s，P=%.3g。',main$increment,hrfmt(main),main$P),
    sprintf('- 逐一剔除 9 个组分后，每个新评分增加 1 SD 的 HR 范围为 %.3f–%.3f；具体置信区间见 leave_one_component_out.csv。',min(loo$HR),max(loo$HR)),
    '- 考虑死亡竞争风险、对主队列全体参与者标准化后的 10 年 IS 风险：',
    sprintf('  - %s：%.3f%%。',rr$group,100*rr$CIF),
    sprintf('- Q4 相对 Q1 的 10 年风险差为 %.3f 个百分点，风险比为 %.3f。',
            100*rr$risk_difference_vs_Q1[4],rr$risk_ratio_vs_Q1[4]),
    if(!nrow(sig)) '- 连续 DOSE 的分层交互检验经 BH 校正后，没有达到 FDR <0.05 的分层变量。' else
      paste('- 连续 DOSE 交互检验 FDR <0.05 的分层变量：',paste(sig$subgroup,collapse='、'),'。'),
    '', '## 使用说明',
    '- Table_sensitivity.csv 和 Table_subgroups.csv 可直接用于核对和整理补充表。完整数值、模型变量及样本筛选信息保存在各 CSV 中。',
    '- Table_S13_Sensitivity.csv 合并连续评分、Q1–Q4 分组、各比较 P 值及趋势检验，并将 9 项逐一剔除分析单列为一个部分。',
    '- 营养素调整使用钠 w.p26052、钾 w.p26024、纤维 w.p26017、饱和脂肪 w.p26014；旧代码 w.p26155 是反式脂肪。其他食物组名单见 other_foodgroup_covariate_QC.csv。',
    '- 组分贡献为评分的描述性分解，不代表某种食物对卒中的因果贡献。组分赋分方向与阈值见 component_food_rules.csv。',
    '- 逐一剔除和替代评分采用各自的 SD/P10–P90；分层及排除分析采用固定的原始 DOSE 增量与分组界值。',
    '- 前 2/5 年排除采用 landmark 设计，结果适用于该时间点仍未发生 IS、仍存活且仍在随访的人群。',
    '- 绝对风险使用死亡、失访、行政截止的最早日期；followup_date_QC.csv 记录旧随访时长的校正，并提供 Earliest_censoring_date_corrected 对照模型。',
    '- 死亡模型按 DOSE 四分位组和年龄带分层基线风险，并保留年龄连续项；原死亡模型的对照估计见 death_model_10year_risk_comparison.csv。',
    if(all(grepl('Not estimated',rr$CI_method))) '- 本次绝对风险为点估计，未计算其置信区间；HR 均提供 95% CI。设置 DOSE_RISK_BOOT=200 后重跑可估计风险及风险差/比的 bootstrap 区间。' else paste('- 绝对风险区间：',rr$CI_method[1]),
    '- 检查 analysis_warnings.csv、subgroup_eligibility.csv 及 PH_P。统计显著性或亚组内 P 值差异不能代替交互检验。',
    '- 这些分析不能替代独立人群对评分的外部验证。')
  added <- c(Additional_sodium='调整钠',Additional_potassium='调整钾',
    Additional_fiber='调整纤维',Additional_saturated_fat='调整饱和脂肪',
    Additional_four_nutrients='同时调整四种营养素',
    Additional_other_foodgroups='同时调整其他食物组',
    Additional_nutrients_and_other_foodgroups='同时调整四种营养素和其他食物组')
  if(all(names(added) %in% sens$Analysis)) {
    z <- sens[match(names(added),sens$Analysis),]
    lines <- c(lines,'','## 新增营养素及其他食物组调整',
      '各模型保持原 DOSE P10–P90 增量（4.5 分）和四分位界值。',
      '| 分析 | P10–P90 HR（95% CI） | Q4 对 Q1 HR（95% CI） |',
      '|---|---|---|',sprintf('| %s | %s | %s |',unname(added),z$HR_P10_P90,z$Q4_vs_Q1),
      '营养素和其他食物组可能反映评分的部分组成或作用路径，因此上述结果属于进一步/过度调整检验，不能直接解释为中介分析或完全独立的因果效应。',
      '饱和脂肪字段为 w.p26014；旧参考代码使用的 w.p26155 是反式脂肪，已纠正。',
      '参考字段：https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=100117')
  }
  writeLines(enc2utf8(lines),file.path(out,'RESULTS_zh.md'),useBytes=TRUE)
}
write_review_tables(OUT)
}))

# ===========================================================================
# score / plot
# ===========================================================================
dose_register_step("score", "plot", expression({
# ======================================================================
# Figures/tables for score.R. Run this script AFTER that analysis.
# References: dose_v1/scripts/03_plot_mutiPattern.R and the DOSE palette.
# Reads completed results and dat.pattern.rds for score correlations and Table 1.
# Does not refit survival models or select cutoffs. The M2 forest recalculates
# BH FDR within the ten manuscript patterns; other figures retain source FDR.
# --forest-10-only regenerates only this forest and its ten-pattern result table.
# --table1-only generates Table 1 in outputs/ukb_table1 without drawing figures.
# All primary-score analyses and plot labels use dose / DOSE.
# Outputs: ukb/04score/pattern/figures/; see README_figures.txt.
# R dependencies: base R, grid, graphics and grDevices. Table 1 XLSX formatting
# optionally uses f/table1_export.mjs and the bundled Node artifact tool.
#
# 与最新 s3_2 输出的对应关系（不再读取旧版 res_std/res_q5/res_4g）：
#   diet = "dose"；原始评分列由 pattern_spec.csv 指定为 dose。
#   评分范围和显示步长读取结果元数据，不写死具体计分版本。
#   exposure/level -> diet/category；m1-m3 -> M1-M3；五分位 -> Q1-Q4。
#   case/total -> Events/N；分组计数使用 Group_Events/Group_N。
#   .std -> 按 scale 区分 "P10-P90" 与 "Per 1-SD"，不可混合。
#   饮食纳入与标签读取 pattern_spec.csv；结局标签读取 analysis_config.rds。
# 旧版各图在下方第 4 节整合；保留第 1-3 节的 DOSE/IS 主分析图表。
# 此脚本运行时只读取已完成的模型结果；不会 source 上游分析脚本。
# ======================================================================

DIR_PROJECT <- dose_output_dir()
DIR_OUT <- file.path(DIR_PROJECT, "ukb/04score/pattern")
DIR_FIG <- file.path(DIR_OUT, "figures")
read_result <- function(name) {
    f <- file.path(DIR_OUT, name)
    if (!file.exists(f)) stop("Missing ", f, ". Run score.R first.")
    read.csv(f, check.names=FALSE, stringsAsFactors=FALSE)
}
cont <- read_result("pattern_continuous.csv")
cat_res <- read_result("pattern_categories.csv")
tests <- read_result("pattern_overall_tests.csv")
rcs <- read_result("DOSE_RCS_curves.csv")
rcs_tests <- read_result("DOSE_RCS_tests.csv")
comparison <- read_result("pattern_IS_M2_comparison.csv")
prediction <- read_result("pattern_IS_prediction_common_M2.csv")
spec <- read_result("pattern_spec.csv")
config_file <- file.path(DIR_PROJECT, "ukb/00data/analysis_config.rds")
if (!file.exists(config_file)) stop("Missing analysis configuration: ", config_file)
cfg_plot <- readRDS(config_file)
PRIMARY_OUTCOME <- "cvd_stroke_i" # Same primary endpoint as score.R.
# Table 1 uses the exact common IS cohort and saved score groups from s3_2.
# --table1-only writes CSV/QC and, when the bundled Node runtime is available,
# the formatted XLSX without regenerating figures or fitting survival models.
make_dose_table1 <- function() {
    out <- file.path(DIR_PROJECT, "outputs/ukb_table1")
    dir.create(out, recursive=TRUE, showWarnings=FALSE)
    d <- readRDS(file.path(DIR_OUT, "dat.pattern.rds"))
    active <- spec[which(spec$include & spec$available), ]
    covs_all <- unique(c(cfg_plot$covs.m1, cfg_plot$covs.alcohol,
                         cfg_plot$covs.m2.add, cfg_plot$covs.m3.add))
    tv <- paste0(PRIMARY_OUTCOME, ".t2e"); ev <- paste0(PRIMARY_OUTCOME, ".Yt2e")
    keep <- complete.cases(d[c(tv, ev, covs_all, active$variable)]) &
        is.finite(d[[tv]]) & d[[tv]] > 0
    d <- d[keep, ]; g <- droplevels(d$diet.dose.q)
    stopifnot(!anyDuplicated(d$eid), identical(levels(g), paste0("Q",1:4)))
    reference <- cat_res[cat_res$diet=="dose" & cat_res$outcome==PRIMARY_OUTCOME &
                             cat_res$model=="M2", ]
    reference <- reference[match(levels(g), reference$category), ]
    stopifnot(nrow(reference)==4L, all(reference$N==nrow(d)),
              all(reference$Events==sum(d[[ev]])),
              all(as.integer(table(g))==reference$Group_N),
              all(as.integer(tapply(d[[ev]],g,sum))==reference$Group_Events))
    cuts <- read_result("pattern_quartile_definitions.csv")
    cuts <- cuts[cuts$diet=="dose", ];cuts <- cuts[match(levels(g),cuts$category), ]
    intervals <- paste0(ifelse(cuts$lower_inclusive,"[","("),cuts$cut_lower,"-",cuts$cut_upper,"]")
    recut <- cut(d$dose, c(cuts$cut_lower[1],cuts$cut_upper),include.lowest=TRUE,labels=levels(g))
    stopifnot(identical(as.character(recut),as.character(g)))
    covs <- setdiff(unique(c(cfg_plot$covs.m1,cfg_plot$covs.m2.add)),cfg_plot$covs.alcohol)
    labels <- c(dose="DOSE diet score",energy.kJ="Energy intake (kJ/day)",age="Age (years)",
        sex.b="Sex",ethnic.c3="Ethnicity",tdi="Townsend deprivation index",edu.c="Education",
        smoke.c="Smoking status",days_pa_mod.c="Moderate physical activity (days/week)",
        bmi="Body mass index (kg/m2)",fhist_stroke="Family history of stroke",
        primary_ahtn="Hypertension",t2dm="Type 2 diabetes")
    stopifnot(all(covs %in% names(labels)))
    order <- c("dose","energy.kJ","age","sex.b","ethnic.c3","tdi","edu.c",
               "smoke.c","days_pa_mod.c","bmi","fhist_stroke","primary_ahtn","t2dm")
    order <- order[order %in% c("dose",covs)]
    rows <- list(); audit <- list(); numeric_stats <- list()
    add <- function(label,values=rep("",5),p=NA_real_,kind="data",variable="") {
        rows[[length(rows)+1L]] <<- data.frame(Characteristics=label,Total=values[1],
            Q1=values[2],Q2=values[3],Q3=values[4],Q4=values[5],P=p,kind=kind,variable=variable)
    }
    groups <- c(list(Total=seq_len(nrow(d))),lapply(levels(g),function(z)which(g==z)))
    names(groups) <- c("Total",levels(g))
    set.seed(20260922)
    for(v in order) {
        x <- d[[v]];stopifnot(!anyNA(x))
        if(is.numeric(x)) {
            skew <- mean((x-mean(x))^3)/mean((x-mean(x))^2)^1.5
            method <- if(abs(skew)>1) "Kruskal-Wallis" else "One-way ANOVA"
            p <- if(v=="dose") NA_real_ else if(abs(skew)>1)
                kruskal.test(x,g)$p.value else oneway.test(x~g,var.equal=TRUE)$p.value
            vals <- vapply(names(groups),function(lev){
                z<-x[groups[[lev]]];m<-mean(z);s<-sd(z)
                numeric_stats[[length(numeric_stats)+1L]] <<- data.frame(variable=v,group=lev,N=length(z),mean=m,SD=s)
                sprintf("%.2f (%.2f)",m,s)
            },character(1))
            add(paste0(labels[v],", mean (SD)"),vals,p,variable=v)
            audit[[length(audit)+1L]] <- data.frame(variable=v,test=if(v=="dose")"Not tested: grouping exposure" else method,P=p,skewness=skew,min_expected=NA_real_)
        } else {
            x<-droplevels(factor(x));tab<-table(x,g)
            chi<-suppressWarnings(chisq.test(tab,correct=FALSE));mn<-min(chi$expected)
            method<-"Pearson chi-squared"
            if(mn<5){chi<-chisq.test(tab,simulate.p.value=TRUE,B=100000);method<-"Pearson chi-squared (100000 Monte Carlo replicates)"}
            add(paste0(labels[v],", n (%)"),p=chi$p.value,kind="section",variable=v)
            for(lev in levels(x)) {
                lab<-switch(lev,low="Low",middle="Middle",high="High",never="Never",previous="Former",current="Current",lev)
                vals<-vapply(groups,function(ix){n<-sum(x[ix]==lev);sprintf("%s (%.1f%%)",format(n,big.mark=",",trim=TRUE),100*n/length(ix))},character(1))
                add(lab,vals,variable=v)
            }
            stopifnot(all(colSums(tab)==as.integer(table(g))))
            audit[[length(audit)+1L]]<-data.frame(variable=v,test=method,P=chi$p.value,skewness=NA_real_,min_expected=mn)
        }
    }
    write.csv(do.call(rbind,rows),file.path(out,"Table1_data.csv"),row.names=FALSE,na="")
    write.csv(do.call(rbind,audit),file.path(out,"Table1_tests.csv"),row.names=FALSE,na="")
    write.csv(do.call(rbind,numeric_stats),file.path(out,"Table1_continuous_statistics.csv"),row.names=FALSE)
    write.csv(data.frame(group=c("Total",levels(g)),N=c(nrow(d),as.integer(table(g))),
        Events=c(sum(d[[ev]]),reference$Group_Events),range=c("",intervals)),file.path(out,"Table1_cohort_QC.csv"),row.names=FALSE)
    notes <- c(
        "Values are mean (SD) or n (%). The cohort and quartiles exactly match the common complete-case incident IS analysis in score.R. Quartile cutpoints are fixed from the saved source-cohort definitions; tied scores remain in the same group.",
        "Characteristics follow DOSE Model 2: age, sex, ethnicity, education, Townsend deprivation index, energy intake, smoking, moderate physical activity, BMI, family history of stroke, hypertension and type 2 diabetes. DOSE includes wine; alcohol is not separately adjusted. Values use the prepared/imputed covariates used in the fitted models.",
        "P values compare the four quartiles: Pearson chi-squared tests for categorical variables; one-way ANOVA for continuous variables with absolute pooled-sample skewness <=1; Kruskal-Wallis tests when absolute skewness >1. A Monte Carlo chi-squared P value (100,000 replicates) is used if any expected count is <5. P values are unadjusted. DOSE itself is not tested because it defines the groups.",
        "Education is the current low/middle/high categorical variable, not years of education. Physical activity is days of moderate activity per week (UKB Field 884), not MET-minutes/week. Unknown is retained as a category for physical activity and family history of stroke.",
        "Sources: ukb/04score/pattern/dat.pattern.rds, pattern_categories.csv, pattern_quartile_definitions.csv, pattern_spec.csv and ukb/00data/analysis_config.rds. Test details and unrounded statistics are saved alongside this workbook.")
    write.csv(data.frame(note=notes),file.path(out,"Table1_notes.csv"),row.names=FALSE)
    # The optional workbook renderer uses the bundled artifact-tool runtime.
    deps<-file.path(path.expand("~"),".cache/codex-runtimes/codex-primary-runtime/dependencies")
    node<-Sys.getenv("DOSE_TABLE1_NODE",file.path(deps,"node/bin/node.exe"))
    if(!file.exists(node)) node<-Sys.which("node")
    renderer <- file.path(dose_script_dir(), "table1_export.mjs")
    if(nzchar(node) && file.exists(node) && file.exists(renderer)) {
        marker<-file.path(path.expand("~"),".codex/plugins/cache/openai-primary-runtime/spreadsheets/26.905.11957/skills/spreadsheets/container_tools/mark_artifact_operation_started.mjs")
        if(file.exists(marker)) {
            status<-system2(node,c(shQuote(marker),"--operation-kind create --expected-output-count 1 --output-format xlsx"))
            if(status!=0)stop("Table 1 artifact registration failed")
        }
        export_started<-Sys.time()
        status<-system2(node,c(shQuote(file.path(dose_script_dir(), "table1_export.mjs")),shQuote(out)))
        # Some desktop runtime versions return status 1 after a successful export.
        xlsx<-file.path(out,"Table 1.xlsx")
        if(!file.exists(xlsx) || file.info(xlsx)$mtime < export_started)
            stop("Table 1 workbook export failed (status ",status,")")
    } else {
        # Keep workbook export available when the optional custom renderer is absent.
        if (!requireNamespace("openxlsx", quietly=TRUE)) stop("Table 1 export requires openxlsx")
        tables <- c("Table1_data.csv", "Table1_tests.csv", "Table1_continuous_statistics.csv",
                    "Table1_cohort_QC.csv", "Table1_notes.csv")
        sheets <- c("Table 1", "Tests", "Continuous statistics", "Cohort QC", "Notes")
        wb <- openxlsx::createWorkbook()
        for (i in seq_along(tables)) {
            openxlsx::addWorksheet(wb, sheets[i])
            openxlsx::writeData(wb, sheets[i], read.csv(file.path(out,tables[i]), check.names=FALSE),
                headerStyle=openxlsx::createStyle(textDecoration="bold",fgFill="#E7EEF5"))
            openxlsx::freezePane(wb, sheets[i], firstRow=TRUE)
            openxlsx::setColWidths(wb, sheets[i], cols=1:10, widths="auto")
        }
        openxlsx::saveWorkbook(wb, file.path(out,"Table 1.xlsx"), overwrite=TRUE)
    }
    message("Table 1: ",nrow(d)," participants; ",sum(d[[ev]])," IS events. Output: ",out)
}
if (!"--forest-10-only" %in% commandArgs(trailingOnly=TRUE)) make_dose_table1()
if ("--table1-only" %in% commandArgs(trailingOnly=TRUE)) quit(save="no",status=0)
models <- c("M1", "M2", "M3")
scales <- c("P10-P90", "Per 1-SD")
quartiles <- paste0("Q", 1:4)

# Validate current column names and complete panels BEFORE creating any PDF.
require_columns <- function(d, cols, context) {
    missing <- setdiff(cols, names(d))
    if (length(missing)) stop(context, ": missing columns: ", paste(missing, collapse=", "))
}
assert_keys <- function(d, keys, context) {
    require_columns(d, keys, context)
    if (!nrow(d) || anyNA(d[keys]) || anyDuplicated(d[keys])) {
        stop("Missing or duplicated keys: ", context)
    }
}
assert_rows <- function(d, keys, context) {
    assert_keys(d, keys, context)
    require_columns(d, c("HR", "LCI", "HCI"), context)
    values <- as.matrix(d[c("HR", "LCI", "HCI")])
    if (!is.numeric(values) || any(!is.finite(values) | values <= 0) ||
        any(d$LCI > d$HR | d$HR > d$HCI)) stop("Invalid HR/CI: ", context)
}
require_columns(spec, c("diet", "label", "variable", "include", "available", "minimum", "maximum", "step"), "pattern_spec.csv (rerun current s3_2 if outdated)")
use <- spec[which(spec$include & spec$available), ]
assert_keys(use, "diet", "included patterns")
if (!all(c("dose", "dash") %in% use$diet)) stop("DOSE and DASH are required.")
if (any(grepl("^dose", use$diet) & use$diet != "dose"))
    stop("Results still contain alternative DOSE versions. Run the current score.R first.")
primary_spec <- use[use$diet == "dose", ]
if (primary_spec$variable != "dose" ||
    any(!is.finite(unlist(primary_spec[c("minimum", "maximum", "step")]))))
    stop("Invalid primary DOSE metadata in pattern_spec.csv.")
assert_rows(cont, c("diet", "outcome", "model", "scale"), "continuous results")
assert_rows(cat_res, c("diet", "outcome", "model", "category"), "quartile results")
require_columns(cont, c("label", "N", "Events", "P", "FDR"), "continuous results")
require_columns(cat_res, c("N", "Events", "Group_N", "Group_Events", "P", "FDR",
    "score_min", "score_max", "Incidence_percent"), "quartile results")
assert_keys(tests, c("diet", "outcome", "model"), "overall tests")
require_columns(tests, "P_trend", "overall tests")
assert_keys(comparison, "diet", "IS comparison")
require_columns(comparison, "rank_HR_P10P90", "IS comparison")
assert_keys(prediction, "diet", "IS prediction")
require_columns(prediction, c("label", "Delta_C", "Delta_C_LCI", "Delta_C_HCI", "Baseline_C"), "IS prediction")
if (!setequal(comparison$diet, use$diet) || !setequal(prediction$diet, use$diet) ||
    anyNA(comparison$rank_HR_P10P90)) stop("Pattern lists differ between s3_2 exports.")
out_map <- cfg_plot$out_map
Ys <- c(PRIMARY_OUTCOME, setdiff(unique(cont$outcome), PRIMARY_OUTCOME))
if (!length(out_map) || anyNA(out_map[Ys]) || !all(Ys %in% cfg_plot$Ys) || !setequal(Ys, cont$outcome)) {
    stop("Outcome configuration differs from s3_2 results; check input versions.")
}
assert_panel_grid <- function(d, expected, context) {
    keys <- names(expected)
    if (nrow(d) != nrow(expected) ||
        nrow(merge(d[keys], expected, by=keys)) != nrow(expected)) {
        stop("Incomplete or unexpected panels: ", context)
    }
}
assert_panel_grid(cont, expand.grid(diet=use$diet, outcome=Ys, model=models,
    scale=scales, stringsAsFactors=FALSE), "continuous results")
assert_panel_grid(cat_res, expand.grid(diet=use$diet, outcome=Ys, model=models,
    category=quartiles, stringsAsFactors=FALSE), "quartile results")
assert_panel_grid(tests, expand.grid(diet=use$diet, outcome=Ys, model=models,
    stringsAsFactors=FALSE), "overall tests")
assert_rows(rcs, c("outcome", "score"), "DOSE RCS curves")
assert_keys(rcs_tests, "outcome", "DOSE RCS tests")
require_columns(rcs_tests, c("reference", "knot1", "knot2", "knot3", "P_overall", "P_nonlinear"), "DOSE RCS tests")
if (!setequal(rcs$outcome, Ys) || !setequal(rcs_tests$outcome, Ys)) stop("Incomplete RCS outcomes.")
if (!all(rcs$diet == "dose") || !all(rcs_tests$diet == "dose")) stop("RCS files must contain primary DOSE only.")
pattern_order <- comparison$diet[order(comparison$rank_HR_P10P90, comparison$diet)]
diet_labels <- setNames(use$label, use$diet)
cont$label <- unname(diet_labels[cont$diet])
prediction$label <- unname(diet_labels[prediction$diet])
forest_diets <- c("dose", "ahei", "medi24", "mind", "dash", "hlfd",
                  "phdi", "digm", "hpdi", "hlcd")
stopifnot(length(forest_diets) == 10L, all(forest_diets %in% use$diet))
# Preserve the existing IS association ordering, restricted to the requested set.
forest_order <- pattern_order[pattern_order %in% forest_diets]
cont_forest <- cont[cont$diet %in% forest_diets, ]
cont_forest$FDR_all_patterns <- cont_forest$FDR
forest_families <- split(seq_len(nrow(cont_forest)),
    interaction(cont_forest[c("outcome", "model", "scale")], drop=TRUE))
for (ix in forest_families) {
    stopifnot(length(ix) == 10L, setequal(cont_forest$diet[ix], forest_diets),
              all(is.finite(cont_forest$P[ix])),
              all(cont_forest$P[ix] >= 0 & cont_forest$P[ix] <= 1))
    cont_forest$FDR[ix] <- p.adjust(cont_forest$P[ix], method="BH")
}
cont_forest$FDR_family_size <- 10L
write.csv(cont_forest, file.path(DIR_OUT, "pattern_continuous_10patterns.csv"), row.names=FALSE)
cor_input <- file.path(DIR_OUT, "dat.pattern.rds")
if (!file.exists(cor_input)) stop("Missing ", cor_input, ". Run score.R first.")
dir.create(DIR_FIG, recursive=TRUE, showWarnings=FALSE)
C_S <- list(dose="#2F5D8A", sig="#4A7EBB", nonsig="#BFC7D1", ci="#B8C9D9", ref="grey50")
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))
fmt_hr <- function(hr, lo, hi) sprintf("%.2f (%.2f-%.2f)", hr, lo, hi)
point_col <- function(diet, p) ifelse(diet == "dose", C_S$dose,
                                     ifelse(!is.na(p) & p < .05, C_S$sig, C_S$nonsig))
open_pdf <- function(path, width, height) {
    if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf(path, width=width, height=height, family="sans")
    else grDevices::pdf(path, width=width, height=height, family="Helvetica")
}

# 1. Forest report: one page per outcome, with the same diet order throughout.
draw_forest_page <- function(Y) {
    r <- cont_forest[cont_forest$outcome == Y & cont_forest$model == "M2" & cont_forest$scale == "P10-P90", ]
    r <- r[match(forest_order, r$diet), ]
    stopifnot(nrow(r) == 10L, !anyNA(r$HR), length(unique(r$N)) == 1)
    cols <- point_col(r$diet, r$P)
    grid::grid.newpage()
    grid::grid.text(out_map[Y], x=.035, y=.955, just="left",
                    gp=grid::gpar(fontsize=17, fontface="bold", col="#21364B"))
    grid::grid.text(sprintf("Model 2 | %s participants | %s events | P10-P90 increment",
                            format(r$N[1], big.mark=","), format(r$Events[1], big.mark=",")),
                    x=.035, y=.908, just="left", gp=grid::gpar(fontsize=10, col="grey35"))
    grid::grid.text("Dietary pattern", x=.035, y=.851, just="left",
                    gp=grid::gpar(fontsize=10, fontface="bold"))
    grid::grid.text(c("HR (95% CI)", "P", "FDR"),
                    x=c(.738, .89, .97), y=.851,
                    gp=grid::gpar(fontsize=10, fontface="bold"))
    # Center the numeric column headings explicitly.
    # The panel itself is on a log-HR axis; all intervals are kept visible.
    xlo <- min(r$LCI, 1)*.95
    xhi <- max(r$HCI, 1)*1.05
    panel_left <- .27; panel_right <- .59
    xcoord <- function(x) panel_left + (log(x)-log(xlo))/(log(xhi)-log(xlo))*(panel_right-panel_left)
    ys <- seq(.795, .225, length.out=nrow(r))
    row_height <- min(.044, diff(range(ys))/(nrow(r)-1)*.94)
    for (i in seq_len(nrow(r))) {
        if (r$diet[i] == "dose") grid::grid.rect(x=.50, y=ys[i], width=.95, height=row_height,
            gp=grid::gpar(fill="#EDF3F8", col=NA))
    }
    ticks <- pretty(c(xlo, xhi), n=5)
    ticks <- ticks[ticks >= xlo & ticks <= xhi & ticks > 0]
    for (v in ticks) grid::grid.segments(x0=xcoord(v), x1=xcoord(v), y0=.197, y1=.824,
        gp=grid::gpar(col=if (v == 1) "grey50" else "grey88", lty=if (v == 1) 2 else 3, lwd=.8))
    if (!1 %in% ticks) grid::grid.segments(x0=xcoord(1), x1=xcoord(1), y0=.197, y1=.824,
        gp=grid::gpar(col="grey50", lty=2, lwd=.8))
    for (i in seq_len(nrow(r))) {
        face <- if (r$diet[i] == "dose") "bold" else "plain"
        grid::grid.text(r$label[i], x=.035, y=ys[i], just="left", gp=grid::gpar(fontsize=11, fontface=face))
        grid::grid.segments(x0=xcoord(r$LCI[i]), x1=xcoord(r$HCI[i]), y0=ys[i], y1=ys[i],
                            gp=grid::gpar(col=cols[i], lwd=1.8))
        grid::grid.segments(x0=xcoord(c(r$LCI[i], r$HCI[i])), x1=xcoord(c(r$LCI[i], r$HCI[i])),
                            y0=ys[i]-.006, y1=ys[i]+.006, gp=grid::gpar(col=cols[i], lwd=1.4))
        grid::grid.points(x=xcoord(r$HR[i]), y=ys[i], default.units="npc", pch=15, size=grid::unit(2.2, "mm"),
                          gp=grid::gpar(col=cols[i]))
        grid::grid.text(fmt_hr(r$HR[i], r$LCI[i], r$HCI[i]), x=.738, y=ys[i],
                        gp=grid::gpar(fontsize=10, fontface=face))
        grid::grid.text(fmt_p(r$P[i]), x=.89, y=ys[i], gp=grid::gpar(fontsize=9.5))
        grid::grid.text(fmt_p(r$FDR[i]), x=.97, y=ys[i], gp=grid::gpar(fontsize=9.5))
    }
    grid::grid.segments(x0=panel_left, x1=panel_right, y0=.19, y1=.19)
    grid::grid.segments(x0=xcoord(ticks), x1=xcoord(ticks), y0=.19, y1=.182)
    grid::grid.text(format(ticks, trim=TRUE), x=xcoord(ticks), y=.16, gp=grid::gpar(fontsize=9))
    grid::grid.text("Hazard ratio", x=(panel_left+panel_right)/2, y=.115, gp=grid::gpar(fontsize=10))
    grid::grid.text(sprintf("Deep blue: DOSE. Blue: P < 0.05. Grey: P >= 0.05. FDR: BH adjustment across %d patterns.", length(forest_order)),
                    x=.035, y=.057, just="left", gp=grid::gpar(fontsize=8, col="grey35"))
    grid::grid.text("Diet order follows the ischemic-stroke M2 association; HR rankings do not establish predictive superiority.",
                    x=.035, y=.025, just="left", gp=grid::gpar(fontsize=8, col="grey35"))
}
forest_file <- file.path(DIR_FIG, "Fig_Pattern_Stroke_Forest_M2.pdf")
open_pdf(forest_file, 8.4, 6.6)
tryCatch(for (Y in Ys) draw_forest_page(Y), finally=bounded_dev_off())
if ("--forest-10-only" %in% commandArgs(trailingOnly=TRUE)) {
    cat("Saved ten-pattern forest:", forest_file, "\n")
    quit(save="no", status=0L)
}

# 2. DOSE quartile HRs and RCS: two outcomes per page, M2 throughout.
draw_quartile <- function(Y) {
    r <- cat_res[cat_res$diet == "dose" & cat_res$model == "M2" & cat_res$outcome == Y, ]
    r <- r[order(r$category), ]
    tr <- tests[tests$diet == "dose" & tests$model == "M2" & tests$outcome == Y, ]
    labels <- ifelse(r$score_min == r$score_max, as.character(r$score_min), paste0(r$score_min, "-", r$score_max))
    span <- diff(range(c(r$LCI, r$HCI, 1)))
    ylim <- c(max(.01, min(r$LCI)-.12*span), max(r$HCI)+.30*span)
    plot(1:4, r$HR, type="n", xlim=c(.7, 4.3), ylim=ylim, xaxt="n",
         xlab="DOSE score", ylab="HR (95% CI)", main=paste0(out_map[Y], ": quartiles"))
    axis(1, at=1:4, labels=labels)
    abline(h=1, lty=2, col="grey60")
    segments(1:4, r$LCI, 1:4, r$HCI, col=C_S$dose, lwd=1.5)
    segments((1:4)-.06, r$LCI, (1:4)+.06, r$LCI, col=C_S$dose, lwd=1.5)
    segments((1:4)-.06, r$HCI, (1:4)+.06, r$HCI, col=C_S$dose, lwd=1.5)
    lines(1:4, r$HR, col=C_S$dose, lwd=1.5); points(1:4, r$HR, pch=16, col=C_S$dose, cex=1.05)
    text(1:4, r$HCI+.07*span, labels=c("Ref.", sprintf("%.2f", r$HR[-1])), cex=.85)
    mtext(paste0("P-trend: ", fmt_p(tr$P_trend)),
          side=3, line=.1, adj=1, cex=.75)
}
draw_rcs <- function(Y) {
    r <- rcs[rcs$outcome == Y, ]; r <- r[order(r$score), ]
    tt <- rcs_tests[rcs_tests$outcome == Y, ]
    span <- diff(range(c(r$LCI, r$HCI)))
    plot(r$score, r$HR, type="n", xlim=c(primary_spec$minimum, primary_spec$maximum), ylim=c(max(.01, min(r$LCI)-.06*span), max(r$HCI)+.20*span),
         xlab="DOSE score", ylab="HR (95% CI)", main=paste0(out_map[Y], ": RCS"), xaxt="n")
    axis(1, at=pretty(c(primary_spec$minimum, primary_spec$maximum), n=9))
    polygon(c(r$score, rev(r$score)), c(r$LCI, rev(r$HCI)), col=adjustcolor(C_S$ci, alpha.f=.6), border=NA)
    abline(h=1, lty=2, col="grey60"); abline(v=tt$reference, lty=3, col="grey70")
    lines(r$score, r$HR, col=C_S$dose, lwd=2)
    legend("topright", legend=c(paste0("P-overall: ", fmt_p(tt$P_overall)),
                                paste0("P-nonlinear: ", fmt_p(tt$P_nonlinear))), bty="n", cex=.8)
    mtext(paste0("Reference = ", tt$reference, "; knots = ", paste(unlist(tt[c("knot1", "knot2", "knot3")]), collapse=", ")),
          side=3, line=.1, adj=0, cex=.72)
}
rcs_file <- file.path(DIR_FIG, "Fig_DOSE_Quartiles_RCS_M2.pdf")
open_pdf(rcs_file, 10.8, if (length(Ys) == 1) 5 else 8.8)
tryCatch({
    par(mfrow=c(if (length(Ys) == 1) 1 else 2, 2), mar=c(4.2, 4.3, 3.2, 1), oma=c(1, 0, 1.6, 0),
        mgp=c(2.5, .8, 0), cex=.9, las=1, bty="l")
    for (i in seq_along(Ys)) {
        draw_quartile(Ys[i]); draw_rcs(Ys[i])
        if (i %% 2 == 0 || i == length(Ys)) {
            mtext("DOSE and incident stroke | Model 2", outer=TRUE, side=3, line=.2, adj=.03, cex=1.15, font=2)
            mtext("Group cutoffs follow quantile() + cut(right=TRUE); spline bands are pointwise 95% confidence intervals.",
                  outer=TRUE, side=1, line=0, cex=.75)
        }
    }
}, finally=bounded_dev_off())

# 3. Association strength and apparent discrimination are different targets.
comparison_file <- file.path(DIR_FIG, "Fig_Pattern_IS_Comparison_M2.pdf")
open_pdf(comparison_file, 12.8, 6.5)
tryCatch({
    a <- cont[cont$outcome == "cvd_stroke_i" & cont$model == "M2" & cont$scale == "Per 1-SD", ]
    a <- a[order(a$HR), ]; b <- prediction[match(a$diet, prediction$diet), ]
    yy <- rev(seq_len(nrow(a))); cols <- point_col(a$diet, a$P)
    par(mfrow=c(1, 2), mar=c(4.4, 7.5, 3, 1.3), oma=c(3, 0, 2, 0), mgp=c(2.6, .8, 0), las=1, cex=.85, bty="l")
    plot(a$HR, yy, type="n", log="x", xlim=range(c(a$LCI*.98, a$HCI*1.02, 1)),
         ylim=c(.4, nrow(a)+.6), yaxt="n", xlab="HR per 1-SD increment (95% CI)", ylab="", main="Association strength")
    axis(2, at=yy, labels=a$label, tick=FALSE, cex.axis=.9)
    abline(v=1, lty=2, col="grey60")
    segments(a$LCI, yy, a$HCI, yy, col=cols, lwd=1.7); points(a$HR, yy, pch=15, col=cols, cex=1)
    limits <- range(c(b$Delta_C_LCI, b$Delta_C_HCI, 0)); pad <- max(diff(limits)*.12, .0001)
    plot(b$Delta_C, yy, type="n", xlim=limits+c(-pad, pad), ylim=c(.4, nrow(a)+.6), yaxt="n",
         xlab="Change in Harrell C versus covariates alone", ylab="", main="Apparent discrimination")
    axis(2, at=yy, labels=b$label, tick=FALSE, cex.axis=.9)
    abline(v=0, lty=2, col="grey60")
    segments(b$Delta_C_LCI, yy, b$Delta_C_HCI, yy, col=cols, lwd=1.7)
    points(b$Delta_C, yy, pch=15, col=cols, cex=1)
    mtext(sprintf("Common M2 covariates | Baseline C = %.4f", b$Baseline_C[1]), side=3, line=.1, cex=.78)
    mtext("Ischemic stroke: dietary-pattern comparison", outer=TRUE, side=3, line=.3, adj=.03, font=2, cex=1.2)
    mtext("Association: pattern-specific alcohol adjustment. Discrimination: the same DOSE-model covariates for every pattern.",
          outer=TRUE, side=1, line=0, cex=.74)
    mtext("Discrimination is evaluated in the development sample; DOSE's food selection and cutoffs have not been independently validated here.",
          outer=TRUE, side=1, line=1.1, cex=.74)
}, finally=bounded_dev_off())

# Publication-style tables, with all model results available for inspection.
tc <- cont[cont$scale == "P10-P90", ]
tc$Outcome <- unname(out_map[tc$outcome]); tc$Pattern <- tc$label
tc$HR_CI <- fmt_hr(tc$HR, tc$LCI, tc$HCI); tc$P_value <- fmt_p(tc$P); tc$FDR_value <- fmt_p(tc$FDR)
wide_cont <- reshape(tc[c("Outcome", "Pattern", "N", "Events", "model", "HR_CI", "P_value", "FDR_value")],
    idvar=c("Outcome", "Pattern", "N", "Events"), timevar="model", direction="wide")
write.csv(wide_cont, file.path(DIR_OUT, "Table_Pattern_Stroke_M1_M3.csv"), row.names=FALSE)
tq <- cat_res[cat_res$diet == "dose", ]
tq$Outcome <- unname(out_map[tq$outcome])
tq$Score <- ifelse(tq$score_min == tq$score_max, as.character(tq$score_min), paste0(tq$score_min, "-", tq$score_max))
tq$HR_CI <- ifelse(tq$category == "Q1", "Ref.", fmt_hr(tq$HR, tq$LCI, tq$HCI))
tq$P_value <- fmt_p(tq$P)
keys <- c("outcome", "model")
tr <- tests[tests$diet == "dose", c(keys, "P_trend")]
tq <- merge(tq, tr, by=keys, all.x=TRUE, sort=FALSE)
tq$P_trend <- fmt_p(tq$P_trend)
wide_cat <- reshape(tq[c("Outcome", "category", "Score", "Group_N", "Group_Events", "Incidence_percent", "model", "HR_CI", "P_value", "P_trend")],
    idvar=c("Outcome", "category", "Score", "Group_N", "Group_Events", "Incidence_percent"), timevar="model", direction="wide")
wide_cat <- wide_cat[order(match(wide_cat$Outcome, out_map), wide_cat$category), ]
write.csv(wide_cat, file.path(DIR_OUT, "Table_DOSE_Quartiles_M1_M3.csv"), row.names=FALSE)

# 4. Legacy multi-pattern figures, adapted to the current result schema.
# Use base/grid so this script remains runnable without package installation.
diet_colors <- c(medi24="#6E8FB2", dash="#7DA494", mind="#D8A0C1",
    hpdi="#EAB67A", ahei="#E5A79A", phdi="#C16E71", digm="#ABC8E5",
    modern="#B8B8B8", hlcd="#C9A227", hlfd="#4DAF4A", redii="#D95F02",
    fds="#7570B3", medito="#5C798E", maha="#9673A6", dose=C_S$dose)
stopifnot(all(pattern_order %in% names(diet_colors)))
sig_marks <- function(p, fdr) paste0(ifelse(!is.na(p) & p < .05, "*", ""),
                                    ifelse(!is.na(fdr) & fdr < .05, "\u2020", ""))
event_label <- function(events, n) paste0(format(events, big.mark=",", trim=TRUE), "/",
                                         format(n, big.mark=",", trim=TRUE))

# Reusable aligned forest table. Each panel has a true log-HR axis; limits
# include every CI, unlike the fixed/clipped axes in the original script.
draw_multi_forest <- function(panels, row_labels, title, subtitle, footer) {
    n <- length(row_labels); k <- length(panels)
    stopifnot(all(vapply(panels, nrow, integer(1)) == n))
    grid::grid.newpage()
    grid::grid.text(title, x=.025, y=.96, just="left",
                    gp=grid::gpar(fontsize=17, fontface="bold", col="#21364B"))
    grid::grid.text(subtitle, x=.025, y=.918, just="left", gp=grid::gpar(fontsize=10))
    ys <- seq(.78, .22, length.out=n)
    row_h <- min(.09, .54/max(1, n-1))
    for (i in seq_len(n)) if (i %% 2 == 0) grid::grid.rect(x=.5, y=ys[i], width=.96,
        height=row_h, gp=grid::gpar(fill="#F2F5F8", col=NA))
    grid::grid.text(row_labels, x=.025, y=ys, just="left", gp=grid::gpar(fontsize=10))
    limits <- range(c(1, unlist(lapply(panels, function(r) c(r$LCI, r$HCI)))))
    limits <- exp(log(limits) + c(-.08, .08))
    ticks <- pretty(limits, n=4); ticks <- ticks[ticks > 0 & ticks >= limits[1] & ticks <= limits[2]]
    panel_w <- .82/k
    for (j in seq_len(k)) {
        r <- panels[[j]]; left <- .17+(j-1)*panel_w
        x0 <- left+.26*panel_w; x1 <- left+.57*panel_w
        xc <- function(x) x0 + log(x/limits[1])/log(limits[2]/limits[1])*(x1-x0)
        grid::grid.text(names(panels)[j], x=left+panel_w/2, y=.875,
                        gp=grid::gpar(fontsize=12, fontface="bold"))
        grid::grid.text(c("Events/N", "HR", "HR (95% CI)", "P"),
            x=left+panel_w*c(.115, .415, .77, .98), y=.829, gp=grid::gpar(fontsize=8.5))
        grid::grid.segments(x0=xc(1), x1=xc(1), y0=.18, y1=.805,
                            gp=grid::gpar(col="grey55", lty=2))
        grid::grid.text(event_label(r$plot_events, r$plot_n), x=left+.115*panel_w,
                        y=ys, gp=grid::gpar(fontsize=8))
        cols <- unname(diet_colors[r$diet])
        grid::grid.segments(x0=xc(r$LCI), x1=xc(r$HCI), y0=ys, y1=ys,
                            gp=grid::gpar(col=cols, lwd=1.6))
        grid::grid.points(x=xc(r$HR), y=ys, default.units="npc", pch=15, size=grid::unit(2, "mm"),
                          gp=grid::gpar(col=cols))
        est <- ifelse(r$is_ref, "Ref.", paste0(fmt_hr(r$HR, r$LCI, r$HCI), sig_marks(r$P, r$FDR)))
        grid::grid.text(est, x=left+.77*panel_w, y=ys, gp=grid::gpar(fontsize=8))
        grid::grid.text(fmt_p(r$P), x=left+.98*panel_w, y=ys, gp=grid::gpar(fontsize=8))
        grid::grid.segments(x0=x0, x1=x1, y0=.175, y1=.175)
        grid::grid.text(format(ticks, trim=TRUE), x=xc(ticks), y=.145, gp=grid::gpar(fontsize=8))
        grid::grid.segments(x0=xc(ticks), x1=xc(ticks), y0=.175, y1=.167)
    }
    grid::grid.text(footer, x=.025, y=.065, just="left", gp=grid::gpar(fontsize=9))
    grid::grid.text("* P < 0.05; \u2020 BH FDR < 0.05 (from the s3 analysis tables).", x=.025, y=.027,
                    just="left", gp=grid::gpar(fontsize=9, col="grey35"))
}
as_plot_rows <- function(r, categorical=FALSE) {
    r$plot_n <- if (categorical) r$Group_N else r$N
    r$plot_events <- if (categorical) r$Group_Events else r$Events
    r$is_ref <- if (categorical) r$category == "Q1" else FALSE
    r
}
multi_file <- file.path(DIR_FIG, "Fig_Pattern_Stroke_Models_M1_M3.pdf")
open_pdf(multi_file, 19, 7.5)
tryCatch({
    for (scale in c("P10-P90", "Per 1-SD")) for (Y in Ys) {
        panels <- setNames(lapply(models, function(m) {
            r <- cont[cont$outcome == Y & cont$model == m & cont$scale == scale, ]
            if (!setequal(r$diet, pattern_order)) stop("Incomplete model panel: ", Y, " ", m)
            as_plot_rows(r[match(pattern_order, r$diet), ])
        }), models)
        draw_multi_forest(panels, unname(diet_labels[pattern_order]), out_map[Y],
            paste0("Dietary patterns | ", scale, " increment | Models 1-3"),
            "Diet order follows the ischemic-stroke M2 P10-P90 association. Alcohol adjustment follows each pattern's specification.")
    }
}, finally=bounded_dev_off())

# One PDF per pattern; actual Q1 counts and fitted Q2-Q4 effects are used.
quartile_files <- character(0)
for (diet in pattern_order) {
    f <- file.path(DIR_FIG, paste0("Fig_Quartiles_", diet, "_M1_M3.pdf"))
    open_pdf(f, 19, 6)
    tryCatch({
        for (Y in Ys) {
            panels <- setNames(lapply(models, function(m) {
                r <- cat_res[cat_res$diet == diet & cat_res$outcome == Y & cat_res$model == m, ]
                if (!setequal(r$category, paste0("Q", 1:4))) stop("Incomplete quartile panel")
                as_plot_rows(r[match(paste0("Q", 1:4), r$category), ], TRUE)
            }), models)
            tt <- tests[tests$diet == diet & tests$outcome == Y, ]
            trend <- paste(models, fmt_p(tt$P_trend[match(models, tt$model)]), sep=": ", collapse=" | ")
            bounds <- panels[[1]]
            score_labels <- ifelse(bounds$score_min == bounds$score_max,
                format(bounds$score_min, trim=TRUE, digits=4),
                paste0(format(bounds$score_min, trim=TRUE, digits=4), " - ",
                       format(bounds$score_max, trim=TRUE, digits=4)))
            draw_multi_forest(panels, paste0(quartiles, " (", score_labels, ")"), paste(diet_labels[diet], "|", out_map[Y]),
                "Quartiles Q2-Q4 versus Q1 | Models 1-3", paste0("P-trend: ", trend,
                ". Identical scores stay together; group sizes need not be equal."))
        }
    }, finally=bounded_dev_off())
    quartile_files <- c(quartile_files, f)
}

# DOSE-DASH: consolidate the legacy quintile and four-group versions into
# current quartiles, with BOTH continuous scales explicitly identified.
dash_file <- file.path(DIR_FIG, "Fig_DOSE_DASH_Quartiles_M2.pdf")
open_pdf(dash_file, 14.5, 6.5)
tryCatch({
    for (Y in Ys) {
        panels <- setNames(lapply(c("dose", "dash"), function(diet) {
            q <- cat_res[cat_res$diet == diet & cat_res$outcome == Y & cat_res$model == "M2", ]
            q <- as_plot_rows(q[match(paste0("Q", 1:4), q$category), ], TRUE)
            s <- cont[cont$diet == diet & cont$outcome == Y & cont$model == "M2", ]
            s <- as_plot_rows(s[match(c("P10-P90", "Per 1-SD"), s$scale), ])
            common <- intersect(names(q), names(s))
            rbind(q[common], s[common])
        }), c("DOSE", "DASH"))
        draw_multi_forest(panels, c(paste0("Q", 1:4), "P10-P90", "Per 1-SD"),
            paste("DOSE and DASH |", out_map[Y]), "Model 2 | Quartiles and continuous associations",
            "Quartile rows show group events/N; continuous rows show the full analysis cohort. Q1 is the reference.")
    }
}, finally=bounded_dev_off())

# Consistency scatter: label numbers and a keyed legend prevent collisions
# for the many nearly identical HRs, without requiring ggrepel.
scatter_file <- file.path(DIR_FIG, "Fig_Consistency_IS_ICH_M2.pdf")
scatter_data <- NULL
if (all(c("cvd_stroke_i", "cvd_stroke_ih") %in% Ys)) {
    open_pdf(scatter_file, 10, 7)
    tryCatch({
        scatter_tables <- list()
        for (scale in c("P10-P90", "Per 1-SD")) {
            a <- cont[cont$outcome == "cvd_stroke_i" & cont$model == "M2" & cont$scale == scale, ]
            a <- a[match(pattern_order, a$diet), ]
            b <- cont[cont$outcome == "cvd_stroke_ih" & cont$model == "M2" & cont$scale == scale, ]
            b <- b[match(a$diet, b$diet), ]
            lim <- range(c(a$HR, b$HR, 1)); lim <- lim+c(-1, 1)*max(.06, diff(lim)*.15)
            par(mar=c(5, 5, 4, 12), xpd=FALSE, las=1, bty="l")
            plot(a$HR, b$HR, type="n", xlim=lim, ylim=lim, asp=1,
                 xlab="Ischemic stroke HR", ylab="Intracerebral hemorrhage HR",
                 main=paste0("Stroke-type consistency | M2 | ", scale))
            abline(h=1, v=1, lty=2, col="grey65"); abline(a=0, b=1, lty=3, col=C_S$dose)
            points(a$HR, b$HR, pch=16, col=diet_colors[a$diet], cex=1.2)
            # Place identifiers on a spaced rail and connect to the points.
            label_y <- seq(lim[1], lim[2], length.out=nrow(a))
            ord <- order(b$HR); xx <- par("usr")[2]
            segments(a$HR[ord], b$HR[ord], xx, label_y, col=adjustcolor(diet_colors[a$diet[ord]], .6))
            text(xx, label_y, labels=ord, pos=4, xpd=NA, cex=.75)
            legend("right", inset=c(-.56, 0), legend=paste(seq_len(nrow(a)), a$label),
                   col=diet_colors[a$diet], pch=16, bty="n", xpd=NA, cex=.85)
            mtext("Above y=x: lower HR for ischemic stroke. Descriptive comparison; no heterogeneity test.",
                  side=1, line=3.8, cex=.77)
            scatter_tables[[scale]] <- data.frame(diet=a$diet, label=a$label, scale=scale,
                HR_IS=a$HR, HR_ICH=b$HR)
        }
        scatter_data <- do.call(rbind, scatter_tables)
    }, finally=bounded_dev_off())
    write.csv(scatter_data, file.path(DIR_OUT, "Plot_Consistency_IS_ICH_M2.csv"), row.names=FALSE)
}

trend_file <- file.path(DIR_FIG, "Fig_Pattern_Quartile_Trends_M2.pdf")
open_pdf(trend_file, 12, if (length(Ys) == 1) 6.5 else 8.5)
tryCatch({
    ncols <- if (length(Ys) == 1) 1L else 2L
    nrows <- ceiling(length(Ys)/ncols)
    panels <- matrix(seq_len(nrows*ncols), nrow=nrows, byrow=TRUE)
    layout(rbind(panels, rep(nrows*ncols+1L, ncols)), heights=c(rep(1,nrows),.28))
    par(mar=c(4, 4, 3, 1), oma=c(1, 0, 0, 0), las=1, bty="l")
    for (Y in Ys) {
        rr <- cat_res[cat_res$outcome == Y & cat_res$model == "M2", ]
        plot(1:4, rep(1, 4), type="n", ylim=range(c(1, rr$HR)), xaxt="n",
             xlab="Quartile", ylab="Hazard ratio (Model 2)", main=out_map[Y])
        axis(1, at=1:4, labels=paste0("Q", 1:4)); abline(h=1, lty=2, col="grey60")
        for (diet in pattern_order) {
            r <- rr[rr$diet == diet, ]; r <- r[match(paste0("Q", 1:4), r$category), ]
            lines(1:4, r$HR, type="b", pch=16, cex=.65, col=diet_colors[diet],
                  lwd=if (diet == "dose") 2.5 else 1)
        }
    }
    if (length(Ys) < nrows*ncols) for (i in seq_len(nrows*ncols-length(Ys))) plot.new()
    par(mar=rep(0, 4))
    plot.new()
    legend("center", legend=unname(diet_labels[pattern_order]), col=diet_colors[pattern_order],
           lty=1, pch=16, ncol=5, bty="n", cex=.8)
    mtext("Points are HR estimates; confidence intervals are provided in the per-pattern forest PDFs.",
          side=1, line=0, outer=TRUE, cex=.75)
}, finally=bounded_dev_off())

# Score correlation and hierarchical clustering. Explicitly load the cohort;
# the legacy script relied on an unspecified global `dat` object.
dat_cor <- as.data.frame(readRDS(cor_input))
score_vars <- use$variable[match(pattern_order, use$diet)]
stopifnot(all(score_vars %in% names(dat_cor)))
dat_cor <- dat_cor[score_vars]
if (any(!vapply(dat_cor, is.numeric, logical(1))) ||
    any(vapply(dat_cor, function(x) any(is.infinite(x)), logical(1)))) {
    stop("Score correlation requires numeric, non-infinite raw scores.")
}
names(dat_cor) <- unname(diet_labels[pattern_order])
M <- cor(dat_cor, use="pairwise.complete.obs", method="spearman")
if (any(!is.finite(M))) stop("Undefined score correlations: inspect missingness/constant scores.")
n <- ncol(M); P <- matrix(NA_real_, n, n, dimnames=dimnames(M)); Npair <- P
diag(Npair) <- vapply(dat_cor, function(x) sum(is.finite(x)), integer(1))
for (i in seq_len(n-1L)) for (j in seq.int(i+1L, n)) {
    ok <- is.finite(dat_cor[[i]]) & is.finite(dat_cor[[j]])
    Npair[i, j] <- Npair[j, i] <- sum(ok)
    if (sum(ok) >= 3L && length(unique(dat_cor[[i]][ok])) > 1L &&
        length(unique(dat_cor[[j]][ok])) > 1L) {
        P[i, j] <- P[j, i] <- suppressWarnings(cor.test(dat_cor[[i]][ok], dat_cor[[j]][ok],
            method="spearman", exact=FALSE)$p.value)
    }
}
F <- P; F[upper.tri(F)] <- p.adjust(P[upper.tri(P)], "BH")
F[lower.tri(F)] <- t(F)[lower.tri(F)]
for (item in c("M", "P", "F", "Npair")) write.csv(get(item),
    file.path(DIR_OUT, paste0("Pattern_Spearman_", item, ".csv")), row.names=TRUE)
# Match the original pheatmap default: Euclidean distances between rows of
# the correlation matrix, complete linkage, identical row/column ordering.
hc <- hclust(dist(M), method="complete"); ord <- hc$order
write.csv(data.frame(order=seq_len(n), pattern=rownames(M)[ord]),
          file.path(DIR_OUT, "Pattern_Correlation_Cluster_Order.csv"), row.names=FALSE)
cor_file <- file.path(DIR_FIG, "Fig_Pattern_Correlation_Cluster.pdf")
open_pdf(cor_file, 10.5, 10)
tryCatch({
    grid::grid.newpage()
    grid::grid.text("Dietary score correlations", x=.5, y=.97,
                    gp=grid::gpar(fontsize=17, fontface="bold"))
    grid::grid.text(sprintf("Spearman correlation | %s participants | complete-linkage clustering", format(nrow(dat_cor), big.mark=",")),
                    x=.5, y=.935, gp=grid::gpar(fontsize=10))
    xleft <- .17; ybottom <- .18; width <- .70; cell <- width/n
    xc <- xleft + (seq_len(n)-.5)*cell
    yc <- ybottom + width - (seq_len(n)-.5)*cell
    # Top dendrogram shares leaf coordinates with the heatmap columns.
    node_x <- node_y <- numeric(n-1)
    leaf_pos <- match(seq_len(n), ord)
    for (i in seq_len(n-1)) {
        kids <- hc$merge[i, ]
        xx <- vapply(kids, function(k) if (k < 0) xc[leaf_pos[-k]] else node_x[k], numeric(1))
        yy <- vapply(kids, function(k) if (k < 0) .885 else node_y[k], numeric(1))
        h <- .885+.035*hc$height[i]/max(max(hc$height), .Machine$double.eps)
        grid::grid.lines(x=c(xx[1], xx[1], xx[2], xx[2]), y=c(yy[1], h, h, yy[2]),
                         gp=grid::gpar(col="grey35", lwd=.7))
        node_x[i] <- mean(xx); node_y[i] <- h
    }
    pal <- colorRampPalette(c("#C16E71", "white", "#2F5D8A"))(201)
    for (i in seq_len(n)) for (j in seq_len(n)) {
        val <- M[ord[i], ord[j]]; p <- P[ord[i], ord[j]]
        stars <- if (i == j || is.na(p)) "" else if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else ""
        grid::grid.rect(x=xc[j], y=yc[i], width=cell, height=cell,
            gp=grid::gpar(fill=pal[round((val+1)*100)+1], col="grey90", lwd=.5))
        grid::grid.text(paste0(sprintf("%.2f", val), if (nzchar(stars)) paste0("\n", stars) else ""),
            x=xc[j], y=yc[i], gp=grid::gpar(fontsize=7.5, col=if (abs(val) > .65) "white" else "black"))
    }
    grid::grid.text(rownames(M)[ord], x=xleft-.012, y=yc, just="right", gp=grid::gpar(fontsize=9))
    for (j in seq_len(n)) grid::grid.text(colnames(M)[ord[j]], x=xc[j], y=ybottom-.018,
        rot=45, just="right", gp=grid::gpar(fontsize=9))
    for (i in seq_len(201)) grid::grid.rect(x=.915, y=ybottom+(i-.5)/201*width,
        width=.02, height=width/201+.0002, gp=grid::gpar(fill=pal[i], col=NA))
    grid::grid.text(c("-1", "-0.5", "0", "0.5", "1"), x=.94,
        y=ybottom+c(0, .25, .5, .75, 1)*width, just="left", gp=grid::gpar(fontsize=9))
    grid::grid.text("* P < 0.05; ** P < 0.01; *** P < 0.001 (unadjusted, two-sided).", x=.03, y=.065,
                    just="left", gp=grid::gpar(fontsize=9))
    grid::grid.text("Pairwise N, raw P, and BH-adjusted P across unique score pairs are exported separately.",
                    x=.03, y=.035, just="left", gp=grid::gpar(fontsize=9))
}, finally=bounded_dev_off())
rm(dat_cor)

# Appendix: the ten manuscript patterns, using the same score correlations.
cor_appendix_files <- local({
    library(grid)
    input <- DIR_OUT
    output <- DIR_FIG
    stem <- "Fig_Pattern_Correlations_10patterns"
    full <- M
    counts <- Npair
    selected <- unique(cont_forest$label)
    stopifnot(length(selected) == 10L, selected[1] == "DOSE",
        all(selected %in% rownames(full)),
        identical(dimnames(full), dimnames(counts)),
        isTRUE(all.equal(full, t(full))))
    m <- full[selected, selected]
    npair <- counts[selected, selected]
    stopifnot(all(is.finite(m)), all(abs(diag(m) - 1) < 1e-12), all(abs(m) <= 1 + 1e-12))
    m[] <- pmax(-1, pmin(1, m)); diag(m) <- 1
    others <- selected[-1]
    nn <- npair[lower.tri(npair)]
    fmt <- function(x) format(x, big.mark=",", scientific=FALSE, trim=TRUE)
    sample <- if (length(unique(nn)) == 1L) paste0("N = ", fmt(nn[1])) else
        paste0("Pairwise N = ", fmt(min(nn)), "-", fmt(max(nn)))
    red <- colorRampPalette(c("#fff4f0", "#ef8b7a", "#bd1728"))(1001)
    orange <- colorRampPalette(c("#fee6ce", "#f9a15b", "#d65b12"))(1001)
    negative <- colorRampPalette(c("#f0f7ff", "#3182bd", "#08519c"))(1001)
    color <- function(pal, rho) ifelse(rho < 0, negative[round(abs(rho) * 1000) + 1L], pal[round(abs(rho) * 1000) + 1L])
    # Use Unicode rho instead of plotmath so the symbol also uses Arial.
    rho_symbol <- "\u03c1"
    draw <- function() {
        grid.newpage()
        pushViewport(viewport(gp=gpar(fontfamily="Arial", col="black")))
        grid.text("Correlations between different dietary patterns", .5, .961,
            just="top", gp=gpar(fontsize=16, fontface="bold"))
        grid.text(paste("Spearman correlation |", sample), .5, .925,
            gp=gpar(fontsize=10, col="#555555"))
        pushViewport(viewport(x=.525, y=.565, width=unit(7.9,"inches"),
            height=unit(7.9 * 9.9 / 11.05,"inches"),
            xscale=c(-1.45,9.6), yscale=c(9,-.9), clip="off"))
        u <- function(x) unit(x, "native")
        txt <- function(label, x, y, size=10, ...) grid.text(label, u(x), u(y),
            gp=gpar(fontsize=size), ...)
        for (i in 1:8) for (j in 0:(i-1)) {
            rho <- m[i+2,j+2]
            grid.rect(u(j),u(i),u(1),u(1),
                gp=gpar(fill="white",col="#b5b5b5",lwd=.65))
            side <- .94 * sqrt(abs(rho))
            grid.rect(u(j),u(i),u(side),u(side),
                gp=gpar(fill=color(red,rho),col="#9b7773",lwd=.6))
            txt(sprintf("%.2f",rho),j,i)
        }
        hub <- c(7.15,1.05)
        # Sample quadratic Bezier curves in score coordinates.
        t <- seq(0,1,length.out=250)
        for (i in 0:8) {
            rho <- m[1,i+2]
            control <- if (i < hub[1]) c(hub[1],i) else c(i,hub[2])
            xx <- (1-t)^2*i + 2*(1-t)*t*control[1] + t^2*hub[1]
            yy <- (1-t)^2*i + 2*(1-t)*t*control[2] + t^2*hub[2]
            grid.lines(u(xx),u(yy),gp=gpar(col=color(orange,rho),lwd=5.6,lineend="round"))
        }
        for (i in 0:8) {
            grid.points(u(i),u(i),pch=21,size=unit(2.1,"mm"),
                gp=gpar(fill="#d65b12",col="white",lwd=.7))
            grid.rect(u(i),u(i-.34),width=unit(.34,"inches"),height=unit(.17,"inches"),
                gp=gpar(fill="white",col=NA))
            txt(sprintf("%.2f",m[1,i+2]),i,i-.34)
            txt(others[i+1],-.67,i,11,just="right")
            txt(others[i+1],i,8.68,11,rot=90,just="right")
        }
        grid.points(u(hub[1]),u(hub[2]),pch=21,size=unit(3,"mm"),
            gp=gpar(fill="#d65b12",col="white",lwd=.9))
        grid.text("DOSE",u(hub[1]+.46),u(hub[2]),just="left",
            gp=gpar(fontsize=13,fontface="bold"))
        popViewport()
        legend <- function(left, pal, title) {
            steps <- 400L
            grid.rect(x=left+(seq_len(steps)-.5)/steps*.28,y=.102,
                width=.28/steps+.0001,height=.020,
                gp=gpar(fill=color(pal,(seq_len(steps)-.5)/steps),col=NA))
            ticks <- seq(0,1,.2)
            grid.segments(left+ticks*.28,.092,left+ticks*.28,.088,gp=gpar(lwd=.6))
            grid.text(sprintf("%.1f",ticks),left+ticks*.28,.078,gp=gpar(fontsize=9))
            grid.text(paste0(title," (",rho_symbol,")"),left+.14,.134,gp=gpar(fontsize=10))
        }
        legend(.15,red,"Between other dietary patterns")
        legend(.59,orange,"With DOSE")
        grid.text("Area shows absolute correlation; blue indicates negative correlation; labels retain the sign.",
            .5,.029,gp=gpar(fontsize=9,col="#555555"))
        popViewport()
    }
    stopifnot(capabilities("cairo"))
    dir.create(output,recursive=TRUE,showWarnings=FALSE)
    # Cairo embeds the installed Arial font and keeps PDF text selectable.
    cairo_pdf(file.path(output,paste0(stem,".pdf")),width=9.3,height=9.4,family="Arial")
    tryCatch(draw(),finally=bounded_dev_off())
    bounded_png(file.path(output,paste0(stem,".png")),width=9.3,height=9.4,
        units="in",res=400,type="cairo",family="Arial")
    tryCatch(draw(),finally=bounded_dev_off())
    pairs <- do.call(rbind,lapply(1:9,function(i) do.call(rbind,lapply((i+1):10,
        function(j) data.frame(pattern_1=selected[i],pattern_2=selected[j],
            spearman_rho=m[i,j],N=npair[i,j])))))
    stopifnot(nrow(pairs) == 45L)
    write.csv(pairs,file.path(input,"Pattern_Spearman_10patterns_pairs.csv"),row.names=FALSE)
    message("Validated 10 patterns and 45 unique pairs; ",sample,". Font: Arial.")
    message(file.path(output,paste0(stem,".pdf")))
    file.path(output, paste0(stem, c(".pdf", ".png")))
})

figure_files <- c(forest_file, rcs_file, comparison_file, multi_file, quartile_files,
                  dash_file, if (!is.null(scatter_data)) scatter_file, trend_file, cor_file,
                  cor_appendix_files)
stopifnot(all(file.exists(figure_files)), all(file.info(figure_files)$size > 0))
writeLines(c(
    "Run f/score.R -> score.R -> score.R.",
    "Legacy source: D:/projects/dose_v1/scripts/03_plot_mutiPattern.R.",
    "Included diets, labels and raw score variables follow pattern_spec.csv; outcomes follow analysis_config.rds.",
    "The primary score column and plot label are dose / DOSE. Range and scoring metadata come from pattern_spec.csv. All grouping figures use fitted Q1-Q4.",
    "Primary outcome is ischemic stroke; other stroke outcomes are supporting analyses.",
    "Models forest and consistency scatter each contain separate P10-P90 and per-SD pages.",
    "Quartile PDFs contain one outcome per page; DOSE-DASH combines quartiles and both continuous scales.",
    "Trend figure shows point estimates only; see quartile forest plots for 95% CIs.",
    "Consistency scatter is descriptive, not a statistical heterogeneity test; above y=x means lower IS HR.",
    "Association significance/FDR is read unchanged from s3_2 tables, not recalculated across outcomes.",
    "Correlation uses raw score columns from dat.pattern.rds, pairwise complete observations.",
    "Correlation clustering: Euclidean distance between correlation rows, complete linkage (legacy pheatmap defaults).",
    "Correlation heatmap stars use raw P; exported F matrix is BH-adjusted over unique off-diagonal pairs.",
    "The ten-pattern appendix follows the correlation heatmap, reuses its Spearman matrix, and exports PDF/PNG with Arial text plus 45 pairwise correlations.",
    "Current association ranking and apparent C-index do not establish independently validated superiority.",
    "", basename(figure_files)), file.path(DIR_FIG, "README_figures.txt"))
cat("Figures and tables complete:\n", paste(figure_files, collapse="\n"), "\n", sep="")
}))

# ===========================================================================
# score / outcomes
# ===========================================================================
dose_register_step("score", "outcomes", expression({
# ======================================================================
# Primary DOSE outcome-wide analysis; adapted from dose_v1/scripts/03_pattern_otherOutcomes.R.
# 运行顺序：foodwas / cohort -> score.R -> 本脚本。
# 不 source 上游分析，不依赖工作区中的 dat，也不修改原始表型文件。
#
# 主评分统一读取 dose；具体计分版本、成分、范围和步长从评分 RDS 的属性读取。
# 计分规则只在 score.R 中维护，其他版本仅保留作后续敏感性分析。
# Q1-Q4，P10-P90 / Per 1-SD，仅 M1；仅分析 DOSE。
# M1 与 score.R 一致；是否调整饮酒由 score_contains_alcohol 决定。
# 当前评分包含 Wine，因此移除单独的饮酒协变量。
# 不读取 M2/M3 协变量，也不因这些额外变量缺失而排除参与者。
# 所有人继承当前 dat.base.imp 的入组限制（包括排除基线总体卒中）；
# 再逐结局排除该疾病在 start_date 当天及以前的既往病。
# 前列腺癌仅纳入 sex.b == Male；此时性别为常数，自动从协变量中移除并记录。
# 其他疾病定义严格读取 icd10_other.lst，不用同名基线 covariate 替代。
# 尤其 lst 的 cvd = I00-I99，与 s0 构建的基线 cvd 定义不同。
# 随访截止 = min(失访, 死亡, cfg$date_follow_end)；诊断晚于截止则为删失。
# 这修正旧 t2e 对截止日期处理的问题，卒中结果可能因此与既有 s3_2 不同。
#
# 输出：ukb/04score/otherOutcomes/（不覆盖 pattern/ 的 IS 结果）。
# 包含 M1 长表、逐结局 QC、模型状态、格式化表格和可选环形柱形图。
# 本文件整合完整分析与两类绘图；默认完整分析，--plot-only 仅重画现有结果。
# 依赖：survival；图形使用基础 R，无需额外绘图包。
# ======================================================================

DIR_PROJECT <- dose_output_dir()
DIR_ROOT <- dose_drive()
DIR_PHE <- dose_ukb_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_SCORE <- file.path(DIR_PROJECT, "ukb/04score")
DIR_OUT <- file.path(DIR_SCORE, "otherOutcomes")
MODELS <- "M1"
MAKE_FIGURES <- TRUE
QUARTILES <- paste0("Q", 1:4)
SCALES <- c("P10-P90", "Per 1-SD")

CIRCULAR_ONLY <- "--circular-only" %in% commandArgs(trailingOnly=TRUE)
FOREST_ONLY <- "--forest-only" %in% commandArgs(trailingOnly=TRUE)
PLOT_ONLY <- FOREST_ONLY || CIRCULAR_ONLY || "--plot-only" %in% commandArgs(trailingOnly=TRUE) ||
    isTRUE(getOption("dose.other.plot_only", FALSE))

# Shared reverse reference palette: Fasting time (stroke) through Age.
pal <- c("Circulatory system disorders"="#DC666E", Infections="#B06D82",
        Cancer="#B58B2A", "Blood & immune"="#FFA834", Endocrine="#9A715E",
        "Mental & behavioural"="#A581D4", "Nervous system"="#B99CDD",
        "Eye & ear"="#BD67ED", Respiratory="#72AEC8", Digestive="#4F9F9A",
        Skin="#52A3CD", Musculoskeletal="#81BCAE", Genitourinary="#ABCE91", Other="#B8B8B8")

# Fig2 / foodwas.R colours sampled across its hue families.
# Circulatory outcomes (including stroke) retain red; dotplot palette is unchanged.
pal_circular <- c("Circulatory system disorders"="#DC666E",
    Infections="#006B66", Cancer="#3FA087", "Blood & immune"="#78C9A5",
    Endocrine="#3E96BA", "Mental & behavioural"="#084E79",
    "Nervous system"="#5D1987", "Eye & ear"="#B26DDE",
    Respiratory="#B87848", Digestive="#6A3317", Skin="#D96A00",
    Musculoskeletal="#FFC13A", Genitourinary="#E44A7D", Other="#B8B8B8")

# Integrated circular and dot-and-whisker plots; no external plot scripts.
plot_dose_other_outcomes <- function(dir_out, pal) {
    input <- file.path(dir_out, "Plot_DOSE_OutcomeWide_M1.csv")
    estimates <- file.path(dir_out, "DOSE_other_continuous.csv")
    d <- read.csv(input, stringsAsFactors=FALSE, check.names=FALSE)
    raw <- read.csv(estimates, stringsAsFactors=FALSE, check.names=FALSE)
    required <- c("outcome", "model", "scale", "HR", "LCI", "HCI", "P", "FDR")
    stopifnot(all(c(required, "label", "Category") %in% names(d)))
    key <- function(x) paste(x$outcome, x$model, x$scale, sep="|")
    stopifnot(!anyDuplicated(key(d)), !anyDuplicated(key(raw)))
    ix <- match(key(d), key(raw))
    stopifnot(!anyNA(ix), setequal(key(d), key(raw[raw$model == "M1", ])),
              isTRUE(all.equal(d[required], raw[ix, required], check.attributes=FALSE)))
    stopifnot(!anyNA(d$Category), all(d$Category %in% names(pal)))
    fig_dir <- file.path(dir_out, "figures")
    dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
    # Preserve existing figures before replacing them on the first redraw.
    old <- list.files(fig_dir, "^Fig_DOSE_OutcomeWide_M1_(P10P90|PerSD)\\.pdf$", full.names=TRUE)
    backup <- file.path(fig_dir, "before_style_redraw")
    if (length(old) && !dir.exists(backup)) {
        dir.create(backup)
        stopifnot(all(file.copy(old, backup, overwrite=FALSE)))
    }
    audit <- list()
    for (s in c("P10-P90", "Per 1-SD")) {
        a <- d[d$model == "M1" & d$scale == s, ]
        good <- with(a, is.finite(HR) & is.finite(LCI) & is.finite(HCI) &
                     LCI > 0 & LCI <= HR & HR <= HCI)
        r <- a[good, ]
        if (!nrow(r)) stop("No estimable outcomes: ", s)
        r$Category <- factor(r$Category, levels=names(pal))
        r <- r[order(r$Category, r$P, abs(r$HR-1)), ]
        sig <- !is.na(r$FDR) & r$FDR < .05
        # Six blank positions at the top; one blank between disease categories.
        gaps <- c(0L, cumsum(diff(as.integer(r$Category)) != 0L))
        xid <- 6 + seq_len(nrow(r)) + gaps
        total <- max(xid)
        axis_x <- 3.5
        angle <- pi/2 - 2*pi*(xid-axis_x)/total
        # Same linear 1-HR transform as the reference. HR > 1 extends inward.
        # Include the complete confidence-interval range, including HR > 1.
        hr_range <- range(c(r$LCI, r$HCI, 1))
        ticks <- sort(unique(c(1, pretty(hr_range, n=5))))
        ticks <- ticks[ticks > 0 & ticks >= hr_range[1]-.04 & ticks <= hr_range[2]+.04]
        span <- max(diff(range(c(hr_range, ticks))), .1)
        upper <- max(c(hr_range, ticks)) + .1*span
        # Compact ring with larger type; limits below account for full label widths.
        label_cex <- .88
        radius <- function(hr) .85*(.65 + 2.4*(upper-hr)/span)
        lo <- radius(r$HCI); hi <- radius(r$LCI); mid <- radius(r$HR)
        base <- radius(1)
        label_r <- pmax(hi + .18, base + .3)
        draw <- function() {
            layout(matrix(c(1, 2), nrow=1), widths=c(11.6, 4.4))
            par(mar=c(2.5, .5, 3.5, .5), family="sans")
            plot.new()
            # Measure labels in plot units and fit all of them inside the panel.
            limit <- max(label_r) + 2
            for (iteration in 1:12) {
                plot.window(c(-limit, limit), c(-limit, limit), asp=1, xaxs="i", yaxs="i")
                label_extent <- label_r + strwidth(r$label, cex=label_cex)
                limit <- max(abs(label_extent*cos(angle)), abs(label_extent*sin(angle))) + .3
            }
            plot.window(c(-limit, limit), c(-limit, limit), asp=1, xaxs="i", yaxs="i")
            circle <- seq(0, 2*pi, length.out=1200)
            for (tick in ticks) lines(radius(tick)*cos(circle), radius(tick)*sin(circle),
                col=if (tick == 1) "#B8BEC5" else "#ECECEC", lwd=if (tick == 1) .9 else .6,
                lty=if (tick == 1) 2 else 1)
            for (i in seq_len(nrow(r))) {
                arc <- seq(angle[i]-.85*pi/total, angle[i]+.85*pi/total, length.out=12)
                polygon(c(base*cos(arc), mid[i]*cos(rev(arc))),
                        c(base*sin(arc), mid[i]*sin(rev(arc))),
                        col=pal[as.character(r$Category[i])], border=NA)
                colour <- if (sig[i]) "#374151" else "#BDBDBD"
                segments(lo[i]*cos(angle[i]), lo[i]*sin(angle[i]),
                         hi[i]*cos(angle[i]), hi[i]*sin(angle[i]), col=colour, lwd=.85)
                for (rad in c(lo[i], hi[i])) {
                    cap <- seq(angle[i]-.3*pi/total, angle[i]+.3*pi/total, length.out=6)
                    lines(rad*cos(cap), rad*sin(cap), col=colour, lwd=.85)
                }
                left <- cos(angle[i]) < 0
                text(label_r[i]*cos(angle[i]), label_r[i]*sin(angle[i]), r$label[i],
                    srt=angle[i]*180/pi + if (left) 180 else 0,
                    adj=c(if (left) 1 else 0, .5), cex=label_cex,
                    col=if (sig[i]) "#2F3742" else "#9CA3AF")
            }
            segments(0, min(radius(ticks)), 0, max(radius(ticks)), col="#6B7280")
            segments(-.04, radius(ticks), .04, radius(ticks), col="#6B7280")
            text(.08, radius(ticks), sprintf("%.2f", ticks), adj=0, cex=.86, col="#4B5563")
            text(0, max(radius(ticks))+.23, "HR", font=2, cex=1, col="#2F3742")
            title(main="DOSE and health outcomes", cex.main=1.3, line=1.6)
            mtext(paste0("M1  |  ", if (s == "P10-P90") "10th-to-90th percentile increment" else "Per 1-SD increment"),
                side=3, line=.3, cex=.95, col="#4B5563")
            par(mar=c(2.5, 0, 3.5, 1))
            plot.new(); plot.window(c(0, 1), c(0, 1))
            used <- names(pal)[names(pal) %in% as.character(r$Category)]
            legend("topleft", legend=used, fill=pal[used], border=NA, bty="n",
                   title="Disease Category", title.font=2, cex=1.04, y.intersp=1.4)
            text(0, .32, "Bars: HR (linear scale)\nWhiskers: 95% CI\nDashed circle: HR = 1\nLower HR extends outward",
                 adj=c(0, 1), cex=.98, col="#4B5563")
            text(0, .19, "Dark labels: BH FDR < 0.05\nGrey labels: BH FDR >= 0.05\n(or unavailable FDR)",
                 adj=c(0, 1), cex=.98, col="#4B5563")
        }
        suffix <- paste0(if (s == "P10-P90") "P10P90" else "PerSD", "_Fig2Palette")
        stem <- file.path(fig_dir, paste0("Fig_DOSE_OutcomeWide_M1_", suffix))
        grDevices::cairo_pdf(paste0(stem, ".pdf"), width=16, height=12)
        tryCatch(draw(), finally=bounded_dev_off())
        bounded_png(paste0(stem, ".png"), width=3200, height=2400, res=200, type="cairo")
        tryCatch(draw(), finally=bounded_dev_off())
        audit[[s]] <- data.frame(scale=s, planned=nrow(a), plotted=nrow(r), omitted=sum(!good),
            FDR_significant=sum(sig), omitted_outcomes=paste(a$outcome[!good], collapse=";"))
        message("Redrawn: ", stem, " [PDF + PNG]")
    }
    write.csv(do.call(rbind, audit), file.path(fig_dir, "redraw_QC.csv"), row.names=FALSE)
    invisible(audit)
}

plot_dose_other_dotplot <- function(dir_out, pal) {
    d <- read.csv(file.path(dir_out, "Plot_DOSE_OutcomeWide_M1.csv"), stringsAsFactors=FALSE)
    raw <- read.csv(file.path(dir_out, "DOSE_other_continuous.csv"), stringsAsFactors=FALSE)
    fields <- c("outcome", "model", "scale", "HR", "LCI", "HCI", "P", "FDR")
    key <- function(z) paste(z$outcome, z$model, z$scale, sep="|")
    stopifnot(!anyDuplicated(key(d)), !anyDuplicated(key(raw)))
    ix <- match(key(d), key(raw))
    stopifnot(!anyNA(ix), setequal(key(d), key(raw[raw$model == "M1", ])),
              isTRUE(all.equal(d[fields], raw[ix, fields], check.attributes=FALSE)))
    stopifnot(!anyNA(d$Category), all(d$Category %in% names(pal)))
    fig_dir <- file.path(dir_out, "figures")
    dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
    audit <- list()
    for (s in c("P10-P90", "Per 1-SD")) {
        a <- d[d$model == "M1" & d$scale == s, ]
        good <- with(a, is.finite(HR) & is.finite(LCI) & is.finite(HCI) &
                     LCI > 0 & LCI <= HR & HR <= HCI)
        r <- a[good, ]
        if (!nrow(r)) stop("No estimable outcomes: ", s)
        # Keep the order within categories from the saved outcome list.
        r$Category <- factor(r$Category, levels=names(pal))
        r <- r[order(r$Category), ]
        stars <- ifelse(is.na(r$FDR), "", ifelse(r$FDR < .001, "***",
                         ifelse(r$FDR < .01, "**", ifelse(r$FDR < .05, "*", ""))))
        r$stars <- stars
        # All outcomes on one row; keep complete disease categories together.
        panels <- list(seq_len(nrow(r)))
        category_labels <- c("Circulatory system disorders"="Circulatory system",
            Infections="Infections", Cancer="Cancer", "Blood & immune"="Blood &\nimmune",
            Endocrine="Endocrine", "Mental & behavioural"="Mental &\nbehavioural",
            "Nervous system"="Nervous\nsystem", "Eye & ear"="Eye & ear",
            Respiratory="Respiratory", Digestive="Digestive", Skin="Skin",
            Musculoskeletal="Musculoskeletal", Genitourinary="Genitourinary", Other="Other")
        # Display-only short labels; exported key retains the original definitions.
        short <- c(cvd_stroke_i="Ischemic stroke", qi_stroke="Overall stroke",
            cvd_stroke_ih="Intracerebral hemorrhage", cvd_stroke_sh="Subarachnoid hemorrhage",
            infect="All infections", cancer="All cancers", ca_melanoma="Melanoma / other skin cancer",
            ca_leukemia="Hematologic malignancy", blood_immune="Blood / immune disorders",
            mental="Mental / behavioural disorders", neuro="Nervous system disorders",
            cvd="Circulatory disorders", resp="Respiratory disorders",
            flu_pneumonia="Influenza / pneumonia", copd="COPD", digestive="Digestive disorders",
            musculo="Musculoskeletal disorders", ra_disorder="Rheumatoid arthritis disorders",
            gu="Genitourinary disorders")
        r$display_label <- r$label
        matched <- r$outcome %in% names(short)
        r$display_label[matched] <- unname(short[r$outcome[matched]])
        hr_range <- range(c(r$LCI, r$HCI, 1))
        span <- max(diff(hr_range), .1)
        lower <- max(0, hr_range[1]-.09*span)
        upper <- hr_range[2]+.16*span
        ticks <- sort(unique(c(1, pretty(c(lower, upper), n=5))))
        ticks <- ticks[ticks >= lower & ticks <= upper]
        draw <- function() {
            par(mfrow=c(length(panels), 1), oma=c(1.2, 0, 2.4, 0),
                mar=c(10.5, 4.8, 3.8, 1.2), family="sans", cex=1)
            for (j in seq_along(panels)) {
                z <- r[panels[[j]], ]
                group <- as.character(z$Category)
                # Visible gaps separate categories; positions are independent of significance.
                x <- seq_len(nrow(z)) + .65*c(0, cumsum(group[-1] != group[-length(group)]))
                colour <- unname(pal[group])
                plot.new(); plot.window(c(.35, max(x)+.65), c(lower, upper), xaxs="i", yaxs="i")
                bounds <- which(group[-1] != group[-length(group)])
                if (length(bounds)) abline(v=(x[bounds]+x[bounds+1])/2, col="#E5E7EB", lwd=.7)
                abline(h=ticks[ticks != 1], col="#F0F1F3", lwd=.5)
                abline(h=1, lty=2, col="#5B6470", lwd=1.1)
                axis(2, at=ticks, las=1, cex.axis=.95)
                mtext("HR (95% CI)", side=2, line=3.1, cex=1)
                segments(x, z$LCI, x, z$HCI, col=adjustcolor(colour, alpha.f=.9), lwd=2)
                points(x, z$HR, pch=16, col=colour, cex=1.05)
                text(x, z$HCI+.03*span, z$stars, cex=.9, col="#111827")
                axis(1, at=x, labels=FALSE, tck=-.012)
                box(bty="l", col="#4B5563")
                text(x, lower-.045*span, z$display_label, srt=55, adj=1, xpd=NA, cex=.8)
                # Category headings directly above their outcomes, with matching short rules.
                for (g in unique(group)) {
                    gx <- x[group == g]
                    segments(min(gx)-.25, upper+.035*span, max(gx)+.25, upper+.035*span,
                        col=pal[g], lwd=3.5, xpd=NA)
                    text(mean(range(gx)), upper+.115*span, category_labels[g],
                        col="#303843", font=2, cex=.86, xpd=NA)
                }
            }
            mtext(paste0("DOSE and health outcomes | ",
                if (s == "P10-P90") "10th-to-90th percentile increment" else "Per 1-SD increment"),
                side=3, outer=TRUE, line=.7, cex=1.15, font=2)
            mtext("BH FDR: * < 0.05   ** < 0.01   *** < 0.001", side=1, outer=TRUE,
                line=.1, cex=.8, col="#4B5563")
        }
        suffix <- paste0(if (s == "P10-P90") "P10P90" else "PerSD", "_Publication")
        stem <- file.path(fig_dir, paste0("Fig_DOSE_OutcomeWide_M1_Dotplot_", suffix))
        grDevices::cairo_pdf(paste0(stem, ".pdf"), width=28, height=8)
        tryCatch(draw(), finally=bounded_dev_off())
        bounded_png(paste0(stem, ".png"), width=5600, height=1600, res=200, type="cairo")
        tryCatch(draw(), finally=bounded_dev_off())
        write.csv(r[c("outcome", "label", "display_label", "Category", "scale", "HR", "LCI", "HCI", "FDR", "stars")],
            paste0(stem, "_data.csv"), row.names=FALSE)
        audit[[s]] <- data.frame(scale=s, planned=nrow(a), plotted=nrow(r), omitted=sum(!good),
            FDR_significant=sum(nzchar(stars)), omitted_outcomes=paste(a$outcome[!good], collapse=";"))
        message("Saved: ", stem, " [PDF + PNG]")
    }
    write.csv(do.call(rbind, audit), file.path(fig_dir, "dotplot_QC.csv"), row.names=FALSE)
    invisible(audit)
}


# Table-style outcome-wide forest, adapted visually from the supplied Fig. 4.
# No re-fitting, subsetting by significance, or recalculation of FDR.
plot_dose_other_chapter_forest <- function(dir_out) {
    d <- read.csv(file.path(dir_out, "DOSE_other_continuous.csv"), stringsAsFactors=FALSE)
    chapters <- list(
        `Infectious diseases`=c("infect", "infect_bact", "infect_viral"),
        Cancer=c("cancer", "ca_breast", "ca_crc", "ca_leukemia", "ca_lung", "ca_melanoma", "ca_esophagus", "ca_prostate", "ca_stomach"),
        `Endocrine disorders`=c("endo", "dm", "obesity"),
        `Circulatory system disorders`=c("cvd", "arrhythmia", "cvd_afib", "cerebro", "hf", "htn", "ihd", "pad", "cvd_vte", "qi_stroke", "cvd_stroke_i", "cvd_stroke_ih", "cvd_stroke_sh"),
        `Respiratory system disorders`=c("resp", "asthma", "copd", "flu_pneumonia"),
        `Digestive system disorders`=c("digestive", "ibd", "liver", "pancreatitis"),
        `Genitourinary system disorders`=c("gu", "renal_failure", "ckd"),
        `Musculoskeletal disorders`=c("musculo", "oa", "oporosis", "ra_disorder", "sciatica", "soft_tissue"),
        `Nervous system disorders`=c("neuro", "epilepsy", "parkinson", "sleep_disorder", "tia"),
        `Mental and behavioural disorders`=c("mental", "dementia", "mood", "neurotic", "psychotic", "substance_abuse"),
        `Eye disorders`=c("eye", "cataract", "glaucoma"), `Ear disorders`="ear",
        `Skin disorders`=c("skin", "skin_infect_eczema"),
        `Blood and immune disorders`=c("blood_immune", "anemia"))
    present <- unique(d$outcome)
    extra <- setdiff(present, unname(unlist(chapters)))
    if (length(extra)) chapters$Other <- extra
    chapters <- lapply(chapters, intersect, present)
    chapters <- chapters[lengths(chapters) > 0L]
    ord <- unname(unlist(chapters))
    stopifnot(!anyDuplicated(ord))
    fig_dir <- file.path(dir_out, "figures")
    dir.create(fig_dir, recursive=TRUE, showWarnings=FALSE)
    blue <- "#3D6388"
    for (s in c("P10-P90", "Per 1-SD")) {
        r <- d[d$model == "M1" & d$scale == s, ]
        stopifnot(!anyDuplicated(r$outcome), setequal(r$outcome, ord))
        r <- r[match(ord, r$outcome), ]
        r$Disease_Category <- rep(names(chapters), lengths(chapters))
        r$significant <- !is.na(r$FDR) & r$FDR < .05
        valid <- with(r, is.finite(HR) & is.finite(LCI) & is.finite(HCI) & LCI > 0 & LCI <= HR & HR <= HCI)
        stopifnot(any(valid))
        r$HR_CI <- ifelse(valid, sprintf("%.2f (%.2f, %.2f)", r$HR, r$LCI, r$HCI), "Not estimable")
        groupid <- rep(seq_along(chapters), lengths(chapters))
        y <- rev(seq_len(nrow(r))) + .6*(length(chapters)-groupid)
        yrange <- c(-2.6, max(y)+2.4)
        # Linear HR axis includes every CI; no interval is clipped.
        hrange <- range(c(1, r$LCI[valid], r$HCI[valid]))
        ticks <- sort(unique(c(1, round(pretty(hrange, n=5), 8))))
        xlim <- range(c(ticks, hrange))
        map <- function(hr) .55 + (hr-xlim[1])/diff(xlim)*.30
        draw <- function() {
            par(mar=c(.6, .6, .6, .6), family="sans", ps=12, cex=1, xaxs="i", yaxs="i")
            plot.new(); plot.window(c(0, 1), yrange)
            text(0, max(y)+2, "DOSE and health outcomes", adj=0, font=2, cex=1)
            text(1, max(y)+2, if (s == "P10-P90") "M1 | 10th-to-90th percentile increment" else "M1 | Per 1-SD increment", adj=1, cex=1)
            text(c(0,.235,.865), max(y)+.95, c("Disease Category", "Outcome", "Hazard ratio (95% CI)"), adj=0, font=2, cex=1)
            for (g in seq_along(chapters)) {
                ix <- which(groupid == g)
                bottom <- min(y[ix])-.48; top <- max(y[ix])+.48
                rect(0, bottom, .23, top, col="#EBEEF2", border=NA)
                # General chapter outcome is the first row, as in the reference.
                rect(.235, y[ix[1]]-.48, .995, y[ix[1]]+.48, col="#EBEEF2", border=NA)
                for (tick in ticks) segments(map(tick), bottom, map(tick), top,
                    col=if (tick == 1) "#7A818A" else "#DCE0E5", lty=if (tick == 1) 2 else 1,
                    lwd=if (tick == 1) 1.2 else .7)
                chapter <- names(chapters)[g]
                if (nchar(chapter) > 28) chapter <- paste(strwrap(chapter, width=27), collapse="\n")
                text(.006, mean(range(y[ix])), chapter, adj=c(0,.5), cex=1, font=2)
            }
            for (i in seq_len(nrow(r))) {
                # All text is 12 pt; the wider canvas fits full outcome names.
                text(.235, y[i], r$label[i], adj=0, cex=1)
                text(.865, y[i], r$HR_CI[i], adj=0, cex=1)
                if (valid[i]) {
                    segments(map(r$LCI[i]), y[i], map(r$HCI[i]), y[i], col="#41474E", lwd=1.25)
                    points(map(r$HR[i]), y[i], pch=21, cex=1, lwd=1.1,
                        col=if (r$significant[i]) blue else "#41474E",
                        bg=if (r$significant[i]) blue else "white")
                }
            }
            for (tick in ticks) {
                segments(map(tick), -.55, map(tick), -.75, col="#707780")
                text(map(tick), -1.05, format(tick, trim=TRUE), cex=1)
            }
            text(.7, -1.8, "HR (95% CI)", font=2, cex=1)
            text(0, -2.45, "Filled blue: BH FDR < 0.05; open: FDR >= 0.05 or unavailable. Lines: 95% CI; dashed reference: HR = 1.", adj=0, cex=1)
        }
        stem <- file.path(fig_dir, paste0("Fig_DOSE_OutcomeWide_M1_ChapterForest_", if (s == "P10-P90") "P10P90" else "PerSD"))
        cairo_pdf(paste0(stem,".pdf"), width=18, height=20, pointsize=12)
        tryCatch(draw(), finally=bounded_dev_off())
        bounded_png(paste0(stem,".png"), width=3600, height=4000, res=200, type="cairo", pointsize=12)
        tryCatch(draw(), finally=bounded_dev_off())
        write.csv(r, paste0(stem,"_data.csv"), row.names=FALSE)
        message("Saved chapter forest: ", stem)
    }
}


if (!PLOT_ONLY) {
if (!requireNamespace("survival", quietly=TRUE)) stop("Package survival is required.")
library(survival)
inputs <- c(cohort=file.path(DIR_DATA, "dat.base.imp.rds"),
    config=file.path(DIR_DATA, "analysis_config.rds"),
    dose=file.path(DIR_SCORE, "dose_score.rds"),
    definitions=Sys.getenv("DOSE_OUTCOME_DEFINITIONS", dose_definition_file("icd10_other")),
    codes=file.path(DIR_PHE, "Rdata/icd10.code0.rds"),
    dates=file.path(DIR_PHE, "Rdata/icd10.date0.rds"))
if (any(!file.exists(inputs))) stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse=", "))
require_cols <- function(d, cols, context) {
    missing <- setdiff(cols, names(d))
    if (length(missing)) stop(context, ": missing ", paste(missing, collapse=", "))
}
check_ids <- function(d, context) {
    require_cols(d, "eid", context)
    d$eid <- as.character(d$eid)
    if (anyNA(d$eid) || any(!nzchar(trimws(d$eid))) || anyDuplicated(d$eid)) {
        stop(context, ": eid must be unique and non-missing.")
    }
    d
}
as_date <- function(x) {
    if (inherits(x, "Date")) return(x)
    if (inherits(x, "POSIXt")) return(as.Date(x))
    if (is.numeric(x)) return(as.Date(x, origin="1970-01-01"))
    z <- as.Date(as.character(x))
    if (any(!is.na(x) & is.na(z))) stop("Unrecognised non-missing date.")
    z
}
bind <- function(x) { z <- do.call(rbind, x); rownames(z) <- NULL; z }
cfg <- readRDS(inputs[["config"]])
dat <- check_ids(as.data.frame(readRDS(inputs[["cohort"]])), "cohort")
dose <- readRDS(inputs[["dose"]])
component_cols <- attr(dose, "component_columns", exact=TRUE)
dose_meta <- list(method=attr(dose, "primary_scoring_method"),
    range=attr(dose, "score_theoretical_range"), step=attr(dose, "score_step"),
    alcohol=attr(dose, "score_contains_alcohol"), component_steps=attr(dose, "component_steps"))
if (!identical(attr(dose, "score_version"), "dose") ||
    !identical(attr(dose, "score_intake_unit"), "servings/day") ||
    length(dose_meta$method) != 1L || is.na(dose_meta$method) ||
    length(dose_meta$range) != 2L || any(!is.finite(dose_meta$range)) || diff(dose_meta$range) <= 0 ||
    length(dose_meta$step) != 1L || !is.finite(dose_meta$step) || dose_meta$step <= 0 ||
    !is.logical(dose_meta$alcohol) || length(dose_meta$alcohol) != 1L || is.na(dose_meta$alcohol) ||
    !length(component_cols) || anyDuplicated(component_cols) ||
    !all(c("eid", "dose", component_cols) %in% names(dose)) ||
    !identical(names(dose_meta$component_steps), component_cols) ||
    any(!is.finite(dose_meta$component_steps) | dose_meta$component_steps <= 0))
    stop("Invalid primary DOSE metadata. Run the current score.R first.")
dose$eid <- as.character(dose$eid)
if (anyNA(dose$eid) || anyDuplicated(dose$eid) || !setequal(dat$eid,dose$eid)) stop("Score/cohort ID mismatch.")
valid_components <- vapply(component_cols, function(v) is.numeric(dose[[v]]) &&
    all(is.na(dose[[v]]) | dose[[v]] %in% seq(0, 1, by=dose_meta$component_steps[[v]])), logical(1))
if (!all(valid_components) || !is.numeric(dose$dose) ||
    !all(is.na(dose$dose) | dose$dose %in% seq(dose_meta$range[1], dose_meta$range[2], by=dose_meta$step)) ||
    !isTRUE(all.equal(dose$dose, rowSums(dose[component_cols]), check.attributes=FALSE)))
    stop("DOSE must match its component sum, range and step metadata.")
if ("diet.dose.sum" %in% names(dose) && !identical(dose$dose, dose$diet.dose.sum))
    stop("The compatibility alias diet.dose.sum differs from dose.")
dose$diet.dose.sum <- dose$dose
dat$dose <- dose$dose[match(dat$eid, dose$eid)]
dat$diet.dose.sum <- dat$dose # Compatibility alias; all models use dose.
rm(dose)
if (!all(c("covs.m1", "covs.alcohol", "date_follow_end", "Ys", "out_map") %in% names(cfg))) {
    stop("analysis_config.rds is not the current s0_2 configuration.")
}
# Same M1 alcohol handling as s3_2; no dependency on a versioned score name.
m1_covs <- if (dose_meta$alcohol) setdiff(unique(cfg$covs.m1), cfg$covs.alcohol) else
    unique(c(cfg$covs.m1, cfg$covs.alcohol))
if (!length(cfg$covs.m1) || !length(m1_covs) ||
    anyNA(m1_covs) || any(!nzchar(m1_covs))) stop("Invalid M1 covariate configuration.")
stroke_outcomes <- c("cvd_stroke_i", "qi_stroke", "cvd_stroke_ih", "cvd_stroke_sh")
if (!all(stroke_outcomes %in% cfg$Ys) || anyNA(cfg$out_map[stroke_outcomes])) stop("Missing stroke configuration.")
stroke_dates <- paste0("fod_icd10_", stroke_outcomes)
require_cols(dat, c(m1_covs, "start_date", "date_lost", "date_death", stroke_dates), "cohort")
for (v in c("start_date", "date_lost", "date_death", stroke_dates)) dat[[v]] <- as_date(dat[[v]])
end_date <- as_date(cfg$date_follow_end)
if (length(end_date) != 1L || is.na(end_date)) stop("Invalid administrative follow-up end date.")
dat$end0 <- pmin(dat$date_lost, dat$date_death, end_date, na.rm=TRUE)

# Same full-cohort scaling and quantile/cut rule as s3_2, fixed across outcomes.
# Do not recalculate quartile cutoffs or SD after excluding each outcome's cases.
x <- dat$dose
if (!is.numeric(x) || any(!is.na(x) & !x %in% seq(dose_meta$range[1], dose_meta$range[2], by=dose_meta$step))) stop("Invalid primary DOSE score.")
br <- unname(quantile(x, c(0, .25, .5, .75, 1), na.rm=TRUE))
p1090 <- unname(quantile(x, c(.1, .9), na.rm=TRUE))
score_sd <- sd(x, na.rm=TRUE)
if (anyDuplicated(br) || !is.finite(score_sd) || score_sd <= 0 || diff(p1090) <= 0) stop("Degenerate DOSE score distribution.")
dat$diet.dose.q <- cut(x, br, labels=QUARTILES, include.lowest=TRUE, right=TRUE)
contrasts(dat$diet.dose.q) <- contr.treatment(QUARTILES, base=1)
increments <- setNames(c(diff(p1090), score_sd), SCALES)
group_def <- data.frame(diet="dose", scoring_method=dose_meta$method, category=QUARTILES, cut_lower=br[1:4], cut_upper=br[2:5],
    lower_inclusive=c(TRUE, FALSE, FALSE, FALSE), upper_inclusive=TRUE)
group_def$N <- as.integer(table(dat$diet.dose.q))
group_def$score_min <- as.numeric(tapply(x, dat$diet.dose.q, min))
group_def$score_max <- as.numeric(tapply(x, dat$diet.dose.q, max))
cat("Primary score: dose (", dose_meta$method, "); ", nrow(dat), " participants; quartiles: ",
    paste(paste0(group_def$score_min, "-", group_def$score_max), collapse=" / "), "\n", sep="")

# Keep the actual regex/name/label list. No assumed aliases such as dm=t2dm.
fn <- read.delim(inputs[["definitions"]], header=FALSE, sep="\t", quote="", comment.char="",
                 fill=TRUE, stringsAsFactors=FALSE)
if (ncol(fn) < 3L) stop("Expected pattern, name and label in icd10_other.lst.")
fn <- fn[1:3]; names(fn) <- c("pattern", "name", "label")
fn[] <- lapply(fn, trimws)
if (anyNA(fn) || any(!nzchar(as.matrix(fn))) || anyDuplicated(fn$name)) stop("Invalid outcome definitions.")
fn <- fn[!fn$name %in% c("stroke", stroke_outcomes), ]
outcomes <- rbind(data.frame(pattern=NA_character_, name=stroke_outcomes,
    label=unname(cfg$out_map[stroke_outcomes]), stringsAsFactors=FALSE), fn)
outcomes$is_stroke <- outcomes$name %in% stroke_outcomes

# Extract first diagnosis dates using paired ICD code/date slots. Match by eid
# before extraction; missing participants must not silently become non-cases.
needed <- unique(c("eid", m1_covs, "sex.b", "dose", "diet.dose.sum", "diet.dose.q",
    "start_date", "date_lost", "date_death", "end0", stroke_dates,
    intersect(unlist(lapply(stroke_outcomes, function(y) paste0(y, c(".t2e", ".Yt2e")))), names(dat))))
dat <- dat[needed]
invisible(gc())
code0 <- check_ids(as.data.frame(readRDS(inputs[["codes"]])), "ICD codes")
if (!all(dat$eid %in% code0$eid)) stop("ICD codes do not cover the prepared cohort.")
code0 <- code0[match(dat$eid, code0$eid), , drop=FALSE]
invisible(gc())
date0 <- check_ids(as.data.frame(readRDS(inputs[["dates"]])), "ICD dates")
if (!all(dat$eid %in% date0$eid)) stop("ICD dates do not cover the prepared cohort.")
date0 <- date0[match(dat$eid, date0$eid), , drop=FALSE]
invisible(gc())
slots <- setdiff(names(code0), "eid")
date_slots <- setdiff(names(date0), "eid")
# Current D:/scripts/ukb/phe.R writes code slots i1,i2,... and dates d1,d2,... .
# Match their numeric suffixes explicitly instead of trusting physical order.
if (!length(slots) || !all(grepl("^i[0-9]+$", slots)) ||
    !all(grepl("^d[0-9]+$", date_slots)) ||
    !setequal(sub("^i", "", slots), sub("^d", "", date_slots))) {
    stop("Unpaired ICD code/date slots; expected iN/dN from the current phenotype pipeline.")
}
slots <- slots[order(as.integer(sub("^i", "", slots)))]
date_slots <- paste0("d", sub("^i", "", slots))
if (!all(dat$eid %in% code0$eid) || !all(dat$eid %in% date0$eid)) stop("ICD inputs do not cover the prepared cohort.")
code0 <- code0[match(dat$eid, code0$eid), slots, drop=FALSE]
date0 <- date0[match(dat$eid, date0$eid), date_slots, drop=FALSE]
date0[] <- lapply(date0, as_date)
fod <- data.frame(eid=dat$eid)
extraction_qc <- vector("list", nrow(fn))
undated <- list()
for (i in seq_len(nrow(fn))) {
    first <- rep(Inf, nrow(dat)); unknown_date <- rep(FALSE, nrow(dat)); hits <- integer(nrow(dat))
    for (j in seq_along(slots)) {
        hit <- !is.na(code0[[j]]) & grepl(fn$pattern[i], as.character(code0[[j]]), perl=TRUE)
        hits <- hits + as.integer(hit)
        unknown_date <- unknown_date | (hit & is.na(date0[[j]]))
        good <- hit & !is.na(date0[[j]])
        first[good] <- pmin(first[good], as.numeric(date0[[j]][good]))
    }
    first[!is.finite(first)] <- NA_real_
    fod[[paste0("fod_icd10_", fn$name[i])]] <- as.Date(first, origin="1970-01-01")
    # Any matching undated record leaves first-onset timing uncertain, even
    # if another matching record is dated. Exclude it from incident analyses.
    undated[[fn$name[i]]] <- unknown_date
    extraction_qc[[i]] <- data.frame(outcome=fn$name[i], matched_N=sum(hits > 0),
        dated_N=sum(!is.na(first)), undated_record_N=sum(unknown_date))
    cat("Extracted: ", fn$name[i], "\n", sep="")
}
rm(code0, date0)
invisible(gc())
# Replace same-named date columns deliberately, never create .x/.y duplicates.
dat[names(fod)[-1]] <- fod[-1]
dir.create(DIR_OUT, recursive=TRUE, showWarnings=FALSE)
dose_save_rds(fod, file.path(DIR_OUT, "fod.icd10.other.rds"))
write.csv(bind(extraction_qc), file.path(DIR_OUT, "ICD_extraction_QC.csv"), row.names=FALSE)
write.csv(outcomes, file.path(DIR_OUT, "outcome_label_list.csv"), row.names=FALSE)
write.csv(group_def, file.path(DIR_OUT, "DOSE_quartile_definitions.csv"), row.names=FALSE)
write.csv(data.frame(diet="dose", scoring_method=dose_meta$method, P10=p1090[1], P90=p1090[2], SD=score_sd,
    P10_P90_increment=diff(p1090)), file.path(DIR_OUT, "DOSE_score_scaling.csv"), row.names=FALSE)
rm(fod)

# Log failures and preserve NA estimates; never report a failed fit as null HR.
fit_log <- list()
safe_fit <- function(d, terms, context) {
    warnings <- character(0); error <- ""
    fit <- tryCatch(withCallingHandlers(
        coxph(reformulate(terms, response="Surv(time, event)"), data=d,
              ties="efron", singular.ok=FALSE, na.action=na.fail),
        warning=function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
        error=function(e) { error <<- conditionMessage(e); NULL })
    if (!is.null(fit) && (any(!is.finite(coef(fit))) || any(!is.finite(vcov(fit))))) {
        error <- "Non-finite coefficients/covariance"; fit <- NULL
    }
    # Convergence warnings can produce finite but unreliable estimates.
    if (any(grepl("converg|infinite|singular|iteration", warnings, ignore.case=TRUE))) {
        error <- paste(c(error, "Unreliable fit: convergence/singularity warning"), collapse=" "); fit <- NULL
    }
    fit_log[[length(fit_log)+1L]] <<- data.frame(context=context,
        status=if (is.null(fit)) "failed" else if (length(warnings)) "warning" else "ok",
        error=error, warnings=paste(unique(warnings), collapse=" | "))
    fit
}
effect <- function(fit, term, increment=1) {
    missing <- data.frame(HR=NA_real_, LCI=NA_real_, HCI=NA_real_, P=NA_real_)
    if (is.null(fit) || !term %in% names(coef(fit))) return(missing)
    b <- unname(coef(fit)[term]); se <- sqrt(vcov(fit)[term, term])
    if (!is.finite(se) || se <= 0) return(missing)
    z <- data.frame(HR=exp(b*increment), LCI=exp((b-1.959964*se)*increment),
        HCI=exp((b+1.959964*se)*increment), P=2*pnorm(-abs(b/se)))
    if (any(!is.finite(as.matrix(z))) || any(z[c("HR", "LCI", "HCI")] <= 0)) return(missing)
    z
}
ph_p <- function(fit, term) {
    if (is.null(fit)) return(NA_real_)
    tryCatch(suppressWarnings(unname(cox.zph(fit, terms=TRUE)$table[term, "p"])), error=function(e) NA_real_)
}
continuous <- categories <- overall <- cohort_qc <- model_spec <- stroke_qc <- list()
k <- 0L
for (Y in outcomes$name) {
    label <- outcomes$label[match(Y, outcomes$name)]
    onset <- dat[[paste0("fod_icd10_", Y)]]
    prior <- !is.na(onset) & !is.na(dat$start_date) & onset <= dat$start_date
    # C61 is analysed in the male risk population. Including women with no
    # prostate-cancer events causes separation in the adjusted sex coefficient.
    sex_restriction <- if (Y == "ca_prostate") "Male" else "None"
    eligible_sex <- if (Y == "ca_prostate") !is.na(dat$sex.b) & as.character(dat$sex.b) == "Male" else rep(TRUE, nrow(dat))
    unknown <- if (Y %in% names(undated)) undated[[Y]] else rep(FALSE, nrow(dat))
    event <- as.integer(!is.na(onset) & onset > dat$start_date & onset <= dat$end0)
    finish <- dat$end0; finish[which(event == 1L)] <- onset[which(event == 1L)]
    time <- as.numeric(finish-dat$start_date)/365.25
    valid_follow <- !is.na(dat$start_date) & is.finite(time) & time > 0 & !is.na(event)
    complete <- complete.cases(dat[c(m1_covs, "dose", "diet.dose.q")])
    for (v in m1_covs) if (is.numeric(dat[[v]])) complete <- complete & is.finite(dat[[v]])
    ok <- !prior & eligible_sex & !unknown & valid_follow & complete
    d <- dat[ok, c(m1_covs, "dose", "diet.dose.q"), drop=FALSE]
    d[m1_covs] <- lapply(d[m1_covs], function(z) if (is.factor(z)) droplevels(z) else z)
    d$time <- time[ok]; d$event <- event[ok]
    d$pheno <- d$dose
    d$group <- factor(d$diet.dose.q, levels=QUARTILES)
    contrasts(d$group) <- contr.treatment(QUARTILES, base=1)
    medians <- tapply(d$pheno, d$group, median)
    d$trend <- unname(medians[as.character(d$group)])
    cohort_qc[[Y]] <- data.frame(outcome=Y, label=label, input_N=nrow(dat), sex_restriction=sex_restriction,
        prior_N=sum(prior), sex_excluded_after_prior_N=sum(!prior & !eligible_sex),
        undated_after_prior_N=sum(!prior & eligible_sex & unknown),
        invalid_follow_after_exclusions_N=sum(!prior & eligible_sex & !unknown & !valid_follow),
        incomplete_after_exclusions_N=sum(!prior & eligible_sex & !unknown & valid_follow & !complete),
        N=nrow(d), Events=sum(d$event), person_years=sum(d$time),
        diagnoses_after_censor_N=sum(!is.na(onset) & onset > dat$end0, na.rm=TRUE))
    if (Y %in% stroke_outcomes && all(paste0(Y, c(".t2e", ".Yt2e")) %in% names(dat))) {
        old_t <- dat[[paste0(Y, ".t2e")]][ok]; old_y <- dat[[paste0(Y, ".Yt2e")]][ok]
        stroke_qc[[Y]] <- data.frame(outcome=Y, N=sum(ok),
            changed_event_N=sum(is.na(old_y) | old_y != d$event),
            changed_time_N=sum(is.na(old_t) | abs(old_t-d$time) > 1e-8))
    }
    for (m in MODELS) {
        k <- k+1L
        # Outcome-specific exclusions may leave a constant M1 covariate.
        # Remove only constants and report every omission.
        requested <- m1_covs
        constant <- requested[vapply(d[requested], function(z) length(unique(z)) < 2L, logical(1))]
        co <- setdiff(requested, constant)
        model_spec[[k]] <- data.frame(outcome=Y, model=m, requested=paste(requested, collapse=";"),
            covariates=paste(co, collapse=";"), omitted_constant=paste(constant, collapse=";"), region_stratified=FALSE)
        context <- paste(Y, m)
        eligible <- nrow(d) > 0L && sum(d$event) > 0L && length(unique(d$pheno)) > 1L
        fl <- fc <- ft <- f0 <- NULL
        if (eligible) {
            fl <- safe_fit(d, c("pheno", co), paste(context, "continuous"))
            if (all(table(d$group) > 0L)) {
                fc <- safe_fit(d, c("group", co), paste(context, "quartiles"))
                ft <- safe_fit(d, c("trend", co), paste(context, "trend"))
                f0 <- safe_fit(d, co, paste(context, "baseline"))
            } else fit_log[[length(fit_log)+1L]] <- data.frame(context=paste(context, "quartiles/trend"),
                status="skipped", error="One or more quartiles are empty", warnings="")
        } else fit_log[[length(fit_log)+1L]] <- data.frame(context=context,
            status="skipped", error="No participants/events or constant exposure", warnings="")
        meta <- data.frame(diet="dose", scoring_method=dose_meta$method, label=label, outcome=Y, model=m, N=nrow(d), Events=sum(d$event))
        ph_linear <- ph_p(fl, "pheno")
        continuous[[k]] <- bind(lapply(SCALES, function(s) cbind(meta, scale=s,
            increment=unname(increments[s]), PH_P=ph_linear, effect(fl, "pheno", increments[s]))))
        categories[[k]] <- bind(lapply(QUARTILES, function(q) {
            z <- d[d$group == q, , drop=FALSE]
            ef <- if (q == "Q1") data.frame(HR=1, LCI=1, HCI=1, P=NA_real_) else effect(fc, paste0("group", q))
            cbind(meta, category=q, reference="Q1", Group_N=nrow(z), Group_Events=sum(z$event),
                Person_years=sum(z$time), Incidence_percent=if (nrow(z)) 100*mean(z$event) else NA_real_,
                score_min=if (nrow(z)) min(z$pheno) else NA_real_,
                score_max=if (nrow(z)) max(z$pheno) else NA_real_,
                score_median=if (nrow(z)) median(z$pheno) else NA_real_, ef)
        }))
        p_global <- if (!is.null(fc) && !is.null(f0)) {
            df <- attr(logLik(fc), "df")-attr(logLik(f0), "df")
            pchisq(max(0, 2*as.numeric(logLik(fc)-logLik(f0))), df=df, lower.tail=FALSE)
        } else NA_real_
        overall[[k]] <- cbind(meta, P_global=p_global, P_trend=effect(ft, "trend")$P, PH_P=ph_p(fc, "group"))
    }
    cat("Completed: ", Y, " | N=", nrow(d), "; events=", sum(d$event), "\n", sep="")
}
continuous <- bind(continuous); categories <- bind(categories); overall <- bind(overall)

# Outcome-wide FDR differs intentionally from s3_2's across-pattern FDR.
# Families: all planned outcomes per model/scale, or per model/quartile.
# Failed outcomes remain in the family size; their P and FDR stay missing.
adjust_fdr <- function(d, p, by) {
    out <- rep(NA_real_, nrow(d))
    for (ix in split(seq_len(nrow(d)), interaction(d[by], drop=TRUE))) {
        out[ix] <- p.adjust(d[[p]][ix], method="BH", n=length(ix))
    }
    out
}
continuous$FDR <- adjust_fdr(continuous, "P", c("model", "scale"))
categories$FDR <- adjust_fdr(categories, "P", c("model", "category"))
overall$FDR_global <- adjust_fdr(overall, "P_global", "model")
overall$FDR_trend <- adjust_fdr(overall, "P_trend", "model")
write.csv(continuous, file.path(DIR_OUT, "DOSE_other_continuous.csv"), row.names=FALSE)
write.csv(categories, file.path(DIR_OUT, "DOSE_other_categories.csv"), row.names=FALSE)
write.csv(overall, file.path(DIR_OUT, "DOSE_other_overall_tests.csv"), row.names=FALSE)
write.csv(bind(cohort_qc), file.path(DIR_OUT, "analysis_cohort_QC.csv"), row.names=FALSE)
write.csv(bind(model_spec), file.path(DIR_OUT, "model_specification.csv"), row.names=FALSE)
write.csv(bind(fit_log), file.path(DIR_OUT, "model_fit_status.csv"), row.names=FALSE)
if (length(stroke_qc)) write.csv(bind(stroke_qc), file.path(DIR_OUT, "stroke_followup_comparison_QC.csv"), row.names=FALSE)

# Formatted tables: both continuous scales, observed quartile ranges and trend.
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "<0.001", sprintf("%.3f", p)))
fmt_hr <- function(h, l, u) ifelse(is.na(h), "Not estimable", sprintf("%.2f (%.2f-%.2f)", h, l, u))
table_rows <- list(); k <- 0L
for (Y in outcomes$name) for (m in MODELS) {
    a <- continuous[continuous$outcome == Y & continuous$model == m, ]
    q <- categories[categories$outcome == Y & categories$model == m, ]
    t <- overall[overall$outcome == Y & overall$model == m, ]
    intervals <- paste0(group_def$score_min, "-", group_def$score_max)
    k <- k+1L
    table_rows[[k]] <- data.frame(outcome=Y, Outcome=unname(outcomes$label[match(Y, outcomes$name)]), model=m,
        DOSE_score=c(a$scale, paste(q$category, intervals[match(q$category, QUARTILES)]), "P for trend"),
        N=c(a$N, q$Group_N, t$N), Events=c(a$Events, q$Group_Events, t$Events),
        HR_CI=c(fmt_hr(a$HR, a$LCI, a$HCI), ifelse(q$category == "Q1", "Ref.", fmt_hr(q$HR, q$LCI, q$HCI)), ""),
        P_value=fmt_p(c(a$P, q$P, t$P_trend)), FDR=fmt_p(c(a$FDR, q$FDR, t$FDR_trend)))
}
tables <- bind(table_rows)
write.csv(tables, file.path(DIR_OUT, "Table_DOSE_OutcomeWide_M1.csv"), row.names=FALSE)
write.csv(tables[!tables$outcome %in% stroke_outcomes, ],
    file.path(DIR_OUT, "Table_DOSE_OtherOutcomes_M1_excluding_stroke.csv"), row.names=FALSE)

# Disease categories retained from the legacy figure; shared reverse reference palette below.
disease_groups <- list(
    `Circulatory system disorders`=c(stroke_outcomes, "cvd", "htn", "ihd", "arrhythmia", "cvd_afib", "hf", "cvd_vte", "cerebro", "pad"),
    Infections=c("infect", "infect_bact", "infect_viral"),
    Cancer=c("cancer", "ca_crc", "ca_esophagus", "ca_stomach", "ca_lung", "ca_melanoma", "ca_breast", "ca_prostate", "ca_leukemia"),
    `Blood & immune`=c("blood_immune", "anemia"), Endocrine=c("endo", "dm", "obesity"),
    `Mental & behavioural`=c("mental", "dementia", "substance_abuse", "psychotic", "mood", "neurotic"),
    `Nervous system`=c("neuro", "parkinson", "epilepsy", "tia", "sleep_disorder"),
    `Eye & ear`=c("eye", "cataract", "glaucoma", "ear"),
    Respiratory=c("resp", "flu_pneumonia", "copd", "asthma"),
    Digestive=c("digestive", "ibd", "liver", "pancreatitis"), Skin=c("skin", "skin_infect_eczema"),
    Musculoskeletal=c("musculo", "ra_disorder", "oa", "sciatica", "soft_tissue", "oporosis"),
    Genitourinary=c("gu", "renal_failure", "ckd"))
plot_data <- continuous[continuous$model == "M1", ]
plot_data$Category <- "Other"
for (g in names(disease_groups)) plot_data$Category[plot_data$outcome %in% disease_groups[[g]]] <- g
plot_data$Category <- factor(plot_data$Category, levels=names(pal))
plot_data$significant <- !is.na(plot_data$FDR) & plot_data$FDR < .05
write.csv(plot_data, file.path(DIR_OUT, "Plot_DOSE_OutcomeWide_M1.csv"), row.names=FALSE)
if (MAKE_FIGURES) {
    if (!FOREST_ONLY) plot_dose_other_outcomes(DIR_OUT, pal_circular)
    if (!CIRCULAR_ONLY && !FOREST_ONLY) plot_dose_other_dotplot(DIR_OUT, pal)
    if (!CIRCULAR_ONLY) plot_dose_other_chapter_forest(DIR_OUT)
}
write.csv(data.frame(input=names(inputs), path=unname(inputs), md5=unname(tools::md5sum(inputs))),
    file.path(DIR_OUT, "input_manifest.csv"), row.names=FALSE)
code_files <- file.path(dose_script_dir(), c("score.R", "score.R"))
write.csv(data.frame(script=basename(code_files), path=code_files, md5=unname(tools::md5sum(code_files))),
    file.path(DIR_OUT, "code_manifest.csv"), row.names=FALSE)
capture.output(sessionInfo(), file=file.path(DIR_OUT, "sessionInfo.txt"))
writeLines(c(
    "Run after foodwas / cohort and score.R. No dependence on global dat or legacy results.",
    paste0("Exposure: dose (compatibility alias diet.dose.sum); scoring method: ", dose_meta$method,
           "; theoretical range ", paste(dose_meta$range, collapse="-"), "; step ", dose_meta$step, "."),
    "Scoring criteria are defined in score.R and stored in dose_score.rds; no alternative scores are analysed here.",
    "Full prepared-cohort quantile cutoffs and SD are fixed across outcomes; Q1 reference; within-outcome median-score trend.",
    paste0("Observed DOSE quartile ranges: ", paste(paste0(group_def$score_min, "-", group_def$score_max), collapse="; "), "."),
    "Cohort inherits all s0 exclusions (including prevalent overall stroke), then excludes prevalent cases per outcome.",
    "Prostate cancer (C61) is analysed among men only (sex.b = Male); the now-constant sex covariate is omitted and recorded. Other outcomes retain both sexes.",
    "ICD regex and labels are taken unchanged from icd10_other.lst, excluding its duplicate stroke entry.",
    "ICD dates replace same-named cohort dates, especially broad cvd; baseline adjustment covariates remain unchanged.",
    "Matching ICD records without dates are excluded from incident analysis as timing-unknown.",
    paste0("Censoring: earliest loss, death or administrative date ", end_date, "; events must be after baseline and on/before censoring."),
    "Stroke follow-up is rebuilt too; inspect stroke_followup_comparison_QC.csv for differences from the existing helper.",
    paste0("M1 only, matching s3_2 for DOSE. Alcohol in score: ", dose_meta$alcohol, "; separate alcohol adjustment is omitted when true. No regional strata or PCs added."),
    paste0("Requested M1 covariates: ", paste(m1_covs, collapse=", "), "."),
    "Complete cases are selected using only M1 covariates and DOSE within each outcome; continuous, quartile and trend models share that cohort.",
    "Constant covariates after exclusions are omitted and recorded; other singularities/convergence failures yield missing estimates.",
    "FDR families are across all planned outcomes within each model/scale or model/category; tests adjusted separately.",
    "Tables excluding stroke retain the same full outcome-wide FDR family; no recalculation on selected subsets.",
    "Q1 is a reference convention; consult Events, Group_N and model_fit_status before interpreting sparse outcomes.",
    "Outcome-wide associations are exploratory. The current cohort was also used for DOSE development.",
    "MAKE_FIGURES controls circular bars and dot-and-whisker plots in PDF and PNG for both scales. Circular plots use the Fig2 MYTH palette with stroke kept coral red; dotplots retain the reverse reference palette. Use --circular-only to refresh circular plots from saved results. All analysis and plotting code is contained in score.R."
), file.path(DIR_OUT, "README_analysis.txt"))
cat("Primary DOSE outcome-wide analysis complete: ", nrow(outcomes), " outcomes; ",
    paste(MODELS, collapse=", "), "; ", sum(continuous$scale == "Per 1-SD" & continuous$FDR < .05, na.rm=TRUE),
    " outcomes with continuous-score BH FDR < 0.05.\n", sep="")
cat("Outcome-wide outputs saved to: ", DIR_OUT, "\n", sep="")

} # End full analysis; plot-only mode never reads participant-level inputs.
if (PLOT_ONLY) {
    if (!FOREST_ONLY) plot_dose_other_outcomes(DIR_OUT, pal_circular)
    if (!CIRCULAR_ONLY && !FOREST_ONLY) plot_dose_other_dotplot(DIR_OUT, pal)
    if (!CIRCULAR_ONLY) plot_dose_other_chapter_forest(DIR_OUT)
    message("Plot-only refresh complete; analysis tables unchanged.")
}
}))

# ===========================================================================
# score / stroke_comparison
# ===========================================================================
dose_register_step("score", "stroke_comparison", expression({
# Four-pattern UKB incident overall stroke, M2. NHANES uses saved survey ORs.
library(survival)
ROOT <- dose_output_dir()
OUT <- file.path(ROOT,'ukb/04score/overall_stroke_comparison')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
cfg <- readRDS(file.path(ROOT,'ukb/00data/analysis_config.rds'))
dat <- as.data.frame(readRDS(file.path(ROOT,'ukb/00data/dat.base.imp.rds')))
score <- readRDS(file.path(ROOT,'ukb/04score/dose_score.rds'))
stopifnot(!anyDuplicated(dat$eid),!anyDuplicated(score$eid),setequal(dat$eid,score$eid))
dat$dose <- score$dose[match(dat$eid,score$eid)]
vars <- c(dose='dose',medi24='diet.medi24.sum',mind='diet.mind.sum',dash='diet.dash.sum')
labels <- c(dose='DOSE',medi24='MedDiet',mind='MIND',dash='DASH')
covs <- lapply(names(vars),function(v){
 co <- unique(c(cfg$covs.m1,cfg$covs.m2.add,if(v=='dash')cfg$covs.alcohol))
 if(v!='dash')co <- setdiff(co,cfg$covs.alcohol)
 co
});names(covs)<-names(vars)
stopifnot(all(c(unname(vars),unlist(covs),'fod_icd10_qi_stroke','start_date','date_death','date_lost') %in% names(dat)))
# Same date-censoring rule as current score.R.
end <- pmin(dat$date_lost,dat$date_death,as.Date(cfg$date_follow_end),na.rm=TRUE)
onset <- dat$fod_icd10_qi_stroke
prior <- !is.na(onset)&!is.na(dat$start_date)&onset<=dat$start_date
event <- as.integer(!is.na(onset)&onset>dat$start_date&onset<=end)
finish <- end;finish[which(event==1)]<-onset[which(event==1)]
time <- as.numeric(finish-dat$start_date)/365.25
needed<-unique(c(unname(vars),unlist(covs)))
ok <- !prior & is.finite(time)&time>0&!is.na(event)&complete.cases(dat[needed])
for(v in needed)if(is.numeric(dat[[v]]))ok<-ok & is.finite(dat[[v]])
d <- dat[ok,needed];d$time<-time[ok];d$event<-event[ok]
message('UKB common M2 cohort: ',nrow(d),'; events: ',sum(d$event))
effect<-function(f,term,increment=1){
 b<-unname(coef(f)[term]);se<-sqrt(vcov(f)[term,term]);stopifnot(is.finite(b),is.finite(se),se>0)
 data.frame(Estimate=exp(b*increment),LCI=exp((b-qnorm(.975)*se)*increment),UCI=exp((b+qnorm(.975)*se)*increment),P=2*pnorm(-abs(b/se)))
}
rows<-list();breaks<-list();layout_rows<-list()
for(v in names(vars)){
 message('Fitting ',v)
 xall<-dat[[vars[v]]];br<-unname(quantile(xall,c(0,.25,.5,.75,1),na.rm=TRUE))
 stopifnot(!anyDuplicated(br))
 d$x<-d[[vars[v]]];d$group<-cut(d$x,br,include.lowest=TRUE,labels=paste0('Q',1:4))
 stopifnot(!anyNA(d$group),all(table(d$group)>0))
 fg<-coxph(reformulate(c('group',covs[[v]]),'Surv(time,event)'),d,ties='efron',singular.ok=FALSE)
 fc<-coxph(reformulate(c('x',covs[[v]]),'Surv(time,event)'),d,ties='efron',singular.ok=FALSE)
 sdscore<-sd(xall,na.rm=TRUE)
 if(v %in% c('dose','dash')) {
  med<-tapply(d$x,d$group,median)
  d$trend<-unname(med[as.character(d$group)])
  ft<-coxph(reformulate(c('trend',covs[[v]]),'Surv(time,event)'),d,ties='efron',singular.ok=FALSE)
  step<-diff(unname(quantile(xall,c(.1,.9),na.rm=TRUE)))
  layout_rows[[v]]<-cbind(data.frame(Diet=v,Increment=step,P_trend=effect(ft,'trend')$P),effect(fc,'x',step))
 }
 for(ct in c(paste0('Q',1:4),'Per SD')){
  e<-if(ct=='Q1')data.frame(Estimate=1,LCI=1,UCI=1,P=NA_real_) else if(ct=='Per SD')effect(fc,'x',sdscore) else effect(fg,paste0('group',ct))
  i<-if(ct=='Per SD')rep(TRUE,nrow(d)) else d$group==ct
  rows[[paste(v,ct)]]<-cbind(data.frame(Cohort='UK Biobank',Diet=v,Label=labels[v],Model='M2',Outcome='Incident overall stroke',Effect='HR',Contrast=ct,N=sum(i),Events=sum(d$event[i]),Increment=if(ct=='Per SD')sdscore else NA_real_),e)
 }
 breaks[[v]]<-data.frame(Diet=v,Q0=br[1],Q25=br[2],Q50=br[3],Q75=br[4],Q100=br[5],SD=sdscore,Covariates=paste(covs[[v]],collapse=' + '))
}
res<-do.call(rbind,rows);rownames(res)<-NULL
res$FDR<-NA_real_
for(ct in c('Q2','Q3','Q4','Per SD')){i<-which(res$Contrast==ct);res$FDR[i]<-p.adjust(res$P[i],'BH')}
write.csv(res,file.path(OUT,'UKB_overall_stroke_M2.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,breaks),file.path(OUT,'UKB_model_definitions.csv'),row.names=FALSE)
write.csv(do.call(rbind,layout_rows),file.path(OUT,'UKB_DOSE_DASH_P10P90_trend.csv'),row.names=FALSE)
nh<-readRDS(file.path(ROOT,'nhanes/NHANES_analysis.rds'))
nz<-nh$cohorts$two_day_2003_2018
wq<-function(x,w,p){o<-order(x);xx<-x[o];cw<-cumsum(w[o])/sum(w);vapply(p,function(q)xx[which(cw>=q)[1]],numeric(1))}
nh_rows<-lapply(c('dose','dash'),function(v){
 a<-nh$audit$two_day_2003_2018[[v]]
 b<-a$cutpoints
 qq<-wq(nz[[v]],nz$weight,c(.1,.9))
 data.frame(Diet=v,Min=b[1],T33=b[2],T67=b[3],Max=b[4],P10=qq[1],P90=qq[2],Increment=diff(qq),SD=a$sd)
})
write.csv(do.call(rbind,nh_rows),file.path(OUT,'NHANES_DOSE_DASH_ranges.csv'),row.names=FALSE)
write.csv(data.frame(N=nrow(d),Events=sum(d$event),person_years=sum(d$time),Prior_excluded=sum(prior),changed_time=sum(abs(time[ok]-dat$qi_stroke.t2e[ok])>1e-8,na.rm=TRUE)),file.path(OUT,'UKB_cohort_QC.csv'),row.names=FALSE)
writeLines(c('Overall-stroke onset is fod_icd10_qi_stroke. Follow-up ends at earliest death, loss to follow-up or administrative end.',
'All four patterns share the complete M2 cohort. Quartile boundaries and SD use the full input cohort, following s3_2 conventions. Ties are retained.',
'FDR is BH across the four patterns within each contrast. The figure reports raw P values.',
'NHANES: existing two-day mean 2003-2018, survey-weighted prevalent any-stroke ORs, weighted tertiles and weighted SD.',
'UKB MedDiet (medi24) and NHANES aMED differ in implementation. NHANES MIND is adapted. No pooled cross-cohort effect is estimated.'),file.path(OUT,'README.txt'))
message('Saved UKB overall stroke results')
}))

if (sys.nframe() == 0L) dose_run_module("score")
