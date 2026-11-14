# Frozen DOSE scoring recipe from source 971a61b44b84500ab46e2f5592478505a0eae50e.
# Pure base-R functions: no file I/O, cohort filtering or disease analysis.

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
