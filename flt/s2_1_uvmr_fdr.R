# Recalculate BH FDR for IVW only, independently within each outcome sheet.
# FDR规则 共同 182 个 FLT 中，按实际有效 IVW 数进行 BH 校正

library(openxlsx); library(dplyr)

dir0 <- "D:"
project <- file.path(dir0, "projects/flt")
input_file  <- file.path(project, "analysis/result.xlsx")
output_file <- file.path(project, "analysis/result_IVW_BH.xlsx")

category_map <- read.xlsx(file.path(dir0, "files/flt.xlsx"), sheet = "item182", check.names = FALSE) %>%
  select(gwas_id, category) %>% mutate(across(everything(), ~ trimws(as.character(.x))))
message(sprintf("Classification: %d GWAS IDs in %d categories", nrow(category_map), n_distinct(category_map$category)))

sheet_config <- data.frame(
  sheet   = c("Sheet2", "Sheet3", "Sheet4", "Sheet5", "Sheet6"),
  outcome = c("dementia_r13", "anystroke_eur", "sleep_disorder_r13", "mdd_2025", "anxiety_2026"),
  stringsAsFactors = FALSE
)

wb <- loadWorkbook(input_file)
for (i in seq_len(nrow(sheet_config))) {
  sheet   <- sheet_config$sheet[i]
  outcome <- sheet_config$outcome[i]
  dat <- read.xlsx(wb, sheet, check.names = FALSE, skipEmptyRows = FALSE, skipEmptyCols = FALSE)

  stopifnot(all(c("Exposure", "Outcome", "Method", "P-value") %in% names(dat)))

  ivw      <- which(toupper(trimws(as.character(dat$Method))) == "IVW")
  exposure <- trimws(as.character(dat$Exposure[ivw]))

  if (anyNA(exposure) || any(exposure == "") || anyDuplicated(exposure)) stop(sheet, ": IVW Exposure must be present and unique")
  if (!isTRUE(all(dat$Outcome[ivw] == outcome))) stop(sheet, ": unexpected IVW Outcome")

  category <- category_map$category[match(exposure, category_map$gwas_id)]
  if (anyNA(category)) stop(sheet, ": unmapped IVW Exposure: ", paste(exposure[is.na(category)], collapse = ", "))

  p     <- suppressWarnings(as.numeric(as.character(dat[["P-value"]][ivw])))
  valid <- is.finite(p) & p >= 0 & p <= 1
  if (any(!valid)) warning(sheet, ": excluded ", sum(!valid), " invalid IVW P-values")

  # 家族 1：本 outcome 全部有效 IVW，n = 有效总数（常数）
  p_ivw     <- rep(NA_real_,    nrow(dat))
  n_ivw     <- rep(NA_integer_, nrow(dat))
  p_ivw[ivw[valid]] <- p.adjust(p[valid], method = "BH")
  n_ivw[ivw[valid]] <- sum(valid)

  # 家族 2：outcome × category，n = 该分类内有效数（随行变化）
  p_cat <- rep(NA_real_,    nrow(dat))
  n_cat <- rep(NA_integer_, nrow(dat))
  for (group in unique(category)) {
    idx <- which(valid & category == group)
    p_cat[ivw[idx]] <- p.adjust(p[idx], method = "BH")
    n_cat[ivw[idx]] <- length(idx)
    message(sprintf("%s | %s: %d valid IVW tests", sheet, group, length(idx)))
  }

  # 若列已存在则覆盖，否则追加；写回 wb 并调宽
  put <- function(name, values, width) {
    col <- match(name, names(dat))
    if (is.na(col)) col <- ncol(dat) + 1L
    writeData(wb, sheet, setNames(data.frame(values), name), startCol = col, startRow = 1, colNames = TRUE, rowNames = FALSE, keepNA = FALSE)
    setColWidths(wb, sheet, cols = col, widths = width)
    dat[[name]] <<- values
  }
  put("P_FDR_IVW",      p_ivw, 18)
  put("N_IVW",          n_ivw, 12)
  put("P_FDR_category", p_cat, 18)
  put("N_category",     n_cat, 12)

  message(sprintf("%s (%s): %d valid IVW tests; %d category-BH P < 0.05",
                  sheet, outcome, sum(valid), sum(p_cat < 0.05, na.rm = TRUE)))
}

# 保留源工作簿及其原有 P-FDR 列，另存为新文件
saveWorkbook(wb, output_file, overwrite = TRUE)
message("Saved: ", output_file)