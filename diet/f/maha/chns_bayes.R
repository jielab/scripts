# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
maha_load_state("assoc")
pub_precision <- make_precision_table_hr(mort)
for (filename in c("Fig6.mortality_validation.out.xlsx", "publication_inputs.out.xlsx")) {
    path <- cohort_file(filename)
    if (!file.exists(path)) stop("Run the CHNS assoc phase first: ", path)
    sheets <- readxl::excel_sheets(path)
    book <- setNames(lapply(sheets, function(s) readxl::read_excel(path, sheet = s)), sheets)
    book$precision <- pub_precision
    writexl::write_xlsx(book, path)
}
write_xlsx(list(precision = pub_precision), "Bayes.out.xlsx")
maha_write_bayes(pub_precision, "chns")
