# Stream baseline assays column by column into a private temporary file.
args <- commandArgs(TRUE)
x <- readRDS(args[1])
stopifnot(is.data.frame(x), "eid" %in% names(x), !anyDuplicated(names(x)),
          !anyNA(x$eid), !anyDuplicated(x$eid))
fields <- setdiff(names(x), "eid")
technical <- grepl("(^|[. _:-])(plate|well|centre|center)([. _:-]|$)", fields, ignore.case=TRUE)
fields <- fields[!technical]
if (!length(fields) || !all(vapply(x[fields], is.numeric, logical(1)))) stop("Expected named numeric baseline assays")
if (any(tolower(fields) %in% c("eid","age","sex","event","time","target","split","y","group","role","date_attend")))
    stop("An outcome or participant attribute was found among the assays")
writeLines(as.character(x$eid), file.path(args[2],"ids.txt"))
writeLines(fields, file.path(args[2],"features.txt"))
con <- file(file.path(args[2],"baseline.f32"),"wb")
for (nm in fields) writeBin(as.numeric(x[[nm]]),con,size=4L,endian="little")
close(con)
