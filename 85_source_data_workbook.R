## 85_source_data_workbook.R
## The source data as a single Excel workbook, because the submission system unpacks a .zip into one item per file.
## Sheet "README" holds 05_source_data/README.txt, sheet "Index" maps each sheet (T01, T02, ...) to its source file,
## and each table follows on its own sheet with its values unchanged. Output: Source_data.xlsx (next to the zip).
suppressPackageStartupMessages({ library(data.table); library(openxlsx) })
S <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission"
SD <- file.path(S, "05_source_data"); OUT <- file.path(S, "00_UPLOAD", "Source_data.xlsx")
files <- sort(setdiff(list.files(SD), "README.txt"))
wb <- createWorkbook()
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(README = readLines(file.path(SD, "README.txt"), encoding = "UTF-8")))
idx <- data.table(sheet = sprintf("T%02d", seq_along(files)), file = files, rows = NA_integer_, columns = NA_integer_)
addWorksheet(wb, "Index")
for (i in seq_along(files)) {
  d <- fread(file.path(SD, files[i]), encoding = "UTF-8", na.strings = c("", "NA"))
  idx[i, `:=`(rows = nrow(d), columns = ncol(d))]
  addWorksheet(wb, idx$sheet[i]); writeData(wb, idx$sheet[i], d) }
writeData(wb, "Index", idx)
setColWidths(wb, "Index", cols = 1:4, widths = c(8, 80, 8, 8)); setColWidths(wb, "README", cols = 1, widths = 120)
saveWorkbook(wb, OUT, overwrite = TRUE)
cat(sprintf("workbook written: %s (%d tables, %.1f MB)\n", OUT, length(files), file.size(OUT) / 1e6))
