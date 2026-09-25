## G1_gnhs_describe.R
## First look at the published supplementary tables of the Guangzhou Nutrition and Health Study serum proteomics paper
## (Cai et al., Cell Reports Medicine 2023;4:101172), to see whether they can serve as an external test of the question
## this manuscript asks: at the same measured body size, does a proteomic read-out of central adiposity separate risk?
## Reads only; writes a description to the GNHS folder.
suppressPackageStartupMessages({ library(data.table); library(readxl) })
G <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920/GNHS"; TB <- file.path(G, "derived")
dir.create(TB, showWarnings = FALSE)
say <- function(...) { cat(sprintf(...), "\n"); flush.console() }
rd <- function(n, sheet) as.data.table(read_excel(file.path(G, "tables", sprintf("tableS%d.xlsx", n)), sheet = sheet, na = c("", "NA", "NaN")))

## ---------------- clinical table ----------------
cl <- rd(1, "Clinical information"); say("S1 clinical: %d participants, %d columns", nrow(cl), ncol(cl))
print(names(cl))
for (v in names(cl)) { x <- cl[[v]]; u <- unique(x[!is.na(x)])
  say("  %-16s class %-9s missing %5d  unique %5d  %s", v, class(x)[1], sum(is.na(x)), length(u),
      if (length(u) <= 6) paste(sort(u), collapse = ", ") else sprintf("range %.2f to %.2f", min(as.numeric(x), na.rm = TRUE), max(as.numeric(x), na.rm = TRUE))) }

## ---------------- proteomics tables ----------------
for (n in c(2, 3)) { a <- rd(n, "A"); b <- rd(n, "B")
  say("S%d meta: %d rows, columns %s", n, nrow(a), paste(names(a), collapse = ", "))
  say("S%d matrix: %d samples x %d proteins", n, nrow(b), ncol(b) - 1)
  tp <- names(a)[2]; print(table(a[[tp]], useNA = "ifany"))
  say("S%d unique participants %d; samples in matrix also in meta: %d", n, uniqueN(a$Patient_ID), sum(b[[1]] %in% a$Sample_ID))
  fwrite(a, file.path(TB, sprintf("S%d_meta.csv", n))); fwrite(b, file.path(TB, sprintf("S%d_matrix.csv", n))) }

## ---------------- how the sample ids relate to the participant ids ----------------
a2 <- fread(file.path(TB, "S2_meta.csv")); b2 <- fread(file.path(TB, "S2_matrix.csv"))
say("S2 sample ids in the matrix, first 5: %s", paste(head(b2[[1]], 5), collapse = ", "))
say("S2 sample ids of the meta, first 5: %s", paste(head(a2$Sample_ID, 5), collapse = ", "))
say("clinical ids, first 5: %s", paste(head(cl$Patient_ID, 5), collapse = ", "))
say("clinical ids matched by S2 participants: %d of %d", sum(unique(a2$Patient_ID) %in% cl$Patient_ID), uniqueN(a2$Patient_ID))
fwrite(cl, file.path(TB, "S1_clinical.csv")); say("DONE")
