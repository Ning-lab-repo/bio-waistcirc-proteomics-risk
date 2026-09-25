## 86_protein_names_all.R
## Full names of all 2,920 proteins, for the source data. Each protein takes the recommended name of its UniProt
## entry (human, reviewed), matched on the gene symbol (primary symbol first, then synonyms); an assay that measures
## two or more genes takes the names of all of them. The 28 proteins named in the manuscript, in Supplementary Tables 3
## and 5 and in the figures keep the curated names written below, which add a common alternative name where it helps.
## The blank name columns of the two 2,920-protein source tables are then filled from the same list.
## Output: T119_protein_names_all.csv; in the source data, source_table_Methods_protein_names.csv (all 2,920 proteins)
##         and the name columns of source_table_Figure3C_protein_axis.csv and
##         source_table_Results_protein_associations_with_fat_depots.csv (other columns unchanged).
suppressPackageStartupMessages(library(data.table))
W  <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260916ProWC/20260920"; TB <- file.path(W, "tables")
SD <- "/home/data/heamei/nmrLR/yuxuan_pca/yuxuan_wc_fig/buchongshuju/20260920_eBioMedicine_submission/05_source_data"
prot <- fread(file.path(SD, "source_table_Methods_imputation_medians.csv"))$protein   ## the 2,920 analysed proteins
stopifnot(length(prot) == 2920, !anyDuplicated(prot))

## UniProt, human reviewed entries (downloaded once and kept with the tables)
up_file <- file.path(TB, "uniprot_human_reviewed.tsv.gz")
if (!file.exists(up_file))
  download.file(paste0("https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cgene_primary%2Cgene_synonym%2Cprotein_name",
                       "&format=tsv&query=%28organism_id%3A9606%29%20AND%20%28reviewed%3Atrue%29"), up_file, mode = "wb", quiet = TRUE)
up <- fread(cmd = paste("zcat", shQuote(up_file)), sep = "\t", quote = "", header = TRUE, colClasses = "character")
setnames(up, c("accession", "gene", "synonyms", "uniprot_name"))
## recommended name: the text before the first alternative name in parentheses and before any cleavage products
up[, name := sub(" [(].*$", "", sub(" [[](Cleaved into|Includes)[^]]*[]].*$", "", uniprot_name))]
## an entry encoded by two or more genes lists them separated by "; ", in the primary and in the synonym field
pri <- up[, .(s = trimws(unlist(strsplit(gene, ";", fixed = TRUE)))), by = .(accession, name)][s != ""]
syn <- up[, .(s = unlist(strsplit(synonyms, "[; ]+"))), by = .(accession, name)][s != ""]
look1 <- function(g) {
  r <- pri[s == g]
  if (nrow(r) == 0) r <- syn[s == g]
  if (nrow(r) == 0) return(NULL)
  r[order(accession)][1, .(name, accession)]
}
res <- rbindlist(lapply(prot, function(p) {
  hits <- rbindlist(lapply(strsplit(p, "_", fixed = TRUE)[[1]], look1))
  if (nrow(hits) == 0) return(data.table(gene_symbol = p, protein_name = NA_character_, uniprot_accession = NA_character_))
  data.table(gene_symbol = p, protein_name = paste(unique(hits$name), collapse = "; "), uniprot_accession = paste(hits$accession, collapse = ";"))
}))

## curated names of the proteins named in the manuscript, its tables and figures
cur <- fread(text = paste(c("gene_symbol|protein_name",
  "ADM|Adrenomedullin", "CA14|Carbonic anhydrase 14", "CDHR2|Cadherin-related family member 2", "CFH|Complement factor H",
  "CHGB|Secretogranin-1 (chromogranin B)", "CKB|Creatine kinase B-type", "CST3|Cystatin C",
  "FABP4|Fatty acid-binding protein, adipocyte (fatty acid-binding protein 4)", "FGF21|Fibroblast growth factor 21", "FURIN|Furin",
  "GH1|Somatotropin (growth hormone 1)", "IGFBP1|Insulin-like growth factor-binding protein 1",
  "IGFBP2|Insulin-like growth factor-binding protein 2", "IGSF3|Immunoglobulin superfamily member 3",
  "IGSF9|Protein turtle homolog A (immunoglobulin superfamily member 9A)", "IL1RN|Interleukin-1 receptor antagonist protein",
  "INHBC|Inhibin beta C chain", "LEP|Leptin", "NCAN|Neurocan core protein", "OPTC|Opticin", "OXT|Oxytocin-neurophysin 1",
  "PON3|Serum paraoxonase/lactonase 3", "PRAP1|Proline-rich acidic protein 1", "RTN4R|Reticulon-4 receptor",
  "SEZ6L|Seizure 6-like protein", "SLITRK1|SLIT and NTRK-like protein 1",
  "SSC4D|Scavenger receptor cysteine-rich domain-containing group B protein",
  "WFIKKN2|WAP, Kazal, immunoglobulin, Kunitz and NTR domain-containing protein 2"), collapse = "\n"), sep = "|", quote = "")
stopifnot(nrow(cur) == 28, all(cur$gene_symbol %in% prot))
res[cur, on = "gene_symbol", protein_name := i.protein_name]
## two Olink assay names that are not UniProt gene symbols
alias <- data.table(gene_symbol = c("NTproBNP", "ANP32C"),
                    protein_name = c("N-terminal pro-B-type natriuretic peptide (NT-proBNP; natriuretic peptides B, NPPB)",
                                     "Acidic leucine-rich nuclear phosphoprotein 32 family member C (putative protein ANP32CP)"),
                    uniprot_accession = c("P16860", "O43423"))
res[alias, on = "gene_symbol", `:=`(protein_name = i.protein_name, uniprot_accession = i.uniprot_accession)]
cat(sprintf("names: %d of %d proteins (%d curated); not found in UniProt: %s\n", sum(!is.na(res$protein_name)), nrow(res),
            nrow(cur), paste(res[is.na(protein_name)]$gene_symbol, collapse = ", ")))
fwrite(res, file.path(TB, "T119_protein_names_all.csv"))
fwrite(res, file.path(SD, "source_table_Methods_protein_names.csv"))

## the name columns of the two 2,920-protein tables, where blank
for (f in c("source_table_Figure3C_protein_axis.csv", "source_table_Results_protein_associations_with_fat_depots.csv")) {
  d <- fread(file.path(SD, f), colClasses = "character")   ## read as text, so that every other value is written back unchanged
  stopifnot(all(c("protein", "name") %in% names(d)), nrow(d) == 2920)
  i <- which(is.na(d$name) | d$name == "")
  set(d, i, "name", res$protein_name[match(d$protein[i], res$gene_symbol)])
  cat(sprintf("%s: %d names filled, %d still blank\n", f, length(i), sum(is.na(d$name) | d$name == "")))
  fwrite(d, file.path(SD, f))
}
cat("DONE\n")
