rm(list = ls())
gc()

library(data.table)
library(Biostrings)

# --- Αρχεία ---
fasta_file <- "PLANiTS_v2020.3-full/ITS1.fasta"
taxonomy_file <- "PLANiTS_v2020.3-full/ITS1_taxonomy"

# --- Διαβάζουμε taxonomy ---
tax <- fread(taxonomy_file, header = FALSE, sep = "\t", na.strings = "NA")
setnames(tax, c("id", "taxonomy"))
# Αντικαθιστούμε NA με κενό
tax[, taxonomy := gsub("NA", "", taxonomy)]

# --- Διαβάζουμε fasta ---
fasta <- readDNAStringSet(fasta_file)
ids <- names(fasta)

# --- Merge ---
dt <- data.table(id = ids, seq = as.character(fasta))
dt <- merge(dt, tax, by = "id", all.x = TRUE)

# --- Συνάρτηση για split σε 80 bp ---
wrap80 <- function(s) {
  paste0(strwrap(s, width = 80), collapse = "\n")
}

# --- Taxonomy FASTA ---
tax_headers <- paste0(">", dt$taxonomy, ";")
tax_seqs <- vapply(dt$seq, wrap80, "")
writeLines(c(rbind(tax_headers, tax_seqs)), "toSpecies.fasta")

# --- Process FASTA ---
species <- vapply(strsplit(dt$taxonomy, ";"), function(x) tail(x, 1), "")
species <- sapply(species, function(s) {
  # αν είναι μόνο μία λέξη με πρώτο κεφαλαίο (πιθανό γένος)
  if (grepl("^[A-Z][a-z]+$", s)) s <- paste("unknown", s)
  s
})
proc_headers <- paste0(">", dt$id, ".ITS1 ", species)
proc_seqs <- vapply(dt$seq, wrap80, "")
writeLines(c(rbind(proc_headers, proc_seqs)), "assignSpecies.fasta")

