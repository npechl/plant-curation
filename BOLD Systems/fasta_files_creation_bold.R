rm(list = ls())
gc()

library(data.table)
library(stringr)

df <- fread("download", sep = "\t", header = TRUE, na.strings = c("", "NA"), fill = Inf)

#filter by marker_code (pick ONE and uncomment)
df <- df[marker_code == "ITS1"]                       # single code
# dt <- dt[marker_code %in% c("ITS1","ITS2")]           # multiple codes

# --- Make headers ---
# toSpecies FASTA
df[, tax_header := paste0(">", 
                          kingdom, ";", phylum, ";", class, ";", order, ";", 
                          family, ";", genus, ";", species, ";")]

# assignSpecies FASTA
df[, process_header := paste0(">", processid, ".", marker_code, " ", species)]

# --- Helper: wrap seqs to 80 bp ---
wrap_seq <- function(seq, width=80) {
  seq <- toupper(seq)  # enforce uppercase
  paste(str_wrap(seq, width=width), collapse="\n")
}

# --- Write taxonomy FASTA ---
writeLines(
  df[, paste0(tax_header, "\n", sapply(nuc, wrap_seq))],
  "toSpecies.fasta"
)

# --- Write process FASTA ---
writeLines(
  df[, paste0(process_header, "\n", sapply(nuc, wrap_seq))],
  "assignSpecies.fasta"
)
