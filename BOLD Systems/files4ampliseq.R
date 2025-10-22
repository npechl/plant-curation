# Load libraries
library(data.table)
library(tidyverse)
library(stringr)
library(Biostrings)

#------------------INPUT FILES----------------
fasta_file <- "test.fasta.gz"
tax_file   <- "test.tax"

#---------------READ TAX FILE-------------------------
tax <- fread(tax_file, header = FALSE, sep = "\t", 
             col.names = c("ID", "taxonomy"))

#-------------READ FASTA---------------------------------
fasta <- readDNAStringSet(fasta_file)
fasta_dt <- data.table(
  header = names(fasta),
  sequence = as.character(fasta)
)

# Extract the ID (before the first dot)
fasta_dt[, ID := str_extract(header, "^[^.]+")]

# Merge with taxonomy info
merged <- merge(fasta_dt, tax, by = "ID", all.x = TRUE)

#-----------------FASTA 1: Taxonomy as header------------------------------
fasta_taxonomy <- merged %>%
  select(header = taxonomy, sequence)

# Write fasta
write_lines(
  paste0(">", fasta_taxonomy$header, "\n", fasta_taxonomy$sequence),
  "taxonomy_headers.fasta"
)
R.utils::gzip("taxonomy_headers.fasta")
#-------------------FASTA 2: Sample ID + Species---------------------
# Extract last taxon name (species)
merged <- merged %>%
  mutate(species = word(taxonomy, -1, sep = fixed(";"))) %>%
  filter(species != "NA") # remove entries with no species

fasta_species <- merged %>%
  mutate(header = paste(ID, species)) %>%
  select(header, sequence)

# Write fasta
write_lines(
  paste0(">", fasta_species$header, "\n", fasta_species$sequence),
  "id_species.fasta"
)
R.utils::gzip("id_species.fasta")
