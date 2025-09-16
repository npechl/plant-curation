

# load libraries --------------------

library(data.table)
library(stringr)
library(seqinr)

# input data ----------------------

df <- fread("download", sep = "\t", header = TRUE, na.strings = c("", "NA"), fill = Inf)

# filtering out --------------------------

df <- df |> subset(!is.na(nuc))

df$nuc <- df$nuc |> str_remove_all("\\-")

# split per marker ------------------

its0 <- df |> subset(marker_code %in% c("ITS", "ITS1", "ITS2"))
trnl <- df |> subset(marker_code == "trnL")

# Make headers -----------------------

its0[, tax_col := paste(kingdom, phylum, class, order, family, genus, species, sep = ";")]
trnl[, tax_col := paste(kingdom, phylum, class, order, family, genus, species, sep = ";")]

# workflow -----------------------

dir.create("./r-curation", showWarnings = FALSE)

write.fasta(sequences = as.list(its0$nuc), names = as.list(its0$processid), file.out = "r-curation/ITS0.fasta", nbchar = 120)
write.fasta(sequences = as.list(trnl$nuc), names = as.list(trnl$processid), file.out = "r-curation/trnL.fasta", nbchar = 120)

R.utils::gzip("./r-curation/ITS0.fasta")
R.utils::gzip("./r-curation/trnL.fasta")

its0_tax <- its0[, .(processid, tax_col)]
trnl_tax <- trnl[, .(processid, tax_col)]

fwrite(its0_tax, "./r-curation/ITS0.tax", sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)
fwrite(trnl_tax, "./r-curation/trnL.tax", sep = "\t", row.names = FALSE, col.names = F, quote = FALSE)








