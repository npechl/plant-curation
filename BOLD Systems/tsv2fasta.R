

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

its  <- df |> subset(marker_code == "ITS")
its1 <- df |> subset(marker_code == "ITS1")
its2 <- df |> subset(marker_code == "ITS2")
trnl <- df |> subset(marker_code == "trnL")

# Make headers -----------------------

its[, tax_header := paste0(record_id, "_", paste(kingdom, phylum, class, order, family, genus, species, sep = ";"))]
its1[, tax_header := paste0(record_id, "_", paste(kingdom, phylum, class, order, family, genus, species, sep = ";"))]
its2[, tax_header := paste0(record_id, "_", paste(kingdom, phylum, class, order, family, genus, species, sep = ";"))]
trnl[, tax_header := paste0(record_id, "_", paste(kingdom, phylum, class, order, family, genus, species, sep = ";"))]

# workflow -----------------------

dir.create("./r-curation", showWarnings = FALSE)

write.fasta(sequences = as.list(its$nuc), names = as.list(its$tax_header), file.out = "r-curation/ITS.fasta", nbchar = 120)
write.fasta(sequences = as.list(its1$nuc), names = as.list(its1$tax_header), file.out = "r-curation/ITS1.fasta", nbchar = 120)
write.fasta(sequences = as.list(its2$nuc), names = as.list(its2$tax_header), file.out = "r-curation/ITS2.fasta", nbchar = 120)
write.fasta(sequences = as.list(trnl$nuc), names = as.list(trnl$tax_header), file.out = "r-curation/trnL.fasta", nbchar = 120)

fwrite(its, "./r-curation/ITS.tax", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(its1, "./r-curation/ITS1.tax", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(its2, "./r-curation/ITS2.tax", sep = "\t", row.names = FALSE, quote = FALSE)
fwrite(trnl, "./r-curation/trnL.tax", sep = "\t", row.names = FALSE, quote = FALSE)








