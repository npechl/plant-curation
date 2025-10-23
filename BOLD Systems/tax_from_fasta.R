# load libraries --------------------
library(data.table)
library(stringr)
library(seqinr)

# read fasta ------------------------------------------------
fasta <- read.fasta("test.fasta", as.string = T,whole.header = T)
headers <- names(fasta)

# parse headers --------------------------------------------------------------
ids <- str_extract(headers, "^[^.]+")
taxa <- str_replace(headers, "^[^.]+\\.[^_]*_", "")

# combine and write ------------------------------------
out <- data.table(id = ids, taxonomy = taxa)
fwrite(out, "test.tax", sep = "\t", col.names = FALSE)

