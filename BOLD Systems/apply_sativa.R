library(data.table)

tax_file <- "test.tax"
mis_file <- "syntest.mis"
out_file <- "test_corrected.tax"

tax <- fread(tax_file, header = F, sep = "\t")
mis <- fread(mis_file, header = T, sep = "\t")

setnames(tax, c("ID", "Taxonomy"))
setnames(mis, old = ";SeqID", new = "ID")

oldtax <- copy(tax)
tax[mis, on = "ID", Taxonomy := ProposedTaxonomyPath]

oldtax[Taxonomy != tax$Taxonomy]

fwrite(tax, out_file, sep = "\t", quote = FALSE)
