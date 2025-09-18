
library(seqinr)
library(stringr)

concatITS <- function(its1, s58, its2) {
    
    s1 <- its1 |> read.fasta()
    s2 <- s58 |> read.fasta()
    s3 <- its2 |> read.fasta()
    
    names(s1) = s1 |> names() |> str_split_i("\\|", 1)
    names(s2) = s2 |> names() |> str_split_i("\\|", 1)
    names(s3) = s3 |> names() |> str_split_i("\\|", 1)
    
    overlap <- intersect(names(s1), names(s2))
    overlap <- intersect(overlap, names(s3))
    
    s1 <- s1[overlap]
    s2 <- s2[overlap]
    s3 <- s3[overlap]
    
    concat <- mapply(function(x, y, z) {
        
        o <- c(x, y, z)
        
    }, s1, s2, s3)
    
    return(concat)
    
}

its1 = "test/250902.ITS1.fasta"
s58  = "test/250902.5_8S.fasta"
its2 = "test/250902.ITS2.fasta"

outdir <- its1 |> dirname()
db <- its1 |> basename() |> str_split_i("\\.", 1)

full_its <- concatITS(its1, s58, its2)

write.fasta(
    full_its, names(full_its) |> as.list(), nbchar = 120, 
    file.out = paste0(outdir, "/", db, ".full.ITS.fasta")
)

