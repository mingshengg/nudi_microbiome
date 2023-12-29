library(ShortRead)
library(Biostrings)

full_ps <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")
ps <- readRDS("./Output/noncontam_ps_object.RDS")

full_ps@tax_table %>% as.data.frame %>% filter(Genus == 'Mycoplasma') <- myco.asv
myco.seq <- ps@refseq[myco.asv]

fasta_dir <- file.path(getwd(), "refs")
outfile <- file.path(dirname(fasta_dir), "myco_fasta.fasta")
writeFasta(myco.seq, outfile, mode = 'a')
