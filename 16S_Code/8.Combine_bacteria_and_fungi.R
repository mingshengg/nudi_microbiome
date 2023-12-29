##  ###################################################  ##
##  Combine Bacterial and fungal ps objects

##  Author: Ming Sheng                                   ##
##                                                       ##

# load packages
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(ShortRead); packageVersion("ShortRead")

# load data objects
fung <- readRDS("./Output/fungal_noncontam_ps_object.RDS")
fung <- fung %>% subset_samples(Species != 'Trinchesia sibogae')
bact <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")
bact <- bact %>% subset_samples(Species != 'Trinchesia sibogae')
bact@sam_data[bact@sam_data$Species == 'Pteraeolidia semperi complex',]$Species <- "Pteraeolidia semperi"


# pull "Structure" and "Location" from SampleID in fungal objects
fung_meta <- fung %>% microbiome::meta()
bact_meta <- bact %>% microbiome::meta()
bact_meta$Site %>% unique()
bact_meta$Species %>% unique()

fung_meta <- fung_meta %>% 
  mutate(Site = case_when(str_detect(SampleID, pattern = "CY") ~ "Cyrene",
                              str_detect(SampleID, pattern = "TS|HW|HN|HG") ~ "Pulau Hantu",
                              str_detect(SampleID, pattern = "JG") ~ "Jong",
                              str_detect(SampleID, pattern = "SS") ~ "Semakau",
                              str_detect(SampleID, pattern = "TK") ~ "Tekukor")) %>% 
  mutate(Study = "Fungi")

# rename fields in bact
bact@sam_data$Study = "Bacteria"

# reassign metadata
fung@sam_data <- sample_data(fung_meta)

# clean fungal ps objects
fung <- fung %>% 
  subset_taxa(Kingdom == "k__Fungi") %>% 
  subset_samples(sample_sums(fung) > 0)

bact <- bact %>% 
  subset_taxa(taxa_sums(bact) > 0) %>% 
  subset_samples(sample_sums(bact) > 0)

# rename fungal taxa
tax_table(fung)[,1] <- tax_table(fung)[,1] %>% str_remove_all("k__")
tax_table(fung)[,2] <- tax_table(fung)[,2] %>% str_remove_all("p__")
tax_table(fung)[,3] <- tax_table(fung)[,3] %>% str_remove_all("c__")
tax_table(fung)[,4] <- tax_table(fung)[,4] %>% str_remove_all("o__")
tax_table(fung)[,5] <- tax_table(fung)[,5] %>% str_remove_all("f__")
tax_table(fung)[,6] <- tax_table(fung)[,6] %>% str_remove_all("g__")
tax_table(fung)[,7] <- tax_table(fung)[,7] %>% str_remove_all("s__")

tax_table(bact)[,1] %>% table
tax_table(fung)[,1] %>% table

# inspect
fung
bact

# merge_taxa at species level
fung <- tax_glom(fung,taxrank = "Species")
bact <- tax_glom(bact,taxrank = "Genus")

# rename samples
fung <- fung %>% 
  subset_samples(!duplicated(fung@sam_data$SampleID))

sample_names(fung) <- fung@sam_data$SampleID

fung <- fung %>% subset_samples(sample_names(fung) %in% sample_names(bact))


# combine and clean up environment
fung@sam_data <- fung@sam_data[,-3]
bact@sam_data <- bact@sam_data[,-3]

full2 <- merge_phyloseq(fung,otu_table(bact),tax_table(bact),sample_data(bact))
tax_table(full2)[,1] %>% 
  table()


sample_names(bact)
# rm(fung_2019);rm(fung_2020);rm(bact)

# add microbe designation (bact vs fung)
full2@sam_data$Microbe <- ifelse(full2@sam_data$Study == "Bacteria",yes = "Bacteria",no="Fungi")

# save full ps object
saveRDS(full2,"./Output/bact_and_fungi_clean_ps_object.RDS")
