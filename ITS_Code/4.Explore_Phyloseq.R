##  ###################################################  ##
##  Investigate and tidy up full phyloseq object         ##
##                                                       ##
##  Author: Ming Sheng - June 16, 2022                   ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  tidyverse v 1.3.1                                    ##
##  phyloseq v 1.40.0                                    ##
##  vegan v 2.6.2                                        ##
##  VennDiagram v 1.7.3                                  ##
##  patchwork v 1.1.1                                    ##
##                                                       ##
##  ###################################################  ##

# Load packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(VennDiagram); packageVersion("VennDiagram")
library(patchwork); packageVersion("patchwork")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585")
# names(pal) <- c("f","l","p","s")

# Load ps object with tree ####
full_ps <- readRDS("./Output/full_cleaned_ps_object_w_tree.RDS")

# agglomerate taxa at genus level ####
full_ps_genus <- tax_glom(full_ps,"Genus")
full_ps_genus <- full_ps_genus %>% subset_taxa(taxa_sums(full_ps_genus)>0)
full_ps_genus <- full_ps_genus %>% subset_samples(sample_sums(full_ps_genus)>0)
full_ps_genus <- full_ps_genus %>% subset_samples(PCR_Negative != "TRUE")
full_ps_genus <- full_ps_genus %>% subset_samples(Location != 'Cyrene Reef')
full_ps_genus <- full_ps_genus %>% subset_samples(SampleID != 'HW01')
full_ps_genus <- full_ps_genus %>% subset_samples(Species == 'Pteraeolidia semperi')
# output genus-level ps_object for convenience
saveRDS(full_ps_genus, "./Output/full_ps_object_w_tree_genus-glom.RDS")

# VennDiagram of overlap for Species ####

A <- full_ps_genus %>% subset_samples(Species == "Pteraeolidia semperi") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Species == "Pteraeolidia semperi"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(Species == "Trinchesia sibogae") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Species == "Trinchesia sibogae"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B))
n12 <- sum(full %in% unique(A) & full %in% unique(B))

venn.plot1 <- draw.pairwise.venn(area1=length((A)),
                               area2=length((B)),
                               n12,
                               category = c("P. semperi","T. sibogae"),
                               fill = pal[1:2])
dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_Species.png")
grid.draw(venn.plot2)
dev.off()

# VennDiagram of overlap by Location (species combined)(Pool Pulau Hantu together) ####

A <- full_ps_genus %>% subset_samples(Location == "Pulau Hantu (North)"|Location == "Pulau Hantu (West)"|Location == "Pulau Hantu") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Location == "Pulau Hantu (North)"|Location == "Pulau Hantu (West)"|Location == "Pulau Hantu"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(Location == "Pulau Jong") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Location == "Pulau Jong"))>0) %>% tax_table() %>% row.names()
C <- full_ps_genus %>% subset_samples(Location == "Pulau Tekukor") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Location == "Pulau Tekukor"))>0) %>% tax_table() %>% row.names()
D <- full_ps_genus %>% subset_samples(Location == "Cyrene Reef") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Location == "Cyrene Reef"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B,C,D))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n14 <- sum(full %in% unique(A) & full %in% unique(D))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n24 <- sum(full %in% unique(B) & full %in% unique(D))
n34 <- sum(full %in% unique(C) & full %in% unique(D))
n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))
n124 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D))
n134 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D))
n234 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D))
n1234 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

venn.plot2 <- draw.quad.venn(area1=length((A)),
                             area2=length((B)),
                             area3=length((C)),
                             area4=length((D)),
                             n12,n13,n14,n23,n24,n34,n123,n124,n134,n234,n1234,
                             category = c("Pulau Hantu","Pulau Jong","Pulau Tekukor","Cyrene Reef"),
                             fill=pal)
dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_Location.png")
grid.draw(venn.plot2)
dev.off()

# VennDiagram of overlap by different parts of Pulau Hantu ####

A <- full_ps_genus %>% subset_samples(Location == "Pulau Hantu (North)") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Location == "Pulau Hantu (North)"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(Location == "Pulau Hantu (West)") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Location == "Pulau Hantu (West)"))>0) %>% tax_table() %>% row.names()
C <- full_ps_genus %>% subset_samples(Location == "Pulau Hantu") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Location == "Pulau Hantu"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B,C))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))

venn.plot2 <- draw.triple.venn(area1=length((A)),
                             area2=length((B)),
                             area3=length((C)),
                             n12,n13,n23,n123,
                             category = c("Pulau Hantu (North)","Pulau Hantu (West)","Pulau Hantu (T. s)"),
                             fill=pal[1:3])
dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_Hantu_Parts.png")
grid.draw(venn.plot2)
dev.off()

# Explore phylogenetic tree ####
full_ps_genus@sam_data$Species %>% unique()
# colored by phylum (blanks removed)
names(full_ps_genus@sam_data) <- c("SampleID","SampleSpecies","Location","PCR_Negative")

full_ps_genus %>% 
  subset_samples(SampleSpecies != "Blank") %>% 
  plot_tree(ladderize="left", color="SampleSpecies") +
  scale_color_viridis_d()

p1 <- full_ps_genus %>% ##
  subset_samples(Location == "Pulau Hantu (North)"|Location == "Pulau Hantu (West)"|Location == "Pulau Hantu") %>% 
  plot_tree(ladderize="left", color="Location") +
  scale_color_manual(values = pal[1:3]) +
  ggtitle("Pulau Hantu")

p2 <- full_ps_genus %>%  ## Pteraeolidia semperi
  subset_samples(SampleSpecies == "Pteraeolidia semperi") %>% 
  plot_tree(ladderize="left", color="SampleSpecies") +
  scale_color_manual(values = pal [2]) + 
  theme(legend.position = "none") +
  ggtitle("Pteraeolidia semperi")

p3 <- full_ps_genus %>% 
  subset_samples(SampleSpecies == "Trinchesia sibogae") %>% 
  plot_tree(ladderize="left", color="SampleSpecies") +
  scale_color_manual(values = pal[3]) +
  theme(legend.position = "none") +
  ggtitle("Trinchesia sibogae")

p4 <- full_ps_genus %>% 
  subset_samples(Structure == "Sediment") %>% 
  plot_tree(ladderize="left", color="Structure") +
  scale_color_manual(values = pal[4]) +
  theme(legend.position = "none") +
  ggtitle("Sediment")

(p1+p2) / (p3+p4)
ggsave("./Output/Figs/Phylogenetic_Dispersion_by_Plant_Structure.png",dpi=300)

rm(p1);rm(p2);rm(p3);rm(p4)

# VennDiagram at ASV level by Location

A <- full_ps %>% subset_samples(Location == "Pulau Hantu (North)"|Location == "Pulau Hantu (West)"|Location == "Pulau Hantu") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Location == "Pulau Hantu (North)"|Location == "Pulau Hantu (West)"|Location == "Pulau Hantu"))>0) %>% tax_table() %>% row.names()
B <- full_ps %>% subset_samples(Location == "Pulau Jong") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Location == "Pulau Jong"))>0) %>% tax_table() %>% row.names()
C <- full_ps %>% subset_samples(Location == "Pulau Tekukor") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Location == "Pulau Tekukor"))>0) %>% tax_table() %>% row.names()
D <- full_ps %>% subset_samples(Location == "Cyrene Reef") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Location == "Cyrene Reef"))>0) %>% tax_table() %>% row.names()


full <- unique(c(A,B,C,D))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n14 <- sum(full %in% unique(A) & full %in% unique(D))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n24 <- sum(full %in% unique(B) & full %in% unique(D))
n34 <- sum(full %in% unique(C) & full %in% unique(D))
n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))
n124 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D))
n134 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D))
n234 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D))
n1234 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

venn.plot2 <- draw.quad.venn(area1=length((A)),
                             area2=length((B)),
                             area3=length((C)),
                             area4=length((D)),
                             n12,n13,n14,n23,n24,n34,n123,n124,n134,n234,n1234,
                             category = c("Pulau Hantu","Pulau Jong","Pulau Tekukor","Cyrene Reef"),
                             fill=pal)
dev.off()
png("./Output/Figs/VennDiagram_Shared_ASV-Level_Taxa_by_Location.png")
grid.draw(venn.plot2)
dev.off()

# VennDiagram at ASV level of overlap for Species ####

A <- full_ps %>% subset_samples(Species == "Pteraeolidia semperi") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Species == "Pteraeolidia semperi"))>0) %>% tax_table() %>% row.names()
B <- full_ps %>% subset_samples(Species == "Trinchesia sibogae") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Species == "Trinchesia sibogae"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B))
n12 <- sum(full %in% unique(A) & full %in% unique(B))

venn.plot1 <- draw.pairwise.venn(area1=length((A)),
                                 area2=length((B)),
                                 n12,
                                 category = c("P. semperi","T. sibogae"),
                                 fill = pal[1:2])
dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_ASV.png")
grid.draw(venn.plot1)
dev.off()

full_ps_spp <- tax_glom(full_ps,c("Genus","Species"))

# VennDiagram at Species level

A <- full_ps_spp %>% subset_samples(Location == "Pulau Hantu (North)"|Location == "Pulau Hantu (West)"|Location == "Pulau Hantu") %>% 
  subset_taxa(taxa_sums(full_ps_spp %>% subset_samples(Location == "Pulau Hantu (North)"|Location == "Pulau Hantu (West)"|Location == "Pulau Hantu"))>0) %>% tax_table() %>% row.names()
B <- full_ps_spp %>% subset_samples(Location == "Pulau Jong") %>% 
  subset_taxa(taxa_sums(full_ps_spp %>% subset_samples(Location == "Pulau Jong"))>0) %>% tax_table() %>% row.names()
C <- full_ps_spp %>% subset_samples(Location == "Pulau Tekukor") %>% 
  subset_taxa(taxa_sums(full_ps_spp %>% subset_samples(Location == "Pulau Tekukor"))>0) %>% tax_table() %>% row.names()
D <- full_ps_spp %>% subset_samples(Location == "Cyrene Reef") %>% 
  subset_taxa(taxa_sums(full_ps_spp %>% subset_samples(Location == "Cyrene Reef"))>0) %>% tax_table() %>% row.names()


full <- unique(c(A,B,C,D))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n14 <- sum(full %in% unique(A) & full %in% unique(D))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n24 <- sum(full %in% unique(B) & full %in% unique(D))
n34 <- sum(full %in% unique(C) & full %in% unique(D))
n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))
n124 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D))
n134 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D))
n234 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D))
n1234 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D))

venn.plot2 <- draw.quad.venn(area1=length((A)),
                             area2=length((B)),
                             area3=length((C)),
                             area4=length((D)),
                             n12,n13,n14,n23,n24,n34,n123,n124,n134,n234,n1234,
                             category = c("Pulau Hantu","Pulau Jong","Pulau Tekukor","Cyrene Reef"),
                             fill=pal)
dev.off()
png("./Output/Figs/VennDiagram_Shared_Spp-Level_Taxa_by_Location.png")
grid.draw(venn.plot2)
dev.off()

# VennDiagram at Species of overlap for Species ####

A <- full_ps_spp %>% subset_samples(Species == "Pteraeolidia semperi") %>% 
  subset_taxa(taxa_sums(full_ps_spp %>% subset_samples(Species == "Pteraeolidia semperi"))>0) %>% tax_table() %>% row.names()
B <- full_ps_spp %>% subset_samples(Species == "Trinchesia sibogae") %>% 
  subset_taxa(taxa_sums(full_ps_spp %>% subset_samples(Species == "Trinchesia sibogae"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B))
n12 <- sum(full %in% unique(A) & full %in% unique(B))

venn.plot1 <- draw.pairwise.venn(area1=length((A)),
                                 area2=length((B)),
                                 n12,
                                 category = c("P. semperi","T. sibogae"),
                                 fill = pal[1:2])
dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_ASV.png")
grid.draw(venn.plot1)
dev.off()


