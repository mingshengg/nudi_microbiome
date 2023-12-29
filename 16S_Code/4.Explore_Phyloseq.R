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
full_ps_genus <- full_ps_genus %>% subset_samples(Site != "Blank")
# output genus-level ps_object for convenience
saveRDS(full_ps_genus, "./Output/full_ps_object_w_tree_genus-glom.RDS")
full_ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")

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
grid.draw(venn.plot1)
dev.off()

# VennDiagram of overlap by Site (species combined)(Pool Pulau Hantu together) ####

A <- full_ps_genus %>% subset_samples(Site == "Pulau Hantu") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Site == "Pulau Hantu"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(Site == "Jong") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Site == "Jong"))>0) %>% tax_table() %>% row.names()
C <- full_ps_genus %>% subset_samples(Site == "Tekukor") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Site == "Tekukor"))>0) %>% tax_table() %>% row.names()
D <- full_ps_genus %>% subset_samples(Site == "Cyrene") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Site == "Cyrene"))>0) %>% tax_table() %>% row.names()
E <- full_ps_genus %>% subset_samples(Site == "Semakau") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Site == "Semakau"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B,C,D,E))

n12 <- sum(full %in% unique(A) & full %in% unique(B))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n14 <- sum(full %in% unique(A) & full %in% unique(D))
n15 <- sum(full %in% unique(A) & full %in% unique(E))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n24 <- sum(full %in% unique(B) & full %in% unique(D))
n25 <- sum(full %in% unique(B) & full %in% unique(E))
n34 <- sum(full %in% unique(C) & full %in% unique(D))
n35 <- sum(full %in% unique(C) & full %in% unique(E))
n45 <- sum(full %in% unique(D) & full %in% unique(E))
n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))
n124 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D))
n125 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(E))
n134 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D))
n135 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(E))
n145 <- sum(full %in% unique(A) & full %in% unique(D) & full %in% unique(E))
n234 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D))
n235 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(E))
n245 <- sum(full %in% unique(B) & full %in% unique(D) & full %in% unique(E))
n345 <- sum(full %in% unique(C) & full %in% unique(D) & full %in% unique(E))
n1234 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D))
n1235 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(E))
n1245 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(D) & full %in% unique(E))
n1345 <- sum(full %in% unique(A) & full %in% unique(C) & full %in% unique(D) & full %in% unique(E))
n2345 <- sum(full %in% unique(B) & full %in% unique(C) & full %in% unique(D) & full %in% unique(E))
n12345 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C) & full %in% unique(D) %in% unique(E))


venn.plot2 <- draw.quintuple.venn(area1=length((A)),
                             area2=length((B)),
                             area3=length((C)),
                             area4=length((D)),
                             area5=length((E)),
                             n12,n13,n14,n15,n23,n24,n25,n34,n35,n45,n123,n124,n125,n134,n135,n145,n234,n235,n245,n345,
                             n1234,n1235,n1245,n1345,n2345,n12345,
                             category = c("Pulau Hantu","Pulau Jong","Pulau Tekukor","Cyrene Reef","Semakau"),
                             fill=c(pal,"#000000"))
dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_Location.png")
grid.draw(venn.plot2)
dev.off()

# VennDiagram of overlap by different parts of Pulau Hantu ####

A <- full_ps_genus %>% subset_samples(Island == "Hantu North") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Island == "Hantu North"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(Island == "HantuWest") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Island == "HantuWest"))>0) %>% tax_table() %>% row.names()


full <- unique(c(A,B))
n12 <- sum(full %in% unique(A) & full %in% unique(B))

venn.plot2 <- draw.pairwise.venn(area1=length((A)),
                             area2=length((B)),
                             n12,
                             category = c("Pulau Hantu (North)","Pulau Hantu (West)"),
                             fill=pal[1:2])
dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_Hantu_Parts.png")
grid.draw(venn.plot2)
dev.off()

# Explore phylogenetic tree ####
full_ps_genus@sam_data$Species %>% unique()
# colored by phylum (blanks removed)
names(full_ps_genus@sam_data) <- c("SampleID","SampleSpecies","Island","PCR_Negative","Site")

p1 <- full_ps_genus %>% 
  subset_samples(SampleSpecies != "NA") %>% 
  plot_tree(ladderize="left", color="SampleSpecies") +
  scale_color_viridis_d()

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
  subset_samples(SampleSpecies == "Pteraeolidia semperi complex") %>% 
  plot_tree(ladderize="left", color="SampleSpecies") +
  scale_color_manual(values = pal[3]) +
  theme(legend.position = "none") +
  ggtitle("Pteraeolidia semperi complex")

(p2)/(p3 + p4)
ggsave("./Output/Figs/Phylogenetic_Dispersion_by_Species.png",dpi=300)


p1 <- full_ps_genus %>% ##
  subset_samples(Site == "Pulau Hantu") %>% 
  plot_tree(ladderize="left", color="Site") +
  scale_color_manual(values = pal[1]) +
  theme(legend.position = "none") +
  ggtitle("Pulau Hantu")


p2 <- full_ps_genus %>% ##
  subset_samples(Site == "Cyrene") %>% 
  plot_tree(ladderize="left", color="Site") +
  scale_color_manual(values = pal[2]) +
  theme(legend.position = "none") +
  ggtitle("Cyrene")

p3 <- full_ps_genus %>% ##
  subset_samples(Site == "Jong") %>% 
  plot_tree(ladderize="left", color="Site") +
  scale_color_manual(values = pal[3]) +
  theme(legend.position = "none") +
  ggtitle("Jong")

p4 <- full_ps_genus %>% ##
  subset_samples(Site == "Semakau") %>% 
  plot_tree(ladderize="left", color="Site") +
  scale_color_manual(values = pal[4]) +
  theme(legend.position = "none") +
  ggtitle("Semakau")

p5 <- full_ps_genus %>% ##
  subset_samples(Site == "Tekukor") %>% 
  plot_tree(ladderize="left", color="Site") +
  scale_color_manual(values = "#ff0000") +
  theme(legend.position = "none") +
  ggtitle("Tekukor")

(p1+p2) / (p3+p4+p5)
ggsave("./Output/Figs/Phylogenetic_Dispersion_by_Site.png",dpi=300)

rm(p1);rm(p2);rm(p3);rm(p4);rm(p5)

# VennDiagram at ASV level by Site

A <- full_ps %>% subset_samples(Site == "Pulau Hantu") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Site == "Pulau Hantu"))>0) %>% tax_table() %>% row.names()
B <- full_ps %>% subset_samples(Site == "Jong") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Site == "Jong"))>0) %>% tax_table() %>% row.names()
C <- full_ps %>% subset_samples(Site == "Tekukor") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Site == "Tekukor"))>0) %>% tax_table() %>% row.names()
D <- full_ps %>% subset_samples(Site == "Cyrene") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Site == "Cyrene"))>0) %>% tax_table() %>% row.names()
E <- full_ps %>% subset_samples(Site == "Semakau") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Site == "Semakau"))>0) %>% tax_table() %>% row.names()

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


