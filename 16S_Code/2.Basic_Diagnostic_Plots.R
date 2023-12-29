##  ###################################################  ##
##  Combine phyloseq objects from all runs into one      ##
##  This script processes HN00158973 w/o pitcher plant   ##
##                                                       ##
##  Nudi 16S microbiome                                  ##
##                                                       ##
##  Author: Ming Sheng - June 15, 2022                   ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  tidyverse v 1.3.1                                    ##
##  phyloseq v 1.40.0                                    ##
##  vegan v 2.6.2                                        ##
##                                                       ##
##  ###################################################  ##

# Load packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
source("./plot_bar2.R")

# Summary plots ####

# Plot of taxon-level assignment efficiency 
ps_sp <- readRDS("./Output/noncontam_ps_object.RDS")
full_ps <- ps_sp
phy <- !is.na(tax_table(ps_sp)[,2])
cla <- !is.na(tax_table(ps_sp)[,3])
ord <- !is.na(tax_table(ps_sp)[,4])
fam <- !is.na(tax_table(ps_sp)[,5])
gen <- !is.na(tax_table(ps_sp)[,6])
spp <- !is.na(tax_table(ps_sp)[,7])
assignments <- data.frame(Phylum=phy, Class=cla,Order=ord,Family=fam,Genus=gen,Species=spp)

assignments %>% pivot_longer(1:6) %>% mutate(name=factor(name,levels = c("Phylum","Class","Order","Family","Genus","Species"))) %>%
  ggplot(aes(x=name,fill=value)) + geom_bar() + scale_fill_manual(values=c("Gray","Black")) +
  labs(x="Taxonomic level",y="Count",fill="Unambiguous\nassignment")

ggsave("./Output/Figs/Silva_Taxonomic_Assignment_Efficiency_at_Each_Taxonomic_Rank.png",dpi=300)
rm(phy,cla,ord,fam,gen,spp,assignments,ps_sp)


# Plots of sample sums and taxon sums
png("./Output/Figs/Overall_ESV_Count_Distribution.png")
plot(taxa_sums(full_ps),main = "ESV read count distribution",ylab="Raw ESV Count")
dev.off()

png("./Output/Figs/Overall_Sample_Sum_Distribution.png")
plot(sort(sample_sums(full_ps),decreasing = TRUE),main = "Sample ESV Sums",ylab="Sample Sums")
dev.off()

# Rarefaction curves
png("./Output/Figs/Overall_Rarefaction_Curve.png")
rarecurve(otu_table(full_ps),step = 300,label = FALSE)
dev.off()

# ESV richness distribution
png("./Output/Figs/Overall_ESV_Richness_Distribution.png")
plot(sort(specnumber(otu_table(full_ps)),decreasing = TRUE),
     main="ESV Richness by Sample",ylab="ESV Richness")
dev.off()

# alert that script is finished with beep
beepr::beep(sound=8)
