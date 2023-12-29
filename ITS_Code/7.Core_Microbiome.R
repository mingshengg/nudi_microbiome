##  ###################################################  ##
##  Core Microbiome Members                              ##
##                                                       ##
##  Author: Ming Sheng - August, 2022                    ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.0.3                                            ##
##  tidyverse v 1.3.0                                    ##
##  phyloseq v 1.32.0                                    ##
##  vegan v 2.5.6                                        ##
##  broom v 0.7.1                                        ##
##  patchwork v 1.0.1                                    ##
##  microbiome v 1.10.0                                  ##
##                                                       ##
##  ###################################################  ##

# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(RColorBrewer)
source("./R/plot_bar2.R")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585", "#1ABC9C")

# Load ps objects
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")

# rename ps_genusal taxa
tax_table(ps_genus)[,1] <- tax_table(ps_genus)[,1] %>% str_remove_all("k__")
tax_table(ps_genus)[,2] <- tax_table(ps_genus)[,2] %>% str_remove_all("p__")
tax_table(ps_genus)[,3] <- tax_table(ps_genus)[,3] %>% str_remove_all("c__")
tax_table(ps_genus)[,4] <- tax_table(ps_genus)[,4] %>% str_remove_all("o__")
tax_table(ps_genus)[,5] <- tax_table(ps_genus)[,5] %>% str_remove_all("f__")
tax_table(ps_genus)[,6] <- tax_table(ps_genus)[,6] %>% str_remove_all("g__")
tax_table(ps_genus)[,7] <- tax_table(ps_genus)[,7] %>% str_remove_all("s__")


pseq.rel_fung <- microbiome::transform(ps_genus, "compositional")


# find detection thresholds
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

pseq.rel_fung %>%
  plot_core(prevalences = prevalences, 
            detections = det, 
            plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

det <- seq(from = 50, to = round(max(abundances(ps_genus))/10, -1), by = 100)
prevalences <- seq(.05, 1, .05)
ps_genus %>%
  plot_core(plot.type = "heatmap",
            prevalences = prevalences,
            detections = det,
            colours = rev(brewer.pal(5, "Spectral")),
            min.prevalence = .01, horizontal = T)

# find core microbiome members overall
overall_core_taxa <- core_members(pseq.rel_fung,
                                  detection = 0,
                                  prevalence = 0.5,
                                  include.lowest = F)
#prevalence = presence/absence
#detection = relative abundance

length(overall_core_taxa)


# subset to just those taxa
ps_core <- pseq.rel_fung %>% 
  subset_taxa(taxa_names(pseq.rel_fung) %in% overall_core_taxa)

overall_core_taxa_names <- corncob::otu_to_taxonomy(data=ps_core,level = c("Order","Family","Genus"),taxa_names(ps_core))
taxa_names(ps_core) <- overall_core_taxa_names

overall_psm <- ps_core %>% 
  psmelt()

site <- ps_core %>% 
  merge_samples("Site",fun="sum")

site@sam_data$Site <- row.names(site@sam_data)

overall_psm <- overall_psm %>%
  arrange(Site,sample_Species) %>% 
  mutate(sample_Species = factor(sample_Species,levels = unique(sample_Species)),
         Site = factor(Site,levels = unique(Site)),
         Host_Site = (paste0(sample_Species,"_",Site)))

psm$Abundance[is.na(psm$Abundance)] <- 0

psm$SampleID <- factor(psm$SampleID,levels=unique(psm$SampleID))

overall_psm <- overall_psm %>% 
  mutate(Str_color = case_when(Site == "Hantu" ~ pal[1],
                               Site == "Jong" ~ pal[2],
                               Site == "Cyrene" ~ pal[3],
                               Site == "Semakau" ~ pal[4],
                               Site == "Tekukor" ~ pal[5]))

gen <- overall_psm
gen$Genus <- as.factor(gen$Genus)
gen$Genus <- reorder(gen$Genus, gen$Abundance, FUN = mean)

detection_50_fung <- gen %>%   
  ggplot(aes(x=Sample,y=Genus,fill=Abundance)) +
  geom_tile() +
  facet_grid(cols = vars(Site),
             scales = 'free') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face="bold",size=14),
        axis.text.y = element_text(size=12)) +
  scale_fill_viridis_c() +
  labs(y="Taxon",fill="Relative\nabundance")

detection_50_fung
ggsave('./Output/Figs/detection_50_fung.png', dpi = 400, width = 10, height = 6) ##change


se <- function(x) sqrt(var(x) / length(x))
summary_core_taxa <- data.frame(Name = character(),
                                Mean = double(),
                                SD = double(),
                                SE = double(),
                                Max = double(),
                                Min = double(),
                                Det = double())
for (i in overall_core_taxa){
  summary_core_taxa[i,]<- c(pseq.rel_fung@tax_table[i,2],
                            pseq.rel_fung@otu_table[,i] %>% mean(),
                            pseq.rel_fung@otu_table[,i] %>% sd(),
                            pseq.rel_fung@otu_table[,i] %>% se(),
                            pseq.rel_fung@otu_table[,i] %>% max(),
                            pseq.rel_fung@otu_table[,i] %>% min(),
                            pseq.rel_fung@otu_table[,i] %>% as.logical() %>% sum()/43)
}

summary_core_taxa
