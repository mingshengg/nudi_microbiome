##  ###################################################  ##
##  Alpha-diversity measures - btwn site and structure   ##
##                                                       ##
##  Author: Ming Sheng - June 17, 2022                   ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  tidyverse v 1.3.1                                    ##
##  phyloseq v 1.40.0                                    ##
##  vegan v 2.6.2                                        ##
##  broom v 0.8.0                                        ##
##  patchwork v 1.1.1                                    ##
##  microbiome v 1.18.0                                  ##
##  emmeans v 1.7.4.1                                    ##
##                                                       ##
##  ###################################################  ##

# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(emmeans); packageVersion("emmeans")
source("./plot_bar2.R")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585", "#FAEBD7")

# Load ps objects
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")


# Model alpha diversity ####
meta <- microbiome::meta(ps_genus)                 
meta$Shannon <- vegan::diversity(otu_table(ps_genus),index = "shannon")
meta$Richness <- vegan::specnumber(otu_table(ps_genus))

# add to ps object
ps_genus@sam_data$Richness <- meta$Richness
ps_genus@sam_data$Shannon <- meta$Shannon

# lme models
shannon_mod <- aov(data = meta,
                              formula = Shannon ~ Site)
richness_mod <- aov(data = meta,
                               formula = Richness ~ Site)
# send to file
sink("./Output/Shannon_Diversity_Model.txt")
summary(shannon_mod)
anova(shannon_mod)
sink(NULL)

emmeans(shannon_mod, pairwise ~ Site)
# violin plots for shannon and richness 

shannon_plot <- meta %>% 
  ggplot(aes(x = Site, y = Shannon, fill = Site)) + 
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.1, fill = 'white') +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  scale_fill_manual(values = pal) +
  labs(y = 'Shannon diversity index')


shannon_plot
ggsave("./Output/Figs/shannon_by_site.pdf",width = 12,height = 6)

richness_plot <- meta %>% 
  ggplot(aes(x = Site, y = Richness, fill = Site)) + 
  geom_violin(scale = 'width') +
  geom_boxplot(width = 0.1, fill = 'white') +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  scale_fill_manual(values = pal)

ggsave("./Output/Figs/richness_by_site.pdf",width = 12,height = 6)
# Stacked phylum-level barcharts (horizontal) by site ####

# merge samples
ps_genus@sam_data %>% names()
ps_genus@otu_table %>% rowSums()
newmergevar <- paste(ps_genus@sam_data$Species,
                     ps_genus@sam_data$Location,
                     sep = "_")
ps_genus@sam_data$newmergevar <- newmergevar

psm_ra <- ps_genus %>%   
  merge_samples(newmergevar,fun = "sum") %>% 
  transform_sample_counts(function(x){x/sum(x)})

# repair metadata
sp <- psm_ra@sam_data %>% row.names() %>% str_split("_") %>% map_chr(1)
lo <- psm_ra@sam_data %>% row.names() %>% str_split("_") %>% map_chr(2)
psm_ra@sam_data$Species <- sp
psm_ra@sam_data$Location <- lo

top_20 <- psm_ra@otu_table %>% colSums() %>% 
  sort(decreasing = T) %>% 
  as.data.frame() %>% 
  row.names() %>%
  .[1:20]

taxtab20 <- cbind(tax_table(psm_ra), genus_19 = NA)
taxtab20[top_20,"genus_19"] <- as(tax_table(psm_ra)[top_20,"Genus"],"character")
tax_table(psm_ra) <- tax_table(taxtab20)

palette <- c("#d43e6e",
             "#7f352d",
             "#d45030",
             "#d28f73",
             "#d2993b",
             "#7d602a",
             "#a8b639",
             "#7a964e",
             "#57bb47",
             "#375a31",
             "#5fb58f",
             "#59a5c3",
             "#455485",
             "#6d75d9",
             "#593093",
             "#8142d5",
             "#ae8fca",
             "#ce4ebf",
             "#713160",
             "#d17ba0")

ps_plot <- psm_ra %>% 
  plot_bar2(x="Location",fill="genus_19") +
  theme_minimal() +
  scale_fill_manual(values = palette)

ps_plot

ggsave("./Output/Figs/horizontal_bar_charts_ITS.pdf",dpi=400,width = 16,height = 8)

ps_plot <- psm_ra %>% 
  subset_samples(Species == "Pteraeolidia semperi") %>% 
  plot_bar2(x="Site",fill="Phylum") +
  coord_flip() +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=8,margin = margin(c(5,10,5,10))),
        axis.text.x = element_text(face="bold",size=8,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face="bold",size=10),
        axis.title = element_text(face="bold",size=12),
        legend.title = element_text(face="bold",size=12),
        legend.text = element_text(size=10),
        plot.title = element_text(face="bold.italic",size=14,hjust=.5),
        legend.position = "none",
        plot.margin = margin(r=20)) +
  scale_fill_viridis_d() +
  scale_y_continuous(breaks=c(0,.5,1)) +
  labs(y="\nRelative abundance",title = "Ptaraeolidia semperi")

sa_plot <- psm_ra %>% 
  subset_samples(Species == "Sonneratia alba") %>% 
  plot_bar2(x="Location",fill="Phylum",facet_grid = ~Structure) +
  coord_flip() +
  theme_minimal() +
  theme(strip.text = element_text(face="bold",size=8,margin = margin(c(5,-20,5,-20))),
        axis.text.x = element_text(face="bold",size=8,angle=90,hjust=1,vjust=.5),
        axis.text.y = element_text(face="bold",size=10),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_blank(),
        legend.title = element_text(face="bold",size=12),
        legend.text = element_text(size=10),
        plot.title = element_text(face="bold.italic",size=14,hjust=.5)) +
  scale_fill_viridis_d() +
  scale_y_continuous(breaks=c(0,.5,1)) +
  labs(y="\nRelative abundance",title = "Sonneratia alba")

aa_plot + sa_plot
ggsave("./Output/Figs/horizontal_bar_charts.png",dpi=400,width = 16,height = 8)
