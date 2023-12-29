##  ###################################################  ##
##  Beta-diversity measures - btwn site and structure    ##
##                                                       ##
##  Author: Ming Sheng - June 16, 2022                   ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  tidyverse v 1.3.1                                    ##
##  phyloseq v 1.40.0                                    ##
##  vegan v 2.6.2                                        ##
##  broom v 0.8.0                                        ##
##  patchwork v 1.1.1                                    ##
##  microbiome v 1.18.0                                  ##
##  purrr v 0.3.4                                        ##
##  corncob v 0.2.0                                      ##
##  indicspecies v 1.7.12                                ##
##                                                       ##
##  ###################################################  ##

# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")
library(corncob); packageVersion("corncob")
library(indicspecies); packageVersion("indicspecies")
library(pairwiseAdonis); packageVersion("pairwiseAdonis")

source("./bbdml_helper.R")

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585", "#FAEBD7")

# Load ps object glom by genus
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")

# Beta-diversity distances and ordinations ####
set.seed(123)
unifrac.dist <- UniFrac(ps_genus,weighted = TRUE,normalized = TRUE,parallel = TRUE)

# Heatmaps of weighted unifrac distances
Psemperi_unifrac.dist <- ps_genus %>% 
  UniFrac(weighted = TRUE,normalized = TRUE,parallel = TRUE)

Psemperi_names_df <- 
  data.frame(sampleid = Psemperi_unifrac.dist %>% 
               as.matrix() %>% 
               colnames(),
             host="P. semperi") %>% 
  mutate(Location = sampleid %>% substr(start=1,stop=2) %>% map_chr(1))

names_df <- Psemperi_names_df %>% 
  mutate(color=case_when(Location == "HG"|Location == "HN"|Location == "HW" ~ pal[1],
                         Location == "TK" ~ pal[2],
                         Location == "JG" ~ pal[3],
                         Location == "SS" ~ pal[4]))
names_df %>% 
  filter(host=="P. semperi") %>% 
  pull(color)

Psemperi_unifrac.dist %>% as.matrix() %>% heatmap(ColSideColors = names_df %>% 
                                                 filter(host=="P. semperi") %>% 
                                                 pull(color),
                                               RowSideColors = names_df %>% 
                                                 filter(host=="P. semperi") %>% 
                                                 pull(color),
                                               keep.dendro = FALSE,
                                               scale = "none",
                                               col=gray.colors(20),main = "P. semperi")

# export
ggsave("./Output/Figs/W-Unifrac_Heatmap_Plot_by_Site.png", dpi = 300)


set.seed(123)

#UniFrac betadisper, PCoA, and PERMANOVA

unifrac.disper <- betadisper(unifrac.dist, group = ps_genus@sam_data$Site)
anova(unifrac.disper)
permutest(unifrac.disper, pairwise = T)
plot(unifrac.disper, main = 'Beta-Dispersion with UniFrac distances')

png("./Output/Figs/W-Unifrac_Beta-Dispersion_by_Site.png")
plot(unifrac.disper, main = 'Beta-Dispersion with UniFrac distances')
dev.off()

ordu <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate("PCoA","unifrac", weighted=TRUE)

uni <- plot_ordination(ps_genus, ordu, color="Site") +
  geom_point(size=5,alpha=.7) + 
  scale_color_manual(values=pal) +
  labs(caption = "MDS/PCoA on weighted-UniFrac distance") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

uni + stat_ellipse(type = 't')

ggsave("./Output/Figs/W-Unifrac_Ordination_Plot_by_Site.pdf")

u.permanova <- adonis2(unifrac.dist ~ ps_genus@sam_data$Site)
u.permanova

pairwise.adonis(unifrac.dist, ps_genus@sam_data$Site)

# beta-dispersion ####
w <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  otu_table() %>% 
  betadiver("w")
w.disper <- betadisper(w,group = meta(ps_genus)$Site)
anova(w.disper)
permutest(w.disper, pairwise = "T")
plot(w.disper,main = "Beta-Dispersion")

png("./Output/Figs/Beta-Dispersion_by_Site.png")
plot(w.disper,main = "Beta-Dispersion")
dev.off()

# PermANOVA ####
ps_ra <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)})

set.seed(123)
permanova <- adonis2(otu_table(ps_ra) ~ ps_ra@sam_data$Site)

sink("./Output/PermANOVA_test_summary.txt")
permanova
sink(NULL)

pairwise.adonis(otu_table(ps_ra), ps_ra@sam_data$Site)

ordu <- ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  ordinate("PCoA","bray")

bray <- plot_ordination(ps_genus, ordu, color="Site") +
  geom_point(size=5,alpha=.7) + 
  scale_color_manual(values=pal) +
  labs(caption = "MDS/PCoA on Bray-Curtis distance") +
  theme_minimal() + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

bray + stat_ellipse(type = 't')

StatChull <- ggproto("StatChull", Stat,
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     },
                     
                     required_aes = c("x", "y")
)


stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
ggsave("./Output/Figs/BrayCurtis_Ordination_Plot_by_Site.pdf")

