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

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585", "#1ABC9C")

# Load ps objects
ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")

sample_data(ps_genus)$Site[1:26] <- "Hantu"
pseq.rel <- microbiome::transform(ps_genus, "compositional")


# find detection thresholds
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

ps_genus %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_core(prevalences = prevalences, 
            detections = det, 
            plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

det <- seq(from = 50, to = round(max(abundances(ps_genus))/10, -1), by = 300)
prevalences <- seq(.05, 1, .05)
ps_genus %>%
  plot_core(plot.type = "heatmap",
            prevalences = prevalences,
            detections = det,
            colours = rev(brewer.pal(5, "Spectral")),
            min.prevalence = .2, horizontal = T)

# find core microbiome members overall
overall_core_taxa <- core_members(pseq.rel,
                                  detection = 0.01, ##change
                                  prevalence = 0.3,
                                  include.lowest = F)

length(overall_core_taxa)

# subset to just those taxa
ps_core <- pseq.rel %>% 
  subset_taxa(taxa_names(pseq.rel) %in% overall_core_taxa)

overall_core_taxa_names <- corncob::otu_to_taxonomy(data=ps_core,level = c("Family","Genus"),taxa_names(ps_core))
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

# Sorting them as factors based on mean abundance so that the most abundant one will appear on top
gen <- overall_psm
gen$Genus <- as.factor(gen$Genus)
gen$Genus <- reorder(gen$Genus, gen$Abundance, FUN = mean)

detection_50 <- gen %>% ##change
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

detection_50 ##change
ggsave('./Output/Figs/bact_detection_50.png', dpi = 400, width = 10, height = 6) ##change

se <- function(x) sqrt(var(x) / length(x))
pseq.rel@otu_table[,'ASV9']   %>% sd()

pseq.rel@otu_table %>% 
  as.data.frame %>% 
  select(ASV1,ASV9,ASV14) %>%
  mutate(comp = ASV1 + ASV9 + ASV14) %>%
  mutate(Site = pseq.rel@sam_data$Site) %>%
  ggplot(aes(x = Site, y = comp, fill = Site)) +
  geom_boxplot() +
  theme_minimal()

pseq.rel@tax_table[,5] %>%
  as.data.frame %>%
  filter(Family == 'Bacillaceae')

pseq.rel@otu_table[,'ASV636']

#Summary stats for core taxa
summary_core_taxa <- data.frame(Name = character(),
                                Mean = double(),
                                SD = double(),
                                SE = double(),
                                Max = double(),
                                Min = double(),
                                Det = double())
for (i in overall_core_taxa){
  summary_core_taxa[i,]<- c(pseq.rel@tax_table[i,6],
                            pseq.rel@otu_table[,i] %>% mean(),
                            pseq.rel@otu_table[,i] %>% sd(),
                            pseq.rel@otu_table[,i] %>% se(),
                            pseq.rel@otu_table[,i] %>% max(),
                            pseq.rel@otu_table[,i] %>% min(),
                            pseq.rel@otu_table[,i] %>% as.logical() %>% count() %>% filter(x== 'TRUE') %>% select(freq)/46)
}

summary_core_taxa
# find core for each plant site
sites <- pseq.rel@sam_data$Site %>% unique()

for(i in sites){
  pseq.rel %>% 
    subset_samples(Site == i) %>% 
    core_members(detection = 0.5,
                 prevalence = 0,
                 include.lowest = F) %>% 
    assign(paste0(i,"_core_members"),value=.,envir = .GlobalEnv)
}

# combine them
all_core_members <- c(Cyrene_core_members,Jong_core_members,
                      Hantu_core_members,Tekukor_core_members,Semakau_core_members) %>% unique()

ps_genus_all_core <- pseq.rel %>% 
  subset_taxa(taxa_names(pseq.rel) %in% all_core_members)

# rename taxa to family/genus
core_taxa_names <- corncob::otu_to_taxonomy(data=ps_genus_all_core,level = c("Family","Genus"),taxa_names(ps_genus_all_core))
taxa_names(ps_genus_all_core) <- core_taxa_names


psm <- ps_genus_all_core %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  psmelt()

site <- ps_genus_all_core %>% 
  merge_samples("Site",fun="sum")

site@sam_data$Site <- row.names(site@sam_data)

psm <- psm %>%
  arrange(Site,sample_Species) %>% 
  mutate(sample_Species = factor(sample_Species,levels = unique(sample_Species)),
         Site = factor(Site,levels = unique(Site)),
         Host_Site = (paste0(sample_Species,"_",Site)))

psm$Abundance[is.na(psm$Abundance)] <- 0

psm$SampleID <- factor(psm$SampleID,levels=unique(psm$SampleID))

psm <- psm %>% 
  mutate(Str_color = case_when(Site == "Hantu" ~ pal[1],
                               Site == "Jong" ~ pal[2],
                               Site == "Cyrene" ~ pal[3],
                               Site == "Semakau" ~ pal[4],
                               Site == "Tekukor" ~ pal[5]))
psm %>%   
  ggplot(aes(x=Sample,y=Genus,fill=Abundance)) +
  geom_tile() +
  facet_grid(cols = vars(Site),
             scales = 'free') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        strip.text = element_text(face="bold",size=14),
        axis.text.y = element_text(size=10)) +
  scale_fill_viridis_c() +
  labs(y="Taxon",fill="Relative\nabundance")
ggsave("./Output/Figs/Core_Heatmap_by_Site.png",dpi=400,height = 8,width = 12)

count(head(psm[,1:5][order(psm$Abundance, decreasing = TRUE),], n = 100)$OTU == 'Mycoplasmataceae_Mycoplasma')

#lowest reads SS05 vs highest reads HN08
SS05 <- ps_genus %>% subset_samples(SampleID == 'SS05')
HN08 <- ps_genus %>% subset_samples(SampleID == 'HN08')

SS05@otu_table %>% table()
HN08@otu_table %>% table()

#Prevalence of taxonomic groups
prevalence(ps_genus, detection = 10, sort = T, count = T)

#Try with compositional data?
pseq.rel <- microbiome::transform(ps_genus, "compositional")
prevalence(pseq.rel, detection = 0.1, sort = T, count = T)
overall_core_taxa <- core_members(pseq.rel,
                                  detection = 0.1,
                                  prevalence = 0.1,
                                  include.lowest = FALSE)

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

pseq.rel %>% 
  plot_core(prevalences = prevalences, 
            detections = det, 
            plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")


