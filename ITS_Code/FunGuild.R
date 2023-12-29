# INstallation
## devtools::install_github("brendanf/FUNGuildR")

library(FUNGuildR)
library(dplyr)
library(tidyr)
library(stringr)

ps_genus <- readRDS("./Output/full_ps_object_w_tree_genus-glom.RDS")

# rename ps_genusal taxa
tax_table(ps_genus)[,1] <- tax_table(ps_genus)[,1] %>% str_remove_all("k__")
tax_table(ps_genus)[,2] <- tax_table(ps_genus)[,2] %>% str_remove_all("p__")
tax_table(ps_genus)[,3] <- tax_table(ps_genus)[,3] %>% str_remove_all("c__")
tax_table(ps_genus)[,4] <- tax_table(ps_genus)[,4] %>% str_remove_all("o__")
tax_table(ps_genus)[,5] <- tax_table(ps_genus)[,5] %>% str_remove_all("f__")
tax_table(ps_genus)[,6] <- tax_table(ps_genus)[,6] %>% str_remove_all("g__")
tax_table(ps_genus)[,7] <- tax_table(ps_genus)[,7] %>% str_remove_all("s__")

# create new dataframe
fungi_data <- ps_genus@tax_table %>% 
  as.data.frame() %>%
  unite('Taxonomy',Kingdom:Genus, sep =';') %>%
  select('Taxonomy')
  
# Main function of FUNGuild
fun_guilds <- funguild_assign(fungi_data)
rownames(fun_guilds) <- rownames(fungi_data)

write.csv(fun_guilds, "./Output/FUNGuild_Output.csv", row.names = T)

# Inspect common taxa
fun_guilds['ASV60',]
  
fun_guilds %>% 
  group_by(trophicMode) %>%
  summarise(no_rows = length(trophicMode))
