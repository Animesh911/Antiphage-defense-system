##########################################
# 
# R programming to analyze CRISPRCasTyper (cctyper) output result
# Genomes are already classified using GTDB-Tk
#
# Author: Animesh Kumar
#
###########################################

library(tidyr)
library(dplyr) 
library(reshape2)
library(ggplot2)
library(tidyverse)

#Read metadata genomes
metadata_genomes <- read.csv("..\data\metadata_genomes.txt", sep="\t", header=T) %>% as_tibble() %>% 
  distinct(accession_genbank, .keep_all = TRUE) %>% 
  filter(classification != "#N/A") %>%                                                  # filter #N/A
  filter(!grepl("d__Archaea", classification)) %>%                                      # filter archaea
  mutate(gtdb_copy = classification) %>%
  separate(col = gtdb_copy, into = paste0("Column_", 1:7), sep = ";") %>%
  mutate(across(starts_with("Column_2"), ~str_replace(., "_A$|_B$|_F$", ""), .names = "Phylum")) %>% 
  mutate(across(starts_with("Phylum"), ~str_replace(., "^p__", ""), .names = "Phylum")) %>% 
  group_by(Phylum) %>% 
  mutate(P_count = n()) %>% 
  mutate(Phylum_count = paste(Phylum, P_count, sep = "\n")) %>% 
  mutate(Phylum_count_accession = paste(Phylum_count, accession_genbank, sep = " ")) %>% 
  mutate(across(starts_with("Column_5"), ~str_replace(., "_A$|_B$|_F$", ""), .names = "Family")) %>% 
  mutate(across(starts_with("Family"), ~str_replace(., "^f__", ""), .names = "Family")) %>% 
  group_by(Family) %>% mutate(F_count = n()) %>% 
  mutate(Family_count = paste(Family, F_count, sep = "\n")) %>% 
  mutate(Family_count_accession = paste(Family_count, accession_genbank, sep = " ")) %>% 
  ungroup()

#Read crisprs_near_cas as predicted by CRISPRCasTyper
crispr_data <- read.csv("..\data\merged_crisprs_near_cas.tab", sep="\t", header=T) %>% as_tibble() 

# Fix the header
fix_header <- colnames(crispr_data) #last header falls out in the blank column

#merge metadata_genome with crispr_data: includes rows with no CRISPR leads to          # AccessionID: 938 bacterial
crispr_metadata_genomes <- crispr_data %>%                                 
  subset(!grepl("AccessionID", AccessionID)) %>% 
  separate(col = "AccessionID", into = c("AccessionID", "Contig1"), sep = " ") %>%      # Seperate accession and contigs which are in same column 1
  setNames(fix_header) %>% 
  select(-ncol(.)) %>% # remove last column due to header mismatch  
  separate(col = "AccessionID", into = c("AccessionID", "version"), sep = "\\.") %>%    # Due to the absence of genome version, it will give warning as "Expected 2 pieces. Missing pieces filled with `NA` in 244 rows [107, 108, ...]"
  select(-version)  %>%                                                                 # Warning due to version
  left_join(., metadata_genomes, join_by(AccessionID == accession_genbank)) %>% 
  filter(classification != "#N/A") %>%                                                  # filter #N/A
  filter(!grepl("d__Archaea", classification)) %>%                                      # filter archaea, if any, all converted to #N/A
  filter(Subtype_probability > 0.75)

# Read a tab-separated CRISPR type file
crispr_type <- read.table("..\data\crispr_type.txt", sep = "\t", header = TRUE)
rownames(crispr_type) <- crispr_type[,1]


######################### Plot Number of systems per prokaryotic phylum
# Count total phylum in prokaryotic datasets
cas0_crisprtype_phylum_count_matrix <- crispr_metadata_genomes %>%
  select(Phylum, Phylum_count_accession, P_count, AccessionID, Subtype) %>% 
  left_join(., crispr_type, join_by(Subtype == Prediction)) %>% 
  ungroup()

cas0_crisprtype_phylum_count_matrix %>% 
  select(-Phylum_count_accession) %>% 
  write.csv(., "../final_figures_and_tables/cas0_crisprtype_phylum_count_matrix.csv", row.names = FALSE)

# count genomes falling under a particular phylum and respective CRISPR type
cas1_crisprtype_phylum_heatmap_matrix <- cas0_crisprtype_phylum_count_matrix %>% 
  count(Phylum_count_accession, Subtype.y) %>% 
  pivot_wider(names_from = Subtype.y, values_from = n, values_fill = 0) %>% #A tibble: 166 × 20
  pivot_longer(cols = -Phylum_count_accession, names_to = "Subtype", values_to = "crispr_count") %>% # A tibble: 3,154 × 3
  separate(col = "Phylum_count_accession", into = c("Phylum_count", "accession"), sep = " ") %>% 
  group_by(Phylum_count, Subtype) %>%  
  mutate(Phylum_count1 = Phylum_count) %>%
  separate(col = "Phylum_count1", into = c("Phylum", "P_count"), sep = "\\n") %>% 
  mutate(Number_of_phylum_containing_each_crisprtype = sum(crispr_count > 0), "Total_number_of_different_crisprtype_per_prokaryotic_genomes" = sum(crispr_count)) %>%
  mutate(Frequency_of_each_crisprtype_in_a_given_phylum = Number_of_phylum_containing_each_crisprtype/as.numeric(P_count)) %>% 
  slice(which.max(crispr_count)) %>%                                            # To remove duplicates created during the calculation and don't use crispr_count
  left_join(., crispr_type %>% distinct(Subtype, .keep_all = TRUE), by = "Subtype") %>%
  select(-accession, -crispr_count) 

# PLOT
library(ggh4x)

cas1_crisprtype_phylum_heatmap_plot <- cas1_crisprtype_phylum_heatmap_matrix %>% 
  ggplot(aes(x = Subtype, y = fct_reorder(Phylum_count, as.numeric(P_count)), fill = Frequency_of_each_crisprtype_in_a_given_phylum )) +
  geom_tile() +
  geom_text(aes(label = Number_of_phylum_containing_each_crisprtype), size = 4) +
  labs(x = "CRISPR-Cas subtype", y = "Phylum (genomes count) in dataset")+                   #Phylum (genomes > 0)
  scale_fill_gradientn(colors = c("#efedf5", "#a6bddb", "#2b8cbe"), #"#fee0d2", "#9ebcda", "#8856a7"
                       limits = c(0, .5),
                       breaks = seq(0, 1, by = 0.1),
                       labels = seq(0, 1, by = 0.1)) +
  facet_nested(. ~ crispr_class + crispr_type, scales = "free", space = "free_x", nest_line = element_line(colour = "black")) +
  theme(axis.text = element_text(angle = 0, size = 13), axis.title =element_text(size=14), legend.title=element_blank(),
        strip.text = element_text(size = 13),
        legend.key.height=grid::unit(2, "cm"), legend.key.width=grid::unit(0.3, "cm")) 

cas1_crisprtype_phylum_heatmap_plot

ggsave("../final_figures_and_tables/cas1_crisprtype_phylum_heatmap_plot.png", plot = cas1_crisprtype_phylum_heatmap_plot, width = 13, height = 5.5, dpi = 600)

############################## Plot Distribution of CRISPR Cas Types in phylum dataset
cas2_crisprtype_phylum_distribution <- cas1_crisprtype_phylum_heatmap_matrix %>% ungroup() %>% 
  select(Subtype, Total_number_of_different_crisprtype_per_prokaryotic_genomes, crispr_class) %>%
  filter(Total_number_of_different_crisprtype_per_prokaryotic_genomes > 0) %>% 
  ggplot(aes(x = Subtype, y = Total_number_of_different_crisprtype_per_prokaryotic_genomes)) +
  geom_bar(stat = "identity", fill="lightblue", width = 0.5)+
  facet_grid(~ crispr_class, scales = "free_x", space = "free_x", switch = "x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        strip.placement = "outside", strip.background.x = element_blank(),
        axis.line.y = element_blank(), axis.ticks.y = element_blank()) + labs(y = "", x = "") + scale_y_continuous(expand = c(0, 0))  

ggsave("../final_figures_and_tables/cas2_crisprtype_phylum_distribution_plot.png", plot = cas2_crisprtype_phylum_distribution, width = 4.39, height = 3.12, dpi = 600)
