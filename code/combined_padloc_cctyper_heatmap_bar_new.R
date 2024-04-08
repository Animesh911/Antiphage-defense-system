##########################################
# 
# R programming to analyze both PADLOC and CRISPRCasTyper (cctyper) output result
# Generates barplot and heatmap
# 
# Author: Animesh Kumar
#
###########################################

library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(tidyverse)

# Read a tab-separated metadata file of defence system
bact_defence_combined <- read.table("..\data\metadata_bact_defence_combined.txt", sep = "\t", header = TRUE)
rownames(bact_defence_combined) <- bact_defence_combined[,1]                                                      # both are CRISPRcas subsystems I-F and I-F_T

# Read metadata of genomes
metadata_genomes <- read.csv("..\data\metadata_genomes.txt", sep="\t", header=T) %>% as_tibble() %>% 
  distinct(accession_genbank, .keep_all = TRUE) %>% 
  filter(classification != "#N/A") %>%                                                                            # filter #N/A
  filter(!grepl("d__Archaea", classification)) %>%                                                                # filter archaea
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


# Read the padloc file
padloc_data <- read.csv("..\data\all_padloc.csv") %>% as_tibble() %>% 
  filter(full.seq.E.value < 0.01) %>% 
  filter(domain.iE.value < 0.01) %>% 
  filter(target.coverage > 0.8) %>%
  filter(hmm.coverage > 0.8) %>%
  distinct(AccessionID, system, system.number) %>% 
  select(AccessionID, system)


# Fix the header
fix_header <- colnames(read.csv("..\data\merged_crisprs_near_cas.tab", sep="\t", header=T) %>% as_tibble())         # last header falls out in the blank column

# Read the crispr file 
crispr_data <- read.csv("..\data\merged_crisprs_near_cas.tab", sep="\t", header=T) %>% as_tibble() %>%
  subset(!grepl("AccessionID", AccessionID)) %>%
  separate(col = "AccessionID", into = c("AccessionID", "Contig1"), sep = " ") %>%                                  # Seperate accession and contigs which are in same column 1
  setNames(fix_header) %>% 
  select(-ncol(.)) %>%                                                                                              # remove last column due to header mismatch  
  separate(col = "AccessionID", into = c("AccessionID", "version"), sep = "\\.") %>%                                # Due to the absence of genome version, it will give warning as "Expected 2 pieces. Missing pieces filled with `NA` in 244 rows [107, 108, ...]"
  filter(Subtype_probability > 0.75) %>%                                                                            # number could be higher due to presence of !bacteria or 231
  rename(system = Subtype) %>% 
  select(AccessionID, system)                                                                                       # warning due to version

# Append padloc_data0 and crispr_data0 and metadata of genomes and metadata of defence system
data0 <- rbind(padloc_data, crispr_data) %>%
  group_by(AccessionID, system) %>%
  summarise(Frequency = n()) %>%
  ungroup() %>% 
  full_join(., metadata_genomes, join_by(AccessionID == accession_genbank)) %>% 
  filter(classification != "#N/A") %>%                                                                            # filter, already removed archaea turns into #N/A
  filter(!grepl("d__Archaea", classification)) %>%                                                                # filter, already removed archaea
  left_join(., bact_defence_combined, join_by(system == system_type))

data0 %>% 
  select(-Phylum_count, -Phylum_count_accession, -Family_count, -Family_count_accession) %>% 
  write.csv(., "../final_figures_and_tables/system0_count_matrix.csv", row.names = FALSE)

#######################################
######################
##                  ##
##  HEATMAP PHYLUM  ##
##                  ##
######################
#######################################
system1_heatmap_phylum <- data0 %>%
  filter(Frequency>0) %>%                                                                                       # filter(system != "N/A") or remove genome having no system
  count(Phylum_count_accession, bact_defence_system) %>% 
  pivot_wider(names_from = bact_defence_system, values_from = n, values_fill = 0) %>% # A tibble: 883 × 53
  pivot_longer(cols = -Phylum_count_accession, names_to = "bact_defence_system", values_to = "system_count") %>% # A tibble: 45,916 × 3
  separate(col = "Phylum_count_accession", into = c("Phylum_count", "accession"), sep = " ") %>% 
  group_by(Phylum_count, bact_defence_system) %>% 
  mutate(Phylum_count1 = Phylum_count) %>%
  separate(col = "Phylum_count1", into = c("Phylum", "P_count"), sep = "\\n") %>% 
  mutate(Number_of_phylum_containing_each_system = sum(system_count > 0), "Total_number_of_different_system_per_prokaryotic_genomes" = sum(system_count)) %>%
  mutate(Frequency_of_each_system_in_a_given_phylum = Number_of_phylum_containing_each_system/as.numeric(P_count)) %>% 
  slice(which.max(system_count)) %>%                                                                            # To remove duplicates created during the calculation and don't use system_count
  select(-accession, -Phylum) %>% 
  filter(bact_defence_system !=  "CRISPR-Cas_padloc") %>% 
  filter(bact_defence_system !=  "DMS") %>% 
  mutate(bact_defence_system = recode(bact_defence_system, "CRISPR-Cas_cctyper" = "CRISPR-Cas"))

system1_heatmap_phylum %>% 
  separate(col = "Phylum_count", into = c("Phylum", "P_count"), sep = "\\n") %>%
  write.csv(., "../final_figures_and_tables/system1_heatmap_phylum.csv", row.names = FALSE)


#library(tidyverse)
label_angle <- function(label_value) {
  if (label_value > 100) {
    return(90)  # Set angle to 90 degrees if label > 100
  } else {
    return(0)   # Set angle to 0 degrees if label <= 100
  }
}


library(forcats)
system1_heatmap_phylum_plot <- system1_heatmap_phylum %>% 
  filter(as.numeric(P_count) > 5) %>%
  ggplot(aes(x = bact_defence_system, y = fct_reorder(Phylum_count, as.numeric(P_count)), fill = Frequency_of_each_system_in_a_given_phylum )) +
  geom_tile() +
  geom_text(aes(label = Number_of_phylum_containing_each_system, angle = ifelse(Number_of_phylum_containing_each_system > 100, 90, 0)), color = "black", size = 3) +  # Add values to heatmap
  labs(x = "system", y = "phylum (genomes count > 5) in dataset") +
  scale_fill_gradientn(colors = c("#efedf5", "#a6bddb", "#2b8cbe"),
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2),
                       labels = seq(0, 1, by = 0.2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title=element_blank(),
        legend.key.height=grid::unit(2, "cm"), legend.key.width=grid::unit(0.3, "cm"))

system1_heatmap_phylum_plot


ggsave("../final_figures_and_tables/system1_heatmap_phylum_plot.png", plot = system1_heatmap_phylum_plot, width = 9, height = 4.6, dpi = 600)


#######################################
#######################
##                   ##
##  HEATMAP FAMILY   ##
##                   ##
#######################
#######################################
system2_heatmap_family <- data0 %>%
  filter(Frequency>0) %>%                                                                                     # filter(system != "N/A") or remove genome having no system
  count(Family_count_accession, bact_defence_system) %>% 
  pivot_wider(names_from = bact_defence_system, values_from = n, values_fill = 0) %>% # A tibble: 883 × 53
  pivot_longer(cols = -Family_count_accession, names_to = "bact_defence_system", values_to = "system_count") %>% # A tibble: 45,916 × 3
  separate(col = "Family_count_accession", into = c("Family_count", "accession"), sep = " ") %>% 
  group_by(Family_count, bact_defence_system) %>% 
  mutate(Family_count1 = Family_count) %>%
  separate(col = "Family_count1", into = c("Family", "F_count"), sep = "\\n") %>% 
  mutate(Number_of_family_containing_each_system = sum(system_count > 0), "Total_number_of_different_system_per_prokaryotic_genomes" = sum(system_count)) %>%
  mutate(Frequency_of_each_system_in_a_given_family = Number_of_family_containing_each_system/as.numeric(F_count)) %>% 
  slice(which.max(system_count)) %>%                                                                          # To remove duplicates created during the calculation and don't use system_count
  select(-accession, -Family) %>% 
  filter(bact_defence_system !=  "CRISPR-Cas_padloc") %>% 
  filter(bact_defence_system !=  "DMS") %>% 
  mutate(bact_defence_system = recode(bact_defence_system, "CRISPR-Cas_cctyper" = "CRISPR-Cas"))

system2_heatmap_family %>% 
  separate(col = "Family_count", into = c("Family", "F_count"), sep = "\\n") %>% 
  write.csv(., "../final_figures_and_tables/system2_heatmap_family.csv", row.names = FALSE)


#library(tidyverse)
label_angle <- function(label_value) {
  if (label_value > 100) {
    return(90)  # Set angle to 90 degrees if label > 100
  } else {
    return(0)   # Set angle to 0 degrees if label <= 100
  }
}


library(forcats)
system2_heatmap_family_plot <- system2_heatmap_family %>% 
  filter(as.numeric(F_count) > 5) %>%
  ggplot(aes(x = bact_defence_system, y = fct_reorder(Family_count, as.numeric(F_count)), fill = Frequency_of_each_system_in_a_given_family )) +
  geom_tile() +
  geom_text(aes(label = Number_of_family_containing_each_system, angle = ifelse(Number_of_family_containing_each_system > 100, 90, 0)), color = "black", size = 3) +  # Add values to heatmap
  labs(x = "system", y = "Family (genomes count >5 ) in dataset") +
  scale_fill_gradientn(colors = c("#efedf5", "#a6bddb", "#2b8cbe"),
                       limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2),
                       labels = seq(0, 1, by = 0.2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title=element_blank(),
        legend.key.height=grid::unit(2, "cm"), legend.key.width=grid::unit(0.3, "cm"))

system2_heatmap_family_plot

ggsave("../final_figures_and_tables/system2_heatmap_family_plot.png", plot = system2_heatmap_family_plot, width = 10, height = 9, dpi = 600)


#######################################
######################
##                  ##
##  HEATMAP GENUS   ##    few genus missing due to classification at higher taxonomic ranks
##                  ##
######################
#######################################
#######################################
#######################################
# Percentages of genomes encoding in the system, here genomes are counted once if any of the subsystem of a system is present.
library(ggh4x)

system3_percentage_of_genomes <- data0 %>% select(AccessionID, bact_defence_system, Frequency) %>% 
  mutate(Frequency = as.integer(Frequency > 0)) %>%                                                           #convert logical vector to an integer, TRUE becomes 1, and FALSE becomes 0
  distinct() %>% 
  pivot_wider(names_from = AccessionID, values_from = Frequency, values_fill = 0) %>% 
  filter(bact_defence_system != "#N/A") %>% 
  column_to_rownames("bact_defence_system") %>% 
  #mutate(across(everything(), ~ifelse(. > 0, 1, 0))) %>%
  mutate(percentage_system_genome = rowSums(. > 0) /ncol(.) * 100) %>%
  arrange(desc(percentage_system_genome)) %>% 
  rownames_to_column("bact_defence_system") %>% 
  select(bact_defence_system, percentage_system_genome) %>% 
  filter(!(bact_defence_system %in% c("DMS", "CRISPR-Cas_padloc"))) %>% 
  mutate(bact_defence_system = recode(bact_defence_system, "CRISPR-Cas_cctyper" = "CRISPR-Cas"))

system3_percentage_of_genomes %>% 
  write.csv(., "../final_figures_and_tables/system3_percentage_of_genomes.csv", row.names = FALSE)


system3_percentage_of_genomes_barplot <- system3_percentage_of_genomes %>% 
  mutate(System = reorder(row.names(.), desc(percentage_system_genome))) %>%
  ggplot(aes(x = fct_reorder(bact_defence_system, as.numeric(System)), y = percentage_system_genome)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Antiphage systems", y = "Prevalence (%)") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,NA)) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = -90, margin = margin(t = 0.1, unit = "cm")), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray90"))


system3_percentage_of_genomes_barplot

ggsave("../final_figures_and_tables/system3_percentage_of_genomes_barplot.png", bg="white", plot = system3_percentage_of_genomes_barplot, width = 6.25, height = 4.6, dpi = 600)

###########################################
# Abundance of CRISPR-Cas and other systems in cold-living bacteria: VENN diagram
venn_data <- data0 %>% select(AccessionID, bact_defence_system, Frequency) %>% 
  mutate(bact_defence_system = recode(bact_defence_system, "CRISPR-Cas_cctyper" = "CRISPR-Cas")) %>% 
  mutate(Frequency = as.integer(Frequency > 0)) %>%                                                               #convert logical vector to an integer, TRUE becomes 1, and FALSE becomes 0
  distinct() %>% 
  pivot_wider(names_from = bact_defence_system, values_from = Frequency, values_fill = 0) %>% 
  column_to_rownames("AccessionID") %>% 
  select(-`NA`, -`CRISPR-Cas_padloc`, -DMS) %>% 
  mutate(other_system = rowSums(select(., Phosphorothioation, Argonaute, Septu, Thoeris, Zorya, Dpd, PARIS, Mza, Hachiman,
                                       Viperin, BREX, ietAS, Wadjet, GAO, Hhe, Retron, DarTG, `Helicase-DUF2290`, Shedu, 
                                       AVAST, Pycsar, Lamassu, Upx, RADAR, Druantia, Kiwa, qatABCD, Gasdermin, Ppl, Tmn,
                                       Dsr, `TIR-NLR`, Hma, Paris_fused, `TerY-P`, Nhi, DISARM,
                                       DUF4238, `DprA-PRTase`, `3HP`, `Hydrolase-TM`, Abi, DRT, CBASS, Gabija))) %>% 
  select(RM, dXTPase, `CRISPR-Cas`, other_system) %>% 
  rownames_to_column(var = "AccessionID") %>% 
  mutate_at(vars(RM:other_system), ~ifelse(. > 0, AccessionID, .)) %>% 
  select(-AccessionID)

venn_data_list <- lapply(venn_data, function(x) x[x != 0])

#
#library("ggVennDiagram")
#ggVennDiagram(venn_data_list, label_alpha = 0)

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

#display_venn(venn_data_list, 
#             category.names = c("RM" ,  "dXTPase", "CRISPR-Cas", "other_system"),
#             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))

#
library(ggvenn)
system4_venn_plot <- ggvenn(
  venn_data_list, 
  fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"), ## "#0073C2FF", "#EFC000FF", "#868686FF", "#009E73", #"#CD534CFF"
  stroke_size = 0.5, set_name_size = 4,
)

system4_venn_plot

ggsave("../final_figures_and_tables/system4_venn_plot.png", bg="white", plot = system4_venn_plot, width = 7, height = 5.5, dpi = 600)

###########################################
# Abundance of CRISPR-Cas and other systems in cold-living bacteria: PIE diagram
# Create a data frame with two rows: one for CRISPR-Cas_cctyper and one for the rest

pie_data <- data.frame(
  bact_defence_system = c("CRISPR-Cas_cctyper", "Other"),
  percentage_system_genome = c(17.6972281, 100 - 17.6972281)
)


cas3_pie_plot <- ggplot(pie_data, aes(x = "", y = percentage_system_genome, fill = bact_defence_system)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 45) +
  theme_void() +
  scale_fill_manual(values = c("#1499aa", "#d9d9d9"), 
                    labels = c("Genomes with CRISPR-Cas", "Genomes without CRISPR-Cas")) +
  labs(title = "Abundance of CRISPR-Cas systems\nin cold-living bacteria", title.position = "plot") +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  geom_text(aes(label = paste0(round(percentage_system_genome, 2), "%")),                                   # round it to 2 digit
            position = position_stack(vjust = 0.5), size = 7) +
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        legend.box = "horizontal",
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.spacing.x = unit(0.5, "cm"))

cas3_pie_plot

ggsave("../final_figures_and_tables/cas3_pie_plot.png", bg="white", plot = cas3_pie_plot, width = 5.4, height = 3.5, dpi = 600)

###########################################
