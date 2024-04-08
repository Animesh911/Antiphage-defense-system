##########################################
# 
# R programming to analyze both PADLOC and CRISPRCasTyper (cctyper) output result
# Plot PCA
# 
# Author: Animesh Kumar
#
###########################################

library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(factoextra)

# Read a tab-separated metadata file of defence system
bact_defence_combined <- read.table("..\data\metadata_bact_defence_combined.txt", sep = "\t", header = TRUE)
rownames(bact_defence_combined) <- bact_defence_combined[,1]                    # both are CRISPRcas I-F and I-F_T

#Read metadata of genomes
metadata_genomes <- read.csv("..\data\metadata_genomes.txt", sep="\t", header=T) %>% as_tibble() %>% 
  distinct(accession_genbank, .keep_all = TRUE) %>% 
  filter(classification != "#N/A") %>%                                          # filter #N/A
  filter(!grepl("d__Archaea", classification)) %>%                              # filter archaea
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

######################## Metadata genomes S1
library("xlsx")

metadata_genomes %>% 
  select(accession_genbank, database, isol_env, host_species, dependent, oxygen_requirement, classification) %>% 
  as.data.frame() %>% 
  write.xlsx(file = "../final_figures_and_tables/Table_supplementary_publication.xlsx", sheetName = "tableS1_metadata_genomes", append = FALSE, row.names = FALSE)

########################
# Read the PADLOC output file
padloc_data <- read.csv("..\data\all_padloc.tab") %>% as_tibble() %>% 
  filter(full.seq.E.value < 0.01) %>% 
  filter(domain.iE.value < 0.01) %>% 
  filter(target.coverage > 0.8) %>%
  filter(hmm.coverage > 0.8) %>%
  distinct(AccessionID, system, system.number) %>% 
  select(AccessionID, system)

#################### table S3
# pAgo proteins in cold-adapted bacteria
read.csv("all_padloc.csv") %>% as_tibble() %>% 
  filter(full.seq.E.value < 0.01) %>% 
  filter(domain.iE.value < 0.01) %>% 
  filter(target.coverage > 0.8) %>%
  filter(hmm.coverage > 0.8) %>% 
  filter(grepl("pAgo_(I|II)", protein.name) | protein.name == "pAgo") %>% 
  distinct(AccessionID, system, system.number, .keep_all=TRUE) %>%
  left_join(., metadata_genomes[, c("accession_genbank", "classification")], join_by(AccessionID == accession_genbank)) %>% 
  select(-gtdb) %>% 
  as.data.frame() %>% 
  write.xlsx(file = "../final_figures_and_tables/Table_supplementary_publication.xlsx", sheetName = "table2_pAgo_protein", append = TRUE, row.names = FALSE)


###################### table S4
# pAgo proteins in cold-adapted bacteria
read.csv("all_padloc.csv") %>% as_tibble() %>% 
  filter(full.seq.E.value < 0.01) %>% 
  filter(domain.iE.value < 0.01) %>% 
  filter(target.coverage > 0.8) %>%
  filter(hmm.coverage > 0.8) %>% 
  filter(grepl("^retron", system)) %>% 
  filter(grepl("^RT", protein.name)) %>% 
  distinct(AccessionID, system, system.number, .keep_all=TRUE) %>%
  left_join(., metadata_genomes[, c("accession_genbank", "classification")], join_by(AccessionID == accession_genbank)) %>% 
  select(-gtdb) %>% 
  as.data.frame() %>% 
  write.xlsx(file = "../final_figures_and_tables/Table_supplementary_publication.xlsx", sheetName = "table3_retron_system", append = TRUE, row.names = FALSE)


####################
# Fix the header
fix_header <- colnames(read.csv("..\data\merged_crisprs_near_cas.tab", sep="\t", header=T) %>% as_tibble())         #fix column header as last header falls out in the blank column

# Read the crispr result file 
crispr_data <- read.csv("..\data\merged_crisprs_near_cas.tab", sep="\t", header=T) %>% as_tibble() %>%
  subset(!grepl("AccessionID", AccessionID)) %>%
  separate(col = "AccessionID", into = c("AccessionID", "Contig1"), sep = " ") %>%                                   # separate accession and contigs which are in same column1
  setNames(fix_header) %>% 
  select(-ncol(.)) %>% # remove last column due to header mismatch  
  separate(col = "AccessionID", into = c("AccessionID", "version"), sep = "\\.") %>%                                 # due to absence of version in some will receive warning: Expected 2 pieces. Missing pieces filled with `NA` in 244 rows [107, 108, ...
  filter(Subtype_probability > 0.75) %>%                                                                             # number could be higher due to the presence of non prokaryotes (!bacteria or 231)
  rename(system = Subtype) %>% 
  select(AccessionID, system)                                                                                        # warning due to genome version

# Fix the header
#fix_header <- colnames(crispr_data)                                                                                 # fix header, as last header falls out in the blank column

# Append padloc_data0 and crispr_data0 and metadata_genomes and bact_defence_combined
data0 <- rbind(padloc_data, crispr_data) %>%
  group_by(AccessionID, system) %>%
  summarise(Frequency = n()) %>%
  ungroup() %>% 
  full_join(., metadata_genomes, join_by(AccessionID == accession_genbank)) %>% 
  filter(classification != "#N/A") %>%                                                                              # already removed archaea turns into #N/A
  filter(!grepl("d__Archaea", classification)) %>%                                                                  # already removed archaea
  left_join(., bact_defence_combined, join_by(system == system_type)) %>% 
  mutate(bact_defence_system = recode(bact_defence_system, "CRISPR-Cas_cctyper" = "CRISPR-Cas"))


# create unitary matrix
system_unitary_matrix <- data0 %>% select(AccessionID, bact_defence_system, Frequency, Phylum, Phylum_count, Family, Family_count) %>% 
  mutate(Frequency = as.integer(Frequency > 0)) %>%                                                                  # convert to unitary matrix
  distinct() %>% 
  pivot_wider(names_from = bact_defence_system, values_from = Frequency, values_fill = 0) %>% 
  select(-`NA`, -`CRISPR-Cas_padloc`, -DMS) 
  #select(-AccessionID, -Phylum, -Family)

################################################# Plot PCA using unitary matrix; select heatmap family genome count > 5, 
fum <- system_unitary_matrix %>% 
  separate(col = "Family_count", into = c("Family", "F_count"), sep = "\n") %>% 
  filter(as.integer(F_count) > 5) %>%                                                                               # filter Family to match heatmap select family genome count > 5
  select(-AccessionID, -Phylum, -Phylum_count, -F_count) %>% 
  mutate(system_sum = rowSums(select(., -1))) %>%                                                                   # 1 sum1 all system for a family row or remove zero row
  bind_rows(summarise(.data = ., across(where(is.numeric), sum, na.rm = TRUE), Family = "fam_sum"))                 # 2 sum2 all system individually, will through error with across()

fum_filter <- fum %>% 
  filter(as.integer(system_sum) > 0) %>%                                                                            # filter sum1, or remove row/genome with zero value
  select(where(~last(.) > 5)) %>%                                                                                   # filter sum2, or remove column/system with more 5 value
  select(-system_sum) %>% 
  filter(Family != "fam_sum") %>%
  mutate(Family = factor(Family))
  
summary(prcomp(fum_filter[-1], center = TRUE))

sys.pca.f <- PCA(fum_filter[-1], scale.unit = FALSE)                                                                # scaling not necessary as same unit, scale.unit = FALSE, graph = TRUE

# scree plot
systerm5_scree_plot.f <- fviz_eig(sys.pca.f, addlabels = TRUE)
ggsave("../final_figures_and_tables/systerm5_scree_plot_f.png", bg="white", plot = systerm5_scree_plot.f, width = 6.25, height = 4.6, dpi = 600)

# PCA family plot
systerm6_pca_family_plot <- fviz_pca_var(sys.pca.f, col.var = "cos2", # col.var = "contrib", #col.var = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.02), title = "") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
             #select.var = list(contrib = 10))

ggsave("../final_figures_and_tables/systerm6_pca_family_plot.png", bg="white", plot = systerm6_pca_family_plot, width = 6.25, height = 4.6, dpi = 600)
