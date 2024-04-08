##########################################
# 
# R programming to analyze both PADLOC and CRISPRCasTyper (cctyper) output result
# Perform MVA and plot PCA
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
fix_header <- colnames(read.csv("..\data\merged_crisprs_near_cas.tab", sep="\t", header=T) %>% as_tibble()) #last header falls out in the blank column

# Read the crispr file 
crispr_data <- read.csv("..\data\merged_crisprs_near_cas.tab", sep="\t", header=T) %>% as_tibble() %>%
  subset(!grepl("AccessionID", AccessionID)) %>%
  separate(col = "AccessionID", into = c("AccessionID", "Contig1"), sep = " ") %>% #accession and contigs were in same column1: seperate
  setNames(fix_header) %>% 
  select(-ncol(.)) %>% # remove last column due to header mismatch  
  separate(col = "AccessionID", into = c("AccessionID", "version"), sep = "\\.") %>%                                                                                       # due to absence of version in some will receive warning: Expected 2 pieces. Missing pieces filled with `NA` in 244 rows [107, 108, ...
  filter(Subtype_probability > 0.75) %>% # number could be higher due to the presence of !bacteria or 231
  rename(system = Subtype) %>% 
  select(AccessionID, system)     # warning due to version

# Fix the header
#fix_header <- colnames(crispr_data) #last header falls out in the blank column

# Append padloc_data0 and crispr_data0 and metadata of genomes and metadata of defence system
data0 <- rbind(padloc_data, crispr_data) %>%
  group_by(AccessionID, system) %>%
  summarise(Frequency = n()) %>%
  ungroup() %>% 
  full_join(., metadata_genomes, join_by(AccessionID == accession_genbank)) %>% 
  filter(classification != "#N/A") %>%                                                                              # already removed archaea turns into #N/A
  filter(!grepl("d__Archaea", classification)) %>%                                                                  # already removed archaea
  left_join(., bact_defence_combined, join_by(system == system_type)) %>% 
  mutate(bact_defence_system = recode(bact_defence_system, "CRISPR-Cas_cctyper" = "CRISPR-Cas"))


# filter data to crate matrix
system_unitary_matrix <- data0 %>% select(AccessionID, bact_defence_system, Frequency, Phylum, Phylum_count, Family, Family_count) %>% 
  mutate(Frequency = as.integer(Frequency > 0)) %>%                                                                  # convert to unitary matrix
  distinct() %>% 
  pivot_wider(names_from = bact_defence_system, values_from = Frequency, values_fill = 0) %>% 
  select(-`NA`, -`CRISPR-Cas_padloc`, -DMS) 
  #select(-AccessionID, -Phylum, -Family)

################################################# TRY 1 keep unitary matrix, to match heatmap select family genome count > 5, 
fum <- system_unitary_matrix %>% 
  separate(col = "Family_count", into = c("Family", "F_count"), sep = "\n") %>% 
  filter(as.integer(F_count) > 5) %>%                                                                               # filter Family to match heatmap select family genome count > 5
  select(-AccessionID, -Phylum, -Phylum_count, -F_count) %>% 
  mutate(system_sum = rowSums(select(., -1))) %>%                                                                   # 1 sum1 all system for a family row or remove zero row
  bind_rows(summarise(.data = ., across(where(is.numeric), sum, na.rm = TRUE), Family = "fam_sum"))  # 2 sum2 all system individually, will through error with across()

fum_filter <- fum %>% 
  filter(as.integer(system_sum) > 0) %>%                                                                            # filter sum1, or remove row/genome with zero value
  select(where(~last(.) > 5)) %>%                                                                                   # filter sum2, or remove column/system with more 5 value
  select(-system_sum) %>% 
  filter(Family != "fam_sum") %>%
  mutate(Family = factor(Family))
  

x <- prcomp(fum_filter[-1], center = TRUE)
summary(x)


sys.pca.f <- PCA(fum_filter[-1], scale.unit = FALSE)  # scaling not necessary as same unit, scale.unit = FALSE, graph = TRUE

systerm5_scree_plot.f <- fviz_eig(sys.pca.f, addlabels = TRUE)

ggsave("../final_figures_and_tables/systerm5_scree_plot_f.png", bg="white", plot = systerm5_scree_plot.f, width = 6.25, height = 4.6, dpi = 600)


systerm6_pca_family_plot <- fviz_pca_var(sys.pca.f, col.var = "cos2", # col.var = "contrib", #col.var = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.02), title = "") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
             #select.var = list(contrib = 10))

ggsave("../final_figures_and_tables/systerm6_pca_family_plot.png", bg="white", plot = systerm6_pca_family_plot, width = 6.25, height = 4.6, dpi = 600)


#################### didn't worked ##########
fviz_pca_ind(sys.pca, pointshape = 21, pointsize = 2.5,#pointsize = "cos2",#label="none", #select.ind = list(cos2 = 100),
             geom.ind = "point", # show points only (nbut not "text")
             fill.ind = fum_filter$Family, #col.var = factor(c("Flavobacteriaceae", "Alteromonadaceae", "Burkholderiaceae", "Planococcaceae")),
             #habillage = fum_filter$Family, 
             #addEllipses=TRUE, ellipse.level=0.95,
             #col.ind = fum_filter$Family, # color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             #addEllipses = TRUE, # Concentration ellipses
             #addEllipses = TRUE, ellipse.type = "convex", col.ind = fum_filter$Family,
             legend.title = "Groups"
             #,select.ind = list(contrib = 30)
)


fviz_pca_biplot(sys.pca, pointshape = 21, pointsize = 2.5,#pointsize = "cos2",
                geom.ind = "point", 
                col.var = "cos2", 
                fill.ind = fum_filter$Family, #col.var = factor(c("Flavobacteriaceae", "Alteromonadaceae", "Burkholderiaceae", "Planococcaceae")),
                #,select.ind = list(contrib = 30)
                label = "var",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE,
                select.ind = list(cos2 = 0.6),
                select.var = list(cos2 = 0.02),
                legend.title = list(title = "Groups", fill = "Family")
                #,select.ind = list(cos2 = 10)
                
)

##################################################################################################################
##################################################################################################################
################################################# TRY 2 keep unitary matrix, to match heatmap select phylum genome count > 5, 
phy <- system_unitary_matrix %>% 
  separate(col = "Phylum_count", into = c("Phylum", "P_count"), sep = "\n") %>% 
  filter(as.integer(P_count) > 5) %>%                                                                               # filter Family to match heatmap select family genome count > 5
  select(-AccessionID, -Family, -Family_count, -P_count) %>% 
  mutate(system_sum = rowSums(select(., -1))) %>%                                                                   # 1 sum1 all system for a family row or remove zero row
  bind_rows(summarise(.data = ., across(where(is.numeric), sum, na.rm = TRUE), Phylum = "phy_sum"))  # 2 sum2 all system individually, will through error with across()

phy_filter <- phy %>% 
  filter(as.integer(system_sum) > 0) %>%                                                                            # filter sum1, or remove row/genome with zero value
  select(where(~last(.) > 5)) %>%                                                                                   # filter sum2, or remove column/system with more 5 value
  select(-system_sum) %>% 
  filter(Phylum != "phy_sum") %>%
  mutate(Phylum = factor(Phylum))


y <- prcomp(phy_filter[-1], center = TRUE)
summary(y)

sys.pca.p <- PCA(phy_filter[-1], scale.unit = FALSE)  # scaling not necessary as same unit, scale.unit = FALSE, graph = TRUE

systerm7_scree_plot.p <- fviz_eig(sys.pca.p, addlabels = TRUE)
systerm7_scree_plot.p

systerm8_pca_phylum_plot <- fviz_pca_var(sys.pca.p, col.var = "cos2", # col.var = "contrib", #col.var = "cos2", 
                                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                         repel = TRUE,
                                         select.var = list(cos2 = 0.02), title = "") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
systerm8_pca_phylum_plot


################################################# TRY 3 keep unitary matrix, to match heatmap select phylum genome count > 5, 
#Q1 Genomes with CRISPR present, which other system present 

phy_filter_cas_p <- phy %>% 
  filter(as.integer(system_sum) > 0) %>%                                                                            # filter sum1, or remove row/genome with zero value
  select(where(~last(.) > 5)) %>%                                                                                   # filter sum2, or remove column/system with more 5 value
  select(-system_sum) %>% 
  filter(Phylum != "phy_sum") %>%
  mutate(Phylum = factor(Phylum)) %>% 
  filter(`CRISPR-Cas` == "1") #%>% 
  #select(-`CRISPR-Cas`)

fum_filter_cas_p <- fum %>% 
  filter(as.integer(system_sum) > 0) %>%                                                                            # filter sum1, or remove row/genome with zero value
  select(where(~last(.) > 5)) %>%                                                                                   # filter sum2, or remove column/system with more 5 value
  select(-system_sum) %>% 
  filter(Family != "fam_sum") %>%
  mutate(Family = factor(Family)) %>% 
  filter(`CRISPR-Cas` == "1") #%>% 
  #select(-`CRISPR-Cas`)

z1 <- prcomp(phy_filter_cas_p[-1], center = TRUE)
summary(z1)

z2 <- prcomp(fum_filter_cas_p[-1], center = TRUE)
summary(z2)

sys.pca.p_p <- PCA(phy_filter_cas_p[-1], scale.unit = FALSE)  # scaling not necessary as same unit, scale.unit = FALSE, graph = TRUE
sys.pca.f_p <- PCA(fum_filter_cas_p[-1], scale.unit = FALSE)  # scaling not necessary as same unit, scale.unit = FALSE, graph = TRUE

fviz_eig(sys.pca.p_p, addlabels = TRUE)
fviz_eig(sys.pca.f_p, addlabels = TRUE)

fviz_pca_var(sys.pca.p_p, col.var = "cos2", # col.var = "contrib", #col.var = "cos2", 
                                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                         repel = TRUE,
                                         select.var = list(cos2 = 0.02), title = "") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 

fviz_pca_var(sys.pca.f_p, col.var = "cos2", # col.var = "contrib", #col.var = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.02), title = "") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 


################################################# TRY 4 keep unitary matrix, to match heatmap select phylum genome count > 5, 
#Q2 Genomes when CRISPR not present, which system present

phy_filter_cas_a <- phy %>% 
  filter(as.integer(system_sum) > 0) %>%                                                                                      # filter sum1, or remove row/genome with zero value
  select(where(~last(.) > 5)) %>%                                                                                             # filter sum2, or remove column/system with more 5 value
  select(-system_sum) %>% 
  filter(Phylum != "phy_sum") %>%
  mutate(Phylum = factor(Phylum)) %>% 
  filter(`CRISPR-Cas` != "1") #%>% 
#select(-`CRISPR-Cas`)

fum_filter_cas_a <- fum %>% 
  filter(as.integer(system_sum) > 0) %>%                                                                                     # filter sum1, or remove row/genome with zero value
  select(where(~last(.) > 5)) %>%                                                                                            # filter sum2, or remove column/system with more 5 value
  select(-system_sum) %>% 
  filter(Family != "fam_sum") %>%
  mutate(Family = factor(Family)) %>% 
  filter(`CRISPR-Cas` != "1") #%>% 
#select(-`CRISPR-Cas`)

z3 <- prcomp(phy_filter_cas_a[-1], center = TRUE)
summary(z3)

z4 <- prcomp(fum_filter_cas_a[-1], center = TRUE)
summary(z4)

sys.pca.p_a <- PCA(phy_filter_cas_a[-1], scale.unit = FALSE)  # scaling not necessary as same unit, scale.unit = FALSE, graph = TRUE
sys.pca.f_a <- PCA(fum_filter_cas_a[-1], scale.unit = FALSE)  # scaling not necessary as same unit, scale.unit = FALSE, graph = TRUE

fviz_eig(sys.pca.p_a, addlabels = TRUE)
fviz_eig(sys.pca.f_a, addlabels = TRUE)

fviz_pca_var(sys.pca.p_a, col.var = "cos2", # col.var = "contrib", #col.var = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.02), title = "") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 

fviz_pca_var(sys.pca.f_a, col.var = "cos2", # col.var = "contrib", #col.var = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             select.var = list(cos2 = 0.02), title = "") +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 

################################################# TRY 5 keep unitary matrix, to match heatmap select phylum genome count > 5, 
# Dataset with binary variables. Dependent variable is binary (0 or 1), logistic regression would be more appropriate than linear regression.
phy_filter$`CRISPR-Cas` <- as.factor(phy_filter$`CRISPR-Cas`)

# Specify the dependent variable (response variable)
dependent_variable <- "`CRISPR-Cas`"

# Extract the columns of interest
independent_variables <- colnames(phy_filter)[colnames(phy_filter) != dependent_variable & colnames(phy_filter) != "Phylum"] 

# Create an empty data frame to store the results
p_results <- data.frame(variable = character(), coef = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through independent variables and fit logistic regression models
for (variable in independent_variables) {
  # Fit logistic regression model
  model <- glm(paste(dependent_variable, "~", variable), data = phy_filter, family = "binomial")
  
  # Extract coefficients and p-value
  coef_value <- coef(summary(model))[2, "Estimate"]
  p_value <- coef(summary(model))[2, "Pr(>|z|)"]
  
  # Store results in the data frame
  result_row <- data.frame(variable = variable, coef = coef_value, p_value = p_value, stringsAsFactors = FALSE)
  p_results <- rbind(p_results, result_row)
}

# View the results
print(p_results)



