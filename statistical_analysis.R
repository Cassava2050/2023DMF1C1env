# Data analysis
rm(list = ls())

## Load libraries
library(asreml)
source("https://raw.githubusercontent.com/Cassava2050/PPD/main/utilities_tidy.R") # become into a 
trial_interest <- "MDEPR"
year_interest <- 2023


# master_data to save the results
master_data <- list()

# Load the tidy data
trial_set_number = 1

# all files in the folder
list_file = list.files(here::here("output"))

# tidy data of the trials interested
sel_file = list_file[str_detect(list_file, "_tidy_data4analysis_")] 

# the data we will use
sel_file_use = sel_file[1]

sel_file_use

trial1_tidy = read.csv(here::here("output", sel_file_use), header=TRUE,
                       stringsAsFactors = FALSE,
                       as.is=T,
                       check.names = FALSE)
if(trial_set_number == 1){
  trial_tidy_all = trial1_tidy
}

# Obtain all the trait information using a cloud file (gitHub) -------
trait_all <-
  read.csv("https://raw.githubusercontent.com/lfdelgadom/standar_col_names_CB/main/standar_col_names.csv") %>%
  select(analysis_col_name) %>%
  filter(str_detect(analysis_col_name, "obs_"))
trait_all_adj <- gsub("obs_", "", trait_all$analysis_col_name)
trait_all_adj = c(trait_all_adj,
                  "harvest_number_plan", "germination_perc",
                  "yield_ha_v2", "DM_yield_ha", "starch_content", "starch_yield_ha")
trait_all_adj <- gsub("-", "_", trait_all_adj)

# Meta info.
meta_all <-
  read.csv("https://raw.githubusercontent.com/lfdelgadom/standar_col_names_CB/main/standar_col_names.csv") %>%
  select(analysis_col_name) %>%
  filter(str_detect(analysis_col_name, "use_"))
meta_all_adj <- gsub("use_", "", meta_all$analysis_col_name)
meta_all_adj <- c(
  meta_all_adj,
  "check_released", "latitude", "longitude",
  "altitude", "department", "country",
  "ag_zone", "location_short"
)

# Select the observations for analysis
names(trial_tidy_all) <- gsub("-", "_", names(trial_tidy_all))
analysis_trait <- names(trial_tidy_all)[names(trial_tidy_all) %in% trait_all_adj]
print("All the traits investigated:")
print(analysis_trait)

# Select the meta information for analysis
meta_col <- names(trial_tidy_all)[names(trial_tidy_all) %in% meta_all_adj]
print("All the meta information:")
print(meta_col)

# Check the SD of each trait
trial_rm_sd <- remove_no_var_tidy(my_dat = trial_tidy_all,
                                  analysis_trait = analysis_trait,
                                  meta_info = meta_col)
master_data[["mean_of_sd"]] = sd_mean

# Trait ideal
no_traits_for_analysis <- c(
  "stake_plant" , "planted_number_plot", 
  "harvest_number", "root_weight_air", 
  "root_weight_water", "harvest_number_plan",
  "yield_ha_v2", "root_rot_perc", "harvest_index",
  "germinated_number_plot"
)

no_variation_traits <- c() # "CAD_5mon", "CAD_7mon", "CAD_3mon", "lodging1_3_6mon"

no_traits_for_analysis <- c(no_variation_traits, no_traits_for_analysis)

trait_ideal <- analysis_trait[!analysis_trait %in% no_traits_for_analysis]
print("the trait ideal is:"); trait_ideal

trait_ideal %>% as.data.frame() %>% 
  write.table("clipboard", sep = "\t", col.names = T, row.names = F)


# Genotypic correlation (Phenotypic values)
correlation <- gg_cor(
  colours = c("red", "white", "blue"),
  data = trial_rm_sd[, trait_ideal],
  label_size = 2
)

ggsave(paste("images\\pheno_corr", trial_interest, Sys.Date(), ".png", sep = "_"),
       plot = correlation, units = "in", dpi = 300, width = 12, height = 8
)

# Check design experimental
my_dat <- trial_rm_sd %>% 
  add_column(block = NA) %>% mutate(block = as.factor(block)) 

# number of trials
length(unique(my_dat$trial_name)) 

results <- check_design_met(
  data = my_dat,
  genotype = "accession_name",
  trial = "trial_name",
  traits = trait_ideal,
  rep = "rep_number",
  col = "col_number",
  row = "row_number",
  block = "block"
)

summary <- results$summ_traits 

p1 <- summary %>% 
  ggplot(aes(x = traits , y = trial_name, label = round(miss_perc,2),  fill = miss_perc ))+
  geom_tile(color = "gray")+
  geom_text(color = "white")+
  theme_minimal(base_size = 13)+
  labs(title = "Percentage of missing values (exp/trait)", x = "", y = "") +
  theme(axis.text.x = element_text(hjust = 1 , angle = 75, size = 16),
        axis.text.y = element_text(size = 16))
p1
ggsave(paste("images\\missing_", trial_interest, Sys.Date(), ".png", sep = "_"),
       plot = p1, units = "in", dpi = 300, width = 15, height = 6
)
master_data[["summ_traits"]] <- summary


## Single trial analysis
obj <- single_trial_analysis(results = results,
                             progress = TRUE,
                             remove_outliers = FALSE)


trials <- unique(my_dat$trial_name)

header_sort = vector()
i = 1
for (i in 1:length(trials)) {
  
  cat("\n_______________")
  cat("\nTRIAL:", trials[i], "\n")
  cat("_______________\n")
  
  for (j in 1:length(trait_ideal)) {
    
    blue_blup <- obj$blues_blups %>% 
      filter(trial == trials[i]) %>% 
      select(-c(trial, seBLUEs, seBLUPs, wt)) %>% 
      pivot_wider(names_from = "trait", values_from = c("BLUEs", "BLUPs"))
    
    header_sort = c(header_sort,
                    grep(trait_ideal[j], sort(names(blue_blup)), value=TRUE))
    blue_blup <- blue_blup %>% dplyr::select(genotype, any_of(header_sort)) %>% 
      mutate(across(where(is.double), round, 1))
  }
  master_data[[paste0("BLUP_BLUE_", trials[i])]] <- blue_blup
}

# Save the spatial correction plots
folder = paste0(here::here("output"), "/")
pdf(paste(folder, "01_", trial_interest, "_spatial_correction_", Sys.Date(), 
          ".pdf", sep = ""), width = 8, height = 6)
plot(obj, type = "spatial") 
dev.off()


## Single heritability
single_h2 <- obj$resum_fitted_model[ ,1:3] %>% 
  group_by(trial) %>%
  spread(trait, value = heritability) #%>% print(width = Inf) 

single_h2 %>% print(width = Inf)

master_data[["single_h2"]] <- single_h2 

single_h2 %>% 
  write.table("clipboard", sep = "\t", col.names = T, row.names = F, na = "")

# ---------------------------------------------------------------------------

# BLUPs gxe
BLUPs_table <- 
  blue_blup %>% 
  select(genotype, starts_with("BLUPs")) %>% 
  rename("accession_name" = genotype) %>% 
  mutate(across(where(is.numeric), round, 2))

#save the BLUPs data
master_data[[paste0("BLUPs_")]] <- BLUPs_table

# genotypic correlation
geno_cor <- gg_cor(
  colours = c("red", "white", "blue"),
  data = BLUPs_table, 
  label_size = 2.5
) + 
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14))

geno_cor
# save corr plot
ggsave(paste("images\\geno_corr", trial_interest, Sys.Date(), ".png", sep = "_"),
       units = "in", dpi = 300, width = 14, height = 8)

## Save the master data results
folder_output <- here::here("output//")
meta_file_name <- paste0(folder_output, paste(year_interest, trial_interest, "master_results", Sys.Date(), ".xlsx", sep = "_"))
write.xlsx(master_data, file = meta_file_name)


## Index selection
list_file <- list.files(folder_output)
sel_file <- list_file[str_detect(list_file, "_master_results_") &
                        str_detect(list_file, trial_interest)]
sel_file

sel_file[1]
blupDF_kp <- read_excel(
  paste(folder_output,
        sel_file[1],
        sep = ""
  ),
  sheet = "BLUP_BLUE_202307DMF1C_tani"
)


## Selection index
colnames(blupDF_kp)

index_traits <- c("BLUPs_starch_content", 
                  "BLUPs_yield_ha", "BLUPs_CMD_harvest")

index_dat <- blupDF_kp %>%
  select("genotype", all_of(index_traits)) %>% 
  drop_na()


# Selection index function
# multi-trait -------------------------------------------------------------
library(FactoMineR)
library(factoextra)

pca_index <- function(data, id, variables = NULL, percentage = 0.20, b) {
  # The data set to be analyzed. It should be in the form of a data frame.
  data <- as.data.frame(data)
  rownames(data) <- data[, id]
  if (is.null(variables)) variables <- names(data)[names(data) != id]
  data <- data[, variables]
  index <- selIndex(Y = as.matrix(data), b = b, scale = T)
  index <- c(index)
  data$index <- index
  data <- data %>% arrange(desc(index))
  data$selected <- NA
  data$selected[1:(round(percentage * nrow(data)))] <- TRUE # select best genos (larger index selection)
  data$selected <- ifelse(is.na(data$selected), FALSE, data$selected) # reject genos with lower index selection
  res.pca <- PCA(data, graph = T, scale.unit = T, quali.sup = ncol(data))
  
  final <- fviz_pca_biplot(res.pca,
                           habillage = data$selected,
                           geom = c("point"),
                           addEllipses = T,
                           col.var = "black",
                           ggtheme = theme_minimal()
  )
  
  
  selection <- data %>% filter(selected == T)
  return(list(res.pca = res.pca, final = final, results = data, selection = selection))
}

selIndex <- function(Y, b, scale = FALSE) {
  if (scale) {
    return(scale(Y) %*% b)
  }
  return(Y %*% b)
}


## Index selection
res.pca <- pca_index(
  data = index_dat, id = "genotype",
  variables = index_traits,
  b = c(10, 10, -5), percentage = 0.25
)
res.pca_final <- res.pca$final
res.pca_final
ggsave(paste("images/selection", Sys.Date(), ".png"),
       plot = res.pca_final, units = "in", dpi = 300, width = 7, height = 6
)
res.pca$selection
selections <- res.pca$results %>% rownames_to_column(var = "accession_name")


selections %>% 
  select(accession_name, index, everything()) %>% 
  write.table("clipboard", sep = "\t", col.names = T, row.names = F)


# Add index column to BLUEs_BLUPs_MET
blue_blup <- 
  master_data$BLUP_BLUE_202307DMF1C_tani %>% 
  left_join(selections[-c(2:4)], by = c("genotype" = "accession_name")) %>% 
  relocate(index, selected, .before = 2)

blue_blup <- blue_blup %>% 
  arrange(is.na(selected)) %>% select(genotype, index, selected, 
                                      contains("starch"), contains("yield_ha"), everything())
master_data[["BLUP_BLUE_202307DMF1C_tani"]] = blue_blup


## Save the master data results
folder_output <- here::here("output//")
meta_file_name <- paste0(folder_output, paste(year_interest, trial_interest, "master_results", Sys.Date(), ".xlsx", sep = "_"))
write.xlsx(master_data, file = meta_file_name)

