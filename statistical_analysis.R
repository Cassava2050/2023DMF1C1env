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


# Rep_number == 3 was removed from 2021103MDAYT_quan, 
# due to there were several missing values
# my_dat_1 <- trial_rm_sd %>% filter(trial_name 
#                                    %in% c("2021103MDAYT_quan"), 
#                                    !rep_number == 3)

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

# shared clones
shared <- plot(results, type = "connectivity")
ggsave(paste('images\\shared_', trial_interest, Sys.Date(), ".png", sep = "_"),
       plot = shared, units = "in", dpi = 300, width = 8, height = 6)

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

plot(obj, type = "spatial") 

## Single heritability
single_h2 <- obj$resum_fitted_model[ ,1:3] %>% 
  group_by(trial) %>%
  spread(trait, value = heritability) #%>% print(width = Inf) 

single_h2 %>% print(width = Inf)

master_data[["single_h2"]] <- single_h2 

single_h2 %>% 
  write.table("clipboard", sep = "\t", col.names = T, row.names = F, na = "")

# ---------------------------------------------------------------------------
traits_to_remove <- c("root_number_commercial",
                      "root_peduncle1_3",
                      "root_rot_number",
                      "root_shape1_6",	
                      "root_type1_5",
                      "root_weight_plot",
                      "starch_content",
                      "starch_yield_ha",
                      "vigor1_5",
                      "yield_ha",
                      "carotenoid1_8",
                      "CMD_1mon",
                      "CMD_6mon",
                      "CMD_harvest",
                      "germination_perc")

met_results <- met_analysis(obj, 
                            filter_traits = 
                            trait_ideal[!trait_ideal %in% c(traits_to_remove)],
                            h2_filter = 0.09,
                            progress = TRUE
)


# h2 gxe
master_data[["h2_gxe"]] <- 
  met_results$heritability %>% 
  arrange(desc(h2)) %>%
  mutate(across(where(is.numeric), round, 2))

master_data$h2_gxe %>%
  write.table("clipboard", col.names = T, row.names = F, sep = "\t")

# BLUPs gxe
BLUPs_table <- 
  met_results$overall_BLUPs %>% 
  select(-c(std.error, status)) %>% 
  group_by(genotype) %>% 
  spread(trait, value = predicted.value) %>% 
  rename("accession_name" = genotype) %>% 
  mutate(across(where(is.numeric), round, 2)) %>% 
  ungroup() 
#save the BLUPs data
master_data[[paste0("BLUPs_", "gxe")]] <- BLUPs_table

library(corrplot)

# Genotypic Correlation: Locations
# CMD_3mon
# Open a PNG device
png(paste0("images\\corrplot_CMD_3mon", trial_interest, Sys.Date(), ".png"), 
    units = "in", res = 300,
    width = 8, height = 10)

# Create the correlation plot
corrplot(met_results$VCOV$CMD_3mon$CORR, method = "color",  
         type = "lower", order = "hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black",
         tl.cex = 1.5, # Text label color and size,
         number.cex = 1.5,
         diag = TRUE)

# Close the device
dev.off()

# branch_number
png(paste0("images\\corrplot_branch_number", trial_interest, Sys.Date(), ".png"), 
    units = "in", res = 300,
    width = 8, height = 10)

# Create the correlation plot
corrplot(met_results$VCOV$branch_number$CORR, method = "color",  
         type = "lower", order = "hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black",
         tl.cex = 1.5, # Text label color and size,
         number.cex = 1.5,
         diag = TRUE)

# Close the device
dev.off()

# plant_type
png(paste0("images\\corrplot_plant_type", trial_interest, Sys.Date(), ".png"), 
    units = "in", res = 300,
    width = 8, height = 10)

# Create the correlation plot
corrplot(met_results$VCOV$plant_type$CORR, method = "color",  
         type = "lower", order = "hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black",
         tl.cex = 1.5, # Text label color and size,
         number.cex = 1.5,
         diag = TRUE)

# Close the device
dev.off()


# height_1st_branch
png(paste0("images\\corrplot_height_1st_branch", trial_interest, Sys.Date(), ".png"), 
    units = "in", res = 300,
    width = 8, height = 10)

# Create the correlation plot
corrplot(met_results$VCOV$height_1st_branch$CORR, method = "color",  
         type = "lower", order = "hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black",
         tl.cex = 1.5, # Text label color and size,
         number.cex = 1.5,
         diag = TRUE)

# Close the device
dev.off()



## Save variance covariance correlation

as.data.frame(do.call(rbind, met_results$VCOV))$CORR

## Save the BLUEs or raw data across the trials
variables <- colnames(BLUPs_table)[!grepl("accession_name", 
                                          colnames(BLUPs_table))]
for (var in variables) {
  
  cat("\n_______________")
  cat("\nTRIAL:", var, "\n")
  cat("_______________\n")
  
  blue_blup <-
    obj$blues_blups %>%
    select(trial, genotype, trait, BLUEs) %>%
    spread(trait, value = BLUEs) %>%
    select(trial, genotype, any_of(var)) %>%
    group_by(trial, genotype) %>%
    pivot_wider(names_from = trial, values_from = any_of(var)) %>%
    right_join(BLUPs_table %>%
                 select(accession_name, any_of(var)), by = c("genotype" = "accession_name")) %>%
    arrange(is.na(across(where(is.numeric))), across(where(is.numeric))) %>%
    mutate(across(where(is.numeric), round, 2))
  # remove all NA columns
  blue_blup <- blue_blup[, colSums(is.na(blue_blup)) < nrow(blue_blup)]
  
  master_data[[paste0("BLUP_BLUE_", var)]] <- blue_blup
}

## Stability analysis
for (var in variables) {
  
  cat("\n_______________")
  cat("\nTRIAL:", var, "\n")
  cat("_______________\n")
  
  stab <- met_results$stability %>% 
    filter(trait == var) %>% 
    arrange(superiority) %>% 
    pivot_wider(names_from = "trait", values_from = c('predicted.value')) 
  
  # Change colname
  colnames(stab)[5] <- paste('BLUPs', colnames(stab)[5], sep = '_') 
  colnames(stab)[c(2, 3, 4)] <- paste(colnames(stab)[c(2, 3, 4)], var, sep = '_') 
  
  master_data[[paste0("stability_", var)]] <- stab
}


ind <- grep("^stability_", names(master_data))

# select elements that met the condition
stab_values <- master_data[ind] %>% 
  reduce(inner_join, by = "genotype") %>% 
  select(!starts_with("BLUPs_")) %>% 
  mutate(across(where(is.numeric), round, 2))

# remove multiple stability sheets
master_data[ind] <- NULL

## BLUE and BLUP data together
BLUES_dona <- c("root_number_commercial",	"root_peduncle1_3",
               "root_rot_number", "root_shape1_6", "root_type1_5", 
               "root_weight_plot", "starch_content", "starch_yield_ha",	
               "yield_ha")

BLUEs_BLUPs <- 
  obj$blues_blups %>%
  select(trait, genotype, trial, BLUEs, seBLUEs) %>%
  filter(trait %in% variables) %>% 
  pivot_wider(names_from = "trait", values_from = c("BLUEs", "seBLUEs")) %>%
  pivot_wider(names_from = trial, values_from = c(
    paste("BLUEs", variables, sep = "_"),
    paste("seBLUEs", variables, sep = "_")
  )) %>%
  left_join(
    met_results$overall_BLUPs %>%
      select(!status) %>%
      rename(
        BLUPs = predicted.value,
        seBLUPs = std.error
      ) %>%
      pivot_wider(names_from = "trait", values_from = c("BLUPs", "seBLUPs")),
    by = "genotype"
  ) %>%
  #arrange(desc(BLUPs_starch_content)) %>% 
  arrange(is.na(across(where(is.numeric))), across(where(is.numeric))) %>%
  mutate(across(where(is.numeric), round, 2))
# remove all NA columns
BLUEs_BLUPs <- BLUEs_BLUPs[, colSums(is.na(BLUEs_BLUPs)) < nrow(BLUEs_BLUPs)]


# put all together stab_values with blues_blups
BLUEs_BLUPs <- 
  BLUEs_BLUPs %>% left_join(stab_values, by = 'genotype')  


# add BLUES of dona location
BLUES_dona_value <- obj$blues_blups %>% 
  select(trait, genotype, trial, BLUEs, seBLUEs) %>% 
  filter(trial == "202304DMPYT_dona", trait %in% BLUES_dona)

BLUES_dona_wider <- BLUES_dona_value %>% 
  pivot_wider(names_from = "trait", values_from = c("BLUEs", "seBLUEs")) %>%
  pivot_wider(names_from = trial, values_from = c(
    paste("BLUEs", BLUES_dona, sep = "_"),
    paste("seBLUEs", BLUES_dona, sep = "_")
  )) 

BLUEs_BLUPs <- BLUEs_BLUPs %>% left_join(BLUES_dona_wider, by = "genotype")

variables_BLUE_dona <- c(variables , BLUES_dona)
header_sort = vector()
for (i in 1:length(variables_BLUE_dona)) {
  
  header_sort = c(header_sort, 
                  grep(variables_BLUE_dona[i], sort(names(BLUEs_BLUPs)), value=TRUE) 
  )
  
}


BLUEs_BLUPs <- BLUEs_BLUPs %>%
  select(genotype, all_of(header_sort), -starts_with("se")) 


# I need to add the BLUES of dona 
BLUEs_BLUPs <- BLUEs_BLUPs %>% 
  relocate(colnames(BLUEs_BLUPs)[str_detect(colnames(BLUEs_BLUPs), "starch_content")], .after = genotype)


master_data[["BLUEs_BLUPs_MET"]] = BLUEs_BLUPs


# Genotypic correlation
# add the BLUES_dona_value  to BLUPs_table 
BLUES_dona_BLUES <- BLUES_dona_value %>% pivot_wider(names_from = "trait", 
                            values_from = c("BLUEs", "seBLUEs")) %>% 
  select(genotype, starts_with("BLUEs"))

# remove BLUES word from BLUES_dona_BLUES colnames
colnames(BLUES_dona_BLUES)[-1] <- gsub("BLUEs_", "", names(BLUES_dona_BLUES)[-1])

# merging both tables
BLUPs_table <- BLUPs_table %>% left_join(BLUES_dona_BLUES, by = c("accession_name" = "genotype"))


# save again the updated BLUPs table
#save the BLUPs data
master_data[[paste0("BLUPs_", "gxe")]] <- BLUPs_table

geno_cor <- gg_cor(
  colours = c("red", "white", "blue"),
  data = BLUPs_table, 
  label_size = 2.5
) + 
  theme(
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14))

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
  sheet = paste0("BLUPs_", "gxe")
)


## Selection index
colnames(blupDF_kp)

index_traits <- c("starch_content", #"height_1st_branch", 
                  "yield_ha", "plant_type")

index_dat <- blupDF_kp %>%
  select("accession_name", all_of(index_traits)) %>% 
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
  data = index_dat, id = "accession_name",
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
BLUEs_BLUPs <- 
  master_data$BLUEs_BLUPs_MET %>% 
  left_join(selections[-c(2:4)], by = c("genotype" = "accession_name")) %>% 
  relocate(index, selected, .before = 2)

BLUEs_BLUPs <- BLUEs_BLUPs %>% 
  arrange(is.na(selected))
master_data[["BLUEs_BLUPs_MET"]] = BLUEs_BLUPs


## Save the master data results
folder_output <- here::here("output//")
meta_file_name <- paste0(folder_output, paste(year_interest, trial_interest, "master_results", Sys.Date(), ".xlsx", sep = "_"))
write.xlsx(master_data, file = meta_file_name)

