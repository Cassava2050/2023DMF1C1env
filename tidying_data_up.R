# Clearing the Environment
rm(list = ls())

# Sourcing External R Script
source("https://raw.githubusercontent.com/Cassava2050/PPD/main/utilities_tidy.R")

# Setting Up File Paths and Parameters
folder <- here::here("data//")  
file <- "phenotype.csv"
skip_col <- 3 # Double check the number of columns skipped
trial_interest <- "MDEPR"
year_interest <- 2022

# Loading Data
# The function 'read_cassavabase' should be defined in the sourced script
sel_data <- read_cassavabase(phenotypeFile = paste0(folder, file))

# Standardizing Column Names
# The function 'change_colname' should also be defined in the sourced script
# Replace 'NA' with the desired standard naming convention if necessary
sel_data_kp <- change_colname(sel_data, NA)

# Changing Column Classes
# Create a vector 'obs_col' containing names of columns related to observations
obs_col <- c(
  names(sel_data_kp)[str_detect(names(sel_data_kp), "obs_")],
  "use_rep_number", "blockNumber",
  "use_plot_number", "use_plot_width",
  "use_plot_length"
)

# mutate all traits into numeric
sel_data_kp <- sel_data_kp %>% 
  mutate(across(all_of(obs_col), as.numeric))


# remove - , replace by _
names(sel_data_kp) = gsub("-", "_", names(sel_data_kp))

# Duplications in row and cols
duplicated_plot <- row_col_dup(sel_data_kp)

# Plot trial layout
trial_layout(sel_data_kp)


## Check the clone name
cloneName_new_old <- check_clone_name(
  clone_list = sel_data_kp$use_accession_name,
  new_names = NA,
  add_check = NULL
)

trial_standard <- sel_data_kp %>%
  left_join(cloneName_new_old,
            by = c("use_accession_name" = "accession_name_ori")
  ) %>%
  select(-use_accession_name) %>%
  rename(use_accession_name = use_accession_name.y)

## Add GIS data
trial_standard <- add_GIS(trial_standard)

# Get world map data for reference
vietnam_map <- map_data("world", region = "Vietnam")

# Create a data frame with your locations
locations <- data.frame(
  Location = c("Tay Ninh", "Dong Nai"),
  Latitude = c(11.3009, 11.1432),
  Longitude = c(106.1107, 107.2742)
)


library(ggrepel)

# Create a ggplot object for the map
static_map <- ggplot() +
  geom_polygon(data = vietnam_map, aes(x = long, y = lat, group = group), 
               fill = "lightgray", col = "black", linewidth = 0.03) +
  geom_point(data = locations, aes(x = Longitude, y = Latitude), size = 2, color = "red") +
  geom_text_repel(data = locations, aes(x = Longitude, y = Latitude, label = Location), 
                  size = 3) +
  geom_text(aes(x = 108.5, y = 10, label = "Author: Luis Fdo. Delgado"), size = 1.5) +
  labs(x =  "Latitude", y = "Longitude") +
  labs(title = "Locations in Vietnam", subtitle = "Tay Ninh and Dong Nai") +
  theme(plot.title = element_text(face = "bold.italic"),
        plot.subtitle = element_text(face = "italic")) +
  coord_fixed(1.3) # This sets the aspect ratio

# Display the map
print(static_map)

ggsave(paste("images\\map", trial_interest, ".png", sep = "_"),
       plot = static_map, units = "in", dpi = 300, width = 6, height = 5
)

# extract checks
accession_rep_ct <- trial_standard %>%
  count(use_trial_name, use_accession_name, use_rep_number)  %>%
  arrange(use_trial_name) %>%
  filter(n>1)
accession_rep_ct 

# Harvesting time
conducted_trials <- 
  trial_standard %>% group_by(use_trial_name, use_plant_date,use_harvest_date, use_location) %>% 
  summarise(n_gen = n_distinct(use_accession_name)) %>% 
  mutate(harvesting_time = 
           interval(ymd(use_plant_date), ymd(use_harvest_date)) %>% as.period,
         harvesting_time = paste0(harvesting_time@month, "month ", harvesting_time@day, "day")) %>% 
  ungroup()

conducted_trials

conducted_trials %>% relocate(harvesting_time, .after = use_harvest_date) %>% 
  write.table("clipboard", sep="\t", col.names = T, row.names = F)

## plot plant number

plants_plot <- trial_standard %>%
  group_by(use_trial_name) %>%
  count(obs_planted_number_plot) 
plants_plot


## Frequency harvest plant number
plants_harvested <- trial_standard %>%
  group_by(use_trial_name) %>%
  count(obs_harvest_number) %>% arrange(desc(obs_harvest_number))

# planted and harvested
plants_plot %>% select(-n) %>% 
  left_join(plants_harvested %>% 
              summarise(harvested_plants = max(obs_harvest_number)), by = "use_trial_name") %>% 
  write.table("clipboard", sep="\t", col.names = T, row.names = F)


# I dont think I am able to calculate yield_V2
plants_to_harvest <- plants_harvested %>% 
  ggplot(aes(x = factor(obs_harvest_number), 
             y = n, fill = factor(obs_harvest_number))) +
  geom_col(col = 'black') +
  #scale_fill_jco() +
  #theme_xiaofei() +
  theme(axis.text.x = element_text(vjust = 1, angle = 65),
        text = element_text(size = 20),
        legend.position="top")+
  labs(x = "Harvest_plant_number", y = "Freq", fill = "Harvest_plant_number") +
  facet_wrap(~ use_trial_name)

plants_to_harvest
ggsave(paste("images\\bar", trial_interest, Sys.Date(), ".png", sep = "_"),
       plot = plants_to_harvest, units = "in", dpi = 300, width = 9, height = 6)


## Compute germination, yield, yield_starch
trial_standard <- trial_standard %>%
  mutate(obs_harvest_number_plan =
           case_when(str_detect(use_trial_name, "202301")  ~ 8,
                     str_detect(use_trial_name, "202304")  ~ 8),
         obs_germination_perc = obs_germinated_number_plot/obs_planted_number_plot * 100,
         # 2) calculate area per plant
         area_plant = (use_plot_length*use_plot_width)/obs_planted_number_plot,
         # 3) calculate the yield_v4 obs_planted_number_plot
         obs_yield_ha = ifelse(obs_yield_ha == 0, NA, obs_yield_ha),
         obs_starch_content = ifelse(obs_starch_content == 0, NA, obs_starch_content),
         obs_yield_ha_v2 = (((obs_root_weight_plot*10000)/(area_plant*obs_harvest_number_plan))/1000),
         obs_starch_yield_ha = obs_starch_content * obs_yield_ha_v2 / 100)


# Plot the yield_v2 vs yield uploaded by Lizbeth
library(plotly)

p1 <- trial_standard_new %>% ggplot() +
  geom_point(aes(x = obs_yield_ha, y = obs_yield_ha_v2, color = use_plot_number), show.legend = F) +
  facet_wrap(~use_trial_name) +
  theme_xiaofei()

ggplotly(p1)
detach("package:plotly", unload = TRUE)

## Is numeric all traits?
is_numeric(trial_data = trial_standard)

# Get the tidy data
meta_info = names(trial_standard )[str_detect(names(trial_standard), "use_")]
meta_info = gsub("use_", "", meta_info)
meta_info
trial_tidy = trial_standard
names(trial_tidy)= gsub("use_", "", names(trial_standard))

# observations
trait_list = names(trial_tidy)[str_detect(names(trial_tidy), "obs_")]
trait_list = gsub("obs_", "", trait_list)
trait_list
names(trial_tidy)= gsub("obs_", "", names(trial_tidy))
trial_tidy = trial_tidy[c(meta_info, trait_list)]


# Boxplots and save in output folder
dev.off()
boxplot_traits(my_dat = trial_tidy,
               trait_wanted = trait_list,
               folder = paste0(here::here("output"), "/"), 
               trial_interest)


# Grouping boxplot
trait_wanted <- trait_list
plot_bxp <- trial_tidy %>%
  pivot_longer(
    cols = all_of(trait_wanted),
    names_to = "var",
    values_to = "values"
  ) %>%
  filter(!var %in% c(
    "stake_plant", "planted_number_plot",
    "harvest_number", "root_weight_air",
    "root_weight_water", "harvest_number_plan",
    "root_rot_perc", "yield_ha_v2"
  )) %>%
  ggplot(aes(x = trial_name, y = values)) +
  facet_wrap(~var,
             ncol = 7, scales = "free_y"
  ) + 
  geom_violin(fill = "gray") +
  geom_boxplot(width = 0.2) +
  labs(x = NULL, y = NULL, title = "") +
  theme_xiaofei() +
  theme(
    axis.text.x = element_text(size = 8, vjust = 1, angle = 65),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(color = "black"),
    strip.text.x = element_text(
      size = 8, face = "bold.italic")
  ) 

plot_bxp

ggsave(paste0("images\\boxplot_", trial_interest, Sys.Date(), ".png"),
       plot = plot_bxp, units = "in", dpi = 300, width = 16, height = 12
)

# Save the tidy data for analysis

write.csv(trial_tidy, here::here("output", paste("01_", year_interest, trial_interest,
                                                   "_tidy_data4analysis_", Sys.Date(), ".csv", sep = "")), row.names = FALSE)


