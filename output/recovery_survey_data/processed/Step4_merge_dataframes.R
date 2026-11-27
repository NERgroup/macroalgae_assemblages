#Predictive model of urchin gonad index
#Step 4: Merge DFs and harmonize predictors
#Sabrina N. Grant
#November 2026

rm(list = ls())

librarian::shelf(
  "tidyverse",
  "here",
  "janitor",
  "GGally"
)

################################################################################
#Set directories

procdir <- here::here("output")

mba_dir <- file.path(
  procdir, "recovery_survey_data", "processed"
)

dissection_dir <- file.path(
  procdir, "dissection_data", "processed"
)

mlpa_dir <- file.path(
  procdir, "mpla_data", "processed"
)

################################################################################
#Load MBA data

quad_data <- read_csv(
  file.path(mba_dir, "quad_data_recovery_clean.csv"),
  show_col_types = FALSE
)

kelp_data <- read_csv(
  file.path(mba_dir, "kelp_data_recovery_clean.csv"),
  show_col_types = FALSE
)

dissection_data <- read_csv(
  file.path(dissection_dir, "dissection_data_recovery_clean.csv"),
  show_col_types = FALSE
)

################################################################################
#Load MLPA data

mlpa_swath <- read_csv(
  file.path(mlpa_dir, "kelp_swath_counts_CB_BS.csv"),
  show_col_types = FALSE
)

mlpa_upc <- read_csv(
  file.path(mlpa_dir, "kelp_upc_cov_CB_BS.csv"),
  show_col_types = FALSE
)

################################################################################
#Define sampling hierarchy

patch_key <- c("year", "site", "site_type", "zone")

mlpa_key <- c("year", "site", "zone")

patch_levels <- c("BAR", "INCIP", "FOR")

################################################################################
#Inspect imported data

glimpse(quad_data)
glimpse(kelp_data)
glimpse(dissection_data)
glimpse(mlpa_swath)
glimpse(mlpa_upc)

#Confirm MBA relief categories
table(quad_data$relief_category, useNA = "ifany")

#Confirm MLPA habitat variables
summary(mlpa_upc$prop_complex_habitat)

#Compare crustose coralline cover units
summary(quad_data$cov_crustose_coralline)
summary(mlpa_upc$coralline_algae_crustose)

################################################################################
#Helper functions

mean_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

get_mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  ux <- unique(x)
  as.character(ux[which.max(tabulate(match(x, ux)))])
}

check_unique <- function(df, keys) {
  df %>%
    count(across(all_of(keys))) %>%
    filter(n > 1)
}

################################################################################
#MBA quadrat data

quad_build1 <- quad_data %>%
  mutate(
    survey_date = as.Date(survey_date),
    site_type = factor(site_type, levels = patch_levels),
    
    prop_bedrock = as.numeric(substrate == "Bedrock"),
    
    prop_flatrel = as.numeric(relief_category == "FLATREL"),
    prop_sltrel = as.numeric(relief_category == "SLTREL"),
    prop_modrel = as.numeric(relief_category == "MODREL"),
    prop_hirel = as.numeric(relief_category == "HIREL")
  )

#Aggregate quadrats to transects

quad_transect <- quad_build1 %>%
  group_by(
    across(all_of(c(patch_key, "survey_date", "transect")))
  ) %>%
  summarise(
    across(
      where(is.numeric) & !any_of("quadrat"),
      mean_na
    ),
    substrate = get_mode(substrate),
    .groups = "drop"
  )

#Aggregate transects to patch-depth surveys

quad_patch <- quad_transect %>%
  group_by(
    across(all_of(c(patch_key, "survey_date")))
  ) %>%
  summarise(
    across(
      where(is.numeric) & !any_of("transect"),
      mean_na
    ),
    substrate = get_mode(substrate),
    .groups = "drop"
  ) %>%
  mutate(
    prop_complex_habitat = prop_modrel + prop_hirel,
    
    prop_purp_concealed = if_else(
      purple_urchin_densitym2 > 0,
      purple_urchin_conceiledm2 / purple_urchin_densitym2,
      0
    ),
    
    prop_red_concealed = if_else(
      red_urchin_densitym2 > 0,
      red_urchin_conceiledm2 / red_urchin_densitym2,
      0
    ),
    
    substrate = factor(substrate)
  ) %>%
  select(
    -any_of(c(
      "purple_urchin_conceiledm2",
      "red_urchin_conceiledm2",
      "latitude",
      "longitude"
    ))
  )

glimpse(quad_patch)

################################################################################
#MBA kelp data

kelp_patch <- kelp_data %>%
  mutate(
    survey_date = as.Date(survey_date),
    site_type = factor(site_type, levels = patch_levels)
  ) %>%
  group_by(
    across(all_of(c(patch_key, "survey_date")))
  ) %>%
  summarise(
    across(
      where(is.numeric) & !any_of("transect"),
      mean_na
    ),
    .groups = "drop"
  ) %>%
  mutate(
    across(contains("20m2"), ~ .x / 20)
  ) %>%
  rename_with(
    ~ str_remove(.x, "20m2_?"),
    contains("20m2")
  )

glimpse(kelp_patch)

################################################################################
#Merge MBA community data

community_key <- c(patch_key, "survey_date")

check_unique(quad_patch, community_key)
check_unique(kelp_patch, community_key)

community_patch <- kelp_patch %>%
  left_join(
    quad_patch,
    by = community_key,
    relationship = "one-to-one"
  )

#Check unmatched surveys

kelp_patch %>%
  anti_join(quad_patch, by = community_key)

quad_patch %>%
  anti_join(kelp_patch, by = community_key)

glimpse(community_patch)

################################################################################
#MBA dissection data

dissection_zone <- dissection_data %>%
  filter(!is.na(zone)) %>%
  mutate(
    site_type = factor(site_type, levels = patch_levels)
  ) %>%
  group_by(
    across(all_of(c(patch_key, "species")))
  ) %>%
  summarise(
    n_urchins = n(),
    n_gonad_index = sum(!is.na(gonad_index)),
    
    #Calculate variability BEFORE calculating the mean
    gonad_index_sd = sd(gonad_index, na.rm = TRUE),
    gonad_index = mean_na(gonad_index),
    
    test_diameter_mm = mean_na(test_diameter_mm),
    animal_24hr_mass_g = mean_na(animal_24hr_mass_g),
    
    .groups = "drop"
  ) %>%
  mutate(
    gonad_index_se = gonad_index_sd / sqrt(n_gonad_index)
  )

glimpse(dissection_zone)

urchin_colors <- c(
  "purple_urchin" = "#8054A3",
  "red_urchin" = "#B94B45"
)

ggplot(
  dissection_zone,
  aes(
    x = site_type,
    y = gonad_index,
    color = species
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    alpha = 0.6,
    size = 2
  ) +
  scale_color_manual(
    values = urchin_colors,
    labels = c(
      "purple_urchin" = "Purple urchin",
      "red_urchin" = "Red urchin"
    )
  ) +
  facet_wrap(~ year) +
  theme_classic() +
  labs(
    x = "Patch type",
    y = "Mean gonad index (%)",
    color = NULL,
    title = "Urchin gonad index across patch types"
  ) +
  theme(legend.position = "top")

################################################################################
#Merge dissection and community data

gonad_model_data <- dissection_zone %>%
  left_join(
    community_patch,
    by = patch_key,
    relationship = "many-to-one"
  ) %>%
  mutate(
    site = factor(site),
    zone = factor(zone),
    species = factor(species),
    
    patch_uid = interaction(
      site, site_type, drop = TRUE
    ),
    
    zone_uid = interaction(
      patch_uid, zone, drop = TRUE
    ),
    
    survey_uid = interaction(
      zone_uid, year, drop = TRUE
    ),
    
    gonad_prop = gonad_index / 100
  )

################################################################################
#Quality control

gonad_model_data %>%
  summarise(
    n_observations = n(),
    n_sites = n_distinct(site),
    n_patches = n_distinct(patch_uid),
    n_zones = n_distinct(zone_uid),
    n_surveys = n_distinct(survey_uid),
    n_unmatched = sum(is.na(survey_date))
  )

################################################################################
#Export merged MBA data

output <- here::here(
  "output", "gonad_model", "processed"
)

dir.create(output, recursive = TRUE, showWarnings = FALSE)

write_csv(
  gonad_model_data,
  file.path(output, "mba_gonad_model_data.csv")
)