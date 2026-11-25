#Predictive model of urchin gonad index from patch-level community data
#Step 1: build the modeling data set (join quad, kelp swath and dissection data)
#Sabrina N. Grant
#November 24, 2026

librarian::shelf(tidyverse, tidymodels, here)

################################################################################
#set directories and load processed data

procdir <- here::here("output")

quad_data <- read_csv(file.path(procdir, "recovery_survey_data", "processed",
                                "quad_data_recovery_clean.csv"), show_col_types = FALSE)
kelp_data <- read_csv(file.path(procdir, "recovery_survey_data", "processed",
                                "kelp_data_recovery_clean.csv"), show_col_types = FALSE)
dissection_data <- read_csv(file.path(procdir, "dissection_data", "processed",
                                      "dissection_data_recovery_clean.csv"),
                            show_col_types = FALSE)

#join key: one row per patch (site x zone x patch type x year)
patch_key <- c("year", "site", "zone", "site_type")
patch_levels <- c("BAR", "INCIP", "FOR")

################################################################################
#summarize quad data to zone level (quadrats -> transects -> zone)

#helper function to get the most common value (mode)
get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

quad_patch <- quad_data %>%
  group_by(across(all_of(c(patch_key, "survey_date", "transect")))) %>%
  summarise(across(where(is.numeric) & !any_of("quadrat"), ~ mean(.x, na.rm = TRUE)),
            substrate = get_mode(substrate),
            .groups = "drop") %>%
  group_by(across(all_of(c(patch_key, "survey_date")))) %>%
  summarise(across(where(is.numeric) & !any_of("transect"), ~ mean(.x, na.rm = TRUE)),
            substrate = get_mode(substrate),
            .groups = "drop") %>%
  #proportion of urchins concealed (from patch-level mean densities)
  mutate(prop_purp_concealed = if_else(purple_urchin_densitym2 > 0,
                                       purple_urchin_conceiledm2 / purple_urchin_densitym2, 0),
         prop_red_concealed = if_else(red_urchin_densitym2 > 0,
                                      red_urchin_conceiledm2 / red_urchin_densitym2, 0)) %>%
  dplyr::select(-purple_urchin_conceiledm2, -red_urchin_conceiledm2,
                -latitude, -longitude)

################################################################################
#summarize kelp swath data to patch level and convert to density / m2

kelp_patch <- kelp_data %>%
  group_by(across(all_of(c(patch_key, "survey_date")))) %>%
  summarise(across(where(is.numeric) & !any_of("transect"), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  #20 m2 swath counts -> per m2
  mutate(across(contains("20m2"), ~ .x / 20)) %>%
  rename_with(~ str_remove(.x, "20m2_?"), contains("20m2"))

################################################################################
#summarize dissection data to the zone level (one row per patch x species)
#and join with the community data

community_patch <- kelp_patch %>%
  left_join(quad_patch, by = c(patch_key, "survey_date"))

dissection_zone <- dissection_data %>%
  #one 2025 REC_05 BAR sample has no zone and cannot be matched to a patch
  filter(!is.na(zone)) %>%
  mutate(site_type = factor(site_type, levels = patch_levels)) %>%
  group_by(across(all_of(c(patch_key, "species")))) %>%
  #sd is computed before gonad_index is overwritten by its mean
  summarise(n_urchins = n(),
            gonad_index_sd = sd(gonad_index, na.rm = TRUE),
            gonad_index = mean(gonad_index, na.rm = TRUE),
            test_diameter_mm = mean(test_diameter_mm, na.rm = TRUE),
            animal_24hr_mass_g = mean(animal_24hr_mass_g, na.rm = TRUE),
            .groups = "drop")

gonad_model_data <- dissection_zone %>%
  left_join(community_patch, by = patch_key)

################################################################################
#checks

#patches in the community data: one row per key?
community_patch %>% count(across(all_of(patch_key))) %>% filter(n > 1)

#dissected urchins that did not match a patch
gonad_model_data %>% filter(is.na(survey_date)) %>% count(across(all_of(patch_key)))

#row count preserved by the join?
nrow(gonad_model_data) == nrow(dissection_zone)

#one row per patch x species?
gonad_model_data %>% count(across(all_of(c(patch_key, "species")))) %>% filter(n > 1)

#missing values in the joined data
colSums(is.na(gonad_model_data))[colSums(is.na(gonad_model_data)) > 0]

#rows and urchins per species and year; small n_urchins means noisy means
gonad_model_data %>% group_by(species, year) %>%
  summarise(zones = n(), urchins = sum(n_urchins), min_n = min(n_urchins),
            median_n = median(n_urchins), .groups = "drop")

################################################################################
#next: define outcome and predictors, then split data (by site)
