#Predictive model of urchin gonad index
#Step 5: Harmonize MBA and MLPA predictors
#Sabrina N. Grant
#November 2026

rm(list = ls())

librarian::shelf(
  "tidyverse",
  "here",
  "GGally"
)

################################################################################
#Set directories

procdir <- here::here("output")

gonad_dir <- file.path(
  procdir, "gonad_model", "processed"
)

mlpa_dir <- file.path(
  procdir, "mpla_data", "processed"
)

output <- gonad_dir

################################################################################
#Import processed data

gonad_model_data <- read_csv(
  file.path(gonad_dir, "mba_gonad_model_data.csv"),
  show_col_types = FALSE
)

mlpa_swath <- read_csv(
  file.path(mlpa_dir, "kelp_swath_counts_CB_BS.csv"),
  show_col_types = FALSE
)

mlpa_upc <- read_csv(
  file.path(mlpa_dir, "kelp_upc_cov_CB_BS.csv"),
  show_col_types = FALSE
)

################################################################################
#Restore MBA variable classes

gonad_model_data <- gonad_model_data %>%
  mutate(
    site_type = factor(
      site_type,
      levels = c("BAR", "INCIP", "FOR")
    ),
    across(
      c(site, zone, species, patch_uid, zone_uid, survey_uid),
      as.factor
    ),
    survey_date = as.Date(survey_date)
  )

################################################################################
#Helper functions

mean_na <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

################################################################################
#Create MBA predictor groups

gonad_model_predictor_groups <- gonad_model_data %>%
  mutate(
    density_understory_algae =
      densityptecal +
      densityeisarb +
      densitylamset,
    
    density_kelp_recruits =
      lamr + macr,
    
    density_kelp_juveniles =
      macj + nerj + ptej + lsetj + eisj,
    
    density_kelp_stumps =
      densitymacstump + densitylamstump
  )

################################################################################
#Define shared environmental predictors

shared_predictors <- c(
  "macro_stipe_density_",
  "densitynerlue",
  "density_understory_algae",
  "purple_urchin_densitym2",
  "red_urchin_densitym2",
  "cov_crustose_coralline",
  "prop_complex_habitat"
)

################################################################################
#Check shared predictors

gonad_model_predictor_groups %>%
  select(all_of(shared_predictors)) %>%
  summary()

gonad_model_predictor_groups %>%
  summarise(
    across(
      all_of(shared_predictors),
      ~ sum(!is.finite(.x))
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "predictor",
    values_to = "n_missing"
  )

################################################################################
#Prepare MBA training data

mba_training_data <- gonad_model_predictor_groups %>%
  select(
    year, site, site_type, zone, survey_date,
    patch_uid, zone_uid, survey_uid,
    species, n_urchins, n_gonad_index,
    gonad_index, gonad_index_sd, gonad_index_se,
    gonad_prop,
    all_of(shared_predictors)
  )

#Unique environmental observations
mba_environment <- mba_training_data %>%
  distinct(survey_uid, .keep_all = TRUE) %>%
  select(
    year, site, site_type, zone, survey_uid,
    all_of(shared_predictors)
  )

glimpse(mba_training_data)
glimpse(mba_environment)

################################################################################
#MLPA sampling identifiers

mlpa_key <- c("year", "site", "zone")

mlpa_transect_key <- c(mlpa_key, "transect")

#Check for duplicate transects

mlpa_swath %>%
  count(across(all_of(mlpa_transect_key))) %>%
  filter(n > 1)

mlpa_upc %>%
  count(across(all_of(mlpa_transect_key))) %>%
  filter(n > 1)

################################################################################
#MLPA UPC and habitat predictors

mlpa_upc_zone <- mlpa_upc %>%
  rename(
    cov_crustose_coralline = coralline_algae_crustose
  ) %>%
  group_by(across(all_of(mlpa_key))) %>%
  summarise(
    n_upc_transects = n_distinct(transect),
    
    across(
      c(
        cov_crustose_coralline,
        prop_flatrel,
        prop_sltrel,
        prop_modrel,
        prop_hirel,
        prop_complex_habitat,
        prop_bedrock,
        prop_boulder,
        prop_cobble,
        prop_sand
      ),
      mean_na
    ),
    
    .groups = "drop"
  )

glimpse(mlpa_upc_zone)

################################################################################
#Quality control

mlpa_upc_zone %>%
  summarise(
    n_sites = n_distinct(site),
    n_years = n_distinct(year),
    first_year = min(year),
    last_year = max(year),
    n_surveys = n(),
    n_missing_cca = sum(is.na(cov_crustose_coralline)),
    n_missing_relief = sum(is.na(prop_complex_habitat))
  )

################################################################################
#Inspect MLPA SWATH predictors

mlpa_swath_predictors <- c(
  "macrocystis_stipes",
  "nereocystis_luetkeana",
  "pterygophora_californica",
  "eisenia_arborea",
  "laminaria_setchellii",
  "strongylocentrotus_purpuratus",
  "mesocentrotus_franciscanus"
)

#Check that all variables exist
setdiff(mlpa_swath_predictors, names(mlpa_swath))

mlpa_swath %>%
  select(
    all_of(mlpa_transect_key),
    all_of(mlpa_swath_predictors)
  ) %>%
  glimpse()

#Summarize counts
mlpa_swath %>%
  select(all_of(mlpa_swath_predictors)) %>%
  summary()

#Check for sampling effort variables
grep(
  "area|width|length|size|effort",
  names(mlpa_swath),
  value = TRUE,
  ignore.case = TRUE
)

################################################################################
#Convert MLPA SWATH counts to densities per m2
#Standard MLPA SWATH transect: 30 m x 2 m = 60 m2

swath_area_m2 <- 60

mlpa_swath_build1 <- mlpa_swath %>%
  mutate(
    macro_stipe_density_ =
      macrocystis_stipes / swath_area_m2,
    
    densitynerlue =
      nereocystis_luetkeana / swath_area_m2,
    
    density_understory_algae =
      (pterygophora_californica +
         eisenia_arborea +
         laminaria_setchellii) / swath_area_m2,
    
    purple_urchin_densitym2 =
      strongylocentrotus_purpuratus / swath_area_m2,
    
    red_urchin_densitym2 =
      mesocentrotus_franciscanus / swath_area_m2
  )

swath_predictors <- c(
  "macro_stipe_density_",
  "densitynerlue",
  "density_understory_algae",
  "purple_urchin_densitym2",
  "red_urchin_densitym2"
)

mlpa_swath_build1 %>%
  select(all_of(swath_predictors)) %>%
  summary()

################################################################################
#Aggregate MLPA SWATH densities to site-zone-year

mlpa_swath_zone <- mlpa_swath_build1 %>%
  group_by(across(all_of(mlpa_key))) %>%
  summarise(
    latitude = mean_na(latitude),
    longitude = mean_na(longitude),
    affiliated_mpa = first(affiliated_mpa),
    mpa_designation = first(mpa_designation),
    
    n_swath_transects = n_distinct(transect),
    
    across(
      all_of(swath_predictors),
      mean_na
    ),
    
    .groups = "drop"
  )

glimpse(mlpa_swath_zone)

################################################################################
#Quality control

mlpa_swath_zone %>%
  summarise(
    n_sites = n_distinct(site),
    n_years = n_distinct(year),
    first_year = min(year),
    last_year = max(year),
    n_surveys = n()
  )

#Transect replication

mlpa_swath_zone %>%
  count(n_swath_transects)


################################################################################
#Merge MLPA environmental datasets

mlpa_prediction_data <- mlpa_swath_zone %>%
  left_join(
    mlpa_upc_zone,
    by = mlpa_key,
    relationship = "one-to-one"
  ) %>%
  mutate(
    zone_uid = interaction(
      site, zone,
      drop = TRUE
    ),
    
    survey_uid = interaction(
      zone_uid, year,
      drop = TRUE
    )
  ) %>%
  select(
    year, site, zone,
    zone_uid, survey_uid,
    latitude, longitude,
    affiliated_mpa, mpa_designation,
    n_swath_transects, n_upc_transects,
    all_of(shared_predictors),
    everything()
  )

glimpse(mlpa_prediction_data)

################################################################################
#Combine environmental datasets

environment_comparison <- bind_rows(
  
  mba_environment %>%
    select(year, site, zone, survey_uid,
           all_of(shared_predictors)) %>%
    mutate(
      dataset = "MBA",
      across(c(site, zone, survey_uid), as.character)
    ),
  
  mlpa_prediction_data %>%
    select(year, site, zone, survey_uid,
           all_of(shared_predictors)) %>%
    mutate(
      dataset = "MLPA",
      across(c(site, zone, survey_uid), as.character)
    )
  
) %>%
  mutate(
    dataset = factor(dataset, levels = c("MBA", "MLPA"))
  )

#Define predictor labels

predictor_labels <- c(
  macro_stipe_density_ = "Giant kelp stipes",
  densitynerlue = "Bull kelp density",
  density_understory_algae = "Understory kelp density",
  purple_urchin_densitym2 = "Purple urchin density",
  red_urchin_densitym2 = "Red urchin density",
  cov_crustose_coralline = "Crustose coralline cover",
  prop_complex_habitat = "Proportion complex habitat"
)

dataset_colors <- c(
  "MBA" = "#6E8B3D",
  "MLPA" = "#4682B4"
)

#Convert to long format

environment_long <- environment_comparison %>%
  pivot_longer(
    cols = all_of(shared_predictors),
    names_to = "predictor",
    values_to = "value"
  ) %>%
  mutate(
    predictor = factor(
      predictor,
      levels = shared_predictors
    )
  )

glimpse(environment_long)

################################################################################
#Compare environmental predictor distributions

environment_summary <- environment_long %>%
  filter(is.finite(value)) %>%
  group_by(predictor, dataset) %>%
  summarise(
    n = n(),
    pct_zero = 100 * mean(value == 0),
    min = min(value),
    q05 = quantile(value, 0.05),
    q25 = quantile(value, 0.25),
    median = median(value),
    mean = mean(value),
    q75 = quantile(value, 0.75),
    q95 = quantile(value, 0.95),
    max = max(value),
    sd = sd(value),
    .groups = "drop"
  )

print(environment_summary, n = Inf, width = Inf)

################################################################################
#Density distributions

p_environment_density <- environment_long %>%
  filter(is.finite(value)) %>%
  ggplot(
    aes(
      x = value,
      color = dataset,
      fill = dataset
    )
  ) +
  geom_density(
    alpha = 0.2,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ predictor,
    scales = "free",
    ncol = 3,
    labeller = as_labeller(predictor_labels)
  ) +
  scale_color_manual(values = dataset_colors) +
  scale_fill_manual(values = dataset_colors) +
  theme_classic(base_size = 12) +
  labs(
    x = "Predictor value",
    y = "Probability density",
    color = NULL,
    fill = NULL,
    title = "Environmental predictor distributions",
    subtitle = "MBA training data vs. historical MLPA data"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_environment_density


################################################################################
#Empirical cumulative distributions

p_environment_ecdf <- environment_long %>%
  filter(is.finite(value)) %>%
  ggplot(
    aes(
      x = value,
      color = dataset
    )
  ) +
  stat_ecdf(
    geom = "step",
    linewidth = 0.9
  ) +
  facet_wrap(
    ~ predictor,
    scales = "free_x",
    ncol = 3,
    labeller = as_labeller(predictor_labels)
  ) +
  scale_color_manual(values = dataset_colors) +
  scale_y_continuous(
    labels = scales::label_percent()
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Predictor value",
    y = "Cumulative proportion of observations",
    color = NULL,
    title = "Cumulative environmental distributions"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_environment_ecdf


################################################################################
#Compare Spearman predictor correlations

mba_cor <- mba_environment %>%
  select(all_of(shared_predictors)) %>%
  cor(
    method = "spearman",
    use = "pairwise.complete.obs"
  )

mlpa_cor <- mlpa_prediction_data %>%
  select(all_of(shared_predictors)) %>%
  cor(
    method = "spearman",
    use = "pairwise.complete.obs"
  )

round(mba_cor, 2)
round(mlpa_cor, 2)

################################################################################
#Convert correlation matrices to long format

cor_to_long <- function(cor_matrix, dataset_name) {
  
  as.data.frame(as.table(cor_matrix)) %>%
    as_tibble() %>%
    rename(
      predictor1 = Var1,
      predictor2 = Var2,
      rho = Freq
    ) %>%
    mutate(dataset = dataset_name)
}

cor_comparison <- bind_rows(
  cor_to_long(mba_cor, "MBA"),
  cor_to_long(mlpa_cor, "MLPA")
) %>%
  mutate(
    dataset = factor(dataset, levels = c("MBA", "MLPA")),
    predictor1 = factor(predictor1, levels = shared_predictors),
    predictor2 = factor(predictor2, levels = rev(shared_predictors))
  )

################################################################################
#Correlation heatmaps

p_cor_comparison <- ggplot(
  cor_comparison,
  aes(
    x = predictor1,
    y = predictor2,
    fill = rho
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = round(rho, 2)),
    size = 3
  ) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_fill_gradient2(
    low = "#4682B4",
    mid = "white",
    high = "#B94B45",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_x_discrete(labels = predictor_labels) +
  scale_y_discrete(labels = predictor_labels) +
  coord_fixed() +
  theme_classic(base_size = 11) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Spearman rho",
    title = "Environmental predictor correlations"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_cor_comparison

################################################################################
#MBA training ranges

mba_ranges <- environment_long %>%
  filter(
    dataset == "MBA",
    is.finite(value)
  ) %>%
  group_by(predictor) %>%
  summarise(
    mba_min = min(value),
    mba_max = max(value),
    mba_q05 = quantile(value, 0.05),
    mba_q95 = quantile(value, 0.95),
    .groups = "drop"
  )

################################################################################
#Quantify extrapolation

extrapolation_long <- environment_long %>%
  filter(dataset == "MLPA") %>%
  left_join(mba_ranges, by = "predictor") %>%
  mutate(
    extrapolation = case_when(
      !is.finite(value) ~ NA_character_,
      value < mba_min ~ "Below",
      value > mba_max ~ "Above",
      TRUE ~ "Within"
    )
  )

extrapolation_summary <- extrapolation_long %>%
  filter(!is.na(extrapolation)) %>%
  group_by(predictor) %>%
  summarise(
    n_mlpa = n(),
    mba_min = first(mba_min),
    mba_max = first(mba_max),
    mlpa_min = min(value),
    mlpa_max = max(value),
    n_below = sum(extrapolation == "Below"),
    n_above = sum(extrapolation == "Above"),
    pct_below = 100 * mean(extrapolation == "Below"),
    pct_above = 100 * mean(extrapolation == "Above"),
    pct_outside = 100 * mean(extrapolation != "Within"),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_outside))

print(extrapolation_summary, width = Inf)

################################################################################
#Flag historical observations requiring extrapolation

extrapolation_flags <- extrapolation_long %>%
  group_by(year, site, zone, survey_uid) %>%
  summarise(
    n_predictors_available = sum(!is.na(extrapolation)),
    
    n_predictors_outside = sum(
      extrapolation %in% c("Below", "Above")
    ),
    
    extrapolated_predictors = paste(
      as.character(
        predictor[extrapolation %in% c("Below", "Above")]
      ),
      collapse = ", "
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    extrapolation_status = case_when(
      n_predictors_available < length(shared_predictors) ~
        "Incomplete predictors",
      n_predictors_outside > 0 ~
        "Outside training range",
      TRUE ~ "Within training range"
    )
  )

extrapolation_flags %>%
  count(extrapolation_status) %>%
  mutate(percent = 100 * n / sum(n))