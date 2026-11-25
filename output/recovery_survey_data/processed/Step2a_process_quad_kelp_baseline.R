#MBA quadrat and kelp swath baseline data processing
#(formatted to match the dissection data: clean names, year column,
# site_type factor, recovery surveys 2024-2025, checks, then export)
#Sabrina N. Grant
#November 24, 2026

librarian::shelf(tidyverse, janitor, here)

################################################################################
#set directories and load data

#data are on the enhydra server (must be mounted)
datin <- "/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/processed/recovery/kelp_recovery_data.rda"
output <- here::here("output", "recovery_survey_data", "processed")

#loads quad_data and kelp_data
load(datin)

#stop early if the data do not have the columns this script expects
expected_quad <- c("survey_type", "region", "site", "site_type", "zone",
                   "survey_date", "transect", "quadrat", "purple_urchin_densitym2")
expected_kelp <- c("site", "site_type", "zone", "survey_date", "transect",
                   "n_macro_plants_20m2")
for (x in list(list(quad_data, expected_quad), list(kelp_data, expected_kelp))) {
  miss <- setdiff(x[[2]], names(janitor::clean_names(x[[1]])))
  if (length(miss) > 0) stop("Missing expected columns: ", paste(miss, collapse = ", "))
}

################################################################################
#inspect structure before filtering

dim(quad_data); dim(kelp_data)
table(quad_data$survey_type, useNA = "ifany")
table(quad_data$site_type, useNA = "ifany")

################################################################################
#process quad data (quadrat-level)

quad_build1 <- quad_data %>%
  janitor::clean_names() %>%
  mutate(survey_date = as.Date(survey_date),
         year = as.numeric(format(survey_date, "%Y"))) %>%
  #48 rows (2025: REC_03 FOR, REC_05 FOR/INCIP shallow) lack survey_type, region
  #and coordinates; fill from other rows of the same site and patch type
  group_by(site, site_type) %>%
  fill(survey_type, region, latitude, longitude, .direction = "downup") %>%
  ungroup() %>%
  filter(survey_type == "Recovery", year %in% c(2024, 2025)) %>%
  mutate(site_type = factor(site_type, levels = c("BAR", "INCIP", "FOR"))) %>%
  dplyr::select(year, survey_date, region, site, site_type, zone, transect,
                quadrat, everything(), -survey_type)

################################################################################
#process kelp data (transect-level)

kelp_build1 <- kelp_data %>%
  janitor::clean_names() %>%
  mutate(survey_date = as.Date(survey_date),
         year = as.numeric(format(survey_date, "%Y"))) %>%
  filter(year %in% c(2024, 2025)) %>%
  mutate(site_type = factor(site_type, levels = c("BAR", "INCIP", "FOR"))) %>%
  dplyr::select(year, survey_date, site, site_type, zone, transect, everything())

################################################################################
#checks (review output before exporting)

#sample size by year and patch type
quad_build1 %>% count(year, site_type)
kelp_build1 %>% count(year, site_type)

#missing values
colSums(is.na(quad_build1))[colSums(is.na(quad_build1)) > 0]
colSums(is.na(kelp_build1))[colSums(is.na(kelp_build1)) > 0]

#duplicate samples
quad_build1 %>% count(year, site, site_type, zone, survey_date, transect, quadrat) %>% filter(n > 1)
kelp_build1 %>% count(year, site, site_type, zone, survey_date, transect) %>% filter(n > 1)

#percent cover should sum to 100 per quadrat
quad_build1 %>%
  mutate(cover_total = rowSums(across(starts_with("cov_")))) %>%
  filter(abs(cover_total - 100) > 0.5) %>%
  dplyr::select(year, site, zone, transect, quadrat, cover_total)

#implausible values
summary(quad_build1[c("relief_cm", "risk_index", "purple_urchin_densitym2")])

#site-zone-type patches shared with the dissection data
quad_build1 %>% distinct(year, site, zone, site_type) %>%
  anti_join(read_csv(file.path(output, "..", "..", "dissection_data", "processed",
                               "dissection_data_recovery_clean.csv"),
                     show_col_types = FALSE) %>%
              distinct(year, site, zone, site_type = factor(site_type, levels = levels(quad_build1$site_type))),
            by = c("year", "site", "zone", "site_type"))

################################################################################
#export

dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(quad_build1, file.path(output, "quad_data_recovery_clean.csv"), row.names = FALSE)
write.csv(kelp_build1, file.path(output, "kelp_data_recovery_clean.csv"), row.names = FALSE)
