#Urchin dissection (gonad index) data processing
#Sabrina N. Grant
#November 24, 2026

librarian::shelf(tidyverse, janitor, here)

################################################################################
#set directories and load data

#data are on the enhydra server (must be mounted)
datin <- "/Volumes/enhydra/data/kelp_recovery/MBA_kelp_forest_database/processed/dissection"
output <- here::here("output", "dissection_data", "processed")

dissection_raw <- read_csv(file.path(datin, "dissection_data_recovery.csv"),
                           show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  #harmonize names with the MBA quadrat data
  rename(site = site_official,
         site_type = site_type_official,
         date_collected = survey_date)

#stop early if the file does not have the columns this script expects
expected_cols <- c("survey_type", "region", "site", "site_type", "zone",
                   "date_collected", "treatment", "species", "sex", "test_height_mm",
                   "test_diameter_mm", "animal_24hr_mass_g", "gonad_mass_g",
                   "sample_number", "gonad_index")
missing_cols <- setdiff(expected_cols, names(dissection_raw))
if (length(missing_cols) > 0) {
  stop("Missing expected columns: ", paste(missing_cols, collapse = ", "),
       "\nColumns found: ", paste(names(dissection_raw), collapse = ", "))
}

################################################################################
#inspect structure before filtering

dim(dissection_raw)
table(dissection_raw$survey_type, useNA = "ifany")
table(dissection_raw$species, useNA = "ifany")
table(dissection_raw$treatment, useNA = "ifany")
table(dissection_raw$site_type, useNA = "ifany")

################################################################################
#process

dissection_build1 <- dissection_raw %>%
  mutate(date_collected = as.Date(date_collected),
         year = as.numeric(format(date_collected, "%Y"))) %>%
  #keep recovery patch surveys, matching the MBA quadrat data (2024-2025)
  filter(survey_type == "Recovery", year %in% c(2024, 2025)) %>%
  mutate(site_type = factor(site_type, levels = c("BAR", "INCIP", "FOR")),
         #gonad index as a percentage of whole-animal mass
         gonad_index_check = gonad_mass_g / animal_24hr_mass_g * 100) %>%
  dplyr::select(year, date_collected, site, site_type, zone, treatment,
                species, sex, test_height_mm, test_diameter_mm, animal_24hr_mass_g,
                gonad_mass_g, gonad_index, gonad_index_check, sample_number)

################################################################################
#checks (review output before exporting)

#sample size by patch type, species and year
dissection_build1 %>% count(year, site_type, species)

#samples per site, zone and patch type
dissection_build1 %>% count(year, site, zone, site_type, species) %>% summary()

#missing values
colSums(is.na(dissection_build1))

#does the stored gonad_index match gonad mass / body mass? (> 0.5 pct. point difference)
dissection_build1 %>%
  filter(abs(gonad_index - gonad_index_check) > 0.5) %>%
  dplyr::select(year, site, zone, species, gonad_mass_g, animal_24hr_mass_g,
                gonad_index, gonad_index_check)

#range of gonad index and implausible values
summary(dissection_build1$gonad_index)
dissection_build1 %>% filter(gonad_index < 0 | gonad_index > 40) %>%
  dplyr::select(year, site, zone, species, gonad_index)

#duplicate samples
dissection_build1 %>%
  count(year, site, zone, species, sample_number) %>%
  filter(n > 1)

#site-zone-type patches that match the MBA quadrat data
#(compare against kelp_quad_data_CC.csv once loaded)

################################################################################
#export

dissection_final <- dissection_build1 %>% dplyr::select(-gonad_index_check)

dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(dissection_final, file.path(output, "dissection_data_recovery_clean.csv"),
          row.names = FALSE)
