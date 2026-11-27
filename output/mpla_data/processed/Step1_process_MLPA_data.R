#MLPA Community data processing
#Sabrina N. Grant (adapted from Joshua G. Smith's script)
#April 10, 2026

rm(list=ls())
librarian::shelf(tidyverse, here)


###subtidal monitoring data were accessed at
###https://opc.dataone.org/view/doi:10.25494/P6/MLPA_kelpforest.12

librarian::shelf(data.table, dplyr, ggplot2, tidyr, here)

#MLPA_kelpforest_swath.10.csv
#MLPA_kelpforest_upc.10.csv
#MLPA_kelpforest_site_table.10.csv
#MLPA_kelpforest_taxon_table.10.csv


################################################################################
#set directories and load data
basedir <- here::here("output","mpla_data","raw")           
output <- here::here("output","mpla_data","processed")

kelp_upc_raw <- read.csv(file.path(basedir, "MLPA_kelpforest_upc.10.csv")) %>%
  janitor::clean_names()

kelp_swath_raw <- read.csv(file.path(basedir, "MLPA_kelpforest_swath.10.csv")) %>%
  janitor::clean_names()

kelp_taxon <- read.csv(file.path(basedir, "MLPA_kelpforest_taxon_table.10.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(classcode, species_definition) %>%
  distinct()

site_table <- read.csv(file.path(basedir, "MLPA_kelpforest_site_table.10.csv")) %>%
  janitor::clean_names() %>%
  dplyr::select(site, latitude, longitude, ca_mpa_name_short, mpa_class=site_designation, mpa_designation=site_status, 
                baseline_region)%>%
  distinct() #remove duplicates


################################################################################
#process kelp swath

#select vars and clean
kelp_swath_build1 <- kelp_swath_raw %>%
  #drop juvenile sea urchins, since these were not consistently recorded 
  dplyr::filter(!(classcode == "STRPURREC" |
                    classcode == "MESFRAREC"))%>%
  dplyr::select(year, site, zone, transect, classcode, count, size)%>%
  group_by(year, site, zone, transect, classcode)%>%
  dplyr::summarize(total_count = sum(count), #counts in raw data are grouped by size class. Take summary across all sizes
                   total_stipes = sum(size * count, na.rm = TRUE)) %>% #stipe counts only: size = stipes per plant, so weight by plants counted
  #CHANGED: now creates two columns — total_count (individuals for all species) and
  #total_count_stipes (stipes for Macrocystis, individuals for all other species)
  mutate(total_count_stipes = ifelse(classcode == "MACPYRAD", total_stipes, total_count)) %>%
  dplyr::select(!(total_stipes))

#join species names by class code 
kelp_swath_build2 <- left_join(kelp_swath_build1, kelp_taxon, by="classcode")

#CHANGED: added total_count_stipes to carry both count columns forward
kelp_swath_build3 <- kelp_swath_build2 %>%
  dplyr::select(year, site, zone, transect, species=species_definition, 
                total_count, total_count_stipes)


#add affiliated_mpa

kelp_swath_build4 <- left_join(kelp_swath_build3, site_table, by="site")

kelp_swath_build5 <- kelp_swath_build4 %>%
  ungroup() %>%
  #CHANGED: added total_count_stipes to select
  dplyr::select(year, baseline_region, site,latitude, longitude, affiliated_mpa=ca_mpa_name_short,mpa_class, mpa_designation, zone, transect, species,
                total_count, total_count_stipes) %>%
  mutate(mpa_designation = ifelse(mpa_designation=="reference","ref",mpa_class)) %>%
  #make sure species are not duplicates by summarizing at the transect level (total counts)
  group_by(year, baseline_region, latitude, longitude, site, affiliated_mpa, mpa_class, mpa_designation, zone,  transect, species) %>%
  #CHANGED: summarises both counts (individuals) and counts_stipes (stipes for Macrocystis)
  dplyr::summarise(counts = sum(total_count),
                   counts_stipes = sum(total_count_stipes),
                   .groups = "drop")

#NEW: pull macrocystis stipe counts before pivot
macro_stipes <- kelp_swath_build5 %>%
  filter(species == "Macrocystis pyrifera") %>%
  dplyr::select(year, baseline_region, latitude, longitude, site, affiliated_mpa, 
                mpa_class, mpa_designation, zone, transect, 
                macrocystis_stipes = counts_stipes)

#reshape to wide format — pivot individuals only, as original
kelp_swath_build6 <- kelp_swath_build5 %>%
  dplyr::select(-counts_stipes) %>%  #drop counts_stipes before pivot
  pivot_wider(names_from = species, values_from = counts, values_fill = 0) %>%
  janitor::clean_names() %>%
  mutate(MHW = ifelse(year >= 2014 & year <= 2016, "during", ifelse(year < 2014, "before", "after"))) %>%
  dplyr::select(-c(no_organisms_present_in_this_sample)) %>%
  filter(baseline_region == "CENTRAL") %>%
  filter(year >= 2001) %>% #NEW: starting from 2001 
  #join macrocystis stipes as single new column
  left_join(macro_stipes, by = c("year", "baseline_region", "latitude", "longitude",
                                 "site", "affiliated_mpa", "mpa_class", 
                                 "mpa_designation", "zone", "transect")) %>%
  #transects with no Macrocystis have no row in macro_stipes; true count is 0, not NA
  mutate(macrocystis_stipes = replace_na(macrocystis_stipes, 0))

unique(kelp_swath_build6$site)


#drop species and combine to higher taxa
#### see reference taxonomy table https://docs.google.com/spreadsheets/d/1vxy0XVOrlNhXD-i9tWL_F9G5h8S3S4OV/edit#gid=2031917236

kelp_swath_build7 <- as.data.frame(kelp_swath_build6) %>%
  #merge species
  dplyr::mutate(urticina_merge = rowSums(dplyr::select(.,'urticina', 
                                                       'urticina_coriacea', 'urticina_crassicornis', 'urticina_piscivora')
  ))%>%
  dplyr::select(!(c('urticina', 
                    'urticina_coriacea', 'urticina_crassicornis', 'urticina_piscivora'))) %>%
  #drop species
  dplyr::select(!(c('anthopleura',
                    'apostichopus',
                    'asteroidea',
                    'haliotis',
                    'laminaria',
                    'laminariales',
                    'octopus',
                    'pisaster',
                    'pugettia_spp',
                    'stylasterias_forreri',
                    'urticina_columbiana',
                    'urticina_columbiana_mcpeaki',
                    "unidentified_mobile_invert_species",
  )))%>%
  #drop string
  rename(urticina=urticina_merge)


#drop species that were never observed on the central coast
kelp_swath_zero_drop <- kelp_swath_build7 %>% filter(baseline_region=='CENTRAL') %>%
  dplyr::select(!(where(~ any(. != 0))))

#NOTE: macrocystis_stipes may be dropped here if it contains zeros — 
#consider pulling it out before this step and rejoining after if needed
kelp_swath_build8 <- kelp_swath_build7 %>% filter(baseline_region=='CENTRAL') %>%
  dplyr::select(where(~ any(. != 0))) %>%
  dplyr::select(year, MHW, everything())


names(kelp_swath_build8)

################################################################################
#process kelp upc

kelp_upc_build1 <- kelp_upc_raw%>%
  filter(category=="COVER") %>%
  #in the taxonomic guide, "ZONSPP" and "DICTYOTALES" are both classcodes for "dictyotales" so I am removing "ZONSPP" for now
  filter(classcode != "ZONSPP") 

#select vars and clean
kelp_upc_build2 <- kelp_upc_build1 %>%
  dplyr::select(year, site, zone, transect, classcode, count, pct_cov)%>%
  group_by(year, site, zone, transect, classcode)

#join species names by class code 
kelp_upc_build3 <- left_join(kelp_upc_build2, kelp_taxon, by="classcode")


kelp_upc_build4 <- kelp_upc_build3 %>%
  dplyr::select(year, site, zone, transect, classcode, species=species_definition, count, pct_cov)

#drop data 
kelp_upc_build5 <- kelp_upc_build4 %>%
  filter(!(species=="Haliotis"|
             species=="Laminaria"|
             species=="laminariales"|
             species=="UCSC 2000 Definition Of 'Other Red' Included All Red Algae Except Rhospp, Bushy, Leafy, Crucor And Erecor"|
             species=="UCSB 2000 Definition Of 'Other Red' Included All Red Algae Except Rhospp, Bushy, Leafy, Crucor And Erecor"|
             species=="Red Alga With Cylindrical Branches. Definition Used In 2000"|
             species=="UCSC and UCSB 1999 Definition Of 'Other Red' Included All Red Algae Except Rhospp, Crucor And Erecor"|
             species=="Laminariales Holdfast (Alive)"|
             species=="Sargassum"|
             species=="Aglaophenia struthionides"|
             species=="Macrosystis Holdfast (Dead)"|
             species=="Bare Rock"|
             species=="Bare Sand"|
             species=="Unidentified Fish"|
             species=="Actiniaria"|
             species=="Shell Debris"|
             species=="Sediment/Mud"|
             species=="Stylantheca papillosa"|
             species=="Diaperoforma californica"|
             species==""))


#recalculate percent cov
kelp_upc_build6 <- kelp_upc_build5 %>%
  group_by(year, site, zone, transect, classcode, species)%>%
  dplyr::summarize(sum_count = sum(count)) 

#total points per transect, taken BEFORE dropping non-biotic/excluded categories
#(bare rock, sand, etc.) so cover is relative to all points. Uses build1 (no taxon join)
#to avoid double counting classcodes that map to multiple taxa (e.g., UNID).
#Most transects should total ~30 points.
check_counts <- kelp_upc_build1 %>%
  group_by(year, site, zone, transect)%>%
  dplyr::summarize(transect_total = sum(count), .groups = "drop")

#join percent cov with transect total
kelp_upc_build7 <- left_join(kelp_upc_build6, check_counts, by=c("year","site","zone","transect"))%>%
  mutate(pct_cov = (sum_count / transect_total)*100)%>%
  mutate_at(vars(pct_cov), round, 1)



#add affiliated_mpa

kelp_upc_build8 <- left_join(kelp_upc_build7, site_table, by="site")

kelp_upc_build9 <- kelp_upc_build8 %>%
  ungroup() %>%
  dplyr::select(year, baseline_region, site,latitude, longitude, affiliated_mpa=ca_mpa_name_short,mpa_class, mpa_designation, zone, transect, species,
                #sum_count, transect_total,
                pct_cov) %>%
  mutate(mpa_designation = ifelse(mpa_designation=="reference","ref",mpa_class)) 



#reshape to wide format
kelp_upc_build10 <- kelp_upc_build9 %>%
  pivot_wider(names_from=species, values_from=pct_cov) %>%
  janitor::clean_names()%>%
  #drop_na(region4) %>%
  mutate(MHW = ifelse(year>=2014 & year<=2016, "during",ifelse(year<2014, "before","after")))%>%
  #dplyr::select(year, region3, region4, affiliated_mpa, mpa_defacto_class, MHW, everything())%>%
  #filter central coast only
  filter(baseline_region == "CENTRAL") %>%
  #filter years >= 2007
  filter(year >= 2007) %>%
  #absent taxa are 0 cover; fill by name rather than column position
  mutate(across(-c(year, baseline_region, site, latitude, longitude, affiliated_mpa,
                   mpa_class, mpa_designation, zone, transect, MHW),
                ~replace_na(.x, 0))) %>%
  dplyr::select(year, MHW, everything())
#filter(region3 == 'central') %>%
#select(where(~ any(. != 0)))
#mutate_all(~ifelse(is.nan(.), NA, .))


################################################################################
#select for carmel bay and big sur sites

nereocystis_sites <- kelp_swath_build8 %>%
  group_by(site, latitude, longitude) %>%
  summarise(
    n_years = n_distinct(year),
    years_with_nereocystis = n_distinct(
      year[nereocystis_luetkeana > 0]
    ),
    max_nereocystis = max(nereocystis_luetkeana, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(years_with_nereocystis > 0) %>%
  arrange(desc(latitude))

print(nereocystis_sites, n = Inf)

library(sf)
land <- rnaturalearth::ne_countries(scale = "large", country = "United States of America",
                                    returnclass = "sf")

nereocystis_sites |>
  # Draw sites with the most years on top
  arrange(years_with_nereocystis) |>
  ggplot() +
  geom_sf(data = land, fill = "grey90", color = "grey50", linewidth = 0.3) +
  geom_point(aes(longitude, latitude, color = years_with_nereocystis),
             size = 3, alpha = 0.85) +
  # Label only the top sites, placed offshore to avoid the dense cluster
  ggrepel::geom_text_repel(
    data = ~ slice_max(.x, years_with_nereocystis, n = 5, with_ties = FALSE),
    aes(longitude, latitude, label = site), size = 3, min.segment.length = 0,
    xlim = c(NA, -122.3), direction = "y", hjust = 1, segment.color = "grey40") +
  scale_color_viridis_c(direction = -1, end = 0.9) +
  coord_sf(xlim = c(-123.6, -120.4), ylim = c(34.6, 37.1), expand = FALSE) +
  theme_classic() +
  theme(panel.grid = element_line(color = "grey92")) +
  labs(x = NULL, y = NULL, color = "Years with\nbull kelp",
       title = "Sites with historical Nereocystis presence")

################################################################################
#select for Carmel Bay and Big Sur MPA sites
# Define geographic extent
################################################################################
# Select Carmel Bay and Big Sur sites

target_mpas <- c(
  "Carmel Bay SMCA",
  "Carmel Pinnacles SMR",
  "Point Lobos SMR",
  "Point Sur SMR",
  "Big Creek SMR"
)

# Filter SWATH data
kelp_swath_CB_BS <- kelp_swath_build8 %>%
  filter(affiliated_mpa %in% target_mpas)

# Filter UPC data
kelp_upc_CB_BS <- kelp_upc_build10 %>%
  filter(affiliated_mpa %in% target_mpas)

################################################################################
# Evaluate selected sites

site_summary <- kelp_swath_CB_BS %>%
  group_by(affiliated_mpa, site, latitude, longitude) %>%
  summarise(
    first_year = min(year),
    last_year = max(year),
    n_years = n_distinct(year),
    n_transects = n(),
    n_years_bull_kelp = n_distinct(
      year[nereocystis_luetkeana > 0]
    ),
    mean_bull_kelp = mean(nereocystis_luetkeana),
    mean_giant_kelp = mean(macrocystis_pyrifera),
    mean_giant_kelp_stipes = mean(macrocystis_stipes),
    .groups = "drop"
  ) %>%
  mutate(
    pct_years_bull_kelp = 100 * n_years_bull_kelp / n_years
  ) %>%
  arrange(desc(latitude))

ggplot(site_summary, aes(x = longitude, y = latitude)) +
  geom_point(
    aes(color = pct_years_bull_kelp, size = n_years)
  ) +
  geom_text_repel(aes(label = site), size = 3) +
  scale_color_gradient(
    low = "grey80",
    high = "#6E8B3D",
    limits = c(0, 100)
  ) +
  coord_quickmap() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    color = "Years with\nbull kelp (%)",
    size = "Years sampled",
    title = "Carmel Bay and Big Sur monitoring sites"
  )
site_order <- site_summary %>%
  arrange(latitude) %>%
  pull(site) %>%
  unique()

sampling_summary <- kelp_swath_CB_BS %>%
  group_by(site, affiliated_mpa, year) %>%
  summarise(
    n_transects = n(),
    bull_kelp_present = any(nereocystis_luetkeana > 0),
    .groups = "drop"
  ) %>%
  mutate(
    site = factor(site, levels = site_order)
  )

ggplot(sampling_summary, aes(x = year, y = site)) +
  geom_tile(
    aes(fill = bull_kelp_present),
    color = "white",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = c("FALSE" = "grey80", "TRUE" = "#6E8B3D"),
    labels = c("Absent", "Present")
  ) +
  scale_x_continuous(
    breaks = seq(2001, max(sampling_summary$year), 1)
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Survey year",
    y = NULL,
    fill = "Bull kelp",
    title = "Temporal sampling coverage",
    subtitle = "Sites ordered north to south; blank cells indicate unsampled years")
################################################################################
# Export

write.csv(
  kelp_swath_CB_BS,
  file.path(output, "kelp_swath_counts_CB_BS.csv"),
  row.names = FALSE)

write.csv(
  kelp_upc_CB_BS,
  file.path(output, "kelp_upc_cov_CB_BS.csv"),
  row.names = FALSE)

################################################################################
#for carmel bay and monterey bay sites
kelp_swath_build9 <- kelp_swath_build8 %>%
  #select sites in Carmel and Monterey Bay only
  dplyr::filter(latitude >= 36.46575 & latitude <= 36.64045) %>%
  filter(!(site == "ASILOMAR_DC" |
             site == "ASILOMAR_UC" |
             site == "CHINA_ROCK" |
             site == "CYPRESS_PT_DC" |
             site == "CYPRESS_PT_UC" |
             site == "PINNACLES_IN" |
             site == "PINNACLES_OUT" |
             site == "PT_JOE" |
             site == "SPANISH_BAY_DC" |
             site == "SPANISH_BAY_UC" |
             site == "BIRD_ROCK"|
             site == "LINGCOD_DC"|
             site == "LINGCOD_UC"|
             site == "CARMEL_DC"|
             site == "CARMEL_UC"))



kelp_upc_build11 <- kelp_upc_build10 %>%
  #select sites in Carmel and Monterey Bay only
  dplyr::filter(latitude >= 36.46575 & latitude <= 36.64045) %>%
  filter(!(site == "ASILOMAR_DC" |
             site == "ASILOMAR_UC" |
             site == "CHINA_ROCK" |
             site == "CYPRESS_PT_DC" |
             site == "CYPRESS_PT_UC" |
             site == "PINNACLES_IN" |
             site == "PINNACLES_OUT" |
             site == "PT_JOE" |
             site == "SPANISH_BAY_DC" |
             site == "SPANISH_BAY_UC" |
             site == "BIRD_ROCK"|
             site == "LINGCOD_DC"|
             site == "LINGCOD_UC"|
             site == "CARMEL_DC"|
             site == "CARMEL_UC"))

################################################################################
#export 


#Export
write.csv(kelp_swath_build9,file.path(output, "kelp_swath_counts_CC.csv"), row.names = FALSE)
#last write April 22, 2026

#Export
write.csv(kelp_upc_build11,file.path(output, "kelp_upc_cov_CC.csv"), row.names = FALSE)
#last write April 22, 2026

