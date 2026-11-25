#Distribution of urchin gonad index (individual-level dissection data)
#Sabrina N. Grant
#November 24, 2026

librarian::shelf(tidyverse, here, patchwork)

################################################################################
#set directories and load data

procdir <- here::here("output")
figures_out <- here::here("figures")

dissection_data <- read_csv(file.path(procdir, "dissection_data", "processed",
                                      "dissection_data_recovery_clean.csv"),
                            show_col_types = FALSE) %>%
  mutate(site_type = factor(site_type, levels = c("BAR", "INCIP", "FOR")),
         species = factor(species, levels = c("purple_urchin", "red_urchin"),
                          labels = c("Purple urchin", "Red urchin")))

################################################################################
#summary statistics

gonad_summary <- dissection_data %>%
  group_by(species, year, site_type) %>%
  summarise(n = n(),
            n_missing = sum(is.na(gonad_index)),
            mean = mean(gonad_index, na.rm = TRUE),
            sd = sd(gonad_index, na.rm = TRUE),
            median = median(gonad_index, na.rm = TRUE),
            q25 = quantile(gonad_index, 0.25, na.rm = TRUE),
            q75 = quantile(gonad_index, 0.75, na.rm = TRUE),
            min = min(gonad_index, na.rm = TRUE),
            max = max(gonad_index, na.rm = TRUE),
            prop_zero = mean(gonad_index == 0, na.rm = TRUE),
            .groups = "drop")
gonad_summary

################################################################################
#plots

#overall distribution by species
p_species <- ggplot(dissection_data, aes(x = gonad_index)) +
  geom_histogram(bins = 40, fill = "grey70", color = "black", linewidth = 0.2) +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Gonad index (%)", y = "Count") +
  theme_classic(base_size = 14)

#distribution by species and patch type (years pooled; see gonad_summary for year)
p_group <- ggplot(dissection_data, aes(x = gonad_index)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black", linewidth = 0.2) +
  facet_grid(species ~ site_type, scales = "free_y") +
  labs(x = "Gonad index (%)", y = "Count") +
  theme_classic(base_size = 14)

combined <- p_species / p_group + plot_layout(heights = c(1, 1.5))
combined

################################################################################
#export

ggsave(file.path(figures_out, "FigS3_gonad_index_distribution.png"), combined,
       bg = "white", width = 8, height = 8, units = "in", dpi = 600)
