#Figure 3
#Sabrina N. Grant
#July 9th, 2026

rm(list=ls())
librarian::shelf(tidyverse, here, ggplot2, MASS, caret, patchwork, grid, RColorBrewer, reshape2, vegan, ggtext)


##require:
#normalized_kelp_quad_data_CC.csv
#kelp_quad_data_CC.csv

################################################################################
#set directories and load data
basedir     <- here::here("output")
figures_out <- here::here("figures")

mba_data <- read_csv(file.path(basedir, "mba_data", "processed", "kelp_quad_data_CC.csv"))%>%
  janitor::clean_names()
norm_mba_data <- read_csv(file.path(basedir, "mba_data", "processed", "normalized_kelp_quad_data_CC.csv"))%>%
  janitor::clean_names()



################################################################################
#renaming dfs

merged_data <- mba_data

norm_data <- norm_mba_data  


#filtering for 2024 data 
year_one_data <- merged_data %>% 
  filter(year == 2024) %>%
  dplyr::select(year, site, zone, site_type, macro_stipe_density, everything()) %>%
  mutate(site_type = as.factor(site_type))


#getting mean of numeric values across grouping categories
year_one_data %>%
  group_by(site, site_type, zone, year) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop")

#filtering for 2024 z-scored data 
z_year_one_data <- norm_data %>% 
  filter(year == 2024) %>%
  dplyr::select(year, site, zone, site_type, macro_stipe_density, everything()) %>%
  mutate(site_type = as.factor(site_type)) 

#getting mean of numeric values across grouping categories
z_year_one_data_means <-  z_year_one_data %>%
  group_by(site, site_type, zone, year) %>%
  summarise(
    across(where(is.numeric), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop")

################################################################################
#LDA on Z-Scored Data

lda_model_z <- lda(site_type ~ macro_stipe_density + 
                     density_ptecal + density_nerlue + 
                     density_lamset + density_eisarb + 
                     density_lamstump + 
                     density_macstump + 
                     lamr + 
                     macr +
                     macj + 
                     nerj + 
                     ptej + 
                     lsetj + 
                     cov_articulated_coralline + 
                     cov_crustose_coralline + 
                     cov_encrusting_red + 
                     cov_fleshy_red + 
                     cov_stephanocystis + 
                     cov_dictyoneurum_spp  + 
                     cov_desmarestia_spp + 
                     cov_lam_holdfast_live + 
                     cov_mac_holdfast_live, 
                   data = z_year_one_data)

#lda summary
lda_model_z

#getting LDA predictions
pred_z <- predict(lda_model_z)
pred_class_z <- pred_z$class

#getting % accuracy
confusionMatrix(pred_class_z, z_year_one_data$site_type)
#predicts at 81% accuracy which is the same as non z scored


################################################################################
#Fig. 7 A: LDA model B: Coefficients 

#A: LDA Model Fig
#getting LDA scores 
lda_scores_z <- as.data.frame(predict(lda_model_z)$x) 
lda_scores_z$site_type <- z_year_one_data$site_type 

#plotting LD1 vs LD2 
a7 <- ggplot(lda_scores_z, aes(x = LD1, y = LD2, color = site_type, shape = site_type)) + geom_point(size = 3.5, alpha = 0.9) + 
  stat_ellipse(level = 0.95, size = 1, alpha = 0.6) + 
  scale_color_brewer(palette = "Set1") + 
  theme_classic(base_size = 16) + 
  labs( title = "", subtitle = "", x = paste0("LD1 (", round(lda_model_z$svd[1]^2 / sum(lda_model_z$svd^2) * 100, 1), "%)"), y = paste0("LD2 (", round(lda_model_z$svd[2]^2 / sum(lda_model_z$svd^2) * 100, 1), "%)"), color = "Site Type", shape = "Site Type" ) + 
  theme( plot.title = element_text(face = "bold", size = 20, hjust = 0.5, margin = margin(b = 10)), plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)), axis.title = element_text(face = "bold", size = 15), axis.text = element_text(size = 13), legend.title = element_text(face = "bold", size = 14), legend.text = element_text(size = 12), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "white", linewidth = 0.3) )

#B: Coefficients of LDA fig
# Get coefficients into tidy format
coef_z <- as.data.frame(lda_model_z$scaling)
coef_z$variable <- rownames(coef_z)
coef_long <- pivot_longer(coef_z, cols = c(LD1, LD2), names_to = "LD", values_to = "coefficient")

# clean up variable names for display
coef_long <- coef_long %>%
  mutate(variable = str_replace_all(variable, "_", " ") %>%
           str_to_title())

b7 <- ggplot(coef_long, aes(x = coefficient, 
                            y = reorder(variable, coefficient), 
                            fill = LD)) +
  geom_col(position = position_dodge(width = 0.7), 
           width = 0.6, 
           alpha = 0.85) +
  geom_vline(xintercept = 0, linewidth = 0.7, color = "black", linetype = "solid") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic(base_size = 16) +
  labs(
    x = "Discriminant Coefficient",
    y = NULL,
    fill = "Linear\nDiscriminant"
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 13),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

a7 
b7
fig7 <- a7 + b7 +
  plot_layout(ncol = 2, widths = c(1, 1.4)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 18)
    )
  )
fig7

################################################################################

ggsave(file.path(figures_out, "Fig3_LDA_and_Discriminants.png"), fig7, bg = "white",
       width = 10, height = 8, units = "in", dpi = 600)

