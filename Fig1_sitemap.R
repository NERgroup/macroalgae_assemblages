##Fig1_site_map
#Sabrina N. Grant (code adapted from Joshua G. Smith's script)
#June 30th, 2026
 
#install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")

rm(list=ls())

librarian::shelf(patchwork, ggplot2, ggrepel, ggspatial, stringr, tidyverse, here, tidyverse, tidync, sf, raster, terra, ggtext, rnaturalearth, cowplot)


################################################################################
#set directories and load data
basedir <- here::here("output")
figdir <- here::here("analyses")

################################################################################
#shared components of both maps (land polygons, theme, CA locator inset)

usa     <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf")
foreign <- rnaturalearth::ne_countries(country = c("Canada", "Mexico"), returnclass = "sf")
ca_counties <- usa %>% filter(name == "California")

bbox <- st_bbox(c(xmin = -122.08, ymin = 36.4, xmax = -121.7, ymax = 36.670574),
                crs = st_crs(ca_counties))

ca_counties_mpen <- st_intersection(ca_counties, st_as_sfc(bbox))

landmarks_a <- tibble(place = "Monterey \nPeninsula", long = -121.927984, lat = 36.585)
landmarks_b <- tibble(place = "Monterey \nPeninsula", long = -121.93, lat = 36.603)

monterey_label <- data.frame(
  x = c(-121.89, -122.03),
  y = c(36.657, 36.54),
  label = c("Monterey \nBay", "Carmel \nBay")
)

my_theme <- theme(
  axis.text = element_text(size = 6),
  axis.text.y = element_text(angle = 90, hjust = 0.5),
  axis.title = element_text(size = 8),
  plot.tag = element_text(size = 8, face = "bold"),
  plot.title = element_text(size = 7, face = "bold"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  legend.key = element_blank(),
  legend.background = element_rect(fill = alpha('blue', 0)),
  legend.key.height = unit(1, "lines"),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 7),
  strip.background = element_blank(),
  strip.text = element_text(size = 6, face = "bold"),
  panel.background = element_rect(fill = "lightblue")
)

g1_inset <- ggplotGrob(
  ggplot() +
    geom_sf(data = foreign, fill = "grey80", color = "white", lwd = 0.3) +
    geom_sf(data = usa, fill = "grey80", color = "white", lwd = 0.3) +
    annotate("rect", xmin = -122.6, xmax = -121, ymin = 36.2, ymax = 37.1,
             color = "black", fill = NA, lwd = 0.8) +
    labs(x = "", y = "") +
    coord_sf(xlim = c(-124.5, -117), ylim = c(32.5, 42)) +
    theme_bw() + my_theme +
    theme(
      plot.margin = unit(rep(0, 4), "null"),
      panel.background = element_rect(fill = 'transparent'),
      axis.ticks = element_blank(),
      axis.ticks.length = unit(0, "null"),
      axis.text = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title = element_blank()
    )
)

################################################################################
#Panel A: MPLA MAP
#real site coordinates pulled from kelp_swath_counts_CC.csv
################################################################################

#site locations 
site_locations_labels_a <- tribble(
  ~site,            ~longitude,   ~latitude,
  "Otter Pt DC",    -121.921162,  36.636484,
  "Otter Pt UC",    -121.918929,  36.634633,
  "Siren",          -121.918590,  36.630733,
  "Lovers DC",      -121.910735,  36.625774,
  "Lovers UC",      -121.909082,  36.624220,
  "Hopkins DC",     -121.904196,  36.623586,
  "Hopkins UC",     -121.900789,  36.621649,
  "Macabee DC",     -121.896835,  36.618184,
  "Macabee UC",     -121.895695,  36.617148,
  "Cannery DC",     -121.896037,  36.614953,
  "Cannery UC",     -121.894574,  36.612639,
  "Lone Tree",      -121.970688,  36.566418,
  "Pescadero UC",   -121.959761,  36.561121,
  "Stillwater UC",  -121.947320,  36.560116,
  "Stillwater DC",  -121.944014,  36.559729,
  "Pescadero DC",   -121.955191,  36.559412,
  "Butterfly UC",   -121.936297,  36.539655,
  "Butterfly DC",   -121.935670,  36.537439,
  "Monastery DC",   -121.933322,  36.525415,
  "Monastery UC",   -121.930507,  36.525267,
  "Bluefish DC",    -121.945094,  36.523813,
  "Bluefish UC",    -121.943542,  36.522225,
  "Weston UC",      -121.948765,  36.512613,
  "Weston DC",      -121.945786,  36.510347
)

site_locations_sf_a <- site_locations_labels_a %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

#mpla sites map
map_a <- ggplot() +
  geom_sf(data = ca_counties_mpen) +
  geom_sf(data = site_locations_sf_a) +
  ggrepel::geom_label_repel(
    data = site_locations_labels_a,
    aes(x = longitude, y = latitude, label = site),
    box.padding = 0.3, point.padding = 0.5, force = 18,
    size = 2, min.segment.length = 0.1, segment.color = "grey50",
    max.overlaps = Inf  #fix: without this, ggrepel was silently dropping
    #Lovers/Hopkins/Macabee DC labels in the dense cluster
  ) +
  geom_text(data = landmarks_a, mapping = aes(x = long, y = lat, label = place),
            size = 4, fontface = "bold") +
  coord_sf(xlim = c(-122.08, -121.85), ylim = c(36.5, 36.66)) +
  theme_bw() + my_theme +
  xlab("Longitude") + ylab("Latitude") +
  annotation_custom(grob = g1_inset, xmin = -122.08, xmax = -121.997, ymin = 36.615, ymax = 36.665) +
  annotation_scale(location = "br", width_hint = 0.25, style = "ticks") +
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(.5, "cm"), pad_y = unit(0.5, "cm"),  #fix: pad_x was 0.5cm, not
    #enough to clear the inset --
    #needs ~8.5cm at this panel width
    style = north_arrow_fancy_orienteering(),
    height = unit(2, "cm"), width = unit(2, "cm")
  ) +
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 4, fontface = "bold") +
  annotate("segment", x = -122, y = 36.54, xend = -121.96, yend = 36.54,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))

################################################################################
#Panel B: MBA Map
##real site coordinates pulled from kelp_quad_data_CC.csv

site_locations_labels_b <- tribble(
  ~site,    ~zone,      ~longitude,   ~latitude,
  "Rec 6",  "Deep",     -121.924711,  36.635513,
  "Rec 6",  "Shallow",  -121.925471,  36.635311,
  "Rec 7",  "Deep",     -121.905628,  36.622197,
  "Rec 7",  "Shallow",  -121.905469,  36.622116,
  "Rec 1",  "Shallow",  -121.968783,  36.568807,
  "Rec 1",  "Deep",     -121.968767,  36.568195,
  "Rec 2",  "Shallow",  -121.955704,  36.561829,
  "Rec 2",  "Deep",     -121.956116,  36.561439,
  "Rec 3",  "Shallow",  -121.943589,  36.561363,
  "Rec 14", "Shallow",  -121.945361,  36.561260,
  "Rec 13", "Deep",     -121.947136,  36.561251,
  "Rec 12", "Shallow",  -121.945425,  36.561096,
  "Rec 12", "Deep",     -121.945838,  36.560686,
  "Rec 3",  "Deep",     -121.943653,  36.560643,
  "Rec 11", "Shallow",  -121.935624,  36.555265,
  "Rec 11", "Deep",     -121.937602,  36.554974,
  "Rec 4",  "Deep",     -121.935201,  36.548574,
  "Rec 4",  "Shallow",  -121.934696,  36.548096,
  "Rec 5",  "Shallow",  -121.934434,  36.540208,
  "Rec 5",  "Deep",     -121.935253,  36.539736,
  "Rec 10", "Shallow",  -121.928712,  36.532849,
  "Rec 10", "Deep",     -121.929212,  36.532527
) %>%
  mutate(site_zone = paste(site, zone))

site_locations_sf_b <- site_locations_labels_b %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

map_b <- ggplot() +
  geom_sf(data = ca_counties_mpen) +
  geom_sf(data = site_locations_sf_b, aes(shape = zone, color = zone), size = 2) +
  scale_shape_manual(values = c("Deep" = 1, "Shallow" = 17)) +
  scale_color_manual(values = c("Deep" = "black", "Shallow" = "grey40")) +
  ggrepel::geom_label_repel(
    data = site_locations_labels_b,
    aes(x = longitude, y = latitude, label = site_zone),
    box.padding = 0.3, point.padding = 0.5, force = 40,  #fix: bumped from 25 --
    #Rec 11-14 cluster near
    #Carmel Bay still crowded,
    #may need a zoomed inset
    #if this isn't enough
    size = 1.8, min.segment.length = 0.1, segment.color = "grey50",
    max.overlaps = Inf, max.iter = 20000
  ) +
  geom_text(data = landmarks_b, mapping = aes(x = long, y = lat, label = place),
            size = 4, fontface = "bold") +
  coord_sf(xlim = c(-122.08, -121.85), ylim = c(36.5, 36.66)) +
  theme_bw() + my_theme +
  xlab("Longitude") + ylab("Latitude") +
  labs(shape = "Zone", color = "Zone") +  #fix: labs(label = "Zone") did nothing --
  #"label" isn't a guide-producing aesthetic
  #here, so the legend title stayed "zone"
  theme(legend.position = c(0.1, 0.9)) + 
  theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9, face = "bold")
  ) +
  annotation_scale(location = "br", width_hint = 0.25, style = "ticks") +
  geom_text(data = monterey_label, mapping = aes(x = x, y = y, label = label),
            size = 4, fontface = "bold") +
  annotate("segment", x = -122, y = 36.54, xend = -121.96, yend = 36.54,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))

################################################################################
#Combine panels
################################################################################

combined <- map_a + map_b


combined <- plot_grid(
  map_a, map_b,
  labels = c("A", "B"),
  label_size = 15,
  label_fontface = "bold",
  ncol = 2,
  align = "h"
)

combined


ggsave(combined, filename=file.path(figdir, "Fig1_site_map.png"), bg = "white",
       width=10, height=8, units="in", dpi=600) 
