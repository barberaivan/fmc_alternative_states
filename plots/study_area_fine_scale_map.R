library(terra)
library(sf)
library(ggplot2)
library(ggspatial)
library(ggmap)
register_google(key = "AIzaSyAcqcfyYbP-7qSeweEf2wnF9blVV_Ir7Hc", write = TRUE)

# argentina map
library(rnaturalearth)
library(rnaturalearthdata)

# Clave maps
# AIzaSyAcqcfyYbP-7qSeweEf2wnF9blVV_Ir7Hc

fires <- vect(
  "/home/ivan/Insync/patagonian_fires/patagonian_fires/patagonian_fires.shp"
)
fire <- fires[fires$name == "Balcon del Gutiérrez", ]

# Convert to sf if it's a terra object
if (inherits(fire, "SpatVector")) {
  fire <- st_as_sf(fire)
}

# Define bounding box manually: area between Lago Nahuel Huapi (N) and Lago Gutiérrez (S)
# Approximate bounding box around Bariloche
bbox <- st_as_sfc(st_bbox(c(
  xmin = -71.35,
  xmax = -71.25,
  ymin = -41.2,
  ymax = -40.95
), crs = 4326))


# Fetch the map (adjust zoom until you're happy)
map <- get_googlemap(
  center = c(lon = -71.33, lat = -41.1),
  zoom = 10,
  maptype = "satellite",#"roadmap", #"terrain",
  style = "feature:all|element:labels|visibility:off"  # removes all labels
)

# Plot
mapa <- 
ggmap(map, extent = "device") +
  
  # Add white transparent layer to fade the map
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -Inf, ymax = Inf,
           fill = "white", alpha = 0.1) +  # adjust alpha to control fading
  
  # Fire
  geom_sf(data = fire, inherit.aes = FALSE, fill = "red", alpha = 0.7) +
  
  # Custom position: Scale bar near bottom center
  annotation_scale(
    location = "br",  # still needed for units
    width_hint = 0.5,
    line_width = 0.5, height = unit(0.15, "cm"),
    bar_cols = c("gray30", "white")  # two segments
    #pad_x = unit(0.25, "npc"), pad_y = unit(0.02, "npc")
  ) +
  
  scale_y_continuous(breaks = c(-41.0, -41.2)) +
  scale_x_continuous(breaks = c(-71.4, -71.2)) +
  
  theme(
    panel.grid.major = element_line(color = "gray80", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 9),
    axis.title = element_blank()
  ) +
  coord_sf(expand = FALSE)
mapa

ggsave(
  "mapa_local.pdf",
  width = 16, height = 19, units = "cm", dpi = 900
)



# Argentina map -----------------------------------------------------------

# Get country borders
arg <- ne_countries(scale = "medium", country = "Argentina", returnclass = "sf")

# Create point for Bariloche area
bariloche <- st_as_sf(data.frame(
  lon = -71.33,
  lat = -41.1
), coords = c("lon", "lat"), crs = 4326)

# Plot}
mapa_arg <- 
ggplot() +
  geom_sf(data = arg, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = bariloche, color = "blue", size = 2) +
  coord_sf(
    xlim = c(-75, -52), 
    ylim = c(-56, -21), 
    expand = FALSE
  ) +
  theme_void()
mapa_arg

ggsave(
  "mapa_arg.jpeg", plot = mapa_arg,
  width = 3, height = 6, units = "cm", dpi = 500
)
