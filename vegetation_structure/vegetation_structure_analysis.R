library(tidyverse)   # everything
library(mgcv)        # gam
library(plyr)        # revalue
library(vegan)       # metaMDS
library(ggpubr)      # ggarrange
library(viridis)     # nice colors
library(desiderata)  # rotate_x_facet_text, to left-aling facet text

# Prepare data ------------------------------------------------------------

es <- read.csv("vegetation_structure/data_vegetation_structure.csv", header = TRUE, sep = ",",
               stringsAsFactors = FALSE) 
 
# remove "Matorral de caña" and transect 5 (we didn't have dataloggers there)
es <- es[es$comm_details != "Matorral de cania" & 
         es$transect != 5, ]

# rename transects and communities
unique(es$transect)
es$transect <- revalue(as.character(es$transect), 
                       replace = c("2" = "T1",
                                   "4" = "T2",
                                   "6" = "T3"))
es$transect <- factor(es$transect, levels = c("T1", "T2", "T3"))

unique(es$Community)
es$Community <- revalue(as.character(es$Community), 
                       replace = c("Coihual no regenerado" = "Open shrubland",
                                   "Coihual regenerando" = "Young forest",
                                   "Matorral" = "Closed shrubland",              
                                   "Bosque no quemado" = "Mature forest"))
es$Community <- factor(es$Community, levels = c("Mature forest", "Young forest", "Closed shrubland", "Open shrubland"))


es$site <- paste(es$Community, es$transect, sep = " / ")


# recode point_id:
es$point_id <- as.numeric(as.factor(es$point_id))


# function to separate site into community and transect
separate_site <- function(d) {
  if("site" %in% names(d)) {
    d <- separate(d, "site", into = c("Community", "transect"), sep = " / ",
                  remove = FALSE)
    d$Community <- factor(d$Community, levels = c("Mature forest", "Young forest", "Closed shrubland", "Open shrubland"))
    d$transect <- factor(d$transect, levels = c("T1", "T2", "T3"))
  } else {stop("No site column in this dataset.")}
  return(d)
}



# Vertical structure ------------------------------------------------------

# necesito primero llenar una matriz por tipo de fuel (vivo, muerto, total)
# que tenga tantas cols como intercepts y tantas rows como points_id.
# Allí pondré 1 o 0 si el string de un punto contiene el número de cada columna.

# OJO, esto de abajo puede estar MAL.
# parece que puse point_id repetidos en distintos sitios. 
# O sea, el point_id = 1 puede estar en varios sitios, por ejemplo.

puntos <- max(es$point_id)
intMat <- vector(mode = "list", length = 3)
for(i in 1:3) intMat[[i]] <- matrix(0, puntos, 16)
names(intMat) <- c("Live", "Dead", "Live or dead")

# Fill live and dead
names(es)
columns <- c("vivo_interc", "muerto_interc")
for(f in 1:2) {
  for(i in 1:puntos) {
    for(int in 1:16) {
      ind <- which(es$point_id == i)
      texto <- paste(es[ind, columns[f]], sep = "<")
      tokes <- as.numeric(do.call("c", strsplit(texto, "<")))
      toca <- as.numeric(int %in% tokes)
      intMat[[f]][i, int] <- toca  
    }
  }
}

# Fill total
for(i in 1:puntos) {
  for(int in 1:16) {
    ind <- which(es$point_id == i)
    toca <- numeric(2)
    for(col in 15:16) {
      texto <- paste(es[ind, col], sep = "<")
      tokes <- as.numeric(do.call("c", strsplit(texto, "<")))
      toca[col - 14] <- as.numeric(int %in% tokes)
    }
    intMat[["Live or dead"]][i, int] <- max(toca)  
  }
}

str(intMat)
# unique(as.numeric(intMat$vivo))
# unique(as.numeric(intMat$muerto))
# unique(as.numeric(intMat$todes))

# Add variables
temp <- aggregate(rep(1, nrow(es)) ~ Community + transect + site + point_id, data = es,
                  FUN = mean)[, 1:4]
for(i in 1:3) {
  intMat[[i]] <- cbind(temp, intMat[[i]]) %>% as.data.frame
  intMat[[i]]$fuel_type <- names(intMat)[i]
}

# paso a long
intMat2 <- do.call("rbind", intMat)

int_names <- as.character(1:16)
intmat_long <- pivot_longer(intMat2, cols = which(names(intMat2) %in% int_names),
                            names_to = "intercept", values_to = "touch")

sum_len <- function(x) return(c(suma = sum(x), largo = length(x)))
int_agg <- aggregate(touch ~ Community + transect + site + fuel_type + intercept,
                     intmat_long, FUN = sum_len)
str(int_agg)
int_agg <- do.call("data.frame", int_agg)
names(int_agg)[(ncol(int_agg)-1) : ncol(int_agg)] <- c("touch", "total")

# define height:
int_agg$h <- as.numeric(int_agg$intercept) / 4 - 0.25 / 2

# define site * fuel_type for the GAM
int_agg$fuel_site <- as.factor(paste(int_agg$fuel_type, int_agg$site, sep = " / "))

# Compute observed proportions to plot
int_agg$prop <- int_agg$touch / int_agg$total


#### GAMs fit 

# VIVO
mgam <- gam(cbind(touch, total - touch) ~ 
              fuel_site + s(h, by = fuel_site, bs = "cr", k = 12),
            method = "REML", family = "binomial", data = int_agg)

# data for prediction
nd <- expand.grid(
  h = seq(0, 4, length.out = 150),
  fuel_site = unique(int_agg$fuel_site)
)

preds <- predict(mgam, newdata = nd, se.fit = TRUE)
nd$p_mle <- plogis(preds$fit)
nd$p_lower <- plogis(preds$fit - qnorm(0.975) * preds$se.fit)
nd$p_upper <- plogis(preds$fit + qnorm(0.975) * preds$se.fit)

nd <- separate(nd, "fuel_site", into = c("fuel_type", "Community", "transect"),
               sep = " / ", remove = F)

nd$fuel_type <- factor(nd$fuel_type, levels = c("Live", "Dead", "Live or dead"))
nd$Community <- factor(nd$Community, levels = levels(es$Community))
int_agg$fuel_type <- factor(int_agg$fuel_type, levels = c("Live", "Dead", "Live or dead"))

# plot
theme_set(theme_bw())
ggplot(data = nd, aes(x = h, y = p_mle, ymin = p_lower, ymax = p_upper, 
                      colour = Community, fill = Community)) +
  geom_line() + 
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_point(data = int_agg, 
             mapping = aes(x = h, y = prop, colour = Community),
             inherit.aes = FALSE,
             shape = 1
             ) +
  facet_grid(rows = vars(fuel_type), cols = vars(transect)) + 
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank()) +
  coord_flip() +
  xlab("Height (m)") +
  ylab("Fuel presence probability or proportion") 
  

# Fuel proportion (flatten vertical dimension) ------------------------------

fuel_prop <- aggregate(cbind(touch, total) ~ Community + transect + fuel_type,
                       int_agg, FUN = sum)
str(fuel_prop)
fuel_prop$p_obs <- fuel_prop$touch / fuel_prop$total

# sd for the proportion, from Gelman and Hill 2007, page 17
fuel_prop$p_sd <- sqrt(fuel_prop$p_obs * (1 - fuel_prop$p_obs) / fuel_prop$total)

# confidence interval:
fuel_prop$p_lower <- fuel_prop$p_obs - qnorm(0.975) * fuel_prop$p_sd
fuel_prop$p_upper <- fuel_prop$p_obs + qnorm(0.975) * fuel_prop$p_sd

# plot
w <- 0.3
ggplot(data = fuel_prop, 
       aes(x = transect, y = p_obs, ymin = p_lower, ymax = p_upper, 
           colour = Community)) +
  geom_point(position = position_dodge2(width = w),
             size = 3) +
  geom_linerange(position = position_dodge2(width = w),
                 size = 0.7) +
  facet_wrap(vars(fuel_type), nrow = 2, scales = "free_y") + 
  theme(legend.position = c(0.75, 0.25),
        panel.grid.minor = element_blank()) +
  xlab("Transect") +
  ylab("Fuel presence proportion") 



# Species composition -----------------------------------------------------

# por sitio, por especie, 
# length(unique(c(texto muerto[>0], texto vivo[>0])))
# hacer matriz de tantas filas como puntos y tantas cols como spp

names(es)

species <- unique(es$sp_code)[!is.na(unique(es$sp_code))]
nsp <- length(species)

es2 <- es[complete.cases(es$sp_code), ]
sites <- unique(es2$site)
nsites <- length(unique(es2$site))

# Create array with number of intercept touches by species
# for live, dead and both.

sp_touch <- array(0, dim = c(nsites, nsp, 3),
                  dimnames = list(
                    site = sites,
                    sp = species,
                    fuel_type = c("Live", "Dead", "Live or dead")
                  ))

for(s in sites) {
  print(s)
  for(sp in species) {
    # s = "Closed shrubland / T1"
    # sp = "Nant"
    filter <- es2$sp_code == sp & es2$site == s
    temp <- es2[filter, c("vivo_interc", "muerto_interc")]
  
    if(nrow(temp) > 0) {
      temp$tot <- 0
      # loop through points where the species was found to count number of 
      # touches in each point
      
      for(i in 1:nrow(temp)) {
        v <- temp[i, "vivo_interc"]
        m <- temp[i, "muerto_interc"]
        tocv <- as.numeric(do.call("c", strsplit(v, "<")))
        tocm <- as.numeric(do.call("c", strsplit(m, "<")))
        total <- unique(c(tocv, tocm))
      
        temp$vivo_interc[i] <- length(tocv[tocv > 0])
        temp$muerto_interc[i] <- length(tocm[tocm > 0])
        temp$tot[i] <- length(total[total > 0])
      }
      #str(temp)
    
      sp_touch[s, sp, "Live"] <- sum(temp$vivo_interc %>% as.numeric)
      sp_touch[s, sp, "Dead"] <- sum(temp$muerto_interc %>% as.numeric)
      sp_touch[s, sp, "Live or dead"] <- sum(temp$tot)
    }
  }
}

# Normalize and get proportions by site

# number o voxels by site (points * 16)
site_points <- aggregate(rep(1, nrow(es)) ~ site, es, length)
dtemp <- data.frame(site = dimnames(sp_touch)$site)
dtemp <- left_join(dtemp, site_points, by = "site")

sp_touch_norm <- sp_touch
sp_touch_prop <- sp_touch
for(i in 1:3) {
  sp_touch_norm[, , i] <- sp_touch[, , i] / rowSums(sp_touch[, , i])
  sp_touch_prop[, , i] <- sp_touch[, , i] / dtemp[, 2]
}

# check
for(i in 1:3) {
  dddd <- rowSums(sp_touch_norm[, , i])
  print(unname(dddd))
}

# Plot relative abundances by site for every fuel type.
sp_norm <- as.data.frame.table(sp_touch_norm)
sp_norm <- separate_site(d = sp_norm)
names(sp_norm)[names(sp_norm) == "Freq"] <- "p"
head(sp_norm)

# get full species name
sp_data <- es2[!duplicated(es2[, c("sp_code", "sp_name")]), c("sp_code", "sp_name")]
names(sp_data) <- c("sp", "Species")
sp_norm <- left_join(sp_norm, sp_data, by = "sp")

ggplot(sp_norm, aes(x = Community, colour = sp, fill = sp, y = p)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_grid(cols = vars(transect), rows = vars(fuel_type))

# Live fuel
ggplot(sp_norm[sp_norm$fuel_type == "Live", ], 
       aes(x = Community, colour = Species, fill = Species, y = p)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(vars(transect), nrow = 3)

# ggsave("species relative abundance plot.png", width = 20, height = 30,
#        units = "cm")

# Choose the k most abundant species by site for live fuel.
k <- 5

sp_most <- sp_touch[, , "Live"]
for(i in 1:nsites) {
  #i = 1
  temp <- sp_touch_norm[i, , "Live"]
  ord <- order(temp, decreasing = TRUE)
  sp_most[i, ] <- species[ord]
}

sp_norm$inside <- FALSE
for(i in sites) {
  #i = 1
  sp_focal <- sp_most[i, 1:k] %>% unname
  filter <- sp_norm$site == i & 
            sp_norm$sp %in% sp_focal &
            sp_norm$fuel_type == "Live"
  sp_norm$inside[filter] <- TRUE
}

# Plot with less species
theme_set(theme_bw())
ggplot(sp_norm[sp_norm$inside, ], 
       aes(x = Community, fill = Species, y = p)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(vars(transect), nrow = 2) + 
  ylab("Relative abundance") + 
  theme(legend.position = c(0.75, 0.25))

# ggsave("species relative abundance plot.png", width = 20, height = 20,
#        units = "cm")

# Export data with normalized abundance by site to compute NMDS
# saveRDS(sp_touch_norm, "vegetation_structure/sp_abundance_norm_by_site_Robject.R")
# view(sp_touch_norm)


# NMDS (performed in the other plot script)

nmds <- metaMDS(sp_touch_norm[, , "Live"], k = 2)
plot(nmds)
nmds$points
dplot <- data.frame(site = rownames(nmds$points))
dplot <- cbind(dplot, nmds$points)
dplot <- separate_site(dplot)

ggplot(data = dplot, 
       aes(x = MDS1, y = MDS2, colour = Community, shape = transect)) +
  geom_point()


# Litter depth and volume and canopy cover --------------------------------

# get litter
lit_raw <- es[complete.cases(es[, c("prof_mantillo_cm", "cob_mantillo_percent")]), ]
variables <- c("prof_mantillo_cm", "cob_mantillo_percent")

lit <- pivot_longer(lit_raw, which(names(lit_raw) %in% variables), names_to = "var",
                    values_to = "response")

# get canopy cover 

cover <- read.csv("vegetation_structure/data_canopy_cover.csv", header = TRUE, sep = ",")

cover$Community <- factor(cover$Community, levels = c("Mature forest", "Young forest", "Closed shrubland", "Open shrubland"))
cover$transect <- factor(cover$transect, levels = c("T1", "T2", "T3"))
names(cover)[3] <- "response"
cover$var <- "canopy_cover"

# merge litter and canopy cover
lit %>% head
ddd <- rbind(lit[, c("transect", "Community", "var", "response")],
             cover[, c("transect", "Community", "var", "response")])

unique(ddd$var)
ddd$var <- factor(ddd$var, levels = unique(ddd$var),
                  labels = c(
                    "Litter depth (cm)", 
                    "Litter cover (%)",
                    "Canopy cover (%)"
                  ))

# Boxplot without extreme points
# https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap
calc_stat <- function(x) {
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), method = 8)
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(stats)
}

# plot
ggplot(ddd, aes(x = transect, colour = Community, y = response)) + 
  stat_summary(fun.data = calc_stat, geom = "boxplot", 
               position = position_dodge2(preserve = "single", padding = 0.05,
                                          width = 0.05)) + 
  facet_wrap(vars(var), nrow = 2, scales = "free_y", strip.position = "left") + 
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = c(0.75, 0.25)
  ) + 
  ylab(NULL) +
  xlab("Transect")


# Species abundance table -------------------------------------------------

# sp_touch_norm[, , "Live"]
# dim(sp_touch_norm[, , "Live"])

spt <- as.data.frame.table(sp_touch_norm[, , "Live"])
names(spt)[2] <- "sp_code"

spcodes <- es[!duplicated(es[, c("sp_code", "sp_name")]), c("sp_code", "sp_name")]
spcodes <- spcodes[complete.cases(spcodes), ]

spt <- separate_site(spt)
spt <- left_join(spt, 
                 spcodes, 
                 by = "sp_code")

names(spt) <- c("site", "Community", "Transect", "sp_code", "Relative abundance", "Species")
spt <- spt[order(spt$Transect, spt$Community), ]
spt$site2 <- paste(spt$Transect, spt$Community, sep = " / ")

sptable <- pivot_wider(spt[, c("Relative abundance", "Species", "site2")], 
                       names_from = "site2", values_from = "Relative abundance")
View(sptable)
sptable[1, -1] %>% as.numeric()
spsums <- sapply(1:nrow(sptable), function(i) sum(sptable[i, -1]))
sptable <- sptable[order(spsums, decreasing = TRUE), ]

# View(sptable)
write.csv(sptable, "vegetation_structure/data_live_fuel_relative_abundance.csv")



# Vegetation structure plots ----------------------------------------------

# Graphic parameters

# Considering A4 format is 21 cm wide
# Units in cm
fig_width <- 21 - 4
fig_height <- 12

a4_length <- 29.7 
fig_height_full <- a4_length - 7

# 2 cm lateral margins in paper
# Font size 11 for almost every text, with axis titles size 12.

# viridis colors, from
# https://www.thinkingondata.com/something-about-viridis-library/ 
colorcines <- c("#000004FF", "#711A6EFF", "#C43C4EFF", "#FCA007FF")
colorcines_diff <- c("#711A6EFF", "#C43C4EFF", "#FCA007FF")

# custom ggplot theme

# from https://rpubs.com/mclaire19/ggplot2-custom-themes

theme_mine <- function() { 
  font <- "Arial"   #assign font family up front
  marg <- 2 # figure margin in mm
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      #panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      #axis.ticks = element_blank(),          #strip axis ticks
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 16,                #set font size
        #face = 'bold',            #bold typeface
        hjust = -0.1,                #left align
        vjust = 1),               
      
      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        size = 12),               
      
      # para separar el eje y de los nros
      axis.title.y = element_text(             
        margin = margin(t = 0, r = 2, b = 0, l = 0, "mm"),
        angle = 90),
      
      axis.text = element_text(              #axis text
        family = font,            #axis family
        size = 9),                #font size
      
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9, family = font),
      
      strip.text = element_text(size = 12, family = font, color = "white"),
      strip.text.x = element_text(margin = margin(1.2,0,1.2,0, "mm")), # tamaño de la cajita
      strip.text.y = element_text(margin = margin(0,1.2,0,1.2, "mm")),
      strip.background = element_rect(fill = "gray10", color = "gray10"),
      
      plot.margin = unit(c(marg, marg, marg, marg), "mm")
    )
}

theme_set(theme_mine())


#### Fig s2: fuel proportion.

w <- 0.3
ggplot(data = fuel_prop, 
       aes(x = transect, y = p_obs, ymin = p_lower, ymax = p_upper, 
           colour = Community)) +
  geom_point(position = position_dodge2(width = w),
             size = 3) +
  geom_linerange(position = position_dodge2(width = w),
                 size = 0.7) +
  facet_wrap(vars(fuel_type), nrow = 2, scales = "free_y") + 
  scale_color_manual(values = colorcines) +
  #scale_fill_manual(values = colorcines) + 
  theme(legend.position = c(0.75, 0.25),
        panel.grid.minor = element_blank()) +
  xlab("Transect") +
  ylab("Fuel presence proportion") 
ggsave(paste("plots/Fig S02_fuel proportion.png", sep = ""),
       width = fig_width, height = fig_height * 1.15, units = "cm")

#### Fig S3: vertical profile

ggplot(data = nd, aes(x = h, y = p_mle, ymin = p_lower, ymax = p_upper, 
                      colour = Community, fill = Community)) +
  geom_line() + 
  geom_ribbon(color = NA, alpha = 0.2) +
  geom_point(data = int_agg, 
             mapping = aes(x = h, y = prop, colour = Community),
             inherit.aes = FALSE,
             shape = 1
  ) +
  facet_grid(rows = vars(fuel_type), cols = vars(transect)) + 
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    legend.text = element_text(margin = margin(t = 2, b = 2, l = -2, r = 2, unit = "mm"))
  ) +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75)) +
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) + 
  xlab("Height (m)") +
  ylab("Fuel presence probability or proportion") +
  rotate_y_facet_text(angle = 270, align = 0.5)
ggsave(paste("plots/Fig S03_fuel vertical profile.png", sep = ""),
       width = fig_width, height = fig_height * 1.5, units = "cm")

#### Litter depth and cover and canopy cover

ggplot(ddd, aes(x = transect, colour = Community, y = response, fill = Community)) + 
  stat_summary(fun.data = calc_stat, geom = "boxplot", 
               position = position_dodge2(preserve = "single", padding = 0.05,
                                          width = 0.05),
               alpha = 0.25) + 
  facet_wrap(vars(var), nrow = 2, scales = "free_y", strip.position = "left") + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 11),
    strip.placement = "outside",
    legend.position = c(0.75, 0.25)
  ) + 
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) + 
  ylab(NULL) +
  xlab("Transect")
ggsave(paste("plots/Fig S04_litter and canopy.png", sep = ""),
       width = fig_width, height = fig_height * 1.2, units = "cm")


#### Species composition plot

# Choose the k most abundant species by site for live fuel.
k <- 6

sp_most <- sp_touch[, , "Live"]
for(i in 1:nsites) {
  #i = 1
  temp <- sp_touch_norm[i, , "Live"]
  ord <- order(temp, decreasing = TRUE)
  sp_most[i, ] <- species[ord]
}

sp_norm$inside <- FALSE
for(i in sites) {
  #i = 1
  sp_focal <- sp_most[i, 1:k] %>% unname
  filter <- sp_norm$site == i & 
    sp_norm$sp %in% sp_focal &
    sp_norm$fuel_type == "Live"
  sp_norm$inside[filter] <- TRUE
}

# order species
sp_norm$Species <- factor(sp_norm$Species, levels = sptable$Species[])

sp_plot <- sp_norm[sp_norm$inside, ]
sps <- unique(sp_plot$Species)
# set.seed(2342)
# sps_ids <- sample(1:length(sps))
# sp_plot$Species <- factor(sp_plot$Species, levels = sps[sps_ids])

nsp <- length(unique(sp_plot$Species)) # substitute nrow  for length if you have rows instead of a list
col_pal <- viridis(nsp, option = "H")
# trying seeds
# seed <- abs(rnorm(1, 1000, 1000)) %>% as.integer 
# good ones: 934, 677, 126, 5, 2003
#set.seed(seed)
set.seed(934)
idsp <- sample(1:nsp)
col_pal2 <- col_pal[idsp]
# 15  5  8 10  6 14 13  4  2 16 11 12  1 17  7  9  3

ggplot(sp_plot, 
       aes(x = Community, fill = Species, y = p)) +
  geom_bar(position = "stack", stat = "identity", color = "black", size = 0.4) +
  facet_wrap(vars(transect), nrow = 3) + 
  ylab("Relative abundance") + 
  theme(legend.position = "right") +
  scale_fill_manual(values = col_pal2)

ggsave(paste("plots/Fig S05_species composition.png", sep = ""),
       width = fig_width, height = fig_height_full, units = "cm")
