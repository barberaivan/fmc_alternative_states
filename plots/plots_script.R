# Plots 

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


# Packages ----------------------------------------------------------------

library(tidyverse)
library(dplyr) # revalue
library(viridis) # colors
library(desiderata) # rotate_x_facet_text, to left-aling facet text
library(gridExtra) # to replace ggarrange
library(ggh4x) # facet_nested
library(vegan) # NMDS


# ggpubr::ggarrange has problems with fonts. gridExtra worked.
# gridExtra tutorial
# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html

# custom ggplot theme -----------------------------------------------------

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

# function to share legend --------------------------------------------------
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

g_ylab<-function(a.gplot){
  a.gplot <- p1
  a.gplot
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "ylab")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

g_ylab(p1)


# Data and predictions ------------------------------------------------------

# Load FMC raw data
fdata <- readRDS("data_fmc_dataframe.R")

# Load models predictions
lll <- list.files("fmc_models/")

listitas <- lll[startsWith(lll, "model_predictions")]

N_fuels <- 7
fuel_names <- c("1 h fuel sticks", "10 h fuel sticks", 
                "Cania", "Laura", "Coihue",
                "Litter", "Live fuel mixture")

flist <- lapply(1:N_fuels, function(f) {
  readRDS(paste("fmc_models/", listitas[f], sep = ""))
}) 

names(flist) <- fuel_names
#str(flist[[7]])

# Table with probabilities of difference <0 (for fig 2)
diff_table <- do.call("rbind", lapply(1:N_fuels, function(f) {
  d <- as.data.frame(flist[[f]]$diff_table)
}))

# Table with FMC as a function of VPD - sites estimates (Fig 3)
vpd_sites <- do.call("rbind", lapply(1:N_fuels, function(f) {
  #f = 3
  d <- as.data.frame(flist[[f]]$site_vpd_prediction) 
}))

# Table with FMC as a function of VPD - prediction (Fig 3)
vpd_pred <- do.call("rbind", lapply(1:N_fuels, function(f) {
  d <- as.data.frame(flist[[f]]$vpd_prediction)
}))

# Functions to rename variables

renamer <- function(data) {
  #data <- fdata
  data_new <- data
  n <- names(data)
  
  # to put italics in facets
  # https://stackoverflow.com/questions/16490331/combining-new-lines-and-italics-in-facet-labels-with-ggplot2
  # fuel type
  if("fuel_type" %in% n) {
    data_new$fuel_type <- factor(as.character(data$fuel_type),
                                 levels = c("1 h fuel sticks",
                                            "10 h fuel sticks",
                                            "Litter",
                                            "Laura",
                                            "Coihue",
                                            "Cania",
                                            "Live fuel mixture"),
                                labels = c(
                                   "textstyle(a.~1~h~fuel~sticks)",
                                   "textstyle(b.~10~h~fuel~sticks)",
                                   "textstyle(c.~Litter)",
                                   "textstyle(d.)~italic('S. patagonicus')",
                                   "textstyle(e.)~italic('N. dombeyi')",
                                   "textstyle(f.)~italic('C. culeou')",
                                   "textstyle(g.~Live~fuel~mixture)"))
    
    # separating live and dead for letters (a, b, c)
    # (used in figures 3 and 4)
    data_new$fuel_type_sep <- factor(as.character(data$fuel_type),
                                 levels = c("1 h fuel sticks",
                                            "10 h fuel sticks",
                                            "Litter",
                                            "Laura",
                                            "Coihue",
                                            "Cania",
                                            "Live fuel mixture"),
                                 labels = c(
                                   "textstyle(a.~1~h~fuel~sticks)",
                                   "textstyle(b.~10~h~fuel~sticks)",
                                   "textstyle(c.~Litter)",
                                   "textstyle(a.)~italic('S. patagonicus')",
                                   "textstyle(b.)~italic('N. dombeyi')",
                                   "textstyle(c.)~italic('C. culeou')",
                                   "textstyle(d.~Live~fuel~mixture)"))
  }
 
  # transect
  if("transect" %in% n) {
    data_new$transect <- factor(as.character(data_new$transect),
                                levels = c("t1", "t2", "t3"),
                                labels = c("T1", "T2", "T3"))
  }
  
  # Community
  if("Community" %in% n) {
    data_new$Community <- factor(as.character(data_new$Community),
                                 # order as a function of cover, not spatial order
                                 levels =  c("U forest", "B forest +reg", 
                                             "B shrubland", "B forest -reg"),
                                 labels = c("Mature forest", "Young forest",
                                            "Closed shrubland", "Open shrubland")
                                 )
  }
  
  return(data_new)
}


# change names

fdata <- renamer(fdata)
diff_table <- renamer(diff_table)
vpd_pred <- renamer(vpd_pred)
vpd_sites <- renamer(vpd_sites)


# Function to make boxplot without extreme points
# https://stackoverflow.com/questions/25124895/no-outliers-in-ggplot-boxplot-with-facet-wrap
calc_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  #stats <- quantile(x, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), method = 8)
  stats <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), method = 8)
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(stats)
}


# Extra: raw data for figure 3, asked by Thomas. Add fuel_type raw column.
# fdata_export <- fdata
# fdata_export$fuel_type <- factor(fdata_export$fuel_type, 
#                                 levels = levels(fdata$fuel_type),
#                            labels = c(
#                              c("1 h fuel sticks",
#                                "10 h fuel sticks",
#                                "Laura",
#                                "Coihue",
#                                "Cania",
#                                "Litter",
#                                "Live fuel mixture")
#                            ))
# head(fdata_export)
# write.csv(fdata_export[, c("fuel_type", "transect", "Community", "point_id",
#                            "date", "moisture")], 
#           "data_fmc_table2.csv")



# Fig 2: VPD ----------------------------------------------------------------


# order VPD data  ---------------------------------------------------------

clim <- readRDS("microclimate/data_microclimate_daily_dataframe.R")
clim <- clim[clim$time <= 85, 
             c("fDate", "time", "Community", "site", "transect", "hum_min", "temp_max", "vpd_max",
               "hum_mean", "temp_mean", "vpd_mean")]

clim2 <- clim
posits <- which(names(clim2) %in% c("hum_min", "temp_max", "vpd_max"))
names(clim2)[posits] <- c("Minimum relative humidity (%)", 
                          "Maximum temperature (°C)", 
                          "Maximum vapour pressure deficit (kPa)")
clim_long <- pivot_longer(clim2, posits, names_to = "variable", values_to = "y")

# ggplot(clim_long, aes(x = time, y = y, colour = Community)) +
#   geom_line() + 
#   facet_grid(rows = vars(variable), cols = vars(transect), 
#              scales = "free_y")
# es demasiado grande

clim$Community <- factor(as.character(clim$Community),
                         levels = c("U forest", "B forest +reg", "B shrubland", "B forest -reg"),
                         labels = levels(fdata$Community))

clim$transect <- factor(as.character(clim$transect), 
                        levels = c("t1", "t2", "t3"),
                        labels = c("T1", "T2", "T3"))


# Extra: raw data for figure 2, asked by Thomas. Add fuel_type raw column.
#write.csv(clim, "microclimate/data_microclimate_daily2.csv")

# add column for ts plot
clim$data_type <- "Time series"

colorcines_air <- c("#000004FF", "#8c2369", "#C43C4EFF", "#FCA007FF")
# colorcines2 <- c("fcffa4", "#ed6925", "#781c6d", "#000004")
# colorcines2 <- c("#000004FF", "#ed6925", "#781c6d", "#000004")
# colorcines_air <- c("#000004FF", "#ed6925", "#C43C4EFF", "#FCA007FF")


# get rain data
raind <- read.csv("microclimate/data_rain.csv")
raind$transect <- "T3"
raind$data_type <- "Time series"
# scale rain
plot(rain_mm ~ time, raind); abline(v = -40)
raind$rain_scaled <- raind$rain_mm / max(clim$vpd_max)
raind$Community <- "Rain"

ggplot(raind, aes(x = time, y = rain_scaled)) +
  geom_line()


# VPD ts and means --------------------------------------------------------
theme_set(theme_mine())
# time series
#colorcines_air_rain <- c(colorcines_air, "gray40")


(vpd_ts <- ggplot(clim, aes(x = time, y = vpd_max, colour = Community)) +
  
  # add rain 
  # geom_bar(data = raind[raind$time > -40 & raind$time <= 85, ],
  #          mapping = aes(x = time, y = rain_scaled),
  #          stat = 'identity', color = NA, fill = "gray50") +
  # 
  #  scale_y_continuous(
  #    
  #    # Add a second axis and specify its features
  #    sec.axis = sec_axis(~ . * max(clim$vpd_max), name = "Pp (mm)")
  #                        #guide = guide_axis(position = "top"))
  #  )  +
  # geom_line(data = raind[raind$time > -40, ], 
  #            mapping = aes(x = time, y = rain_scaled),
  #            color = "gray50") + 
  #  
  geom_line(size = 0.6) + #0.6
   
  facet_grid(cols = vars(data_type),
             rows = vars(transect)) +
  scale_color_manual(values = colorcines_air) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(margin = margin(1,0,1,0, "mm"), angle = 0), # tamaño de la cajita
        # strip.text.y = element_text(margin = margin(0,1,0,1, "mm"), angle = 0),
        panel.spacing = unit(2, "mm"),
        legend.box.margin = margin(-4, 0, 0, 0, "mm"),  # y esto igual
        legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")),  # mejora el asunto de separación de símbolos
        strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        #panel.border = element_rect(fill = "gray93", color = "gray93")  
  ) + 
  ylab("max VPD (kPa)") + 
  xlab("Days since summer solstice")
)

# colorcines_air <- c("#000004FF", "#8c2369", "#C43C4EFF", "#FCA007FF")
colorcines_air <- viridis(4, option = "B", end = 0.8)

(vpd_ts2 <- ggplot(clim, aes(x = time, y = vpd_max, colour = Community)) +
    
    geom_line(size = 0.6) + #0.6
    
    facet_grid(cols = vars(data_type),
               rows = vars(transect)) +
    scale_color_manual(values = colorcines_air) +
    theme(legend.position = "bottom",
          strip.text.x = element_text(margin = margin(1,0,1,0, "mm"), angle = 0), # tamaño de la cajita
          # strip.text.y = element_text(margin = margin(0,1,0,1, "mm"), angle = 0),
          panel.spacing = unit(2, "mm"),
          legend.box.margin = margin(-4, 0, 0, 0, "mm"),  # y esto igual
          legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")),  # mejora el asunto de separación de símbolos
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          #panel.border = element_rect(fill = "gray93", color = "gray93")  
          
          # remove x axis
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(color = "white")
          
    ) + 
    ylab("max VPD (kPa)") + 
    xlab("Days since summer solstice") +
    scale_x_continuous(breaks = c(-30, 0, 30, 60), limits = c(-45, max(clim$time)))
)

# precipitation separate plot
(rain_plot <- ggplot(raind[raind$time > -40 & raind$time <= 85, ], 
       aes(x = time, y = rain_mm/2)) + 
  geom_bar(stat = 'identity', color = "gray50", fill = "gray50") +
  xlab("Days since summer solstice") + 
  ylab("Precipitation (mm)") + 
  #xlim(-45, max(clim$time)) +
  scale_x_continuous(breaks = c(-30, 0, 30, 60), limits = c(-45, max(clim$time)))
)
  


# means 

air_pred <- readRDS("microclimate/predictions_vpd.R")

# means by site:
airmeans <- renamer(air_pred$by_site)

diff_air <- air_pred$diff[air_pred$diff$transect != "average", ]
diff_air <- renamer(diff_air)

# get 0.95 percentiles to place the probability above
# percents <- aggregate(vpd_max ~ Community + transect, 
#                       clim, quantile, probs = 0.95, method = 8)
#percents$vpd_max <- percents$vpd_max #* 1.05
diff_air <- left_join(diff_air, airmeans[, c("Community", "transect", "vpd_upper")], 
                      by = c("Community", "transect"))
diff_air$p <- as.character(round(diff_air$p, digits = 2))

# Add blank p for Mature forest
diff_air2 <- diff_air[c(1, 4, 7), ]
diff_air2$Community <- "Mature forest"
diff_air2$p <- "   "

diff_air3 <- rbind(diff_air, diff_air2)

# add column for facet_grid
diff_air3$data_type <- "Season average"
airmeans$data_type <- "Season average"

(vpd_means <- 
ggplot(airmeans, 
       aes(x = Community, y = vpd_mean, ymin = vpd_lower, ymax = vpd_upper, 
           color = Community, #fill = Community,
           shape = Community)) + 
  
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) + 
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  
  facet_grid(cols = vars(data_type), 
             rows = vars(transect)) +
  
  # estimated means:
  geom_point(size = 4) +
  geom_linerange(size = 2, alpha = 0.4) +
  
  # probabilities
  geom_text(diff_air3, 
            mapping = aes(label = p, y = vpd_upper * 1.05, x = Community,
                          colour = Community),
            vjust = -0.1, color = "gray20", size = 3,
            inherit.aes = F) +

  # observed means:
  # ggnewscale::new_scale_fill() + 
  # #scale_shape_manual(values = c(0, 1, 2, 5)) +
  # scale_fill_manual(values = colorcines) +
  # 
  # geom_point(airmeans, mapping = aes(y = vpd_max), 
  #            size = 2.8, stroke = 0.5) +
  
  theme(legend.position = "bottom",
        axis.text.x = element_text(color = "white"),
        axis.ticks.x = element_line(color = "white"),
        panel.spacing = unit(2, "mm"),
        strip.text.x = element_text(margin = margin(1,0,1,0, "mm"), angle = 0), # tamaño de la cajita
        legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")),  # mejora el asunto de separación de símbolos
        ) +#,
  #theme(legend.position = "right") +#,
  #strip.text.left = element_text(angle = 0, vjust = 1)) + 
  ylab(NULL) + 
  xlab("Community") + 
  ylim(1, 6.5)
)


(vpd_means2 <- 
    ggplot(airmeans, 
           aes(x = Community, y = vpd_mean, ymin = vpd_lower, ymax = vpd_upper, 
               color = Community, #fill = Community,
               shape = Community)) + 
    
    scale_color_manual(values = colorcines) +
    scale_fill_manual(values = colorcines) + 
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    
    facet_grid(cols = vars(data_type), 
               rows = vars(transect)) +
    
    # estimated means:
    geom_point(size = 4) +
    geom_linerange(size = 2, alpha = 0.4) +
    
    # probabilities
    geom_text(diff_air3, 
              mapping = aes(label = p, y = vpd_upper * 1.05, x = Community,
                            colour = Community),
              vjust = -0.1, color = "gray20", size = 3,
              inherit.aes = F) +
 
    theme(legend.position = "right",
          # axis.text.x = element_text(color = "white"),
          # axis.ticks.x = element_line(color = "white"),
          panel.spacing = unit(2, "mm"),
          strip.text.x = element_text(margin = margin(1,0,1,0, "mm"), angle = 0), # tamaño de la cajita
          legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")),  # mejora el asunto de separación de símbolos

          # remove x axis
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
           
    ) +#,
    #theme(legend.position = "right") +#,
    #strip.text.left = element_text(angle = 0, vjust = 1)) + 
    ylab(NULL) + 
    xlab("Community") + 
    ylim(1, 6.5)
)


# VPD ts and means 

vpd_legend <- g_legend(vpd_means)

vpd_all <- grid.arrange(arrangeGrob(vpd_ts + theme(legend.position = "none"),
                                    vpd_means + theme(legend.position = "none"),
                                    nrow = 1, widths = c(2.6, 1)),
                        vpd_legend, nrow = 2, heights = c(10, 1))

# ggsave("plots/Fig 02_vpd_ts_means.png", plot = vpd_all, 
#        width = fig_width, height = fig_height, # fig_height + 7,
#        units = "cm") # fig_height_full leaves 4 cm margins and 3 cm text


# plot including rain

# join vpd plots
vpd_ts_mean <- grid.arrange(vpd_ts2 + theme(legend.position = "none"),
                      vpd_means2 + theme(legend.position = "none"),
                      nrow = 1, widths = c(2.6, 1))

# get legend and merge it with a white rectangle at the bottom and at the left
vpd_leg_tmp1 <- g_legend(vpd_means2)# + 
                #theme(plot.margin = margin(t = 0, l = 5, b = 5, r = 2, unit = 'mm')))

blankPlot <- ggplot() + geom_blank(aes(1,1)) + cowplot::theme_nothing()
vpd_leg_tmp2 <- grid.arrange(vpd_leg_tmp1, blankPlot,
                             nrow = 2, heights = c(4, 1.3))
vpd_leg <- grid.arrange(blankPlot, vpd_leg_tmp2,
                        nrow = 1, widths = c(0.5, 4))
# merge legend with pp (the plot.margin aligns the vpd ts with rain ts)
rain_leg <- grid.arrange(rain_plot + theme(plot.margin = margin(l = 4.7, r = 2, unit = "mm")), 
                         vpd_leg,
                         nrow = 1, widths = c(2.6, 1))

# merge all
vpd_rain <- grid.arrange(vpd_ts_mean, rain_leg, nrow = 2, heights = c(2.8, 1))

ggsave("plots/Fig 02_vpd_ts_means_rain.png", plot = vpd_rain, 
       width = fig_width, height = fig_height + 4, # fig_height + 7,
       units = "cm") # fig_height_full leaves 4 cm margins and 3 cm text

# Fig 3 Dead fuel ts and means --------------------------------------------

# Decrease font size 
theme_mine_small <- function() {
  theme_mine() %+replace%
    theme(
      # decrease text sizes
      strip.text = element_text(size = 11, color = "white"), # 12
      strip.text.y = element_text(size = 8.5, color = "white"), # 12
      axis.title = element_text(size = 11), # 12              
      axis.text = element_text(size = 8)    # 9
    )
}
theme_set(theme_mine_small())

# prepare data

fdata_ts <- aggregate(moisture ~ fuel_type + fuel_type_sep + Community + transect + time, fdata, mean)

levels_live <- levels(fdata$fuel_type)[4:7]
levels_dead <- levels(fdata$fuel_type)[1:3]
live_filter <- fdata_ts$fuel_type %in% levels_live
dead_filter <- fdata_ts$fuel_type %in% levels_dead

fdata_ts$data_type <- "textstyle(Time~series)"

# manage the probabilities
diff_table_temp <- diff_table[, c("fuel_type", "fuel_type_sep", "transect", "Community", "p")]
diff_extra <- diff_table_temp[!duplicated(diff_table_temp[, c("fuel_type", "transect")]), ]
diff_extra$Community <- factor("Mature forest", levels = levels(diff_table$Community))
diff_extra$p <- 0

diff_table_temp2 <- rbind(diff_table_temp, diff_extra)
diff_table_p <- left_join(diff_table_temp2, vpd_sites[, c("fuel_type", "Community", "transect", "fmc_upper")], 
                          by = c("fuel_type", "Community", "transect"))
diff_table_p$p <- as.character(round(diff_table_p$p, digits = 2))
diff_table_p$p[diff_table_p$Community == "Mature forest"] <- "   "

# a few column names

fdata_ts2 <- fdata_ts
fdata_ts2$data_type <- "textstyle(Dates~average)"
diff_table_p$data_type <- "textstyle(Dates~average)"
vpd_sites$data_type <- "textstyle(Dates~average)"

www <- 0.7
lll <- 2.5
psize <- 3


# Dead fuels TS -----------------------------------------------------------

(dead_ts <- ggplot(fdata_ts[dead_filter, ], 
                   aes(x = time, y = moisture, colour = Community,
                       shape = Community, fill = Community)) + 
    geom_point(size = 2) + #, key_glyph = "rect") + 
    geom_line(alpha = 0.5) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    #scale_shape_manual(values = c(0, 1, 2, 5)) +
    scale_color_manual(values = colorcines) +
    scale_fill_manual(values = colorcines) +
    
    facet_nested(rows = vars(fuel_type_sep), 
                 cols = vars(data_type, transect),
                 scales = "free_y", labeller = label_parsed) +
    rotate_y_facet_text(angle = 270, align = 0.05) + 
    ylab("Moisture content (%)") + 
    xlab("Days since summer solstice") + 
    theme(legend.position = "bottom",
          legend.box.margin = margin(-4, 0, 0, 0, "mm"),  # y esto igual
          legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")), # mejora el asunto de separación de símbolos
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
          strip.background = element_rect(fill = "gray10", color = "gray80",
                                          size = 0.1),
          #legend.key.size = unit(3,"mm")
    ) 
)

# ?draw_key
# Dead fuel means ---------------------------------------------------------

www <- 0.7
lll <- 2.5
psize <- 3

# to have 2 shape scales:
# https://github.com/eliocamp/ggnewscale

(dead_means <- ggplot(vpd_sites[vpd_sites$fuel_type %in% levels_dead, ], 
                      mapping = aes(y = fmc_mean, x = Community, 
                                    color = Community, 
                                    shape = Community)) + 
    # transparent points to increase y range
    # geom_point(diff_table_p[!is.na(diff_table_p$fmc_upper) & 
    #                              diff_table_p$fuel_type %in% levels_live, , ], 
    #            mapping = aes(x = transect, y = posit * 1.05),
    #            alpha = 0) +
    
    scale_color_manual(values = colorcines) +
    scale_fill_manual(values = colorcines) + 
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    
    # probabilities:
    geom_text(diff_table_p[!is.na(diff_table_p$fmc_upper) &
                                 diff_table_p$fuel_type %in% levels_dead, ],
              mapping = aes(label = p, y = fmc_upper * 1.04, x = Community, color = Community),
              # position = position_dodge2(preserve = "single", padding = 0.1,
              #                            width = www),
              vjust = -0.3, size = 2.6) + # , color = "gray20" # size 2.6 before
    
    # estimated means:
    geom_point(size = psize,
               position = position_dodge2(preserve = "single", padding = 0.1,
                                          width = www)) +
    geom_linerange(vpd_sites[vpd_sites$fuel_type %in% levels_dead, ], 
                   mapping = aes(ymin = fmc_lower, ymax = fmc_upper),
                   position = position_dodge2(preserve = "single", padding = 0.1,
                                              width = www), lwd = lll, alpha = 0.5) +
    # observed means
    
    #scale_color_manual(values = colorcines) +
    ggnewscale::new_scale_fill() + 
    # #scale_shape_manual(values = c(0, 1, 2, 5)) +
    scale_fill_manual(values = colorcines) +
    
    geom_point(vpd_sites[vpd_sites$fuel_type %in% levels_dead, ],
               mapping = aes(y = moisture, fill = Community), 
               size = psize, stroke = 0.5,
               position = position_dodge2(preserve = "single", padding = 0.1,
                                          width = www)) +
    
    # match ts means :
    geom_point(data = fdata_ts2[fdata_ts2$fuel_type %in% levels_dead, ], 
               mapping = aes(x = Community, y = moisture, colour = Community),
               alpha = 0) +
    
    facet_nested(rows = vars(fuel_type_sep), 
                 cols = vars(data_type, transect), scales = "free_y",
                 labeller = label_parsed) +
    
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(color = "white"),
          axis.ticks.x = element_line(color = "white"),
          panel.spacing.x = unit(0.0, units = "mm"),
          panel.grid.major.x = element_blank(),
          strip.text.x = element_text(margin = margin(0.2,0,0.2,0, "mm")),
          strip.text.y = element_text(margin = margin(0.5,1,0.5,1, "mm"),
                                      size = 9),
          strip.background = element_rect(fill = "gray10", color = "gray80",
                                          size = 0.1)
    ) +
    # el strip para este plot tiene que tener menor margen que el strip
    # para la ts. Allá usé 0.5, acá uso 0.1, Eso genera tamaños muy similares
    # que alinea los paneles.
    
    ylab(NULL) + 
    xlab("Transect - community") + 
    rotate_y_facet_text(angle = 270, align = 0.05) +
    scale_x_discrete(labels = c("a", "b", "c", "d"))
  
)

# Dead TS and means -----------------------------------------------------

dead_legend <- g_legend(dead_ts)

p5 <- grid.arrange(arrangeGrob(dead_ts + theme(legend.position ="none"),
                               dead_means + theme(legend.position ="none"),
                               nrow = 1, widths = c(1.5, 1)),
                   dead_legend, nrow = 2, heights = c(10, 0.8))

ggsave("plots/Fig 03_dead_ts_means.png", plot = p5, 
       width = fig_width, height = fig_height + 0,#fig_height_full - 8, 
       units = "cm") # fig_height_full leaves 4 cm margins and 3 cm text


# Fig 4 Live fuel ts and means --------------------------------------------

# Live fuels TS -----------------------------------------------------------


(live_ts <- ggplot(fdata_ts[live_filter, ], 
                   aes(x = time, y = moisture, colour = Community,
                       shape = Community, fill = Community)) + 
    geom_point(size = 2) + 
    geom_line(alpha = 0.5) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    #scale_shape_manual(values = c(0, 1, 2, 5)) +
    scale_color_manual(values = colorcines) +
    scale_fill_manual(values = colorcines) +
    
    facet_nested(rows = vars(fuel_type_sep), 
                 cols = vars(data_type, transect),
                 scales = "free_y", labeller = label_parsed) +
    rotate_y_facet_text(angle = 270, align = 0.05) + 
    ylab("Moisture content (%)") + 
    xlab("Days since summer solstice") + 
    theme(legend.position = "bottom",
          legend.box.margin = margin(-4, 0, 0, 0, "mm"),  # y esto igual
          legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")), # mejora el asunto de separación de símbolos
          strip.background.y = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
          strip.background = element_rect(fill = "gray10", color = "gray80",
                                          size = 0.1)
    ) +
    scale_x_continuous(breaks = c(-30, 0, 30, 60))
)


# Live fuel means ---------------------------------------------------------

# move probability value for T2, live fuel mixture, closed shrubland
names(diff_table_p)
filt <- with(diff_table_p, 
             fuel_type_sep == levels(fuel_type_sep)[7] &
               transect == "T2" & 
               Community == "Closed shrubland" &
               !is.na(fmc_upper))
diff_table_p$fmc_upper[filt] <- diff_table_p$fmc_upper[filt] + 20

# to have 2 shape scales:
# https://github.com/eliocamp/ggnewscale

(live_means <- ggplot(vpd_sites[vpd_sites$fuel_type %in% levels_live, ], 
                      mapping = aes(y = fmc_mean, x = Community, 
                                    color = Community, 
                                    shape = Community)) + 
    # transparent points to increase y range
    # geom_point(diff_table_p[!is.na(diff_table_p$fmc_upper) & 
    #                              diff_table_p$fuel_type %in% levels_live, , ], 
    #            mapping = aes(x = transect, y = posit * 1.05),
    #            alpha = 0) +
    
    scale_color_manual(values = colorcines) +
    scale_fill_manual(values = colorcines) + 
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    
    # probabilities:
    geom_text(diff_table_p[!is.na(diff_table_p$fmc_upper) &
                             diff_table_p$fuel_type %in% levels_live, ],
              mapping = aes(label = p, y = fmc_upper * 1.04, x = Community, color = Community),
              # position = position_dodge2(preserve = "single", padding = 0.1,
              #                            width = www),
              vjust = -0.3, size = 2.6) + # , color = "gray20" # size 2.6 before
    
    # estimated means:
    geom_point(size = psize,
               position = position_dodge2(preserve = "single", padding = 0.1,
                                          width = www)) +
    geom_linerange(vpd_sites[vpd_sites$fuel_type %in% levels_live, ], 
                   mapping = aes(ymin = fmc_lower, ymax = fmc_upper),
                   position = position_dodge2(preserve = "single", padding = 0.1,
                                              width = www), lwd = lll, alpha = 0.5) +
    # observed means
    
    #scale_color_manual(values = colorcines) +
    ggnewscale::new_scale_fill() + 
    # #scale_shape_manual(values = c(0, 1, 2, 5)) +
    scale_fill_manual(values = colorcines) +
    
    geom_point(vpd_sites[vpd_sites$fuel_type %in% levels_live, ],
               mapping = aes(y = moisture, fill = Community), 
               size = psize, stroke = 0.5,
               position = position_dodge2(preserve = "single", padding = 0.1,
                                          width = www)) +
    
    # match ts means :
    geom_point(data = fdata_ts2[live_filter, ], 
               mapping = aes(x = Community, y = moisture, colour = Community),
               alpha = 0) +
    
    facet_nested(rows = vars(fuel_type_sep), 
                 cols = vars(data_type, transect), scales = "free_y",
                 labeller = label_parsed) +
    
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(color = "white"),
          axis.ticks.x = element_line(color = "white"),
          panel.spacing.x = unit(0.0, units = "mm"),
          panel.grid.major.x = element_blank(),
          strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "mm")),
          strip.text.y = element_text(margin = margin(0.5,1.5,0.5,1.5, "mm")),
          strip.background = element_rect(fill = "gray10", color = "gray80",
                                          size = 0.1)
    ) +
    # el strip para este plot tiene que tener menor margen que el strip
    # para la ts. Allá usé 0.5, acá uso 0.1, Eso genera tamaños muy similares
    # que alinea los paneles.
    
    ylab(NULL) + 
    xlab("Transect - community") + 
    rotate_y_facet_text(angle = 270, align = 0.05) + 
    scale_x_discrete(labels = c("a", "b", "c", "d"))
  
)


# Live TS and means -----------------------------------------------------

live_legend <- g_legend(live_ts)

p4 <- grid.arrange(arrangeGrob(live_ts + theme(legend.position = "none"),
                               live_means + theme(legend.position = "none"),
                               nrow = 1, widths = c(1.5, 1)),
                   live_legend, nrow = 2, heights = c(10, 1))

ggsave("plots/Fig 04_live_ts_means.png", plot = p4, 
       width = fig_width, height =  fig_height + 3, # fig_height + 7,
       units = "cm") # fig_height_full leaves 4 cm margins and 3 cm text



# Fig 3: mc ts and means for all fuel types together ----------------------

# Decrease font size 
theme_mine_small <- function() {
  theme_mine() %+replace%
    theme(
      # decrease text sizes
      strip.text = element_text(size = 11, color = "white"), # 12
      strip.text.y = element_text(size = 8, color = "white"), # 12
      axis.title = element_text(size = 11), # 12              
      axis.text = element_text(size = 8)    # 9
    )
}
theme_set(theme_mine_small())

# prepare data

fdata_ts <- aggregate(moisture ~ fuel_type + fuel_type_sep + Community + transect + time, fdata, mean)

levels_live <- levels(fdata$fuel_type)[4:7]
levels_dead <- levels(fdata$fuel_type)[1:3]
live_filter <- fdata_ts$fuel_type %in% levels_live
dead_filter <- fdata_ts$fuel_type %in% levels_dead

fdata_ts$data_type <- "textstyle(Time~series)"

# manage the probabilities
diff_table_temp <- diff_table[, c("fuel_type", "fuel_type_sep", "transect", "Community", "p")]
diff_extra <- diff_table_temp[!duplicated(diff_table_temp[, c("fuel_type", "transect")]), ]
diff_extra$Community <- factor("Mature forest", levels = levels(diff_table$Community))
diff_extra$p <- 0

diff_table_temp2 <- rbind(diff_table_temp, diff_extra)
diff_table_p <- left_join(diff_table_temp2, vpd_sites[, c("fuel_type", "Community", "transect", "fmc_upper")], 
                          by = c("fuel_type", "Community", "transect"))
diff_table_p$p <- as.character(round(diff_table_p$p, digits = 2))
diff_table_p$p[diff_table_p$Community == "Mature forest"] <- "   "

# a few column names

fdata_ts2 <- fdata_ts
fdata_ts2$data_type <- "textstyle(Dates~average)"
diff_table_p$data_type <- "textstyle(Dates~average)"
vpd_sites$data_type <- "textstyle(Dates~average)"

www <- 0.7
lll <- 2.5
psize <- 3


# ts ----------------------------------------------------------------------

(ld_ts <- ggplot(fdata_ts, 
                   aes(x = time, y = moisture, colour = Community,
                       shape = Community, fill = Community)) + 
   geom_point(size = 2) + #, key_glyph = "rect") + 
   geom_line(alpha = 0.5) +
   scale_shape_manual(values = c(21, 22, 23, 24)) +
   #scale_shape_manual(values = c(0, 1, 2, 5)) +
   scale_color_manual(values = colorcines) +
   scale_fill_manual(values = colorcines) +
   
   facet_nested(rows = vars(fuel_type), 
                cols = vars(data_type, transect),
                scales = "free_y", labeller = label_parsed) +
   rotate_y_facet_text(angle = 270, align = 0.05) + 
   ylab("Moisture content (%)") + 
   xlab("Days since summer solstice") + 
   theme(legend.position = "bottom",
         legend.box.margin = margin(-4, 0, 0, 0, "mm"),  # y esto igual
         legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")), # mejora el asunto de separación de símbolos
         strip.background.y = element_blank(),
         strip.text.y = element_blank(),
         strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")),
         strip.background = element_rect(fill = "gray10", color = "gray80",
                                         size = 0.1),
         #legend.key.size = unit(3,"mm")
   ) 
)


# means -------------------------------------------------------------------

www <- 0.7
lll <- 2.5
psize <- 3

# to have 2 shape scales:
# https://github.com/eliocamp/ggnewscale

(ld_means <- ggplot(vpd_sites, 
                      mapping = aes(y = fmc_mean, x = Community, 
                                    color = Community, 
                                    shape = Community)) + 
    # transparent points to increase y range
    # geom_point(diff_table_p[!is.na(diff_table_p$fmc_upper) & 
    #                              diff_table_p$fuel_type %in% levels_live, , ], 
    #            mapping = aes(x = transect, y = posit * 1.05),
    #            alpha = 0) +
    
    scale_color_manual(values = colorcines) +
    scale_fill_manual(values = colorcines) + 
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    
    # probabilities:
    geom_text(diff_table_p[!is.na(diff_table_p$fmc_upper), ],
              mapping = aes(label = p, y = fmc_upper * 1.04, x = Community, color = Community),
              # position = position_dodge2(preserve = "single", padding = 0.1,
              #                            width = www),
              vjust = -0.3, size = 2.6) + # , color = "gray20" # size 2.6 before
    
    # estimated means:
    geom_point(size = psize,
               position = position_dodge2(preserve = "single", padding = 0.1,
                                          width = www)) +
    geom_linerange(vpd_sites, 
                   mapping = aes(ymin = fmc_lower, ymax = fmc_upper),
                   position = position_dodge2(preserve = "single", padding = 0.1,
                                              width = www), lwd = lll, alpha = 0.5) +
    # observed means
    
    #scale_color_manual(values = colorcines) +
    ggnewscale::new_scale_fill() + 
    # #scale_shape_manual(values = c(0, 1, 2, 5)) +
    scale_fill_manual(values = colorcines) +
    
    geom_point(vpd_sites,
               mapping = aes(y = moisture, fill = Community), 
               size = psize, stroke = 0.5,
               position = position_dodge2(preserve = "single", padding = 0.1,
                                          width = www)) +
    
    # match ts means :
    geom_point(data = fdata_ts2, 
               mapping = aes(x = Community, y = moisture, colour = Community),
               alpha = 0) +
    
    facet_nested(rows = vars(fuel_type), 
                 cols = vars(data_type, transect), scales = "free_y",
                 labeller = label_parsed) +
    
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(color = "white"),
          axis.ticks.x = element_line(color = "white"),
          panel.spacing.x = unit(0.0, units = "mm"),
          panel.grid.major.x = element_blank(),
          strip.text.x = element_text(margin = margin(0.1,0,0.1,0, "mm")),
          strip.text.y = element_text(margin = margin(0.5,1,0.5,1, "mm"),
                                      size = 9),
          strip.background = element_rect(fill = "gray10", color = "gray80",
                                          size = 0.1)
    ) +
    # el strip para este plot tiene que tener menor margen que el strip
    # para la ts. Allá usé 0.5, acá uso 0.1, Eso genera tamaños muy similares
    # que alinea los paneles.
    
    ylab(NULL) + 
    xlab("Transect - community") + 
    rotate_y_facet_text(angle = 270, align = 0.05) +
    scale_x_discrete(labels = c("a", "b", "c", "d"))
  
)

# ts and means -----------------------------------------------------

ld_legend <- g_legend(ld_ts)

p5 <- grid.arrange(arrangeGrob(ld_ts + theme(legend.position ="none"),
                               ld_means + theme(legend.position ="none"),
                               nrow = 1, widths = c(1.5, 1)),
                   ld_legend, nrow = 2, heights = c(10, 0.45))

ggsave("plots/Fig 03_mc_ts_means.png", plot = p5, 
       width = fig_width, height = fig_height_full + 1.5, 
       units = "cm") # fig_height_full leaves 4 cm margins and 3 cm text



# Fig 5: FMC as a function of VPD -----------------------------------------

theme_set(theme_mine())

# change vpd_mean por vpd_max if observed vpd is used

# shorten vpd range for points in vpd_pred:
vpd_from <- min(vpd_sites$vpd_max) - 0.1
vpd_to <- max(vpd_sites$vpd_max) + 0.1
filter_vpd <- !(vpd_pred$vpd < vpd_from | vpd_pred$vpd > vpd_to)

# and shorten within fuels
vpd_pred$show <- TRUE
for(f in unique(vpd_pred$fuel_type)) {
  f = "textstyle(f.)~italic('C. culeou')"
  vpd_from <- min(vpd_sites$vpd_max[vpd_sites$fuel_type == f]) - 0.1
  vpd_to <- max(vpd_sites$vpd_max[vpd_sites$fuel_type == f]) + 0.1
  vpd_pred$show[vpd_pred$fuel_type == f & (vpd_pred$vpd < vpd_from | vpd_pred$vpd > vpd_to)] <- FALSE
  #vpd_pred[vpd_pred$fuel_type == f, ] %>% View
}

ggplot(vpd_sites,
       aes(x = vpd_max, #xmin = vpd_lower, xmax = vpd_upper,
           y = fmc_mean, ymin = fmc_lower, ymax = fmc_upper,
           color = Community, shape = transect)) +
  geom_line(data = vpd_pred[vpd_pred$show, ], linetype = "dotted",
            mapping = aes(x = vpd, y = fmc_mean),
            inherit.aes = FALSE) + 
  geom_ribbon(data = vpd_pred[vpd_pred$show, ], alpha = 0.1, 
              mapping = aes(x = vpd, y = fmc_mean, ymin = fmc_lower, 
                            ymax = fmc_upper),
              inherit.aes = FALSE) + 
  geom_point(size = 2.75) +
  # geom_point(mapping = aes(y = moisture), 
  #            size = 2.75, shape = 1) +
  geom_errorbar(alpha = 0.7, size = 0.9, width = 0) +
  #geom_errorbarh(alpha = 0.7, size = 0.9, height = 0) +
  facet_wrap(vars(fuel_type), ncol = 2, scales = "free_y",
             labeller = label_parsed) +
  rotate_x_facet_text(angle = 0, align = 0.05) +
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) + 
  ylab("Moisture content (%)") + 
  xlab("max VPD (kPa)") +
  theme(
    legend.box = "horizontal",    # una legend arriba de la otra
    #legend.spacing.y = unit(2, "cm"),
    legend.text = element_text(margin = margin(t = 2, b = 2, l = 0, r = 0, unit = "mm")), # mejora el asunto de separación de símbolosc("#FFFFFF", "#FFFFFF", "#FFFFFF")
    legend.box.margin = margin(-3, 0, 0, 0, "mm"),
    legend.position = c(0.75, 0.12),
    panel.spacing = unit(3, "mm"),
  ) +
  xlim(min(vpd_pred$vpd[filter_vpd]), max(vpd_pred$vpd[filter_vpd]))


ggsave("plots/Fig 05_fmc and vpd.png",
       width = fig_width, height = fig_height_full, units = "cm")


# Fig 6: Ordination --------------------------------------------------------

# structure data: normalized abundance by site to compute NMDS
sp_touch_norm <- readRDS("vegetation_structure/sp_abundance_norm_by_site_Robject.R")
# view(sp_touch_norm)

# NMDS 

nmds <- metaMDS(sp_touch_norm[, , "Live"], k = 2)

# check nmds fit
# goodness(nmds) %>% plot(ylim = c(0, max(goodness(nmds))))
# I think it's good and has little variation
# stressplot(nmds) # very nice # Linear r2 = 0.958
# round(nmds$stress, digits = 3) # 0.079 


#nmds$points

# scale to [0, 1],
extreme <- max(abs(c(nmds$points %>% range, nmds$species %>% range)))
# so correlation among axis and covariates is displayed well
nmds$points <- nmds$points / extreme
spp_coord <- (nmds$species / extreme) %>% as.data.frame
spp_coord$spp <- rownames(spp_coord)

# create df to plot
dplot <- data.frame(site = rownames(nmds$points))
dplot <- cbind(dplot, nmds$points)

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
dplot <- separate_site(dplot)

# merge with fuel moisture data
fuel_lfm <- flist[["Live fuel mixture"]]$site_vpd_prediction %>% renamer
dplot <- left_join(dplot, fuel_lfm[, c("Community", "transect", "fmc_mean", 
                                       "vpd_mean", "fmc_lower", "fmc_upper", "moisture")], 
                   by = c("Community", "transect"))


# correlations between nmds axes and fmc
cors <- cor(dplot[, c("fmc_mean", "MDS1", "MDS2")])
#cors <- cor(dplot[, c("moisture", "MDS1", "MDS2")])
cor_data <- data.frame(MDS1 = c(0, cors[2, 1]),
                       MDS2 = c(0, cors[3, 1]))

# get most abundant spp
sp_sums <- sp_touch_norm[, , "Live"] %>% colSums
sp_sums <- sp_sums[order(sp_sums, decreasing = TRUE)]
plot(sp_sums)
out_idx <- length(sp_sums) : (length(sp_sums) - 3)
in_idx <- 1:10
sp_abundant <- names(sp_sums)[-out_idx]
sp_abundant <- names(sp_sums)[in_idx]
sp_show <- spp_coord[spp_coord$spp %in% sp_abundant, ]
names(dplot)

# Plot 2
ggplot(data = dplot, 
       aes(x = MDS1, y = MDS2, colour = Community)) +
  #geom_point(mapping = aes(shape = transect), size = 3, stroke = 1.5) +
  geom_point(mapping = aes(size = moisture), stroke = 1.5) + 
  geom_line(data = cor_data, mapping = aes(x = MDS1, y = MDS2),
            color = "black",
            arrow = arrow(length = unit(0.30, "cm"), 
                          ends = "first", type = "open"),
            inherit.aes = FALSE) +
  geom_text(sp_show, 
            mapping = aes(x = MDS1, y = MDS2, label = spp),
            color = "gray55",
            inherit.aes = FALSE) +
  geom_text(mapping = aes(x = MDS1, y = MDS2, label = transect),
            color = "black", vjust = -1) +
  scale_color_manual(values = colorcines) +
  #scale_shape_manual(values = c(1, 2, 0)) +
  scale_size(breaks = c(120, 150, 180, 220),
              labels = as.character(c(120, 150, 180, 220))) +
  theme(
    legend.position = "right",
    legend.title = element_text(),
    legend.box = "vertical",    # una legend arriba de la otra
    legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")), # mejora el asunto de separación de símbolos
  ) + 
  guides(colour = guide_legend(title = "Community", title.hjust = 0),
         shape = guide_legend(title = "Transect", title.hjust = 0),
         size = guide_legend(title = "Moisture content (%)", title.hjust = 0)) +
  ylim(range(c(dplot$MDS2, sp_show$MDS2)) * 1.1) +
  xlim(range(c(dplot$MDS1, sp_show$MDS1)) * 1.1) +
  # labs(colour = "Community",
  #      shape = "Transect") + 
  coord_fixed() +  # para resaltar la correlación
  ylab("NMDS axis 2") + xlab("NMDS axis 1") 

ggsave("plots/Fig 06_ordination.png",
       width = fig_width, height = fig_height, units = "cm")

# table with full names for image caption:
veg <- read.csv("vegetation_structure/data_vegetation_structure.csv")
veg_translate <- veg[!duplicated(veg$sp_name), c("sp_name", "sp_code")]

sporder <- sp_touch_norm[, , "Live"] %>% colSums()
spp_order <- names(sporder)[order(sporder, decreasing = TRUE)]
dnames <- data.frame(sp_code = spp_order)

left_join(dnames, veg_translate, by = "sp_code")

#    sp_code                 sp_name
# 1     Ndom      Nothofagus dombeyi
# 2     Herb              Herbaceous
# 3       Sp     Schinus patagonicus
# 4       Cc         Chusquea culeou
# 5       Ms         Mutisia spinosa
# 6       Lh         Lomatia hirsuta
# 7       Fi       Fabiana imbricata
# 8      Ari   Aristotelia chilensis
# 9     Nant   Nothofagus antarctica
# 10      Mb         Maytenus boaria
# 11     Aus  Austrocedrus chilensis
# 12      Bm    Berberis microphylla
# 13      Dj          Diostea juncea
# 14      Ec    Embothrium coccineum
# 15      Dc        Discaria chacaye
# 16     Bac             Bacaris sp.
# 17      Rm      Ribes magellanicum
# 18    Gamu    Gaultheria mucronata
# 19      Mc    Maytenus chubutensis
# 20    Ades             Adesmia sp.
# 21      Bs Berberis serratodentata
# 22    Myob     Myoschilos oblongum
# 23      Mm    Maytenus magellanica
# 24      Rr         Rosa rubiginosa


# A stackoverflow response explaining why nMDS doesn't show 
# explained variability of each axis:
# https://stackoverflow.com/questions/49223740/cumulative-variance-explained-for-nmds-in-r

# Fig 7: Ignition probability ----------------------------------------------

p_ig <- rbind(
  flist[["Litter"]]$p_ig,
  flist[["Live fuel mixture"]]$p_ig
)
p_ig_rename <- renamer(p_ig)
p_ig_rename$fuel_type <- factor(as.character(p_ig_rename$fuel_type),
                                 levels = c("textstyle(c.~Litter)",
                                            "textstyle(g.~Live~fuel~mixture)"),
                                 labels = c("Litter",
                                            "Live fuel mixture"))

# It's not possible to use facet_grid with scale = free_x, so I have to make 2
# plots and use ggarrange.

alpha_bar <- 0.6

(p1 <- 
  ggplot(p_ig_rename[p_ig_rename$fuel_type == "Litter", ], 
         aes(x = fmc_mean, xmin = fmc_lower, xmax = fmc_upper,
             y = ig_mean, ymin = ig_lower, ymax = ig_upper,
             colour = Community, shape = Community)) +
  
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
    
  geom_point(size = 3, alpha = 1, stroke = 0.7) + 
  geom_linerange(alpha = alpha_bar) +
  geom_linerange(alpha = alpha_bar, orientation = "y") +
  geom_point(size = 3,
    mapping = aes(y = ig_obs, x = moisture, 
                           fill = Community)) + 
    
  facet_grid(rows = vars(fuel_type), cols = vars(transect)) +
  # facet_wrap(cols = vars(transect), rows = vars(fuel_type),
  #            labeller = label_parsed) +  
    ylab("Ignition probability") + 
    xlab(NULL) +
  theme(
    legend.box = "vertical",    # una legend arriba de la otra
    legend.spacing.y = unit(-2, "mm"),
    legend.text = element_text(margin = margin(t = 2, b = 2, l = -2, r = 2, unit = "mm")), # mejora el asunto de separación de símbolosc("#FFFFFF", "#FFFFFF", "#FFFFFF")
    legend.box.margin = margin(-2, 0, 0, 0, "mm"),
    legend.position = "none",
    panel.spacing = unit(3, "mm"),
    strip.text.y = element_text(angle = 270, hjust = 0.5, vjust = 0.5)
  )
)


(p2 <- 
  ggplot(p_ig_rename[p_ig_rename$fuel_type == "Live fuel mixture", ], 
         aes(x = fmc_mean, xmin = fmc_lower, xmax = fmc_upper,
             y = ig_mean, ymin = ig_lower, ymax = ig_upper,
             colour = Community, shape = Community)) +
  
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
  
  geom_point(size = 3, alpha = 1, stroke = 0.7) + 
  geom_linerange(alpha = alpha_bar) +
  geom_linerange(alpha = alpha_bar, orientation = "y") +
  geom_point(size = 3,
             mapping = aes(y = ig_obs, x = moisture, 
                           fill = Community)) + 
  
  facet_grid(rows = vars(fuel_type), cols = vars(transect)) +
  # facet_wrap(cols = vars(transect), rows = vars(fuel_type),
  #            labeller = label_parsed) +  
  ylab("Ignition probability") + 
  xlab("Moisture content (%)") +
  
  theme(
    legend.box = "vertical",    # una legend arriba de la otra
    legend.spacing.y = unit(-2, "mm"),
    legend.text = element_text(margin = margin(t = 2, b = 2, l = -2, r = 2, unit = "mm")), # mejora el asunto de separación de símbolosc("#FFFFFF", "#FFFFFF", "#FFFFFF")
    legend.box.margin = margin(-2, 0, 0, 0, "mm"), # era -2
    legend.position = "bottom",
    panel.spacing = unit(3, "mm"),
    strip.text.y = element_text(angle = 270, hjust = 0.5, vjust = 0.5),
    strip.text.x = element_text(size = 0, margin = margin(0,0,0,0, "mm")), # tamaño de la cajita
    strip.background = element_rect(fill = "gray10", color = "gray10"),
    #plot.margin = unit(c(2,1,1,1),"mm")
  )
)

# ggsave doen't work 
png(
  "plots/Fig 07_ignition probability_means.png",
  height = fig_height * 1.15, width = fig_width,
  units = "cm", res = 300,
)
grid.arrange(arrangeGrob(p1 + theme(legend.position = "none") + ylab(NULL),
                         p2 + ylab(NULL),
                         nrow = 2, heights = c(1, 1.18),
                         left = grid::textGrob("Ignition probability", 
                                               rot = 90, vjust = 1)))
dev.off()
# Since the packages were updated, ggsave doen't work well for plots made with 
# ggarrange, adding a gray grided background to images. 


# Fig S6 and S7: air humidity and temperature -------------------------------------

colorcines_air <- viridis(4, option = "B", end = 0.8)


# Humedad

ggplot(clim, aes(x = time, y = hum_min, colour = Community)) +
  geom_line(size = 0.6) + 
  facet_wrap(vars(transect), ncol = 1, strip.position = "right") +
  scale_color_manual(values = colorcines_air) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")), # tamaño de la cajita
        strip.text.y = element_text(margin = margin(0,1,0,1, "mm")),
        panel.spacing = unit(3, "mm"),
        legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")), # mejora el asunto de separación de símbolos
        legend.box.margin = margin(-4, 0, 0, 0, "mm")) + # y esto igual
  ylab("Minimum relative humidity (%)") + 
  xlab("Days since summer solstice")

ggsave("plots/Fig S06_air humidity through time.png", 
       width = fig_width, height = fig_height, units = "cm")

# Temperatura

ggplot(clim, aes(x = time, y = temp_max, colour = Community)) +
  geom_line(size = 0.6) + 
  facet_wrap(vars(transect), ncol = 1, strip.position = "right") +
  scale_color_manual(values = colorcines_air) +
  theme(legend.position = "bottom",
        strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm")), # tamaño de la cajita
        strip.text.y = element_text(margin = margin(0,1,0,1, "mm")),
        panel.spacing = unit(3, "mm"),
        legend.text = element_text(margin = margin(l = - 2, r = 3, unit = "mm")), # mejora el asunto de separación de símbolos
        legend.box.margin = margin(-4, 0, 0, 0, "mm")) + # y esto igual
  ylab("Maximum temperature (°C)") + 
  xlab("Days since summer solstice")

ggsave("plots/Fig S07_air temperature through time.png", 
       width = fig_width, height = fig_height, units = "cm")


# Fig S01: FWI ts ---------------------------------------------------------

fwi_data <- read.table("data_fwi.txt", header = TRUE)
fwi_data$date <- as.Date(fwi_data$date, format = "%Y-%m-%d")

# get 31 day average
fwi_data$mean <- NA
for(d in 16:(nrow(fwi_data) - 15)) {
  fwi_data$mean[d] <- mean(fwi_data$fwi[(d-15):(d+15)])
}

colnames(fwi_data) <- c("date", "daily", "mean")
fwi_long <- pivot_longer(fwi_data, cols = 2:3, names_to = "var", values_to = "fwi")

ggplot(fwi_long, aes(x = date, y = fwi, colour = var)) +
  geom_line(size = 0.5) + 
  scale_color_viridis(discrete = TRUE, labels = c("Daily", "31 days average"),
                      option = "D", end = 0.5) 
# not replaced.

# Fig S8 to 14: FMC through time, data + predictions ------------------------

# flist contains data and predictions for every fuel type
f <- 1
for(f in 1:N_fuels) {
  flist[[f]]$ds_predictions$Community <- factor(flist[[f]]$ds_predictions$Community,
                                            levels = levels(flist[[f]]$ds_predictions$Community),
                                            labels = c("Mature forest", "Young forest",
                                                       "Open shrubland", "Closed shrubland"
                                            ))

  flist[[f]]$ds_predictions$transect <- factor(flist[[f]]$ds_predictions$transect,
                                levels = levels(flist[[f]]$ds_predictions$transect),
                                labels = gsub("t", "T", levels(flist[[f]]$ds_predictions$transect)))

  flist[[f]]$fdata$Community <- factor(flist[[f]]$fdata$Community,
                                            levels = levels(flist[[f]]$fdata$Community),
                                            labels = c("Mature forest", "Young forest",
                                                       "Open shrubland", "Closed shrubland"
                                            ))

  flist[[f]]$fdata$transect <- factor(flist[[f]]$fdata$transect,
                                           levels = levels(flist[[f]]$fdata$transect),
                                           labels = gsub("t", "T", levels(flist[[f]]$fdata$transect)))

}

# relevel shrublands.
for(f in 1:N_fuels) {
  flist[[f]]$ds_predictions$Community <- factor(
    as.character(flist[[f]]$ds_predictions$Community),
    levels = c("Mature forest", "Young forest",
               "Closed shrubland", "Open shrubland")
    )
  
  flist[[f]]$fdata$Community <- factor(
    as.character(flist[[f]]$fdata$Community),
    levels = c("Mature forest", "Young forest",
               "Closed shrubland", "Open shrubland")
  )
}



# 1 h sticks
num <- 1
f <- flist[[num]]$fdata$fuel_type[1]
ggplot(flist[[num]]$ds_predictions, aes(x = time, y = mu_mean, colour = Community)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = Community), 
              colour = NA, alpha = 0.2) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper, fill = Community), 
              colour = NA, alpha = 0.4) +
  geom_line(size = 0.8) +
  geom_point(data = flist[[num]]$fdata, mapping = aes(x = time, y = moisture),
             shape = 1) +
  # geom_point(data = fagg_gam, mapping = aes(x = time, y = moisture),
  #            shape = 1) +
  facet_grid(rows = vars(transect), cols = vars(Community)) +
  #facet_grid(cols = vars(transect)) + 
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
  theme(plot.title = element_text(hjust = 0, vjust = 1)) + 
  ggtitle(f) +
  ylab("Moisture content (%)") +
  xlab("Days since summer solstice") +
  theme(legend.position = "none")

ggsave(paste("plots/Fig S08_fuel moisture time_", f, ".png", sep = ""),
       width = a4_length - 4, height = fig_width, units = "cm")

# 10 h sticks
num <- 2
f <- flist[[num]]$fdata$fuel_type[1]
ggplot(flist[[num]]$ds_predictions, aes(x = time, y = mu_mean, colour = Community)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = Community), 
              colour = NA, alpha = 0.2) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper, fill = Community), 
              colour = NA, alpha = 0.4) +
  geom_line(size = 0.8) +
  geom_point(data = flist[[num]]$fdata, mapping = aes(x = time, y = moisture),
             shape = 1) +
  # geom_point(data = fagg_gam, mapping = aes(x = time, y = moisture),
  #            shape = 1) +
  facet_grid(rows = vars(transect), cols = vars(Community)) +
  #facet_grid(cols = vars(transect)) + 
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
  theme(plot.title = element_text(hjust = 0, vjust = 1)) + 
  ggtitle(f) +
  ylab("Moisture content (%)") +
  xlab("Days since summer solstice") +
  theme(legend.position = "none")

ggsave(paste("plots/Fig S09_fuel moisture time_", f, ".png", sep = ""),
       width = a4_length - 4, height = fig_width, units = "cm")



# Litter
num <- 6
f <- flist[[num]]$fdata$fuel_type[1]; f
ggplot(flist[[num]]$ds_predictions, aes(x = time, y = mu_mean, colour = Community)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = Community), 
              colour = NA, alpha = 0.2) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper, fill = Community), 
              colour = NA, alpha = 0.4) +
  geom_line(size = 0.8) +
  geom_point(data = flist[[num]]$fdata, mapping = aes(x = time, y = moisture),
             shape = 1) +
  # geom_point(data = fagg_gam, mapping = aes(x = time, y = moisture),
  #            shape = 1) +
  facet_grid(rows = vars(transect), cols = vars(Community)) +
  #facet_grid(cols = vars(transect)) + 
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
  theme(plot.title = element_text(hjust = 0, vjust = 1)) + 
  ggtitle(f) +
  ylab("Moisture content (%)") +
  xlab("Days since summer solstice") +
  theme(legend.position = "none")

ggsave(paste("plots/Fig S10_fuel moisture time_", f, ".png", sep = ""),
       width = a4_length - 4, height = fig_width, units = "cm")

# Laura
num <- 5
f <- flist[[num]]$fdata$fuel_type[1]
ggplot(flist[[num]]$ds_predictions, aes(x = time, y = mu_mean, colour = Community)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = Community), 
              colour = NA, alpha = 0.2) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper, fill = Community), 
              colour = NA, alpha = 0.4) +
  geom_line(size = 0.8) +
  geom_point(data = flist[[num]]$fdata, mapping = aes(x = time, y = moisture),
             shape = 1) +
  # geom_point(data = fagg_gam, mapping = aes(x = time, y = moisture),
  #            shape = 1) +
  facet_grid(rows = vars(transect), cols = vars(Community)) +
  #facet_grid(cols = vars(transect)) + 
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
  theme(plot.title = element_text(hjust = 0, vjust = 1)) + 
  ggtitle(expression(italic("S. patagonicus"))) +
  ylab("Moisture content (%)") +
  xlab("Days since summer solstice") +
  theme(legend.position = "none")

ggsave(paste("plots/Fig S11_fuel moisture time_", f, ".png", sep = ""),
       width = a4_length - 4, height = fig_width, units = "cm")


# Coihue
num <- 4
f <- flist[[num]]$fdata$fuel_type[1]
ggplot(flist[[num]]$ds_predictions, aes(x = time, y = mu_mean, colour = Community)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = Community), 
              colour = NA, alpha = 0.2) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper, fill = Community), 
              colour = NA, alpha = 0.4) +
  geom_line(size = 0.8) +
  geom_point(data = flist[[num]]$fdata, mapping = aes(x = time, y = moisture),
             shape = 1) +
  # geom_point(data = fagg_gam, mapping = aes(x = time, y = moisture),
  #            shape = 1) +
  facet_grid(rows = vars(transect), cols = vars(Community)) +
  #facet_grid(cols = vars(transect)) + 
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
  theme(plot.title = element_text(hjust = 0, vjust = 1)) + 
  ggtitle(expression(italic("N. dombeyi"))) +
  ylab("Moisture content (%)") +
  xlab("Days since summer solstice") +
  theme(legend.position = "none")

ggsave(paste("plots/Fig S12_fuel moisture time_", f, ".png", sep = ""),
       width = a4_length - 4, height = fig_width, units = "cm")

# Cania
num <- 3
f <- flist[[num]]$fdata$fuel_type[1]
ggplot(flist[[num]]$ds_predictions, aes(x = time, y = mu_mean, colour = Community)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = Community), 
              colour = NA, alpha = 0.2) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper, fill = Community), 
              colour = NA, alpha = 0.4) +
  geom_line(size = 0.8) +
  geom_point(data = flist[[num]]$fdata, mapping = aes(x = time, y = moisture),
             shape = 1) +
  # geom_point(data = fagg_gam, mapping = aes(x = time, y = moisture),
  #            shape = 1) +
  facet_grid(rows = vars(transect), cols = vars(Community)) +
  #facet_grid(cols = vars(transect)) + 
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
  theme(plot.title = element_text(hjust = 0, vjust = 1)) + 
  ggtitle(expression(italic("C. culeou"))) +
  ylab("Moisture content (%)") +
  xlab("Days since summer solstice") +
  theme(legend.position = "none")

ggsave(paste("plots/Fig S13_fuel moisture time_", f, ".png", sep = ""),
       width = a4_length - 4, height = fig_width, units = "cm")

# Live mixture
num <- 7
f <- flist[[num]]$fdata$fuel_type[1]; f
ggplot(flist[[num]]$ds_predictions, aes(x = time, y = mu_mean, colour = Community)) +
  geom_ribbon(aes(ymin = y_lower, ymax = y_upper, fill = Community), 
              colour = NA, alpha = 0.2) +
  geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper, fill = Community), 
              colour = NA, alpha = 0.4) +
  geom_line(size = 0.8) +
  geom_point(data = flist[[num]]$fdata, mapping = aes(x = time, y = moisture),
             shape = 1) +
  # geom_point(data = fagg_gam, mapping = aes(x = time, y = moisture),
  #            shape = 1) +
  facet_grid(rows = vars(transect), cols = vars(Community)) +
  #facet_grid(cols = vars(transect)) + 
  scale_color_manual(values = colorcines) +
  scale_fill_manual(values = colorcines) +
  theme(plot.title = element_text(hjust = 0, vjust = 1)) + 
  ggtitle(f) +
  ylab("Moisture content (%)") +
  xlab("Days since summer solstice") +
  theme(legend.position = "none")

ggsave(paste("plots/Fig S14_fuel moisture time_", f, ".png", sep = ""),
       width = a4_length - 4, height = fig_width, units = "cm")


# Ns por fuel  ------------------------------------------------------------


fd <- readRDS("/home/ivan/Insync/Humedad y EEA/Datos y scripts/fuel_moisture_data_complete_goodnames_object.R")
aggregate(moisture ~ fuel_type, data = fd[complete.cases(fd$moisture) &
                                          fd$transect != "t4", ], 
          FUN = length)
aggregate(moisture ~ fuel_type, data = fd, FUN = length)


