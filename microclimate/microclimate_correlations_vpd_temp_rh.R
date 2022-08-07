micro_daily <- read.csv("microclimate/data_microclimate_daily.csv")
names(micro_daily)[which(names(micro_daily) == "hum_min")] <- "rh_min"

micro_avg <- aggregate(cbind(vpd_max, temp_max, rh_min) ~ 
                         transect + Community + site, 
                      micro_daily, mean)



# Correlation of daily extremes
pairs(micro_daily[, c("vpd_max", "temp_max", "rh_min")],
      col = rgb(0, 0, 0, 0.1), pch = 19, 
      main = "Daily values")

# Correlation of averages across sites
pairs(micro_avg[, c("vpd_max", "temp_max", "rh_min")],
      col = rgb(0, 0, 0, 0.7), pch = 19, cex = 2, 
      main = "Sites means")
