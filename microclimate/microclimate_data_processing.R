# Preprocesamiento de datos microclimáticos.
library(plantecophys) # para calcular VPD
library(tidyverse)

clim <- read.csv("data_microclimate_raw.csv", header = TRUE, sep = ";")
head(clim)

summary(clim)
clim$date <- as.Date(clim$date, format = "%Y-%m-%d")
clim$vpd <- with(clim, RHtoVPD(RH = hum, TdegC = temp)) # VPD is in kPa 
# Agrego columna de "día funcional": si son más de las 18, cuenta para el día siguiente.
clim$fDate <- clim$date
past18 <- clim$hour > 18
clim$fDate[past18] <- clim$date[past18] + 1

clim <- clim[, c("datalog", "date", "hour", "fDate", "temp", "hum","vpd")]
hist(clim$vpd, breaks = 50)
hist(log(clim$vpd), breaks = 30)
summary(clim$vpd)

#### agrego por día ####

summaries <- c("min", "median", "mean", "max")
#day_data <- matrix(NA, nrow = 1680, ncol = 2 + length(summaries) * 2)
resumer <- function(x) {
  return(c(min = min(x),
           median = median(x),
           mean = mean(x),
           max = max(x))
         )
}

day_data <- aggregate(cbind(temp, hum, vpd) ~ fDate + datalog, data = clim, 
                      FUN = resumer)
day_data <- do.call("data.frame", day_data)
str(day_data)

colnames(day_data) <- c("fDate", "datalogger", 
                        paste("temp", summaries, sep = "_"),
                        paste("hum", summaries, sep = "_"),
                        paste("vpd", summaries, sep = "_"))

#day_data[, 2] <- as.character(temper$Group.2)
#for(i in 3:ncol(day_data)) day_data[, i] <- round(as.numeric(as.character(day_data[, i])), digits = 1)


# Agrego columnas de sitio y comunidad

day_data$transect <- rep(c("t1", "t1", "t1", "t1", "t2", "t2", "t2", "t2", "t3", "t3"), each = 168)
day_data$Community <- rep(c("B shrubland", "B forest -reg", 
                            "B forest +reg", "U forest", 
                            "B shrubland", "B forest -reg", 
                            "B forest +reg", "U forest", 
                            "B forest +reg", "U forest"), each = 168)
day_data$site <- paste(day_data$Community, day_data$transect, sep = " / ")

day_data$site <- factor(day_data$site, 
                          levels = c("U forest / t1", "U forest / t2", "U forest / t3",
                                     "B forest +reg / t1", "B forest +reg / t2", "B forest +reg / t3",
                                     "B forest -reg / t1", "B forest -reg / t2",
                                     "B shrubland / t1", "B shrubland / t2"))


View(day_data)
str(day_data)
summary(day_data)

ggplot(day_data, aes(x = fDate, y = hum_min, colour = Community)) +
  geom_line() + 
  facet_wrap(vars(transect), nrow = 3)

ggplot(day_data, aes(x = fDate, y = vpd_max, colour = Community)) +
  geom_line() + 
  facet_wrap(vars(transect), nrow = 3)

ggplot(clim, aes(x = hum, y = vpd)) +
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(datalog))

ggplot(clim, aes(x = temp, y = vpd)) +
  geom_point(alpha = 0.3) + 
  facet_wrap(vars(datalog))
# vpd tiene una forma exponencial negativa con rh y positiva con temp,
# y está más asociada a temp (menos ruido)

pairs(day_data[, 3:(ncol(day_data) - 3)])

ori <- as.Date("2019-12-21")
#endnum <- as.numeric(sdates$live_date[nrow(sdates)] - as.Date(ori))
day_data$time <- as.numeric(day_data$fDate - as.Date(ori))

View(day_data)
saveRDS(day_data, "data_microclimate_daily_dataframe.R")
saveRDS(day_data, "data_microclimate_daily.csv")
