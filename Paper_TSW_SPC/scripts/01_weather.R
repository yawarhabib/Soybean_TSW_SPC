# ==================================================================================================
#                 Script 01 - Growing-season weather at the Khabarovsk trial site
#              Field trials (2017–2024): daily precipitation and air temperature
# --------------------------------------------------------------------------------------------------
#                           Author :-  Yawar Habib, habib.yawr.10@gmail.com
#                              Center for Bio and Medical Technologies
#                    Skolkovo Institute of Science and Technology, Moscow, Russia
# --------------------------------------------------------------------------------------------------
# Input : daily weather of station WMO 31735 (Khabarovsk airport), 1 May - 30 September, from rp5.ru
#         data/raw/weather/rp5_daily.csv
# Steps : season summary of temperature and precipitation
#         per-season climographs with the crop cycle (Figure 1)
# Output: results/01_weather/season_summary.csv and Fig_weather_by_season.tiff
# ==================================================================================================


### Clean the working environment
rm(list = ls())
graphics.off()
gc()

### Libraries
library(dplyr)
library(ggplot2)
library(zoo)
library(here)


### Load the data from working environment
rp5_daily <- read.csv(here("Paper_TSW_SPC/data/raw/weather/rp5_daily.csv"), header = TRUE)
View(rp5_daily)

str(rp5_daily)


################ Daily weather by season

# Season as the year; day of year on a common calendar axis so that the seasons can be faceted
weather <- rp5_daily %>%
  mutate(date = as.Date(date),
         Year = as.integer(format(date, "%Y")),
         doy  = as.Date(paste0("2000-", format(date, "%m-%d")))) %>%
  arrange(date)

# 7-day centred running means of the daily temperatures, within each season
weather <- weather %>%
  group_by(Year) %>%
  mutate(t_mean_7d = rollmean(t_mean_c, 7, fill = NA, align = "center"),
         t_min_7d  = rollmean(t_min_c,  7, fill = NA, align = "center"),
         t_max_7d  = rollmean(t_max_c,  7, fill = NA, align = "center")) %>%
  ungroup()

# Sowing and emergence of each season, from the field records, on the common calendar axis
phenology <- tribble(~Year, ~sowing, ~emergence,
                     2017,  "05-19", "06-07",
                     2018,  "05-17", "05-30",
                     2019,  "05-27", "06-10",
                     2020,  "05-13", "05-28",
                     2021,  "05-22", "06-11",
                     2022,  "05-18", "06-06",
                     2023,  "05-15", "06-03",
                     2024,  "05-20", "06-10") %>%
  mutate(sowing    = as.Date(paste0("2000-", sowing)),
         emergence = as.Date(paste0("2000-", emergence)))


################ Season summary

# Rain days are days with at least 1 mm; 2017 has one day without a precipitation record
season_summary <- weather %>%
  group_by(Year) %>%
  summarise(Tmean_C     = round(mean(t_mean_c), 1),
            Tmax_C      = max(t_max_c),
            Tmin_C      = min(t_min_c),
            Precip_mm   = round(sum(precip_mm, na.rm = TRUE)),
            RainDays    = sum(precip_mm >= 1, na.rm = TRUE),
            MaxDaily_mm = round(max(precip_mm, na.rm = TRUE), 1),
            .groups = "drop") %>%
  left_join(phenology %>%
              transmute(Year, Sowing = format(sowing, "%d.%m"), Emergence = format(emergence, "%d.%m")),
            by = "Year")
season_summary

write.csv(season_summary, here("Paper_TSW_SPC/results/01_weather/season_summary.csv"), row.names = FALSE)


################ Figure: growing-season weather by season

# Colours as in the trait figure: blue precipitation, red temperature
weather_colors <- c(precip = "#1F77B4", temp = "#CC0000")

# Crop-cycle bands under each panel: sowing to 30 September, and sowing to emergence
bands <- bind_rows(phenology %>%
                     transmute(Year, xmin = sowing, xmax = as.Date("2000-09-30"), phase = "Growing season"),
                   phenology %>%
                     transmute(Year, xmin = sowing, xmax = emergence, phase = "Sowing to emergence")) %>%
  mutate(phase = factor(phase, levels = c("Growing season", "Sowing to emergence")))

band_colors <- c("Growing season" = "#BFE3BC", "Sowing to emergence" = "#4FAE51")

# Season precipitation total and mean temperature printed in each panel
season_labels <- season_summary %>%
  mutate(label = sprintf("%d mm  ·  %.1f°C", Precip_mm, Tmean_C))

# Walter-Lieth scaling: temperature drawn on the precipitation axis at 2 mm per °C
temp_scale <- 2

weather_plot <- ggplot(weather, aes(doy)) +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -9.5, ymax = -3.5, fill = phase),
            inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = t_min_7d * temp_scale, ymax = t_max_7d * temp_scale),
              fill = weather_colors[["temp"]], alpha = 0.16) +
  geom_line(aes(y = t_mean_7d * temp_scale), colour = weather_colors[["temp"]], linewidth = 0.9) +
  geom_col(aes(y = precip_mm), fill = weather_colors[["precip"]], width = 1) +
  geom_text(data = season_labels, aes(as.Date("2000-05-02"), 78, label = label),
            hjust = 0, size = 3.1, colour = "grey25", fontface = "italic", inherit.aes = FALSE) +
  facet_wrap(~ Year, nrow = 2) +
  scale_fill_manual(name = NULL, values = band_colors) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month", expand = expansion(add = 2)) +
  scale_y_continuous(name = "Precipitation (mm/day)", limits = c(-13, 80), breaks = seq(0, 80, 20),
                     expand = expansion(mult = c(0, 0.02)),
                     sec.axis = sec_axis(~ . / temp_scale, name = "Temperature (°C)", breaks = seq(0, 40, 10))) +
  labs(x = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position    = "bottom",
        legend.key.size    = unit(0.9, "lines"),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 13),
        panel.spacing      = unit(1.1, "lines"),
        plot.margin        = margin(6, 12, 4, 12),
        axis.title.y.left  = element_text(face = "bold", margin = margin(r = 6)),
        axis.title.y.right = element_text(face = "bold", margin = margin(l = 6)),
        axis.text.x        = element_text(size = 9),
        axis.line          = element_line(colour = "black"),
        axis.ticks         = element_line(colour = "black"),
        axis.text          = element_text(colour = "black"),
        axis.title         = element_text(colour = "black"))

ggsave(here("Paper_TSW_SPC/results/01_weather/Fig_weather_by_season.tiff"), weather_plot,
       width = 14, height = 7.4, units = "in", dpi = 600, compression = "lzw")
