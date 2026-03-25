library(tidyverse)
library(here)
library(readxl)
source(here::here("R/S2-TPC_equations.R"))
source(here::here("R/S1-Functions.R"))
library(metafor)

# 1. Load -----------------------------------------------------------------
int_rate_data <- read_excel(here("data/data_source/IRTemp.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)

int_rate_data_info <- int_rate_data |> 
  group_by(id_pop) |> 
  slice(1) |> 
  select(reference, species, order, family, lat, lon, id_pop)

thermal_limits <- read_rds(here("data/data_sink/thermal_limits_all.rds")) |> 
  mutate(id_pop = as_factor(id_pop)) |>
  drop_na() |> 
  group_by(id_pop) |>
  mutate(n_boots = n()) |> 
  summarise(tmin_est = mean(tmin),
            tmax_est = mean(tmax),
            n_boots_est = mean(n_boots),
            tmin_se = sd(tmin),
            tmax_se = sd(tmax),
            topt_est = mean(topt),
            topt_se = sd(topt))

thermal_limits_intrate <- inner_join(int_rate_data_info, thermal_limits) 

# 2. Tmin ~ Lat -----------------------------------------------------------------

## a) metafor fit ----------------------------------------------------------
fit_tmin_lat_meta <- rma.mv(
  yi = tmin_est,
  V  = tmin_se^2,
  mods = ~ abs(lat),
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate,
  method = "REML")

summary(fit_tmin_lat_meta)

fit_tmin_lat_lmer <- lmerTest::lmer(tmin_est ~ abs(lat)+ (1|species) + (1|reference),
               data = thermal_limits_intrate)
summary(fit_tmin_lat_lmer)

## b) plot  -----------------------------------------------------------------
tmin_lat_plot <- sim_and_plot_linears(model_object = fit_tmin_lat_meta,
                                      var_x = abs(thermal_limits_intrate$lat),
                                      var_y = thermal_limits_intrate$tmin_est,
                                      n_sims = 1000,
                                      your_title = "Minimum",
                                      your_subtitle = "N = 238",
                                      lab_x = "absolute latitude (º)",
                                      lab_y = "thermal limit (ºC)",
                                      color_points = "#0a9396",
                                      color_central = "#005f73",
                                      color_uncertainty = "#94d2bd",
                                      model_type = "metafor")
print(tmin_lat_plot)
ggsave(here("data/data_sink/figures/tmin_lat.png"), height = 2600, width = 2600,
       units = "px")


sim_and_plot_linears(model_object = fit_tmin_lat_lmer,
                     var_x = abs(thermal_limits_intrate$lat),
                     var_y = thermal_limits_intrate$tmin_est,
                     n_sims = 1000,
                     your_title = "Minimum",
                     your_subtitle = "N = 238",
                     lab_x = "absolute latitude (º)",
                     lab_y = "thermal limit (ºC)",
                     color_points = "#0a9396",
                     color_central = "#005f73",
                     color_uncertainty = "#94d2bd",
                     model_type = "lmer")


brm_tmin <- brm(
  tmin_est | se(tmin_se) ~ abs(lat) + (1 | species) + (1 | reference),
  data = thermal_limits_intrate,
  family = gaussian(),
  prior = c(
    prior(normal(0, 2), class = "b"),             # slope prior (degrees per latitude)
    prior(normal(0, 10), class = "Intercept"),
    prior(exponential(1), class = "sd")           # random effect SDs
  ),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.95)
)
# 3. tmax ~ Lat -----------------------------------------------------------------

## a) metafor fit ----------------------------------------------------------
fit_tmax_lat_meta <- rma.mv(
  yi = tmax_est,
  V  = tmax_se^2,
  mods = ~ abs(lat),
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate,
  method = "REML")
summary(fit_tmax_lat_meta)

## b) plot  -----------------------------------------------------------------
tmax_lat_plot <- sim_and_plot_linears(model_object = fit_tmax_lat_meta,
                                      var_x = abs(thermal_limits_intrate$lat),
                                      var_y = thermal_limits_intrate$tmax_est,
                                      n_sims = 1000,
                                      your_title = "Maximum",
                                      your_subtitle = "N = 236",
                                      lab_x = "absolute latitude (º)",
                                      lab_y = "thermal limit (ºC)",
                                      color_points = "#bb3e03",
                                      color_central = "#9b2226",
                                      color_uncertainty = "#ee9b00",
                                      model_type = "metafor")
print(tmax_lat_plot)
ggsave(here("data/data_sink/figures/tmax_lat.png"), height = 2600, width = 2600,
       units = "px")

## g) topt  ----
## metafor ----
fit_topt_lat_meta <- rma.mv(
  yi = topt_est,
  V  = topt_se^2,
  mods = ~ abs(lat),
  random = list(~1|species, ~1|reference),
  method = "REML",
  data = thermal_limits_intrate)
summary(fit_topt_lat_meta)

# 4. Composed plots --------------------------------------------------------

## a) pool -----------------------------------------------------------------
tmin_tmax_lat <- cowplot::plot_grid(tmin_lat_plot, tmax_lat_plot)
tmin_tmax_lat
ggsave(here("data/data_sink/figures/tmin_tmax_lat.png"), height = 2600, width = 3200,
       units = "px")
ggsave(here("data/data_sink/figures/tmin_tmax_lat.svg"), height = 1500, width = 1900,
       units = "px")

therm_lims_vert <- thermal_limits_lats |> 
  rename(tmin = tmin_est,
         tmax = tmax_est) |> 
  pivot_longer(cols = c("tmin", "tmax"),
               values_to = "estimates",
               names_to = "thermal_limit")
