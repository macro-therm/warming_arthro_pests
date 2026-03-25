library(tidyverse)
library(here)
library(readxl)
library(metafor)
source(here::here("R/S2-TPC_equations.R"))
source(here::here("R/S1-Functions.R"))
# 1. Load -----------------------------------------------------------------
int_rate_data <- read_excel(here("source data/IRTemp.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)

int_rate_data_info <- int_rate_data |> 
  group_by(id_pop) |> 
  slice(1) |> 
  select(reference, species, order, family, lat, lon, id_pop)

thermal_limits <- read_rds(here("export data/thermal_limits_all.rds")) |> 
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

# 2. Estimates  ----
## a) tmin  ----

fit_tmin_meta <- rma.mv(
  yi = tmin_est,
  V = tmin_se^2,
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate,
  control = list(optimizer = "optim", optmethod = "BFGS")
  
)
summary(fit_tmin_meta)

## b) tmax  ----

## metafor 
fit_tmax_meta <- rma.mv(
  yi = tmax_est,
  V = tmax_se^2,
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate,
  control = list(optimizer = "optim", optmethod = "BFGS")
) 
summary(fit_tmax_meta)

## c) topt  ----

fit_topt_meta <- rma.mv(
  yi = topt_est,
  V = topt_se^2,
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate,
  control = list(optimizer = "optim", optmethod = "BFGS")
) 
summary(fit_topt_meta) 

# 3. Low vs. high latitudes  ----
thermal_limits_intrate_latitudes <- thermal_limits_intrate |> 
  mutate(lat_group = case_when(abs(lat) <= 35 ~ "low latitudes",
                               abs(lat) > 35  ~ "high latitudes"))
  ## a) tmin  ----

fit_tmin_meta_latgroups <- rma.mv(
  yi = tmin_est,
  V = tmin_se^2,
  mods = ~ as_factor(lat_group),
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate_latitudes,
  control = list(optimizer = "optim", optmethod = "BFGS")
  
)
summary(fit_tmin_meta_latgroups)

## b) tmax  ----

fit_tmax_meta_latgroups <- rma.mv(
  yi = tmax_est,
  V = tmax_se^2,
  mods = ~ as_factor(lat_group),
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate_latitudes,
  control = list(optimizer = "optim", optmethod = "BFGS")
  
)
summary(fit_tmax_meta_latgroups)

## c) topt  ----

fit_topt_meta_latgroups <- rma.mv(
  yi = topt_est,
  V = topt_se^2,
  mods = ~ as_factor(lat_group),
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate_latitudes,
  control = list(optimizer = "optim", optmethod = "BFGS")
  
)
summary(fit_topt_meta_latgroups)

# 4. Orders  ----
## a) tmin  ----
fit_tmin_order_meta <- rma.mv(
  yi = tmin_est,
  V  = tmin_se^2,
  mods = ~ as_factor(order),
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate,
  control = list(optimizer = "optim", optmethod = "BFGS")
)
summary(fit_tmin_order_meta) 
## b) tmax  ----
fit_tmax_order_meta <- rma.mv(
  yi = tmax_est,
  V  = tmax_se^2,
  mods = ~ as_factor(order),
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate,
  control = list(optimizer = "optim", optmethod = "BFGS")
)
summary(fit_tmax_order_meta) 
## b) topt  ----
fit_topt_order_meta <- rma.mv(
  yi = topt_est,
  V  = topt_se^2,
  mods = ~ as_factor(order),
  random = list(~1|species, ~1|reference),
  data = thermal_limits_intrate,
  control = list(optimizer = "optim", optmethod = "BFGS")
)
summary(fit_topt_order_meta) 



# 4.  visual inspection ----

trait_clouds <- ggplot(thermal_limits_intrate_latitudes |> 
                         filter(is.finite(tmax_est)), aes(x = 1.5)) + 
  ggdist::stat_halfeye(aes(y = tmin_est),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#0a9396",
                       fill = "#94d2bd"
  ) + 
  ggdist::stat_dots(aes(y = tmin_est),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#0a9396",
                    fill= "#94d2bd",
                    alpha = 1
  )+
  coord_cartesian(xlim = c(1, NA))+
  ggdist::stat_halfeye(aes(y = tmax_est),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#bb3e03",
                       fill = "#ee9b00"
  ) + 
  ggdist::stat_dots(aes(y = tmax_est),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#bb3e03",
                    fill= "#ee9b00",
                    alpha = 1
  )+
  ggthemes::theme_few()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = element_blank(),
       y = "Temperature (ºC)")
trait_clouds
ggsave(here("data/data_sink/figures/clouds.svg"), height = 1400, width = 1400,
       units = "px")

# orders
trait_clouds+
  facet_wrap(~order)
ggsave(here("data/data_sink/figures/clouds_tlims_order.png"), height = 1400, width = 1400,
         units = "px")

# latitudial groups
ggplot(thermal_limits_intrate_latitudes |> 
         filter(is.finite(tmax_est)) |> 
         filter(!is.na(lat_group)), aes(x = 1.5)) + 
  ggdist::stat_halfeye(aes(y = tmin_est),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#0a9396",
                       fill = "#94d2bd"
  ) + 
  ggdist::stat_dots(aes(y = tmin_est),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#0a9396",
                    fill= "#94d2bd",
                    alpha = 1
  )+
  coord_cartesian(xlim = c(1, NA))+
  ggdist::stat_halfeye(aes(y = tmax_est),
                       adjust = .5,
                       width = .5, ## set slab interval to show IQR and 95% data range
                       color = "#bb3e03",
                       fill = "#ee9b00"
  ) + 
  ggdist::stat_dots(aes(y = tmax_est),
                    side = "left",
                    dotsize = 1.2,
                    justification = 1.05,
                    binwidth = .5,
                    color = "#bb3e03",
                    fill= "#ee9b00",
                    alpha = 1
  )+
  ggthemes::theme_few()+
  facet_wrap(~lat_group)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  labs(x = element_blank(),
       y = "Temperature (ºC)")

ggsave(here("figures/clouds_tlims_latgroups.png"), height = 1400, width = 1400,
       units = "px")
