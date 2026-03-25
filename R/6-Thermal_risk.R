library(tidyverse)
library(here)
library(readxl)
library(sf)
library(terra)
library(leaflet)
source(here::here("R/S2-TPC_equations.R"))
source(here::here("R/S1-Functions.R"))
library(metafor)
library(patchwork)

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

# 2. Climate data -----------------------------------------------------------------

source_dataset_points_ids <- int_rate_data |> 
  select(lat, lon) |> 
  unique() |> 
  mutate(id_location = 1:n())

list_locations <- int_rate_data |> 
  inner_join(source_dataset_points_ids) |>
  group_by(id_location) |> 
  slice(1) |> 
  filter(!is.na(lon)) 
  

points_intrates <-st_as_sf(x = list_locations,
                           coords = c("lon", "lat")) |> 
  st_set_crs(4326)

#points_vect_intrates <- terra::vect(points_intrates)
present_tmin_climate_wc <-  geodata::worldclim_global(var = "tmin",
                                                      res = 2.5,
                                                      path = tempdir())

points_tmin_present <- terra::extract(present_tmin_climate_wc,
                                      points_intrates, 
                                      id = points_intrates$id_location) |> 
  as_tibble() |> 
  pivot_longer(cols = -ID,
               names_to = "month",
               values_to = "tmin_avg") |> 
  mutate(month = str_sub(month, -2),
         month = as_factor(month))

#### a) Historical ----
present_tmax_climate_wc <-  geodata::worldclim_global(var = "tmax",
                                                      res = 2.5,
                                                      path = tempdir())

points_tmax_present <- terra::extract(present_tmax_climate_wc,
                                      points_intrates, 
                                      id = points_intrates$id_location) |> 
  as_tibble() |> 
  pivot_longer(cols = -ID,
               names_to = "month",
               values_to = "tmax_avg") |> 
  mutate(month = str_sub(month, -2),
         month = as_factor(month))

points_tavg_present <- inner_join(points_tmin_present,
                                  points_tmax_present) |> 
  mutate(tavg = map2_dbl(.x = tmin_avg,
                         .y = tmax_avg,
                         .f = ~mean(c(.x, .y))),
         model = as_factor("WorldClim_1970-2000_v21")) |> 
  rename(tmin = tmin_avg,
         tmax = tmax_avg) |> 
  pivot_longer(cols = 3:5, 
               names_to = "var",
               values_to = "temp_value") |> 
  mutate(time_scenario = as_factor("Present")) |>
  filter(var == "tavg")

### b) Future CMIP6 RCP 4.5 ----

cmip6_avg_tbl <- read_rds(here("source data/cmip6_tavg_2041-2060_ssp245_res25.rds"))

points_tavg_future <- cmip6_avg_tbl |>
  select(ID, month, model, tavg) |> 
  pivot_longer(tavg, names_to = "var", values_to = "temp_value") |> 
  mutate(time_scenario = "Future")
  
### c) Joined (Pres-Future) ----
points_tavg_worldclim <- points_tavg_present |> 
  bind_rows(points_tavg_future)

points_tavg_synth <- points_tavg_worldclim |> 
  group_by(ID, month, time_scenario) |> 
  summarise(temp_value = mean(temp_value)) |> 
  rename(id_location = ID)

int_rate_locations <- int_rate_data |> 
  select(lat, lon) |> 
  unique() |> 
  mutate(id_location = 1:n())

int_rate_data_locations <- int_rate_data |> 
  group_by(lat, lon) |> 
  inner_join(int_rate_locations) |> 
  filter(!is.na(lat)) |> 
  ungroup() |> 
  group_by(id_location) |> 
  mutate(id_location = cur_group_id())

extracted_tavg_intrate <- int_rate_data_locations |> 
  group_by(id_pop, id_location) |> 
  slice(1) |> 
  inner_join(points_tavg_synth) |> 
  select(reference, species, order, family, lat, lon, id_pop, id_location, month, time_scenario, temp_value)


# d) Temperature visualization -------------------------------------------

ggplot(extracted_tavg_intrate |> 
         mutate(time_scenario = if_else(time_scenario == "Present",
                                        "Historic",
                                        "Future")),
       aes(x = temp_value,
           y = as_factor(month),
           fill = after_stat(x))) +
  facet_wrap(~fct_rev(as_factor(time_scenario)))+
  ggridges::geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01,
                                         linewidth = .5, color = "gray76") +
  viridis::scale_fill_viridis(name = "Temperature (ºC)", 
                              option = "B") + 
  labs(
    title = "Temperature Distributions Across Months",
    x = "Temperature (°C)", y = "Month") +
  ggthemes::theme_clean() +
  theme(legend.position = "none")

ggsave(here("data/data_sink/figures/supplementary_figs/ridges.png"),
       width = 1500,
       height = 1500,
       units = "px")

# 3. Project rates -----------------------------------------------------------------
tpcs_selected <- readxl::read_excel(here("export data/tpcs_selection_filters_completed.xlsx")) |> 
  filter(step3_uncertainties == "y")

tpcs_boots <- tibble()

pb <- progress::progress_bar$new(
  format = "Loading Bootstrapped TPCs :percent",
  total = nrow(tpcs_selected),
  clear = F)

for(population_id in tpcs_selected$id_pop) {
  tpc_selected_i <- tpcs_selected |> 
    filter(id_pop == population_id) |> 
    pull(tpc_model)
  int_rate_i <- int_rate_data |> 
    filter(!is.na(int_rate)) |> 
    filter(id_pop == population_id)
  boots_tpc_id_i <- read_rds(paste0(here("export data/boots_tpcs/boots_tpc_"), population_id, ".rds")) |>
    filter(model_name_iter == tpc_selected_i) |> 
    mutate(id_pop = population_id)
  calc_parameters_all <- get_therm_lims(boots_tpc = boots_tpc_id_i,
                              temp = int_rate_i$temperature,
                              int_rate = int_rate_i$int_rate,
                              tpc_model = tpc_selected_i,
                              epsilon = 1e-04)
  calc_parameters <- calc_parameters_all |> 
    drop_na() |> 
    summarise(across(where(is.numeric),
                     ~mean(.x[is.finite(.x)], na.rm = TRUE)))
  
  boots_tpc_id_i_params <- boots_tpc_id_i |> 
    mutate(tmin = calc_parameters$tmin,
           tmax = calc_parameters$tmax,
           topt = calc_parameters$topt,
           t_L50 = calc_parameters$t_L50,
           t_R50 = calc_parameters$t_R50,
           rmax = calc_parameters$rmax,
    )
  tpcs_boots <- bind_rows(tpcs_boots, boots_tpc_id_i_params)
  pb$tick()
}
tpcs_boots

#predict performance based on tpcs_boots and extracted monthly tavg

predict_r <- tibble()
pb <- progress::progress_bar$new(
  format = "Predicting population's performance (historical and future) :percent",
  total = 310,
  clear = F)
tpcs_boots_index <- tpcs_boots |> 
  group_by(id_pop) |> 
  slice(1)
iterations_both <- inner_join(tpcs_boots_index, extracted_tavg_intrate)

for(i in unique(iterations_both$id_pop)) {
  tpcs_boots_i <- tpcs_boots |> 
    filter(id_pop == i)  
  tpc_model_i <- tpcs_boots_i$model_name_iter[1]
  extracted_tavg_i <- extracted_tavg_intrate |> 
    filter(id_pop == i)
  reference_i <- as.character(extracted_tavg_i$reference[1])
  species_i <- as.character(extracted_tavg_i$species[1])
  rmax_i <- tpcs_boots_i$rmax[1]
  tmin_i <- tpcs_boots_i$tmin[1]
  t_L50 <- tpcs_boots_i$t_L50[1]
  topt_i <- tpcs_boots_i$topt[1]
  t_R50 <- tpcs_boots_i$t_R50[1]
  tmax_i <- tpcs_boots_i$tmax[1]
  rmax_i <- tpcs_boots_i$rmax[1]
  if(any(is.na(extracted_tavg_i$temp_value))) {
    warning(paste0("id_pop ", i, " had no climate data and was discarded"))
    next
  }
  if(all(is.na(tpcs_boots_i$tmin))) {
    warning(paste0("id_pop ", i, " yield unreliable estimates of Tmin. Please consider selecting a different TPC"))
    
    next}
  
  if(all(is.na(tpcs_boots_i$tmax))) {
    warning(paste0("id_pop ", i, " yield unreliable estimates of Tmax"))
    
    next}
  #loop for each month
  predict_r_pop <- tibble()
  for(month_i in unique(extracted_tavg_intrate$month)) {
    monthly_ext_tavg_i <- extracted_tavg_i |> 
      filter(month == month_i)
    #loop for each time scenario
    predict_r_month <- tibble()
    for(time_i in c("Present", "Future")){
      monthly_ext_tavg_time_i <- monthly_ext_tavg_i |> 
        filter(time_scenario == time_i)
      monthly_time_tvalue <- monthly_ext_tavg_time_i$temp_value
      predicted_r_time_i <- tpcs_boots_i |>
        filter(abs(preds) < 1) |> 
        group_by(id_pop, temp) |> 
        summarise(preds = mean(preds, na.rm = TRUE), .groups = "drop") |> 
        slice_min(temp <= monthly_time_tvalue) |> 
        slice(1) |> 
        select(id_pop, preds)
      if (monthly_time_tvalue < tmin_i || monthly_time_tvalue > tmax_i) {
        predicted_r_time_i$preds <- 0
      }
      predict_r_month_i <- inner_join(monthly_ext_tavg_time_i, predicted_r_time_i, by = "id_pop")
      predict_r_month <- bind_rows(predict_r_month, predict_r_month_i)
    }
    predict_r_month <- predict_r_month |> 
      mutate(rmax = tpcs_boots_i$rmax[1],
             tmin = tpcs_boots_i$tmin[1],
             t_L50 = tpcs_boots_i$t_L50[1],
             topt = tpcs_boots_i$topt[1],
             t_R50 = tpcs_boots_i$t_R50[1],
             tmax = tpcs_boots_i$tmax[1],
             rmax = tpcs_boots_i$rmax[1])
    
    predict_r_pop <- bind_rows(predict_r_pop, predict_r_month)
  }
  predict_r <- bind_rows(predict_r, predict_r_pop)
  pb$tick()
}

predict_r_shift <- predict_r |>
  group_by(reference, species, order, family, lat, lon, id_pop, id_location, month) |> 
  summarise(temp_present = temp_value[time_scenario == "Present"],
            temp_future = temp_value[time_scenario == "Future"],
            preds_present = preds[time_scenario == "Present"],
            preds_future  = preds[time_scenario == "Future"],
            preds_diff = preds_future - preds_present,
            psi = preds_diff/rmax,
            r_max = rmax,
            curve_zone_present = case_when(temp_present < tmin ~ "CE",
                                           temp_present < t_L50 &
                                             temp_present >= tmin ~ "CLP",
                                           temp_present >= t_L50 &
                                             temp_present < topt ~ "OPS",
                                           temp_present >= topt & 
                                             temp_present < t_R50 ~ "OPD",
                                           temp_present >= t_R50 &
                                             temp_present < tmax ~ "HLP",
                                           temp_present > tmax ~ "HE"),
            curve_zone_future = case_when(temp_future < tmin ~ "CE",
                                          temp_future < t_L50 &
                                            temp_future >= tmin ~ "CLP",
                                          temp_future >= t_L50 &
                                            temp_future < topt ~ "OPS",
                                          temp_future >= topt & 
                                            temp_future < t_R50 ~ "OPD",
                                          temp_future >= t_R50 &
                                            temp_future < tmax ~ "HLP",
                                          temp_future >= tmax ~ "HE")
           ) |> 
  slice(1) |> 
  ungroup()

r_max_prop <- predict_r_shift |> 
  group_by(id_pop) |> 
  summarise(preds_present = mean(preds_present, na.rm = TRUE),
            preds_future = mean(preds_future, na.rm = TRUE),
            r_max = mean(r_max, na.rm = TRUE),
            r_hist_perc = (100*preds_present)/r_max,
            r_future_perc = (100*preds_future)/r_max,)
r_max_prop_overall <- r_max_prop |> 
  summarise(r_hist_perc_est = mean(r_hist_perc),
            r_hist_perc_se = sd(r_hist_perc)/sqrt(nrow(r_max_prop)),
            r_future_perc_est = mean(r_future_perc),
            r_future_perc_se = sd(r_future_perc, na.rm = TRUE)/sqrt(nrow(r_max_prop)))

r_count_zones <- predict_r_shift |> 
  select(curve_zone_present, 
         curve_zone_future) |> 
  pivot_longer(cols = 1:2,
               names_to = "time_scenario",
               values_to = "curve_zone") |> 
  mutate(time_scenario = str_sub(time_scenario, 12, -1L)) |> 
  group_by(time_scenario) |> 
  count(curve_zone) |> 
  mutate(n = 100*n/2772) |> 
  mutate(time_scenario = ifelse(time_scenario == "present",
                                "historical",
                                "RCP 4.5 (2041-2060)")) |> 
  mutate(sort = case_when(curve_zone == "CE" ~ 1,
                          curve_zone == "CLP" ~ 2,
                          curve_zone == "OPS" ~ 3,
                          curve_zone == "OPD" ~ 4,
                          curve_zone == "HLP" ~ 5,
                          curve_zone == "HE" ~ 6),
         colorado = case_when(curve_zone == "CE" ~ "#B2F2FD",
                              curve_zone == "CLP" ~ "#61DCA9",
                              curve_zone == "OPS" ~ "#8DB63B",
                              curve_zone == "OPD" ~ "#9B7424",
                              curve_zone == "HLP" ~ "#944046",
                              curve_zone == "HE" ~ "#8C0172"))

ggplot(r_count_zones,
       aes(x = time_scenario, 
           stratum = fct_reorder(curve_zone, sort),
           alluvium = fct_reorder(curve_zone, sort),
           y = n,
           fill = colorado)) +
  ggalluvial::geom_flow(stat = "alluvium", 
                        lode.guidance = "frontback",
                        curve_type = "cubic",
                        alpha = 0.4) +
  ggalluvial::geom_stratum(width = 0.35, color = "gray76") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_identity()+
  labs(x = "Scenario",
       y = element_blank())

ggsave(here("data/data_sink/figures/supplementary_figs/flowchart_curvezones.svg"),
       width = 2800,
       height = 280,
       units = "px")



write_rds(predict_r_shift, here("export data/predictions_psi.rds"))


# 4. Analysis risk -----------------------------------------------------------------
monthly_r_shift <- read_rds(here("export data/predictions_psi.rds")) |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri") 

summary(monthly_r_shift$psi)

##plot the distribution of psi's
ggplot(monthly_r_shift, aes(y = 1.5)) + 
  ggdist::stat_slab(
    aes(x = psi, fill = after_stat(x)),
    adjust = .5,
    fill_type = "segments",
    color = "gray69",
    linewidth = .5
  )+
  coord_cartesian(ylim = c(1.5, NA))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(y = element_blank(),
       x = "Temperature (ºC)")+
  scale_fill_gradientn(colours = khroma::colour("vik")(256),
                       values = scales::rescale(c(-0.946331, 0, 0.433069)),
                       limits = c(-0.946331, 0.433069))

predict_r_shift <- monthly_r_shift |> 
  group_by(reference, species, order, family, id_pop, id_location) |> 
  summarise(psi = mean(psi)) |> 
  mutate(win_or_lose = case_when(psi < 0 ~ as_factor("losers"),
                                 psi >= 0 ~ as_factor("winners")))
summary(predict_r_shift$psi)


predict_r_shift
summary(predict_r_shift$psi)
predict_r_shift |> 
  ungroup() |> 
  count(win_or_lose) |> 
  mutate(n = n*100/231)

predict_r_shift |> 
  ungroup() |> 
  group_by(order) |> 
  summarise(psi = mean(psi))

predict_r_shift <- read_rds(here("export data/predictions_psi.rds")) |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri")

predict_r_shift_preds <- predict_r_shift |>
  group_by(family) |> 
  summarise(r_hist = mean(preds_present),
            r_fut = mean(preds_future)) |> 
  mutate(dt_hist = log(2)/r_hist,
         dt_fut = log(2)/r_fut) |> 
  mutate(dt_magn = dt_fut - dt_hist)


join_psi_lat <- read_rds(here("export data/predictions_psi.rds")) |> 
  group_by(reference, species, order, family, id_pop, id_location) |> 
  summarise(psi = mean(psi),
            lat = mean(lat))  


## winners/losers lat groups
psi_lat_groups_df <- read_rds(here("export data/predictions_psi.rds")) |> 
  group_by(reference, species, order, family, id_pop, id_location, month) |> 
  summarise(psi = mean(psi),
            lat = mean(lat)) |> 
  mutate(win_or_lose = case_when(psi < 0 ~ as_factor("losers"),
                                 psi >= 0 ~ as_factor("winners")),
         lat_group = case_when(abs(lat) <=35 ~ "Low-Latitude",
                               abs(lat) > 35 ~ "High-Latitude")) |> 
  group_by(lat_group) |> 
  count(win_or_lose))


# 5. Visualization maps -----------------------------------------------------------------


## a) Static all ---------------------------------------------------------------


predict_r_shift <- read_rds(here("export data/predictions_psi.rds"))


predict_r_shift_sf <- predict_r_shift |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri") |> 
  group_by(reference, species, order, family, id_pop, id_location) |> 
  summarise(psi = mean(psi)) |> 
  mutate(win_or_lose = case_when(psi < 0 ~ as_factor("losers"),
                                 psi >= 0 ~ as_factor("winners")))

worldmap <- rnaturalearth::ne_countries(scale = 50) |> 
  st_transform(crs = "+proj=robin")


ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom")+
  khroma::scale_color_bam(reverse = T, midpoint = 0)+
  facet_wrap(~forcats::fct_rev(win_or_lose), ncol = 1)+
  labs(color = expression(italic(psi)[italic(r)]))



ggsave(here("data/data_sink/figures/supplementary_figs/winners_losers.png"),
       width = 2800,
       height = 3600,
       units = "px")
palette_rdyblu <- hcl.colors(100, palette = "RdYlBu", rev = TRUE)


ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray16"),
        panel.grid = element_line(colour = "gray2"))+
  scale_color_gradientn(
    colors = palette_rdyblu,
    values = scales::rescale(c(min(predict_r_shift_sf$psi, na.rm=TRUE),
                               0,
                               max(predict_r_shift_sf$psi, na.rm=TRUE)))
    )+
  labs(color = expression(italic(psi)[italic(r)]))


ggsave(here("data/data_sink/figures/supplementary_figs/psi.svg"),
       width = 2800,
       height = 2800,
       units = "px")


#small sample
ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray16"),
        panel.grid = element_line(colour = "gray2"))+
  khroma::scale_color_bam(reverse = T, midpoint = 0)+
  labs(color = expression(italic(psi)[italic(r)]))


ggsave(here("data/data_sink/figures/supplementary_figs/psi_sample.png"),
       width = 21000,
       height = 21000,
       units = "px", dpi = 600)

## b) Static no 0s ---------------------------------------------------------------

predict_r_shift_sf_0s <- predict_r_shift_0s |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri") |> 
  group_by(reference, species, order, family, id_pop, id_location) |> 
  summarise(psi = mean(psi)) |> 
  mutate(win_or_lose = case_when(psi < 0 ~ as_factor("losers"),
                                 psi >= 0 ~ as_factor("winners")))

worldmap <- rnaturalearth::ne_countries(scale = 50) |> 
  st_transform(crs = "+proj=robin")


ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf_0s,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray16"),
        panel.grid = element_line(colour = "gray2"))+
  khroma::scale_color_bam(reverse = T, midpoint = 0)+
  facet_wrap(~forcats::fct_rev(win_or_lose), ncol = 1)+
  labs(color = expression(italic(psi)[italic(r)]))



ggsave(here("figures/supporting information figures/winners_losers_0s.png"),
       width = 2800,
       height = 3600,
       units = "px")
palette_rdyblu <- hcl.colors(100, palette = "RdYlBu", rev = TRUE)


ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf_0s,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.background = element_rect(fill = "gray16"),
        panel.grid = element_line(colour = "gray2"))+
  scale_color_gradientn(
    colors = palette_rdyblu,
    values = scales::rescale(c(min(predict_r_shift_sf$psi, na.rm=TRUE),
                               0,
                               max(predict_r_shift_sf$psi, na.rm=TRUE)))
  )+
  labs(color = expression(italic(psi)[italic(r)]))


ggsave(here("figures/supporting information figures/psi_0s.svg"),
       width = 2800,
       height = 2800,
       units = "px")


#small sample
ggplot()+
  geom_sf(data = worldmap, fill = "gray45", color = NA)+
  geom_sf(data = predict_r_shift_sf_0s,
          aes(color = psi))+
  theme_bw()+
  theme(legend.position = "bottom")+
  khroma::scale_color_bam(reverse = T, midpoint = 0)+
  labs(color = expression(italic(psi)[italic(r)]))


ggsave(here("figures/supporting information figures/psi_sample_0s.png"),
       width = 21000,
       height = 21000,
       units = "px", dpi = 600)


## c) Interactive ---------------------------------------------------------------

#map it with leaflet
predict_r_shift <- read_rds(here("export data/predictions_psi.rds"))

for(i in unique(predict_r_shift$id_pop)){
  predict_r_shift_i <- predict_r_shift  |> 
    filter(id_pop == i) 
  predict_r_shift_curvezone_future <- predict_r_shift_i |> 
    group_by(reference, species, order, family, lat, lon, id_pop, id_location) |> 
    count(curve_zone_future) |> 
    mutate(sort = case_when(curve_zone_future == "CE" ~ 1,
                            curve_zone_future == "CLP" ~ 2,
                            curve_zone_future == "OPS" ~ 3,
                            curve_zone_future == "OPD" ~ 4,
                            curve_zone_future == "HLP" ~ 5,
                            curve_zone_future == "HE" ~ 6),
           colorado = case_when(curve_zone_future == "CE" ~ "#B2F2FD",
                                curve_zone_future == "CLP" ~ "#61DCA9",
                                curve_zone_future == "OPS" ~ "#8DB63B",
                                curve_zone_future == "OPD" ~ "#9B7424",
                                curve_zone_future == "HLP" ~ "#944046",
                                curve_zone_future == "HE" ~ "#8C0172"))
  predict_r_shift_curvezone_present <- predict_r_shift_i |> 
    group_by(reference, species, order, family, lat, lon, id_pop, id_location) |> 
    count(curve_zone_present) |> 
    mutate(sort = case_when(curve_zone_present == "CE" ~ 1,
                            curve_zone_present == "CLP" ~ 2,
                            curve_zone_present == "OPS" ~ 3,
                            curve_zone_present == "OPD" ~ 4,
                            curve_zone_present == "HLP" ~ 5,
                            curve_zone_present == "HE" ~ 6),
           colorado = case_when(curve_zone_present == "CE" ~ "#B2F2FD",
                                curve_zone_present == "CLP" ~ "#61DCA9",
                                curve_zone_present == "OPS" ~ "#8DB63B",
                                curve_zone_present == "OPD" ~ "#9B7424",
                                curve_zone_present == "HLP" ~ "#944046",
                                curve_zone_present == "HE" ~ "#8C0172"))
  predict_r_shift_curvezone <- bind_rows(predict_r_shift_curvezone_present,
                                         predict_r_shift_curvezone_future) |> 
    pivot_longer(cols = c(curve_zone_present, curve_zone_future),
                 values_to = "curve_zone",
                 names_to = "time_scenario") |> 
    mutate(time_scenario = str_sub(time_scenario, 12, -1L)) |> 
    mutate(time_scenario = ifelse(time_scenario == "present",
                                  "a) historical",
                                  "b) future")) |> 
    tidyr::drop_na()
  
  species_name <- parse_character(predict_r_shift_curvezone$species)[1]
  reference_name <- as.character(predict_r_shift_curvezone$reference)[1]
  
  barplot_risk_study <- ggplot(data = predict_r_shift_curvezone)+
    geom_bar(aes(x = fct_reorder(as_factor(curve_zone), sort), 
                 y = n, 
                 fill =  as_factor(colorado)),
             stat = "identity"
    )+
    scale_fill_identity()+
    theme_minimal()+
    theme(legend.position = "bottom")+
    labs(title = species_name,
         x = element_blank(),
         y = "#months",
         caption = NULL)+
    theme(plot.title = element_text(face = "italic"))+
    facet_wrap(~time_scenario)+
    theme(panel.spacing = unit(3, "lines"))  # default is usually 0.5 lines

  ##psi barplots
  
  y_title <- expression(psi[r])
  fill_title <- expression(psi[r])
  barplot_psi_study <- ggplot(data = predict_r_shift_i)+
    geom_bar(aes(x = month, 
                 y = psi, 
                 fill = psi),
             stat = "identity")+
    scale_fill_gradientn(colours = khroma::colour("vik")(256),
                         values = scales::rescale(c(-0.946331, 0, 0.433069)),
                         limits = c(-0.946331, 0.433069))+
    theme_minimal()+
    theme(legend.position = "bottom")+
    labs(x = "Month",
         y = y_title,
         fill = fill_title,
         caption = reference_name)+
    theme(plot.title = element_text(face = "italic"))+
    theme(panel.spacing = unit(3, "lines"))  # default is usually 0.5 lines
  
  barplot_combined_study <- (barplot_risk_study / barplot_psi_study)+
    plot_annotation(tag_levels = "a")
  save(barplot_combined_study, 
       file = here(paste0("figures/supporting information figures/riskpoints/barplots_combined_study",i,".RData")))
  ggsave(here(paste0("figures/supporting information figures/riskpoints/barplots_combined_study",i,".jpg")),
         width = 12,
         height = 12,
         units = "cm")
}

list_ggbarplots <- tibble(id = NULL,
                          ggbarplot = NULL)
for(pop_i in unique(predict_r_shift$id_pop)){
  load(here(paste0("figures/supporting information figures/riskpoints/barplots_combined_study", pop_i, ".RData")))
  ggbarplot_study <- tibble(id_pop = pop_i, ggbarplot = list(barplot_combined_study)) 
  list_ggbarplots <- bind_rows(list_ggbarplots, ggbarplot_study)
}

leaflet_risk_shift <- predict_r_shift |>  
  inner_join(list_ggbarplots, by = "id_pop") %>% 
  mutate(icons = case_when(curve_zone_present == "CE" ~ "#B2F2FD",
                           curve_zone_present == "CLP" ~ "#61DCA9",
                           curve_zone_present == "OPS" ~ "#8DB63B",
                           curve_zone_present == "OPD" ~ "#9B7424",
                           curve_zone_present == "HLP" ~ "#944046",
                           curve_zone_present == "HE" ~ "#8C0172")
  ) |> 
  group_by(id_pop) |> 
  mutate(psi = mean(psi))

mytext <- paste(
  "<i>Click to see monthly distribution of performance risk</i>", "<br/>",
  "<i>Species</i>: ", leaflet_risk_shift$species,"<br/>", 
  "<i>Order</i>: ", leaflet_risk_shift$order,"<br/>", 
  "<i>Family</i>: ", leaflet_risk_shift$family,"<br/>", 
  "<i>Reference</i>: ", leaflet_risk_shift$reference, "<br/>", 
  "<i>Performance Shift Index </i>: ", leaflet_risk_shift$psi, "<br/>",
  "<i>Performance raw difference </i>: ", leaflet_risk_shift$preds_diff, "<br/>",
  "<i>Performance historic </i>: ", leaflet_risk_shift$preds_present, "<br/>",
  "<i>Performance historic (%)</i>: ", leaflet_risk_shift$preds_future, "<br/>",
  
  
  sep="") %>%
  lapply(htmltools::HTML)
colors_riskcurve <- c("#B2F2FD", "#61DCA9", "#8DB63B", "#9B7424", "#944046", "#8C0172")

colors_shift <- palette_rdyblu

palette_riskcurve <- colorFactor(palette = colors_riskcurve,
                                 domain = unique(leaflet_risk_shift$curve_zone_future))
palette_shift <- colorNumeric(palette = colors_shift,
                              domain = unique(leaflet_risk_shift$psi))

intrapests_future_risk_leaflet <- leaflet(data = leaflet_risk_shift) %>% 
  addTiles() %>% 
  addCircleMarkers(lng = ~lon, 
                   lat = ~lat, 
                   stroke = FALSE,
                   fillColor = ~palette_shift(psi), 
                   fillOpacity = 0.8, 
                   color = ~palette_shift(psi), 
                   popup = ~reference,
                   label = mytext,
                   group = "reference",
                   labelOptions = labelOptions( 
                     style = list("font-weight" = "normal", 
                                  padding = "3px 8px"), 
                     textsize = "13px", 
                     direction = "auto")) %>% 
  leafpop::addPopupGraphs(leaflet_risk_shift$ggbarplot , 
                          group = "reference", 
                          width = 400, height = 450) %>% 
  addProviderTiles('CartoDB.DarkMatterNoLabels') #%>%
#addLegend(pal= palette_shift, 
#          values= ~psi, 
#          opacity = 0.9,
#          title = "Performance shift",
#          position = "bottomleft")
intrapests_future_risk_leaflet <- intrapests_future_risk_leaflet %>% 
  addLegend(pal= palette_shift, 
            values= ~psi, 
            opacity = 0.9,
            title = "Performance shift",
            position = "bottomleft")

save(intrapests_future_risk_leaflet, file = here("figures/supporting information figures/intrapests_future_risk_leaflet.RData"))
htmlwidgets::saveWidget(intrapests_future_risk_leaflet,
           file = here("index.html"))


# 6. Examples Figure insets TPC -------------------------------------------


tpcs_selected_inset <- readxl::read_excel(here("data/data_sink/tpcs_selection_filters_completed.xlsx")) |> 
  filter(step3_uncertainties == "y") |> 
  filter(id_pop %in% c(91, 93, 223))

simulated_tpcs <- tibble()


for(i in c(223, 93, 91)) {
  cat(paste0("Beginning population ", i, " of 313\n"))
  int_rate_i <-int_rate_data |> 
    filter(id_pop == i) |> 
    filter(!is.na(int_rate))
  species_i <- int_rate_i |> 
    slice(1) |> 
    pull(species)
  reference_i <- int_rate_i |> 
    slice(1) |> 
    pull(reference) |> 
    as.character()
  tpcs_selected_inset_i <- tpcs_selected_inset |> 
    filter(id_pop == i)
  model_name_i <- tpcs_selected_inset_i$tpc_model
  fitted_tpc_equations_i <- fit_tpcs(temp = int_rate_i$temperature,
                                     int_rate = int_rate_i$int_rate,
                                     model_name = model_name_i)
  possible_error <- tryCatch(expr = {
    sim_tpcs_i <- predict_curves(temp = int_rate_i$temperature,
                                 int_rate = int_rate_i$int_rate,
                                 fitted_parameters = fitted_tpc_i,
                                 model_name_2boot = fitted_tpc_equations_i,n_boots_samples = 100)
    plot_uncertainties(bootstrap_uncertainties_tpcs = sim_tpcs_i,
                       temp = int_rate_i$temperature,
                       int_rate = int_rate_i$int_rate,
                       species = species_i,
                       reference = reference_i,
                       pop_id = i)
    ggsave(paste0(here("data/data_sink/figures/supplementary_figs/inset"), i, ".svg"),
           width = 2100,
           height = 2100,
           units = "px")
  }, # <- inside tryCatch
  error = function(e) e)
  if (inherits(possible_error, "error")) {
    fit_nls <- NULL
    warning(paste0("Reference ID", i, "(",reference_i, ") could not fit bootstrapped TPCs"))
  }
}

# 7. Th. safety margins -------------------------------------------
predict_r_shift_lats <- read_rds(here("export data/predictions_psi.rds"))
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


# a) Third coldest months -----------------------------------------------------
predict_r_coldest <- predict_r_shift_lats |> 
  group_by(id_pop) |> 
  slice_min(temp_present, n = 3) |> 
  summarise(temp_present = mean(temp_present),
            species = unique(species),
            reference = unique(reference),
            lat = mean(lat)) |> 
  inner_join(thermal_limits) |> 
  mutate(csm = temp_present-tmin_est)

fit_csm_lat_meta_quarter <- rma.mv(
  yi = csm,
  V  = tmin_se^2,
  mods = ~ abs(lat),
  random = list(~1|species, ~1|reference),
  data = predict_r_coldest,
  control = list(optimizer = "optim", optmethod = "BFGS")
)
summary(fit_csm_lat_meta_quarter)
csm_lat_plot_coldest <- sim_and_plot_linears(model_object = fit_csm_lat_meta_quarter,
                                     var_x = abs(predict_r_coldest$lat),
                                     var_y = predict_r_coldest$csm,
                                     n_sims = 1000,
                                     your_title = "Cold",
                                     your_subtitle = "Three coldest months",
                                     lab_x = "absolute latitude (º)",
                                     lab_y = "Safety Margin (ºC)",
                                     color_points = "#9d4edd",
                                     color_central = "#3c096c",
                                     color_uncertainty = "#e0aaff",
                                     model_type = "metafor")
print(csm_lat_plot_coldest)
ggsave(here("data/data_sink/figures/csm_lat_quarter.png"), height = 1500, width = 1500,
       units = "px")

# b) Third warmest months -----------------------------------------------------
predict_r_warmest <- predict_r_shift_lats |> 
  group_by(id_pop) |> 
  slice_max(temp_present, n = 3) |> 
  summarise(temp_present = mean(temp_present),
            species = unique(species),
            reference = unique(reference),
            lat = mean(lat)) |> 
  inner_join(thermal_limits) |> 
  mutate(wsm = tmax_est-temp_present)

fit_wsm_lat_meta_quarter <- rma.mv(
  yi = wsm,
  V  = tmax_se^2,
  mods = ~ abs(lat),
  random = list(~1|species, ~1|reference),
  data = predict_r_warmest,
  control = list(optimizer = "optim", optmethod = "BFGS")
)
summary(fit_wsm_lat_meta_quarter)

wsm_lat_plot_warmest <- sim_and_plot_linears(model_object = fit_wsm_lat_meta_quarter,
                                     var_x = abs(predict_r_warmest$lat),
                                     var_y = predict_r_warmest$wsm,
                                     n_sims = 1000,
                                     your_title = "Warm",
                                     your_subtitle = "Three warmest months",
                                     lab_x = "absolute latitude (º)",
                                     lab_y = NULL,
                                     color_points = "#e85d04",
                                     color_central = "#9b2226",
                                     color_uncertainty = "#faa307",
                                     model_type = "metafor")
print(wsm_lat_plot_warmest)

ggsave(here("data/data_sink/figures/wsm_lat_quarter.png"), height = 2600, width = 2600,
       units = "px")
tsm_lat_plot_quarters <- cowplot::plot_grid(csm_lat_plot_coldest, wsm_lat_plot_warmest)
print(tsm_lat_plot_quarters)
ggsave(here("data/data_sink/figures/tsm_lat_plot_quarters.png"), height = 1500, width = 1500,
       units = "px")
ggsave(here("data/data_sink/figures/tsm_lat_plot_quarters.svg"), height = 1500, width = 1900,
       units = "px")


# c) average year -----------------------------------------------------
predict_r_avg_tsm <- predict_r_shift_lats |> 
  group_by(id_pop) |> 
  summarise(temp_present = mean(temp_present),
            species = unique(species),
            reference = unique(reference),
            lat = mean(lat)) |> 
  inner_join(thermal_limits) |> 
  mutate(csm = temp_present-tmin_est,
         wsm = tmax_est-temp_present)

fit_wsm_lat_meta <- rma.mv(
  yi = wsm,
  V  = tmax_se^2,
  mods = ~ abs(lat),
  random = list(~1|species, ~1|reference),
  data = predict_r_avg_tsm,
  control = list(optimizer = "optim", optmethod = "BFGS")
)
summary(fit_wsm_lat_meta)
wsm_lat_plot <- sim_and_plot_linears(model_object = fit_wsm_lat_meta,
                                      var_x = abs(predict_r_avg_tsm$lat),
                                      var_y = predict_r_avg_tsm$wsm,
                                      n_sims = 1000,
                                      your_title = "Warm",
                                      your_subtitle = "Yearly average",
                                      lab_x = "absolute latitude (º)",
                                      lab_y = NULL,
                                     color_points = "#e85d04",
                                     color_central = "#9b2226",
                                     color_uncertainty = "#faa307",
                                      model_type = "metafor")
print(wsm_lat_plot)

fit_csm_lat_meta <- rma.mv(
  yi = csm,
  V  = tmin_se^2,
  mods = ~ abs(lat),
  random = list(~1|species, ~1|reference),
  data = predict_r_avg_tsm,
  control = list(optimizer = "optim", optmethod = "BFGS")
)
summary(fit_csm_lat_meta)
csm_lat_plot <- sim_and_plot_linears(model_object = fit_csm_lat_meta,
                                      var_x = abs(predict_r_avg_tsm$lat),
                                      var_y = predict_r_avg_tsm$csm,
                                      n_sims = 1000,
                                      your_title = "Cold",
                                      your_subtitle = "Yearly average",
                                      lab_x = "absolute latitude (º)",
                                      lab_y = "Safety Margin (ºC)",
                                     color_points = "#9d4edd",
                                     color_central = "#3c096c",
                                     color_uncertainty = "#e0aaff",
                                      model_type = "metafor")
print(csm_lat_plot)

tsm_lat_plot <- cowplot::plot_grid(csm_lat_plot, wsm_lat_plot)
print(tsm_lat_plot)
ggsave(here("data/data_sink/figures/tsm_lat_plot.png"), height = 1500, width = 1500,
       units = "px")
ggsave(here("data/data_sink/figures/tsm_lat_plot.svg"), height = 1500, width = 1900,
       units = "px")

# 8. Transitions -------------------------------------------
predict_r_shift_trans <- predict_r_shift_lats |> 
  mutate(trans = case_when(curve_zone_present == "CE" &
                             curve_zone_future == "CLP" ~ "cold_release",
                           curve_zone_present == "CLP" &
                             curve_zone_future == "OPS" ~ "cold_release",
                           curve_zone_present == "CE" &
                             curve_zone_future == "OPS" ~ "cold_release",
                           curve_zone_present == "OPS" &
                             curve_zone_future == "OPD" ~ "thermal_risk_entry",
                           curve_zone_present == "OPS" &
                             curve_zone_future == "HLP" ~ "thermal_risk_entry",
                           curve_zone_present == "CLP" &
                             curve_zone_future == "OPD" ~ "thermal_risk_entry",
                           curve_zone_present == "CLP" &
                             curve_zone_future == "HLP" ~ "thermal_risk_entry",
                           .default = NA)) |> 
  drop_na() |> 
  group_by(id_pop, lat, species, reference) |> 
  count(trans) |> 
  mutate(trans = as_factor(trans))
print(predict_r_shift_trans)

cold_test <- glmmTMB::glmmTMB(n ~ abs(lat) + (1|species) + (1|reference),
                              data = predict_r_shift_trans |> 
                                filter(trans == "cold_release"),
                              family = poisson)
  
  
summary(cold_test)

warm_test <- glmmTMB::glmmTMB(n ~ abs(lat) + (1|species) + (1|reference),
                              data = predict_r_shift_trans |> 
                                filter(trans == "thermal_risk_entry"),
                              family = poisson)


summary(warm_test)

### tropical vs. temperate

trans_lat_groups <- predict_r_shift_trans |> 
  group_by(id_pop) |> 
  mutate(lat_group = case_when(abs(lat) <= 23.5 ~ "tropical",
                               abs(lat) > 23.5 ~ "temperate")) 
trans_lat_summary <- trans_lat_groups |> 
  group_by(lat_group, trans) |> 
  summarise(total_n = sum(n),
            n_pops = n_distinct(id_pop),
            .groups = "drop")

# Create contingency table
cont_table <- trans_lat_groups |> 
  filter(trans %in% c("cold_release", "thermal_risk_entry")) |> 
  group_by(lat_group, trans) |> 
  summarise(total_n = sum(n), 
            .groups = "drop")  |> 
  tidyr::pivot_wider(names_from = trans, 
                     values_from = total_n, 
                     values_fill = 0)  |> 
  column_to_rownames("lat_group")

cont_table
chisq.test(cont_table)


cold_glm <- glmmTMB(
  n ~ lat_group + (1|species) + (1|reference),
  data = trans_lat_groups  |>  filter(trans == "cold_release"),
  family = poisson
)

summary(cold_glm)


warm_glm <- glmmTMB::glmmTMB(
  n ~ lat_group + (1|species) + (1|reference),
  data = trans_lat_groups  |>  filter(trans == "thermal_risk_entry"),
  family = poisson
)

summary(warm_glm)


# 9. Occupied TPC zones ---------------------------------------------------
### a) pooled across latitudes ----
count_zones <- predict_r_shift_lats |> 
  pivot_longer(c(curve_zone_present, curve_zone_future),
               names_to = "time_scenario",
               values_to = "curve_zone")  |> 
  group_by(time_scenario) |> 
  count(curve_zone) |> 
  mutate(time_scenario = if_else(
    time_scenario == "curve_zone_future",
    "future",
    "historical")) |> 
  mutate(sort = case_when(curve_zone == "CE" ~ 1,
                          curve_zone == "CLP" ~ 2,
                          curve_zone == "OPS" ~ 3,
                          curve_zone == "OPD" ~ 4,
                          curve_zone == "HLP" ~ 5,
                          curve_zone == "HE" ~ 6),
         colorado = case_when(curve_zone == "CE" ~ "#B2F2FD",
                              curve_zone == "CLP" ~ "#61DCA9",
                              curve_zone == "OPS" ~ "#8DB63B",
                              curve_zone == "OPD" ~ "#9B7424",
                              curve_zone == "HLP" ~ "#944046",
                              curve_zone == "HE" ~ "#8C0172"))

ggplot(count_zones,
       aes(x = curve_zone,
           y = n))+
  geom_bar(aes(x = fct_reorder(as_factor(curve_zone), sort), 
               y = n, 
               fill =  as_factor(colorado)),
           stat = "identity"
  )+
  scale_fill_identity()+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(title = NULL,
       x = element_blank(),
       y = "#months")+
  theme(plot.title = element_text(face = "italic"))+
  facet_wrap(~time_scenario)+
  theme(panel.spacing = unit(3, "lines"))  # default is usually 0.5 lines

ggsave(here("data/data_sink/figures/supplementary_figs/curvezones.png"),
       width = 1500,
       height = 1500,
       units = "px")

### b) lat groups ----
count_zones_lat <- predict_r_shift_lats |> 
    mutate(group_lat = ifelse(abs(lat) < 35,
                            "low-latitude", 
                            "high-latitude")) |> 
   mutate(trans = case_when(curve_zone_present == "CE" &
                             curve_zone_future == "CLP" ~ "cold_release",
                           curve_zone_present == "CLP" &
                             curve_zone_future == "OPS" ~ "cold_release",
                           curve_zone_present == "CE" &
                             curve_zone_future == "OPS" ~ "cold_release",
                           curve_zone_present == "OPS" &
                             curve_zone_future == "OPD" ~ "thermal_risk_entry",
                           curve_zone_present == "OPS" &
                             curve_zone_future == "HLP" ~ "thermal_risk_entry",
                           curve_zone_present == "CLP" &
                             curve_zone_future == "OPD" ~ "thermal_risk_entry",
                           curve_zone_present == "CLP" &
                             curve_zone_future == "HLP" ~ "thermal_risk_entry",
                           curve_zone_present != "HE" &
                             curve_zone_future == "HE" ~ "thermal_risk_exclusion",
                           .default = NA)) |> 
  drop_na() |> 
  group_by(group_lat) |> 
  count(trans) |> 
  mutate(trans = as_factor(trans)) |> 
  filter(trans != "cold_release")

ggplot(count_zones_lat,
       aes(x = curve_zone,
           y = n))+
  geom_bar(aes(x = trans, 
               y = n), fill = "gray",
           stat = "identity"
  )+
  scale_fill_identity()+
  theme_bw()+
  theme(legend.position = "bottom")+
  labs(title = NULL,
       x = element_blank(),
       y = "#months")+
  theme(plot.title = element_text(face = "italic"))+
  facet_wrap(~group_lat)+
  theme(panel.spacing = unit(3, "lines"))  # default is usually 0.5 lines

ggsave(here("data/data_sink/figures/supplementary_figs/curvezones.png"),
       width = 1500,
       height = 1500,
       units = "px")
# 10. PSI ~ covars -------------------------------------------
### we assigned 0's to r_hist and r_fut whenever they lie below tmin or above
### tmax (line 265 in this script). Negative values 
### of r, which we cannot adequately describe in our modelling framework, 
### negatively affect fitness (Kingsolver et al. 2013). Next, we exclude these
### dual 0's to better assess the analysis of PSI within the TPC range.

### a) psi distribution ----

monthly_r_shift_0s <- read_rds(here("export data/predictions_psi.rds")) |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri") |> 
  mutate(exclude_logic = if_else(preds_present == 0 &
                                   preds_future == 0,
                                 "drops",
                                 "stays")) |> 
  filter(exclude_logic == "stays")

## some comprobations for the sensitivity analysis
length(unique(monthly_r_shift$id_pop)) - length(unique(monthly_r_shift_0s$id_pop)) # only one location is entirely dropped

nrow(monthly_r_shift) - nrow(monthly_r_shift_0s) #345 months discarded (12.45%)

summary(monthly_r_shift_0s$psi) #mean = 0.073, median = 0.095, same range as before

##plot the distribution of psi's
ggplot(monthly_r_shift_0s, aes(y = 1.5)) + 
  ggdist::stat_slab(
    aes(x = psi, fill = after_stat(x)),
    adjust = .5,
    fill_type = "segments",
    color = "gray69",
    linewidth = .5
  )+
  geom_vline(xintercept = mean(monthly_r_shift_0s$psi), 
             linetype = "dashed",
             color = "gray25")+
  coord_cartesian(ylim = c(1.5, NA))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(y = element_blank(),
       x = expression(psi[italic(r)]),
       fill = expression(psi[italic(r)]))+
  scale_fill_gradientn(colours = khroma::colour("vik")(256),
                       values = scales::rescale(c(-0.946331, 0, 0.433069)),
                       limits = c(-0.946331, 0.433069))
ggsave(here("figures/supporting information figures/psi_dist.png"), height = 1500, width = 1500,
       units = "px")

predict_r_shift_0s <- monthly_r_shift_0s |> 
  group_by(reference, species, order, family, id_pop, id_location) |> 
  summarise(psi = mean(psi)) |> 
  mutate(win_or_lose = case_when(psi < 0 ~ as_factor("losers"),
                                 psi >= 0 ~ as_factor("winners")))

predict_r_shift_0s |> 
  ungroup() |> 
  count(win_or_lose) |> 
  mutate(n = n*100/length(unique(predict_r_shift_0s$id_pop)))

predict_r_shift_0s |> 
  ungroup() |> 
  group_by(order) |> 
  summarise(psi = mean(psi))

predict_r_shift_0s_spatial <- monthly_r_shift_0s |> 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, na.fail = FALSE) |> 
  st_transform(crs = "+proj=wintri")

predict_r_shift_preds_0s <- predict_r_shift_0s_spatial |>
  group_by(family) |> 
  summarise(r_hist = mean(preds_present),
            r_fut = mean(preds_future)) |> 
  mutate(dt_hist = log(2)/r_hist,
         dt_fut = log(2)/r_fut) |> 
  mutate(dt_magn = dt_fut - dt_hist)


join_psi_lat_0s <- read_rds(here("export data/predictions_psi.rds")) |>  
  mutate(exclude_logic = if_else(preds_present == 0 &
                                   preds_future == 0,
                                 "drops",
                                 "stays")) |> 
  filter(exclude_logic == "stays")


### b) lmer latitude ----

#### -> no pooling ----
lat_psi_lm_0s <- lmerTest::lmer(psi~ abs(lat) + (1|reference) + (1|species) + (1|month),
                                data = join_psi_lat_0s)
summary(lat_psi_lm_0s)

psi_lat_plot <- sim_and_plot_linears(model_object = lat_psi_lm_0s,
                                             var_x = abs(join_psi_lat_0s$lat),
                                             var_y = join_psi_lat_0s$psi,
                                             n_sims = 1000,
                                             your_title = "Fitness across latitudes",
                                             your_subtitle = NULL,
                                             lab_x = "absolute latitude (º)",
                                             lab_y = expression(psi(r)),
                                             color_points = "#506610",
                                             color_central = "#11220D",
                                             color_uncertainty = "#DAA614",
                                             model_type = "lmer")
print(psi_lat_plot)
ggsave(here("figures/supporting information figures/psi_lat_month_linear.png"), height = 1500, width = 1500,
       units = "px")

ggplot(join_psi_lat_0s, aes(x = lat, y = psi))+
  geom_point(color = "#506610", alpha = .45)+
  geom_smooth(method = "gam",
              color = "#11220D",
              fill = "#DAA614")+
  theme_bw()+
  labs(x = "absolute latitude (º)",
       y = expression(psi[italic(r)]))+
  facet_wrap(~month)

ggsave(here("figures/supporting information figures/psi_lat_months.png"), height = 1500, width = 1500,
       units = "px")

## order
order_psi_lm <- lmerTest::lmer(psi~ as_factor(order) + (1|reference) + (1|species),
                               data = predict_r_shift)
summary(order_psi_lm)
anova(order_psi_lm)


psi_lat_groups  <- read_rds(here("export data/predictions_psi.rds")) |> 
  mutate(lat_group = case_when(abs(lat) <= 23.4366 ~ "tropical",
                               abs(lat) > 23.4366 ~ "temperate")
  ) |> 
  group_by(reference, species, order, family, id_pop, id_location, lat_group) |> 
  summarise(psi = mean(psi),
            lat = mean(lat))

#### ->  pooling ----

psi_lat_02_pool <- join_psi_lat_0s |> 
  group_by(id_pop, species, reference, order) |> 
  summarise(lat = mean(lat),
            psi = mean(psi),
            preds_diff = mean(preds_diff))

lat_psi_lm_0s_pool <- lmerTest::lmer(psi~ abs(lat) + (1|reference) + (1|species),
                                data = psi_lat_02_pool)
summary(lat_psi_lm_0s_pool)

psi_lat_plot_pool <- sim_and_plot_linears(model_object = lat_psi_lm_0s_pool,
                                     var_x = abs(psi_lat_02_pool$lat),
                                     var_y = psi_lat_02_pool$psi,
                                     n_sims = 1000,
                                     your_title = NULL,
                                     your_subtitle = NULL,
                                     lab_x = "absolute latitude (º)",
                                     lab_y = expression(psi[italic(r)]),
                                     color_points = "#506610",
                                     color_central = "#11220D",
                                     color_uncertainty = "#DAA614",
                                     model_type = "lmer")
print(psi_lat_plot_pool)

ggsave(here("figures/supporting information figures/psi_lat_pool.png"), height = 1500, width = 1500,
       units = "px")

ggplot(psi_lat_02_pool, aes(x = lat, y = psi))+
  geom_point(color = "#506610")+
  geom_smooth(method = "gam",
              color = "#11220D",
              fill = "#DAA614")+
  theme_bw()+
  labs(x = "absolute latitude (º)",
       y = expression(psi[italic(r)]))
  

ggsave(here("figures/supporting information figures/psi_lat_pool_smooth.png"), height = 1500, width = 1500,
       units = "px")


## order
order_psi_lm <- lmerTest::lmer(psi~ as_factor(order) + (1|reference) + (1|species),
                               data = join_psi_lat_0s)
summary(order_psi_lm)
anova(order_psi_lm)


psi_lat_groups  <- read_rds(here("export data/predictions_psi.rds")) |> 
  mutate(lat_group = case_when(abs(lat) <= 23.4366 ~ "tropical",
                               abs(lat) > 23.4366 ~ "temperate")
  ) |> 
  group_by(reference, species, order, family, id_pop, id_location, lat_group) |> 
  summarise(psi = mean(psi),
            lat = mean(lat))


lat_psi_lm <- lmerTest::lmer(psi~ abs(lat) + (1|reference) + (1|species),
                             data = join_psi_lat)
summary(lat_psi_lm)

order_psi_lm <- lmerTest::lmer(psi~ as_factor(order) + (1|reference) + (1|species),
                               data = predict_r_shift)
summary(order_psi_lm)
anova(order_psi_lm)


psi_lat_groups  <- read_rds(here("export data/predictions_psi.rds")) |> 
  mutate(lat_group = case_when(abs(lat) <= 23.4366 ~ "tropical",
                               abs(lat) > 23.4366 ~ "temperate")
  ) |> 
  group_by(reference, species, order, family, id_pop, id_location, lat_group) |> 
  summarise(psi = mean(psi),
            lat = mean(lat))

### c) lmer order ----

psi_order_lmer <- lmerTest::lmer(psi ~ as_factor(order) + (1|species) + (1|reference),
                                 data = join_psi_lat_0s)
summary(psi_order_lmer)
anova(psi_order_lmer)
