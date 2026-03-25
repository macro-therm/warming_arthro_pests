#### script for assemblying the Future CMIP6 models (RCP 8.5, 2041-2060)
## data downloaded from https://geodata.ucdavis.edu/cmip6 for 2.5 resolution

library(tidyverse)
library(terra)
library(readxl)
library(here)
library(sf)

int_rate_data <- read_excel(here("data/data_source/int_rate_dataset_new.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)
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

## a. cmip6_tmax ssp245 ----
# cmip6_tmin
cmip6_models <-  c("ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2", "EC-Earth3-Veg", 
                   "FIO-ESM-2-0", "GISS-E2-1-G", "HadGEM3-GC31-LL", 
                   "INM-CM5-0", "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR",
                   "MRI-ESM2-0", "UKESM1-0-LL")

cmip6_models_tmin_tbl <- tibble(ID = NULL,
                                month = NULL,
                                tmin = NULL,
                                model = NULL)

pb <- progress::progress_bar$new(
  format = "Extracting CMIP6 data for tmin :percent",
  total = length(cmip6_models),
  clear = F)

for(model_cmip6 in cmip6_models){
  raster_model_i <- terra::rast(paste0("E:/cmip6_ssp245/wc2.1_2.5m_tmin_", model_cmip6, "_ssp245_2041-2060.tif"))
  points_tmin_cmip6_i<- terra::extract(raster_model_i,
                                       points_intrates, 
                                       id = points_intrates$id_location) |> 
    as_tibble() |> 
    pivot_longer(cols = -ID,
                 names_to = "month",
                 values_to = "tmin") |> 
    mutate(month = as_factor(str_sub(month, -2)),
           model = as_factor(model_cmip6))
    cmip6_models_tmin_tbl <- bind_rows(cmip6_models_tmin_tbl, points_tmin_cmip6_i)
    pb$tick()
    
}

## b. cmip6_tmax ----
# cmip6_tmax
cmip6_models_tmax_tbl <- tibble(ID = NULL,
                                month = NULL,
                                tmax = NULL,
                                model = NULL)

pb <- progress::progress_bar$new(
  format = "Extracting CMIP6 data for tmax :percent",
  total = length(cmip6_models),
  clear = F)

for(model_cmip6 in cmip6_models){
  raster_model_i <- terra::rast(paste0("E:/cmip6_ssp245/wc2.1_2.5m_tmax_", model_cmip6, "_ssp245_2041-2060.tif"))
  points_tmax_cmip6_i<- terra::extract(raster_model_i,
                                       points_intrates, 
                                       id = points_intrates$id_location) |> 
    as_tibble() |> 
    pivot_longer(cols = -ID,
                 names_to = "month",
                 values_to = "tmax") |> 
    mutate(month = as_factor(str_sub(month, -2)),
           model = as_factor(model_cmip6))
  cmip6_models_tmax_tbl <- bind_rows(cmip6_models_tmax_tbl, points_tmax_cmip6_i)
  pb$tick()
}

## c. cmip6_tavg ----

cmip6_models_tavg <- inner_join(cmip6_models_tmin_tbl, cmip6_models_tmax_tbl) |> 
  mutate(tavg = map2_dbl(.x = tmin, 
                         .y = tmax,
                         .f = ~mean(c(.x, .y)))
  )
saveRDS(object = cmip6_models_tavg,
        file = here("data/data_source/cmip6_tavg_2041-2060_ssp245_res25.rds"))

## e. visualization ----
cmip6_models_tavg_long <- cmip6_models_tavg |> 
  pivot_longer(cols = c(3, 5, 6),
               names_to = "var",
               values_to = "temp")  
  
ggplot(cmip6_models_tavg_long, aes(y = model, x = temp,
                              fill = model))+ 
  ggdist::stat_dist_halfeye(color = "gray24")+
  khroma::scale_fill_batlow(discrete = TRUE)+
  facet_wrap(~var)+
  theme_bw()+
  theme(legend.position = "none")+
  labs(x = "Temperature (ºC)",
       y = "CMIP6 model")

ggsave(here("data/data_sink/figures/tmin_tmax_equation_plot.png"), height = 1900,
       width = 1900,
       units = "px")
