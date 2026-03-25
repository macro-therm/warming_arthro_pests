library(tidyverse)
library(here)
library(readxl)
source(here::here("R/S2-TPC_equations.R"))
source(here::here("R/S1-Functions.R"))

# 1. Load -----------------------------------------------------------------
int_rate_data <- read_excel(here("source data/IRTemp.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)


# 2. Thermal limits for selected TPCs -----------------------------------------------------------------

## a) import selected tpcs ----
tpcs_selected <- readxl::read_excel(here("export data/tpcs_selection_filters_completed.xlsx")) |> 
  filter(step3_uncertainties == "y")

## b) import bootstrapped tpcs ----

tpcs_boots_all <- tibble()
pb <- progress::progress_bar$new(
  format = "Binding Bootstrapped TPCs [:bar] :percent",
  total = 246,
  clear = F)
for(i in unique(tpcs_selected$id_pop)) {
  equation_i <- tpcs_selected |> 
    filter(id_pop == 1) |> 
    pull(tpc_model)
  tpc_boots_i <- read_rds(here(paste0("export data/boots_tpcs/boots_tpc_", i, ".rds"))) |> 
    mutate(id_pop = i) |> 
    filter(model_name_iter == equation_i)
  tpcs_boots_all <- bind_rows(tpcs_boots_all, tpc_boots_i)
  pb$tick()
}

head(tpcs_boots_all)


## b) join bootstraps for each selected tpc ----
thermal_limits_all <- tibble()

for(i in unique(tpcs_boots_all$id_pop)){
  set.seed(2025)
  cat(paste0("Calculating thermal limits for population ", which(unique(tpcs_boots_all$id_pop) == i), " of 246\n" ))
  tpcs_boots_i <- tpcs_boots_all |> 
    filter(id_pop == i)
  int_rate_i <- int_rate_data |> 
    filter(id_pop == i) |> 
    filter(!is.na(int_rate))
  tpc_equation_i <- tpcs_boots_i$model_name_iter[1]
  ## apply our function
  thermal_limits_i <- get_therm_lims(boots_tpc = tpcs_boots_i,
                                     temp = int_rate_i$temperature,
                                     int_rate = int_rate_i$int_rate,
                                     tpc_model = tpc_equation_i,
                                     epsilon = 1e-04) |> 
    mutate(id_pop = i)
  thermal_limits_all <- bind_rows(thermal_limits_all, thermal_limits_i)
  pb$tick
}
write_rds(thermal_limits_all, here("export data/thermal_limits_all.rds"))


