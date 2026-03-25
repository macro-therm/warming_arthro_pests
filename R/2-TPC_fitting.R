library(here)
library(readxl)
library(MuMIn)
source(here::here("R/S2-TPC_equations.R"))
source(here::here("R/S1-Functions.R"))

# 1. Load -----------------------------------------------------------------
int_rate_data <- read_excel(here("source data/IRTemp.xlsx")) |> 
  mutate(reference = as_factor(reference),
         id_pop = as_factor(id_pop)) |> 
  separate(order, into = c("order", "family"), sep = ": ", remove = TRUE)

# 2. Fitting TPCs  -----------------------------------------------------------------
tpcs_AICs_params <- tibble()

pb <- progress::progress_bar$new(
  format = "Fitting  TPCs [:bar] :percent",
  total = 312,
  clear = F)
for(i in unique(int_rate_data$id_pop)) {
  int_rate_i <-int_rate_data |> 
    filter(id_pop == i) |> 
    filter(!is.na(int_rate))
  species_i <- int_rate_i |> 
    slice(1) |> 
    pull(species)
  reference_i <- int_rate_i |> 
    slice(1) |> 
    pull(reference)
  fitted_tpc_i <- suppressMessages(fit_tpcs(temp = int_rate_i$temperature,
                                            int_rate = int_rate_i$int_rate,
                                            model_name = "all"))
  aics_params_i <- fitted_tpc_i |> 
    group_by(model_name) |> 
    mutate(n_params = n_distinct(param_name)) |> 
    slice(1) |> 
    select(model_name, model_AIC, n_params) |> 
    mutate(id = i)
  tpcs_AICs_params <- bind_rows(tpcs_AICs_params, aics_params_i)
  pb$tick()
  }
tpcs_AICs_params
writexl::write_xlsx(tpcs_AICs_params, here("export data/filter_for_aics.xlsx"))


# 3. Bootstrapping TPCs -----------------------------------------------------------------
simulated_tpcs <- tibble()

for(i in unique(int_rate_data$id_pop)) {
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
  
  ## save predictions for manual model selection later
  possible_error <- tryCatch(expr = {
  fitted_tpc_i <- fit_tpcs(temp = int_rate_i$temperature,
                           int_rate = int_rate_i$int_rate,
                           model_name = "all")
  plot_tpcs(temp = int_rate_i$temperature,
            int_rate = int_rate_i$int_rate,
            fitted_parameters = fitted_tpc_i,
            species = species_i,
            reference = paste0(reference_i,"; ID: ", i))
  ggsave(paste0(here("figures/supplementary_figs/fit_tpcs/id"), i, ".png"),
         width = 2100,
         height = 2100,
         units = "px")
  }, error = function(e) e)
  
  fitted_tpc_equations_i <- unique(fitted_tpc_i$model_name)
  
  ## save bootstrapped predictions for propagating uncertainty later
  possible_error <- tryCatch(expr = {
    sim_tpcs_i <- predict_curves(temp = int_rate_i$temperature,
                               int_rate = int_rate_i$int_rate,
                               fitted_parameters = fitted_tpc_i,
                               model_name_2boot = fitted_tpc_equations_i,n_boots_samples = 100)
  write_rds(sim_tpcs_i, file = paste0(here("data/data_sink/boots_tpc_"), i, ".rds"))
  
  ## plot bootstrapped predictions for model selection later
  plot_uncertainties(bootstrap_uncertainties_tpcs = sim_tpcs_i,
                     temp = int_rate_i$temperature,
                     int_rate = int_rate_i$int_rate,
                     species = species_i,
                     reference = reference_i,
                     pop_id = i)
  ggsave(paste0(here("data/data_sink/figures/supplementary_figs/sim_tpcs/id"), i, ".png"),
         width = 2100,
         height = 2100,
         units = "px")
  }, # <- inside tryCatch
  error = function(e) e)
if (inherits(possible_error, "error")) {
  fit_nls <- NULL
  warning(paste0("Reference ID", i, "(",reference_i, ") could not fit bootstrapped TPCs"))
}
if (is.null(fit_nls)) {
  simulated_tpcs <- simulated_tpcs
} else {
  simulated_tpcs <- bind_rows(simulated_tpcs, sim_tpcs_i) |> 
  drop_na()
  }
  cat(paste0("Ending population ", i, " of 313\n"))
 }

# 4. Fit TPCs for AIC filtering -----------------------------------------------------------------
filtered_tpcs <- readxl::read_excel(here("data/data_sink/filter_for_aics.xlsx")) |> 
  filter(tpc_model != "simp_briere1") |> # discarded due to biological reasoning
  rename(id_pop = id) |> 
  mutate(id_pop = as_factor(id_pop))
joined_int_rate_selected <- left_join(filtered_tpcs, int_rate_data)

tpcs_AICs_params <- tibble()

pb <- progress::progress_bar$new(
  format = "Fitting  TPCs [:bar] :percent",
  total = 312,
  clear = F)

for(i in unique(filtered_tpcs$id_pop)) {
  int_rate_i <-int_rate_data |> 
    filter(id_pop == i) |> 
    filter(!is.na(int_rate))
  species_i <- int_rate_i |> 
    slice(1) |> 
    pull(species)
  reference_i <- int_rate_i |> 
    slice(1) |> 
    pull(reference)
  models_i <- filtered_tpcs |> 
    filter(id_pop == i) |> 
    pull(tpc_model)
  
  fitted_tpc_i <- suppressMessages(fit_tpcs(temp = int_rate_i$temperature,
                                            int_rate = int_rate_i$int_rate,
                                            model_name = models_i))
  aics_params_i <- fitted_tpc_i |> 
    group_by(model_name) |> 
    mutate(n_params = n_distinct(param_name)) |> 
    slice(1) |> 
    select(model_name, model_AIC, n_params) |> 
    mutate(id = i)
  tpcs_AICs_params <- bind_rows(tpcs_AICs_params, aics_params_i)
  pb$tick()
}

#export list of fits with AICs
writexl::write_xlsx(tpcs_AICs_params, here("data/data_sink/tpcs_aics_params.xlsx"))


# 5. Select deltas AICs -----------------------------------------------------------------
int_rate_data_chrs <- int_rate_data |> 
  select(id_pop, reference, species) |> 
  group_by(id_pop) |> 
  slice(1)

tpc_modelselection_aics <- readxl::read_excel(here("export data/aics_for_deltas.xlsx")) |>
  rename(id_pop = id) |> 
  mutate(id_pop = as_factor(id_pop)) |> 
  select(-species, -reference) |> 
  left_join(int_rate_data_chrs, by = "id_pop") |> 
  relocate(c(1, 9, 10, 2:8)) |> 
  mutate(AIC = as.numeric(AIC)) |> 
  group_by(id_pop) |> 
  mutate(
    min_AIC = min(AIC, na.rm = TRUE),
    delta_AIC = AIC - min_AIC,
    within_delta2 = delta_AIC < 2
  ) |> 
  ungroup()  
  
#to help identify models with similar AICs

writexl::write_xlsx(tpc_modelselection_aics,
                    here("data/data_sink/tpcs_selection_filters.xlsx"))



# 6. plot selected TPCs ------------------------------------------------------
tpcs_selected <- readxl::read_excel(here("export data/tpcs_selection_filters_completed.xlsx")) |> 
  filter(step3_uncertainties == "y")

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

tpcs_boots_all

int_rate_max <- int_rate_data |> 
  group_by(id_pop) |> 
  summarise(int_rate_max = max(int_rate, na.rm = TRUE),
            species = unique(species),
            reference = unique(reference))

tpcs_boots_intrate <- inner_join(tpcs_boots_all, int_rate_max) |> 
  filter(preds <= 1.5*int_rate_max) |> 
  mutate(label_pops = paste0("Population ", id_pop)) 

int_rate_sub <- int_rate_data |> 
  filter(id_pop %in% unique(tpcs_boots_intrate$id_pop))

label_vec <- setNames(tpcs_boots_intrate$label_pops, 
                      tpcs_boots_intrate$id_pop)

# Get unique IDs and split into chunks of 24
id_chunks <- split(unique(tpcs_boots_intrate$id_pop),
                   ceiling(seq_along(unique(tpcs_boots_intrate$id_pop)) / 41))


for (i in seq_along(id_chunks)) {
  ids <- id_chunks[[i]]
  tpcs_boots_intrate_sub <- tpcs_boots_intrate|> 
    filter(id_pop %in% ids) |> 
    mutate(id_pop = as.numeric(id_pop))
  int_rate_sub <- int_rate_data |> 
    filter(id_pop %in% ids)
  my_title <- NULL
  life_stage <- NULL
  plot_boot_tpcs <- ggplot2::ggplot() +
    ggplot2::geom_line(data = tpcs_boots_intrate_sub |> 
                         filter(curvetype == "uncertainty"), 
                       aes(x = temp, y = preds, 
                           group = as_factor(boots_iter)), 
                       col = "#0E4D62", linewidth = 0.32, alpha = .16) + 
    ggplot2::geom_line(data = tpcs_boots_intrate_sub |> 
                         filter(curvetype == "estimate"), 
                       aes(x = temp, y = preds), 
                       col = "#CF8143", linewidth = 0.85) + 
    ggplot2::geom_point(data = int_rate_sub, 
                        ggplot2::aes(temperature, int_rate), size = 2)+
    ggplot2::scale_x_continuous(limits = c(0, 50)) +
    ggplot2::theme_bw(base_size = 12) + 
    facet_wrap(~id_pop, scales = "free", 
               labeller = as_labeller(label_vec),
               ncol = 6)+
    ggplot2::labs(x = "Temperature (ºC)", y = italic(r)[m] ~ 
                    (d^-1), 
                  title = my_title, 
                  subtitle = life_stage, 
                  caption = "Bootstrapping with residual resampling, see `rTPC` package vignettes")
  ggsave(
    filename = here(paste0("figures/supporting information figures/tpcs_facet_page_", i, ".png")),
    plot = plot_boot_tpcs,
    width = 8.27, height = 11.69, units = "in", dpi = 300  # A4 size in inches
  )
}



# 8. TPCs supporting ------------------------------------------------------

tpcs_selected_counts <- tpcs_selected |>  
  count(tpc_model, n_params) |> 
  mutate(n_params = if_else(tpc_model %in% c("briere1", 
                                             "simp_beta",
                                             "taylor_sexton",
                                             "analytis_kontodimas"),
                            3, 
                            n_params))

ggplot(tpcs_selected_counts, aes(y = fct_reorder(tpc_model, -n),
                          x = n,
                          fill = as_factor(n_params)))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#CFB15B","#72987F", "#33444B" ))+
  theme_bw()+
  labs(x = "Number of populations in the dataset",
       y = "TPC model equation",
       fill = "Number of parameters")

ggsave(here("data/data_sink/figures/supplementary_figs/tpcs_selected_counts.png"),
       width = 1600,
       height = 1600,
       units = "px")
