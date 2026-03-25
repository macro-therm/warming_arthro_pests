#remotes::install_github("padpadpadpad/rTPC")
library(tidyverse)
library(rTPC)
library(devRate)
library(here)
source(here::here("R/S2-TPC_equations.R"))



# 1. auxiliary functions for fitting -----------------------------------------------------------------
model_name_translate <- function(user_model_name) {
  if (!all(user_model_name %in% c("all", available_models$model_name))) {
    stop("model name not available. Please check ?available_models")
  }
  
  model_eq <- available_models |>
    dplyr::filter(model_name == user_model_name) |>
    dplyr::select(source_model_name) |>
    dplyr::pull()
  
  return(model_eq)
  
}



# take names from a fitted model to assign them as names for start values later
extract_param_names <- function(nls_object){
  parameter_est <- coef(nls_object)
  param_names <- names(parameter_est)
  return(param_names)
}

## change names of parameters in devRate to easier argument names of parameters
startvals_names_translate_devrate <- function(start_vals, model_name){
  start_vals_names <- if(model_name == "briere1"){
    c("a", "tmin", "tmax")
  } else if (model_name == "lactin1"){
    c("a", "tmax", "delta_t")
  } else if (model_name == "lactin2"){
    c("a", "tmax", "delta_t", "b")
  } else if (model_name == "janisch"){
    c("dmin", "topt", "a", "b")
  } else if (model_name == "linear_campbell"){
    c("intercept", "slope")
  } else if (model_name == "wang"){
    c("k", "r", "topt", "tmin", "tmax", "a")
  } else if (model_name == "mod_polynomial") {
    c("a_0", "a_1", "a_2", "a_3", "a_4")
  } else if (model_name == "ssi") {
    c("p25", "a", "b", "c", "d", "e")
  } else if (model_name == "regniere") {
    c("tmin", "tmax", "phi", "delta_b", "delta_m", "b")
  }  else (NA)
  
  return(start_vals_names)
}

sim_tpc_gridparams <- function(grid_parameters, temperature, model_name){
  params_i <- grid_parameters
  model_i <- model_name
  model_eq <- available_models |>
    dplyr::filter(model_name == model_i)
  tpc_sim_i <- purrr::map(.x = temperature,
                          .f = reformulate(termlabels = unique(model_eq$params_formula))
  )
  tpc_sim_tbl <- dplyr::tibble(temperature,
                               pred_int_rate = tpc_sim_i) |>
    dplyr::mutate(pred_int_rate = unlist(pred_int_rate))
  return(tpc_sim_tbl)
}

start_vals_devRate <- function (model_name_2fit, temperature, dev_rate) {
  check_data(temp = temperature,
             dev_rate)
  
  if (model_name_2fit$model_name == "briere1") {
    start_vals_explore <- c(a = 2e-04, tmin = 8, tmax = 32)
    message("Poorly informative start values for briere1 model")
  } else {
    model_name_devrate <- model_name_2fit$source_model_name
    popdata <- dplyr::tibble(temp = temperature,
                             rate_development = dev_rate)
    start_vals_prev <- devRate::devRateEqStartVal[[model_name_devrate]]
    names(start_vals_prev) <- startvals_names_translate_devrate(start_vals_prev,
                                                                model_name = model_name_2fit$model_name)
    start_upper_vals <- purrr::map(.x = start_vals_prev,
                                   .f = ~.x + abs(.x/2))
    start_lower_vals <- purrr::map(.x = start_vals_prev,
                                   .f = ~.x - abs(.x/2))
    
    multstart_vals_fit <- nls.multstart::nls_multstart(formula = reformulate(response = "rate_development",
                                                                             termlabels = model_name_2fit |>
                                                                               dplyr::pull(formula)),
                                                       data = popdata,
                                                       iter = 500,
                                                       start_lower = start_lower_vals,
                                                       start_upper = start_upper_vals,
                                                       supp_errors = "Y")
    sum_start_vals_fit <- summary(multstart_vals_fit)
    
    if (is.null(multstart_vals_fit)) {
      start_vals_explore <- dplyr::tibble(param_name = names(start_vals_prev),
                                          start_value = unlist(start_vals_prev),
                                          model_name = model_name_2fit$model_name) |>
        dplyr::pull(start_value)
      
      message("generic starting values")
    } else { start_vals_names <- extract_param_names(multstart_vals_fit)
    start_vals <- sum_start_vals_fit$parameters[1:length(start_vals_names), 1]
    start_vals_explore <- dplyr::tibble(param_name = start_vals_names,
                                        start_value = start_vals,
                                        model_name = model_name_2fit$model_name) |>
      dplyr::pull(start_value)
    }
  }
  return(start_vals_explore)
}

check_data <- function(temp = NULL, int_rate = NULL) {
  
  if (any(is.na(int_rate))) {
    stop("intrinsic rate of increase  data have NAs; please consider removing them or fixing them")
  }
  if (any(is.na(temp))) {
    stop("temperature data have NAs; please consider removing them or fixing them")
  }
  if (!is.numeric(temp)) {
    stop("temperature data is not numeric. Please check it.")
  }
  if (!is.numeric(int_rate)) {
    stop("intrinsic rate of increase data is not numeric. Please check it.")
  }
  if (length(temp) != length(int_rate)) {
    stop("intrinsic rate of increase and temperature inputs are not of same length. Please check it.")
  }
  
  if (any(int_rate > 10)) {
    stop("Extremely high values of int_rate intrinsic rate of increase data might contain a typo error. Please check it.")
  }
  
  if (any(temp < -10) | any(temp > 56)) {
    stop("experienced temperatures by active organisms are usually between 0 and 50 degrees centigrades")
  }
  
}


# 2. TPC model fitting -----------------------------------------------------------------

fit_tpcs <- function(temp = NULL,
                     int_rate = NULL,
                     model_name = NULL){
  
  n_params <- AIC <- BIC <- NULL
  
  check_data(temp, int_rate)
  
  if (!is.character(model_name)){
    stop("model_name must be a string in ?available_models")}
  
  if (!all(model_name %in% c("all", available_models$model_name))) {
    stop("model not available. For available model names, see `available_models`")
  }
  
  if (any(model_name == "all")) {
    models_2fit <- available_models |>
      dplyr::filter(n_params <= dplyr::n_distinct(temp)) |>
      dplyr::pull(model_name)
  } else {
    models_2fit <- model_name
  }
  
  list_param <- dplyr::tibble(model_name = NULL,
                              param_name = NULL,
                              start_vals = NULL,
                              param_est = NULL,
                              param_se = NULL,
                              model_AIC = NULL,
                              model_BIC = NULL,
                              model_fit = NULL)
  
  for (i in models_2fit) {
    message(paste0("fitting model ", i)) # to let people know that the function is working and R is not crashing
    model_i <- available_models |>
      dplyr::filter(model_name == i)
    
    if (available_models$package[available_models$model_name == i] == "devRate") {
      
      if (available_models$n_params[available_models$model_name == i] >= length(temp)) {
        fit_nls <- NULL
      } else {
        start_vals <- start_vals_devRate(model_name_2fit = model_i,
                                         temperature = temp,
                                         int_rate = int_rate)
        
        possible_error <- tryCatch(expr = {
          popdata <- dplyr::tibble(temp, int_rate)
          start_upper_vals <- purrr::map(.x = start_vals,
                                         .f = ~.x + abs(.x/2))
          start_lower_vals <- purrr::map(.x = start_vals,
                                         .f = ~.x - abs(.x/2))
          fit_nls <- nls.multstart::nls_multstart(
            formula = stats::reformulate(response = "int_rate",
                                         termlabels = unique(model_i$formula)),
            data = popdata,
            iter = 500,
            start_lower = start_lower_vals,
            start_upper = start_upper_vals,
            supp_errors = "Y")
          
          sum_fit_nls <- summary(fit_nls)
          list_param_tbl <- dplyr::tibble(model_name = i,
                                          param_name = extract_param_names(fit_nls),
                                          start_vals = tidyr::replace_na(start_vals, 0),
                                          param_est = sum_fit_nls$parameters[1:model_i$n_params, 1],
                                          param_se = sum_fit_nls$parameters[1:model_i$n_params, 2],
                                          model_AIC = AIC(fit_nls),
                                          model_BIC = BIC(fit_nls),
                                          model_fit = list(fit_nls))
        }, # <- inside tryCatch
        error = function(e) e)
        if (inherits(possible_error, "error")) {
          fit_nls <- NULL
        }
        if (is.null(fit_nls)) {
          list_param <- list_param
        } else {
          list_param <- list_param |>
            dplyr::bind_rows(list_param_tbl)
        }
      }
    }
    # end of devRate
    
    if (available_models$package[available_models$model_name == i] == "rTPC") {
      possible_error <- tryCatch(expr = {start_vals <- rTPC::get_start_vals(x = temp,
                                                                            y = int_rate,
                                                                            model_name = model_name_translate(i))
      popdata <- dplyr::tibble(temp, int_rate)
      start_upper_vals <- purrr::map(.x = start_vals,
                                     .f = ~.x + abs(.x/2))
      start_lower_vals <- purrr::map(.x = start_vals,
                                     .f = ~.x - abs(.x/2))
      fit_nls <- nls.multstart::nls_multstart(formula = stats::reformulate(response = "int_rate",
                                                                           termlabels = unique(model_i$formula)),
                                              data = popdata,
                                              iter = 500,
                                              start_lower = start_lower_vals,
                                              start_upper = start_upper_vals,
                                              lower = rTPC::get_lower_lims(popdata$temp,
                                                                           popdata$int_rate,
                                                                           model_name = model_name_translate(i)),
                                              upper = rTPC::get_upper_lims(popdata$temp,
                                                                           popdata$int_rate,
                                                                           model_name = model_name_translate(i)),
                                              supp_errors = "Y")
      sum_fit_nls <- summary(fit_nls)
      list_param_tbl <- dplyr::tibble(model_name = i,
                                      param_name = extract_param_names(fit_nls),
                                      start_vals = tidyr::replace_na(start_vals, 0),
                                      param_est = sum_fit_nls$parameters[1:model_i$n_params, 1],
                                      param_se = sum_fit_nls$parameters[1:model_i$n_params, 2],
                                      model_AIC = AIC(fit_nls),
                                      model_BIC = BIC(fit_nls),
                                      model_fit = list(fit_nls))
      }, # <- inside tryCatch
      error = function(e) e)
      if (inherits(possible_error, "error")) {
        fit_nls <- NULL
      }
      if (is.null(fit_nls)) {
        list_param <- list_param
      } else {list_param <- list_param |>
        dplyr::bind_rows(list_param_tbl)}
    }
    # end of rTPC processing
    
  } # <- loop ends
  if (length(list_param) == 0) {
    warning("no model converged adequately for fitting your data")
  } else {
    list_param <-   list_param |>
      dplyr::group_by(model_name) |>
      dplyr::mutate(model_fit = dplyr::if_else(dplyr::row_number() == 1,
                                               model_fit,
                                               list(NULL)))  |>
      dplyr::ungroup()
  }
  return(list_param)
}


# 3. Explore fitted TPCs  -----------------------------------------------------------------

plot_tpcs <- function(temp = NULL,
                           int_rate = NULL,
                           fitted_parameters = NULL,
                           species = NULL,
                           reference = NULL) {
  
  check_data(temp, int_rate)

  if (is.null(fitted_parameters)) {
    stop("`fitted_parameters` is NULL; use `mappestRisk::fit_tpcs()` to check that at least one model converged")
  }
  
  if (!inherits(fitted_parameters, "data.frame")) {
    stop("The argument `fitted_parameters` must be a tibble or data.frame as produced by `mappestRisk::fit_tpcs()` function. No modifications of columns of the fitted_parameters are allowed, but you can subset observations by filtering or subsetting by rows if desired.")
  }
  
  if (suppressWarnings(any(!c("param_name", "start_vals", "param_est",
                              "param_se", "model_name", "model_AIC",
                              "model_BIC") %in% names(fitted_parameters)))) {
    stop("The argument `fitted_parameters` must be a tibble or data.frame as produced by `mappestRisk::fit_tpcs()` function. No modifications of columns of the fitted_parameters are allowed, but you can subset observations by filtering or subsetting by rows if desired.")
  }
  
  if (nrow(fitted_parameters) == 0) {
    stop("no model has converged in your `fitted_parameters` data.frame. Is it the appropriate object coming from converged `fit_tpcs()`?")
  }
  
  if (typeof(species) != "character" && !is.null(species)) {
    stop("`species` must be a character or NULL")
  }
  
  if (typeof(reference) != "character" && !is.null(reference)) {
    stop("`reference` must be a character or NULL")
  }
  
  popdata <- dplyr::tibble(temperature = temp,
                           intrinsic_rate = int_rate)
  
  predict2fill <- dplyr::tibble(temp = NULL,
                                int_rate = NULL,
                                model_name = NULL,
                                model_AIC = NULL)
  
  model_names2plot <- fitted_parameters |>
    dplyr::distinct(model_name) |>
    dplyr::pull(model_name)
  
  for (i in model_names2plot) {
    fitted_parameters_i <- fitted_parameters |>
      dplyr::filter(model_name == i)
    model_AIC_i <- fitted_parameters_i |>
      dplyr::pull(model_AIC)
    params_i <- fitted_parameters_i |>
      dplyr::pull(param_est)
    formula_i <- available_models |>
      dplyr::filter(model_name == i) |>
      dplyr::pull(params_formula)
    ##predict based on parameters
    explore_preds <- dplyr::tibble(temp = seq(min(popdata$temperature) - 15,
                                              max(popdata$temperature) + 15,
                                              .01),
                                   model_name = i,
                                   model_AIC = model_AIC_i[1],
                                   preds = NULL,
                                   n_params = length(params_i))
    fit_vals_tbl <- explore_preds |>
      dplyr::select(temp, model_name, model_AIC, n_params) |>
      dplyr::mutate(formula = formula_i) |>
      dplyr::mutate(preds = purrr::map_dbl(.x = temp,
                                           .f = stats::reformulate(unique(formula_i)))) |>
      dplyr::filter(preds >= 0) |>
      dplyr::select(-formula) |>
      dplyr::mutate(preds = dplyr::case_when(model_name == "ratkowsky" & temp > params_i[2] ~ NA_real_,
                                             model_name == "ratkowsky" & temp < params_i[1] ~ NA_real_,
                                             model_name == "briere1" & temp < params_i[1] ~ NA_real_,
                                             model_name == "briere2" & temp < params_i[1] ~ NA_real_,
                                             TRUE ~ preds)
      ) # to exclude biological non-sense predictions due to model mathematical properties
    predict2fill <- predict2fill |>
      dplyr::bind_rows(fit_vals_tbl)
  }
  
  aic_text <-  predict2fill  |>
    dplyr::group_by(model_name)  |>
    dplyr::summarise(aic = mean(model_AIC),
                     n_params = paste(mean(n_params), "parameters"))  |>
    dplyr::arrange(aic)
  aic_order <- aic_text  |>
    dplyr::pull(model_name)
  aic_values <- aic_text |>
    dplyr::mutate(aic =   paste("AIC =", round(aic, 2)),
                  temp = min(popdata$temperature),
                  preds = 1.5*max(popdata$intrinsic_rate))
  my_title <- substitute(italic(paste(especie)), list(especie = species))
  ggplot_models <- ggplot2::ggplot() +
    ggplot2::geom_line(data = predict2fill |>
                         dplyr::filter(preds < (1.5*max(popdata$intrinsic_rate))),
                       ggplot2::aes(x = temp, y = preds),
                       color = "#CF8143",
                       linewidth = 1.3) +
    ggplot2::geom_point(data = popdata, ggplot2::aes(x = temperature,
                                                     y = intrinsic_rate),
                        color = "#0E4D62",
                        alpha = .8,
                        size = 1.5) +
    ggplot2::facet_wrap(~factor(model_name)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = "Temperature",
                  y = expression(italic(R(T))~(d^-1)),
                  title = my_title,
                  subtitle = reference) +
    ggplot2::geom_label(data = aic_values,
                        ggplot2::aes(label = aic,
                                     x = temp,
                                     y = preds),
                                     fill = "#0E4D62",
                        color = "white",
                        size = 3) +
    ggplot2::geom_label(data = aic_values,
                        ggplot2::aes(x = temp,
                                     y = preds-preds/8,
                                     label = n_params),
                                     fill = "#0E4D62",
                        color = "white",
                        size = 3)
  
  return(ggplot_models)
}

# 4. Bootstrapping TPCs for parameter uncertainty propagation  -----------------------------------------------------------------
predict_curves <- function(temp = NULL,
                           int_rate = NULL,
                           fitted_parameters = NULL,
                           model_name_2boot = NULL,
                           propagate_uncertainty = TRUE,
                           n_boots_samples = 100) {
  check_data(temp, int_rate)
  if(is.null(fitted_parameters)) {
    stop("`fitted_parameters` must be provided.")
  }
  
  if (!all(model_name_2boot %in% fitted_parameters$model_name)) {
    message(paste0("Models available: ", paste0(unique(fitted_parameters$model_name), collapse = ", ")))
    stop("model not available. Check the models that converged in `fitted_parameters`")
  }
  

  if (n_boots_samples < 100){
    warning("100 iterations might be desirable. Consider increasing `n_boots_samples` if possible")
  }
  
  if (!is.logical(propagate_uncertainty)) {
    stop("`propagate_uncertainty` must be `TRUE` or `FALSE` (def. `TRUE`)")
  }
  
  popdata <- dplyr::tibble(temp,
                           int_rate)
  predict2fill <- dplyr::tibble(temp = NULL,
                                int_rate = NULL,
                                model_name = NULL,
                                model_fit = NULL,
                                model_AIC = NULL)
  
  for (model_name_i in model_name_2boot){
    fitted_parameters_i <- fitted_parameters |>
      dplyr::filter(model_name == model_name_i)
    model_AIC_i <- fitted_parameters_i |>
      dplyr::pull(model_AIC)
    model_fit_i <- fitted_parameters_i |>
      dplyr::pull(model_fit)
    params_i <- fitted_parameters_i |>
      dplyr::pull(param_est)
    formula_i <- available_models |>
      dplyr::filter(model_name == model_name_i) |>
      dplyr::pull(params_formula)
    ##predict based on parameters
    explore_preds <- dplyr::tibble(temp = seq(min(popdata$temp) - 20,
                                              max(popdata$temp) + 15,
                                              .1),
                                   model_name = model_name_i,
                                   model_fit = model_fit_i[1],
                                   model_AIC = model_AIC_i[1],
                                   preds = NULL,
                                   n_params = length(params_i))
    fit_vals_tbl <- explore_preds |>
      dplyr::select(temp, model_name, model_AIC, n_params, model_fit) |>
      dplyr::mutate(formula = formula_i) |>
      dplyr::mutate(preds = purrr::map_dbl(.x = temp,
                                           .f = stats::reformulate(unique(formula_i)))) |>
      dplyr::filter(preds >= -0.2) |>
      dplyr::select(-formula)
    # to exclude biological non-sense predictions due to model mathematical properties
    predict2fill <- predict2fill |>
      dplyr::bind_rows(fit_vals_tbl)
  }
  predict2fill_complete <- predict2fill |>
    dplyr::mutate(tidyr::nest(popdata)) |>
    dplyr::mutate(coefs = purrr::map(.x = model_fit,
                                     .f = stats::coef))
  if(nrow(predict2fill_complete) == 0) {stop("No bootstrap was attempted. Check your model(s)")}
  if (propagate_uncertainty == FALSE) {
    tpc_estimate <- dplyr::tibble(model_name = predict2fill_complete$model_name,
                                  iter = rep(NA, nrow(predict2fill_complete)),
                                  temp = predict2fill_complete$temp,
                                  pred = predict2fill_complete$preds,
                                  curvetype = rep("estimate", nrow(predict2fill_complete))
    )
    warning("No bootstrap was performed. We strongly recommend to propagate uncertainty.")
    return(tpc_estimate)
    
  } else {
    cat("\nADVISE: the simulation of new bootstrapped curves takes some time. Await patiently or reduce your `n_boots_samples`\n")
    sim_boots_tpcs <- dplyr::tibble()
    for (model_i in model_name_2boot){
      predict_model_i <- predict2fill_complete |>
        dplyr::filter(model_name == model_i) |>
        dplyr::filter(preds >= 0)
      coefs_i <- purrr::as_vector(unique(predict_model_i$coefs))
      temp_data_i <- popdata
      formula_i <- available_models |>
        dplyr::filter(model_name == model_i) |>
        dplyr::pull(params_formula)
      model_fit_i <- predict_model_i$model_fit[[1]]
      # extract the residuals and the fitted values
      resids_i <- residuals(model_fit_i)
      fit_vals_i <- fitted(model_fit_i)
      
      ## residual resampling
      resampled_data_resid <- dplyr::tibble()
      pb <- progress::progress_bar$new(
        format = paste0(model_i,": Predicting bootstrapped TPCs [:bar] :percent"),
        total = n_boots_samples,
        clear = F)
      
      for(n_boot_sample in 1:n_boots_samples) {
        resampled_resids_i <- sample(resids_i,
                                     size = length(resids_i),
                                     replace = TRUE)
        resampled_obs_i <- fit_vals_i + resampled_resids_i
        resampled_data_i <- dplyr::tibble(
          predict_var = popdata$temp,
          response_var = resampled_obs_i,
          boot_sample_id = n_boot_sample,
          model_name = model_i
        )
        resampled_data_resid <- dplyr::bind_rows(resampled_data_resid, resampled_data_i)
      }
      # Then fit TPC to each bootstrapped iterated with residual resampling data set and get predictions similarly as before
      predicted_boots_resid <- dplyr::tibble()
      
      for (iter in 1:n_boots_samples) {
        resid_resampled_i <- resampled_data_resid |>
          dplyr::filter(boot_sample_id == iter) |>
          dplyr::filter(response_var >= 0)
        
        resid_fitted_tpc_iter <- suppressMessages(
          fit_tpcs(temp = resid_resampled_i$predict_var,
                        int_rate = resid_resampled_i$response_var,
                        model_name = model_i))
        
        if (nrow(resid_fitted_tpc_iter) == 0) {
          next  # skip this iteration
        }
        
        resid_predictions_temp_iter <- seq(min(popdata$temp, na.rm = TRUE) - 20,
                                           max(popdata$temp, na.rm = TRUE) + 15,
                                           0.1)
        model_fit_boot_iter <- resid_fitted_tpc_iter$model_fit[[1]]
        params_i <- stats::coef(model_fit_boot_iter)
        resid_predictive_tbl_iter <- dplyr::tibble(
          resid_predictions_temp_iter,
          preds_rate =  purrr::map_dbl(.x = resid_predictions_temp_iter,
                                       .f = stats::reformulate(formula_i))
        ) |>
          dplyr::filter(preds_rate >= 0) |>
          dplyr::mutate(boots_iter = iter) |>
          dplyr::mutate(model_name_iter = model_i) |>
          dplyr::mutate(curvetype = "uncertainty") |>
          dplyr::rename(temp = resid_predictions_temp_iter,
                        preds = preds_rate)
        
        predicted_boots_resid <- dplyr::bind_rows(predicted_boots_resid, resid_predictive_tbl_iter)
        pb$tick()
      } # <- bootstrapped simulations for each model equation
      cat(paste0("\n Bootstrapping simulations completed for ", model_i, " \n"))
      estimated_tpc_i <- predict_model_i |>
        dplyr::select(temp, preds, model_name) |>
        dplyr::mutate(boots_iter = NULL,
                      curvetype = "estimate") |>
        dplyr::rename(model_name_iter = model_name)
      
      predicted_boots_resid <- dplyr::bind_rows(predicted_boots_resid, estimated_tpc_i)
      
      ## end of loop for each model TPC simulations
      sim_boots_tpcs <- dplyr::bind_rows(sim_boots_tpcs, predicted_boots_resid)
    } ## end of model equation simulation loop
    if(!any(sim_boots_tpcs$curvetype == "uncertainty")){
      warning("No bootstrap was accomplished. Your model might not be suitable for bootstrapping
due to convergence problems")
    }
  } # end of condition for `uncertainty == TRUE`
  return(sim_boots_tpcs)
}


plot_uncertainties <- function(temp = NULL, 
                               int_rate = NULL, 
                               bootstrap_tpcs = NULL,
                               species = NULL, 
                               life_stage = NULL, 
                               alpha = 0.2) 
{

  devdata <- dplyr::tibble(temp, int_rate)
  central_curve <- dplyr::filter(bootstrap_tpcs, curvetype == 
                                   "estimate")
  uncertainty_curves <- dplyr::filter(bootstrap_tpcs, curvetype == 
                                        "uncertainty")
  my_title <- substitute(italic(paste(x)), list(x = species))
  plot_boot_tpcs <- ggplot2::ggplot() +
    ggplot2::geom_line(data = uncertainty_curves, 
                       aes(x = temp, y = preds, 
                           group = as_factor(boots_iter)), 
                       col = "#0E4D62", linewidth = 0.32) + 
    ggplot2::geom_line(data = central_curve, 
                       aes(x = temp, y = preds), 
                       col = "#CF8143", linewidth = 0.85) + 
    ggplot2::geom_point(data = devdata, ggplot2::aes(temp, 
                                                     int_rate), size = 2)+
    ggplot2::scale_x_continuous(limits = c(0, 50)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.5 *max(devdata$int_rate, na.rm = TRUE))) + 
    ggplot2::theme_bw(base_size = 12) + 
    facet_wrap(~id_pop, scales = "free_y")+
    ggplot2::labs(x = "Temperature", y = italic(r)[m](T) ~ 
                    (d^-1), title = my_title, subtitle = life_stage, 
                  caption = "Bootstrapping with residual resampling, see `rTPC` package vignettes")
  return(plot_boot_tpcs)
  }



# 6. Calculate Thermal Limits  -----------------------------------------------------------------
get_therm_lims <- function(boots_tpc,
                           temp,
                           int_rate,
                           tpc_model,
                           epsilon = 1e-4
                           ) {

  # Find where the predicted values cross zero (or near-zero, defined by epsilon)
  fitted_tpc_i <- suppressMessages(
    fit_tpcs(temp = int_rate_i$temperature,
             int_rate = int_rate_i$int_rate,
             model_name = tpc_model))
  
  topt_i <- as.numeric(rTPC::calc_params(model = fitted_tpc_i$model_fit[[1]])[2])
  rmax_i <- as.numeric(rTPC::calc_params(model = fitted_tpc_i$model_fit[[1]])[1])
  ##left-side
  therm_lims <- tibble()
  for (iter in unique(boots_tpc$boots_iter)) {

    if (is.na(iter)) {
      boots_near_zero <- boots_tpc |> 
        filter(is.na(boots_iter))|> 
        mutate(is_near_zero = ifelse(abs(preds <= epsilon),
                                     TRUE,
                                     FALSE)) 
    } else {
      boots_near_zero <- boots_tpc |> 
      filter(boots_iter == iter) |> 
        mutate(is_near_zero = ifelse(abs(preds <= epsilon),
                                     TRUE,
                                     FALSE)) 
    }

    boots_tpc_left <- boots_near_zero |> 
      filter(temp < topt_i) |> 
      filter(abs(preds) < 1)
    if(!any(boots_tpc_left$is_near_zero)) {
      therm_lim_minimum <- min(boots_tpc_left$temp, na.rm = TRUE)
      left_semimaximum_index <- max(which(boots_tpc_left$preds <= (rmax_i/2)))
      left_semimaximum <- boots_tpc_left$temp[left_semimaximum_index]
    } else {
      therm_lim_minimum_index <- max(which(boots_tpc_left$is_near_zero), na.rm = TRUE)
      therm_lim_minimum <- boots_tpc_left$temp[therm_lim_minimum_index]
      
      left_semimaximum_index <- max(which(boots_tpc_left$preds <= (rmax_i/2)))
      left_semimaximum <- boots_tpc_left$temp[left_semimaximum_index]
    }
    
    ##right-side
    boots_tpc_right <- boots_near_zero |> 
      filter(temp >= topt_i) |> 
      filter(abs(preds) < 1)
    if(!any(boots_tpc_right$is_near_zero)) {
      therm_lim_maximum <- max(boots_tpc_right$temp, na.rm = TRUE)
      right_semimaximum_index <- min(which(boots_tpc_right$preds <= (rmax_i/2)), na.rm = TRUE)
      right_semimaximum <- boots_tpc_right$temp[right_semimaximum_index]
    } else {
      therm_lim_maximum_index <- min(which(boots_tpc_right$is_near_zero), na.rm = TRUE)
      therm_lim_maximum <- boots_tpc_right$temp[therm_lim_maximum_index]
      right_semimaximum_index <- min(which(boots_tpc_right$preds <= (rmax_i/2)), na.rm = TRUE)
      right_semimaximum <- boots_tpc_right$temp[right_semimaximum_index]
      
    }
    if(therm_lim_minimum < 0 || is.na(therm_lim_minimum)) { therm_lim_minimum <- NA}
    if(therm_lim_maximum > 42 || #2 degrees above the maximum measured positive intrinsic rate in our data set 
       is.na(therm_lim_maximum)) {therm_lim_maximum <- NA} 
    
    valid_range <- boots_near_zero |>
      filter(temp >= therm_lim_minimum,
             temp <= therm_lim_maximum)
    if (nrow(valid_range) == 0) {
      topt_iter <- NA
      r_max_iter <- NA
    } else {
      topt_iter <- valid_range$temp[which.max(valid_range$preds)]
      r_max_iter <-  max(valid_range$preds, na.rm = TRUE)
    }
    therm_lims_i <- c(therm_lim_minimum,
                      therm_lim_maximum,
                      topt_iter,
                      left_semimaximum,
                      right_semimaximum,  
                      r_max_iter) 
    names(therm_lims_i) <- c("tmin", "tmax", "topt", "t_L50", "t_R50", "rmax")
    therm_lims <- bind_rows(therm_lims, therm_lims_i)
  }
  return(therm_lims)
}


# 7. Sim&Plot Hierarchical Models  -----------------------------------------------------------------

sim_and_plot_linears <- function(model_object,
                                 var_x,
                                 var_y,
                                 n_sims,
                                 your_title,
                                 your_subtitle,
                                 lab_x,
                                 lab_y,
                                 color_points,
                                 color_central,
                                 color_uncertainty,
                                 model_type) {
  
  modelled_relation <- tibble(indep_var = var_x, 
                              dep_var = var_y)
  sum_model <- summary(model_object)
  if(model_type == "lmer") {
    model_params <- tibble(intercept = sum_model$coefficients[1,1],
                           slope =  sum_model$coefficients[2,1],
                           intercept_se = sum_model$coefficients[1,2],
                           slope_se =  sum_model$coefficients[2,2],
                           sample_size = nrow(modelled_relation)) |> 
      mutate(intercept_sd = intercept_se * sqrt(sample_size),
             slope_sd = slope_se * sqrt(sample_size))  
  } else if (model_type == "metafor") {
    model_params <- tibble(intercept = sum_model$beta[1,1],
                           slope =  sum_model$beta[2,1],
                           intercept_se = sum_model$se[1],
                           slope_se =  sum_model$se[2],
                           sample_size = nrow(modelled_relation))
  }
  set.seed(2023)
  sim_slope_model <- rnorm(n_sims, 
                           mean = model_params |> select(slope) |> as_vector(),
                           sd = model_params |> select(slope_se) |> as_vector())
  
  sim_intercept_model <- rnorm(n_sims, 
                               mean = model_params |> select(intercept) |> as_vector(),
                               sd = model_params |> select(intercept_se) |> as_vector())
  
  sim_model_lm <- tibble(slope = sim_slope_model,
                         intercept = sim_intercept_model)
  
  n_draws <- n_sims
  alpha_level <- .15
  col_draw <- color_uncertainty
  col_median <-  color_central
  
  model_plot <- ggplot(modelled_relation,
                       aes(x = indep_var, y = dep_var))+
    labs(title= your_title,
         subtitle= your_subtitle,
         x = lab_x,
         y = lab_y)+
    ggthemes::theme_few()+ 
    geom_abline(aes(slope = slope,
                    intercept =intercept),
                data = slice_sample(sim_model_lm, n = n_draws),
                color = col_draw,
                alpha = alpha_level)+
    geom_abline(aes(slope = model_params$slope,
                    intercept = model_params$intercept),
                color = color_central,
                linewidth = 1.4)+
    geom_point(alpha=0.82,
               color= color_points)+
    theme(plot.title = element_text(face = "bold"))
  
  model_plot
  
}
