
model_names <- c("analytis_kontodimas", "ashrafi2", "atkin",
                 "beta", "simp_beta", "briere1", 
                 "simp_briere1", "deutsch",
                 "johnk", "pawar", "rezende", "ssi",
                 "taylor_sexton", "thomas",
                 "tomlinson_phillips", "weibull")

pkg_model_names <- c("analytiskontodimas_2004", "ashrafi2_2018", "atkin_2005",
                     "beta_2012", "betatypesimplified_2008", "briere1_1999", 
                     "briere1simplified_1999", "deutsch_2008",
                     "joehnk_2008", "pawar_2018", "rezende_2019", "sharpeschoolfull_1981",
                     "taylorsexton_1972", "thomas_2012",
                     "tomlinsonphillips_2015", "weibull_1995")

available_models <- dplyr::tibble(model_name = model_names) |>
  dplyr::mutate(package = "rTPC",
                source_model_name = pkg_model_names) |>
  dplyr::mutate(formula = dplyr::case_when(
    model_name == "analytis_kontodimas" ~ "rTPC::analytiskontodimas_2004(temp, a, tmin, tmax)",
    model_name == "ashrafi2" ~ "rTPC::ashrafi2_2018(temp, a, b, c)",
    model_name == "atkin" ~ "rTPC::atkin_2005(temp, r0, a, b)",
    model_name == "beta" ~ "rTPC::beta_2012(temp, a, b, c, d, e)",
    model_name == "simp_beta" ~ "rTPC::betatypesimplified_2008(temp, rho, alpha, beta)",
    model_name == "briere1" ~ "rTPC::briere1_1999(temp, tmin, tmax, a)",
    model_name == "simp_briere1" ~ "rTPC::briere1simplified_1999(temp, tmin, tmax, a)",
    model_name == "deutsch" ~ "rTPC::deutsch_2008(temp, rmax, topt, ctmax, a)",
    model_name == "johnk" ~ "rTPC::joehnk_2008(temp, rmax, topt, a, b, c)",
    model_name == "pawar" ~ "rTPC::pawar_2018(temp, r_tref, e, eh, topt, tref = 25)",
    model_name == "rezende" ~ "rTPC::rezende_2019(temp, q10, a, b, c)",
    model_name == "ssi" ~ "rTPC::sharpeschoolfull_1981(temp, r_tref, e, el, tl, eh, th, tref)",
    model_name == "taylor_sexton" ~ "rTPC::taylorsexton_1972(temp, rmax, tmin, topt)",
    model_name == "thomas" ~ "rTPC::thomas_2012(temp, a, b, c, t_ref)",
    model_name == "tomlinson_phillips" ~ "rTPC::tomlinsonphillips_2015(temp, a, b, c)",
    model_name == "weibull" ~ "rTPC::weibull_1995(temp, a, topt, b, c)")
    ) |>
  dplyr::mutate(params_formula = dplyr::case_when(
    model_name == "analytis_kontodimas" ~ "rTPC::analytiskontodimas_2004(.x, params_i[1], params_i[2], params_i[3])",
    model_name == "ashrafi2" ~ "rTPC::ashrafi2_2018(.x, params_i[1], params_i[2], params_i[3])",
    model_name == "atkin" ~ "rTPC::atkin_2005(.x, params_i[1], params_i[2], params_i[3])",
    model_name == "beta" ~ "rTPC::beta_2012(.x, params_i[1], params_i[2], params_i[3], params_i[4], params_i[5])",
    model_name == "simp_beta" ~ "rTPC::betatypesimplified_2008(.x, params_i[1], params_i[2], params_i[3])",
    model_name == "briere1" ~ "rTPC::briere1_1999(.x, params_i[1], params_i[2], params_i[3])",
    model_name == "simp_briere1" ~ "rTPC::briere1simplified_1999(.x, params_i[1], params_i[2], params_i[3])",
    model_name == "deutsch" ~ "rTPC::deutsch_2008(.x, params_i[1], params_i[2], params_i[3], params_i[4])",
    model_name == "johnk" ~ "rTPC::joehnk_2008(.x, params_i[1], params_i[2], params_i[3], params_i[4], params_i[5])",
    model_name == "pawar" ~ "rTPC::pawar_2018(.x, params_i[1], params_i[2], params_i[3], params_i[4], tref = 25)",
    model_name == "rezende" ~ "rTPC::rezende_2019(.x, params_i[1], params_i[2], params_i[3], params_i[4])",
    model_name == "ssi" ~ "rTPC::sharpeschoolfull_1981(.x, params_i[1], params_i[2], params_i[3], params_i[4], params_i[5], params_i[6], tref = 25)",
    model_name == "taylor_sexton" ~ "rTPC::taylorsexton_1972(.x, params_i[1], params_i[2], params_i[3])",
    model_name == "thomas" ~ "rTPC::thomas_2012(.x, params_i[1], params_i[2], params_i[3], params_i[4])",
    model_name == "tomlinson_phillips" ~ "rTPC::tomlinsonphillips_2015(.x, params_i[1], params_i[2], params_i[3])",
    model_name == "weibull" ~ "rTPC::weibull_1995(.x, params_i[1], params_i[2], params_i[3], params_i[4])"
    )) |>
  dplyr::mutate(working_formula = dplyr::case_when(
    model_name == "analytis_kontodimas" ~ "rTPC::analytiskontodimas_2004(.x, start_vals[1], start_vals[2], start_vals[3])",
    model_name == "ashrafi2" ~ "rTPC::ashrafi2_2018(.x, start_vals[1], start_vals[2], start_vals[3])",
    model_name == "atkin" ~ "rTPC::atkin_2005(.x, start_vals[1], start_vals[2], start_vals[3])",
    model_name == "beta" ~ "rTPC::beta_2012(.x, start_vals[1], start_vals[2], start_vals[3], start_vals[4], start_vals[5])",
    model_name == "simp_beta" ~ "rTPC::betatypesimplified_2008(.x, start_vals[1], start_vals[2], start_vals[3])",
    model_name == "briere1" ~ "rTPC::briere1_1999(.x, start_vals[1], start_vals[2], start_vals[3])",
    model_name == "simp_briere1" ~ "rTPC::briere1simplified_1999(.x, start_vals[1], start_vals[2], start_vals[3])",
    model_name == "deutsch" ~ "rTPC::deutsch_2008(.x, start_vals[1], start_vals[2], start_vals[3], start_vals[4])",
    model_name == "johnk" ~ "rTPC::joehnk_2008(.x, start_vals[1], start_vals[2], start_vals[3], start_vals[4], start_vals[5])",
    model_name == "pawar" ~ "rTPC::pawar_2018(.x, start_vals[1], start_vals[2], start_vals[3], start_vals[4], tref = 25)",
    model_name == "rezende" ~ "rTPC::rezende_2019(.x, start_vals[1], start_vals[2], start_vals[3], start_vals[4])",
    model_name == "ssi" ~ "rTPC::sharpeschoolfull_1981(.x, start_vals[1], start_vals[2], start_vals[3], start_vals[4], start_vals[5], start_vals[6], tref = 25)",
    model_name == "taylor_sexton" ~ "rTPC::taylorsexton_1972(.x, start_vals[1], start_vals[2], start_vals[3])",
    model_name == "thomas" ~ "rTPC::thomas_2012(.x, start_vals[1], start_vals[2], start_vals[3], start_vals[4])",
    model_name == "tomlinson_phillips" ~ "rTPC::tomlinsonphillips_2015(.x, start_vals[1], start_vals[2], start_vals[3])",
    model_name == "weibull" ~ "rTPC::weibull_1995(.x, start_vals[1], start_vals[2], start_vals[3], start_vals[4])"
  )) |>
  dplyr::mutate(n_params = dplyr::case_when(
    model_name == "analytis_kontodimas" ~ 3,
    model_name == "ashrafi2" ~ 3,
    model_name == "atkin" ~ 3,
    model_name == "beta" ~ 5,
    model_name == "simp_beta" ~ 3,
    model_name == "briere1" ~ 3,
    model_name == "simp_briere1" ~ 3,
    model_name == "deutsch" ~ 4,
    model_name == "johnk" ~ 5,
    model_name == "pawar" ~ 4,
    model_name == "rezende" ~ 4,
    model_name == "ssi" ~ 6,
    model_name == "taylor_sexton" ~ 3,
    model_name == "thomas" ~ 4,
    model_name == "tomlinson_phillips" ~ 3,
    model_name == "weibull" ~ 4
    ))

