load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

pkgs <- c(
  "tidyverse",
  "brms",
  "posterior",
  "patchwork",
  "ggplot2"
)

invisible(lapply(pkgs, load_or_install))


# 1) POWER

##  simulate basket
simulate_basket <- function(id, n, mu_treat, mu_control = 0, sd = 1) {
  y <- c(
    rnorm(n, mu_control, sd),
    rnorm(n, mu_treat,   sd)
  )
  y_scaled <- (y - mu_control) / sd
  tibble(
    basket    = factor(id),
    treatment = rep(c(0, 1), each = n),
    y         = y_scaled
  )
}


# fit model for heterogeneity 

fit_basket_model <- function(data, heterogeneity = c("high", "low", "moderate"),
                             iter = 5000, chains = 4, seed = 123) {
  heterogeneity <- match.arg(heterogeneity)
  
  cauchy_scale <- switch(heterogeneity,
                         high     = 0.8,
                         low      = 0.2,
                         moderate = 0.4
  )
  
  priors <- c(
    prior(student_t(3, 0, 5), class = Intercept),
    prior(normal(0, 1),        class = b),
    set_prior(paste0("cauchy(0, ", cauchy_scale, ")"),
              class = "sd", group = "basket", coef = "treatment"),
    prior(cauchy(0, 0.5), class = sigma)
  )
  
  brm(
    y ~ treatment + (1 | basket) + (0 + treatment | basket),
    data = data, prior = priors,
    iter = iter, warmup = iter / 2, chains = chains,
    seed = seed, cores = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 18),
    refresh = 0
  )
}

## Extract decisions for MCMC

extract_decisions <- function(fit, true_effects, mu_control = 0,
                              prob_threshold = 0.975) {
  draws <- as_draws_df(fit)
  map_dfr(seq_along(true_effects), function(j) {
    eff <- draws$b_treatment +
      draws[[paste0("r_basket[", j, ",treatment]")]]
    tibble(
      basket      = j,
      true_effect = true_effects[j],
      post_mean   = mean(eff),
      prob_eff    = mean(eff > 0),
      decision    = mean(eff > 0) > prob_threshold
    )
  })
}

extract_decisions_fdr <- function(fit, true_effects, mu_control = 0.8,
                                  prob_threshold = 0.975) {
  draws <- as_draws_df(fit)
  map_dfr(seq_along(true_effects), function(j) {
    eff <- draws$b_treatment +
      draws[[paste0("r_basket[", j, ",treatment]")]]
    tibble(
      basket      = j,
      true_effect = true_effects[j] - mu_control,
      post_mean   = mean(eff),
      prob_eff    = mean(eff > 0),
      decision    = mean(eff > 0) > prob_threshold
    )
  })
}



## run simulation for operating characteristics basket

run_simulation_generic <- function(n_sim, n_baskets, n, mu_treat, sd, scenario,
                                   mu_control = 0,
                                   sim_fn     = simulate_basket,
                                   extract_fn = extract_decisions) {
  if (length(n) == 1L) n <- rep(n, n_baskets)
  stopifnot(length(n) == n_baskets)
  
  results <- vector("list", n_sim)
  
  for (s in seq_len(n_sim)) {
    data <- map_dfr(seq_len(n_baskets),
                    ~ sim_fn(.x, n[.x], mu_treat[.x], mu_control = mu_control, sd = sd)
    )
    fit <- fit_basket_model(
      data,
      heterogeneity = switch(scenario, P1 = "high", P2 = "low", P3 = "moderate"),
      seed = 1000 + s
    )
    results[[s]] <- extract_fn(fit, mu_treat, mu_control = mu_control)
  }
  
  bind_rows(results, .id = "sim")
}

run_simulation <- function(n_sim, n_baskets, n, mu_treat, sd, scenario) {
  run_simulation_generic(
    n_sim      = n_sim,
    n_baskets  = n_baskets,
    n          = n,
    mu_treat   = mu_treat,
    sd         = sd,
    scenario   = scenario,
    mu_control = 0,
    sim_fn     = simulate_basket,
    extract_fn = extract_decisions
  )
}

run_simulation_fdr <- function(n_sim, n_baskets, n, mu_treat, sd, scenario,
                               mu_control = 0.8) {
  run_simulation_generic(
    n_sim      = n_sim,
    n_baskets  = n_baskets,
    n          = n,
    mu_treat   = mu_treat,
    sd         = sd,
    scenario   = scenario,
    mu_control = mu_control,
    sim_fn     = simulate_basket,      # stessa funzione, mu_control diverso
    extract_fn = extract_decisions_fdr
  )
}

## performance calculations
evaluate_performance <- function(res) {
  res %>%
    group_by(basket) %>%
    summarise(
      power = mean(decision[true_effect != 0]),
      .groups = "drop"
    )
}

evaluate_performance_fdr <- function(res) {
  res %>%
    group_by(sim) %>%
    summarise(
      any_fp = any(decision & true_effect == 0),
      fp     = sum(decision & true_effect == 0),
      rej    = sum(decision),
      .groups = "drop"
    ) %>%
    summarise(
      FWER = mean(any_fp),
      FDR  = mean(ifelse(rej > 0, fp / rej, 0)),
      .groups = "drop"
    )
}


## power curve 
find_power_curve <- function(
    n_grid        = seq.int(from = 6, to = 30, by = 2),
    n_sim         = 500,
    n_baskets     = 5,
    mu_treat,
    sd            = 1,
    scenario_code) {
  
  map_dfr(n_grid, function(n) {
    sim <- run_simulation(
      n_sim     = n_sim,
      n_baskets = n_baskets,
      n         = n,
      mu_treat  = mu_treat,
      sd        = sd,
      scenario  = scenario_code
    )
    power <- evaluate_performance(sim) %>% pull(power) %>% mean()
    tibble(
      scenario    = scenario_code,
      n           = n,
      power       = power,
      effect_size = mean(mu_treat)
    )
  })
}



find_error_curve_fdr <- function(
    n_grid        = seq.int(from = 6, to = 30, by = 2),
    n_sim         = 500,
    n_baskets     = 5,
    mu_treat,
    mu_control    = 0.8,
    sd            = 1,
    scenario_code) {
  
  map_dfr(n_grid, function(n) {
    sim <- run_simulation_fdr(
      n_sim      = n_sim,
      n_baskets  = n_baskets,
      n          = rep(n, n_baskets),
      mu_treat   = mu_treat,
      mu_control = mu_control,
      sd         = sd,
      scenario   = scenario_code
    )
    perf <- evaluate_performance_fdr(sim)
    tibble(
      scenario    = scenario_code,
      n           = n,
      FDR         = perf$FDR,
      FWER        = perf$FWER,
      effect_size = mu_treat[1] - mu_control
    )
  })
}