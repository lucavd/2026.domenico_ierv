load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

invisible(lapply(
  c("tidyverse", "brms", "posterior", "patchwork", "ggplot2"),
  load_or_install
))


## Data generation
simulate_umbrella <- function(id, n_control, n_treat,
                              mu_treat, mu_control = NULL, sd = 1,
                              standardise = TRUE) {
  if (is.null(mu_control)) mu_control <- if (standardise) 0 else 0.8
  
  y_raw <- c(rnorm(n_control, mu_control, sd),
             rnorm(n_treat,   mu_treat,   sd))
  
  tibble(
    arm       = factor(id),
    treatment = rep(c(0L, 1L), times = c(n_control, n_treat)),
    y         = if (standardise) (y_raw - mu_control) / sd else y_raw
  )
}


## Model fit for Heterogeneity
.fit_umbrella <- function(data, priors, iter = 5000, chains = 4, seed = 123) {
  brm(
    y ~ treatment + (1 | arm) + (0 + treatment | arm),
    data    = data,
    prior   = priors,
    iter    = iter,
    warmup  = iter / 2,
    chains  = chains,
    seed    = seed,
    cores   = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 18),
    refresh = 0
  )
}

fit_umbrella_model <- function(data, scenario, iter = 5000, chains = 4, seed = 123) {
  priors <- switch(scenario,
                   P1 = c(prior(student_t(3, 0, 10), class = Intercept),
                          prior(normal(0, 1),        class = b),
                          prior(cauchy(0, 0.32),     class = sd, group = arm, coef = treatment),
                          prior(cauchy(0, 0.1),      class = sigma)),
                   P2 = c(prior(student_t(3, 0, 2),  class = Intercept),
                          prior(normal(0, 1),        class = b),
                          prior(cauchy(0, 0.4),      class = sd, group = arm, coef = treatment),
                          prior(cauchy(0, 1.2),      class = sigma)),
                   P3 = c(prior(student_t(3, 0, 5),  class = Intercept),
                          prior(normal(0, 1),        class = b),
                          prior(cauchy(0, 0.4),      class = sd, group = arm, coef = treatment),
                          prior(cauchy(0, 0.6),      class = sigma)),
                   stop("Unknown scenario: ", scenario)
  )
  .fit_umbrella(data, priors, iter, chains, seed)
}


## Decisions extraction for MCMC
extract_decisions_umbrella <- function(fit, true_effects, prob_threshold = 0.975) {
  draws     <- as_draws_df(fit)
  arm_names <- levels(fit$data$arm)
  
  map_dfr(seq_along(true_effects), function(j) {
    eff <- draws[[paste0("b_arm", arm_names[j], ":treatment")]]
    tibble(
      arm         = j,
      true_effect = true_effects[j],
      post_mean   = mean(eff),
      prob_eff    = mean(eff > 0),
      decision    = mean(eff > 0) > prob_threshold
    )
  })
}

extract_decisions_fdr_umb <- function(fit, true_effects,
                                      mu_control = 0.8, prob_threshold = 0.975) {
  draws <- as_draws_df(fit)
  
  map_dfr(seq_along(true_effects), function(j) {
    eff <- draws$b_treatment + draws[[paste0("r_arm[", j, ",treatment]")]]
    tibble(
      arm         = j,
      true_effect = true_effects[j] - mu_control,
      post_mean   = mean(eff),
      prob_eff    = mean(eff > 0),
      decision    = mean(eff > 0) > prob_threshold
    )
  })
}


## Simulations
run_simulation_umbrella <- function(n_sim, n_arms, n_control, n_treat,
                                    mu_treat, sd, scenario) {
  map_dfr(seq_len(n_sim), function(s) {
    data <- map_dfr(seq_len(n_arms),
                    ~ simulate_umbrella(.x, n_control[.x], n_treat[.x],
                                        mu_treat[.x], sd = sd, standardise = TRUE))
    fit  <- fit_umbrella_model(data, scenario, seed = 1000 + s)
    extract_decisions_umbrella(fit, mu_treat)
  }, .id = "sim")
}

run_simulation_fdr_umb <- function(n_sim, n_arms, n_control, n_treat,
                                   mu_treat, mu_control = 0.8, sd, scenario) {
  map_dfr(seq_len(n_sim), function(s) {
    data <- map_dfr(seq_len(n_arms),
                    ~ simulate_umbrella(.x, n_control[.x], n_treat[.x],
                                        mu_treat[.x], mu_control, sd, standardise = FALSE))
    fit  <- fit_umbrella_model(data, scenario, seed = 1000 + s)
    extract_decisions_fdr_umb(fit, mu_treat, mu_control)
  }, .id = "sim")
}


## Metrics' Calculations
evaluate_performance_umbrella <- function(res) {
  res |>
    group_by(arm) |>
    summarise(power = mean(decision, na.rm = TRUE), .groups = "drop")
}

evaluate_performance_fdr_umb <- function(res) {
  res |>
    group_by(sim) |>
    summarise(
      any_fp = any(decision & true_effect == 0),
      fp     = sum(decision & true_effect == 0),
      rej    = sum(decision),
      .groups = "drop"
    ) |>
    summarise(
      FWER = mean(any_fp),
      FDR  = mean(ifelse(rej > 0, fp / rej, 0))
    )
}


## Summary
run_all <- function(scenarios) {
  map(scenarios, function(sc) {
    sim <- run_simulation_umbrella(
      n_sim     = n_sim_real,
      n_arms    = sc$n_arms,
      n_control = sc$n_control,
      n_treat   = sc$n_treat,
      mu_treat  = sc$mu_treat,
      sd        = sc$sd,
      scenario  = sc$scenario_code
    )
    list(scenario    = sc$name,
         performance = evaluate_performance_umbrella(sim))
  })
}

build_summary <- function(results) {
  map_dfr(results, ~ .x$performance |> mutate(scenario = .x$scenario))
}


## Power and FDR curves
find_power_curve_umbrella <- function(
    n_grid         = seq.int(6, 30, by = 2),
    n_sim          = 500,
    n_arms         = 5,
    mu_treat,
    sd             = 1,
    scenario_code,
    scenario_label = scenario_code) {
  
  map_dfr(n_grid, function(n) {
    sim   <- run_simulation_umbrella(
      n_sim     = n_sim,          n_arms    = n_arms,
      n_control = rep(n, n_arms), n_treat   = rep(n, n_arms),
      mu_treat  = mu_treat,       sd        = sd,
      scenario  = scenario_code
    )
    power <- evaluate_performance_umbrella(sim) |> pull(power) |> mean(na.rm = TRUE)
    tibble(scenario = scenario_label, n = n, power = power, effect_size = mu_treat[1])
  })
}

find_error_curve_fdr_umb <- function(
    n_grid        = seq.int(6, 30, by = 2),
    n_sim         = 500,
    n_arms        = 5,
    mu_treat,
    mu_control    = 0.8,
    sd            = 1,
    scenario_code) {
  
  map_dfr(n_grid, function(n) {
    sim  <- run_simulation_fdr_umb(
      n_sim     = n_sim,          n_arms    = n_arms,
      n_control = rep(n, n_arms), n_treat   = rep(n, n_arms),
      mu_treat  = mu_treat,       mu_control = mu_control,
      sd        = sd,             scenario   = scenario_code
    )
    perf <- evaluate_performance_fdr_umb(sim)
    tibble(scenario    = scenario_code, n = n,
           FDR         = perf$FDR,      FWER = perf$FWER,
           effect_size = mu_treat[1] - mu_control)
  })
}


## Plot Functions
plot_power_by_arm <- function(df, title_suffix) {
  df |>
    filter(!is.nan(power)) |>
    ggplot(aes(x = arm, y = power, color = scenario, group = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = 1:5) +
    labs(title = paste("Power by Arm –", title_suffix),
         x = "Arm", y = "Power") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_power_curve <- function(df, title_suffix) {
  ggplot(df, aes(x = n, y = power, color = scenario)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    facet_grid(~ effect_size,
               labeller = labeller(effect_size = \(x) paste0("Effect size = ", x))) +
    geom_hline(yintercept = target_power, linetype = "dashed", color = "black") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = paste("Power Curve –", title_suffix),
         x = "Sample Size per Arm", y = "Power", color = "Prior scenario") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}