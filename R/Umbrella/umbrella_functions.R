load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

invisible(lapply(
  c("tidyverse", "brms", "posterior", "patchwork", "ggplot2", "lme4", "lmerTest"),
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


## Decisions extraction for MCMC and Credible Intervals
extract_decisions_umbrella <- function(fit, true_effects,
                                       prob_threshold = 0.975,
                                       ci_level = 0.95) {
  
  draws <- as_draws_df(fit)
  alpha <- (1 - ci_level) / 2
  
  map_dfr(seq_along(true_effects), function(j) {
    
    eff <- draws$b_treatment +
      draws[[paste0("r_arm[", j, ",treatment]")]]
    
    ci_lower <- quantile(eff, alpha)
    ci_upper <- quantile(eff, 1 - alpha)
    
    tibble(
      arm         = j,
      true_effect = true_effects[j],
      post_mean   = mean(eff),
      prob_eff    = mean(eff > 0),
      decision    = mean(eff > 0) > prob_threshold,
      ci_lower    = ci_lower,
      ci_upper    = ci_upper
    )
  })
}

extract_decisions_fdr_umb <- function(fit, true_effects,
                                      mu_control = 0.8,
                                      prob_threshold = 0.975,
                                      ci_level = 0.95) {
  
  draws <- as_draws_df(fit)
  alpha <- (1 - ci_level) / 2
  
  map_dfr(seq_along(true_effects), function(j) {
    
    eff <- draws$b_treatment +
      draws[[paste0("r_arm[", j, ",treatment]")]]
    
    true_eff <- true_effects[j] - mu_control
    
    ci_lower <- quantile(eff, alpha)
    ci_upper <- quantile(eff, 1 - alpha)
    
    tibble(
      arm         = j,
      true_effect = true_eff,
      post_mean   = mean(eff),
      prob_eff    = mean(eff > 0),
      decision    = mean(eff > 0) > prob_threshold,
      ci_lower    = ci_lower,
      ci_upper    = ci_upper
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
    summarise(
      power = mean(decision, na.rm = TRUE),
      
      MAE  = mean(abs(post_mean - true_effect), na.rm = TRUE),
      MAPE = mean(
        ifelse(true_effect != 0,
               abs((post_mean - true_effect) / true_effect),
               NA_real_),
        na.rm = TRUE
      ),
      
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      ci_lower_mean = mean(ci_lower, na.rm = TRUE),
      ci_upper_mean = mean(ci_upper, na.rm = TRUE),
      coverage      = mean(
        true_effect >= ci_lower & true_effect <= ci_upper,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
}

evaluate_performance_fdr_umb <- function(res) {
  
  # Error rates
  error_part <- res |>
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
  
  # Estimation metrics
  estimation_part <- res |>
    group_by(arm) |>
    summarise(
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      ci_lower_mean = mean(ci_lower, na.rm = TRUE),
      ci_upper_mean = mean(ci_upper, na.rm = TRUE),
      coverage      = mean(
        true_effect >= ci_lower & true_effect <= ci_upper,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  list(
    error_rates = error_part,
    estimation  = estimation_part
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
    
    sim <- run_simulation_umbrella(
      n_sim     = n_sim,
      n_arms    = n_arms,
      n_control = rep(n, n_arms),
      n_treat   = rep(n, n_arms),
      mu_treat  = mu_treat,
      sd        = sd,
      scenario  = scenario_code
    )
    
    perf <- evaluate_performance_umbrella(sim)
    
    tibble(
      scenario    = scenario_label,
      n           = n,
      power       = mean(perf$power, na.rm = TRUE),
      effect_size = mu_treat[1]
    )
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
    
    sim <- run_simulation_fdr_umb(
      n_sim      = n_sim,
      n_arms     = n_arms,
      n_control  = rep(n, n_arms),
      n_treat    = rep(n, n_arms),
      mu_treat   = mu_treat,
      mu_control = mu_control,
      sd         = sd,
      scenario   = scenario_code
    )
    
    perf <- evaluate_performance_fdr_umb(sim)
    
    tibble(
      scenario    = scenario_code,
      n           = n,
      FDR         = perf$error_rates$FDR,
      FWER        = perf$error_rates$FWER,
      effect_size = mu_treat[1] - mu_control
    )
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




# ==========================================================
# 1) Simulazione arms (frequentista)
# ==========================================================
simulate_umbrella_freq <- function(id, n_control, n_treat, mu_treat, mu_control = 0, sd = 1) {
  y <- c(
    rnorm(n_control, mu_control, sd),
    rnorm(n_treat, mu_treat, sd)
  )
  treatment <- rep(c(0,1), times = c(n_control, n_treat))
  tibble(
    arm       = factor(id),
    treatment = treatment,
    y         = y
  )
}

# ==========================================================
# 2) Fit frequentista con borrowing tra arms
# ==========================================================
fit_umbrella_freq <- function(data, heterogeneity = c("high","moderate","low")) {
  heterogeneity <- match.arg(heterogeneity)
  
  arm_stats <- data %>%
    group_by(arm, treatment) %>%
    summarise(
      n    = n(),
      mean = mean(y),
      var  = var(y),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = treatment, values_from = c(n, mean, var), names_prefix = "treat") %>%
    mutate(
      eff_raw = mean_treat1 - mean_treat0,
      se_raw  = sqrt(var_treat0 / n_treat0 + var_treat1 / n_treat1)
    )
  
  weight <- switch(heterogeneity,
                   high     = 0.1,
                   moderate = 0.5,
                   low      = 0.9)
  
  mu_global <- sum(arm_stats$eff_raw / arm_stats$se_raw^2) / sum(1 / arm_stats$se_raw^2)
  
  arm_stats %>%
    mutate(
      eff_borrow = (1 - weight) * eff_raw + weight * mu_global,
      se_borrow  = se_raw * sqrt(1 - weight + weight * 0.1)
    )
}

# ==========================================================
# 3) Estrarre decisioni (CI, alpha)
# ==========================================================
extract_decisions_umbrella_freq <- function(fit, true_effects = NULL, alpha = 0.05) {
  zcrit <- qnorm(1 - alpha / 2)
  if(is.null(true_effects)) true_effects <- rep(0, nrow(fit))
  
  fit %>%
    mutate(
      ci_lower = eff_borrow - zcrit * se_borrow,
      ci_upper = eff_borrow + zcrit * se_borrow,
      decision = !(ci_lower <= 0 & ci_upper >= 0),
      true_effect = true_effects
    ) %>%
    select(arm, eff_borrow, ci_lower, ci_upper, decision, true_effect)
}

# ==========================================================
# 4) Simulazione completa per scenario
# ==========================================================
run_simulation_umbrella_freq <- function(n_sim, n_arms, n_control, n_treat, mu_treat, sd = 1,
                                         heterogeneity = "high") {
  results <- vector("list", n_sim)
  for(s in seq_len(n_sim)) {
    data <- map_dfr(seq_len(n_arms),
                    ~ simulate_umbrella_freq(.x,
                                             n_control = n_control[.x],
                                             n_treat   = n_treat[.x],
                                             mu_treat  = mu_treat[.x],
                                             sd = sd))
    
    fit <- fit_umbrella_freq(data, heterogeneity)
    results[[s]] <- extract_decisions_umbrella_freq(fit, true_effects = mu_treat - 0.8)
  }
  bind_rows(results, .id = "sim")
}

# ==========================================================
# 5) Performance evaluation
# ==========================================================
evaluate_performance_umbrella_freq <- function(res) {
  res %>%
    group_by(arm) %>%
    summarise(
      power = mean(decision),
      MAE   = mean(abs(eff_borrow - true_effect), na.rm = TRUE),
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      mean_ci_lower = mean(ci_lower, na.rm = TRUE),
      mean_ci_upper = mean(ci_upper, na.rm = TRUE),
      coverage      = mean(ci_lower <= true_effect & ci_upper >= true_effect, na.rm = TRUE),
      .groups = "drop"
    )
}

# ==========================================================
# 6) Wrapper run_all / build_summary
# ==========================================================
run_all_freq <- function(scenarios) {
  map(scenarios, function(sc) {
    sim <- run_simulation_umbrella_freq(
      n_sim     = n_sim_real,
      n_arms    = sc$n_arms,
      n_control = sc$n_control,
      n_treat   = sc$n_treat,
      mu_treat  = sc$mu_treat,
      sd        = sc$sd,
      heterogeneity = switch(sc$scenario_code,
                             P1 = "high",
                             P2 = "low",
                             P3 = "moderate")
    )
    list(
      scenario = sc$name,
      performance = evaluate_performance_umbrella_freq(sim)
    )
  })
}

build_summary_freq <- function(results) {
  map_dfr(results, ~ .x$performance %>% mutate(scenario = .x$scenario))
}

# ==========================================================
# 7) Power curve (frequentista)
# ==========================================================
find_power_curve_umbrella_freq <- function(n_grid = seq(6,30,2), n_sim = 100, n_arms = 5,
                                           mu_treat, sd = 1, heterogeneity = "high",
                                           scenario_label = "") {
  map_dfr(n_grid, function(n) {
    sim <- run_simulation_umbrella_freq(
      n_sim     = n_sim,
      n_arms    = n_arms,
      n_control = rep(n, n_arms),
      n_treat   = rep(n, n_arms),
      mu_treat  = mu_treat,
      sd        = sd,
      heterogeneity = heterogeneity
    )
    
    power <- mean(evaluate_performance_umbrella_freq(sim)$power)
    
    tibble(
      scenario = scenario_label,
      n        = n,
      power    = power,
      effect_size = mean(mu_treat - 0.8)
    )
  })
}

# ==========================================================
# 8) Plot power curve
# ==========================================================
plot_power_curve_freq <- function(df, title = "") {
  ggplot(df, aes(x = n, y = power, color = scenario, group = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(title = title, x = "Sample size per arm", y = "Power") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

# ==========================================================
# 9) FDR extraction / simulation / evaluation
# ==========================================================
extract_decisions_fdr_umb_freq <- function(fit, mu_treat, mu_control = 0.8, alpha = 0.05) {
  zcrit <- qnorm(1 - alpha / 2)
  true_effect <- mu_treat - mu_control
  fit %>%
    mutate(
      ci_lower = eff_borrow - zcrit * se_borrow,
      ci_upper = eff_borrow + zcrit * se_borrow,
      decision = !(ci_lower <= 0 & ci_upper >= 0),
      true_effect = true_effect
    ) %>%
    select(arm, eff_borrow, ci_lower, ci_upper, decision, true_effect)
}

run_simulation_fdr_umbrella_freq <- function(n_sim, n_arms, n_control, n_treat,
                                             mu_treat, mu_control = 0.8, sd = 1,
                                             heterogeneity = "high") {
  results <- vector("list", n_sim)
  for(s in seq_len(n_sim)) {
    data <- map_dfr(seq_len(n_arms),
                    ~ simulate_umbrella_freq(.x,
                                             n_control = n_control[.x],
                                             n_treat   = n_treat[.x],
                                             mu_treat  = mu_treat[.x],
                                             mu_control = mu_control,
                                             sd = sd))
    fit <- fit_umbrella_freq(data, heterogeneity)
    results[[s]] <- extract_decisions_fdr_umb_freq(fit, mu_treat, mu_control)
  }
  bind_rows(results, .id = "sim")
}

evaluate_performance_fdr_umbrella_freq <- function(res) {
  error_part <- res %>%
    group_by(sim) %>%
    summarise(
      any_fp = any(decision & true_effect == 0),
      fp     = sum(decision & true_effect == 0),
      rej    = sum(decision),
      .groups = "drop"
    ) %>%
    summarise(
      FWER = mean(any_fp),
      FDR  = mean(ifelse(rej > 0, fp / rej, 0))
    )
  
  estimation_part <- res %>%
    group_by(arm) %>%
    summarise(
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      ci_lower_mean = mean(ci_lower, na.rm = TRUE),
      ci_upper_mean = mean(ci_upper, na.rm = TRUE),
      coverage      = mean(true_effect >= ci_lower & true_effect <= ci_upper, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    error_rates = error_part,
    estimation  = estimation_part
  )
}

find_error_curve_fdr_umbrella_freq <- function(n_grid = seq(6,30,2), n_sim = 100, n_arms = 5,
                                               mu_treat, mu_control = 0.8, sd = 1,
                                               heterogeneity = "high") {
  map_dfr(n_grid, function(n) {
    sim <- run_simulation_fdr_umbrella_freq(
      n_sim      = n_sim,
      n_arms     = n_arms,
      n_control  = rep(n, n_arms),
      n_treat    = rep(n, n_arms),
      mu_treat   = mu_treat,
      mu_control = mu_control,
      sd         = sd,
      heterogeneity = heterogeneity
    )
    perf <- evaluate_performance_fdr_umbrella_freq(sim)
    tibble(
      n      = n,
      FDR    = perf$error_rates$FDR,
      FWER   = perf$error_rates$FWER,
      effect_size = mu_treat[1] - mu_control
    )
  })
}

run_all_fdr_freq <- function(scenarios) {
  map(scenarios, function(sc) {
    sim <- run_simulation_fdr_umbrella_freq(
      n_sim     = n_sim_real,
      n_arms    = sc$n_arms,
      n_control = sc$n_control,
      n_treat   = sc$n_treat,
      mu_treat  = sc$mu_treat,
      mu_control = 0.8,
      sd        = sc$sd,
      heterogeneity = switch(sc$scenario_code,
                             P1 = "high",
                             P2 = "low",
                             P3 = "moderate")
    )
    list(scenario    = sc$name,
         performance = evaluate_performance_fdr_umbrella_freq(sim))
  })
}

build_summary_fdr_freq <- function(results) {
  map_dfr(results, function(x) {
    perf <- x$performance$estimation
    perf$scenario <- x$scenario
    perf
  })
}