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
  "ggplot2",
  "lme4",
  "lmerTest"
)

invisible(lapply(pkgs, load_or_install))


# ==========================================================
# 1) BAYESIANO
# ==========================================================

## Simulate basket
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


## Fit model for heterogeneity
fit_basket_model <- function(data, heterogeneity = c("high", "low", "moderate"),
                             iter = 5000, chains = 4, seed = 123) {
  heterogeneity <- match.arg(heterogeneity)
  
  cauchy_scale <- switch(heterogeneity,
                         high     = 0.8,
                         low      = 0.2,
                         moderate = 0.4)
  
  priors <- c(
    prior(student_t(3, 0, 5), class = Intercept),
    prior(normal(0, 1),        class = b),
    set_prior(paste0("cauchy(0, ", cauchy_scale, ")"),
              class = "sd", group = "basket", coef = "treatment"),
    prior(cauchy(0, 0.5), class = sigma)
  )
  
  brm(
    y ~ treatment + (1 | basket) + (0 + treatment | basket),
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


## Extract decisions for MCMC with CI
extract_decisions <- function(fit, true_effects, mu_control = 0,
                              prob_threshold = 0.975, ci_level = 0.95) {
  draws <- as_draws_df(fit)
  alpha <- (1 - ci_level) / 2
  
  map_dfr(seq_along(true_effects), function(j) {
    eff <- draws$b_treatment +
      draws[[paste0("r_basket[", j, ",treatment]")]]
    
    ci_lower <- quantile(eff, alpha)
    ci_upper <- quantile(eff, 1 - alpha)
    
    tibble(
      basket      = j,
      true_effect = true_effects[j],
      post_mean   = mean(eff),
      prob_eff    = mean(eff > 0),
      decision    = mean(eff > 0) > prob_threshold,
      ci_lower    = ci_lower,
      ci_upper    = ci_upper
    )
  })
}

extract_decisions_fdr <- function(fit, true_effects, mu_control = 0.8,
                                  prob_threshold = 0.975, ci_level = 0.95) {
  draws <- as_draws_df(fit)
  alpha <- (1 - ci_level) / 2
  
  map_dfr(seq_along(true_effects), function(j) {
    eff <- draws$b_treatment +
      draws[[paste0("r_basket[", j, ",treatment]")]]
    
    true_eff <- true_effects[j] - mu_control
    
    ci_lower <- quantile(eff, alpha)
    ci_upper <- quantile(eff, 1 - alpha)
    
    tibble(
      basket      = j,
      true_effect = true_eff,
      post_mean   = mean(eff),
      prob_eff    = mean(eff > 0),
      decision    = mean(eff > 0) > prob_threshold,
      ci_lower    = ci_lower,
      ci_upper    = ci_upper
    )
  })
}


# ----------------------------------------------------------
# FIX: aggiunta extract_tau_draws() — estrae draws di τ e σ
# dal fit brms e calcola I² per ogni campione MCMC.
# Viene chiamata internamente da run_simulation_generic().
# ----------------------------------------------------------
extract_tau_draws <- function(fit, scenario, sim_id) {
  draws <- as_draws_df(fit)
  tau_draws   <- draws$sd_basket__treatment
  sigma_draws <- draws$sigma
  I2_draws    <- tau_draws^2 / (tau_draws^2 + sigma_draws^2)
  
  tibble(
    sim      = sim_id,
    scenario = scenario,
    tau      = tau_draws,
    sigma    = sigma_draws,
    I2       = I2_draws
  )
}


## Run simulation — FIX: ora restituisce lista con $decisions e $tau_df
run_simulation_generic <- function(n_sim, n_baskets, n, mu_treat, sd, scenario,
                                   mu_control = 0,
                                   sim_fn     = simulate_basket,
                                   extract_fn = extract_decisions) {
  if (length(n) == 1L) n <- rep(n, n_baskets)
  stopifnot(length(n) == n_baskets)
  
  decisions_list <- vector("list", n_sim)
  tau_list       <- vector("list", n_sim)
  
  for (s in seq_len(n_sim)) {
    data <- map_dfr(seq_len(n_baskets),
                    ~ sim_fn(.x, n[.x], mu_treat[.x], mu_control = mu_control, sd = sd))
    
    fit <- fit_basket_model(
      data,
      heterogeneity = switch(scenario, P1 = "high", P2 = "low", P3 = "moderate"),
      seed = 1000 + s
    )
    
    decisions_list[[s]] <- extract_fn(fit, mu_treat, mu_control = mu_control)
    
    # FIX: salva draws di τ, σ e I² per ogni simulazione
    tau_list[[s]] <- extract_tau_draws(fit, scenario, sim_id = s)
  }
  
  list(
    decisions = bind_rows(decisions_list, .id = "sim"),
    tau_df    = bind_rows(tau_list)
  )
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
    sim_fn     = simulate_basket,
    extract_fn = extract_decisions_fdr
  )
}


## Performance calculations
# FIX: aggiunto MAPE con gestione true_effect == 0
evaluate_performance <- function(res) {
  res %>%
    group_by(basket) %>%
    summarise(
      power = mean(decision[true_effect != 0]),
      
      MAE  = mean(abs(post_mean - true_effect), na.rm = TRUE),
      MAPE = {
        vals <- ifelse(
          true_effect != 0,
          abs((post_mean - true_effect) / true_effect),
          NA_real_
        )
        if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
      },
      
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      ci_lower_mean = mean(ci_lower, na.rm = TRUE),
      ci_upper_mean = mean(ci_upper, na.rm = TRUE),
      coverage      = mean(true_effect >= ci_lower & true_effect <= ci_upper, na.rm = TRUE),
      
      .groups = "drop"
    )
}

# FIX: aggiunto MAPE in evaluate_performance_fdr
evaluate_performance_fdr <- function(res) {
  
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
    group_by(basket) %>%
    summarise(
      MAE  = mean(abs(post_mean - true_effect), na.rm = TRUE),
      MAPE = {
        vals <- ifelse(
          true_effect != 0,
          abs((post_mean - true_effect) / true_effect),
          NA_real_
        )
        if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
      },
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


# ----------------------------------------------------------
# FIX: summarise_tau() — calcola media, mediana e CI al 95%
# di τ e I² aggregando tutti i draws di tutte le simulazioni.
# Da chiamare sul $tau_df restituito da run_simulation_generic().
# ----------------------------------------------------------
summarise_tau <- function(tau_df) {
  tau_df |>
    group_by(scenario) |>
    summarise(
      tau_mean   = mean(tau),
      tau_median = median(tau),
      tau_lower  = quantile(tau, 0.025),
      tau_upper  = quantile(tau, 0.975),
      I2_mean    = mean(I2) * 100,
      I2_median  = median(I2) * 100,
      I2_lower   = quantile(I2, 0.025) * 100,
      I2_upper   = quantile(I2, 0.975) * 100,
      .groups    = "drop"
    )
}


# ----------------------------------------------------------
# FIX: print_tau_ci() — stampa una tabella leggibile con
# τ e I² (media ± 95% CI) per ogni scenario.
# Utile da chiamare direttamente nel paper/report.
# ----------------------------------------------------------
print_tau_ci <- function(tau_df) {
  summarise_tau(tau_df) |>
    mutate(
      `τ  [95% CI]`  = sprintf("%.3f  [%.3f – %.3f]", tau_mean,  tau_lower,  tau_upper),
      `I² [95% CI]`  = sprintf("%.1f%%  [%.1f%% – %.1f%%]", I2_mean, I2_lower, I2_upper)
    ) |>
    select(scenario, `τ  [95% CI]`, `I² [95% CI]`) |>
    print()
}


## Power curve
# FIX: ora restituisce anche tau_df aggregato per n
find_power_curve <- function(
    n_grid        = seq.int(from = 6, to = 30, by = 2),
    n_sim         = 500,
    n_baskets     = 5,
    mu_treat,
    sd            = 1,
    scenario_code) {
  
  results <- map(n_grid, function(n) {
    sim <- run_simulation(
      n_sim     = n_sim,
      n_baskets = n_baskets,
      n         = n,
      mu_treat  = mu_treat,
      sd        = sd,
      scenario  = scenario_code
    )
    
    power <- evaluate_performance(sim$decisions) %>% pull(power) %>% mean()
    
    list(
      power_row = tibble(
        scenario    = scenario_code,
        n           = n,
        power       = power,
        effect_size = mean(mu_treat)
      ),
      tau_df = sim$tau_df |> mutate(n = n)
    )
  })
  
  list(
    power_curve = map_dfr(results, ~ .x$power_row),
    tau_df      = map_dfr(results, ~ .x$tau_df)
  )
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
    # FIX: accesso corretto a $decisions dopo la modifica di run_simulation_generic
    perf <- evaluate_performance_fdr(sim$decisions)
    tibble(
      scenario    = scenario_code,
      n           = n,
      FDR         = perf$error_rates$FDR,
      FWER        = perf$error_rates$FWER,
      effect_size = mu_treat[1] - mu_control
    )
  })
}


# ==========================================================
# 2) FREQUENTISTA
# ==========================================================

simulate_basket_freq <- function(id, n_control, n_treat, mu_treat, mu_control = 0, sd = 1) {
  y <- c(
    rnorm(n_control, mu_control, sd),
    rnorm(n_treat,   mu_treat,   sd)
  )
  treatment <- rep(c(0, 1), times = c(n_control, n_treat))
  tibble(
    basket    = factor(id),
    treatment = treatment,
    y         = y
  )
}

fit_basket_freq <- function(data, heterogeneity = c("high", "moderate", "low")) {
  heterogeneity <- match.arg(heterogeneity)
  
  basket_stats <- data %>%
    group_by(basket, treatment) %>%
    summarise(
      n    = n(),
      mean = mean(y),
      var  = var(y),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = treatment, values_from = c(n, mean, var),
                names_prefix = "treat") %>%
    mutate(
      eff_raw = mean_treat1 - mean_treat0,
      se_raw  = sqrt(var_treat0 / n_treat0 + var_treat1 / n_treat1)
    )
  
  weight <- switch(heterogeneity,
                   high     = 0.1,
                   moderate = 0.5,
                   low      = 0.9)
  
  mu_global <- sum(basket_stats$eff_raw / basket_stats$se_raw^2) /
    sum(1 / basket_stats$se_raw^2)
  
  basket_stats %>%
    mutate(
      eff_borrow = (1 - weight) * eff_raw + weight * mu_global,
      se_borrow  = se_raw * sqrt(1 - weight + weight * 0.1)
    )
}

extract_decisions_basket_freq <- function(fit, true_effects = NULL, alpha = 0.05) {
  zcrit <- qnorm(1 - alpha / 2)
  if (is.null(true_effects)) true_effects <- rep(0, nrow(fit))
  
  fit %>%
    mutate(
      ci_lower    = eff_borrow - zcrit * se_borrow,
      ci_upper    = eff_borrow + zcrit * se_borrow,
      decision    = !(ci_lower <= 0 & ci_upper >= 0),
      true_effect = true_effects
    ) %>%
    select(basket, eff_borrow, ci_lower, ci_upper, decision, true_effect)
}

run_simulation_basket_freq <- function(n_sim, n_baskets, n_control, n_treat,
                                       mu_treat, sd = 1, heterogeneity = "high") {
  results <- vector("list", n_sim)
  for (s in seq_len(n_sim)) {
    data <- map_dfr(seq_len(n_baskets),
                    ~ simulate_basket_freq(.x,
                                           n_control = n_control[.x],
                                           n_treat   = n_treat[.x],
                                           mu_treat  = mu_treat[.x],
                                           sd        = sd))
    fit <- fit_basket_freq(data, heterogeneity)
    results[[s]] <- extract_decisions_basket_freq(fit, true_effects = mu_treat - 0)
  }
  bind_rows(results, .id = "sim")
}

# FIX: aggiunto MAPE con gestione true_effect == 0
evaluate_performance_basket_freq <- function(res) {
  res %>%
    group_by(basket) %>%
    summarise(
      power         = mean(decision),
      MAE           = mean(abs(eff_borrow - true_effect), na.rm = TRUE),
      MAPE          = {
        vals <- ifelse(
          true_effect != 0,
          abs((eff_borrow - true_effect) / true_effect),
          NA_real_
        )
        if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
      },
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      mean_ci_lower = mean(ci_lower, na.rm = TRUE),
      mean_ci_upper = mean(ci_upper, na.rm = TRUE),
      coverage      = mean(ci_lower <= true_effect & ci_upper >= true_effect, na.rm = TRUE),
      .groups       = "drop"
    )
}

run_all_basket_freq <- function(scenarios) {
  map(scenarios, function(sc) {
    sim <- run_simulation_basket_freq(
      n_sim         = n_sim_real,
      n_baskets     = sc$n_baskets,
      n_control     = sc$n,
      n_treat       = sc$n,
      mu_treat      = sc$mu_treat,
      sd            = sc$sd,
      heterogeneity = switch(sc$scenario_code,
                             P1 = "high",
                             P2 = "low",
                             P3 = "moderate")
    )
    list(
      scenario    = sc$name,
      performance = evaluate_performance_basket_freq(sim)
    )
  })
}

build_summary_basket_freq <- function(results) {
  map_dfr(results, ~ .x$performance %>%
            select(basket, power, MAE, MAPE, mean_ci_width,
                   mean_ci_lower, mean_ci_upper, coverage) %>%
            mutate(scenario = .x$scenario))
}

find_power_curve_basket_freq <- function(n_grid = seq(6, 30, 2), n_sim = 500,
                                         n_baskets = 5, mu_treat, sd = 1,
                                         heterogeneity = "high", scenario_label = "") {
  map_dfr(n_grid, function(n) {
    sim <- run_simulation_basket_freq(
      n_sim         = n_sim,
      n_baskets     = n_baskets,
      n_control     = rep(n, n_baskets),
      n_treat       = rep(n, n_baskets),
      mu_treat      = mu_treat,
      sd            = sd,
      heterogeneity = heterogeneity
    )
    tibble(
      scenario    = scenario_label,
      n           = n,
      power       = mean(evaluate_performance_basket_freq(sim)$power),
      effect_size = mean(mu_treat)
    )
  })
}

plot_power_curve_basket_freq <- function(df, title = "") {
  ggplot(df, aes(x = n, y = power, color = scenario, group = scenario)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    labs(title = title, x = "Sample size per basket", y = "Power") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
}

extract_decisions_fdr_basket_freq <- function(fit, mu_treat, mu_control = 0.8, alpha = 0.05) {
  zcrit       <- qnorm(1 - alpha / 2)
  true_effect <- mu_treat - mu_control
  fit %>%
    mutate(
      ci_lower    = eff_borrow - zcrit * se_borrow,
      ci_upper    = eff_borrow + zcrit * se_borrow,
      decision    = !(ci_lower <= 0 & ci_upper >= 0),
      true_effect = true_effect
    ) %>%
    select(basket, eff_borrow, ci_lower, ci_upper, decision, true_effect)
}

run_simulation_fdr_basket_freq <- function(n_sim, n_baskets, n_control, n_treat,
                                           mu_treat, mu_control = 0.8, sd = 1,
                                           heterogeneity = "high") {
  results <- vector("list", n_sim)
  for (s in seq_len(n_sim)) {
    data <- map_dfr(seq_len(n_baskets),
                    ~ simulate_basket_freq(.x,
                                           n_control  = n_control[.x],
                                           n_treat    = n_treat[.x],
                                           mu_treat   = mu_treat[.x],
                                           mu_control = mu_control,
                                           sd         = sd))
    fit <- fit_basket_freq(data, heterogeneity)
    results[[s]] <- extract_decisions_fdr_basket_freq(fit, mu_treat, mu_control)
  }
  bind_rows(results, .id = "sim")
}

# FIX: aggiunto MAPE in evaluate_performance_fdr_basket_freq
evaluate_performance_fdr_basket_freq <- function(res) {
  
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
    group_by(basket) %>%
    summarise(
      MAE  = mean(abs(eff_borrow - true_effect), na.rm = TRUE),
      MAPE = {
        vals <- ifelse(
          true_effect != 0,
          abs((eff_borrow - true_effect) / true_effect),
          NA_real_
        )
        if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
      },
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

run_all_fdr_basket_freq <- function(scenarios) {
  map(scenarios, function(sc) {
    sim <- run_simulation_fdr_basket_freq(
      n_sim         = n_sim_real,
      n_baskets     = sc$n_baskets,
      n_control     = sc$n,
      n_treat       = sc$n,
      mu_treat      = sc$mu_treat,
      mu_control    = 0.8,
      sd            = sc$sd,
      heterogeneity = switch(sc$scenario_code,
                             P1 = "high",
                             P2 = "low",
                             P3 = "moderate")
    )
    list(scenario    = sc$name,
         performance = evaluate_performance_fdr_basket_freq(sim))
  })
}

# FIX: build_summary_fdr_basket_freq ora restituisce lista separata
# con $error_rates (FWER/FDR) e $estimation (MAE, MAPE, CI, coverage),
# invece di cbind() che produceva colonne duplicate e silenziose
build_summary_fdr_basket_freq <- function(results) {
  
  error_df <- map_dfr(results, ~ .x$performance$error_rates |>
                        mutate(scenario = .x$scenario))
  
  estimation_df <- map_dfr(results, ~ .x$performance$estimation |>
                             select(basket, MAE, MAPE, mean_ci_width,
                                    ci_lower_mean, ci_upper_mean, coverage) |>
                             mutate(scenario = .x$scenario))
  
  list(
    error_rates = error_df,
    estimation  = estimation_df
  )
}

find_error_curve_fdr_basket_freq <- function(n_grid = seq(6, 30, 2), n_sim = 500,
                                             n_baskets = 5, mu_treat, mu_control = 0.8,
                                             sd = 1, heterogeneity = "high") {
  map_dfr(n_grid, function(n) {
    sim <- run_simulation_fdr_basket_freq(
      n_sim         = n_sim,
      n_baskets     = n_baskets,
      n_control     = rep(n, n_baskets),
      n_treat       = rep(n, n_baskets),
      mu_treat      = mu_treat,
      mu_control    = mu_control,
      sd            = sd,
      heterogeneity = heterogeneity
    )
    perf <- evaluate_performance_fdr_basket_freq(sim)
    tibble(
      n           = n,
      FDR         = perf$error_rates$FDR,
      FWER        = perf$error_rates$FWER,
      effect_size = mu_treat[1] - mu_control
    )
  })
}