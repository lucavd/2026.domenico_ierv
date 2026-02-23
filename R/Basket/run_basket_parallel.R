##############################################################################
# Basket Trial Parallel Simulation Runner (Power + FDR)
# Runs all scenarios sequentially, each parallelized across 115 cores
##############################################################################

# ===========================
# Libraries
# ===========================
library(tidyverse)
library(brms)
library(posterior)
library(future)
library(future.apply)
library(telegram.bot)

# ===========================
# Setup
# ===========================
cat("[BASKET] Avvio script...\n"); flush.console()
script_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])))
csv_dir    <- file.path(script_dir, "../..")

readRenviron(file.path(script_dir, "../../.env"))
bot     <- Bot(token = Sys.getenv("TELEGRAM_BOT_TOKEN"))
chat_id <- Sys.getenv("TELEGRAM_CHAT_ID")
cat("[BASKET] Telegram OK\n"); flush.console()

# ===========================
# Source basket functions
# ===========================
source(file.path(script_dir, "basket_functions.R"))
cat("[BASKET] Funzioni caricate\n"); flush.console()

# ===========================
# Pre-compile brms models (once, before spawning workers)
# ===========================
cat("[BASKET] Pre-compilazione modelli Stan...\n"); flush.console()
dummy_data <- map_dfr(1:5, ~ simulate_basket(.x, 10, 0.5, mu_control = 0, sd = 1))

precompiled <- list()
for (het in c("high", "low", "moderate")) {
  cauchy_scale <- switch(het, high = 0.8, low = 0.2, moderate = 0.4)
  priors <- c(
    prior(student_t(3, 0, 5), class = Intercept),
    prior(normal(0, 1), class = b),
    set_prior(paste0("cauchy(0, ", cauchy_scale, ")"),
              class = "sd", group = "basket", coef = "treatment"),
    prior(cauchy(0, 0.5), class = sigma)
  )
  precompiled[[het]] <- brm(
    y ~ treatment + (1 | basket) + (0 + treatment | basket),
    data = dummy_data, prior = priors, chains = 0,
    control = list(adapt_delta = 0.999, max_treedepth = 18)
  )
  cat("  Compilato:", het, "\n"); flush.console()
}

# Override fit_basket_model: reuse pre-compiled model (no recompilation in workers)
fit_basket_model <- function(data, heterogeneity = c("high", "low", "moderate"),
                             iter = 5000, chains = 4, seed = 123) {
  heterogeneity <- match.arg(heterogeneity)
  update(precompiled[[heterogeneity]], newdata = data,
         iter = iter, warmup = iter / 2, chains = chains,
         seed = seed, cores = 1, refresh = 0)
}

# ===========================
# Parallel setup (after pre-compilation)
# ===========================
options(mc.cores = 128)
n_workers <- 115
plan(multisession, workers = n_workers)
cat("[BASKET] Plan multisession:", n_workers, "workers\n"); flush.console()

# ===========================
# Settings
# ===========================
set.seed(2025)
n_sim   <- 500
n_grid  <- seq.int(from = 6, to = 30, by = 2)

# Helper: restart workers to free /tmp space
restart_workers <- function() {
  plan(sequential)
  gc(verbose = FALSE)
  invisible(file.remove(list.files(tempdir(), full.names = TRUE, recursive = TRUE)))
  plan(multisession, workers = n_workers)
  cat("[BASKET] Workers riavviati, temp pulito\n"); flush.console()
}

# ===========================
# Parallel wrappers
# ===========================

run_simulation_par <- function(n_sim, n_baskets, n, mu_treat, sd, scenario) {
  if (length(n) == 1L) n <- rep(n, n_baskets)
  stopifnot(length(n) == n_baskets)
  results <- future_lapply(seq_len(n_sim), function(s) {
    data <- purrr::map_dfr(seq_len(n_baskets),
              ~ simulate_basket(.x, n[.x], mu_treat[.x], mu_control = 0, sd = sd))
    fit <- fit_basket_model(data,
      heterogeneity = switch(scenario, P1 = "high", P2 = "low", P3 = "moderate"),
      seed = 1000 + s)
    extract_decisions(fit, mu_treat, mu_control = 0)
  }, future.seed = 123L,
     future.packages = c("brms", "posterior", "purrr", "dplyr", "tibble"))
  bind_rows(results, .id = "sim")
}

run_simulation_fdr_par <- function(n_sim, n_baskets, n, mu_treat,
                                   mu_control = 0.8, sd, scenario) {
  if (length(n) == 1L) n <- rep(n, n_baskets)
  stopifnot(length(n) == n_baskets)
  results <- future_lapply(seq_len(n_sim), function(s) {
    data <- purrr::map_dfr(seq_len(n_baskets),
              ~ simulate_basket(.x, n[.x], mu_treat[.x], mu_control = mu_control, sd = sd))
    fit <- fit_basket_model(data,
      heterogeneity = switch(scenario, P1 = "high", P2 = "low", P3 = "moderate"),
      seed = 1000 + s)
    extract_decisions_fdr(fit, mu_treat, mu_control = mu_control)
  }, future.seed = 123L,
     future.packages = c("brms", "posterior", "purrr", "dplyr", "tibble"))
  bind_rows(results, .id = "sim")
}

##############################################################################
# ======================
# POWER SIMULATIONS
# ======================
##############################################################################

bot$sendMessage(chat_id = chat_id,
                text = paste0("[BASKET] Inizio POWER - ",
                              n_sim, " sim su ", n_workers, " core"))
power_start <- Sys.time()

# --- Fixed-n scenarios (3 effects x 3 heterogeneity = 9) ---
if (file.exists(file.path(csv_dir, "basket_power_fixed_n.csv"))) {
  cat("[BASKET] Power fixed-n gia' esistente, salto\n"); flush.console()
  bot$sendMessage(chat_id = chat_id, text = "[BASKET] Power fixed-n gia' fatto, salto")
} else {
scenarios_power <- list(
  list(name = "P1_effect08", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.8, 5), sd = 1, sc = "P1"),
  list(name = "P2_effect08", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.8, 5), sd = 1, sc = "P2"),
  list(name = "P3_effect08", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.8, 5), sd = 1, sc = "P3"),
  list(name = "P1_effect05", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.5, 5), sd = 1, sc = "P1"),
  list(name = "P2_effect05", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.5, 5), sd = 1, sc = "P2"),
  list(name = "P3_effect05", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.5, 5), sd = 1, sc = "P3"),
  list(name = "P1_effect02", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.2, 5), sd = 1, sc = "P1"),
  list(name = "P2_effect02", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.2, 5), sd = 1, sc = "P2"),
  list(name = "P3_effect02", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.2, 5), sd = 1, sc = "P3")
)

power_fixed <- map_dfr(scenarios_power, function(s) {
  bot$sendMessage(chat_id = chat_id,
    text = paste("[BASKET] Power:", s$name, "-", format(Sys.time(), "%H:%M:%S")))
  sim  <- run_simulation_par(n_sim, s$n_baskets, s$n, s$mu_treat, s$sd, s$sc)
  perf <- evaluate_performance(sim)
  perf %>% mutate(scenario = s$name, effect_size = s$mu_treat[1], heterogeneity = s$sc)
})

write_csv(power_fixed, file.path(csv_dir, "basket_power_fixed_n.csv"))
bot$sendMessage(chat_id = chat_id, text = "[BASKET] Power fixed-n salvato")
restart_workers()
} # end if !file.exists

# --- Power curves (3 heterogeneity x 3 effects x n_grid) ---
effects    <- list(rep(0.8, 5), rep(0.5, 5), rep(0.2, 5))
sc_codes   <- c("P1", "P2", "P3")
sc_labels  <- c("P1 - High heterogeneity",
                 "P2 - Low heterogeneity",
                 "P3 - Moderate heterogeneity")

if (file.exists(file.path(csv_dir, "basket_power_curves.csv"))) {
  cat("[BASKET] Power curves gia' esistente, salto\n"); flush.console()
  bot$sendMessage(chat_id = chat_id, text = "[BASKET] Power curves gia' fatto, salto")
} else {
power_curves_list <- list()
for (i in seq_along(sc_codes)) {
  for (eff in effects) {
    bot$sendMessage(chat_id = chat_id,
      text = paste("[BASKET] Power curve:", sc_codes[i],
                   "effect =", eff[1], "-", format(Sys.time(), "%H:%M:%S")))
    chunk <- map_dfr(n_grid, function(n_val) {
      sim   <- run_simulation_par(n_sim, 5, n_val, eff, 1, sc_codes[i])
      power <- evaluate_performance(sim) %>% pull(power) %>% mean(na.rm = TRUE)
      tibble(scenario = sc_labels[i], n = n_val,
             power = power, effect_size = eff[1])
    })
    power_curves_list <- c(power_curves_list, list(chunk))
    restart_workers()
  }
}
power_curves <- bind_rows(power_curves_list)

write_csv(power_curves, file.path(csv_dir, "basket_power_curves.csv"))
} # end if !file.exists

power_end <- Sys.time()
bot$sendMessage(chat_id = chat_id,
  text = paste0("[BASKET] POWER completata in ",
                round(difftime(power_end, power_start, units = "hours"), 2), " ore"))


##############################################################################
# ======================
# FDR SIMULATIONS
# ======================
##############################################################################

bot$sendMessage(chat_id = chat_id,
                text = paste0("[BASKET] Inizio FDR - ",
                              n_sim, " sim su ", n_workers, " core"))
fdr_start <- Sys.time()

# --- Fixed-n FDR scenarios (3 heterogeneity under H0) ---
if (file.exists(file.path(csv_dir, "basket_fdr_fixed_n.csv"))) {
  cat("[BASKET] FDR fixed-n gia' esistente, salto\n"); flush.console()
  bot$sendMessage(chat_id = chat_id, text = "[BASKET] FDR fixed-n gia' fatto, salto")
} else {
scenarios_fdr <- list(
  list(name = "P1_null", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.8, 5), mu_control = 0.8, sd = 1, sc = "P1"),
  list(name = "P2_null", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.8, 5), mu_control = 0.8, sd = 1, sc = "P2"),
  list(name = "P3_null", n_baskets = 5, n = c(17,15,18,12,23),
       mu_treat = rep(0.8, 5), mu_control = 0.8, sd = 1, sc = "P3")
)

fdr_fixed <- map_dfr(scenarios_fdr, function(s) {
  bot$sendMessage(chat_id = chat_id,
    text = paste("[BASKET] FDR:", s$name, "-", format(Sys.time(), "%H:%M:%S")))
  sim  <- run_simulation_fdr_par(n_sim, s$n_baskets, s$n,
                                  s$mu_treat, s$mu_control, s$sd, s$sc)
  perf <- evaluate_performance_fdr(sim)
  perf %>% mutate(scenario = s$name)
})

write_csv(fdr_fixed, file.path(csv_dir, "basket_fdr_fixed_n.csv"))
bot$sendMessage(chat_id = chat_id, text = "[BASKET] FDR fixed-n salvato")
restart_workers()
} # end if !file.exists

# --- FDR error curves (3 heterogeneity x n_grid) ---
no_effect <- rep(0.8, 5)

if (file.exists(file.path(csv_dir, "basket_fdr_curves.csv"))) {
  cat("[BASKET] FDR curves gia' esistente, salto\n"); flush.console()
  bot$sendMessage(chat_id = chat_id, text = "[BASKET] FDR curves gia' fatto, salto")
} else {
fdr_curves_list <- list()
for (sc in sc_codes) {
  bot$sendMessage(chat_id = chat_id,
    text = paste("[BASKET] FDR curve:", sc, "-", format(Sys.time(), "%H:%M:%S")))
  chunk <- map_dfr(n_grid, function(n_val) {
    sim  <- run_simulation_fdr_par(n_sim, 5, rep(n_val, 5),
                                    no_effect, 0.8, 1, sc)
    perf <- evaluate_performance_fdr(sim)
    tibble(scenario = sc, n = n_val, FDR = perf$FDR, FWER = perf$FWER)
  })
  fdr_curves_list <- c(fdr_curves_list, list(chunk))
  restart_workers()
}
fdr_curves <- bind_rows(fdr_curves_list)

write_csv(fdr_curves, file.path(csv_dir, "basket_fdr_curves.csv"))
} # end if !file.exists

fdr_end <- Sys.time()
bot$sendMessage(chat_id = chat_id,
  text = paste0("[BASKET] FDR completata in ",
                round(difftime(fdr_end, fdr_start, units = "hours"), 2), " ore"))


##############################################################################
# Chiusura
##############################################################################
plan(sequential)

bot$sendMessage(chat_id = chat_id,
  text = paste0("[BASKET] TUTTE LE SIMULAZIONI COMPLETATE\nTempo totale: ",
                round(difftime(fdr_end, power_start, units = "hours"), 2), " ore"))
