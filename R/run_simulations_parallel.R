##############################################################################
# SMART-Vent Parallel Simulation Runner (FDR + Power)
# Runs both simulation scenarios sequentially, each parallelized
# across 115 cores using future_lapply.
##############################################################################

# ===========================
# Libraries
# ===========================
library(rstanarm)
library(dplyr)
library(tidyverse)
library(rstan)
library(janitor)
library(gsDesign)
library(HDInterval)
library(future)
library(future.apply)
library(telegram.bot)

# ===========================
# Telegram setup
# ===========================
readRenviron(".env")
bot <- Bot(token = Sys.getenv("TELEGRAM_BOT_TOKEN"))
chat_id <- Sys.getenv("TELEGRAM_CHAT_ID")

# ===========================
# Parallel setup
# ===========================
options(mc.cores = 128)
rstan_options(auto_write = TRUE)
n_workers <- 115
plan(multisession, workers = n_workers)

# ===========================
# Settings
# ===========================
set.seed(123)
n_simulations <- 1000
n_sample <- seq(from = 200, to = 1200, by = 100)
subgroup_probs <- c("ARDS" = 0.1, "non-ARDS" = 0.9)
eps <- 0.95

# Gamma (response model) - shared between FDR and Power
gamma0_par <- 0.2
gamma_strategy_par <- c("A" = 0.1, "B" = 0.3)
gamma_subgroup_par <- c("ARDS" = -0.2, "non-ARDS" = 0.1)

# ===========================
# Helper functions
# ===========================

get_response_prob <- function(group, strategy) {
  base <- substr(strategy, 1, 1)
  linpred <- gamma0_par + gamma_strategy_par[base] + gamma_subgroup_par[group]
  return(1 / (1 + exp(-linpred)))
}

assign_by_blocks <- function(response,
                             block_size = 8,
                             labels = c("C", "D"),
                             n_labels_each = block_size / 2,
                             seed = NULL) {
  stopifnot(is.numeric(response) || is.logical(response))
  stopifnot(block_size %% length(labels) == 0,
            n_labels_each * length(labels) == block_size)
  if (!is.null(seed)) set.seed(seed)
  elig <- response != 1
  n_elig <- sum(elig)
  n_blocks <- ceiling(n_elig / block_size)
  make_block <- function() sample(rep(labels, each = n_labels_each), block_size)
  blocks <- replicate(n_blocks, make_block(), simplify = FALSE)
  treat_seq <- unlist(blocks, use.names = FALSE)
  treatment <- character(length(response))
  treatment[elig] <- treat_seq[seq_len(n_elig)]
  treatment[!elig] <- "Continue"
  return(treatment)
}

prob_all <- function(p_post, alfa = 0.5) {
  p_alfa <- p_post^alfa
  p_t <- sum(p_alfa)
  p <- p_alfa / p_t
  return(p)
}

stopping_function <- function(v1, stopper, v2, alfa) {
  if (length(v1) == 2 && (all(v1 == c(1, 0)) || all(v1 == c(0, 1)))) {
    return(v1)
  }
  nm <- names(v1)
  if (!is.null(nm) && length(nm) >= 2 && !is.na(stopper)) {
    if (identical(stopper, nm[1])) {
      res <- c(1, 0); names(res) <- nm; return(res)
    }
    if (identical(stopper, nm[2])) {
      res <- c(0, 1); names(res) <- nm; return(res)
    }
  }
  return(prob_all(v2, alfa = alfa))
}

winner <- function(x, eps) {
  idx <- which(x > eps)
  if (length(idx) == 1) {
    return(names(x)[idx])
  } else {
    return('NA')
  }
}

decision_function <- function(prob, res, esp) {
  if (all(prob == c(1, 0))) return(names(prob)[1])
  if (all(prob == c(0, 1))) return(names(prob)[2])
  return(winner(res, esp))
}

calc_best_prob <- function(df, cols, default_val = 0.5) {
  if (all(cols %in% colnames(df))) {
    m <- df[, cols, drop = FALSE]
    best <- apply(m, 1, function(row) cols[which.max(row)])
    tab <- table(factor(best, levels = cols))
    prob <- as.numeric(prop.table(tab))
    names(prob) <- cols
    list(best = best, prob = prob)
  } else {
    prob <- setNames(rep(default_val, length(cols)), cols)
    list(best = NULL, prob = prob)
  }
}

# ===========================
# Single simulation function
# ===========================
# sim_type: "fdr" or "power"
run_one_sim <- function(sim_id, n_t,
                        beta0, beta_strategy, beta_subgroup, beta_interaction,
                        model_formula_str, sim_type) {

  options(mc.cores = 1)

  batch_size <- 50
  n_batches <- (n_t - 100) / batch_size
  x <- gsDesign::gsDesign(k = n_batches + 1, test.type = 1, alpha = 0.05)
  bound <- c(pnorm(x$upper$bound), 1)

  # --- Local utility generators (depend on betas) ---
  gen_util <- function(strategy, group, stage1) {
    if (grepl("A", strategy)) base <- "A"
    if (grepl("B", strategy)) base <- "B"
    if (grepl("C", strategy)) base <- "C"
    if (grepl("D", strategy)) base <- "D"
    mu <- beta0 + beta_strategy[base] + beta_subgroup[group] + beta_interaction[[base]][group]
    return(rnorm(1, mu, 1))
  }

  gen_util_ref <- function(strategy, group, stage1) {
    if (grepl("A", strategy)) base <- "A"
    if (grepl("B", strategy)) base <- "B"
    if (grepl("C", strategy)) base <- "C"
    if (grepl("D", strategy)) base <- "D"
    mu <- beta0 + beta_strategy[base] + beta_subgroup[group] + beta_interaction[[base]][group]
    return(mu)
  }

  model_formula <- as.formula(model_formula_str)

  # --- Initial equal randomization ---
  p_r_ards <- c(A = 0.5, B = 0.5)
  p_r_non_ards <- c(A = 0.5, B = 0.5)
  p_non_r_ards <- c(C = 0.5, D = 0.5)
  p_non_r_non_ards <- c(C = 0.5, D = 0.5)

  trial_data <- data.frame()

  # --- First 100 patients ---
  subgroup_1 <- sample(c("ARDS", "non-ARDS"), size = 100, replace = TRUE, prob = subgroup_probs)
  stage1_1 <- sample(rep(c('A', 'B'), 50), size = 100, replace = FALSE)
  response_prob_1 <- mapply(get_response_prob, subgroup_1, stage1_1)
  response_1 <- rbinom(100, 1, prob = response_prob_1)
  stage2_1 <- assign_by_blocks(response_1)
  strategy_1 <- ifelse(stage2_1 == "Continue", stage1_1, stage2_1)

  outcome_1 <- mapply(gen_util, strategy_1, subgroup_1, stage1_1)
  outcome_ref_1 <- mapply(gen_util_ref, strategy_1, subgroup_1, stage1_1)

  trial_data <- data.frame(
    PatientID = 1:100,
    Subgroup = subgroup_1,
    Stage1 = stage1_1,
    Response = response_1,
    Stage2 = stage2_1,
    Strategy = strategy_1,
    Outcome = outcome_1,
    Outcome_ref = outcome_ref_1
  )

  # --- Initial model fit ---
  model <- stan_glmer(
    model_formula,
    data = trial_data,
    family = gaussian(),
    chains = 2, iter = 1000, refresh = 0
  )

  newdat <- expand.grid(
    Strategy = unique(trial_data$Strategy),
    Subgroup = unique(trial_data$Subgroup)
  )

  post <- posterior_linpred(model, transform = TRUE, newdata = newdat)
  colnames(post) <- paste(newdat$Strategy, newdat$Subgroup, sep = "_")
  post_db <- as.data.frame(post)

  diff_r_ards <- post_db$A_ARDS - post_db$B_ARDS
  diff_r_non_ards <- post_db$`A_non-ARDS` - post_db$`B_non-ARDS`
  diff_non_r_ards <- post_db$C_ARDS - post_db$D_ARDS
  diff_non_r_non_ards <- post_db$`C_non-ARDS` - post_db$`D_non-ARDS`

  res_r_ards <- calc_best_prob(post, c("A_ARDS", "B_ARDS"))
  prob_best_strategy_r_ards <- res_r_ards$prob

  res_r_non_ards <- calc_best_prob(post, c("A_non-ARDS", "B_non-ARDS"))
  prob_best_strategy_r_non_ards <- res_r_non_ards$prob

  res_non_r_ards <- calc_best_prob(post, c("C_ARDS", "D_ARDS"))
  prob_best_strategy_non_r_ards <- res_non_r_ards$prob

  res_non_r_non_ards <- calc_best_prob(post, c("C_non-ARDS", "D_non-ARDS"))
  prob_best_strategy_non_r_non_ards <- res_non_r_non_ards$prob

  alfa <- length(trial_data$PatientID) / (2 * n_t)

  # --- Initial stopping checks ---
  cred_r_ards <- hdi(diff_r_ards, credMass = bound[1])
  check_r_ards <- ifelse(cred_r_ards[2] < 0, 'B_ARDS',
                         ifelse(cred_r_ards[1] > 0, 'A_ARDS', 'No winner'))
  p_r_ards <- stopping_function(v1 = c('A_ARDS' = 0.5, 'B_ARDS' = 0.5),
                                stopper = check_r_ards,
                                v2 = prob_best_strategy_r_ards, alfa = alfa)

  cred_r_non_ards <- hdi(diff_r_non_ards, credMass = bound[1])
  check_r_non_ards <- ifelse(cred_r_non_ards[2] < 0, 'B_non-ARDS',
                             ifelse(cred_r_non_ards[1] > 0, 'A_non-ARDS', 'No winner'))
  p_r_non_ards <- stopping_function(v1 = c('A_non-ARDS' = 0.5, 'B_non-ARDS' = 0.5),
                                    stopper = check_r_non_ards,
                                    v2 = prob_best_strategy_r_non_ards, alfa = alfa)

  cred_non_r_ards <- hdi(diff_non_r_ards, credMass = bound[1])
  check_non_r_ards <- ifelse(cred_non_r_ards[2] < 0, 'D_ARDS',
                             ifelse(cred_non_r_ards[1] > 0, 'C_ARDS', 'No winner'))
  p_non_r_ards <- stopping_function(v1 = c('C_ARDS' = 0.5, 'D_ARDS' = 0.5),
                                    stopper = check_non_r_ards,
                                    v2 = prob_best_strategy_non_r_ards, alfa = alfa)

  cred_non_r_non_ards <- hdi(diff_non_r_non_ards, credMass = bound[1])
  check_non_r_non_ards <- ifelse(cred_non_r_non_ards[2] < 0, 'D_non-ARDS',
                                 ifelse(cred_non_r_non_ards[1] > 0, 'C_non-ARDS', 'No winner'))
  p_non_r_non_ards <- stopping_function(v1 = c('C_non-ARDS' = 0.5, 'D_non-ARDS' = 0.5),
                                        stopper = check_non_r_non_ards,
                                        v2 = prob_best_strategy_non_r_non_ards, alfa = alfa)

  # =============================
  # Adaptive trial loop
  # =============================
  for (batch in 1:n_batches) {

    subgroup <- sample(c("ARDS", "non-ARDS"), size = batch_size, replace = TRUE, prob = subgroup_probs)

    stage1 <- ifelse(subgroup == 'ARDS',
                     sample(c("A", "B"), size = batch_size, replace = TRUE, prob = p_r_ards),
                     sample(c("A", "B"), size = batch_size, replace = TRUE, prob = p_r_non_ards))

    response_prob <- mapply(get_response_prob, subgroup, stage1)
    response <- rbinom(batch_size, 1, prob = response_prob)

    stage2 <- ifelse(response == 1,
                     "Continue",
                     ifelse(subgroup == 'ARDS',
                            sample(c("C", "D"), size = batch_size, replace = TRUE, prob = p_non_r_ards),
                            sample(c("C", "D"), size = batch_size, replace = TRUE, prob = p_non_r_non_ards)))

    strategy <- ifelse(stage2 == "Continue", stage1, stage2)

    outcome <- mapply(gen_util, strategy, subgroup, stage1)
    outcome_ref <- mapply(gen_util_ref, strategy, subgroup, stage1)

    batch_data <- data.frame(
      PatientID = (nrow(trial_data) + 1):(nrow(trial_data) + batch_size),
      Subgroup = subgroup,
      Stage1 = stage1,
      Response = response,
      Stage2 = stage2,
      Strategy = strategy,
      Outcome = outcome,
      Outcome_ref = outcome_ref
    )
    trial_data <- rbind(trial_data, batch_data)

    model <- stan_glmer(
      model_formula,
      data = trial_data,
      family = gaussian(),
      chains = 2, iter = 1000, refresh = 0
    )

    newdat <- expand.grid(
      Strategy = unique(trial_data$Strategy),
      Subgroup = unique(trial_data$Subgroup)
    )

    post <- posterior_linpred(model, transform = TRUE, newdata = newdat)
    colnames(post) <- paste(newdat$Strategy, newdat$Subgroup, sep = "_")

    res_r_ards <- calc_best_prob(post, c("A_ARDS", "B_ARDS"))
    prob_best_strategy_r_ards <- res_r_ards$prob

    res_r_non_ards <- calc_best_prob(post, c("A_non-ARDS", "B_non-ARDS"))
    prob_best_strategy_r_non_ards <- res_r_non_ards$prob

    res_non_r_ards <- calc_best_prob(post, c("C_ARDS", "D_ARDS"))
    prob_best_strategy_non_r_ards <- res_non_r_ards$prob

    res_non_r_non_ards <- calc_best_prob(post, c("C_non-ARDS", "D_non-ARDS"))
    prob_best_strategy_non_r_non_ards <- res_non_r_non_ards$prob

    alfa <- length(trial_data$PatientID) / (2 * n_t)

    cred_r_ards <- hdi(diff_r_ards, credMass = bound[batch + 1])
    check_r_ards <- ifelse(cred_r_ards[2] < 0, 'B_ARDS',
                           ifelse(cred_r_ards[1] > 0, 'A_ARDS', 'No winner'))
    p_r_ards <- stopping_function(v1 = p_r_ards, stopper = check_r_ards,
                                  v2 = prob_best_strategy_r_ards, alfa = alfa)

    cred_r_non_ards <- hdi(diff_r_non_ards, credMass = bound[batch + 1])
    check_r_non_ards <- ifelse(cred_r_non_ards[2] < 0, 'B_non-ARDS',
                               ifelse(cred_r_non_ards[1] > 0, 'A_non-ARDS', 'No winner'))
    p_r_non_ards <- stopping_function(v1 = p_r_non_ards, stopper = check_r_non_ards,
                                      v2 = prob_best_strategy_r_non_ards, alfa = alfa)

    cred_non_r_ards <- hdi(diff_non_r_ards, credMass = bound[batch + 1])
    check_non_r_ards <- ifelse(cred_non_r_ards[2] < 0, 'D_ARDS',
                               ifelse(cred_non_r_ards[1] > 0, 'C_ARDS', 'No winner'))
    p_non_r_ards <- stopping_function(v1 = p_non_r_ards, stopper = check_non_r_ards,
                                      v2 = prob_best_strategy_non_r_ards, alfa = alfa)

    cred_non_r_non_ards <- hdi(diff_non_r_non_ards, credMass = bound[batch + 1])
    check_non_r_non_ards <- ifelse(cred_non_r_non_ards[2] < 0, 'D_non-ARDS',
                                   ifelse(cred_non_r_non_ards[1] > 0, 'C_non-ARDS', 'No winner'))
    p_non_r_non_ards <- stopping_function(v1 = p_non_r_non_ards, stopper = check_non_r_non_ards,
                                          v2 = prob_best_strategy_non_r_non_ards, alfa = alfa)
  }

  # =============================
  # Final decision
  # =============================
  w_r_ards <- decision_function(prob = p_r_ards, res = prob_best_strategy_r_ards, esp = eps)
  w_r_non_ards <- decision_function(prob = p_r_non_ards, res = prob_best_strategy_r_non_ards, esp = eps)
  w_non_r_ards <- decision_function(prob = p_non_r_ards, res = prob_best_strategy_non_r_ards, esp = eps)
  w_non_r_non_ards <- decision_function(prob = p_non_r_non_ards, res = prob_best_strategy_non_r_non_ards, esp = eps)

  if (sim_type == "fdr") {
    c_r_ards <- ifelse(w_r_ards == 'B_ARDS' | w_r_ards == 'A_ARDS', 1, 0)
    c_r_non_ards <- ifelse(w_r_non_ards == 'A_non-ARDS' | w_r_non_ards == 'B_non-ARDS', 1, 0)
    c_non_r_ards <- ifelse(w_non_r_ards == 'D_ARDS' | w_non_r_ards == 'C_ARDS', 1, 0)
    c_non_r_non_ards <- ifelse(w_non_r_non_ards == 'D_non-ARDS' | w_non_r_non_ards == 'C_non-ARDS', 1, 0)
  } else {
    c_r_ards <- ifelse(w_r_ards == 'B_ARDS', 1, 0)
    c_r_non_ards <- ifelse(w_r_non_ards == 'A_non-ARDS', 1, 0)
    c_non_r_ards <- ifelse(w_non_r_ards == 'D_ARDS', 1, 0)
    c_non_r_non_ards <- ifelse(w_non_r_non_ards == 'D_non-ARDS', 1, 0)
  }

  # =============================
  # Compute metrics
  # =============================
  diff_resp_ARDS <- mean(trial_data$Outcome_ref[trial_data$Strategy == 'A' & trial_data$Subgroup == 'ARDS']) -
    mean(trial_data$Outcome_ref[trial_data$Strategy == 'B' & trial_data$Subgroup == 'ARDS'])
  diff_resp_non_ARDS <- mean(trial_data$Outcome_ref[trial_data$Strategy == 'A' & trial_data$Subgroup == 'non-ARDS']) -
    mean(trial_data$Outcome_ref[trial_data$Strategy == 'B' & trial_data$Subgroup == 'non-ARDS'])
  diff_non_resp_ARDS <- mean(trial_data$Outcome_ref[trial_data$Strategy == 'C' & trial_data$Subgroup == 'ARDS']) -
    mean(trial_data$Outcome_ref[trial_data$Strategy == 'D' & trial_data$Subgroup == 'ARDS'])
  diff_non_resp_non_ARDS <- mean(trial_data$Outcome_ref[trial_data$Strategy == 'C' & trial_data$Subgroup == 'non-ARDS']) -
    mean(trial_data$Outcome_ref[trial_data$Strategy == 'D' & trial_data$Subgroup == 'non-ARDS'])

  tbl_sub_resp <- matrix(table(trial_data$Subgroup, trial_data$Response))
  tbl_strat <- matrix(table(trial_data$Strategy, trial_data$Subgroup, trial_data$Response))

  list(
    result_resp_ARDS = c_r_ards,
    result_resp_non_ARDS = c_r_non_ards,
    result_non_resp_ARDS = c_non_r_ards,
    result_non_resp_non_ARDS = c_non_r_non_ards,
    n_r_ARDS = tbl_sub_resp[3],
    n_non_r_ARDS = tbl_sub_resp[1],
    n_r_non_ARDS = tbl_sub_resp[4],
    n_non_r_non_ARDS = tbl_sub_resp[2],
    n_a_n_r_ARDS = tbl_strat[9],
    n_a_n_r_non_ARDS = tbl_strat[13],
    n_c_n_non_r_ARDS = tbl_strat[3],
    n_c_n_non_r_non_ARDS = tbl_strat[7],
    stop_r_ARDS = ifelse(p_r_ards[1] == 1 | p_r_ards[2] == 1, 1, 0),
    stop_r_non_ARDS = ifelse(p_r_non_ards[1] == 1 | p_r_non_ards[2] == 1, 1, 0),
    stop_non_r_ARDS = ifelse(p_non_r_ards[1] == 1 | p_non_r_ards[2] == 1, 1, 0),
    stop_non_r_non_ARDS = ifelse(p_non_r_non_ards[1] == 1 | p_non_r_non_ards[2] == 1, 1, 0),
    error_resp_ARDS = (mean(diff_r_ards) - diff_resp_ARDS) * 100 / diff_resp_ARDS,
    error_resp_non_ARDS = (mean(diff_r_non_ards) - diff_resp_non_ARDS) * 100 / diff_resp_non_ARDS,
    error_non_resp_ARDS = (mean(diff_non_r_ards) - diff_non_resp_ARDS) * 100 / diff_non_resp_ARDS,
    error_non_resp_non_ARDS = (mean(diff_non_r_non_ards) - diff_non_resp_non_ARDS) * 100 / diff_non_resp_non_ARDS,
    ic_size_resp_ARDS = abs(cred_r_ards[2] - cred_r_ards[1]),
    ic_size_resp_non_ARDS = abs(cred_r_non_ards[2] - cred_r_non_ards[1]),
    ic_size_non_resp_ARDS = abs(cred_non_r_ards[2] - cred_non_r_ards[1]),
    ic_size_non_resp_non_ARDS = abs(cred_non_r_non_ards[2] - cred_non_r_non_ards[1])
  )
}

# ===========================
# Parallel simulator wrapper
# ===========================
simulatore_parallel <- function(n_t, beta0, beta_strategy, beta_subgroup,
                                beta_interaction, model_formula_str, sim_type) {

  results_list <- future_lapply(1:n_simulations, function(i) {
    run_one_sim(i, n_t,
                beta0, beta_strategy, beta_subgroup, beta_interaction,
                model_formula_str, sim_type)
  }, future.seed = 123L)

  df <- bind_rows(results_list)

  power_res <- data.frame(
    power = c(mean(df$result_resp_ARDS) * 100,
              mean(df$result_resp_non_ARDS) * 100,
              mean(df$result_non_resp_ARDS) * 100,
              mean(df$result_non_resp_non_ARDS) * 100),
    sample_size = rep(n_t, 4),
    ic_size = c(mean(df$ic_size_resp_ARDS),
                mean(df$ic_size_resp_non_ARDS),
                mean(df$ic_size_non_resp_ARDS),
                mean(df$ic_size_non_resp_non_ARDS)),
    mape = c(mean(df$error_resp_ARDS),
             mean(df$error_non_resp_ARDS),
             mean(df$error_resp_non_ARDS),
             mean(df$error_non_resp_non_ARDS)),
    early_stop = c(mean(df$stop_r_ARDS) * 100,
                   mean(df$stop_non_r_ARDS) * 100,
                   mean(df$stop_r_non_ARDS) * 100,
                   mean(df$stop_non_r_non_ARDS) * 100),
    n = c(mean(df$n_r_ARDS),
          mean(df$n_r_non_ARDS),
          mean(df$n_non_r_ARDS),
          mean(df$n_non_r_non_ARDS)),
    n_strate = c(mean(df$n_a_n_r_ARDS),
                 mean(df$n_a_n_r_non_ARDS),
                 mean(df$n_c_n_non_r_ARDS),
                 mean(df$n_c_n_non_r_non_ARDS)),
    strategy = c('A', 'A', 'C', 'C'),
    group = c('Responder ARDS', 'Responder non_ARDS',
              'Non Responder ARDS', 'Non Responder non_ARDS')
  )
  return(power_res)
}


##############################################################################
# =====================
# FDR SCENARIO (alpha)
# =====================
##############################################################################

fdr_beta0 <- 0
fdr_beta_strategy <- c("A" = 0, "B" = 0, "C" = 0, "D" = 0)
fdr_beta_subgroup <- c("ARDS" = 0, "non-ARDS" = 0)
fdr_beta_interaction <- list(
  "A" = c("ARDS" = 0, "non-ARDS" = 0),
  "B" = c("ARDS" = 0, "non-ARDS" = 0),
  "C" = c("ARDS" = 0, "non-ARDS" = 0),
  "D" = c("ARDS" = 0, "non-ARDS" = 0)
)
fdr_formula <- "Outcome ~ 0 + Strategy + (0 + Strategy | Subgroup)"

bot$sendMessage(chat_id = chat_id,
                text = paste0("Inizio simulazione FDR (alpha) - ",
                              n_simulations, " sim x ", length(n_sample),
                              " scenari su ", n_workers, " core"))

fdr_start <- Sys.time()

fdr_results <- n_sample %>%
  map_df(~ {
    msg <- paste("FDR: inizio n =", .x, "-", format(Sys.time(), "%H:%M:%S"))
    bot$sendMessage(chat_id = chat_id, text = msg)
    res <- simulatore_parallel(.x, fdr_beta0, fdr_beta_strategy, fdr_beta_subgroup,
                               fdr_beta_interaction, fdr_formula, "fdr")
    msg <- paste("FDR: fine n =", .x, "-", format(Sys.time(), "%H:%M:%S"))
    bot$sendMessage(chat_id = chat_id, text = msg)
    res
  })

write_csv(fdr_results, "risultati_simulazione_alpha_bayes.csv")

fdr_end <- Sys.time()
bot$sendMessage(chat_id = chat_id,
                text = paste0("FDR completata in ",
                              round(difftime(fdr_end, fdr_start, units = "hours"), 2),
                              " ore"))


##############################################################################
# =====================
# POWER SCENARIO
# =====================
##############################################################################

pow_beta0 <- 0
pow_beta_strategy <- c("A" = 0, "B" = .8, "C" = 0, "D" = .5)
pow_beta_subgroup <- c("ARDS" = -1, "non-ARDS" = 0)
pow_beta_interaction <- list(
  "A" = c("ARDS" = 0, "non-ARDS" = .9),
  "B" = c("ARDS" = 0, "non-ARDS" = -.4),
  "C" = c("ARDS" = 0, "non-ARDS" = 0),
  "D" = c("ARDS" = 0, "non-ARDS" = -0.3)
)
pow_formula <- "Outcome ~ 0 + Strategy*Subgroup + (0 + Strategy | Subgroup)"

bot$sendMessage(chat_id = chat_id,
                text = paste0("Inizio simulazione POWER - ",
                              n_simulations, " sim x ", length(n_sample),
                              " scenari su ", n_workers, " core"))

pow_start <- Sys.time()

pow_results <- n_sample %>%
  map_df(~ {
    msg <- paste("POWER: inizio n =", .x, "-", format(Sys.time(), "%H:%M:%S"))
    bot$sendMessage(chat_id = chat_id, text = msg)
    res <- simulatore_parallel(.x, pow_beta0, pow_beta_strategy, pow_beta_subgroup,
                               pow_beta_interaction, pow_formula, "power")
    msg <- paste("POWER: fine n =", .x, "-", format(Sys.time(), "%H:%M:%S"))
    bot$sendMessage(chat_id = chat_id, text = msg)
    res
  })

write_csv(pow_results, "risultati_simulazione_power_bayes.csv")

pow_end <- Sys.time()
bot$sendMessage(chat_id = chat_id,
                text = paste0("POWER completata in ",
                              round(difftime(pow_end, pow_start, units = "hours"), 2),
                              " ore"))

##############################################################################
# Chiusura
##############################################################################
plan(sequential)

bot$sendMessage(chat_id = chat_id,
                text = paste0("TUTTE LE SIMULAZIONI COMPLETATE\n",
                              "Tempo totale: ",
                              round(difftime(pow_end, fdr_start, units = "hours"), 2),
                              " ore"))
