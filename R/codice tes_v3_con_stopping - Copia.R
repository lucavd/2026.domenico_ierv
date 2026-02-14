
# =========================================================
# SMART-Vent Trial Design (Bayesian Adaptive Platform + SMART)
# =========================================================
# - Bayesian adaptive platform trial with sequential randomization
# - Subgroups: ARDS vs non-ARDS
# - Response Adaptive Randomization (RAR): allocation probabilities updated dynamically
# - Borrowing across subgroups: hierarchical Bayesian model (partial pooling)

# ------------------------
# Stage 1: Initial Allocation
# ------------------------
# Patients randomized to one of two initial strategies:
#   Strategy A = Lung-Protective Low Tidal Volume (LTV) Ventilation
#   Strategy B = Personalized Pressure-Controlled Ventilation (PCV)
# Allocation starts 50/50, then adapts with RAR.

# ------------------------
# Response Assessment (Day 3)
# ------------------------
# Patients classified as:
#   - Responders   → continue initial treatment
#   - Non-responders → re-randomized to alternative strategies

# ------------------------
# Stage 2: Adaptive Re-Randomization for Non-Responders
# ------------------------
# Responders:
#   A → Continue A (LTV continued)
#   B → Continue B (PCV continued)
#
# Non-Responders → adaptive randomization between:
#   Strategy C = Rescue High PEEP Ventilation
#   Strategy D = Diaphragm-Sparing Support Ventilation (e.g. NAVA, PAV)

# ------------------------
# Primary Endpoint
# ------------------------
# Utility-weighted composite:
#   - Ventilator-free days
#   - Survival
#   - Post-discharge recovery (weighted by PROMIS-29 / EQ-5D preferences)

# ------------------------
# Adaptive Features
# ------------------------
# - RAR: allocation shifts toward promising strategies batch-by-batch
# - Dynamic Borrowing: hierarchical Bayesian model pools information across ARDS vs non-ARDS
# - Stopping Rules:
#     * Early superiority if P(best) > 0.95
#     * Drop arms if P(best) < 0.01



library(rstanarm)
library(dplyr)

library(rstan)
library(janitor)
library(gsDesign)
library(HDInterval)
library(future)
library(future.apply)
library(progressr)
library(tidyverse)
set.seed(123)



# Fondamentale: impedisce che ogni singolo worker cerchi di usare tutti i core
#options(mc.cores = 4) 
rstan_options(auto_write = TRUE)

# 2. IMPOSTA IL PIANO PARALLELO "SAFE"
# Con 16GB di RAM, 4 worker sono il punto ideale (4GB RAM ciascuno).
# Non salire sopra i 4 o rischi il crash della memoria.
plan(multisession, workers = 4)


# Settings
n_patients_total <- 200
batch_size <- 50
n_batches <- (n_patients_total-100) / batch_size
n_simulations<-1000
subgroup_probs <- c("ARDS" = 0.1, "non-ARDS" = 0.9)
eps<-0.95 # per la scelta di epslion vedere pagina 200 del pdf o 187 del testo adapt
eps1<-0.02
eps2<-0.05
n_sample<-seq(from = 200, to = 1200, by = 100)

# -----------------------------
# Define betas for outcome model
# -----------------------------
beta0 <- 0  # baseline utility
beta_strategy <- c("A" = 0, "B" = .8, "C" = 0, "D" = .5)   # main effect of strategies
beta_subgroup <- c("ARDS" = -1, "non-ARDS" = 0)           # subgroup penalty
beta_stage1 <- c("A" = 1.5, "B" = 2)
beta_interaction <- list(
  "A" = c("ARDS" = 0, "non-ARDS" = .9),
  "B" = c("ARDS" = 0, "non-ARDS" = -.4),
  "C" = c("ARDS" = 0, "non-ARDS" = 0),
  "D" = c("ARDS" = 0, "non-ARDS" = -0.3)
)



# -----------------------------
# Define gammas for response model
# -----------------------------
gamma0 <- 0.2  # baseline logit
gamma_strategy <- c("A" = 0.1, "B" = 0.3)     # only stage1 relevant
gamma_subgroup <- c("ARDS" = -0.2, "non-ARDS" = 0.1)

# Function: utility outcome generation
generate_utility <- function(strategy, group, stage1, sampling = TRUE) {
  
  # Logica per determinare la base
  if (grepl("A", strategy)) base <- "A"
  if (grepl("B", strategy)) base <- "B"
  if (grepl("C", strategy)) base <- "C"
  if (grepl("D", strategy)) base <- "D"
  
  # Calcolo del valore atteso (mu)
  mu <- beta0 + beta_strategy[base] + beta_subgroup[group] + beta_interaction[[base]][group]
  
  # Controllo condizionale per il return
  if (sampling) {
    # Se TRUE (default): restituisce il valore campionato con rumore (SD=6)
    return(rnorm(1, mu, 1)) 
  } else {
    # Se FALSE: restituisce il valore mu pulito deterministico
    return(mu)
  }
}

generate_utility_2 <- function(strategy, group, stage1, sampling = F) {
  
  # Logica per determinare la base
  if (grepl("A", strategy)) base <- "A"
  if (grepl("B", strategy)) base <- "B"
  if (grepl("C", strategy)) base <- "C"
  if (grepl("D", strategy)) base <- "D"
  
  # Calcolo del valore atteso (mu)
  mu <- beta0 + beta_strategy[base] + beta_subgroup[group] + beta_interaction[[base]][group]
  
  # Controllo condizionale per il return
  if (sampling) {
    # Se TRUE (default): restituisce il valore campionato con rumore (SD=6)
    return(rnorm(1, mu, 1)) 
  } else {
    # Se FALSE: restituisce il valore mu pulito deterministico
    return(mu)
  }
}


generate_utility_tester <- function(strategy, group) {
  # parse strategy (e.g., "A → C" → take last letter if switched, else Continue)
  if (grepl("A", strategy)) base <- "A"
  if (grepl("B", strategy)) base <- "B"
  if (grepl("C", strategy)) base <- "C"
  if (grepl("D", strategy)) base <- "D"
  
  mu <- beta0 + beta_strategy[base] + beta_subgroup[group] + beta_interaction[[base]][group]
  #return(rnorm(1, mu, 10)) # noise SD=10
  return(mu)
}
test_data<-data.frame(
  group=rep(c('ARDS','non-ARDS'),each=4),
  strategy=c('A','B','C','D',
             'A','B','C','D'))

test_data$outcome<-mapply(generate_utility_tester,test_data$strategy, test_data$group)
test_data$summary<-paste(test_data$group,test_data$strategy)
test_data
table(test_data$outcome,test_data$strategy)

# Function: subgroup-specific response probabilities (logit model)
get_response_prob <- function(group, strategy) {
  # only Stage1 strategy matters for early response
  base <- substr(strategy, 1, 1)  # "A" or "B"
  linpred <- gamma0 + gamma_strategy[base] + gamma_subgroup[group]
  return(1 / (1 + exp(-linpred)))  # inverse logit
}

assign_by_blocks <- function(response,
                             block_size = 8,
                             labels = c("C", "D"),
                             n_labels_each = block_size / 2,
                             seed = NULL) {
  # Controlli di consistenza
  stopifnot(is.numeric(response) || is.logical(response))
  stopifnot(block_size %% length(labels) == 0,
            n_labels_each * length(labels) == block_size)
  
  if (!is.null(seed)) set.seed(seed)
  
  # Identifica posizioni eleggibili (response != 1)
  elig <- response != 1
  n_elig <- sum(elig)
  
  # Numero di blocchi necessari
  n_blocks <- ceiling(n_elig / block_size)
  
  # Genera blocchi bilanciati e randomizzati
  make_block <- function() sample(rep(labels, each = n_labels_each), block_size)
  blocks <- replicate(n_blocks, make_block(), simplify = FALSE)
  treat_seq <- unlist(blocks, use.names = FALSE)
  
  # Crea vettore di output
  treatment <- character(length(response))
  treatment[elig] <- treat_seq[seq_len(n_elig)]
  treatment[!elig] <- "Continue"
  
  return(treatment)
}
prob_all<-function(p_post,alfa=0.5){
  p_alfa<-p_post^alfa
  p_t<-sum(p_alfa)
  p<-p_alfa/p_t
  return(p)
}#funzione che aggiorna la RAR in modo proporzionale in base ad un parametro alfa

stopping_function <- function(v1, stopper, v2, alfa) {
  # v1: vettore numerico di lunghezza 2 (con o senza names)
  # stopper: stringa da confrontare con i nomi di v1
  # v2: vettore (input per prob_all)
  # alfa: parametro passato a prob_all
  
  ## 1. Se v1 è (1,0) o (0,1) → restituisci v1 così com'è
  if (length(v1) == 2 && (all(v1 == c(1, 0)) || all(v1 == c(0, 1)))) {
    return(v1)
  }
  
  ## 2. Controlla se "stopper" coincide con uno dei names(v1)
  nm <- names(v1)
  
  if (!is.null(nm) && length(nm) >= 2 && !is.na(stopper)) {
    # se coincide col nome del primo elemento → c(1,0)
    if (identical(stopper, nm[1])) {
      res <- c(1, 0)
      names(res) <- nm
      return(res)
    }
    # se coincide col nome del secondo elemento → c(0,1)
    if (identical(stopper, nm[2])) {
      res <- c(0, 1)
      names(res) <- nm
      return(res)
    }
  }
  
  ## 3. Nessuna delle condizioni precedenti → applica prob_all a v2 e alfa
  return(prob_all(v2, alfa = alfa))
}

winner <- function(x, eps) {
  # x: vettore numerico con nomi, es: c(A_ARDS = 0.029, B_ARDS = 0.971)
  # eps: soglia
  
  # quale valore supera eps?
  idx <- which(x > eps)
  
  if (length(idx) == 1) {
    # ritorna il nome della colonna vincente
    return(names(x)[idx])
  } else {
    # nessuno supera eps, oppure (impossibile nel tuo caso) più di uno
    return('NA')
  }
}

decision_function <- function(prob, res, esp) {
  # prob: vettore numerico di lunghezza 2 con names
  # res: oggetto passato alla funzione winner()
  # esp: parametro numerico per winner()
  
  ## 1. Se prob è c(1,0) → restituisci il nome del primo elemento
  if (all(prob == c(1, 0))) {
    return(names(prob)[1])
  }
  
  ## 2. Se prob è c(0,1) → restituisci il nome del secondo elemento
  if (all(prob == c(0, 1))) {
    return(names(prob)[2])
  }
  
  ## 3. Nessuna delle condizioni precedenti → applica winner(res, esp)
  return(winner(res, esp))
}


calc_best_prob <- function(df, cols, default_val = 0.5) {
  if (all(cols %in% colnames(df))) {
    m <- df[, cols, drop = FALSE]
    best <- apply(m, 1, function(row) cols[which.max(row)])
    # garantisce entrambe le categorie in output
    tab <- table(factor(best, levels = cols))
    prob <- as.numeric(prop.table(tab))
    names(prob) <- cols
    list(best = best, prob = prob)
  } else {
    prob <- setNames(rep(default_val, length(cols)), cols)
    list(best = NULL, prob = prob)
  }
}# questa funzione estrae dal vettore post la 
#probabilità di essere il migliore, in caso alcuni gruppi manchino restiruisce un vettore c(0.5,0.5)



# storage -----------------------------------------------------------------
n_r_ARDS<-numeric(n_simulations)
n_non_r_ARDS<-numeric(n_simulations)
n_r_non_ARDS<-numeric(n_simulations)
n_non_r_non_ARDS<-numeric(n_simulations)

n_a_n_r_ARDS<-numeric(n_simulations)
n_a_n_r_non_ARDS<-numeric(n_simulations)
n_c_n_non_r_ARDS<-numeric(n_simulations)
n_c_n_non_r_non_ARDS<-numeric(n_simulations)


stop_r_ARDS<-numeric(n_simulations)
stop_r_non_ARDS<-numeric(n_simulations)
stop_non_r_ARDS<-numeric(n_simulations)
stop_non_r_non_ARDS<-numeric(n_simulations)

error_resp_ARDS<-numeric(n_simulations)
error_resp_non_ARDS<-numeric(n_simulations)
error_non_resp_ARDS<-numeric(n_simulations)
error_non_resp_non_ARDS<-numeric(n_simulations)

results_non_resp_non_ARDS <- numeric(n_simulations)
results_non_resp_ARDS <- numeric(n_simulations)
results_resp_non_ARDS <- numeric(n_simulations)
results_resp_ARDS<- numeric(n_simulations)

ic_size_resp_ARDS<-numeric(n_simulations)
ic_size_resp_non_ARDS<-numeric(n_simulations)
ic_size_non_resp_ARDS<-numeric(n_simulations)
ic_size_non_resp_non_ARDS<-numeric(n_simulations)

simulatore<-function(n_t){
  n_r_ARDS<-numeric(n_simulations)
  n_non_r_ARDS<-numeric(n_simulations)
  n_r_non_ARDS<-numeric(n_simulations)
  n_non_r_non_ARDS<-numeric(n_simulations)
  
  n_a_n_r_ARDS<-numeric(n_simulations)
  n_a_n_r_non_ARDS<-numeric(n_simulations)
  n_c_n_non_r_ARDS<-numeric(n_simulations)
  n_c_n_non_r_non_ARDS<-numeric(n_simulations)
  
  
  stop_r_ARDS<-numeric(n_simulations)
  stop_r_non_ARDS<-numeric(n_simulations)
  stop_non_r_ARDS<-numeric(n_simulations)
  stop_non_r_non_ARDS<-numeric(n_simulations)
  
  error_resp_ARDS<-numeric(n_simulations)
  error_resp_non_ARDS<-numeric(n_simulations)
  error_non_resp_ARDS<-numeric(n_simulations)
  error_non_resp_non_ARDS<-numeric(n_simulations)
  
  results_non_resp_non_ARDS <- numeric(n_simulations)
  results_non_resp_ARDS <- numeric(n_simulations)
  results_resp_non_ARDS <- numeric(n_simulations)
  results_resp_ARDS<- numeric(n_simulations)
  
  ic_size_resp_ARDS<-numeric(n_simulations)
  ic_size_resp_non_ARDS<-numeric(n_simulations)
  ic_size_non_resp_ARDS<-numeric(n_simulations)
  ic_size_non_resp_non_ARDS<-numeric(n_simulations)
  
  for (i in 1:n_simulations) {
    
    n_patients_total <- n_t
    batch_size <- 50
    n_batches <- (n_patients_total -100)/ batch_size
    x <- gsDesign::gsDesign(k = n_batches+1, test.type = 1,alpha = 0.05)
    bound<-c(pnorm(x$upper$bound),1)
    


# Initial equal randomization
p_r_ards<-c(A=0.5,B=0.5)
p_r_non_ards<-c(A=0.5,B=0.5)

p_non_r_ards<-c(C=0.5,D=0.5)

p_non_r_non_ards<-c(C=0.5,D=0.5)

stage1_probs <- c(A = 0.5, B = 0.5)
stage2_probs <- c(C = 0.5, D = 0.5)


# -----------------------------
# Storage
# -----------------------------
trial_data <- data.frame()


# Credo che la function debba partire da qui ------------------------------
# Simulate subgroup membership primi 100
subgroup_1 <- sample(c("ARDS", "non-ARDS"), size = 100, replace = TRUE, prob = subgroup_probs) 
# Stage 1 primi 100 allocation (RAR)
stage1_1 <- sample(rep(c('A','B'),50), size = 100, replace = F)
# Simulate response primi 100
response_prob_1 <- mapply(get_response_prob, subgroup_1, stage1_1)
response_1 <- rbinom(100, 1, prob = response_prob_1)

# Stage 2 allocation (RAR for non-responders)
random_block<-unlist(replicate(13, sample(rep(c('C','D'), 4), 8), simplify = FALSE))
stage2_1 <- assign_by_blocks(response_1)

# Final strategy label

strategy_1 <- ifelse(stage2_1 == "Continue", stage1_1,stage2_1)


outcome_1 <- mapply(generate_utility, strategy_1, subgroup_1,stage1_1)
outcome_ref_1<-mapply(generate_utility_2, strategy_1, subgroup_1,stage1_1)

trial_data<- data.frame(
  PatientID = c(1:100),
  Subgroup = subgroup_1,
  Stage1 = stage1_1,
  Response = response_1,
  Stage2 = stage2_1,
  Strategy = strategy_1,
  Outcome = outcome_1,
  Outcome_ref=outcome_ref_1
)

model <- stan_glmer(
  Outcome ~ 0 + Strategy*Subgroup + (0 + Strategy | Subgroup),
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
post_db<-as.data.frame(post)
diff_r_ards<-post_db$A_ARDS-post_db$B_ARDS
diff_r_non_ards<-post_db$`A_non-ARDS`-post_db$`B_non-ARDS`

diff_non_r_ards<-post_db$C_ARDS-post_db$D_ARDS
diff_non_r_non_ards<-post_db$`C_non-ARDS`-post_db$`D_non-ARDS`


# --- R ARDS ---
res_r_ards <- calc_best_prob(post, c("A_ARDS", "B_ARDS"))
best_strategy_r_ards <- res_r_ards$best
prob_best_strategy_r_ards <- res_r_ards$prob

# --- R non-ARDS ---
res_r_non_ards <- calc_best_prob(post, c("A_non-ARDS", "B_non-ARDS"))
best_strategy_r_non_ards <- res_r_non_ards$best
prob_best_strategy_r_non_ards <- res_r_non_ards$prob

# --- non-R ARDS ---
res_non_r_ards <- calc_best_prob(post, c("C_ARDS","D_ARDS"))
best_strategy_non_r_ards <- res_non_r_ards$best
prob_best_strategy_non_r_ards <- res_non_r_ards$prob


# --- non-R non-ARDS ---
res_non_r_non_ards <- calc_best_prob(post, c("C_non-ARDS","D_non-ARDS"))
best_strategy_non_r_non_ards <- res_non_r_non_ards$best
prob_best_strategy_non_r_non_ards <- res_non_r_non_ards$prob
 

# Compute best strategy within each subgroup
best_strategy <- apply(post, 1, function(row) names(which.max(row)))
prob_best <- prop.table(table(best_strategy))


alfa<-length(trial_data$PatientID)/(2*n_patients_total)



# otteniamo le probabilità di allocazione per i responder ards ------------

cred_r_ards<-hdi(diff_r_ards,credMass = bound[1])
check_r_ards<-ifelse(cred_r_ards[2]<0,'B_ARDS',
       ifelse(cred_r_ards[1]>0,'A_ARDS','No winner'))
p_r_ards<-stopping_function(v1=c('A_ARDS'=0.5,'B_ARDS'=0.5),stopper = check_r_ards,
                  v2=prob_best_strategy_r_ards,alfa=alfa)

# otteniamo le probabilità di allocazione per i responder non ards ------------

cred_r_non_ards<-hdi(diff_r_non_ards,credMass = bound[1])
check_r_non_ards<-ifelse(cred_r_non_ards[2]<0,'B_non-ARDS',
                     ifelse(cred_r_non_ards[1]>0,'A_non-ARDS','No winner'))
p_r_non_ards<-stopping_function(v1=c('A_non-ARDS'=0.5,'B_non-ARDS'=0.5),stopper = check_r_non_ards,
                            v2=prob_best_strategy_r_non_ards,alfa=alfa)

# otteniamo le probabilità di allocazione per i non responder ards ------------

cred_non_r_ards<-hdi(diff_non_r_ards,credMass = bound[1])
check_non_r_ards<-ifelse(cred_non_r_ards[2]<0,'D_ARDS',
                         ifelse(cred_non_r_ards[1]>0,'C_ARDS','No winner'))
p_non_r_ards<-stopping_function(v1=c('C_ARDS'=0.5,'D_ARDS'=0.5),stopper = check_non_r_ards,
                                v2=prob_best_strategy_non_r_ards,alfa=alfa)
# otteniamo le probabilità di allocazione per i non responder non ards ------------

cred_non_r_non_ards<-hdi(diff_non_r_non_ards,credMass = bound[1])
check_non_r_non_ards<-ifelse(cred_non_r_non_ards[2]<0,'D_non-ARDS',
                         ifelse(cred_non_r_non_ards[1]>0,'C_non-ARDS','No winner'))
p_non_r_non_ards<-stopping_function(v1=c('C_non-ARDS'=0.5,'D_non-ARDS'=0.5),stopper = check_non_r_non_ards,
                                v2=prob_best_strategy_non_r_non_ards,alfa=alfa)


# -----------------------------
# Adaptive trial loop
# -----------------------------
for (batch in 1:n_batches) {
  

  # Simulate subgroup membership
  subgroup <- sample(c("ARDS", "non-ARDS"), size = batch_size, replace = TRUE, prob = subgroup_probs)
  
  # Stage 1 allocation (RAR)
  stage1 <- ifelse(subgroup=='ARDS',
                   sample(c("A", "B"), size = batch_size, replace = TRUE, prob = p_r_ards),
                   sample(c("A", "B"), size = batch_size, replace = TRUE, prob = p_r_non_ards))

  
  # Simulate response
  response_prob <- mapply(get_response_prob, subgroup, stage1)
  response <- rbinom(batch_size, 1, prob = response_prob)
  
  # Stage 2 allocation (RAR for non-responders)
  stage2 <- ifelse(response == 1,
                   "Continue",
                   ifelse(subgroup=='ARDS',
                          sample(c("C", "D"), size = batch_size, replace = TRUE, prob = p_non_r_ards),
                          sample(c("C", "D"), size = batch_size, replace = TRUE, prob = p_non_r_non_ards)
                          )
                   )
  # Final strategy label
  strategy <- ifelse(stage2 == "Continue",stage1, stage2)
  
  # Simulate outcomes
  outcome <- mapply(generate_utility, strategy, subgroup,stage1)
  outcome_ref<-mapply(generate_utility_2, strategy, subgroup,stage1)
  
  # Save batch data
  batch_data <- data.frame(
    PatientID = (nrow(trial_data) + 1):(nrow(trial_data) + batch_size),
    Subgroup = subgroup,
    Stage1 = stage1,
    Response = response,
    Stage2 = stage2,
    Strategy = strategy,
    Outcome = outcome,
    Outcome_ref=outcome_ref
  )
  trial_data <- rbind(trial_data, batch_data)
  
  # Bayesian hierarchical update (after ≥ 2 batches)
    
    # Hierarchical model with borrowing across subgroups
    model <- stan_glmer(
      Outcome ~ 0 + Strategy*Subgroup + (0 + Strategy | Subgroup),
      data = trial_data,
      family = gaussian(),
      chains = 2, iter = 1000, refresh = 0
    )
    
    # Posterior predictions for each Strategy × Subgroup
    newdat <- expand.grid(
      Strategy = unique(trial_data$Strategy),
      Subgroup = unique(trial_data$Subgroup)
    )
    
    post <- posterior_linpred(model, transform = TRUE, newdata = newdat)
    colnames(post) <- paste(newdat$Strategy, newdat$Subgroup, sep = "_")
    
    # Compute best strategy within each subgroup
    best_strategy <- apply(post, 1, function(row) names(which.max(row)))
    prob_best <- prop.table(table(best_strategy))
    
    
    # Compute the best strategy in subgroup -----------------------------------
    # --- R ARDS ---
    res_r_ards <- calc_best_prob(post, c("A_ARDS", "B_ARDS"))
    best_strategy_r_ards <- res_r_ards$best
    prob_best_strategy_r_ards <- res_r_ards$prob
    
    # --- R non-ARDS ---
    res_r_non_ards <- calc_best_prob(post, c("A_non-ARDS", "B_non-ARDS"))
    best_strategy_r_non_ards <- res_r_non_ards$best
    prob_best_strategy_r_non_ards <- res_r_non_ards$prob
    
    # --- non-R ARDS ---
    res_non_r_ards <- calc_best_prob(post, c("C_ARDS","D_ARDS"))
    best_strategy_non_r_ards <- res_non_r_ards$best
    prob_best_strategy_non_r_ards <- res_non_r_ards$prob
   
    # --- non-R non-ARDS ---
    res_non_r_non_ards <- calc_best_prob(post, c("C_non-ARDS","D_non-ARDS"))
    best_strategy_non_r_non_ards <- res_non_r_non_ards$best
    prob_best_strategy_non_r_non_ards <- res_non_r_non_ards$prob
    
    
    
    alfa<-length(trial_data$PatientID)/(2*n_patients_total)
    
    # otteniamo le probabilità di allocazione per i responder ards ------------
    
    cred_r_ards<-hdi(diff_r_ards,credMass = bound[batch+1])
    check_r_ards<-ifelse(cred_r_ards[2]<0,'B_ARDS',
                         ifelse(cred_r_ards[1]>0,'A_ARDS','No winner'))
    p_r_ards<-stopping_function(v1=p_r_ards,stopper = check_r_ards,
                                v2=prob_best_strategy_r_ards,alfa=alfa)
    
    # otteniamo le probabilità di allocazione per i responder non ards ------------
    
    cred_r_non_ards<-hdi(diff_r_non_ards,credMass = bound[batch+1])
    check_r_non_ards<-ifelse(cred_r_non_ards[2]<0,'B_non-ARDS',
                             ifelse(cred_r_non_ards[1]>0,'A_non-ARDS','No winner'))
    p_r_non_ards<-stopping_function(v1=p_r_non_ards,stopper = check_r_non_ards,
                                    v2=prob_best_strategy_r_non_ards,alfa=alfa)
    
    # otteniamo le probabilità di allocazione per i non responder ards ------------
    
    cred_non_r_ards<-hdi(diff_non_r_ards,credMass = bound[batch+1])
    check_non_r_ards<-ifelse(cred_non_r_ards[2]<0,'D_ARDS',
                             ifelse(cred_non_r_ards[1]>0,'C_ARDS','No winner'))
    p_non_r_ards<-stopping_function(v1=p_non_r_ards,stopper = check_non_r_ards,
                                    v2=prob_best_strategy_non_r_ards,alfa=alfa)
    # otteniamo le probabilità di allocazione per i non responder non ards ------------
    
    cred_non_r_non_ards<-hdi(diff_non_r_non_ards,credMass = bound[batch+1])
    check_non_r_non_ards<-ifelse(cred_non_r_non_ards[2]<0,'D_non-ARDS',
                                 ifelse(cred_non_r_non_ards[1]>0,'C_non-ARDS','No winner'))
    p_non_r_non_ards<-stopping_function(v1=p_non_r_non_ards,stopper = check_non_r_non_ards,
                                        v2=prob_best_strategy_non_r_non_ards,alfa=alfa)
    

}



w_r_ards<-decision_function(prob =p_r_ards,res = prob_best_strategy_r_ards,esp = eps)
c_r_ards<-ifelse(w_r_ards=='B_ARDS',1,0)

w_r_non_ards<-decision_function(prob =p_r_non_ards,res = prob_best_strategy_r_non_ards,esp = eps)
c_r_non_ards<-ifelse(w_r_non_ards=='A_non-ARDS',1,0)

#w_r_ards<-winner(prob_best_strategy_r_ards,eps = eps)
#c_r_ards<-ifelse(w_r_ards=='B_ARDS',1,0)
#w_r_non_ards<-winner(prob_best_strategy_r_non_ards,eps=eps)
#c_r_non_ards<-ifelse(w_r_non_ards=='A_non-ARDS',1,0)

w_non_r_ards<-decision_function(prob =p_non_r_ards,res = prob_best_strategy_non_r_ards,esp = eps)
c_non_r_ards<-ifelse(w_non_r_ards=='D_ARDS',1,0)

w_non_r_non_ards<-decision_function(prob =p_non_r_non_ards,res = prob_best_strategy_non_r_non_ards,esp = eps )
c_non_r_non_ards<-ifelse(w_non_r_non_ards=='D_non-ARDS',1,0)



#w_non_r_ards<-winner(prob_best_strategy_non_r_ards,eps=eps)
#c_non_r_ards<-ifelse(w_non_r_ards=='D_ARDS',1,0)
#w_non_r_non_ards<-winner(prob_best_strategy_non_r_non_ards,eps=eps)
#c_non_r_non_ards<-ifelse(w_non_r_non_ards=='D_non-ARDS',1,0)

diff_resp_ARDS<-mean(trial_data$Outcome_ref[trial_data$Strategy=='A'&trial_data$Subgroup=='ARDS'])-
  mean(trial_data$Outcome_ref[trial_data$Strategy=='B'&trial_data$Subgroup=='ARDS'])
diff_resp_non_ARDS<-mean(trial_data$Outcome_ref[trial_data$Strategy=='A'&trial_data$Subgroup=='non-ARDS'])-
  mean(trial_data$Outcome_ref[trial_data$Strategy=='B'&trial_data$Subgroup=='non-ARDS'])

diff_non_resp_ARDS<-mean(trial_data$Outcome_ref[trial_data$Strategy=='C'&trial_data$Subgroup=='ARDS'])-
  mean(trial_data$Outcome_ref[trial_data$Strategy=='D'&trial_data$Subgroup=='ARDS'])
diff_non_resp_non_ARDS<-mean(trial_data$Outcome_ref[trial_data$Strategy=='C'&trial_data$Subgroup=='non-ARDS'])-
  mean(trial_data$Outcome_ref[trial_data$Strategy=='D'&trial_data$Subgroup=='non-ARDS'])

# Risultati simulazione iesima --------------------------------------------


n_r_ARDS[i]<-matrix(table(trial_data$Subgroup,trial_data$Response))[3]
n_non_r_ARDS[i]<-matrix(table(trial_data$Subgroup,trial_data$Response))[1]
n_r_non_ARDS[i]<-matrix(table(trial_data$Subgroup,trial_data$Response))[4]
n_non_r_non_ARDS[i]<-matrix(table(trial_data$Subgroup,trial_data$Response))[2]

n_a_n_r_ARDS[i]<-matrix(table(trial_data$Strategy,trial_data$Subgroup,trial_data$Response))[9]
n_a_n_r_non_ARDS[i]<-matrix(table(trial_data$Strategy,trial_data$Subgroup,trial_data$Response))[13]
n_c_n_non_r_ARDS[i]<-matrix(table(trial_data$Strategy,trial_data$Subgroup,trial_data$Response))[3]
n_c_n_non_r_non_ARDS[i]<-matrix(table(trial_data$Strategy,trial_data$Subgroup,trial_data$Response))[7]


stop_r_ARDS[i]<-ifelse(p_r_ards[1]==1|p_r_ards[2]==1,1,0)
stop_r_non_ARDS[i]<-ifelse(p_r_non_ards[1]==1|p_r_non_ards[2]==1,1,0)
stop_non_r_ARDS[i]<-ifelse(p_non_r_ards[1]==1|p_non_r_ards[2]==1,1,0)
stop_non_r_non_ARDS[i]<-ifelse(p_non_r_non_ards[1]==1|p_non_r_non_ards[2]==1,1,0)

error_resp_ARDS[i]<-(mean(diff_r_ards)-diff_resp_ARDS)*100/diff_resp_ARDS
error_resp_non_ARDS[i]<-(mean(diff_r_non_ards)-diff_resp_non_ARDS)*100/diff_resp_non_ARDS
error_non_resp_ARDS[i]<-(mean(diff_non_r_ards)-diff_non_resp_ARDS)*100/diff_non_resp_ARDS
error_non_resp_non_ARDS[i]<-(mean(diff_non_r_non_ards)-diff_non_resp_non_ARDS)*100/diff_non_resp_non_ARDS

ic_size_resp_ARDS[i]<-abs(cred_r_ards[2]-cred_r_ards[1])
ic_size_resp_non_ARDS[i]<-abs(cred_r_non_ards[2]-cred_r_non_ards[1])
ic_size_non_resp_ARDS[i]<-abs(cred_non_r_ards[2]-cred_non_r_ards[1])
ic_size_non_resp_non_ARDS[i]<-abs(cred_non_r_non_ards[2]-cred_non_r_non_ards[1])

results_resp_ARDS[i] <- c_r_ards
results_resp_non_ARDS[i] <- c_r_non_ards
results_non_resp_ARDS[i] <- c_non_r_ards
results_non_resp_non_ARDS[i]<-c_non_r_non_ards
}



power_res<-data.frame(
  power=c(sum(results_resp_ARDS)/length(results_resp_ARDS)*100,
          sum(results_resp_non_ARDS)/length(results_resp_non_ARDS)*100,
          sum(results_non_resp_ARDS)/length(results_non_resp_ARDS)*100,
          sum(results_non_resp_non_ARDS)/length(results_non_resp_non_ARDS)*100),
  sample_size=rep(n_t,4),
  ic_size=c(mean(ic_size_resp_ARDS),
             mean(ic_size_resp_non_ARDS),
             mean(ic_size_non_resp_ARDS),
             mean(ic_size_non_resp_non_ARDS)
  ),
  mape=c(mean(error_resp_ARDS),
         mean(error_non_resp_ARDS),
         mean(error_resp_non_ARDS),
         mean(error_non_resp_non_ARDS)),
  early_stop=c(sum(stop_r_ARDS)/length(stop_r_ARDS)*100,
               sum(stop_non_r_ARDS)/length(stop_non_r_ARDS)*100,
               sum(stop_r_non_ARDS)/length(stop_r_non_ARDS)*100,
               sum(stop_non_r_non_ARDS)/length(stop_non_r_non_ARDS)*100
    
  ),
  n=c(mean(n_r_ARDS),
            mean(n_r_non_ARDS),
            mean(n_non_r_ARDS),
            mean(n_non_r_non_ARDS)),
  n_strate=c(mean(n_a_n_r_ARDS),
      mean(n_a_n_r_non_ARDS),
      mean(n_c_n_non_r_ARDS),
      mean(n_c_n_non_r_non_ARDS)),
  strategy=c('A','A','C','C'),
  group=c('Responder ARDS','Responder non_ARDS','Non Responder ARDS','Non Responder non_ARDS')
)
return(power_res)
}
# Definiamo i campioni da testare


# Esecuzione iterativa
final_results <- n_sample %>% 
  map_df(~ {
    message(paste("Inizio simulazione per n =", .x))
    simulatore(.x)
  })

# Visualizza i risultati finali
print(final_results)
df_finale<-final_results
write_csv(df_finale, "risultati_simulazione_power_bayes.csv")
