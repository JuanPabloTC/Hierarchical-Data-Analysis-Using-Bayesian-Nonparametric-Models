#@ simulations data: 3. Model: M2

# to run this script first execute the functions of the model BNP

# ==============================================================================================
#                                    PREPARACION
# ==============================================================================================

#------------------------------------------------------------------------------
# Install and load the required packages
#
#------------------------------------------------------------------------------

#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("corrplot")
#install.packages("tictoc")
#install.packages("coda")
#install.packages("posterior")
#install.packages("mclust")               # For ARI

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(dplyr))) 
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(corrplot)))
suppressMessages(suppressWarnings(library(tictoc)))
suppressMessages(suppressWarnings(library(coda)))
suppressMessages(suppressWarnings(library(posterior)))
suppressMessages(suppressWarnings(library(mclust)))

#---------------------
# Path specification
#---------------------
#Main path
ruta <- setwd('D:/Actualizado/Maestria estadistica/Tesis/Modelos BNP/20-04-26/Simulaciones/sim_data3_M2')

# Path for text files
path_text_files <- paste0(ruta,'/text_files')
# Create dir (folder) if does not exist
if (!dir.exists(path_text_files)) dir.create(path_text_files, recursive = TRUE, showWarnings = FALSE)

# Path for results
path_to_results <- paste0(ruta,'/Resultados')
# Create dir (folder) if does not exist
if (!dir.exists(path_to_results)) dir.create(path_to_results, recursive = TRUE, showWarnings = FALSE)

#=============================================

#             DATABASE CREATION

#=============================================

#-------------------------------------------------------------------------------
# This function generates the simulated database for a model with three 
#Dirichlet Process, keeping the hierarchies of individual, municipality 
#and department (3 hierarchies).Generates clusters of individuals by 
#their differentiated beta coefficients
#-------------------------------------------------------------------------------



generate_hierarchical_data_3DP <- function(
    m,                           # number of departments
    n_mun_per_dept,              # vector of length m: municipalities per department
    n_obs_range,                 # c(min, max): range for observations per municipality
    p_E,                         # number of individual-level covariates
    p_M,                         # number of municipality-level covariates
    p_D,                         # number of department-level covariates
    k_true_E,                    # true number of E-level clusters (individual)
    k_true_M,                    # true number of M-level clusters (municipal)
    k_true_D,                    # true number of D-level clusters (departmental)
    mix_probs_E,                 # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
    mix_probs_M,                 # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
    mix_probs_D,                 # vector of length k_true_D: D-level mixing probabilities (must sum to 1)
    mu_centers_E,                # vector of length k_true_E: centers for beta_k_E atoms
    mu_centers_M,                # vector of length k_true_M: centers for beta_k_M atoms
    mu_centers_D,                # vector of length k_true_D: centers for beta_k_D atoms
    sigma_k_E,                   # sd for sampling beta_k_E coefficients around centers
    sigma_k_M,                   # sd for sampling beta_k_M coefficients around centers
    sigma_k_D,                   # sd for sampling beta_k_D coefficients around centers
    kappa_min,                   # minimum municipality variance
    kappa_max                    # maximum municipality variance
) {
  
  # ===================================================================
  #                        Input validation
  # ===================================================================
  
  # Check that n_mun_per_dept has length m
  if (length(n_mun_per_dept) != m) {
    stop("Error: length(n_mun_per_dept) must equal m. ",
         "You provided m = ", m, " but n_mun_per_dept has length ",
         length(n_mun_per_dept))
  }
  
  # Check that mix_probs_E has length k_true_E and sums to 1
  if (length(mix_probs_E) != k_true_E)
    stop("Error: length(mix_probs_E) must equal k_true_E = ", k_true_E)
  if (abs(sum(mix_probs_E) - 1.0) > 1e-10)
    stop("Error: mix_probs_E must sum to 1. Current sum = ", sum(mix_probs_E))
  
  # Check that mix_probs_M has length k_true_M and sums to 1
  if (length(mix_probs_M) != k_true_M)
    stop("Error: length(mix_probs_M) must equal k_true_M = ", k_true_M)
  if (abs(sum(mix_probs_M) - 1.0) > 1e-10)
    stop("Error: mix_probs_M must sum to 1. Current sum = ", sum(mix_probs_M))
  
  # Check that mix_probs_D has length k_true_D and sums to 1
  if (length(mix_probs_D) != k_true_D)
    stop("Error: length(mix_probs_D) must equal k_true_D = ", k_true_D)
  if (abs(sum(mix_probs_D) - 1.0) > 1e-10)
    stop("Error: mix_probs_D must sum to 1. Current sum = ", sum(mix_probs_D))
  
  # Check that mu_centers have correct lengths
  if (length(mu_centers_E) != k_true_E)
    stop("Error: length(mu_centers_E) must equal k_true_E = ", k_true_E)
  if (length(mu_centers_M) != k_true_M)
    stop("Error: length(mu_centers_M) must equal k_true_M = ", k_true_M)
  if (length(mu_centers_D) != k_true_D)
    stop("Error: length(mu_centers_D) must equal k_true_D = ", k_true_D)
  
  # Check that n_obs_range has length 2 and minimum >= 1
  if (length(n_obs_range) != 2)
    stop("Error: n_obs_range must be a vector of length 2: c(min, max)")
  if (n_obs_range[1] < 1)
    stop("Error: minimum observations per municipality must be >= 1")
  
  set.seed(777)
  
  # ===================================================================
  #                    Derived quantities
  # ===================================================================
  
  n_mun_total <- sum(n_mun_per_dept)   # total municipalities
  
  # Sample observations per municipality from discrete uniform
  n_obs_per_mun_vec <- sample(n_obs_range[1]:n_obs_range[2],
                              size    = n_mun_total,
                              replace = TRUE)
  
  n <- sum(n_obs_per_mun_vec)          # total observations
  
  cat("Generating 3-DP dataset with:\n")
  cat("  Departments (m)         :", m,          "\n")
  cat("  Total municipalities    :", n_mun_total, "\n")
  cat("  Total observations (n)  :", n,           "\n")
  cat("  Covariates: p_E =", p_E, ", p_M =", p_M, ", p_D =", p_D, "\n")
  cat("  True E-level clusters (k_true_E):", k_true_E, "\n")
  cat("  True M-level clusters (k_true_M):", k_true_M, "\n")
  cat("  True D-level clusters (k_true_D):", k_true_D, "\n\n")
  
  # ===================================================================
  #                   Build hierarchical structure
  # ===================================================================
  
  dept_id    <- integer(n)
  mun_id     <- integer(n)
  mun_map    <- vector("list", n_mun_total)
  dept_to_mun <- vector("list", m)
  obs_to_mun <- integer(n)
  
  # Build mapping (order by department -> municipality -> observations)
  obs_counter <- 1L
  mun_counter <- 1L
  
  for (q in seq_len(m)) {
    dept_to_mun[[q]] <- integer(0)
    n_mun_q <- n_mun_per_dept[q]
    
    for (j in seq_len(n_mun_q)) {
      n_jq    <- n_obs_per_mun_vec[mun_counter]   # observations for this municipality
      start_i <- obs_counter
      end_i   <- obs_counter + n_jq - 1L
      
      mun_map[[mun_counter]] <- list(
        q     = q,
        j     = j,
        start = start_i,
        end   = end_i,
        n     = n_jq
      )
      dept_to_mun[[q]] <- c(dept_to_mun[[q]], mun_counter)
      
      # Fill IDs
      dept_id[start_i:end_i]    <- q
      mun_id[start_i:end_i]     <- mun_counter
      obs_to_mun[start_i:end_i] <- mun_counter
      
      obs_counter <- obs_counter + n_jq
      mun_counter <- mun_counter + 1L
    }
  }
  
  # ===================================================================
  #                      Generate covariates
  # ===================================================================
  
  # 1) Individual-level covariates: X_E_full (n x p_E)
  #    Varies per observation — drawn independently from N(0,1)
  if (p_E > 0) {
    X_E_full <- matrix(rnorm(n * p_E, mean = 0, sd = 1), nrow = n, ncol = p_E)
  } else {
    X_E_full <- matrix(nrow = n, ncol = 0)
  }
  colnames(X_E_full) <- paste0("Individual_", seq_len(p_E))
  X_E_names          <- colnames(X_E_full)
  
  # 2) Municipality-level covariates: Z_full (n x p_M)
  #    Constant per municipality — one value per municipality, expanded to observations
  if (p_M > 0) {
    mun_cov_vals <- matrix(rnorm(n_mun_total * p_M, mean = 0, sd = 1),
                           nrow = n_mun_total, ncol = p_M)
    Z_full <- mun_cov_vals[mun_id, , drop = FALSE]  # expand to observations
  } else {
    Z_full <- matrix(nrow = n, ncol = 0)
  }
  colnames(Z_full) <- paste0("Municipal_", seq_len(p_M))
  Z_M_names        <- colnames(Z_full)   
  
  # 3) Department-level covariates: W_full (n x p_D)
  #    Constant per department — one value per department, expanded to observations
  if (p_D > 0) {
    dept_cov_vals <- matrix(rnorm(m * p_D, mean = 0, sd = 1),
                            nrow = m, ncol = p_D)
    W_full <- dept_cov_vals[dept_id, , drop = FALSE]  # expand to observations
  } else {
    W_full <- matrix(nrow = n, ncol = 0)
  }
  colnames(W_full) <- paste0("Departmental_", seq_len(p_D))
  W_D_names        <- colnames(W_full)
  
  # Combined flat covariate matrix for one-DP model compatibility (n x (p_E + p_M + p_D))
  X_full  <- cbind(X_E_full, Z_full, W_full)
  x_names <- colnames(X_full)
  
  # Precompute squared norms (used by samplers)
  xTx_vec <- rowSums(X_E_full^2)    # x_i^T x_i per observation (E-level)
  zTz_vec <- rowSums(Z_full^2)      # z_i^T z_i per observation (M-level)
  wTw_vec <- rowSums(W_full^2)      # w_i^T w_i per observation (D-level)
  
  # With the full matrix
  xTx_full_vec <- rowSums(X_full^2) 
  
  # ===================================================================
  #                   Generate true parameters
  # ===================================================================
  
  # True intercept from N(0, 1)
  beta_int_true <- rnorm(1, mean = 0, sd = 1)
  
  # True E-level atoms: each cluster k has p_E coefficients centered at mu_centers_E[k]
  beta_k_E_true <- vector("list", k_true_E)
  for (k in seq_len(k_true_E)) {
    beta_k_E_true[[k]] <- rnorm(p_E, mean = mu_centers_E[k], sd = sigma_k_E)
  }
  
  # True M-level atoms: each cluster k has p_M coefficients centered at mu_centers_M[k]
  beta_k_M_true <- vector("list", k_true_M)
  for (k in seq_len(k_true_M)) {
    beta_k_M_true[[k]] <- rnorm(p_M, mean = mu_centers_M[k], sd = sigma_k_M)
  }
  
  # True D-level atoms: each cluster k has p_D coefficients centered at mu_centers_D[k]
  beta_k_D_true <- vector("list", k_true_D)
  for (k in seq_len(k_true_D)) {
    beta_k_D_true[[k]] <- rnorm(p_D, mean = mu_centers_D[k], sd = sigma_k_D)
  }
  
  # True municipality variances: uniform on [kappa_min, kappa_max]
  mun_kappa  <- runif(n_mun_total, min = kappa_min, max = kappa_max)
  sqrt_kappa <- sqrt(mun_kappa)
  
  # ===================================================================
  #                   Assign clusters and generate Y
  # ===================================================================
  
  # Assign E-level cluster for each observation based on mix_probs_E
  xi_E_true <- sample.int(k_true_E, size = n, replace = TRUE, prob = mix_probs_E)
  
  # Assign M-level cluster for each observation based on mix_probs_M
  xi_M_true <- sample.int(k_true_M, size = n, replace = TRUE, prob = mix_probs_M)
  
  # Assign D-level cluster for each observation based on mix_probs_D
  xi_D_true <- sample.int(k_true_D, size = n, replace = TRUE, prob = mix_probs_D)
  
  # Build beta per observation matrices (K_true x p)
  beta_E_mat     <- do.call(rbind, beta_k_E_true)        # k_true_E x p_E
  beta_M_mat     <- do.call(rbind, beta_k_M_true)        # k_true_M x p_M
  beta_D_mat     <- do.call(rbind, beta_k_D_true)        # k_true_D x p_D
  
  # Rearange beta_E_mat according to cluster assignments
  beta_E_per_obs <- beta_E_mat[xi_E_true, , drop = FALSE]  # n x p_E
  beta_M_per_obs <- beta_M_mat[xi_M_true, , drop = FALSE]  # n x p_M
  beta_D_per_obs <- beta_D_mat[xi_D_true, , drop = FALSE]  # n x p_D
  
  # Full linear predictor:
  # vartheta_i = beta_int + x_i^T beta_k_E[xi_E[i]] + z_i^T beta_k_M[xi_M[i]] + w_i^T beta_k_D[xi_D[i]]
  lin_pred <- beta_int_true +                                    #length = n
    rowSums(X_E_full * beta_E_per_obs) +   # x_i^T beta_k_E
    rowSums(Z_full   * beta_M_per_obs) +   # z_i^T beta_k_M
    rowSums(W_full   * beta_D_per_obs)     # w_i^T beta_k_D
  
  # Add observation-level noise with variance = mun_kappa[obs_to_mun[i]]
  y <- lin_pred + rnorm(n, mean = 0, sd = sqrt(mun_kappa[obs_to_mun]))
  
  # ===================================================================
  #                   Create nested Y structure
  # ===================================================================
  
  Y <- vector("list", m)
  for (q in seq_len(m)) {
    n_mun_q <- length(dept_to_mun[[q]])
    Y[[q]] <- vector("list", length = n_mun_q)
    mun_indices_q <- dept_to_mun[[q]]
    
    for (idx_local in seq_along(mun_indices_q)) {
      mun_idx <- mun_indices_q[idx_local]
      s <- mun_map[[mun_idx]]$start
      e <- mun_map[[mun_idx]]$end
      Y[[q]][[idx_local]] <- y[s:e]
    }
  }
  
  # ===================================================================
  #                   Create data frame
  # ===================================================================
  
  datos <- data.frame(
    obs = seq_len(n), 
    dept = dept_id, 
    mun = mun_id, 
    y = y)
  
  datos <- cbind(datos, X_full)
  
  # ===================================================================
  #                   Print summary
  # ===================================================================
  
  cat("Summary:\n")
  cat("  Municipality sizes (min, median, max):",
      min(n_obs_per_mun_vec), median(n_obs_per_mun_vec), max(n_obs_per_mun_vec), "\n")
  cat("  True intercept (beta_int_true):", round(beta_int_true, 3), "\n")
  cat("  mun_kappa range:", round(range(mun_kappa), 2), "\n")
  cat("  E-level cluster sizes:\n"); print(table(xi_E_true))
  cat("  M-level cluster sizes:\n"); print(table(xi_M_true))
  cat("  D-level cluster sizes:\n"); print(table(xi_D_true))
  cat("\n")
  
  # ===================================================================
  #                  Initialization objects — 3-DP model
  # ===================================================================
  
  # --- xi_E_false: random cluster assignment for E-level (initialization) ---
  xi_E_false <- sample.int(k_true_E, size = n, replace = TRUE, prob = mix_probs_E)
  # In case some cluster is empty, it deletes it from the betas and from xi
  used_E     <- sort(unique(xi_E_false))
  xi_E_false <- as.integer(factor(xi_E_false, levels = used_E, labels = seq_along(used_E)))
  
  # --- beta_k_E_false: initial E-level atoms (drop atoms for empty clusters) ---
  beta_k_E_false <- vector("list", k_true_E)
  for (k in seq_len(k_true_E)){
    beta_k_E_false[[k]] <- rnorm(p_E, mean = mu_centers_E[k], sd = sigma_k_E)
  }
  # in case of empty clusters, deletes its beta
  beta_k_E_false <- beta_k_E_false[used_E]
  
  
  # --- xi_M_false: random cluster assignment for M-level (initialization) ---
  xi_M_false <- sample.int(k_true_M, size = n, replace = TRUE, prob = mix_probs_M)
  used_M     <- sort(unique(xi_M_false))
  xi_M_false <- as.integer(factor(xi_M_false, levels = used_M, labels = seq_along(used_M)))
  
  # --- beta_k_M_false: initial M-level atoms (drop atoms for empty clusters) ---
  beta_k_M_false <- vector("list", k_true_M)
  for (k in seq_len(k_true_M)){
    beta_k_M_false[[k]] <- rnorm(p_M, mean = mu_centers_M[k], sd = sigma_k_M)
  }
  beta_k_M_false <- beta_k_M_false[used_M]
  
  
  # --- xi_D_false: random cluster assignment for D-level (initialization) ---
  xi_D_false <- sample.int(k_true_D, size = n, replace = TRUE, prob = mix_probs_D)
  used_D     <- sort(unique(xi_D_false))
  xi_D_false <- as.integer(factor(xi_D_false, levels = used_D, labels = seq_along(used_D)))
  
  # --- beta_k_D_false: initial D-level atoms (drop atoms for empty clusters) ---
  beta_k_D_false <- vector("list", k_true_D)
  for (k in seq_len(k_true_D)){
    beta_k_D_false[[k]] <- rnorm(p_D, mean = mu_centers_D[k], sd = sigma_k_D)
  }
  beta_k_D_false <- beta_k_D_false[used_D]
  
  
  # ===================================================================
  #         Initialization objects — one-DP model (model 2 compatibility)
  #         Uses combined X_full (p_E + p_M + p_D columns) and a single xi
  # ===================================================================
  
  p       <- p_E + p_M + p_D       # total covariates for one-DP model
  k_true  <- k_true_E              # use E-level cluster count as reference for one-DP
  
  # --- xi_false: single cluster assignment for one-DP model ---
  # Uses E level xi 
  xi_false <- sample.int(k_true, size = n, replace = TRUE, prob = mix_probs_E)
  used_1dp <- sort(unique(xi_false))
  xi_false <- as.integer(factor(xi_false, levels = used_1dp, labels = seq_along(used_1dp)))
  
  # --- beta_k_false: initial atoms for one-DP model (length p = p_E + p_M + p_D) ---
  # uses beta_K_E
  beta_k_false <- vector("list", k_true)
  for (k in seq_len(k_true)){
    beta_k_false[[k]] <- rnorm(p, mean = mu_centers_E[k], sd = sigma_k_E)
  }
  beta_k_false <- beta_k_false[used_1dp]
  
  # --- beta_int_false: initial intercept for one-DP model ---
  beta_int_false <- rnorm(1, mean = 0, sd = 1)
  
  # ===============================================
  #              kappa2_q initialization 
  # ===============================================
  
  kappa2_q <- sapply(seq_len(m), function(q) {
    mun_indices_q <- dept_to_mun[[q]]
    mean(mun_kappa[mun_indices_q])   # mean of true kappa2_{j,q} within dept q
  })
  
  # ===================================================================
  #                     Return objects
  # ===================================================================
  
  return(list(
    
    # --- Response ---
    y = y,
    
    # --- Three-level covariate matrices ---
    X_E_full = X_E_full,      # (n x p_E) individual-level
    Z_full   = Z_full,        # (n x p_M) municipality-level
    W_full   = W_full,        # (n x p_D) department-level
    
    X_E_names  = X_E_names,
    Z_M_names  = Z_M_names,
    W_D_names  = W_D_names,
    
    X_full   = X_full,        # (n x p)   combined (for one-DP model)
    x_names  = x_names,
    
    # --- Dimensions ---
    n           = n,
    m           = m,
    p_E         = p_E,
    p_M         = p_M,
    p_D         = p_D,
    p           = p_E + p_M + p_D,
    k_true_E    = k_true_E,
    k_true_M    = k_true_M,
    k_true_D    = k_true_D,
    n_mun_total = n_mun_total,
    
    k_true      = k_true,
    
    # --- Hierarchical structure ---
    dept_id     = dept_id,
    mun_id      = mun_id,
    obs_to_mun  = obs_to_mun,
    mun_map     = mun_map,
    dept_to_mun = dept_to_mun,
    
    # --- Precomputed values (used by 3-DP samplers) ---
    xTx_vec      = xTx_vec,
    zTz_vec      = zTz_vec,
    wTw_vec      = wTw_vec,
    
    xTx_full_vec = xTx_full_vec,
    
    sqrt_kappa = sqrt_kappa,
    
    # --- True parameters ---
    beta_int_true  = beta_int_true,
    beta_k_E_true  = beta_k_E_true,
    beta_k_M_true  = beta_k_M_true,
    beta_k_D_true  = beta_k_D_true,
    mun_kappa      = mun_kappa,
    xi_E_true      = xi_E_true,
    xi_M_true      = xi_M_true,
    xi_D_true      = xi_D_true,
    mix_probs_E    = mix_probs_E,
    mix_probs_M    = mix_probs_M,
    mix_probs_D    = mix_probs_D,
    
    # --- Nested Y structure ---
    Y = Y,
    
    # --- Initialization objects: 3-DP model ---
    xi_E_false     = xi_E_false,
    xi_M_false     = xi_M_false,
    xi_D_false     = xi_D_false,
    beta_k_E_false = beta_k_E_false,
    beta_k_M_false = beta_k_M_false,
    beta_k_D_false = beta_k_D_false,
    
    # --- Initialization objects: one-DP model (model 2 compatibility) ---
    xi_false       = xi_false,
    beta_k_false   = beta_k_false,
    
    beta_int_false = beta_int_false,
    
    # --- Shared initialization ---
    kappa2_q = kappa2_q,
    
    # --- Data frame ---
    datos = datos
  ))
}


# ==============================================================================================
#                                  SIMULATED SCENARIOS
# ==============================================================================================

#-----------------------------------------------------------------------------
# Generate 15 simulated scenarios using the function 
# generate_hierarchical_data_3DP
#-----------------------------------------------------------------------------

# =============================
# 1. VARIANCE SCENARIOS
# =============================

# In this variance scenarios we keep the number of departments moderated, the number 
#of municipalities per department and the number of observations per municipality homogeneous.
#Number of covariates is also moderated 12. Mixing probabilities homogeneous in the 4 clusters

# The only change along scenarios is the variance. We process three scenarios: low, moderated and
#high variance:

# The variance is relative to the cluster centers

# Low: ~ 4
# Moderated: ~ 10 
# High: ~ 20

# Constant parameters
m_1              = 32                          # number of departments
n_mun_per_dept_1 = rep(5, m_1)                 # vector of length m: municipalities per department
n_obs_range_1    = c(8, 11)                   # c(min, max): range for observations per municipality 
p_E_1            = 3                           # number of individual-level covariates
p_M_1            = 3                           # number of municipality-level covariates
p_D_1            = 3                           # number of department-level covariates
k_true_E_1       = 2                           # true number of E-level clusters (individual)
k_true_M_1       = 2                           # true number of M-level clusters (municipal)
k_true_D_1       = 2                           # true number of D-level clusters (departmental)
mix_probs_E_1    = rep(1/2, 2)                 # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_1    = rep(1/2, 2)                 # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_1    = rep(1/2, 2)                 # vector of length k_true_D: D-level mixing probabilities (must sum to 1)
mu_centers_E_1   = c(1, 8)                     # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_1   = c(1, 8)                     # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_1   = c(1, 8)                     # vector of length k_true_D: centers for beta_k_D atoms
sigma_k_E_1      = 1.0                         # sd for sampling beta_k_E coefficients around centers
sigma_k_M_1      = 1.0                         # sd for sampling beta_k_M coefficients around centers
sigma_k_D_1      = 1.0                         # sd for sampling beta_k_D coefficients around centers


# ----------------------------------
# 1. Low municipal variance scenario
# ----------------------------------
kappa_min_1_1     = 2                           # min value of municipal variance (runif)              
kappa_max_1_1     = 4                           # max value of municipal variance (runif)


simulated_data_1 <- generate_hierarchical_data_3DP(
  m                   = m_1,                          
  n_mun_per_dept      = n_mun_per_dept_1,             
  n_obs_range         = n_obs_range_1,                
  p_E                 = p_E_1,                        
  p_M                 = p_M_1,                        
  p_D                 = p_D_1,                       
  k_true_E            = k_true_E_1,                 
  k_true_M            = k_true_M_1,                   
  k_true_D            = k_true_D_1,                   
  mix_probs_E         = mix_probs_E_1,                 
  mix_probs_M         = mix_probs_M_1,                 
  mix_probs_D         = mix_probs_D_1,               
  mu_centers_E        = mu_centers_E_1,               
  mu_centers_M        = mu_centers_M_1,                
  mu_centers_D        = mu_centers_D_1,             
  sigma_k_E           = sigma_k_E_1,             
  sigma_k_M           = sigma_k_M_1,             
  sigma_k_D           = sigma_k_D_1,             
  kappa_min           = kappa_min_1_1,                    
  kappa_max           = kappa_max_1_1                    
)


# ----------------------------------
# 2. Moderated municipal variance scenario
# ----------------------------------
kappa_min_1_2     = 8                            # min value of municipal variance (runif)              
kappa_max_1_2     = 10                           # max value of municipal variance (runif)

simulated_data_2 <- generate_hierarchical_data_3DP(
  m                   = m_1,                          
  n_mun_per_dept      = n_mun_per_dept_1,             
  n_obs_range         = n_obs_range_1,                
  p_E                 = p_E_1,                        
  p_M                 = p_M_1,                        
  p_D                 = p_D_1,                       
  k_true_E            = k_true_E_1,                 
  k_true_M            = k_true_M_1,                   
  k_true_D            = k_true_D_1,                   
  mix_probs_E         = mix_probs_E_1,                 
  mix_probs_M         = mix_probs_M_1,                 
  mix_probs_D         = mix_probs_D_1,               
  mu_centers_E        = mu_centers_E_1,               
  mu_centers_M        = mu_centers_M_1,                
  mu_centers_D        = mu_centers_D_1,             
  sigma_k_E           = sigma_k_E_1,             
  sigma_k_M           = sigma_k_M_1,             
  sigma_k_D           = sigma_k_D_1,             
  kappa_min           = kappa_min_1_2,                    
  kappa_max           = kappa_max_1_2                    
)


# ----------------------------------
# 3. High municipal variance scenario
# ----------------------------------
kappa_min_1_3     = 12                            # min value of municipal variance (runif)              
kappa_max_1_3     = 14                            # max value of municipal variance (runif)

simulated_data_3 <- generate_hierarchical_data_3DP(
  m                   = m_1,                          
  n_mun_per_dept      = n_mun_per_dept_1,             
  n_obs_range         = n_obs_range_1,                
  p_E                 = p_E_1,                        
  p_M                 = p_M_1,                        
  p_D                 = p_D_1,                       
  k_true_E            = k_true_E_1,                 
  k_true_M            = k_true_M_1,                   
  k_true_D            = k_true_D_1,                   
  mix_probs_E         = mix_probs_E_1,                 
  mix_probs_M         = mix_probs_M_1,                 
  mix_probs_D         = mix_probs_D_1,               
  mu_centers_E        = mu_centers_E_1,               
  mu_centers_M        = mu_centers_M_1,                
  mu_centers_D        = mu_centers_D_1,             
  sigma_k_E           = sigma_k_E_1,             
  sigma_k_M           = sigma_k_M_1,             
  sigma_k_D           = sigma_k_D_1,             
  kappa_min           = kappa_min_1_3,                    
  kappa_max           = kappa_max_1_3                    
)



# =============================
# 2. NUMBER OF CLUSTER K SCENARIOS
# =============================

# In this three scenarios we maintain number of departments moderated, homogeneous
#number of municipalities per department and observations per municipality. Also, 
#moderated number of covariates (15), and homogeneous mixing probabilities, 
#and low variance.

# The only change is the number of clusters:
# - Low clusters K = 3
# - Moderated    K = 10
# - High         K = 20


# Constant parameters
m_2              = 32                          # number of departments
n_mun_per_dept_2 = rep(4, m_2)                 # vector of length m: municipalities per department
n_obs_range_2    = c(17, 18)                   # c(min, max): range for observations per municipality
p_E_2            = 3                           # number of individual-level covariates
p_M_2            = 3                           # number of municipality-level covariates
p_D_2            = 3                           # number of department-level covariates
sigma_k_E_2      = 1.0                         # sd for sampling beta_k_E coefficients around centers
sigma_k_M_2      = 1.0                         # sd for sampling beta_k_M coefficients around centers
sigma_k_D_2      = 1.0                         # sd for sampling beta_k_D coefficients around centers
kappa_min_2      = 2                           # min value of municipal variance (runif)              
kappa_max_2      = 4                           # max value of municipal variance (runif)


# ----------------------------------
# 4. Low number of clusters k scenario
# ----------------------------------

k_true_E_2_1     = 2                           # true number of E-level clusters (individual)
k_true_M_2_1     = 2                           # true number of M-level clusters (municipal)
k_true_D_2_1     = 2                           # true number of D-level clusters (departmental)

mu_centers_E_2_1 = c(1, 8)                   # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_2_1 = c(1, 8)                   # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_2_1 = c(1, 8)                   # vector of length k_true_D: centers for beta_k_D atoms

mix_probs_E_2_1  = rep(1/2, 2)                 # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_2_1  = rep(1/2, 2)                 # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_2_1  = rep(1/2, 2)                 # vector of length k_true_D: D-level mixing probabilities (must sum to 1)



simulated_data_4 <- generate_hierarchical_data_3DP(
  m                   = m_2,                          
  n_mun_per_dept      = n_mun_per_dept_2,             
  n_obs_range         = n_obs_range_2,                
  p_E                 = p_E_2,                        
  p_M                 = p_M_2,                        
  p_D                 = p_D_2,                       
  k_true_E            = k_true_E_2_1,                 
  k_true_M            = k_true_M_2_1,                   
  k_true_D            = k_true_D_2_1,                   
  mix_probs_E         = mix_probs_E_2_1,                 
  mix_probs_M         = mix_probs_M_2_1,                 
  mix_probs_D         = mix_probs_D_2_1,               
  mu_centers_E        = mu_centers_E_2_1,               
  mu_centers_M        = mu_centers_M_2_1,                
  mu_centers_D        = mu_centers_D_2_1,             
  sigma_k_E           = sigma_k_E_2,             
  sigma_k_M           = sigma_k_M_2,             
  sigma_k_D           = sigma_k_D_2,             
  kappa_min           = kappa_min_2,                    
  kappa_max           = kappa_max_2                    
)


# ----------------------------------
# 5. Moderated number of clusters k scenario
# ----------------------------------

k_true_E_2_2     = 4                          # true number of E-level clusters (individual)
k_true_M_2_2     = 4                          # true number of M-level clusters (municipal)
k_true_D_2_2     = 4                          # true number of D-level clusters (departmental)

mu_centers_E_2_2 = seq(from = 1, by = 6, 
                       length.out = k_true_E_2_2)   # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_2_2 = seq(from = 1, by = 6, 
                       length.out = k_true_M_2_2)   # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_2_2 = seq(from = 1, by = 6, 
                       length.out = k_true_D_2_2)   # vector of length k_true_D: centers for beta_k_D atoms

mix_probs_E_2_2  = rep(0.25, k_true_E_2_2)      # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_2_2  = rep(0.25, k_true_M_2_2)      # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_2_2  = rep(0.25, k_true_D_2_2)      # vector of length k_true_D: D-level mixing probabilities (must sum to 1)


simulated_data_5 <- generate_hierarchical_data_3DP(
  m                   = m_2,                          
  n_mun_per_dept      = n_mun_per_dept_2,             
  n_obs_range         = n_obs_range_2,                
  p_E                 = p_E_2,                        
  p_M                 = p_M_2,                        
  p_D                 = p_D_2,                       
  k_true_E            = k_true_E_2_2,                 
  k_true_M            = k_true_M_2_2,                   
  k_true_D            = k_true_D_2_2,                   
  mix_probs_E         = mix_probs_E_2_2,                 
  mix_probs_M         = mix_probs_M_2_2,                 
  mix_probs_D         = mix_probs_D_2_2,               
  mu_centers_E        = mu_centers_E_2_2,               
  mu_centers_M        = mu_centers_M_2_2,                
  mu_centers_D        = mu_centers_D_2_2,             
  sigma_k_E           = sigma_k_E_2,             
  sigma_k_M           = sigma_k_M_2,             
  sigma_k_D           = sigma_k_D_2,             
  kappa_min           = kappa_min_2,                    
  kappa_max           = kappa_max_2                    
)

# ----------------------------------
# 6. High number of clusters k scenario
# ----------------------------------

k_true_E_2_3     = 6                          # true number of E-level clusters (individual)
k_true_M_2_3     = 6                          # true number of M-level clusters (municipal)
k_true_D_2_3     = 6                          # true number of D-level clusters (departmental)

mu_centers_E_2_3 = seq(from = 1, by = 5, 
                       length.out = k_true_E_2_3)   # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_2_3 = seq(from = 1, by = 5, 
                       length.out = k_true_M_2_3)   # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_2_3 = seq(from = 1, by = 5, 
                       length.out = k_true_D_2_3)   # vector of length k_true_D: centers for beta_k_D atoms

mix_probs_E_2_3  = rep(1/k_true_E_2_3, k_true_E_2_3)  # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_2_3  = rep(1/k_true_M_2_3, k_true_M_2_3)  # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_2_3  = rep(1/k_true_D_2_3, k_true_D_2_3)  # vector of length k_true_D: D-level mixing probabilities (must sum to 1)


simulated_data_6 <- generate_hierarchical_data_3DP(
  m                   = m_2,                          
  n_mun_per_dept      = n_mun_per_dept_2,             
  n_obs_range         = n_obs_range_2,                
  p_E                 = p_E_2,                        
  p_M                 = p_M_2,                        
  p_D                 = p_D_2,                       
  k_true_E            = k_true_E_2_3,                 
  k_true_M            = k_true_M_2_3,                   
  k_true_D            = k_true_D_2_3,                   
  mix_probs_E         = mix_probs_E_2_3,                 
  mix_probs_M         = mix_probs_M_2_3,                 
  mix_probs_D         = mix_probs_D_2_3,               
  mu_centers_E        = mu_centers_E_2_3,               
  mu_centers_M        = mu_centers_M_2_3,                
  mu_centers_D        = mu_centers_D_2_3,             
  sigma_k_E           = sigma_k_E_2,             
  sigma_k_M           = sigma_k_M_2,             
  sigma_k_D           = sigma_k_D_2,             
  kappa_min           = kappa_min_2,                    
  kappa_max           = kappa_max_2                    
)



# =============================
# 3. NUMBER OF COVARIATES p
# =============================

# In this three scenarios we maintain number of departments moderated, 
#moderated number of clusters, and homogeneous number of municipalities 
#per department and observations per municipality. 
#Also homogeneous mixing probabilities, and low variance.


# The only change is the number of covariates:
# - Low          p_l = 3
# - Moderated    p_l = 10
# - High         p_l = 25


# Constant parameters
m_3              = 32                          # number of departments
n_mun_per_dept_3 = rep(5, m_3)                 # vector of length m: municipalities per department
n_obs_range_3    = c(16, 18)                   # c(min, max): range for observations per municipality
k_true_E_3       = 2                           # true number of E-level clusters (individual)
k_true_M_3       = 2                           # true number of M-level clusters (municipal)
k_true_D_3       = 2                           # true number of D-level clusters (departmental)
mix_probs_E_3    = rep(1/2, 2)                 # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_3    = rep(1/2, 2)                 # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_3    = rep(1/2, 2)                 # vector of length k_true_D: D-level mixing probabilities (must sum to 1)
mu_centers_E_3   = c(1, 8)                   # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_3   = c(1, 8)                   # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_3   = c(1, 8)                   # vector of length k_true_D: centers for beta_k_D atoms
sigma_k_E_3      = 1.0                         # sd for sampling beta_k_E coefficients around centers
sigma_k_M_3      = 1.0                         # sd for sampling beta_k_M coefficients around centers
sigma_k_D_3      = 1.0                         # sd for sampling beta_k_D coefficients around centers
kappa_min_3      = 2                           # min value of municipal variance (runif)              
kappa_max_3      = 4                           # max value of municipal variance (runif)


# ----------------------------------
# 7. Low number of covariates p scenario
# ----------------------------------

p_E_3_1          = 2                           # number of individual-level covariates
p_M_3_1          = 2                           # number of municipality-level covariates
p_D_3_1          = 2                           # number of department-level covariates


simulated_data_7 <- generate_hierarchical_data_3DP(
  m                   = m_3,                          
  n_mun_per_dept      = n_mun_per_dept_3,             
  n_obs_range         = n_obs_range_3,                
  p_E                 = p_E_3_1,                        
  p_M                 = p_M_3_1,                        
  p_D                 = p_D_3_1,                       
  k_true_E            = k_true_E_3,                 
  k_true_M            = k_true_M_3,                   
  k_true_D            = k_true_D_3,                   
  mix_probs_E         = mix_probs_E_3,                 
  mix_probs_M         = mix_probs_M_3,                 
  mix_probs_D         = mix_probs_D_3,               
  mu_centers_E        = mu_centers_E_3,               
  mu_centers_M        = mu_centers_M_3,                
  mu_centers_D        = mu_centers_D_3,             
  sigma_k_E           = sigma_k_E_3,             
  sigma_k_M           = sigma_k_M_3,             
  sigma_k_D           = sigma_k_D_3,             
  kappa_min           = kappa_min_3,                    
  kappa_max           = kappa_max_3                    
)

# ----------------------------------
# 8. Moderated number of covariates p scenario
# ----------------------------------

p_E_3_2          = 5                           # number of individual-level covariates
p_M_3_2          = 5                           # number of municipality-level covariates
p_D_3_2          = 5                           # number of department-level covariates

simulated_data_8 <- generate_hierarchical_data_3DP(
  m                   = m_3,                          
  n_mun_per_dept      = n_mun_per_dept_3,             
  n_obs_range         = n_obs_range_3,                
  p_E                 = p_E_3_2,                        
  p_M                 = p_M_3_2,                        
  p_D                 = p_D_3_2,                       
  k_true_E            = k_true_E_3,                 
  k_true_M            = k_true_M_3,                   
  k_true_D            = k_true_D_3,                   
  mix_probs_E         = mix_probs_E_3,                 
  mix_probs_M         = mix_probs_M_3,                 
  mix_probs_D         = mix_probs_D_3,               
  mu_centers_E        = mu_centers_E_3,               
  mu_centers_M        = mu_centers_M_3,                
  mu_centers_D        = mu_centers_D_3,             
  sigma_k_E           = sigma_k_E_3,             
  sigma_k_M           = sigma_k_M_3,             
  sigma_k_D           = sigma_k_D_3,             
  kappa_min           = kappa_min_3,                    
  kappa_max           = kappa_max_3                    
)


# ----------------------------------
# 9. High number of covariates p scenario
# ----------------------------------

p_E_3_3          = 7                           # number of individual-level covariates
p_M_3_3          = 7                           # number of municipality-level covariates
p_D_3_3          = 5                            # number of department-level covariates

simulated_data_9 <- generate_hierarchical_data_3DP(
  m                   = m_3,                          
  n_mun_per_dept      = n_mun_per_dept_3,             
  n_obs_range         = n_obs_range_3,                
  p_E                 = p_E_3_3,                        
  p_M                 = p_M_3_3,                        
  p_D                 = p_D_3_3,                       
  k_true_E            = k_true_E_3,                 
  k_true_M            = k_true_M_3,                   
  k_true_D            = k_true_D_3,                   
  mix_probs_E         = mix_probs_E_3,                 
  mix_probs_M         = mix_probs_M_3,                 
  mix_probs_D         = mix_probs_D_3,               
  mu_centers_E        = mu_centers_E_3,               
  mu_centers_M        = mu_centers_M_3,                
  mu_centers_D        = mu_centers_D_3,             
  sigma_k_E           = sigma_k_E_3,             
  sigma_k_M           = sigma_k_M_3,             
  sigma_k_D           = sigma_k_D_3,             
  kappa_min           = kappa_min_3,                    
  kappa_max           = kappa_max_3                    
)


# =============================
# 4. MIXING PROBABILITIES
# =============================

# In this scenario we keep a moderated number departments (5) and of covariates
#we keep an homogeneous number of municipalities per department, and number observations per municipality
#. Also, keep a moderated number of clusters (6). 

# The only change is the mixing probabilities:
# - 1 cluster with high probability (0.9, 0.05, 0.05)
# - 1 cluster with low probability (0.1, 0.45, 0.45)
# - 2 clusters with equal probability, 2 with low (0.33, 0.33, 0.33)


# Constant parameters
m_4              = 32                          # number of departments
n_mun_per_dept_4 = rep(4, m_4)                 # vector of length m: municipalities per department
n_obs_range_4    = c(14, 15)                   # c(min, max): range for observations per municipality
p_E_4            = 4                          # number of individual-level covariates
p_M_4            = 4                          # number of municipality-level covariates
p_D_4            = 4                          # number of department-level covariates
sigma_k_E_4      = 1.0                        # sd for sampling beta_k_E coefficients around centers
sigma_k_M_4      = 1.0                        # sd for sampling beta_k_M coefficients around centers
sigma_k_D_4      = 1.0                        # sd for sampling beta_k_D coefficients around centers
kappa_min_4      = 2                          # min value of municipal variance (runif)              
kappa_max_4      = 4                          # max value of municipal variance (runif)


# -------------------------------------------------
# 10. One cluster with high probability
# -------------------------------------------------
k_true_E_4_1       = 3                                      # true number of E-level clusters (individual)
k_true_M_4_1       = 3                                      # true number of M-level clusters (municipal)
k_true_D_4_1       = 3                                      # true number of D-level clusters (departmental)
mix_probs_E_4_1    = c(0.8, 0.1, 0.1)                     # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_4_1    = c(0.8, 0.1, 0.1)                     # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_4_1    = c(0.8, 0.1, 0.1)                     # vector of length k_true_D: D-level mixing probabilities (must sum to 1)
mu_centers_E_4_1   = seq(from = 1, by = 7, length.out = k_true_E_4_1)  # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_4_1   = seq(from = 1, by = 7, length.out = k_true_M_4_1)                # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_4_1   = seq(from = 1, by = 7, length.out = k_true_D_4_1)                # vector of length k_true_D: centers for beta_k_D atoms


simulated_data_10 <- generate_hierarchical_data_3DP(
  m                   = m_4,                          
  n_mun_per_dept      = n_mun_per_dept_4,             
  n_obs_range         = n_obs_range_4,                
  p_E                 = p_E_4,                        
  p_M                 = p_M_4,                        
  p_D                 = p_D_4,                       
  k_true_E            = k_true_E_4_1,                 
  k_true_M            = k_true_M_4_1,                   
  k_true_D            = k_true_D_4_1,                   
  mix_probs_E         = mix_probs_E_4_1,                 
  mix_probs_M         = mix_probs_M_4_1,                 
  mix_probs_D         = mix_probs_D_4_1,               
  mu_centers_E        = mu_centers_E_4_1,               
  mu_centers_M        = mu_centers_M_4_1,                
  mu_centers_D        = mu_centers_D_4_1,             
  sigma_k_E           = sigma_k_E_4,             
  sigma_k_M           = sigma_k_M_4,             
  sigma_k_D           = sigma_k_D_4,             
  kappa_min           = kappa_min_4,                    
  kappa_max           = kappa_max_4                    
)


# -------------------------------------------------
# 11. One cluster with low probability
# -------------------------------------------------
k_true_E_4_2       = 3                                            # true number of E-level clusters (individual)
k_true_M_4_2       = 3                                            # true number of M-level clusters (municipal)
k_true_D_4_2       = 3                                            # true number of D-level clusters (departmental)
mix_probs_E_4_2    = c(0.1, 0.45, 0.45)                           # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_4_2    = c(0.1, 0.45, 0.45)                           # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_4_2    = c(0.1, 0.45, 0.45)                           # vector of length k_true_D: D-level mixing probabilities (must sum to 1)
mu_centers_E_4_2   = seq(from = 1, by = 7, length.out = k_true_E_4_2)  # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_4_2   = seq(from = 1, by = 7, length.out = k_true_M_4_2)  # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_4_2   = seq(from = 1, by = 7, length.out = k_true_D_4_2)  # vector of length k_true_D: centers for beta_k_D atoms


simulated_data_11 <- generate_hierarchical_data_3DP(
  m                   = m_4,                          
  n_mun_per_dept      = n_mun_per_dept_4,             
  n_obs_range         = n_obs_range_4,                
  p_E                 = p_E_4,                        
  p_M                 = p_M_4,                        
  p_D                 = p_D_4,                       
  k_true_E            = k_true_E_4_2,                 
  k_true_M            = k_true_M_4_2,                   
  k_true_D            = k_true_D_4_2,                   
  mix_probs_E         = mix_probs_E_4_2,                 
  mix_probs_M         = mix_probs_M_4_2,                 
  mix_probs_D         = mix_probs_D_4_2,               
  mu_centers_E        = mu_centers_E_4_2,               
  mu_centers_M        = mu_centers_M_4_2,                
  mu_centers_D        = mu_centers_D_4_2,             
  sigma_k_E           = sigma_k_E_4,             
  sigma_k_M           = sigma_k_M_4,             
  sigma_k_D           = sigma_k_D_4,             
  kappa_min           = kappa_min_4,                    
  kappa_max           = kappa_max_4                    
)

# -------------------------------------------------
# 12.Three clusters with similar probability
# -------------------------------------------------
k_true_E_4_3       = 3                                            # true number of E-level clusters (individual)
k_true_M_4_3       = 3                                            # true number of M-level clusters (municipal)
k_true_D_4_3       = 3                                            # true number of D-level clusters (departmental)
mix_probs_E_4_3    = c(0.4, 0.25, 0.35)                           # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_4_3    = c(0.4, 0.25, 0.35)                           # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_4_3    = c(0.4, 0.25, 0.35)                           # vector of length k_true_D: D-level mixing probabilities (must sum to 1)
mu_centers_E_4_3   = seq(from = 1, by = 7, length.out = k_true_E_4_3)  # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_4_3   = seq(from = 1, by = 7, length.out = k_true_M_4_3)  # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_4_3   = seq(from = 1, by = 7, length.out = k_true_D_4_3)  # vector of length k_true_D: centers for beta_k_D atoms


simulated_data_12 <- generate_hierarchical_data_3DP(
  m                   = m_4,                          
  n_mun_per_dept      = n_mun_per_dept_4,             
  n_obs_range         = n_obs_range_4,                
  p_E                 = p_E_4,                        
  p_M                 = p_M_4,                        
  p_D                 = p_D_4,                       
  k_true_E            = k_true_E_4_3,                 
  k_true_M            = k_true_M_4_3,                   
  k_true_D            = k_true_D_4_3,                   
  mix_probs_E         = mix_probs_E_4_3,                 
  mix_probs_M         = mix_probs_M_4_3,                 
  mix_probs_D         = mix_probs_D_4_3,               
  mu_centers_E        = mu_centers_E_4_3,               
  mu_centers_M        = mu_centers_M_4_3,                
  mu_centers_D        = mu_centers_D_4_3,             
  sigma_k_E           = sigma_k_E_4,             
  sigma_k_M           = sigma_k_M_4,             
  sigma_k_D           = sigma_k_D_4,             
  kappa_min           = kappa_min_4,                    
  kappa_max           = kappa_max_4                    
)


# =============================
# 5. Number of observations n
# =============================

# In this scenario we keep a moderated number departments and of covariates
#we keep an homogeneous number of municipalities per department, and number observations per municipality
#. Also, keep a moderated number of clusters (6). Mixing probabilities are homogeneous across clusters

# The only change is the number of observations per municipality n_obs_range:
# - Low: 5-8
# - Moderated: 30-35
# - High: 80-100


# Constant parameters
m_5              = 32                          # number of departments
n_mun_per_dept_5 = rep(5, m_5)                 # vector of length m: municipalities per department
p_E_5            = 4                           # number of individual-level covariates
p_M_5            = 4                           # number of municipality-level covariates
p_D_5            = 4                           # number of department-level covariates
k_true_E_5       = 2                           # true number of E-level clusters (individual)
k_true_M_5       = 2                           # true number of M-level clusters (municipal)
k_true_D_5       = 2                           # true number of D-level clusters (departmental)
mix_probs_E_5    = rep(1/2, 2)                 # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_5    = rep(1/2, 2)                 # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_5    = rep(1/2, 2)                 # vector of length k_true_D: D-level mixing probabilities (must sum to 1)
mu_centers_E_5   = c(1, 7)                   # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_5   = c(1, 7)                   # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_5   = c(1, 7)                   # vector of length k_true_D: centers for beta_k_D atoms
sigma_k_E_5      = 1.0                         # sd for sampling beta_k_E coefficients around centers
sigma_k_M_5      = 1.0                         # sd for sampling beta_k_M coefficients around centers
sigma_k_D_5      = 1.0                         # sd for sampling beta_k_D coefficients around centers
kappa_min_5      = 2                           # min value of municipal variance (runif)              
kappa_max_5      = 5                           # max value of municipal variance (runif)




# -------------------------------------------------
# 13. Low number of observations n
# -------------------------------------------------
n_obs_range_5_1  = c(6, 8)                   # c(min, max): range for observations per municipality

simulated_data_13 <- generate_hierarchical_data_3DP(
  m                   = m_5,                          
  n_mun_per_dept      = n_mun_per_dept_5,             
  n_obs_range         = n_obs_range_5_1,                
  p_E                 = p_E_5,                        
  p_M                 = p_M_5,                        
  p_D                 = p_D_5,                       
  k_true_E            = k_true_E_5,                 
  k_true_M            = k_true_M_5,                   
  k_true_D            = k_true_D_5,                   
  mix_probs_E         = mix_probs_E_5,                 
  mix_probs_M         = mix_probs_M_5,                 
  mix_probs_D         = mix_probs_D_5,               
  mu_centers_E        = mu_centers_E_5,               
  mu_centers_M        = mu_centers_M_5,                
  mu_centers_D        = mu_centers_D_5,             
  sigma_k_E           = sigma_k_E_5,             
  sigma_k_M           = sigma_k_M_5,             
  sigma_k_D           = sigma_k_D_5,             
  kappa_min           = kappa_min_5,                    
  kappa_max           = kappa_max_5                    
)


# -------------------------------------------------
# 14.Moderated number of observations n
# -------------------------------------------------
n_obs_range_5_2    = c(12, 13)                   # c(min, max): range for observations per municipality

simulated_data_14 <- generate_hierarchical_data_3DP(
  m                   = m_5,                          
  n_mun_per_dept      = n_mun_per_dept_5,             
  n_obs_range         = n_obs_range_5_2,                
  p_E                 = p_E_5,                        
  p_M                 = p_M_5,                        
  p_D                 = p_D_5,                       
  k_true_E            = k_true_E_5,                 
  k_true_M            = k_true_M_5,                   
  k_true_D            = k_true_D_5,                   
  mix_probs_E         = mix_probs_E_5,                 
  mix_probs_M         = mix_probs_M_5,                 
  mix_probs_D         = mix_probs_D_5,               
  mu_centers_E        = mu_centers_E_5,               
  mu_centers_M        = mu_centers_M_5,                
  mu_centers_D        = mu_centers_D_5,             
  sigma_k_E           = sigma_k_E_5,             
  sigma_k_M           = sigma_k_M_5,             
  sigma_k_D           = sigma_k_D_5,             
  kappa_min           = kappa_min_5,                    
  kappa_max           = kappa_max_5                    
)


# -------------------------------------------------
# 15. High number of observations n
# -------------------------------------------------
n_obs_range_5_3    = c(16, 18)                   # c(min, max): range for observations per municipality

simulated_data_15 <- generate_hierarchical_data_3DP(
  m                   = m_5,                          
  n_mun_per_dept      = n_mun_per_dept_5,             
  n_obs_range         = n_obs_range_5_3,                
  p_E                 = p_E_5,                        
  p_M                 = p_M_5,                        
  p_D                 = p_D_5,                       
  k_true_E            = k_true_E_5,                 
  k_true_M            = k_true_M_5,                   
  k_true_D            = k_true_D_5,                   
  mix_probs_E         = mix_probs_E_5,                 
  mix_probs_M         = mix_probs_M_5,                 
  mix_probs_D         = mix_probs_D_5,               
  mu_centers_E        = mu_centers_E_5,               
  mu_centers_M        = mu_centers_M_5,                
  mu_centers_D        = mu_centers_D_5,             
  sigma_k_E           = sigma_k_E_5,             
  sigma_k_M           = sigma_k_M_5,             
  sigma_k_D           = sigma_k_D_5,             
  kappa_min           = kappa_min_5,                    
  kappa_max           = kappa_max_5                    
)


# =============================
# 6. NUMBER OF DEPARTMENTS
# =============================

# In this scenario we keep a moderated number of covariates and homogeneous number of
#municipalities per department and observations per municipality. Also, keep a moderated 
#number of clusters with homogeneous mixing probabilities. 

# The only change is the number of departments:
# - Low: 6
# - Moderated: 20
# - High: 40


# Constant parameters
p_E_6            = 4                           # number of individual-level covariates
p_M_6            = 4                           # number of municipality-level covariates
p_D_6            = 3                           # number of department-level covariates
k_true_E_6       = 2                           # true number of E-level clusters (individual)
k_true_M_6       = 2                           # true number of M-level clusters (municipal)
k_true_D_6       = 2                           # true number of D-level clusters (departmental)
mix_probs_E_6    = rep(1/2, 2)                 # vector of length k_true_E: E-level mixing probabilities (must sum to 1)
mix_probs_M_6    = rep(1/2, 2)                 # vector of length k_true_M: M-level mixing probabilities (must sum to 1)
mix_probs_D_6    = rep(1/2, 2)                 # vector of length k_true_D: D-level mixing probabilities (must sum to 1)
mu_centers_E_6   = c(1, 7)                   # vector of length k_true_E: centers for beta_k_E atoms
mu_centers_M_6   = c(1, 7)                   # vector of length k_true_M: centers for beta_k_M atoms
mu_centers_D_6   = c(1, 7)                   # vector of length k_true_D: centers for beta_k_D atoms
sigma_k_E_6      = 1.0                         # sd for sampling beta_k_E coefficients around centers
sigma_k_M_6      = 1.0                         # sd for sampling beta_k_M coefficients around centers
sigma_k_D_6      = 1.0                         # sd for sampling beta_k_D coefficients around centers
kappa_min_6      = 2                           # min value of municipal variance (runif)              
kappa_max_6      = 4                           # max value of municipal variance (runif)


# ----------------------------------
# 16. Low number of departments m scenario
# ----------------------------------
m_6_1              = 6                           # number of departments
n_mun_per_dept_6_1 = rep(5, m_6_1)                  # vector of length m: municipalities per department
n_obs_range_6_1    = c(20, 28)                   # c(min, max): range for observations per municipality


simulated_data_16 <- generate_hierarchical_data_3DP(
  m                   = m_6_1,                          
  n_mun_per_dept      = n_mun_per_dept_6_1,             
  n_obs_range         = n_obs_range_6_1,                
  p_E                 = p_E_6,                        
  p_M                 = p_M_6,                        
  p_D                 = p_D_6,                       
  k_true_E            = k_true_E_6,                 
  k_true_M            = k_true_M_6,                   
  k_true_D            = k_true_D_6,                   
  mix_probs_E         = mix_probs_E_6,                 
  mix_probs_M         = mix_probs_M_6,                 
  mix_probs_D         = mix_probs_D_6,               
  mu_centers_E        = mu_centers_E_6,               
  mu_centers_M        = mu_centers_M_6,                
  mu_centers_D        = mu_centers_D_6,             
  sigma_k_E           = sigma_k_E_6,             
  sigma_k_M           = sigma_k_M_6,             
  sigma_k_D           = sigma_k_D_6,             
  kappa_min           = kappa_min_6,                    
  kappa_max           = kappa_max_6                    
)

# ----------------------------------
# 17. Moderated number of departments m scenario
# ----------------------------------
m_6_2              = 30                           # number of departments
n_mun_per_dept_6_2 = rep(7, m_6_2)                # vector of length m: municipalities per department
n_obs_range_6_2    = c(7, 9)                   # c(min, max): range for observations per municipality

simulated_data_17 <- generate_hierarchical_data_3DP(
  m                   = m_6_2,                          
  n_mun_per_dept      = n_mun_per_dept_6_2,             
  n_obs_range         = n_obs_range_6_2,                
  p_E                 = p_E_6,                        
  p_M                 = p_M_6,                        
  p_D                 = p_D_6,                       
  k_true_E            = k_true_E_6,                 
  k_true_M            = k_true_M_6,                   
  k_true_D            = k_true_D_6,                   
  mix_probs_E         = mix_probs_E_6,                 
  mix_probs_M         = mix_probs_M_6,                 
  mix_probs_D         = mix_probs_D_6,               
  mu_centers_E        = mu_centers_E_6,               
  mu_centers_M        = mu_centers_M_6,                
  mu_centers_D        = mu_centers_D_6,             
  sigma_k_E           = sigma_k_E_6,             
  sigma_k_M           = sigma_k_M_6,             
  sigma_k_D           = sigma_k_D_6,             
  kappa_min           = kappa_min_6,                    
  kappa_max           = kappa_max_6                    
)

# ----------------------------------
# 18. High number of departments m scenario
# ----------------------------------
m_6_3              = 50                           # number of departments
n_mun_per_dept_6_3 = rep(4, m_6_3)                # vector of length m: municipalities per department
n_obs_range_6_3    = c(12, 18)                   # c(min, max): range for observations per municipality


simulated_data_18 <- generate_hierarchical_data_3DP(
  m                   = m_6_3,                          
  n_mun_per_dept      = n_mun_per_dept_6_3,             
  n_obs_range         = n_obs_range_6_3,                
  p_E                 = p_E_6,                        
  p_M                 = p_M_6,                        
  p_D                 = p_D_6,                       
  k_true_E            = k_true_E_6,                 
  k_true_M            = k_true_M_6,                   
  k_true_D            = k_true_D_6,                   
  mix_probs_E         = mix_probs_E_6,                 
  mix_probs_M         = mix_probs_M_6,                 
  mix_probs_D         = mix_probs_D_6,               
  mu_centers_E        = mu_centers_E_6,               
  mu_centers_M        = mu_centers_M_6,                
  mu_centers_D        = mu_centers_D_6,             
  sigma_k_E           = sigma_k_E_6,             
  sigma_k_M           = sigma_k_M_6,             
  sigma_k_D           = sigma_k_D_6,             
  kappa_min           = kappa_min_6,                    
  kappa_max           = kappa_max_6                    
)



# ==============================================================================================
#                     TEST SIMULATED SCENARIOS (~30% of training size)
# ==============================================================================================


#-----------------------------------------------------------------------------
# Multiply n_obs_range by 0.3 and ceiling, keeping all other
# parameters identical to the training datasets.
#
# Training total obs = n_mun_total × mean(n_obs_range)
# Test    total obs = n_mun_total × mean(n_obs_range × 0.3)

#-----------------------------------------------------------------------------

# =============================
# 1. VARIANCE SCENARIOS
# =============================

# Training: n_obs_range_1    = c(14, 16)

# ----------------------------------
# Test 1. Low municipal variance
# ----------------------------------

test_data_1 <- generate_hierarchical_data_3DP(
  m                   = m_1,                          
  n_mun_per_dept      = n_mun_per_dept_1,             
  n_obs_range         = ceiling(n_obs_range_1 * 0.3),                
  p_E                 = p_E_1,                        
  p_M                 = p_M_1,                        
  p_D                 = p_D_1,                       
  k_true_E            = k_true_E_1,                 
  k_true_M            = k_true_M_1,                   
  k_true_D            = k_true_D_1,                   
  mix_probs_E         = mix_probs_E_1,                 
  mix_probs_M         = mix_probs_M_1,                 
  mix_probs_D         = mix_probs_D_1,               
  mu_centers_E        = mu_centers_E_1,               
  mu_centers_M        = mu_centers_M_1,                
  mu_centers_D        = mu_centers_D_1,             
  sigma_k_E           = sigma_k_E_1,             
  sigma_k_M           = sigma_k_M_1,             
  sigma_k_D           = sigma_k_D_1,             
  kappa_min           = kappa_min_1_1,                    
  kappa_max           = kappa_max_1_1                    
)

# ----------------------------------
# Test 2. Moderated municipal variance scenario
# ----------------------------------

test_data_2 <- generate_hierarchical_data_3DP(
  m                   = m_1,                          
  n_mun_per_dept      = n_mun_per_dept_1,             
  n_obs_range         = ceiling(n_obs_range_1 * 0.3),                
  p_E                 = p_E_1,                        
  p_M                 = p_M_1,                        
  p_D                 = p_D_1,                       
  k_true_E            = k_true_E_1,                 
  k_true_M            = k_true_M_1,                   
  k_true_D            = k_true_D_1,                   
  mix_probs_E         = mix_probs_E_1,                 
  mix_probs_M         = mix_probs_M_1,                 
  mix_probs_D         = mix_probs_D_1,               
  mu_centers_E        = mu_centers_E_1,               
  mu_centers_M        = mu_centers_M_1,                
  mu_centers_D        = mu_centers_D_1,             
  sigma_k_E           = sigma_k_E_1,             
  sigma_k_M           = sigma_k_M_1,             
  sigma_k_D           = sigma_k_D_1,             
  kappa_min           = kappa_min_1_2,                    
  kappa_max           = kappa_max_1_2                    
)


# ----------------------------------
# Test 3. High municipal variance scenario
# ----------------------------------

test_data_3 <- generate_hierarchical_data_3DP(
  m                   = m_1,                          
  n_mun_per_dept      = n_mun_per_dept_1,             
  n_obs_range         = ceiling(n_obs_range_1 * 0.3),                
  p_E                 = p_E_1,                        
  p_M                 = p_M_1,                        
  p_D                 = p_D_1,                       
  k_true_E            = k_true_E_1,                 
  k_true_M            = k_true_M_1,                   
  k_true_D            = k_true_D_1,                   
  mix_probs_E         = mix_probs_E_1,                 
  mix_probs_M         = mix_probs_M_1,                 
  mix_probs_D         = mix_probs_D_1,               
  mu_centers_E        = mu_centers_E_1,               
  mu_centers_M        = mu_centers_M_1,                
  mu_centers_D        = mu_centers_D_1,             
  sigma_k_E           = sigma_k_E_1,             
  sigma_k_M           = sigma_k_M_1,             
  sigma_k_D           = sigma_k_D_1,             
  kappa_min           = kappa_min_1_3,                    
  kappa_max           = kappa_max_1_3                    
)



# =============================
# 2. NUMBER OF CLUSTER K SCENARIOS
# =============================

# Training: n_obs_range_2    = c(15, 20)

# ----------------------------------
# Test 4. Low number of clusters k 
# ----------------------------------

test_data_4 <- generate_hierarchical_data_3DP(
  m                   = m_2,                          
  n_mun_per_dept      = n_mun_per_dept_2,             
  n_obs_range         = ceiling(n_obs_range_2 * 0.3),                
  p_E                 = p_E_2,                        
  p_M                 = p_M_2,                        
  p_D                 = p_D_2,                       
  k_true_E            = k_true_E_2_1,                 
  k_true_M            = k_true_M_2_1,                   
  k_true_D            = k_true_D_2_1,                   
  mix_probs_E         = mix_probs_E_2_1,                 
  mix_probs_M         = mix_probs_M_2_1,                 
  mix_probs_D         = mix_probs_D_2_1,               
  mu_centers_E        = mu_centers_E_2_1,               
  mu_centers_M        = mu_centers_M_2_1,                
  mu_centers_D        = mu_centers_D_2_1,             
  sigma_k_E           = sigma_k_E_2,             
  sigma_k_M           = sigma_k_M_2,             
  sigma_k_D           = sigma_k_D_2,             
  kappa_min           = kappa_min_2,                    
  kappa_max           = kappa_max_2                    
)


# ----------------------------------
# Test 5. Moderated number of clusters k 
# ----------------------------------

test_data_5 <- generate_hierarchical_data_3DP(
  m                   = m_2,                          
  n_mun_per_dept      = n_mun_per_dept_2,             
  n_obs_range         = ceiling(n_obs_range_2 * 0.3),                
  p_E                 = p_E_2,                        
  p_M                 = p_M_2,                        
  p_D                 = p_D_2,                       
  k_true_E            = k_true_E_2_2,                 
  k_true_M            = k_true_M_2_2,                   
  k_true_D            = k_true_D_2_2,                   
  mix_probs_E         = mix_probs_E_2_2,                 
  mix_probs_M         = mix_probs_M_2_2,                 
  mix_probs_D         = mix_probs_D_2_2,               
  mu_centers_E        = mu_centers_E_2_2,               
  mu_centers_M        = mu_centers_M_2_2,                
  mu_centers_D        = mu_centers_D_2_2,             
  sigma_k_E           = sigma_k_E_2,             
  sigma_k_M           = sigma_k_M_2,             
  sigma_k_D           = sigma_k_D_2,             
  kappa_min           = kappa_min_2,                    
  kappa_max           = kappa_max_2                    
)

# ----------------------------------
# Test 6. High number of clusters k 
# ----------------------------------

test_data_6 <- generate_hierarchical_data_3DP(
  m                   = m_2,                          
  n_mun_per_dept      = n_mun_per_dept_2,             
  n_obs_range         = ceiling(n_obs_range_2 * 0.3),                
  p_E                 = p_E_2,                        
  p_M                 = p_M_2,                        
  p_D                 = p_D_2,                       
  k_true_E            = k_true_E_2_3,                 
  k_true_M            = k_true_M_2_3,                   
  k_true_D            = k_true_D_2_3,                   
  mix_probs_E         = mix_probs_E_2_3,                 
  mix_probs_M         = mix_probs_M_2_3,                 
  mix_probs_D         = mix_probs_D_2_3,               
  mu_centers_E        = mu_centers_E_2_3,               
  mu_centers_M        = mu_centers_M_2_3,                
  mu_centers_D        = mu_centers_D_2_3,             
  sigma_k_E           = sigma_k_E_2,             
  sigma_k_M           = sigma_k_M_2,             
  sigma_k_D           = sigma_k_D_2,             
  kappa_min           = kappa_min_2,                    
  kappa_max           = kappa_max_2                    
)



# =============================
# 3. NUMBER OF COVARIATES p
# =============================

#Training: n_obs_range_3    = c(15, 18)                 

# ----------------------------------
# Test 7. Low number of covariates p 
# ----------------------------------

test_data_7 <- generate_hierarchical_data_3DP(
  m                   = m_3,                          
  n_mun_per_dept      = n_mun_per_dept_3,             
  n_obs_range         = ceiling(n_obs_range_3 * 0.3),                
  p_E                 = p_E_3_1,                        
  p_M                 = p_M_3_1,                        
  p_D                 = p_D_3_1,                       
  k_true_E            = k_true_E_3,                 
  k_true_M            = k_true_M_3,                   
  k_true_D            = k_true_D_3,                   
  mix_probs_E         = mix_probs_E_3,                 
  mix_probs_M         = mix_probs_M_3,                 
  mix_probs_D         = mix_probs_D_3,               
  mu_centers_E        = mu_centers_E_3,               
  mu_centers_M        = mu_centers_M_3,                
  mu_centers_D        = mu_centers_D_3,             
  sigma_k_E           = sigma_k_E_3,             
  sigma_k_M           = sigma_k_M_3,             
  sigma_k_D           = sigma_k_D_3,             
  kappa_min           = kappa_min_3,                    
  kappa_max           = kappa_max_3                    
)

# ----------------------------------
# Test 8. Moderated number of covariates p 
# ----------------------------------

test_data_8 <- generate_hierarchical_data_3DP(
  m                   = m_3,                          
  n_mun_per_dept      = n_mun_per_dept_3,             
  n_obs_range         = ceiling(n_obs_range_3 * 0.3),                
  p_E                 = p_E_3_2,                        
  p_M                 = p_M_3_2,                        
  p_D                 = p_D_3_2,                       
  k_true_E            = k_true_E_3,                 
  k_true_M            = k_true_M_3,                   
  k_true_D            = k_true_D_3,                   
  mix_probs_E         = mix_probs_E_3,                 
  mix_probs_M         = mix_probs_M_3,                 
  mix_probs_D         = mix_probs_D_3,               
  mu_centers_E        = mu_centers_E_3,               
  mu_centers_M        = mu_centers_M_3,                
  mu_centers_D        = mu_centers_D_3,             
  sigma_k_E           = sigma_k_E_3,             
  sigma_k_M           = sigma_k_M_3,             
  sigma_k_D           = sigma_k_D_3,             
  kappa_min           = kappa_min_3,                    
  kappa_max           = kappa_max_3                    
)


# ----------------------------------
# Test 9. High number of covariates p 
# ----------------------------------

test_data_9 <- generate_hierarchical_data_3DP(
  m                   = m_3,                          
  n_mun_per_dept      = n_mun_per_dept_3,             
  n_obs_range         = ceiling(n_obs_range_3 * 0.3),                
  p_E                 = p_E_3_3,                        
  p_M                 = p_M_3_3,                        
  p_D                 = p_D_3_3,                       
  k_true_E            = k_true_E_3,                 
  k_true_M            = k_true_M_3,                   
  k_true_D            = k_true_D_3,                   
  mix_probs_E         = mix_probs_E_3,                 
  mix_probs_M         = mix_probs_M_3,                 
  mix_probs_D         = mix_probs_D_3,               
  mu_centers_E        = mu_centers_E_3,               
  mu_centers_M        = mu_centers_M_3,                
  mu_centers_D        = mu_centers_D_3,             
  sigma_k_E           = sigma_k_E_3,             
  sigma_k_M           = sigma_k_M_3,             
  sigma_k_D           = sigma_k_D_3,             
  kappa_min           = kappa_min_3,                    
  kappa_max           = kappa_max_3                    
)


# =============================
# 4. MIXING PROBABILITIES
# =============================

# Training: n_obs_range_4    = c(15, 18)

# -------------------------------------------------
# Test 10. One cluster with high probability
# -------------------------------------------------

test_data_10 <- generate_hierarchical_data_3DP(
  m                   = m_4,                          
  n_mun_per_dept      = n_mun_per_dept_4,             
  n_obs_range         = ceiling(n_obs_range_4 * 0.3),                
  p_E                 = p_E_4,                        
  p_M                 = p_M_4,                        
  p_D                 = p_D_4,                       
  k_true_E            = k_true_E_4_1,                 
  k_true_M            = k_true_M_4_1,                   
  k_true_D            = k_true_D_4_1,                   
  mix_probs_E         = mix_probs_E_4_1,                 
  mix_probs_M         = mix_probs_M_4_1,                 
  mix_probs_D         = mix_probs_D_4_1,               
  mu_centers_E        = mu_centers_E_4_1,               
  mu_centers_M        = mu_centers_M_4_1,                
  mu_centers_D        = mu_centers_D_4_1,             
  sigma_k_E           = sigma_k_E_4,             
  sigma_k_M           = sigma_k_M_4,             
  sigma_k_D           = sigma_k_D_4,             
  kappa_min           = kappa_min_4,                    
  kappa_max           = kappa_max_4                    
)


# -------------------------------------------------
# Test 11. One cluster with low probability
# -------------------------------------------------

test_data_11 <- generate_hierarchical_data_3DP(
  m                   = m_4,                          
  n_mun_per_dept      = n_mun_per_dept_4,             
  n_obs_range         = ceiling(n_obs_range_4 * 0.3),                
  p_E                 = p_E_4,                        
  p_M                 = p_M_4,                        
  p_D                 = p_D_4,                       
  k_true_E            = k_true_E_4_2,                 
  k_true_M            = k_true_M_4_2,                   
  k_true_D            = k_true_D_4_2,                   
  mix_probs_E         = mix_probs_E_4_2,                 
  mix_probs_M         = mix_probs_M_4_2,                 
  mix_probs_D         = mix_probs_D_4_2,               
  mu_centers_E        = mu_centers_E_4_2,               
  mu_centers_M        = mu_centers_M_4_2,                
  mu_centers_D        = mu_centers_D_4_2,             
  sigma_k_E           = sigma_k_E_4,             
  sigma_k_M           = sigma_k_M_4,             
  sigma_k_D           = sigma_k_D_4,             
  kappa_min           = kappa_min_4,                    
  kappa_max           = kappa_max_4                    
)

# -------------------------------------------------
# Test 12.Two clusters with high and two with low probability
# -------------------------------------------------

test_data_12 <- generate_hierarchical_data_3DP(
  m                   = m_4,                          
  n_mun_per_dept      = n_mun_per_dept_4,             
  n_obs_range         = ceiling(n_obs_range_4 * 0.3),                
  p_E                 = p_E_4,                        
  p_M                 = p_M_4,                        
  p_D                 = p_D_4,                       
  k_true_E            = k_true_E_4_3,                 
  k_true_M            = k_true_M_4_3,                   
  k_true_D            = k_true_D_4_3,                   
  mix_probs_E         = mix_probs_E_4_3,                 
  mix_probs_M         = mix_probs_M_4_3,                 
  mix_probs_D         = mix_probs_D_4_3,               
  mu_centers_E        = mu_centers_E_4_3,               
  mu_centers_M        = mu_centers_M_4_3,                
  mu_centers_D        = mu_centers_D_4_3,             
  sigma_k_E           = sigma_k_E_4,             
  sigma_k_M           = sigma_k_M_4,             
  sigma_k_D           = sigma_k_D_4,             
  kappa_min           = kappa_min_4,                    
  kappa_max           = kappa_max_4                    
)


# =============================
# 5. Number of observations n
# =============================


# -------------------------------------------------
# Test 13. Low number of observations n
# -------------------------------------------------
# Training: n_obs_range_5_1  = c(5, 8)  

test_data_13 <- generate_hierarchical_data_3DP(
  m                   = m_5,                          
  n_mun_per_dept      = n_mun_per_dept_5,             
  n_obs_range         = ceiling(n_obs_range_5_1 * 0.3),                
  p_E                 = p_E_5,                        
  p_M                 = p_M_5,                        
  p_D                 = p_D_5,                       
  k_true_E            = k_true_E_5,                 
  k_true_M            = k_true_M_5,                   
  k_true_D            = k_true_D_5,                   
  mix_probs_E         = mix_probs_E_5,                 
  mix_probs_M         = mix_probs_M_5,                 
  mix_probs_D         = mix_probs_D_5,               
  mu_centers_E        = mu_centers_E_5,               
  mu_centers_M        = mu_centers_M_5,                
  mu_centers_D        = mu_centers_D_5,             
  sigma_k_E           = sigma_k_E_5,             
  sigma_k_M           = sigma_k_M_5,             
  sigma_k_D           = sigma_k_D_5,             
  kappa_min           = kappa_min_5,                    
  kappa_max           = kappa_max_5                    
)


# -------------------------------------------------
# Test 14.Moderated number of observations n
# -------------------------------------------------
# Training: n_obs_range_5_2    = c(30, 35)

test_data_14 <- generate_hierarchical_data_3DP(
  m                   = m_5,                          
  n_mun_per_dept      = n_mun_per_dept_5,             
  n_obs_range         = ceiling(n_obs_range_5_2 * 0.3),                
  p_E                 = p_E_5,                        
  p_M                 = p_M_5,                        
  p_D                 = p_D_5,                       
  k_true_E            = k_true_E_5,                 
  k_true_M            = k_true_M_5,                   
  k_true_D            = k_true_D_5,                   
  mix_probs_E         = mix_probs_E_5,                 
  mix_probs_M         = mix_probs_M_5,                 
  mix_probs_D         = mix_probs_D_5,               
  mu_centers_E        = mu_centers_E_5,               
  mu_centers_M        = mu_centers_M_5,                
  mu_centers_D        = mu_centers_D_5,             
  sigma_k_E           = sigma_k_E_5,             
  sigma_k_M           = sigma_k_M_5,             
  sigma_k_D           = sigma_k_D_5,             
  kappa_min           = kappa_min_5,                    
  kappa_max           = kappa_max_5                    
)


# -------------------------------------------------
# Test 15. High number of observations n
# -------------------------------------------------
# Training: n_obs_range_5_3    = c(80, 100)                  

test_data_15 <- generate_hierarchical_data_3DP(
  m                   = m_5,                          
  n_mun_per_dept      = n_mun_per_dept_5,             
  n_obs_range         = ceiling(n_obs_range_5_3 * 0.3),                
  p_E                 = p_E_5,                        
  p_M                 = p_M_5,                        
  p_D                 = p_D_5,                       
  k_true_E            = k_true_E_5,                 
  k_true_M            = k_true_M_5,                   
  k_true_D            = k_true_D_5,                   
  mix_probs_E         = mix_probs_E_5,                 
  mix_probs_M         = mix_probs_M_5,                 
  mix_probs_D         = mix_probs_D_5,               
  mu_centers_E        = mu_centers_E_5,               
  mu_centers_M        = mu_centers_M_5,                
  mu_centers_D        = mu_centers_D_5,             
  sigma_k_E           = sigma_k_E_5,             
  sigma_k_M           = sigma_k_M_5,             
  sigma_k_D           = sigma_k_D_5,             
  kappa_min           = kappa_min_5,                    
  kappa_max           = kappa_max_5                    
)



# =============================
# 6. NUMBER OF DEPARTMENTS
# =============================

# ----------------------------------
# Test 16. Low number of departments m scenario
# ----------------------------------

test_data_16 <- generate_hierarchical_data_3DP(
  m                   = m_6_1,                          
  n_mun_per_dept      = n_mun_per_dept_6_1,             
  n_obs_range         = ceiling(n_obs_range_6_1 * 0.3),                
  p_E                 = p_E_6,                        
  p_M                 = p_M_6,                        
  p_D                 = p_D_6,                       
  k_true_E            = k_true_E_6,                 
  k_true_M            = k_true_M_6,                   
  k_true_D            = k_true_D_6,                   
  mix_probs_E         = mix_probs_E_6,                 
  mix_probs_M         = mix_probs_M_6,                 
  mix_probs_D         = mix_probs_D_6,               
  mu_centers_E        = mu_centers_E_6,               
  mu_centers_M        = mu_centers_M_6,                
  mu_centers_D        = mu_centers_D_6,             
  sigma_k_E           = sigma_k_E_6,             
  sigma_k_M           = sigma_k_M_6,             
  sigma_k_D           = sigma_k_D_6,             
  kappa_min           = kappa_min_6,                    
  kappa_max           = kappa_max_6                    
)

# ----------------------------------
# Test 17. Moderated number of departments m scenario
# ----------------------------------


test_data_17 <- generate_hierarchical_data_3DP(
  m                   = m_6_2,                          
  n_mun_per_dept      = n_mun_per_dept_6_2,             
  n_obs_range         = ceiling(n_obs_range_6_2 * 0.3),                
  p_E                 = p_E_6,                        
  p_M                 = p_M_6,                        
  p_D                 = p_D_6,                       
  k_true_E            = k_true_E_6,                 
  k_true_M            = k_true_M_6,                   
  k_true_D            = k_true_D_6,                   
  mix_probs_E         = mix_probs_E_6,                 
  mix_probs_M         = mix_probs_M_6,                 
  mix_probs_D         = mix_probs_D_6,               
  mu_centers_E        = mu_centers_E_6,               
  mu_centers_M        = mu_centers_M_6,                
  mu_centers_D        = mu_centers_D_6,             
  sigma_k_E           = sigma_k_E_6,             
  sigma_k_M           = sigma_k_M_6,             
  sigma_k_D           = sigma_k_D_6,             
  kappa_min           = kappa_min_6,                    
  kappa_max           = kappa_max_6                    
)

# ----------------------------------
# Test 18. High number of departments m scenario
# ----------------------------------

test_data_18 <- generate_hierarchical_data_3DP(
  m                   = m_6_3,                          
  n_mun_per_dept      = n_mun_per_dept_6_3,             
  n_obs_range         = ceiling(n_obs_range_6_3 * 0.3),                
  p_E                 = p_E_6,                        
  p_M                 = p_M_6,                        
  p_D                 = p_D_6,                       
  k_true_E            = k_true_E_6,                 
  k_true_M            = k_true_M_6,                   
  k_true_D            = k_true_D_6,                   
  mix_probs_E         = mix_probs_E_6,                 
  mix_probs_M         = mix_probs_M_6,                 
  mix_probs_D         = mix_probs_D_6,               
  mu_centers_E        = mu_centers_E_6,               
  mu_centers_M        = mu_centers_M_6,                
  mu_centers_D        = mu_centers_D_6,             
  sigma_k_E           = sigma_k_E_6,             
  sigma_k_M           = sigma_k_M_6,             
  sigma_k_D           = sigma_k_D_6,             
  kappa_min           = kappa_min_6,                    
  kappa_max           = kappa_max_6                    
)


# ==============================================================================================
#@                            PREPARATION FOR THE GIBBS
# ==============================================================================================
#------------------------------------------------
#      Parameters for the GIBBS
#------------------------------------------------

B               = 20000 

nu_beta         = 1

a_alpha         = 1
b_alpha         = 1
alpha           = 1


sigma2_base     = 5
a_sigma2        = 1


for (s in 1:18) {
  
  simulated_df_name <- sprintf("simulated_data_%d", s)
  obj <- get(simulated_df_name, envir = .GlobalEnv)
  
  #------------------------------------------------
  #      Retrieve information from the database
  #------------------------------------------------
  
  #Dataframe extraction
  df          <- obj$datos
  x_names     <- obj$x_names
  Y           <- obj$Y
  y           <- obj$y
  m           <- obj$m
  X_full      <- obj$X_full
  mun_kappa   <- obj$mun_kappa
  sqrt_kappa  <- obj$sqrt_kappa
  xTx_vec     <- obj$xTx_full_vec        
  xi          <- obj$xi_false            #xi to initialize (Not the real xi cluster assignments)
  beta_k      <- obj$beta_k_false        #beta_k to initialize (Not the real beta_k atoms)
  beta_int    <- obj$beta_int_false
  
  
  p           <- obj$p
  obs_to_mun  <- obj$obs_to_mun
  mun_map     <- obj$mun_map
  kappa2_q    <- obj$kappa2_q
  dept_to_mun <- obj$dept_to_mun
  n_mun_total <- obj$n_mun_total
  
  mu_beta     <- beta_int
  
  # Unitary priors hiperparameter specification
  lm_previa_int       <-  lm(y ~ 1)
  lm_previa_unitaria  <-  lm(y ~ X_full -1)
  
  sigma2_beta         <-  summary(lm_previa_int)$sigma
  gamma_beta          <-  summary(lm_previa_int)$sigma
  b_sigma2            <-  summary(lm_previa_int)$sigma
  
  a_b_sigma2      = a_sigma2*b_sigma2
  
  sigma2_beta         <-  summary(lm_previa_unitaria)$sigma
  
  
  nu2_mu              <- summary(lm_previa_unitaria)$sigma
  nu2_mu_inv          <- 1 / nu2_mu
  
  eta_mu              <- coef(lm_previa_unitaria)
  mu_vec              <- coef(lm_previa_unitaria)
  
  # VARIANCE COMPONENT HIPERPARAMETERS
  # using CV (1) and var(y)
  alpha_kappa      <-   max(0.1, 2 / (1^2))
  beta_kappa       <-   max(0.1, alpha_kappa / var(y))
  
  nu_kappa         <-   4 + (2 / 1^2)
  
  #------------------------------------------------
  #      Execute the GIBBS sampler
  #------------------------------------------------
  tictoc::tic()
  cadena <- MCMC2_BNP(
    B = B,
    nu_beta = nu_beta,
    gamma_beta = gamma_beta,
    nu_kappa = nu_kappa,
    alpha_kappa = alpha_kappa,
    beta_kappa = beta_kappa
  )
  
  tiempo = tictoc::toc()
  cadena$info$Tiempo = tiempo$callback_msg
  
  # ======================================
  # LOAD TEXT FILES
  # ======================================
  
  #-------------------------------
  # load xi samples
  #-------------------------------
  
  xi_file <- file.path(path_text_files, "xi_samples.txt")
  
  xi_samples <- as.matrix(read.table(xi_file, header = FALSE))
  
  cat("Xi loaded:", nrow(xi_samples), "samples ×", ncol(xi_samples), "observations\n\n")
  
  #-------------------------------
  # Load beta_k samples line by line
  #-------------------------------
  
  beta_k_file <- file.path(path_text_files, "beta_k_samples.txt")
  
  # Count lines first: Functions to create, open and close connections, i.e., “generalized files”,
  #such as possibly compressed files, URLs, pipes, etc.
  
  con <- file(beta_k_file, "r")
  lines <- readLines(con)
  close(con)
  
  n_samples <- length(lines)
  cat("Found", n_samples, "samples\n")
  
  # Read each line and parse
  beta_k_samples <- vector("list", n_samples)
  
  for (i in 1:n_samples) {
    # Split line by spaces and convert to numeric
    values <- as.numeric(strsplit(lines[i], "\\s+")[[1]])
    
    K <- as.integer(values[1])          # First element is number of clusters
    betas <- values[2:(K * p + 1)]      # Remaining elements are beta vectors
    
    # Reconstruct list of K vectors
    beta_k_samples[[i]] <- vector("list", K)
    for (k in 1:K) {
      start <- (k - 1) * p + 1
      end   <- k * p
      beta_k_samples[[i]][[k]] <- betas[start:end]
    }
    
    # Progress indicator
    #if (i %% 100 == 0) cat("  Processed", i, "/", n_samples, "samples\n")
  }
  
  cat("beta_k loaded:", length(beta_k_samples), "samples\n")
  
  
  # ==================================
  # Add to the chain results
  # ==================================
  
  cadena$xi <- xi_samples
  cadena$beta_k <- beta_k_samples
  
  #-------------------------------
  # Verify
  #-------------------------------
  cat("Summary:\n")
  cat("  Number of samples:", nrow(xi_samples), "\n")
  cat("  Number of observations:", ncol(xi_samples), "\n")
  cat("  Number of covariates:", p, "\n")
  cat("  Sample 1 has", length(beta_k_samples[[1]]), "clusters\n")
  cat("  Sample", n_samples, "has", length(beta_k_samples[[n_samples]]), "clusters\n\n")
  
  
  # ==============================================================================
  #@                            LOG-LIKELIHOOD COMPUTATION
  # ==============================================================================
  
  cat("Computing log-likelihood for", B, "samples...\n")
  cat("Validating data structures...\n")
  
  # ===================================================================
  #    VALIDATE DATA STRUCTURES
  # ===================================================================
  
  # Check that we have the right dimensions
  stopifnot(nrow(cadena$xi) == B)
  stopifnot(ncol(cadena$xi) == length(y))
  stopifnot(length(cadena$beta_k) == B)
  stopifnot(length(cadena$beta_int) == B)
  
  cat("xi dimensions:", nrow(cadena$xi), "samples ×", ncol(cadena$xi), "observations\n")
  cat("beta_k samples:", length(cadena$beta_k), "\n")
  cat("beta_int samples:", length(cadena$beta_int), "\n\n")
  
  # ============================================
  # Compute Log-likelihood for each sample
  # ============================================
  
  log_lik_vec <- numeric(B)
  
  for (sample_count in seq_len(B)) {
    
    # Progress message
    if (sample_count %% 100 == 0 || sample_count == 1 || sample_count == B) {
      cat("Iteración", sample_count, "de", B, "\n")
      flush.console()
    }
    
    # ------------------------------------------
    # Extract parameters for the sample 
    # ------------------------------------------
    
    xi_sample <- cadena$xi[sample_count, ]          # Cluster assignments (length n)
    beta_k_sample <- cadena$beta_k[[sample_count]]  # List of K beta vectors
    beta_int_sample <- cadena$beta_int[sample_count] # Global intercept beta_0
    
    K_sample <- length(beta_k_sample)  # Number of clusters in this sample
    
    # ----------------------------------------------------
    # Compute fitted values for al observations as vector
    # ----------------------------------------------------
    
    # Stack all beta_k vectors into a matrix: K x p
    # Row k contains the coefficients β_k for cluster k
    
    # do.call function: do.call constructs and executes a function call from a name or 
    #a function and a list of arguments to be passed to it.
    beta_mat <- do.call(rbind, beta_k_sample)
    
    
    # Extract beta for each observation based on cluster assignment
    # beta_per_obs[i, ] = beta_{k(i)}: the coefficients for observation i's cluster
    beta_per_obs <- beta_mat[xi_sample, , drop = FALSE]  # n x p matrix
    
    # Compute fitted values: vartheta_i = beta_0 + x_i^T beta_{k(i)} for all observations
    # rowSums(X_full * beta_per_obs) computes x_i^T beta_{k(i)} for each i
    fitted_values <- beta_int_sample + rowSums(X_full * beta_per_obs)    #lenght = n
    
    # ---------------------------------------
    # Compute log-likelihood by municipality
    # ---------------------------------------
    
    total_ll <- 0  # Initialize total log-likelihood for this sample
    
    for (q in 1:m) {
      # Number of municipalities in department q
      n_mun_q <- length(Y[[q]])
      
      for (j_local in 1:n_mun_q) {   #iterates over municipalties in q
        
        # Get global municipality index
        mun_idx <- dept_to_mun[[q]][j_local]
        
        # Get observation i indices for this municipality
        idx_jq <- mun_map[[mun_idx]]$start:mun_map[[mun_idx]]$end
        
        # Extract response and fitted values for this municipality j
        y_jq <- y[idx_jq]
        fitted_jq <- fitted_values[idx_jq]
        
        # Get kappa2_{j,q} for this sample (From the list)
        kappa2_jq <- cadena$kappa2_jq[[q]][[j_local]][sample_count]
        
        # Compute log-likelihood contribution for this municipality:
        # it is computed as the sum of individual ll of the municipality
        # ll_{j,q} = sum_{i ∈ I_{j,q}} log N(y_i | beta0 + x_i^T beta_{k(i)}, kappa2_{j,q})
        ll_jq <- sum(dnorm(x = y_jq,
                           mean = fitted_jq,
                           sd = sqrt(kappa2_jq),
                           log = TRUE))           #log scale
        
        # Add to total. Sum of all loglikelihood
        total_ll <- total_ll + ll_jq
      }
    }
    
    # Store log-likelihood for this sample as vector
    log_lik_vec[sample_count] <- total_ll
  }
  
  # ===================
  # Store to the chain
  # ==================
  
  cadena$log_likelihood <- log_lik_vec
  
  
  # Name of the stored file
  results_file_name <- paste0("/Chain_Simulation_", s ,".RData")
  
  save(cadena, file = paste0(path_to_results, results_file_name))
  
  cat("----- Simulacion", s, "finalizada ------", "\n")
}



# ==============================================================================================
#                                     LOAD MCMC CHAINS
# ==============================================================================================

chains_executed <- 18


for (s in 1:chains_executed) {
  
  results_file_name <- file.path(path_to_results, paste0("Chain_Simulation_", s, ".RData"))
  
  # Load into a temporary environment to avoid overwriting
  temp_env <- new.env()
  load(file = results_file_name, envir = temp_env)
  
  # Assign each chain to a uniquely named object
  assign(paste0("cadena_", s), temp_env$cadena)
  
  cat("Loaded:", results_file_name, "\n")
}


# ==============================================================================================
#@                                     PLOTS
# ==============================================================================================

cadenas <- (1:18)

for (s in cadenas) {
  
  #the current MCMC chain
  cadena <- get(paste0("cadena_",s))
  data <- get(paste0("simulated_data_",s))
  
  # ===============================================
  #@             LOGLIKELIHOOD PLOTS
  # ===============================================
  
  B <- nrow(cadena$xi)
  
  cat("Summary statistics:\n")
  cat("  Min log-lik:  ", sprintf("%.2f", min(cadena$log_likelihood)), "\n")
  cat("  Max log-lik:  ", sprintf("%.2f", max(cadena$log_likelihood)), "\n")
  cat("  Mean log-lik: ", sprintf("%.2f", mean(cadena$log_likelihood)), "\n")
  cat("  Median log-lik:", sprintf("%.2f", median(cadena$log_likelihood)), "\n")
  cat("  SD log-lik:   ", sprintf("%.2f", sd(cadena$log_likelihood)), "\n\n")
  
  # Traceplot 1
  
  df_loglik <- data.frame(
    iteration = 1:B,
    log_likelihood = cadena$log_likelihood
  )
  
  ll_trace <- ggplot(df_loglik, aes(x = iteration, y = log_likelihood)) +
    geom_line(color = "steelblue", alpha = 0.4) +
    geom_point(size = 1.1, alpha = 0.55) +
    labs(title = paste0("Gráfico de traza log-verosimilitud: Sim ",s),
         x = "Iteración MCMC",
         y = "Log-verosimilitud") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  x11()
  print(ll_trace)
}

# ==============================================================================
#@                 INFERENCE ON THE NUMBER OF CLUSTERS
# ==============================================================================
#Posterior number of clusters


for (s in cadenas) {  
  
  #the current MCMC chain
  cadena <- get(paste0("cadena_",s))
  data <- get(paste0("simulated_data_",s))
  
  # True number of clusters:
  k_true <- data$k_true
  
  # Compute the number of clusters K at each iteration. 
  k_values <- apply(cadena$xi, 1, function(x) length(unique(x))) #note: The 1 applies for rows
  
  # Compute the posterior distribution of K
  k_table <- table(k_values) / length(k_values)
  
  # Plot
  x11()
  par(mfrow = c(1,1), mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
  plot(as.numeric(names(k_table)), k_table, type = "h", lwd = 3, col = "dodgerblue4",
       main = paste0("Distribución posterior del número de clusters: Sim", s),
       xlab = "Número de clusters", ylab = "Densidad",
       yaxt = "n")  # Suppress default y-axis
  
  axis(2, at = pretty(k_table), labels = pretty(k_table))
}




# ==============================================================================================
#@                             ADJUSTED RAND INDEX (ARI)
# ==============================================================================================

# Compute the ARI for all simulations (all iterations within a simulation)
ARI_mean <- numeric(length(cadenas))

for (s in cadenas) {
  
  current_chain <- get(paste0("cadena_", s))
  
  current_sampled <- get(paste0("simulated_data_", s))
  
  # Create xi_E_M_D: joint cluster assignment from the intersection of E, M, and D levels
  # Each unique combination of (xi_E, xi_M, xi_D) gets a consecutive integer label
  
  # Paste the three assignments into a single string key per observation
  joint_key <- paste(current_sampled$xi_E_true,
                     current_sampled$xi_M_true,
                     current_sampled$xi_D_true,
                     sep = "_")
  
  # Relabel unique combinations consecutively: 1, 2, 3, ...
  xi_E_M_D <- as.integer(factor(joint_key, levels = unique(joint_key)))
  
  # Store in the list
  current_sampled$xi_E_M_D <- xi_E_M_D
  
  
  # ----------------------------------------------------------------
  # Compute ARI per iteration comparing true joint clusters vs
  # sampled joint clusters from the chain
  # ----------------------------------------------------------------
  
  n_iter <- nrow(current_chain$xi)
  
  ARI <- numeric(n_iter)
  
  # Compute the ARI per iteration
  for (i in seq(n_iter)) {
    ARI[i] <- adjustedRandIndex(
      current_sampled$xi_E_M_D, 
      as.vector(current_chain$xi[i, ])
    )
  }
  # Mean of ARI iterations
  ARI_mean[s] <- mean(ARI)
  
  assign(
    paste0("simulated_data_", s),
    current_sampled
  )
}
ARI_mean
round(ARI_mean,2)


# ==============================================================================================
#@                                    INFERENCE
# ==============================================================================================

cadenas <- c(1:18)

# ===============================================
#@     FILTER THE CHAIN BY k_mode ITERATIONS
# ===============================================

for (s in cadenas) {
  
  cat("--- Filtering simulation", s, " --- \n")
  
  # The current MCMC chain
  cadena <- get(paste0("cadena_",s))
  
  # Name of the filtered chain
  
  filtered_chain_name <- paste0("cadena_", s, "_filtered")
  
  #------------------------
  # Compute the mode
  #------------------------
  
  # Compute the number of clusters K at each iteration. 
  k_values <- apply(cadena$xi, 1, function(x) length(unique(x))) #note: The 1 applies for rows
  
  # Compute the posterior distribution of K (proportion of each k)
  k_table <- table(k_values) / length(k_values)
  
  # Find the mode of k (most repeated k) 
  k_mode <- as.numeric(names(k_table)[which.max(k_table)])
  
  #--------------------------------------------------
  # Filter the chain by iterations which k = k_mode
  #--------------------------------------------------
  
  # 1) Execute the function, then 2) store the filtered chain as cadena_1_filtered
  # Function filter_cadena_by_k is stored in funciones.R
  
  assign(
    filtered_chain_name,
    filter_cadena_by_k(cadena, k_mode)
  )
  
  cat("--- Completed --- \n")
}



# ===============================================
#@          SOLVE LABEL SWITCHING
# ===============================================

#-----------------------------------------------------------------------------
# Label switching will be solved by finding the permutation of betas that minimise the MSE
#Find the correct beta vector for each cluster
#-----------------------------------------------------------------------------


for (s in 8:18) {
  
  cat("--- PERMUTATING SIMULATION", s, " --- \n")
  
  # Get the cadena_s_filtered
  filtered_chain_name    <- paste0("cadena_", s, "_filtered")      
  current_filtered_chain <- get(filtered_chain_name)
  
  # Get its corresponding simulated data
  simulated_data_name    <- paste0("simulated_data_",s)
  current_simulated_data <- get(simulated_data_name)
  
  result <- solve_label_switching_cpp(
    cadena_filtered = current_filtered_chain,                   # Filtered MCMC chain
    y               = current_simulated_data$y,                 # The real y data
    X_full          = current_simulated_data$X_full,            # The real X data
    k_mode          = current_filtered_chain$filter_info$k_mode # The mode of k clusters
  )
  
  # Assigns the ordered betas to the cadena_s_filtered
  #assign(filtered_chain_name, result)     
  
  # Pint execution
  cat("--- Permutation for simulation", s, "Completed --- \n")
  cat("beta_k_correct_order has been added to the chain \n")
  
  # Store the results
  results_to_save <- paste0("/",filtered_chain_name, ".RData")
  
  save(result, file = paste0(path_to_results, results_to_save))
}


#--------------------------------------
# LOAD THE CHAINS (ALREADY PERMUTATED)
#--------------------------------------

for (s in cadenas) {
  
  # Filtered chain name (permuted)
  results_file_name <- file.path(path_to_results, paste0("cadena_", s, "_filtered", ".RData"))
  #Load chain
  temp_env <- new.env(parent = emptyenv())
  load(file = results_file_name, envir = temp_env)
  
  # Get the chain loaded
  loaded_obj <- get(ls(temp_env)[1], envir = temp_env)
  
  # Assign it with the name "cadena_s_filtered"
  assign(paste0("cadena_", s, "_filtered"), loaded_obj)
  
  
  cat("Loaded:", results_file_name, "\n")
}



# ===============================================
# COMPUTE EACH CLUSTER BETAs MEAN AND QUANTILES
# ===============================================

#-----------------------------------------------------------------------------
# Compute the mean and quantiles of each beta (in beta_k) for every cluster, and
# also compute the mean and quantiles for the intercept
#-----------------------------------------------------------------------------

for (s in cadenas) {
  
  cat("--- COMPUTING BETAS MEAN AND QUANTILES FOR SIMULATION", s, " --- \n")
  
  # Get the cadena_s_filtered
  filtered_chain_name    <- paste0("cadena_", s, "_filtered")      
  current_filtered_chain <- get(filtered_chain_name)
  
  
  result <- compute_posterior_summaries(
    cadena_filtered      <- current_filtered_chain,
    k_mode               <- current_filtered_chain$filter_info$k_mode
  )
  
  # Assigns the ordered betas to the cadena_s_filtered
  assign(filtered_chain_name, result)     
  
  # Pint execution
  cat("--- Computation for simulation", s, "Completed --- \n")
  cat("cadena_filtered$beta_k_summary has been added to the chain \n")
  cat("cadena_filtered$beta_int_summary has been added to the chain \n")
}

# ==============================================================================================
#@                                EXTERNAL VALIDATION (PREDICTION)
# ==============================================================================================

#---------------------------------------------------------------------
# Using the testing samples, predictive capacity metrics will be 
#computed (MSE, MAE, R2)
#---------------------------------------------------------------------


for (s in cadenas) {
  
  cat("--- COMPUTING BETAS MEAN AND QUANTILES FOR SIMULATION", s, " --- \n")
  
  # Get the cadena_s_filtered
  filtered_chain_name    <- paste0("cadena_", s, "_filtered")      
  current_filtered_chain <- get(filtered_chain_name)
  
  # Get its corresponding testing data
  testing_data_name    <- paste0("test_data_",s)
  current_testing_data <- get(testing_data_name)
  
  result <- compute_test_metrics(
    cadena_filtered      <- current_filtered_chain,
    test_data            <- current_testing_data
  )
  
  # Name of a list to store results of prediction metrics.
  list_metrics_name <- paste0("externa_validation_results_", s)
  
  
  # Assigns the results to the cadena_s_filtered
  assign(list_metrics_name, result)     
  
  # Pint execution
  cat("--- Metrics computation for simulation", s, "Completed --- \n")
  cat("MSE, MAE, R2 Stored in: ",list_metrics_name, "\n")
  cat("cadena_filtered$beta_int_summary has been added to the chain \n")
}



# ===================================================================
#                         EXECUTE WAIC — MODEL 2
# ===================================================================

#---------------------------------------------------------------------
# Compute WAIC
#---------------------------------------------------------------------

for (s in cadenas) {
  
  cat("--- COMPUTING WAIC FOR SIMULATION", s, " --- \n")
  
  # Get the cadena_s_filtered
  filtered_chain_name    <- paste0("cadena_", s, "_filtered")      
  current_filtered_chain <- get(filtered_chain_name)
  
  # Get its corresponding simulated data
  simulated_data_name    <- paste0("simulated_data_",s)
  current_simulated_data <- get(simulated_data_name)
  
  # Retrieve data from the simulated data
  y      <- current_simulated_data$y
  X_full <- current_simulated_data$X_full
  obs_to_mun <- current_simulated_data$obs_to_mun
  mun_map <- current_simulated_data$mun_map
  
  # Compute WAIC
  result_waic <- compute_WAIC_model2(cadena_filtered = current_filtered_chain,
                                     y               = y,
                                     X_full          = X_full,
                                     obs_to_mun      = obs_to_mun,
                                     mun_map         = mun_map)
  
  # Store WAIC inside the filtered chain 
  current_filtered_chain$WAIC_summary <- result_waic
  
  # Assigns the waic to the cadena_s_filtered
  assign(filtered_chain_name, current_filtered_chain)
  
  # Pint execution
  cat("--- WAIC computation for simulation", s, "Completed --- \n")
  
}

















