#@ utils

# ===================================================================
#@                             PREPARATION
# ===================================================================

#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("corrplot")
#install.packages("tictoc")
#install.packages("coda")
#install.packages("posterior")

suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(dplyr))) 
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(corrplot)))
suppressMessages(suppressWarnings(library(tictoc)))
suppressMessages(suppressWarnings(library(coda)))
suppressMessages(suppressWarnings(library(posterior)))


# Load xi samplers. C++ functions

Rcpp::sourceCpp("path/samplers_M3_rcpp.cpp")



#============================================================================================


#@--------                    3. NON-PARAMETRIC MODEL: THREE DP               ---------------


#============================================================================================




# =====================================================================================

#@                                3.1  GIBBS

# =====================================================================================

#---------------------------------------------------------------------------------
# The following functions are the Full Conditional Distributions used in the GIBBS
#---------------------------------------------------------------------------------

# ---------------------------------------
# 1) Update E cluster assignments (xi_E)
# ---------------------------------------

# Update the individual level cluster asignments (xi_E)

# Replaced by its equivalent C++ function
sample_xi_E <- function(y, X_E_full, xi_E, beta_k_E, sigma2_E_base, alpha_E, mu_E,
                        Z_full, W_full, xi_M, xi_D, beta_k_M, beta_k_D,
                        mun_kappa, obs_to_mun, xTx_vec, sqrt_kappa, beta_int) {
  # Sample FCD of xi_i^E for the E-level Dirichlet Process
  #
  # INPUTS:
  # y                : numeric vector (n)       - response variable
  # X_E_full         : numeric matrix (n x p_E) - individual-level covariates
  # xi_E             : integer vector (n)       - current E-level cluster assignments
  # beta_k_E         : list of vectors          - E-level cluster coefficients
  # sigma2_E_base    : scalar                   - base measure variance (sigma^2_E)
  # alpha_E          : scalar                   - E-level DP concentration parameter
  # mu_E             : numeric vector (p_E)     - base measure mean
  # Z_full           : numeric matrix (n x p_M) - municipality-level covariates (flat)
  # W_full           : numeric matrix (n x p_D) - department-level covariates (flat)
  # xi_M             : integer vector (n)       - current M-level cluster assignments
  # xi_D             : integer vector (n)       - current D-level cluster assignments
  # beta_k_M         : list of vectors          - M-level cluster coefficients
  # beta_k_D         : list of vectors          - D-level cluster coefficients
  # mun_kappa        : numeric vector           - kappa^2_{j,q} per municipality (flat)
  # obs_to_mun       : integer vector (n)       - pre-computed municipality indices
  # xTx_vec          : numeric vector (n)       - pre-computed x_i^T x_i per observation
  # sqrt_kappa       : numeric vector           - pre-computed sqrt(kappa^2_{j,q})
  # beta_int         : scalar                   - global intercept beta
  #
  # OUTPUTS:
  # list with:
  #   xi_E    : updated E-level cluster assignments (relabelled 1,...,K)
  #   beta_k_E: updated list of E-level cluster coefficients
  
  n   <- length(y)
  p_E <- ncol(X_E_full)
  
  # Pre-compute constants for new cluster beta sampling
  sigma2_inv  <- 1 / sigma2_E_base
  sigma2_diag <- sigma2_inv * diag(p_E)
  
  if (n > 1000) cat("Starting sample_xi_E for", n, "observations...\n")
  
  for (i in seq_len(n)) {
    
    if (n > 1000 && i %% 1000 == 0) {
      cat("  Processed", i, "/", n, "observations (",
          round(100 * i / n, 1), "%)\n", sep = "")
    }
    
    # ----------------------------------------------------------------
    # Pre-computed values for observation i
    # ----------------------------------------------------------------
    mun_idx      <- obs_to_mun[i]               # Lookup
    kappa2_i     <- mun_kappa[mun_idx]          # Variance 
    sqrt_kappa_i <- sqrt_kappa[mun_idx]
    kappa2_inv   <- 1 / kappa2_i
    x_i          <- X_E_full[i, ]               # row of individual covariates
    xTx_i        <- xTx_vec[i] 
    
    # ----------------------------------------------------------------
    # Compute M-level and D-level contributions for observation i:
    # z_{j,q}^T beta_k_M  +  w_q^T beta_k_D
    # These are fixed while we sample xi_E[i]
    # ----------------------------------------------------------------
    # Use municipal covariates of i and its corresponding cluster beta_k_M vector
    contrib_M <- sum(Z_full[i, ] * beta_k_M[[xi_M[i]]])         # z_{j,q}^T beta_k_M  
    # Use departamental covariates of i and its corresponding cluster beta_k_D vector
    contrib_D <- sum(W_full[i, ] * beta_k_D[[xi_D[i]]])         # w_q^T beta_k_D
    
    # ----------------------------------------------------------------
    # Residual r_i^(E) = y_i - (beta + z^T beta_k_M + w^T beta_k_D)
    # Used for: (1) new cluster marginal likelihood
    #           (2) new beta_k_E posterior sampling
    # ----------------------------------------------------------------
    r_i_E <- y[i] - (beta_int + contrib_M + contrib_D)
    
    # ----------------------------------------------------------------
    # Cluster counts excluding observation i
    # ----------------------------------------------------------------
    xi_minus_i    <- xi_E[-i]                           # xi_E excluding observation i
    cluster_table <- table(xi_minus_i)
    unique_k      <- as.integer(names(cluster_table))
    counts_k      <- as.integer(cluster_table)          # Counts of observations per cluster
    K_new_E       <- length(unique_k)
    
    # ----------------------------------------------------------------
    # Log-probability vector: slots 1..K_new_E = existing, K_new_E+1 = new
    # ----------------------------------------------------------------
    log_probs <- numeric(K_new_E + 1)                   # Last slot = new cluster
    
    # --- Probabilities for EXISTING clusters ---
    
    # p(xi_i^E = k_E | -) propto n_{-i,k} * N(y_i | beta + x^T beta_k_E + contrib_M + contrib_D,  kappa2)
    for (k_idx in seq_len(K_new_E)) {
      cluster_id  <- unique_k[k_idx]                    # Extracts cluster ID
      n_k         <- counts_k[k_idx]                    # Exctracts n_{-i,K_E}
      # Compute the mean beta + x^T beta_k_E + z_{j,q}^T beta_k_M + w_q^T beta_k_D
      mean_k      <- beta_int + sum(x_i * beta_k_E[[cluster_id]]) + contrib_M + contrib_D
      
      # Compute p(xi_i^E = k_E | -) propto n_{-i,k} * N(y_i | beta + x^T beta_k_E 
      #                                               + contrib_M + contrib_D,  kappa2)
      log_probs[k_idx] <- log(n_k) +
        dnorm(y[i], mean = mean_k, sd = sqrt_kappa_i, log = TRUE)
    }
    
    # --- Probability for NEW cluster ---
    
    # p(xi_i^E = K+1 | -) propto alpha * N(r_i^E | x_i^T mu_E, kappa2 + sigma2_E_base * x_i^T x_i)
    var_new  <- kappa2_i + (sigma2_E_base * xTx_i)     # kappa2 + sigma2_E_base * x_i^T x_i
    mean_new <- sum(x_i * mu_E)                        # x_i^T mu_E  (evaluated on r_i^E)
    
    # compute the log of p(xi_i^E = K+1 | -) = log(alpha_E) + log N(r^E_i | mean_new, var_new)
    log_probs[K_new_E + 1] <- log(alpha_E) +
      dnorm(r_i_E, mean = mean_new, sd = sqrt(var_new), log = TRUE)
    
    # --- Sample new cluster assignment ---
    # Normalize probabilities (log-sum-exp trick for numerical stability)
    max_log   <- max(log_probs)
    probs     <- exp(log_probs - max_log)              # subtracts max value for stability
    probs     <- probs / sum(probs)
    sampled_k <- sample.int(K_new_E + 1, size = 1, prob = probs)
    
    # ----------------------------------------------------------------
    # Update assignment
    # ----------------------------------------------------------------
    if (sampled_k <= K_new_E) {
      
      # Assign to existing cluster (Use original cluster ID)
      xi_E[i] <- unique_k[sampled_k]
      
    } else {
      
      
      # ---- New cluster: sample beta_{K+1}^E from its posterior ----
      # Create new cluster
      
      # V = (kappa2_inv * x_i x_i^T + sigma2_inv * I)^{-1}
      # M = V * (kappa2_inv * r_i^E * x_i + sigma2_inv * mu_E)
      
      A <- kappa2_inv * tcrossprod(x_i) + sigma2_diag
      
      # OPTIMIZATION: Use Cholesky throughout instead of solve()
      # Cholesky decomposition: A = R^T R
      R <- chol(A)                                           # A = R^T R
      
      # Mean vector: M = A^{-1} b
      # where b = (1/kappa^2) t_i x_i + (1/sigma^2) mu
      # Note: notice that the sum of kappa2_inv does not appear because we are sampling
      # based on observation i, the only observation in that cluster)
      b       <- kappa2_inv * r_i_E * x_i + sigma2_inv * mu_E
      
      # Solve A M = b using backsolve: M = R^{-1} (R^{-T} b)
      M_post  <- backsolve(R, backsolve(R, b, transpose = TRUE))
      
      # Sample: beta_new ~ N(M, A^{-1}) = M + R^{-1} z
      z <- rnorm(p_E)
      
      beta_new <- as.numeric(M_post + backsolve(R, z))
      
      # Find next available cluster ID
      # Find smallest available cluster ID (fills gaps if any)
      used_ids      <- unique(xi_E)
      max_id        <- max(used_ids)
      # Find first gap or use max+1
      all_ids       <- seq_len(max_id + 1)
      available_ids <- setdiff(all_ids, used_ids)  # Looks for gapps, if exists uses min(gap)
      #If gap exists, uses the smallest missing label otherwise use max_id + 1
      new_id <- if (length(available_ids) > 0) min(available_ids) else (max_id + 1L)
      
      beta_k_E[[new_id]] <- beta_new
      xi_E[i]            <- new_id
    }
  }
  # ----------------------------------------------------------------
  # Final cleanup: relabel clusters consecutively 1, 2, ..., K
  # ----------------------------------------------------------------
  used_clusters <- sort(unique(xi_E))
  K_final       <- length(used_clusters)
  
  # Extract only used beta_k
  beta_k_E_clean <- vector("list", K_final)
  for (idx in seq_along(used_clusters)) {
    beta_k_E_clean[[idx]] <- beta_k_E[[used_clusters[idx]]]
  }
  # Relabel xi to be 1, 2, 3, ..., K
  xi_E_final <- as.integer(factor(xi_E,
                                  levels = used_clusters,
                                  labels = seq_along(used_clusters)))
  
  if (n > 1000) cat("Completed! Final E-level clusters:", K_final, "\n")
  
  return(list(xi_E    = xi_E_final,
              beta_k_E = beta_k_E_clean))
}

#result_xi_E <- sample_xi_E(y, X_E_full, xi_E, beta_k_E, sigma2_E_base, alpha_E, mu_E,
#                    Z_full, W_full, xi_M, xi_D, beta_k_M, beta_k_D,
#                    mun_kappa, obs_to_mun, xTx_vec, sqrt_kappa, beta_int)

#xi_E     <- result_xi_E$xi_E
#beta_k_E <- result_xi_E$beta_k_E

# ---------------------------------------
# 2) Update M cluster assignments (xi_M)
# ---------------------------------------

# Update the municipal level cluster assignments (xi_M)

# Replaced by its equivalent C++ function
sample_xi_M <- function(y, Z_full, xi_M, beta_k_M, sigma2_M_base, alpha_M, mu_M,
                        X_E_full, W_full, xi_E, xi_D, beta_k_E, beta_k_D,
                        mun_kappa, obs_to_mun, zTz_vec, sqrt_kappa, beta_int) {
  # Sample FCD of xi_i^E for the E-level Dirichlet Process
  #
  # INPUTS:
  # y           : numeric vector (n)       - response variable
  # Z_full         : numeric matrix (n x p_M) - Municipal-level covariates
  # xi_M        : integer vector (n)       - current M-level cluster assignments
  # beta_k_M    : list of vectors          - M-level cluster coefficients
  # sigma2_E_base    : scalar                   - base measure variance (sigma^2_M)
  # alpha_M     : scalar                   - M-level DP concentration parameter
  # mu_E        : numeric vector (p_E)     - base measure mean
  # X_E_full      : numeric matrix (n x p_E) - Individual-level covariates (flat)
  # W_full      : numeric matrix (n x p_D) - department-level covariates (flat)
  # xi_E        : integer vector (n)       - current E-level cluster assignments
  # xi_D        : integer vector (n)       - current D-level cluster assignments
  # beta_k_E    : list of vectors          - E-level cluster coefficients
  # beta_k_D    : list of vectors          - D-level cluster coefficients
  # mun_kappa   : numeric vector           - kappa^2_{j,q} per municipality (flat)
  # obs_to_mun  : integer vector (n)       - pre-computed municipality indices
  # xTx_vec     : numeric vector (n)       - pre-computed x_i^T x_i per observation
  # sqrt_kappa  : numeric vector           - pre-computed sqrt(kappa^2_{j,q})
  # beta_int    : scalar                   - global intercept beta
  #
  # OUTPUTS:
  # list with:
  #   xi_M    : updated M-level cluster assignments (relabelled 1,...,K_M)
  #   beta_k_M: updated list of M-level cluster coefficients
  
  n   <- length(y)
  p_M <- ncol(Z_full)
  
  # Pre-compute constants for new cluster beta sampling
  sigma2_inv  <- 1 / sigma2_M_base
  sigma2_diag <- sigma2_inv * diag(p_M)
  
  if (n > 1000) cat("Starting sample_xi_M for", n, "observations...\n")
  
  for (i in seq_len(n)) {
    
    if (n > 1000 && i %% 1000 == 0) {
      cat("  Processed", i, "/", n, "observations (",
          round(100 * i / n, 1), "%)\n", sep = "")
    }
    
    # ----------------------------------------------------------------
    # Pre-computed values for observation i
    # ----------------------------------------------------------------
    
    mun_idx      <- obs_to_mun[i]               # Lookup
    kappa2_i     <- mun_kappa[mun_idx]          # Variance 
    sqrt_kappa_i <- sqrt_kappa[mun_idx]
    kappa2_inv   <- 1 / kappa2_i
    z_i          <- Z_full[i, ]               # row of individual covariates
    zTz_i        <- zTz_vec[i] 
    
    # ----------------------------------------------------------------
    # Compute M-level and D-level contributions for observation i:
    # z_{j,q}^T beta_k_M  +  w_q^T beta_k_D
    # These are fixed while we sample xi_M[i]
    # ----------------------------------------------------------------
    # Use municipal covariates of i and its corresponding cluster beta_k_M vector
    contrib_E <- sum(X_E_full[i, ] * beta_k_E[[xi_E[i]]])         # x_{i,j,q}^T beta_k_E  
    # Use departamental covariates of i and its corresponding cluster beta_k_D vector
    contrib_D <- sum(W_full[i, ] * beta_k_D[[xi_D[i]]])         # w_q^T beta_k_D
    
    # ----------------------------------------------------------------
    # Residual r_i^(E) = y_i - (beta + x_i^T beta_k_E + w_q^T beta_k_D)
    # Used for: (1) new cluster marginal likelihood
    #           (2) new beta_k_M posterior sampling
    # ----------------------------------------------------------------
    r_i_M <- y[i] - (beta_int + contrib_E + contrib_D)           # Residual
    
    # ----------------------------------------------------------------
    # Cluster counts excluding observation i
    # ----------------------------------------------------------------
    xi_minus_i    <- xi_M[-i]                           # xi_M excluding observation i
    cluster_table <- table(xi_minus_i)
    unique_k      <- as.integer(names(cluster_table))
    counts_k      <- as.integer(cluster_table)          # Counts of observations per cluster
    K_new_M       <- length(unique_k)
    
    # ----------------------------------------------------------------
    # Log-probability vector: slots 1..K_new_M = existing, K_new_M+1 = new
    # ----------------------------------------------------------------
    log_probs <- numeric(K_new_M + 1)                   # Last slot = new cluster
    
    # --- Probabilities for EXISTING clusters ---
    
    # p(xi_i^M = k | -) propto n_{-i,k} * N(y_i | beta + contrib_E + z^T beta_k_M + contrib_D, kappa2)
    for (k_idx in seq_len(K_new_M)) {
      cluster_id  <- unique_k[k_idx]                    # Extracts cluster ID
      n_k         <- counts_k[k_idx]                    # Exctracts n_{-i,K_M}
      # Compute the mean beta + x^T beta_k_E + z_{j,q}^T beta_k_M + w_q^T beta_k_D
      mean_k      <- beta_int + contrib_E + sum(z_i * beta_k_M[[cluster_id]]) + contrib_D
      
      # Compute p(xi_i^E = k_E | -) propto n_{-i,k} * N(y_i | beta + contrib_E 
      #                                               + z^T beta_k_M + contrib_D, kappa2)
      log_probs[k_idx] <- log(n_k) +
        dnorm(y[i], mean = mean_k, sd = sqrt_kappa_i, log = TRUE)
    }
    
    # --- Probability for NEW cluster ---
    
    # p(xi_i^M = K+1 | -) propto alpha * N(r_i^M | z_i^T mu_M, kappa2 + sigma2_M * z_i^T z_i)
    var_new  <- kappa2_i + (sigma2_M_base * zTz_i)     # kappa2 + sigma2_M_base * x_i^T x_i
    mean_new <- sum(z_i * mu_M)                        # z_i^T mu_M  (evaluated on r_i^M)
    
    # compute the log of p(xi_i^E = K+1 | -) = log(alpha_M) + log N(r^M_i | mean_new, var_new)
    log_probs[K_new_M + 1] <- log(alpha_M) +
      dnorm(r_i_M, mean = mean_new, sd = sqrt(var_new), log = TRUE)
    
    # --- Sample new cluster assignment ---
    # Normalize probabilities (log-sum-exp trick for numerical stability)
    max_log   <- max(log_probs)
    probs     <- exp(log_probs - max_log)              # subtracts max value for stability
    probs     <- probs / sum(probs)
    sampled_k <- sample.int(K_new_M + 1, size = 1, prob = probs)
    
    # ----------------------------------------------------------------
    # Update assignment
    # ----------------------------------------------------------------
    if (sampled_k <= K_new_M) {
      
      # Assign to existing cluster (Use original cluster ID)
      xi_M[i] <- unique_k[sampled_k]
      
    } else {
      
      
      # ---- New cluster: sample beta_{K+1}^E from its posterior ----
      # Create new cluster
      
      # V = (kappa2_inv * z_i z_i^T + sigma2_inv * I)^{-1}
      # M = V * (kappa2_inv * r_i^M * z_i + sigma2_inv * mu_M)
      
      A <- kappa2_inv * tcrossprod(z_i) + sigma2_diag
      
      # OPTIMIZATION: Use Cholesky throughout instead of solve()
      # Cholesky decomposition: A = R^T R
      R <- chol(A)                                           # A = R^T R
      
      # Mean vector: M = A^{-1} b
      # where b = (1/kappa^2) t_i z_i + (1/sigma^2) mu_Z
      # Note: notice that the sum of kappa2_inv does not appear because we are sampling
      # based on observation i, the only observation in that cluster)
      b       <- kappa2_inv * r_i_M * z_i + sigma2_inv * mu_M
      
      # Solve A M = b using backsolve: M = R^{-1} (R^{-T} b)
      M_post  <- backsolve(R, backsolve(R, b, transpose = TRUE))
      
      # Sample: beta_new ~ N(M, A^{-1}) = M + R^{-1} z_rand
      z_rand <- rnorm(p_M)
      
      beta_new <- as.numeric(M_post + backsolve(R, z_rand))
      
      # Find next available cluster ID
      # Find smallest available cluster ID (fills gaps if any)
      used_ids      <- unique(xi_M)
      max_id        <- max(used_ids)
      # Find first gap or use max+1
      all_ids       <- seq_len(max_id + 1)
      available_ids <- setdiff(all_ids, used_ids)  # Looks for gapps, if exists uses min(gap)
      #If gap exists, uses the smallest missing label otherwise use max_id + 1
      new_id <- if (length(available_ids) > 0) min(available_ids) else (max_id + 1L)
      
      beta_k_M[[new_id]] <- beta_new
      xi_M[i]            <- new_id
    }
  }
  
  # ----------------------------------------------------------------
  # Final cleanup: relabel clusters consecutively 1, 2, ..., K
  # ----------------------------------------------------------------
  used_clusters <- sort(unique(xi_M))
  K_final       <- length(used_clusters)
  
  # Extract only used beta_k
  beta_k_M_clean <- vector("list", K_final)
  for (idx in seq_along(used_clusters)) {
    beta_k_M_clean[[idx]] <- beta_k_M[[used_clusters[idx]]]
  }
  # Relabel xi to be 1, 2, 3, ..., K
  xi_M_final <- as.integer(factor(xi_M,
                                  levels = used_clusters,
                                  labels = seq_along(used_clusters)))
  
  if (n > 1000) cat("Completed! Final M-level clusters:", K_final, "\n")
  
  return(list(xi_M    = xi_M_final,
              beta_k_M = beta_k_M_clean))
}


#result_xi_M <- sample_xi_M(y, Z_full, xi_M, beta_k_M, sigma2_M_base, alpha_M, mu_M,
#                    X_E_full, W_full, xi_E, xi_D, beta_k_E, beta_k_D,
#                    mun_kappa, obs_to_mun, zTz_vec, sqrt_kappa, beta_int)

#xi_M <- result_xi_M$xi_M
#beta_k_M <- result_xi_M$beta_k_M



# ---------------------------------------
# 3) Update D cluster assignments (xi_D)
# ---------------------------------------

# Update the municipal level cluster assignments (xi_D)

# Replaced by its equivalent C++ function
sample_xi_D <- function(y, W_full, xi_D, beta_k_D, sigma2_D_base, alpha_D, mu_D,
                        X_E_full, Z_full, xi_E, xi_M, beta_k_E, beta_k_M,
                        mun_kappa, obs_to_mun, wTw_vec, sqrt_kappa, beta_int) {
  # Sample FCD of xi_i^D for the D-level Dirichlet Process
  #
  # INPUTS:
  # y           : numeric vector (n)       - response variable
  # W_full      : numeric matrix (n x p_D) - department-level covariates (flat)
  # xi_D        : integer vector (n)       - current D-level cluster assignments
  # beta_k_D    : list of vectors          - D-level cluster coefficients
  # sigma2_D_base: scalar                  - base measure variance (sigma^2_D)
  # alpha_D     : scalar                   - D-level DP concentration parameter
  # mu_D        : numeric vector (p_D)     - base measure mean
  # X_E_full    : numeric matrix (n x p_E) - individual-level covariates (flat)
  # Z_full      : numeric matrix (n x p_M) - municipality-level covariates (flat)
  # xi_E        : integer vector (n)       - current E-level cluster assignments
  # xi_M        : integer vector (n)       - current M-level cluster assignments
  # beta_k_E    : list of vectors          - E-level cluster coefficients
  # beta_k_M    : list of vectors          - M-level cluster coefficients
  # mun_kappa   : numeric vector           - kappa^2_{j,q} per municipality (flat)
  # obs_to_mun  : integer vector (n)       - pre-computed municipality indices
  # wTw_vec     : numeric vector (n)       - pre-computed w_q^T w_q per observation
  # sqrt_kappa  : numeric vector           - pre-computed sqrt(kappa^2_{j,q})
  # beta_int    : scalar                   - global intercept beta
  #
  # OUTPUTS:
  # list with:
  #   xi_D    : updated D-level cluster assignments (relabelled 1,...,K_D)
  #   beta_k_D: updated list of D-level cluster coefficients
  
  n   <- length(y)
  p_D <- ncol(W_full)
  
  # Pre-compute constants for new cluster beta sampling
  sigma2_inv  <- 1 / sigma2_D_base
  sigma2_diag <- sigma2_inv * diag(p_D)
  
  if (n > 1000) cat("Starting sample_xi_D for", n, "observations...\n")
  
  for (i in seq_len(n)) {
    
    if (n > 1000 && i %% 1000 == 0) {
      cat("  Processed", i, "/", n, "observations (",
          round(100 * i / n, 1), "%)\n", sep = "")
    }
    
    # ----------------------------------------------------------------
    # Pre-computed values for observation i
    # ----------------------------------------------------------------
    
    mun_idx      <- obs_to_mun[i]               # Lookup
    kappa2_i     <- mun_kappa[mun_idx]          # Variance 
    sqrt_kappa_i <- sqrt_kappa[mun_idx]
    kappa2_inv   <- 1 / kappa2_i
    w_i          <- W_full[i, ]               # row of departamental-level covariates
    wTw_i        <- wTw_vec[i] 
    
    # ----------------------------------------------------------------
    # Compute E-level and M-level contributions for observation i:
    # x_i^T beta_k_E  +  z_{j,q}^T beta_k_M
    # These are fixed while we sample xi_D[i]
    # ----------------------------------------------------------------
    # Use individual covariates of i and its corresponding cluster beta_k_M vector
    contrib_E <- sum(X_E_full[i, ] * beta_k_E[[xi_E[i]]])        # x_{i,j,q}^T beta_k_E  
    # Use municipal covariates of i and its corresponding cluster beta_k_D vector
    contrib_M <- sum(Z_full[i, ] * beta_k_M[[xi_M[i]]])          # z_{j,q}^T beta_k_M
    
    # ----------------------------------------------------------------
    # Residual r_i^(D) = y_i - (beta + x_i^T beta_k_E + z_{j,q}^T beta_k_M)
    # Used for: (1) new cluster marginal likelihood
    #           (2) new beta_k_D posterior sampling
    # ----------------------------------------------------------------
    r_i_D <- y[i] - (beta_int + contrib_E + contrib_M)           # Residual
    
    # ----------------------------------------------------------------
    # Cluster counts excluding observation i
    # ----------------------------------------------------------------
    xi_minus_i    <- xi_D[-i]                           # xi_D excluding observation i
    cluster_table <- table(xi_minus_i)
    unique_k      <- as.integer(names(cluster_table))
    counts_k      <- as.integer(cluster_table)          # Counts of observations per cluster
    K_new_D       <- length(unique_k)
    
    # ----------------------------------------------------------------
    # Log-probability vector: slots 1..K_new_M = existing, K_new_M+1 = new
    # ----------------------------------------------------------------
    log_probs <- numeric(K_new_D + 1)                   # Last slot = new cluster
    
    # --- Probabilities for EXISTING clusters ---
    
    # p(xi_i^D = k | -) propto n_{-i,k} * N(y_i | beta + contrib_E + contrib_M + w^T beta_k_D, kappa2)
    for (k_idx in seq_len(K_new_D)) {
      cluster_id  <- unique_k[k_idx]                    # Extracts cluster ID
      n_k         <- counts_k[k_idx]                    # Exctracts n_{-i,K_D}
      # Compute the mean beta + x^T beta_k_E + z_{j,q}^T beta_k_M + w_q^T beta_k_D
      mean_k      <- beta_int + contrib_E + contrib_M + sum(w_i * beta_k_D[[cluster_id]])
      
      # Compute p(xi_i^D = k_D | -) propto n_{-i,k} * N(y_i | beta + contrib_E
      #                                               + contrib_M + w^T beta_k_D, kappa2)
      log_probs[k_idx] <- log(n_k) +
        dnorm(y[i], mean = mean_k, sd = sqrt_kappa_i, log = TRUE)
    }
    
    # --- Probability for NEW cluster ---
    
    # p(xi_i^D = K+1 | -) propto alpha * N(r_i^D | w_i^T mu_D, kappa2 + sigma2_D * w_i^T w_i)
    var_new  <- kappa2_i + (sigma2_D_base * wTw_i)     # kappa2 + sigma2_D * w_i^T w_i
    mean_new <- sum(w_i * mu_D)                        # w_i^T mu_D  (evaluated on r_i^D)
    
    # Compute the log of p(xi_i^D = K+1 | -) = log(alpha_D) + log N(r^D_i | mean_new, var_new)
    log_probs[K_new_D + 1] <- log(alpha_D) +
      dnorm(r_i_D, mean = mean_new, sd = sqrt(var_new), log = TRUE)
    
    # --- Sample new cluster assignment ---
    # Normalize probabilities (log-sum-exp trick for numerical stability)
    max_log   <- max(log_probs)
    probs     <- exp(log_probs - max_log)              # Subtracts max value for stability
    probs     <- probs / sum(probs)
    sampled_k <- sample.int(K_new_D + 1, size = 1, prob = probs)
    
    # ----------------------------------------------------------------
    # Update assignment
    # ----------------------------------------------------------------
    if (sampled_k <= K_new_D) {
      
      # Assign to existing cluster (Use original cluster ID)
      xi_D[i] <- unique_k[sampled_k]
      
    } else {
      
      
      # ---- New cluster: sample beta_{K+1}^E from its posterior ----
      # Create new cluster
      
      # V = (kappa2_inv * w_i w_i^T + sigma2_inv * I)^{-1}
      # M = V * (kappa2_inv * r_i^D * w_i + sigma2_inv * mu_D)
      
      A <- kappa2_inv * tcrossprod(w_i) + sigma2_diag
      
      # OPTIMIZATION: Use Cholesky throughout instead of solve()
      # Cholesky decomposition: A = R^T R
      R <- chol(A)                                           # A = R^T R
      
      # Mean vector: M = A^{-1} b
      # where b = (1/kappa^2) r_i^D * w_i + (1/sigma^2) mu_D
      # Note: notice that the sum of kappa2_inv does not appear because we are sampling
      # based on observation i, the only observation in that cluster
      b       <- kappa2_inv * r_i_D * w_i + sigma2_inv * mu_D
      
      # Solve A M = b using backsolve: M = R^{-1} (R^{-T} b)
      M_post  <- backsolve(R, backsolve(R, b, transpose = TRUE))
      
      # Sample: beta_new ~ N(M, A^{-1}) = M + R^{-1} z_rand
      z_rand <- rnorm(p_D)
      
      beta_new <- as.numeric(M_post + backsolve(R, z_rand))
      
      # Find next available cluster ID
      # Find smallest available cluster ID (fills gaps if any)
      used_ids      <- unique(xi_D)
      max_id        <- max(used_ids)
      # Find first gap or use max+1
      all_ids       <- seq_len(max_id + 1)
      available_ids <- setdiff(all_ids, used_ids)  # Looks for gaps, if exists uses min(gap)
      # If gap exists, uses the smallest missing label otherwise use max_id + 1
      new_id <- if (length(available_ids) > 0) min(available_ids) else (max_id + 1L)
      
      beta_k_D[[new_id]] <- beta_new
      xi_D[i]            <- new_id
    }
  }
  
  # ----------------------------------------------------------------
  # Final cleanup: relabel clusters consecutively 1, 2, ..., K
  # ----------------------------------------------------------------
  used_clusters <- sort(unique(xi_D))
  K_final       <- length(used_clusters)
  
  # Extract only used beta_k
  beta_k_D_clean <- vector("list", K_final)
  for (idx in seq_along(used_clusters)) {
    beta_k_D_clean[[idx]] <- beta_k_D[[used_clusters[idx]]]
  }
  # Relabel xi to be 1, 2, 3, ..., K
  xi_D_final <- as.integer(factor(xi_D,
                                  levels = used_clusters,
                                  labels = seq_along(used_clusters)))
  
  if (n > 1000) cat("Completed! Final D-level clusters:", K_final, "\n")
  
  return(list(xi_D    = xi_D_final,
              beta_k_D = beta_k_D_clean))
}


#result_xi_D <- sample_xi_D(y, W_full, xi_D, beta_k_D, sigma2_D_base, alpha_D, mu_D,
#                           X_E_full, Z_full, xi_E, xi_M, beta_k_E, beta_k_M,
#                           mun_kappa, obs_to_mun, wTw_vec, sqrt_kappa, beta_int)

#xi_D <- result_xi_D$xi_D
#beta_k_D <- result_xi_D$beta_k_D



# --------------------------------------------
# 4) Update individual level atoms (beta_k_E)
# --------------------------------------------

sample_beta_k_E <- function(y, X_E_full, xi_E, beta_k_E, beta_int, sigma2_E_base, mu_E,
                            Z_full, W_full, xi_M, xi_D, beta_k_M, beta_k_D,
                            kappa_inv, obs_to_mun) {
  # Update E-level cluster-specific coefficient vectors beta_k_E
  #
  # INPUTS:
  # y            : numeric vector (n)       - response variable
  # X_E_full     : numeric matrix (n x p_E) - individual-level covariates matrix
  # xi_E         : integer vector (n)       - current E-level cluster assignments
  # beta_k_E     : list of vectors          - current E-level cluster coefficients
  # beta_int     : scalar                   - global intercept beta
  # sigma2_E_base: scalar                   - base measure variance (sigma^2_E)
  # mu_E         : numeric vector (p_E)     - base measure mean
  # Z_full       : numeric matrix (n x p_M) - municipality-level covariates matrix
  # W_full       : numeric matrix (n x p_D) - department-level covariates matrix
  # xi_M         : integer vector (n)       - current M-level cluster assignments
  # xi_D         : integer vector (n)       - current D-level cluster assignments
  # beta_k_M     : list of vectors          - M-level cluster coefficients
  # beta_k_D     : list of vectors          - D-level cluster coefficients
  # kappa_inv    : numeric vector           - pre-computed 1/kappa^2_{j,q} per municipality
  # obs_to_mun   : integer vector (n)       - pre-computed municipality indices
  #
  # OUTPUTS:
  # beta_k_E     : list of updated E-level cluster coefficient vectors
  
  K_E <- length(beta_k_E)
  p_E <- ncol(X_E_full)
  
  # Pre-compute constants
  sigma2_inv <- 1 / sigma2_E_base
  sigma2_I   <- sigma2_inv * diag(p_E)
  sigma2_mu  <- sigma2_inv * mu_E
  
  for (k in seq_len(K_E)) {
    
    # Find observations in E-level cluster k
    idx_k <- which(xi_E == k)
    n_k   <- length(idx_k)          # n of observations in cluster k
    
    # Skip empty clusters (shouldn't happen)
    if (n_k == 0) next
    
    # Get covariates and response for cluster k
    X_k <- X_E_full[idx_k, , drop = FALSE]      #dim (n_k x p_E)
    y_k <- y[idx_k]                             #length = n_k
    
    # ----------------------------------------------------------------
    # Compute residual r_i^(E) = y_i - (beta + z_{j,q}^T beta_k_M + w_q^T beta_k_D)
    # for each observation i in cluster k
    # ----------------------------------------------------------------
    
    # Compute z_i ^T Beta_k_M_i
    
    # z_{j,q}^T beta_k_M. Subsets in cluster k(n_k x p_E)
    # sapply(xi_M[idx_k], function(m) beta_k_M[[m]]) Returns a matrix (p_m x n_k)
    contrib_M_k <- rowSums(Z_full[idx_k, , drop = FALSE] *                      #length = n_k
                             t(sapply(xi_M[idx_k], function(m) beta_k_M[[m]])))
    contrib_D_k <- rowSums(W_full[idx_k, , drop = FALSE] *   # w_q^T beta_k_D
                             t(sapply(xi_D[idx_k], function(d) beta_k_D[[d]])))
    
    # Compute the residuals
    r_k <- y_k - (beta_int + contrib_M_k + contrib_D_k)      # r_i^(E) length n_k
    
    # Get weights: 1/kappa^2_{j,q} for each observation in cluster k
    kappa_inv_k <- kappa_inv[obs_to_mun[idx_k]]              # Subsets for cluster k
    
    # ----------------------------------------------------------------
    # Precision matrix:
    # A = sum_{i in C_k^E} (1/kappa^2_{j,q}) x_i x_i^T + (1/sigma^2_E) I_{p_E}
    # ----------------------------------------------------------------
    
    # Weighted X^TX: sum_i (1/kappa^2_i) x_i x_i^T = X^T diag(weights) X
    X_weighted   <- X_k * kappa_inv_k  # Each row of X_k scaled by 1/kappa^2_i (Element-wise multiply each row by weight)
    XtX_weighted <- t(X_weighted) %*% X_k
    
    # Weighted X^T r: sum_i (1/kappa^2_i) x_i r_i^(E)
    Xtr_weighted <- t(X_k) %*% (kappa_inv_k * r_k)
    
    # Precision matrix: A = sum_i (1/kappa^2_i) x_i x_i^T + (1/sigma^2_E) I_{p_E}
    A <- XtX_weighted + sigma2_I
    
    # Optimisation: Use Cholesky decomposition instead of solve()
    # Cholesky decomposition: A = R^T R
    R <- chol(A)
    
    # ----------------------------------------------------------------
    # Mean vector:
    # M = V * b
    # where b = sum_i (1/kappa^2_i) x_i r_i^(E) + (1/sigma^2_E) mu_E
    # ----------------------------------------------------------------
    b <- Xtr_weighted + sigma2_mu
    
    # Solve A M = b using backsolve: M = R^{-1} (R^{-T} b)
    M <- backsolve(R, backsolve(R, b, transpose = TRUE))
    
    # Sample: beta_k_E ~ N(M, V) where V = A^{-1}
    # Using: beta_k_E = M + R^{-1} * z_rand where z_rand ~ N(0, I)
    z_rand <- rnorm(p_E)
    beta_k_E[[k]] <- as.numeric(M + backsolve(R, z_rand))
  }
  
  return(beta_k_E)
}


#beta_k_E <- sample_beta_k_E(y, X_E_full, xi_E, beta_k_E, beta_int, sigma2_E_base, mu_E,
#                            Z_full, W_full, xi_M, xi_D, beta_k_M, beta_k_D,
#                            kappa_inv, obs_to_mun)


# --------------------------------------------
# 5) Update municipal level atoms (beta_k_M)
# --------------------------------------------

sample_beta_k_M <- function(y, Z_full, xi_M, beta_k_M, beta_int, sigma2_M_base, mu_M,
                            X_E_full, W_full, xi_E, xi_D, beta_k_E, beta_k_D,
                            kappa_inv, obs_to_mun) {
  # Update M-level cluster-specific coefficient vectors beta_k_M
  #
  # INPUTS:
  # y            : numeric vector (n)       - response variable
  # Z_full       : numeric matrix (n x p_M) - municipality-level covariates matrix
  # xi_M         : integer vector (n)       - current M-level cluster assignments
  # beta_k_M     : list of vectors          - current M-level cluster coefficients
  # beta_int     : scalar                   - global intercept beta
  # sigma2_M_base: scalar                   - base measure variance (sigma^2_M)
  # mu_M         : numeric vector (p_M)     - base measure mean
  # X_E_full     : numeric matrix (n x p_E) - individual-level covariates matrix
  # W_full       : numeric matrix (n x p_D) - department-level covariates matrix
  # xi_E         : integer vector (n)       - current E-level cluster assignments
  # xi_D         : integer vector (n)       - current D-level cluster assignments
  # beta_k_E     : list of vectors          - E-level cluster coefficients
  # beta_k_D     : list of vectors          - D-level cluster coefficients
  # kappa_inv    : numeric vector           - pre-computed 1/kappa^2_{j,q} per municipality
  # obs_to_mun   : integer vector (n)       - pre-computed municipality indices
  #
  # OUTPUTS:
  # beta_k_M     : list of updated M-level cluster coefficient vectors
  
  K_M <- length(beta_k_M)
  p_M <- ncol(Z_full)
  
  # Pre-compute constants
  sigma2_inv <- 1 / sigma2_M_base
  sigma2_I   <- sigma2_inv * diag(p_M)
  sigma2_mu  <- sigma2_inv * mu_M
  
  for (k in seq_len(K_M)) {
    
    # Find observations in M-level cluster k
    idx_k <- which(xi_M == k)
    n_k   <- length(idx_k)          # n of observations in cluster k
    
    # Skip empty clusters (shouldn't happen)
    if (n_k == 0) next
    
    # Get covariates and response for cluster k
    Z_k <- Z_full[idx_k, , drop = FALSE]        # dim (n_k x p_M)
    y_k <- y[idx_k]                             # length = n_k
    
    # ----------------------------------------------------------------
    # Compute residual r_i^(M) = y_i - (beta + x_i^T beta_k_E + w_q^T beta_k_D)
    # for each observation i in cluster k
    # ----------------------------------------------------------------
    
    # Compute x_i^T beta_k_E_i
    # x_i^T beta_k_E. Subsets in cluster k (n_k x p_E)
    # sapply(xi_E[idx_k], function(e) beta_k_E[[e]]) Returns a matrix (p_E x n_k)
    contrib_E_k <- rowSums(X_E_full[idx_k, , drop = FALSE] *                    # length = n_k
                             t(sapply(xi_E[idx_k], function(e) beta_k_E[[e]])))
    contrib_D_k <- rowSums(W_full[idx_k, , drop = FALSE] *   # w_q^T beta_k_D
                             t(sapply(xi_D[idx_k], function(d) beta_k_D[[d]])))
    
    # Compute the residuals
    r_k <- y_k - (beta_int + contrib_E_k + contrib_D_k)      # r_i^(M) length n_k
    
    # Get weights: 1/kappa^2_{j,q} for each observation in cluster k
    kappa_inv_k <- kappa_inv[obs_to_mun[idx_k]]              # Subsets for cluster k
    
    # ----------------------------------------------------------------
    # Precision matrix:
    # A = sum_{i in C_k^M} (1/kappa^2_{j,q}) z_i z_i^T + (1/sigma^2_M) I_{p_M}
    # ----------------------------------------------------------------
    
    # Weighted Z^TZ: sum_i (1/kappa^2_i) z_i z_i^T = Z^T diag(weights) Z
    Z_weighted   <- Z_k * kappa_inv_k  # Each row of Z_k scaled by 1/kappa^2_i (Element-wise multiply each row by weight)
    ZtZ_weighted <- t(Z_weighted) %*% Z_k
    
    # Weighted Z^T r: sum_i (1/kappa^2_i) z_i r_i^(M)
    Ztr_weighted <- t(Z_k) %*% (kappa_inv_k * r_k)
    
    # Precision matrix: A = sum_i (1/kappa^2_i) z_i z_i^T + (1/sigma^2_M) I_{p_M}
    A <- ZtZ_weighted + sigma2_I
    
    # Optimisation: Use Cholesky decomposition instead of solve()
    # Cholesky decomposition: A = R^T R
    R <- chol(A)
    
    # ----------------------------------------------------------------
    # Mean vector:
    # M = V * b
    # where b = sum_i (1/kappa^2_i) z_i r_i^(M) + (1/sigma^2_M) mu_M
    # ----------------------------------------------------------------
    b <- Ztr_weighted + sigma2_mu
    
    # Solve A M = b using backsolve: M = R^{-1} (R^{-T} b)
    M <- backsolve(R, backsolve(R, b, transpose = TRUE))
    
    # Sample: beta_k_M ~ N(M, V) where V = A^{-1}
    # Using: beta_k_M = M + R^{-1} * z_rand where z_rand ~ N(0, I)
    z_rand <- rnorm(p_M)
    beta_k_M[[k]] <- as.numeric(M + backsolve(R, z_rand))
  }
  
  return(beta_k_M)
}


#beta_k_M <- sample_beta_k_M(y, Z_full, xi_M, beta_k_M, beta_int, sigma2_M_base, mu_M,
#                           X_E_full, W_full, xi_E, xi_D, beta_k_E, beta_k_D,
#                           kappa_inv, obs_to_mun)


# --------------------------------------------
# 6) Update departament level atoms (beta_k_D)
# --------------------------------------------


sample_beta_k_D <- function(y, W_full, xi_D, beta_k_D, beta_int, sigma2_D_base, mu_D,
                            X_E_full, Z_full, xi_E, xi_M, beta_k_E, beta_k_M,
                            kappa_inv, obs_to_mun) {
  # Update D-level cluster-specific coefficient vectors beta_k_D
  #
  # INPUTS:
  # y            : numeric vector (n)       - response variable
  # W_full       : numeric matrix (n x p_D) - department-level covariates matrix
  # xi_D         : integer vector (n)       - current D-level cluster assignments
  # beta_k_D     : list of vectors          - current D-level cluster coefficients
  # beta_int     : scalar                   - global intercept beta
  # sigma2_D_base: scalar                   - base measure variance (sigma^2_D)
  # mu_D         : numeric vector (p_D)     - base measure mean
  # X_E_full     : numeric matrix (n x p_E) - individual-level covariates matrix
  # Z_full       : numeric matrix (n x p_M) - municipality-level covariates matrix
  # xi_E         : integer vector (n)       - current E-level cluster assignments
  # xi_M         : integer vector (n)       - current M-level cluster assignments
  # beta_k_E     : list of vectors          - E-level cluster coefficients
  # beta_k_M     : list of vectors          - M-level cluster coefficients
  # kappa_inv    : numeric vector           - pre-computed 1/kappa^2_{j,q} per municipality
  # obs_to_mun   : integer vector (n)       - pre-computed municipality indices
  #
  # OUTPUTS:
  # beta_k_D     : list of updated D-level cluster coefficient vectors
  
  K_D <- length(beta_k_D)
  p_D <- ncol(W_full)
  
  # Pre-compute constants
  sigma2_inv <- 1 / sigma2_D_base
  sigma2_I   <- sigma2_inv * diag(p_D)
  sigma2_mu  <- sigma2_inv * mu_D
  
  for (k in seq_len(K_D)) {
    
    # Find observations in D-level cluster k
    idx_k <- which(xi_D == k)
    n_k   <- length(idx_k)          # n of observations in cluster k
    
    # Skip empty clusters (shouldn't happen)
    if (n_k == 0) next
    
    # Get covariates and response for cluster k
    W_k <- W_full[idx_k, , drop = FALSE]        # dim (n_k x p_D)
    y_k <- y[idx_k]                             # length = n_k
    
    # ----------------------------------------------------------------
    # Compute residual r_i^(D) = y_i - (beta + x_i^T beta_k_E + z_{j,q}^T beta_k_M)
    # for each observation i in cluster k
    # ----------------------------------------------------------------
    
    # Compute x_i^T beta_k_E_i
    # x_i^T beta_k_E. Subsets in cluster k (n_k x p_E)
    # sapply(xi_E[idx_k], function(e) beta_k_E[[e]]) Returns a matrix (p_E x n_k)
    contrib_E_k <- rowSums(X_E_full[idx_k, , drop = FALSE] *                    # length = n_k
                             t(sapply(xi_E[idx_k], function(e) beta_k_E[[e]])))
    contrib_M_k <- rowSums(Z_full[idx_k, , drop = FALSE] *   # z_{j,q}^T beta_k_M
                             t(sapply(xi_M[idx_k], function(m) beta_k_M[[m]])))
    
    # Compute the residuals
    r_k <- y_k - (beta_int + contrib_E_k + contrib_M_k)      # r_i^(D) length n_k
    
    # Get weights: 1/kappa^2_{j,q} for each observation in cluster k
    kappa_inv_k <- kappa_inv[obs_to_mun[idx_k]]              # Subsets for cluster k
    
    # ----------------------------------------------------------------
    # Precision matrix:
    # A = sum_{i in C_k^D} (1/kappa^2_{j,q}) w_q w_q^T + (1/sigma^2_D) I_{p_D}
    # ----------------------------------------------------------------
    
    # Weighted W^TW: sum_i (1/kappa^2_i) w_i w_i^T = W^T diag(weights) W
    W_weighted   <- W_k * kappa_inv_k  # Each row of W_k scaled by 1/kappa^2_i (Element-wise multiply each row by weight)
    WtW_weighted <- t(W_weighted) %*% W_k
    
    # Weighted W^T r: sum_i (1/kappa^2_i) w_i r_i^(D)
    Wtr_weighted <- t(W_k) %*% (kappa_inv_k * r_k)
    
    # Precision matrix: A = sum_i (1/kappa^2_i) w_i w_i^T + (1/sigma^2_D) I_{p_D}
    A <- WtW_weighted + sigma2_I
    
    # Optimisation: Use Cholesky decomposition instead of solve()
    # Cholesky decomposition: A = R^T R
    R <- chol(A)
    
    # ----------------------------------------------------------------
    # Mean vector:
    # M = V * b
    # where b = sum_i (1/kappa^2_i) w_i r_i^(D) + (1/sigma^2_D) mu_D
    # ----------------------------------------------------------------
    b <- Wtr_weighted + sigma2_mu
    
    # Solve A M = b using backsolve: M = R^{-1} (R^{-T} b)
    M <- backsolve(R, backsolve(R, b, transpose = TRUE))
    
    # Sample: beta_k_D ~ N(M, V) where V = A^{-1}
    # Using: beta_k_D = M + R^{-1} * z_rand where z_rand ~ N(0, I)
    z_rand <- rnorm(p_D)
    beta_k_D[[k]] <- as.numeric(M + backsolve(R, z_rand))
  }
  
  return(beta_k_D)
}


#beta_k_D <- sample_beta_k_D(y, W_full, xi_D, beta_k_D, beta_int, sigma2_D_base, mu_D,
#                            X_E_full, Z_full, xi_E, xi_M, beta_k_E, beta_k_M,
#                            kappa_inv, obs_to_mun)




# -------------------------------------------------------
# 7) Update concentration parameter (E-level) (alpha_E)
# -------------------------------------------------------

sample_alpha_E <- function(alpha_E, xi_E, a_alpha_E, b_alpha_E) {
  
  # Update E-level concentration parameter using Escobar-West method
  #
  # INPUTS:
  # alpha_E : current E-level concentration parameter
  # xi_E    : integer vector (n) - current E-level cluster assignments
  # a_E     : scalar - shape hyperparameter of prior alpha_E ~ Gamma(a_E, b_E)
  # b_E     : scalar - rate hyperparameter of prior alpha_E ~ Gamma(a_E, b_E)
  #
  # OUTPUTS:
  # alpha_E : updated E-level concentration parameter
  
  K_E <- length(unique(xi_E))   # Number of occupied E-level clusters
  n   <- length(xi_E)           # Total number of observations
  
  # Sample auxiliary variable eta_E ~ Beta(alpha_E + 1  ,  n)
  eta_E    <- rbeta(1, shape1 = alpha_E + 1, shape2 = n)
  
  # Compute mixture probability pi_eta_E
  pi_eta_E <- (a_alpha_E + K_E - 1) / (a_alpha_E + K_E - 1 + n*(b_alpha_E - log(eta_E)))
  
  # Sample alpha_E as mixture of two Gamma distributions
  if (runif(1) < pi_eta_E) {
    return(rgamma(1, shape = a_alpha_E + K_E, rate = b_alpha_E - log(eta_E)))
  } else {
    return(rgamma(1, shape = a_alpha_E + K_E - 1, rate = b_alpha_E - log(eta_E)))
  }
}

#alpha_E <- sample_alpha_E(alpha_E, xi_E, a_alpha_E, b_alpha_E)


# -------------------------------------------------------
# 8) Update concentration parameter (M-level) (alpha_M)
# -------------------------------------------------------

sample_alpha_M <- function(alpha_M, xi_M, a_alpha_M, b_alpha_M) {
  
  # Update M-level concentration parameter using Escobar-West method
  #
  # INPUTS:
  # alpha_M : current M-level concentration parameter
  # xi_M    : integer vector (n) - current M-level cluster assignments
  # a_M     : scalar - shape hyperparameter of prior alpha_M ~ Gamma(a_M, b_M)
  # b_M     : scalar - rate hyperparameter of prior alpha_M ~ Gamma(a_M, b_M)
  #
  # OUTPUTS:
  # alpha_M : updated M-level concentration parameter
  
  K_M <- length(unique(xi_M))   # Number of occupied M-level clusters
  n   <- length(xi_M)           # Total number of observations
  
  # Sample auxiliary variable eta_M ~ Beta(alpha_M + 1  ,  n)
  eta_M    <- rbeta(1, shape1 = alpha_M + 1, shape2 = n)
  
  # Compute mixture probability pi_eta_M
  pi_eta_M <- (a_alpha_M + K_M - 1) / (a_alpha_M + K_M - 1 + n*(b_alpha_M - log(eta_M)))
  
  # Sample alpha_M as mixture of two Gamma distributions
  if (runif(1) < pi_eta_M) {
    return(rgamma(1, shape = a_alpha_M + K_M, rate = b_alpha_M - log(eta_M)))
  } else {
    return(rgamma(1, shape = a_alpha_M + K_M - 1, rate = b_alpha_M - log(eta_M)))
  }
}

#alpha_M <- sample_alpha_M(alpha_M, xi_M, a_alpha_M, b_alpha_M)


# -------------------------------------------------------
# 9) Update concentration parameter (D-level) (alpha_D)
# -------------------------------------------------------

sample_alpha_D <- function(alpha_D, xi_D, a_alpha_D, b_alpha_D) {
  
  # Update D-level concentration parameter using Escobar-West method
  #
  # INPUTS:
  # alpha_D : current D-level concentration parameter
  # xi_D    : integer vector (n) - current D-level cluster assignments
  # a_D     : scalar - shape hyperparameter of prior alpha_D ~ Gamma(a_D, b_D)
  # b_D     : scalar - rate hyperparameter of prior alpha_D ~ Gamma(a_D, b_D)
  #
  # OUTPUTS:
  # alpha_D : updated D-level concentration parameter
  
  K_D <- length(unique(xi_D))   # Number of occupied D-level clusters
  n   <- length(xi_D)           # Total number of observations
  
  # Sample auxiliary variable eta_D ~ Beta(alpha_D + 1  ,  n)
  eta_D    <- rbeta(1, shape1 = alpha_D + 1, shape2 = n)
  
  # Compute mixture probability pi_eta_D
  pi_eta_D <- (a_alpha_D + K_D - 1) / (a_alpha_D + K_D - 1 + n*(b_alpha_D - log(eta_D)))
  
  # Sample alpha_D as mixture of two Gamma distributions
  if (runif(1) < pi_eta_D) {
    return(rgamma(1, shape = a_alpha_D + K_D, rate = b_alpha_D - log(eta_D)))
  } else {
    return(rgamma(1, shape = a_alpha_D + K_D - 1, rate = b_alpha_D - log(eta_D)))
  }
}

#alpha_D <- sample_alpha_D(alpha_D, xi_D, a_alpha_D, b_alpha_D)



# ---------------------------------------------------
# 10) Update E-level mean of base distribution (mu_E)
# ---------------------------------------------------

sample_mu_E <- function(beta_k_E, sigma2_E_base, eta_mu_E, nu2_mu_E_inv) {
  
  # Update mean (\mu_E) of E-level base distribution mu_E vector
  #
  # INPUTS:
  # beta_k_E        : coefficient vectors
  # sigma2_E_base   : variance of beta_k_M ~ N_p_E(\mu_E, \sigma_E^2 I)
  # eta_mu_E        : mean of mu_E ~ N_p_E(eta_mu_E, nu2_mu_E)
  # nu2_mu_E_inv    : inverse of variance of nu2_mu_E  in mu_E ~ N_p_E(eta_mu_E, nu2_mu_E)
  #
  # OUTPUTS:
  # mu_E            : updated means vector 
  
  K_E <- length(beta_k_E)               # Number of occupied E-level clusters
  
  # Stack all beta_k_E into matrix (K_E x p_E) and sum columns. Stack by rows
  beta_E_matrix <- do.call(rbind, beta_k_E)
  sum_beta      <- colSums(beta_E_matrix)
  
  # Compute constants
  sigma2_E_inv  <- 1 / sigma2_E_base
  
  # Compute V and M matrix
  v <- 1 / (K_E * sigma2_E_inv + nu2_mu_E_inv)
  # Posterior mean: M = v * (1/sigma_E² * sum beta_K_E + 1/nu_E² * eta_E)
  M <- v * (sigma2_E_inv * sum_beta + nu2_mu_E_inv * eta_mu_E)
  
  # Sample from N_p_E(M, v*I_p)
  # Since V = v*I_p is diagonal, we can sample each component independently:
  # mu_E = M + sqrt(v) * z where z ~ N(0, I_p)
  mu_vec_E <- M + sqrt(v) * rnorm(length(M))
  
  return(mu_vec_E)
}

#mu_E <- sample_mu_E(beta_k_E, sigma2_E_base, eta_mu_E, nu2_mu_E_inv)

# ---------------------------------------------------
# 11) Update M-level mean of base distribution (mu_M)
# ---------------------------------------------------

sample_mu_M <- function(beta_k_M, sigma2_M_base, eta_mu_M, nu2_mu_M_inv) {
  
  # Update mean (\mu_M) of M-level base distribution mu_M vector
  #
  # INPUTS:
  # beta_k_M        : coefficient vectors
  # sigma2_M_base   : variance of beta_k_M ~ N_p_M(\mu_M, \sigma_M^2 I)
  # eta_mu_M        : mean of mu_M ~ N_p_M(eta_mu_M, nu2_mu_M)
  # nu2_mu_M_inv    : inverse of variance of nu2_mu_M  in mu_M ~ N_p_M(eta_mu_M, nu2_mu_M)
  #
  # OUTPUTS:
  # mu_M            : updated means vector 
  
  K_M <- length(beta_k_M)               # Number of occupied M-level clusters
  
  # Stack all beta_k_M into matrix (K_M x p_M) and sum columns. Stack by rows
  beta_M_matrix <- do.call(rbind, beta_k_M)
  sum_beta      <- colSums(beta_M_matrix)
  
  # Compute constants
  sigma2_M_inv  <- 1 / sigma2_M_base
  
  # Compute V and M matrix
  v <- 1 / (K_M * sigma2_M_inv + nu2_mu_M_inv)
  # Posterior mean: M = v * (1/sigma_M² * sum beta_K_M + 1/nu_M² * eta_M)
  M <- v * (sigma2_M_inv * sum_beta + nu2_mu_M_inv * eta_mu_M)
  
  # Sample from N_p_M(M, v*I_p)
  # Since V = v*I_p_M is diagonal, we can sample each component independently:
  # mu_M = M + sqrt(v) * z where z ~ N(0, I_p)
  mu_vec_M <- M + sqrt(v) * rnorm(length(M))
  
  return(mu_vec_M)
}

#mu_M <- sample_mu_M(beta_k_M, sigma2_M_base, eta_mu_M, nu2_mu_M_inv)



# ---------------------------------------------------
# 12) Update D-level mean of base distribution (mu_D)
# ---------------------------------------------------

sample_mu_D <- function(beta_k_D, sigma2_D_base, eta_mu_D, nu2_mu_D_inv) {
  
  # Update mean (\mu_D) of M-level base distribution mu_D vector
  #
  # INPUTS:
  # beta_k_D        : coefficient vectors
  # sigma2_D_base   : variance of beta_k_D ~ N_p_D(\mu_D, \sigma_D^2 I)
  # eta_mu_D        : mean of mu_D ~ N_p_D(eta_mu_D, nu2_mu_D)
  # nu2_mu_D_inv    : inverse of variance of nu2_mu_D  in mu_D ~ N_p_D(eta_mu_D, nu2_mu_D)
  #
  # OUTPUTS:
  # mu_D            : updated means vector 
  
  K_D <- length(beta_k_D)               # Number of occupied D-level clusters
  
  # Stack all beta_k_D into matrix (K_D x p_D) and sum columns. Stack by rows
  beta_D_matrix <- do.call(rbind, beta_k_D)
  sum_beta      <- colSums(beta_D_matrix)
  
  # Compute constants
  sigma2_D_inv  <- 1 / sigma2_D_base
  
  # Compute V and M matrix
  v <- 1 / (K_D * sigma2_D_inv + nu2_mu_D_inv)
  # Posterior mean: M = v * (1/sigma_D² * sum beta_K_D + 1/nu_D² * eta_D)
  M <- v * (sigma2_D_inv * sum_beta + nu2_mu_D_inv * eta_mu_D)
  
  # Sample from N_p_D(M, v*I_p)
  # Since V = v*I_p_D is diagonal, we can sample each component independently:
  # mu_D = M + sqrt(v) * z where z ~ N(0, I_p)
  mu_vec_D <- M + sqrt(v) * rnorm(length(M))
  
  return(mu_vec_D)
}

#mu_D <- sample_mu_D(beta_k_D, sigma2_D_base, eta_mu_D, nu2_mu_D_inv)



# -------------------------------------------------------------
# 13) Update E-level variance of base distribution (sigma2_E_base)
# -------------------------------------------------------------

sample_sigma2_E_base <- function(beta_k_E, mu_E, p_E, nu_E, nu_gamma2_E) {
  # Update variance (\sigma_E^2) of the individual base distribution
  #
  # INPUTS:
  # beta_k_E         : list of K_E vectors - individual cluster-specific coefficients
  # mu_E             : numeric vector (p_E) - mean of base distribution
  # p_E              : integer - dimension of individual coefficient vectors
  # nu_E             : scalar - prior shape hyperparameter
  # nu_gamma2_E      : scalar - pre-computed nu_E * gamma2
  #
  # OUTPUTS:
  # sigma2_E_base    : scalar - updated variance of base distribution
  
  K_E <- length(beta_k_E)
  
  # Posterior shape: (K_E p_E + nu_E)/2
  shape_post <- (K_E * p_E + nu_E) / 2
  
  # Optimisation: Vectorized computation of sum of squared deviations
  # Method 1: Stack into matrix and compute efficiently
  beta_matrix <- do.call(rbind, beta_k_E)       # K_E x p_E matrix
  
  # Center: subtract mu from each row
  centered <- sweep(beta_matrix, 2, mu_E, "-")  # K_E x p_E matrix
  
  # Sum of (beta_k_E - mu_E)^T(beta_k_E - mu_E) = sum of all squared elements
  sum_sq_dev <- sum(centered^2)
  
  # Posterior rate: (nu_E gamma2_E + Σ(beta_k_E - mu_E)^T(beta_k_E - mu_E)) / 2
  rate_post <- (nu_gamma2_E + sum_sq_dev) / 2
  
  # Sample from Inverse Gamma
  # In R: rgamma gives Gamma, so IG(a,b) = 1/Gamma(a,b)
  sigma2 <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
  
  return(sigma2)
}


#sigma2_E_base <- sample_sigma2_E_base(beta_k_E, mu_E, p_E, nu_E, nu_gamma2_E)



# -------------------------------------------------------------
# 14) Update M-level variance of base distribution (sigma2_M_base)
# -------------------------------------------------------------

sample_sigma2_M_base <- function(beta_k_M, mu_M, p_M, nu_M, nu_gamma2_M) {
  # Update variance (\sigma_M^2) of the municipal base distribution
  #
  # INPUTS:
  # beta_k_M         : list of K_M vectors - municipal cluster-specific coefficients
  # mu_M             : numeric vector (p_M) - mean of base distribution
  # p_M              : integer - dimension of municipal coefficient vectors
  # nu_M             : scalar - prior shape hyperparameter
  # nu_gamma2_M      : scalar - pre-computed nu_M * gamma2_M
  #
  # OUTPUTS:
  # sigma2_M_base    : scalar - updated variance of base distribution
  
  K_M <- length(beta_k_M)
  
  # Posterior shape: (K_M p_M + nu_M)/2
  shape_post <- (K_M * p_M + nu_M) / 2
  
  # Optimisation: Vectorized computation of sum of squared deviations
  # Method 1: Stack into matrix and compute efficiently
  beta_matrix <- do.call(rbind, beta_k_M)       # K_M x p_M matrix
  
  # Center: subtract mu from each row
  centered <- sweep(beta_matrix, 2, mu_M, "-")  # K_M x p_M matrix
  
  # Sum of (beta_k_M - mu_M)^T(beta_k_M - mu_M) = sum of all squared elements
  sum_sq_dev <- sum(centered^2)
  
  # Posterior rate: (nu_M gamma_M + Σ(beta_k_M - mu_M)^T(beta_k_M - mu_M)) / 2
  rate_post <- (nu_gamma2_M + sum_sq_dev) / 2
  
  # Sample from Inverse Gamma
  # In R: rgamma gives Gamma, so IG(a,b) = 1/Gamma(a,b)
  sigma2 <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
  
  return(sigma2)
}


#sigma2_M_base <- sample_sigma2_M_base(beta_k_M, mu_M, p_M, nu_M, nu_gamma2_M)


# -------------------------------------------------------------
# 15) Update D-level variance of base distribution (sigma2_D_base)
# -------------------------------------------------------------

sample_sigma2_D_base <- function(beta_k_D, mu_D, p_D, nu_D, nu_gamma2_D) {
  # Update variance (\sigma_D^2) of the departamental base distribution
  #
  # INPUTS:
  # beta_k_D         : list of K_D vectors - departament cluster-specific coefficients
  # mu_D             : numeric vector (p_D) - mean of base distribution
  # p_D              : integer - dimension of departament coefficient vectors
  # nu_D             : scalar - prior shape hyperparameter
  # nu_gamma2_D      : scalar - pre-computed nu_D * gamma2_D
  #
  # OUTPUTS:
  # sigma2_D_base    : scalar - updated variance of base distribution
  
  K_D <- length(beta_k_D)
  
  # Posterior shape: (K_D p_D + nu_D)/2
  shape_post <- (K_D * p_D + nu_D) / 2
  
  # Optimisation: Vectorized computation of sum of squared deviations
  # Method 1: Stack into matrix and compute efficiently
  beta_matrix <- do.call(rbind, beta_k_D)       # K_D x p_D matrix
  
  # Center: subtract mu from each row
  centered <- sweep(beta_matrix, 2, mu_D, "-")  # K_D x p_D matrix
  
  # Sum of (beta_k_D - mu_D)^T(beta_k_D - mu_D) = sum of all squared elements
  sum_sq_dev <- sum(centered^2)
  
  # Posterior rate: (nu_D gamma2_D + Σ(beta_k_D - mu_D)^T(beta_k_D - mu_D)) / 2
  rate_post <- (nu_gamma2_D + sum_sq_dev) / 2
  
  # Sample from Inverse Gamma
  # In R: rgamma gives Gamma, so IG(a,b) = 1/Gamma(a,b)
  sigma2 <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
  
  return(sigma2)
}


#sigma2_D_base <- sample_sigma2_D_base(beta_k_D, mu_D, p_D, nu_D, nu_gamma2_D)




# ------------------------------
# 16) Update intercept \beta (beta_int)
# ------------------------------
sample_beta_int <- function(y, X_E_full, Z_full, W_full,
                            xi_E, xi_M, xi_D,
                            beta_k_E, beta_k_M, beta_k_D,
                            mu_beta, sigma2_beta,
                            kappa_inv, obs_to_mun) {
  
  # Update global intercept beta_int
  #
  # INPUTS:
  # y            : numeric vector (n)       - response variable
  # X_E_full     : numeric matrix (n x p_E) - individual-level covariates (flat)
  # Z_full       : numeric matrix (n x p_M) - municipality-level covariates (flat)
  # W_full       : numeric matrix (n x p_D) - department-level covariates (flat)
  # xi_E         : integer vector (n)       - current E-level cluster assignments
  # xi_M         : integer vector (n)       - current M-level cluster assignments
  # xi_D         : integer vector (n)       - current D-level cluster assignments
  # beta_k_E     : list of vectors          - E-level cluster coefficients
  # beta_k_M     : list of vectors          - M-level cluster coefficients
  # beta_k_D     : list of vectors          - D-level cluster coefficients
  # mu_beta      : scalar                   - prior mean for beta ~ N(mu_beta, sigma2_beta)
  # sigma2_beta  : scalar                   - prior variance for beta ~ N(mu_beta, sigma2_beta)
  # kappa_inv    : numeric vector           - pre-computed 1/kappa^2_{j,q} per municipality
  # obs_to_mun   : integer vector (n)       - municipality index for each observation
  #
  # OUTPUTS:
  # beta_int     : scalar - updated global intercept
  
  n <- length(y)
  
  # ----------------------------------------------------------------
  # Compute residual r_i^(beta) = y_i - (x_i^T beta_k_E + z_{j,q}^T beta_k_M + w_q^T beta_k_D)
  # for each observation i (beta is excluded since we are sampling it)
  # ----------------------------------------------------------------
  
  # E-level contribution: x_i^T beta_k_E for each observation i
  # Compute x_i^T beta_k_E_i
  # x_i^T beta_k_E. Subsets in cluster k (n_k x p_E)
  # sapply(xi_E[idx_k], function(e) beta_k_E[[e]]) Returns a matrix (p_E x n_k)
  contrib_E <- rowSums(X_E_full *
                         t(sapply(xi_E, function(e) beta_k_E[[e]])))  # length n
  
  # M-level contribution: z_{j,q}^T beta_k_M for each observation i
  contrib_M <- rowSums(Z_full *
                         t(sapply(xi_M, function(m) beta_k_M[[m]])))  # length n
  
  # D-level contribution: w_q^T beta_k_D for each observation i
  contrib_D <- rowSums(W_full *
                         t(sapply(xi_D, function(d) beta_k_D[[d]])))  # length n
  
  # Residuals: r_i^(beta) = y_i - (x_i^T beta_k_E + z_{j,q}^T beta_k_M + w_q^T beta_k_D)
  residuals <- y - (contrib_E + contrib_M + contrib_D)
  
  # ----------------------------------------------------------------
  # Get 1/kappa2_{j(i),q(i)} for each observation i
  # ----------------------------------------------------------------
  kappa_inv_i <- kappa_inv[obs_to_mun]    # length n
  
  # ----------------------------------------------------------------
  # Posterior variance:
  # V = (sum_{q,j} n_{j,q}/kappa^2_{j,q} + 1/sigma^2_beta)^{-1}
  #   = (sum_i 1/kappa^2_{j(i),q(i)} + 1/sigma^2_beta)^{-1}
  # ----------------------------------------------------------------
  precision_post <- sum(kappa_inv_i) + 1 / sigma2_beta   # scalar
  V              <- 1 / precision_post                   # scalar
  
  # ----------------------------------------------------------------
  # Posterior mean:
  # M = V * (sum_{q,j,i} r_i^(beta)/kappa^2_{j,q} + mu_beta/sigma^2_beta)
  # ----------------------------------------------------------------
  weighted_sum_r <- sum(residuals * kappa_inv_i)         # sum_i r_i/kappa^2_i
  M              <- V * ( weighted_sum_r + (mu_beta / sigma2_beta) )
  
  # Sample beta ~ N(M, V)
  beta_int <- rnorm(1, mean = M, sd = sqrt(V))
  
  return(beta_int)
}

#beta_int <- sample_beta_int(y, X_E_full, Z_full, W_full,
#                            xi_E, xi_M, xi_D,
#                            beta_k_E, beta_k_M, beta_k_D,
#                            mu_beta, sigma2_beta,
#                            kappa_inv, obs_to_mun)


# ------------------------------
# 17) Update sigma2_beta 
# ------------------------------


sample_sigma2_beta <- function(mu_beta, beta_int, nu_beta, gamma2_beta){
  
  # Update beta_int variance sigma2_beta
  #
  # INPUTS:
  # mu_beta       : scalar - prior mean for \beta (beta_int)
  # beta_int      : scalar - global intercept 
  # nu_beta       : shape of GI in sigma2_beta
  # gamma2_beta   : rate of GI in sigma2_beta
  #
  # OUTPUTS:
  # sigma2_beta   : Variance of global intercept
  
  a_sig2_beta <- 0.5 * (nu_beta + 1)
  b_sig2_beta <- 0.5 * ( (nu_beta * gamma2_beta) + (beta_int - mu_beta)^2 )
  sigma2_beta <- 1 / rgamma(1, shape = a_sig2_beta, rate = b_sig2_beta)
  
  return(sigma2_beta)
}

#sigma2_beta <- sample_sigma2_beta(mu_beta, beta_int, nu_beta, gamma2_beta)



# -------------------------------------------------------
# 18) Update kappa2_{j,q} (municipality-level variances)
# -------------------------------------------------------
# Replaced by its C++ equivalent
sample_kappa2_jq <- function(y, X_E_full, Z_full, W_full,
                             xi_E, xi_M, xi_D,
                             beta_k_E, beta_k_M, beta_k_D,
                             beta_int, mun_map, kappa2_q, nu_kappa) {
  # Update municipality-level variances kappa^2_{j,q}
  #
  # INPUTS:
  # y            : numeric vector (n)       - response variable
  # X_E_full     : numeric matrix (n x p_E) - individual-level covariates 
  # Z_full       : numeric matrix (n x p_M) - municipality-level covariates
  # W_full       : numeric matrix (n x p_D) - department-level covariates 
  # xi_E         : integer vector (n)       - current E-level cluster assignments
  # xi_M         : integer vector (n)       - current M-level cluster assignments
  # xi_D         : integer vector (n)       - current D-level cluster assignments
  # beta_k_E     : list of vectors          - E-level cluster coefficients
  # beta_k_M     : list of vectors          - M-level cluster coefficients
  # beta_k_D     : list of vectors          - D-level cluster coefficients
  # beta_int     : scalar                   - global intercept beta
  # mun_map      : list                     - municipality mapping (q, j, start, end, n)
  # mun_kappa_q  : numeric vector (m)       - department-level hyperparameters kappa^2_q
  # nu_kappa     : scalar                   - hyperparameter nu_kappa
  #
  # OUTPUTS:
  # kappa2_jq    : numeric vector (n_mun_total) - updated municipality variances
  
  n_mun_total <- length(mun_map)
  kappa2_jq   <- numeric(n_mun_total)            # list with municipal variances
  
  for (mun_idx in seq_len(n_mun_total)) {
    
    # Get municipality info
    q    <- mun_map[[mun_idx]]$q       # Department index
    j    <- mun_map[[mun_idx]]$j       # Municipality index within department
    n_jq <- mun_map[[mun_idx]]$n       # Number of observations in municipality j,q
    
    # Get observation indices for this municipality
    idx_jq <- mun_map[[mun_idx]]$start:mun_map[[mun_idx]]$end
    
    # Extract data for municipality j,q
    y_jq <- y[idx_jq]                                    # Responses (length n_jq)
    
    # ----------------------------------------------------------------
    # Compute full fitted values: vartheta_{i,j,q} = beta + x_i^T beta_k_E
    #                                               + z_{j,q}^T beta_k_M
    #                                               + w_q^T beta_k_D
    # ----------------------------------------------------------------
    
    # E-level contribution: x_i^T beta_k_E for each i in municipality j,q
    # do.call(rbind, ...) stacks beta vectors -> (n_jq x p_E) matrix
    
    # Matrix with each observation beta_k_E
    beta_E_mat <- do.call(rbind, beta_k_E[xi_E[idx_jq]])               # n_jq x p_E
    # Compute x_i^T beta_k_E for the municipality mun_idx
    contrib_E  <- rowSums(X_E_full[idx_jq, , drop=FALSE] * beta_E_mat) # length n_jq
    
    # M-level contribution: z_{j,q}^T beta_k_M for each i in municipality j,q
    beta_M_mat <- do.call(rbind, beta_k_M[xi_M[idx_jq]])               # n_jq x p_M
    contrib_M  <- rowSums(Z_full[idx_jq, , drop=FALSE] * beta_M_mat)   # length n_jq
    
    # D-level contribution: w_q^T beta_k_D for each i in municipality j,q
    beta_D_mat <- do.call(rbind, beta_k_D[xi_D[idx_jq]])               # n_jq x p_D
    contrib_D  <- rowSums(W_full[idx_jq, , drop=FALSE] * beta_D_mat)   # length n_jq
    
    # Full fitted values: vartheta_{i,j,q} 
    fitted_jq <- beta_int + contrib_E + contrib_M + contrib_D   # length n_jq
    
    # Residuals: (y_{j,q} - psi_{j,q})
    residuals_jq <- y_jq - fitted_jq                            # length n_jq
    
    # Sum of squared residuals: (y_{j,q} - psi_{j,q})^T (y_{j,q} - psi_{j,q})
    sum_sq_resid <- sum(residuals_jq^2)
    
    # ----------------------------------------------------------------
    # Posterior parameters for Inverse Gamma:
    # shape: (nu_kappa + n_{j,q}) / 2
    # rate : (nu_kappa * kappa2_q + sum_sq_resid) / 2
    # ----------------------------------------------------------------
    shape_post <- (nu_kappa + n_jq) / 2
    rate_post  <- (nu_kappa * kappa2_q[q] + sum_sq_resid) / 2
    
    # Sample from Inverse Gamma: IG(a, b) = 1/Gamma(a, rate=b)
    kappa2_jq[mun_idx] <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
  }
  
  return(kappa2_jq)
}


#kappa2_jq <- sample_kappa2_jq(y, X_E_full, Z_full, W_full,
#                             xi_E, xi_M, xi_D,
#                             beta_k_E, beta_k_M, beta_k_D,
#                             beta_int, mun_map, kappa2_q, nu_kappa)

# -----------------------------------------------
# 19) Update kappa2_q (department-level variance)
# -----------------------------------------------

# Replaced by its C++ equivalent 

sample_kappa2_q <- function(mun_kappa, dept_to_mun, nu_kappa,
                            alpha_kappa, beta_kappa, m) {
  # Update department-level variance hyperparameters kappa^2_q
  #
  # INPUTS:
  # mun_kappa   : numeric vector     - current municipality variances kappa^2_{j,q}
  # dept_to_mun : list (m)           - dept_to_mun[[q]] = flat municipality indices in dept q
  # nu_kappa    : scalar             - hyperparameter nu_kappa
  # alpha_kappa : scalar             - prior shape alpha_kappa
  # beta_kappa  : scalar             - prior rate beta_kappa
  # m           : integer            - number of departments
  #
  # OUTPUTS:
  # kappa2_q    : numeric vector (m) - updated department-level hyperparameters
  
  kappa2_q <- numeric(m)
  
  for (q in seq_len(m)) {
    
    # Get pre-computed municipality indices for department q
    mun_indices_q <- dept_to_mun[[q]]
    
    # Number of municipalities n_q in department q
    n_q <- length(mun_indices_q)
    
    # Sum of 1/kappa^2_{j,q} for all municipalities j in department q
    sum_inv_kappa <- sum(1 / mun_kappa[mun_indices_q])
    
    # ----------------------------------------------------------------
    # Posterior parameters for Gamma:
    # shape: (nu_kappa * n_q + alpha_kappa) / 2
    # rate : beta_kappa/2 + (nu_kappa/2) * sum_{j=1}^{n_q} 1/kappa^2_{j,q}
    # ----------------------------------------------------------------
    shape_post <- (nu_kappa * n_q + alpha_kappa) / 2
    rate_post  <- (beta_kappa / 2) + (nu_kappa / 2) * sum_inv_kappa
    
    # Sample from Gamma(shape, rate)
    kappa2_q[q] <- rgamma(1, shape = shape_post, rate = rate_post)
  }
  
  return(kappa2_q)
}

#kappa2_q <- sample_kappa2_q(mun_kappa, dept_to_mun, nu_kappa,
#                            alpha_kappa, beta_kappa, m)


# ----------------------------------------------------------

# Helper functions to read back xi and beta_k from files 

#----------------------------------------------------------

# Read all xi samples from file
# Returns: matrix (n_samples x n_observations)
read_xi_samples <- function(file_path) {
  xi_matrix <- read.table(file_path, header = FALSE)
  return(as.matrix(xi_matrix))
}

# Read all beta_k samples from file
# Returns: list of length n_samples, each element is a list of K vectors of length p
read_beta_k_samples <- function(file_path, p) {
  
  # Read all lines as raw strings (handles variable-length rows)
  lines     <- readLines(file_path)
  n_samples <- length(lines)
  
  beta_k_all <- vector("list", n_samples)
  
  for (i in seq_len(n_samples) ) {
    # Analise row i into numeric vector
    row <- as.numeric(strsplit(lines[i], " ")[[1]])
    
    K <- as.integer(row[1])          # First element is K
    betas <- row[2:(K * p + 1)]      # Remaining elements are the beta vectors
    
    # Reconstruct list of K vectors of length p
    beta_k_all[[i]] <- vector("list", K)
    for (k in seq_len(K) ) {
      start <- (k - 1) * p + 1
      end   <- k * p
      beta_k_all[[i]][[k]] <- betas[start:end]
    }
  }
  
  return(beta_k_all)
}





#===============================================================================

#                    3.2  GIBBS SAMPLER FUNCTION (Model 3)

#===============================================================================


#-------------------------------------------------------------------------------
# This function runs the GIBBS sampler using as input the functions created
# above as the Full Conditional Distributions of each parameter in Model 3.
# Three Dirichlet Processes: E-level (individual), M-level (municipal),
# D-level (departamental).
#
# Data objects and initial values are read from the global environment.
#-------------------------------------------------------------------------------

MCMC3_BNP <- function(B,
                      # --- E-level DP hyperparameters ---
                      a_alpha_E, b_alpha_E,
                      eta_mu_E, nu2_mu_E_inv,
                      nu_E, nu_gamma2_E,
                      # --- M-level DP hyperparameters ---
                      a_alpha_M, b_alpha_M,
                      eta_mu_M, nu2_mu_M_inv,
                      nu_M, nu_gamma2_M,
                      # --- D-level DP hyperparameters ---
                      a_alpha_D, b_alpha_D,
                      eta_mu_D, nu2_mu_D_inv,
                      nu_D, nu_gamma2_D,
                      # --- Global intercept hyperparameters ---
                      mu_beta, nu_beta, gamma2_beta,
                      # --- Variance hyperparameters ---
                      nu_kappa, alpha_kappa, beta_kappa) {
  
  # --- Sampling configuration ---
  
  # Configuración de burn-in y thinning
  burn_in    <- ceiling(B * 1.5)      # 150% iterations as burn-in
  thin       <- 10                      # Store every 10 iterations
  total_iter <- burn_in + (B * thin)   # Total iterations to execute
  n_samples  <- B                      # Number of samples to store
  
  #=====================     FILE CONNECTIONS     =====================
  
  # Define file paths (one per xi level, one per beta_k level = 6 files)
  xi_E_file     <- file.path(path_text_files, "xi_E_samples.txt")           # xi_E file
  xi_M_file     <- file.path(path_text_files, "xi_M_samples.txt")           # xi_M file
  xi_D_file     <- file.path(path_text_files, "xi_D_samples.txt")           # xi_D file
  beta_k_E_file <- file.path(path_text_files, "beta_k_E_samples.txt")       # beta_k_E file
  beta_k_M_file <- file.path(path_text_files, "beta_k_M_samples.txt")       # beta_k_M file
  beta_k_D_file <- file.path(path_text_files, "beta_k_D_samples.txt")       # beta_k_D file
  
  # Open connections (open = "w" truncates file if it exists)
  xi_E_con     <- file(xi_E_file,     open = "w")
  xi_M_con     <- file(xi_M_file,     open = "w")
  xi_D_con     <- file(xi_D_file,     open = "w")
  beta_k_E_con <- file(beta_k_E_file, open = "w")
  beta_k_M_con <- file(beta_k_M_file, open = "w")
  beta_k_D_con <- file(beta_k_D_file, open = "w")
  
  cat("Files opened:\n")
  cat("  xi_E:     ", xi_E_file,     "\n")
  cat("  xi_M:     ", xi_M_file,     "\n")
  cat("  xi_D:     ", xi_D_file,     "\n")
  cat("  beta_k_E: ", beta_k_E_file, "\n")
  cat("  beta_k_M: ", beta_k_M_file, "\n")
  cat("  beta_k_D: ", beta_k_D_file, "\n")
  
  #=====================     STORAGE STRUCTURES     =====================
  
  # Pre-allocated list storage structure. Only for final samples. Parameter samples will be stored
  # as list except from parameters xi and beta_k (for all 3 levels), which will be stored as text files.
  
  cadena <- list()
  
  # 1-3. xi_E, xi_M, xi_D ---> stored in text files
  
  # 4-6. beta_k_E, beta_k_M, beta_k_D ---> stored in text files
  
  # 7. alpha_E, alpha_M, alpha_D
  cadena$alpha_E <- rep(NA, n_samples)
  cadena$alpha_M <- rep(NA, n_samples)
  cadena$alpha_D <- rep(NA, n_samples)
  
  # 8. mu_E, mu_M, mu_D (matrices: rows = samples, cols = covariates)
  cadena$mu_E <- matrix(NA, nrow = n_samples, ncol = p_E)
  cadena$mu_M <- matrix(NA, nrow = n_samples, ncol = p_M)
  cadena$mu_D <- matrix(NA, nrow = n_samples, ncol = p_D)
  
  colnames(cadena$mu_E) <- X_E_names     # Assignate colnames
  colnames(cadena$mu_M) <- Z_M_names
  colnames(cadena$mu_D) <- W_D_names
  
  # 9. sigma2_E_base, sigma2_M_base, sigma2_D_base
  cadena$sigma2_E_base <- rep(NA, n_samples)
  cadena$sigma2_M_base <- rep(NA, n_samples)
  cadena$sigma2_D_base <- rep(NA, n_samples)
  
  # 10. beta_int
  cadena$beta_int <- rep(NA, n_samples)
  
  # 11. sigma2_beta
  cadena$sigma2_beta <- rep(NA, n_samples)
  
  # 12-13. kappa2_jq (municipality) and kappa2_q (department)
  cadena$kappa2_jq <- vector("list", length = m)
  cadena$kappa2_q  <- vector("list", length = m)
  
  # create nested storage lists for kappa: depto -> municipality 
  for (q in 1:m) {
    n_municipios_q <- length(Y[[q]])
    #municipality variance
    cadena$kappa2_jq[[q]] <- vector("list", length = n_municipios_q)
    #departamental variance
    cadena$kappa2_q[[q]]  <- rep(NA, n_samples)
    for (j in 1:n_municipios_q) {
      cadena$kappa2_jq[[q]][[j]] <- rep(NA, n_samples)
    }
  }
  
  # Sample counter for stored samples
  sample_count <- 0L
  
  #=====================     PARAMETER UPDATING     =====================
  
  # GIBBS ITERATIONS
  
  for (b in 1:total_iter) {       
    
    #----------------------------------------------------------
    #      NON-PARAMETRIC PARAMETERS (3 DP levels)
    #----------------------------------------------------------
    
    # 1. Update E-level cluster assignments (xi_E) and E-level atoms (beta_k_E)
    result_xi_E <- sample_xi_E_cpp(y, X_E_full, xi_E, beta_k_E, sigma2_E_base,
                                   alpha_E, mu_E, Z_full, W_full, xi_M, xi_D,
                                   beta_k_M, beta_k_D, mun_kappa, obs_to_mun,
                                   xTx_vec, sqrt_kappa, beta_int)
    xi_E     <- result_xi_E$xi_E          # xi_E
    beta_k_E <- result_xi_E$beta_k_E      # beta_k_E for new clusters
    
    # 2. Update M-level cluster assignments (xi_M) and M-level atoms (beta_k_M)
    result_xi_M <- sample_xi_M_cpp(y, Z_full, xi_M, beta_k_M, sigma2_M_base,
                                   alpha_M, mu_M, X_E_full, W_full, xi_E, xi_D,
                                   beta_k_E, beta_k_D, mun_kappa, obs_to_mun,
                                   zTz_vec, sqrt_kappa, beta_int)
    xi_M     <- result_xi_M$xi_M          # xi_M
    beta_k_M <- result_xi_M$beta_k_M      # beta_k_M for new clusters
    
    # 3. Update D-level cluster assignments (xi_D) and D-level atoms (beta_k_D)
    result_xi_D <- sample_xi_D_cpp(y, W_full, xi_D, beta_k_D, sigma2_D_base,
                                   alpha_D, mu_D, X_E_full, Z_full, xi_E, xi_M,
                                   beta_k_E, beta_k_M, mun_kappa, obs_to_mun,
                                   wTw_vec, sqrt_kappa, beta_int)
    xi_D     <- result_xi_D$xi_D          # xi_D
    beta_k_D <- result_xi_D$beta_k_D      # beta_k_D for new clusters
    
    # 4. Update E-level cluster coefficients (beta_k_E)
    
    kappa_inv <- 1 / mun_kappa           # Avoid repeated division in the loop
    
    beta_k_E  <- sample_beta_k_E(y, X_E_full, xi_E, beta_k_E, beta_int,
                                 sigma2_E_base, mu_E, Z_full, W_full,
                                 xi_M, xi_D, beta_k_M, beta_k_D,
                                 kappa_inv, obs_to_mun)
    
    # 5. Update M-level cluster coefficients (beta_k_M)
    beta_k_M  <- sample_beta_k_M(y, Z_full, xi_M, beta_k_M, beta_int,
                                 sigma2_M_base, mu_M, X_E_full, W_full,
                                 xi_E, xi_D, beta_k_E, beta_k_D,
                                 kappa_inv, obs_to_mun)
    
    # 6. Update D-level cluster coefficients (beta_k_D)
    beta_k_D  <- sample_beta_k_D(y, W_full, xi_D, beta_k_D, beta_int,
                                 sigma2_D_base, mu_D, X_E_full, Z_full,
                                 xi_E, xi_M, beta_k_E, beta_k_M,
                                 kappa_inv, obs_to_mun)
    
    # 7. Update E-level concentration parameter (alpha_E)
    alpha_E <- sample_alpha_E(alpha_E, xi_E, a_alpha_E, b_alpha_E)
    
    # 8. Update M-level concentration parameter (alpha_M)
    alpha_M <- sample_alpha_M(alpha_M, xi_M, a_alpha_M, b_alpha_M)
    
    # 9. Update D-level concentration parameter (alpha_D)
    alpha_D <- sample_alpha_D(alpha_D, xi_D, a_alpha_D, b_alpha_D)
    
    # 10. Update E-level base measure mean (mu_E)
    mu_E <- sample_mu_E(beta_k_E, sigma2_E_base, eta_mu_E, nu2_mu_E_inv)
    
    # 11. Update M-level base measure mean (mu_M)
    mu_M <- sample_mu_M(beta_k_M, sigma2_M_base, eta_mu_M, nu2_mu_M_inv)
    
    # 12. Update D-level base measure mean (mu_D)
    mu_D <- sample_mu_D(beta_k_D, sigma2_D_base, eta_mu_D, nu2_mu_D_inv)
    
    # 13. Update E-level base measure variance (sigma2_E_base)
    sigma2_E_base <- sample_sigma2_E_base(beta_k_E, mu_E, p_E, nu_E, nu_gamma2_E)
    
    # 14. Update M-level base measure variance (sigma2_M_base)
    sigma2_M_base <- sample_sigma2_M_base(beta_k_M, mu_M, p_M, nu_M, nu_gamma2_M)
    
    # 15. Update D-level base measure variance (sigma2_D_base)
    sigma2_D_base <- sample_sigma2_D_base(beta_k_D, mu_D, p_D, nu_D, nu_gamma2_D)
    
    #----------------------------------------------------------
    #            PARAMETRIC PARAMETERS
    #----------------------------------------------------------
    
    # 16. Update global intercept (beta_int)
    beta_int <- sample_beta_int(y, X_E_full, Z_full, W_full,
                                xi_E, xi_M, xi_D,
                                beta_k_E, beta_k_M, beta_k_D,
                                mu_beta, sigma2_beta,
                                kappa_inv, obs_to_mun)
    
    # 17. Update global intercept variance (sigma2_beta)
    sigma2_beta <- sample_sigma2_beta(mu_beta, beta_int, nu_beta, gamma2_beta)
    
    # 18. Update municipality-level variances (kappa2_{j,q} = mun_kappa)
    mun_kappa <- sample_kappa2_jq_cpp(y, X_E_full, Z_full, W_full,
                                  xi_E, xi_M, xi_D,
                                  beta_k_E, beta_k_M, beta_k_D,
                                  beta_int, mun_map, kappa2_q, nu_kappa)
    
    mun_kappa  <- pmax(mun_kappa, 1e-4)  # This prevents the vector to go to cero
    
    # Recompute derived quantities after updating mun_kappa
    # These are used in steps 1-6 of the next iteration
    kappa_inv  <- 1 / mun_kappa
    # These are used in sample_xi and sample_beta_k next iteration
    sqrt_kappa <- sqrt(mun_kappa)
    
    # 19. Update department-level variance hyperparameters (kappa2_q)
    kappa2_q <- sample_kappa2_q_cpp(mun_kappa, dept_to_mun, nu_kappa,
                                alpha_kappa, beta_kappa, m)
    
    #=====================     SAMPLE STORAGE     =====================
    
    # ------------------------------
    # Save results after burn-in & thinning
    # ------------------------------
    
    if (b > burn_in && (b - burn_in) %% thin == 0) {
      sample_count <- sample_count + 1L
      
      # 1. Store xi_E (one row: n integers separated by spaces)
      writeLines(paste(xi_E, collapse = " "), con = xi_E_con)
      
      # 2. Store xi_M
      writeLines(paste(xi_M, collapse = " "), con = xi_M_con)
      
      # 3. Store xi_D
      writeLines(paste(xi_D, collapse = " "), con = xi_D_con)
      
      # 4. Store beta_k_E
      # One row: K_E followed by all beta vectors concatenated
      # |K|  beta_k[[1]]   beta_k[[2]]   beta_k[[3]]
      # Example row: "3 0.5 1.2 ... 0.8 -0.3 ... 1.1 0.4 ..."
      K_E_current  <- length(beta_k_E) 
      beta_k_E_flat <- c(K_E_current, unlist(beta_k_E))      # First element K_E_current, then the vector
      writeLines(paste(beta_k_E_flat, collapse = " "), con = beta_k_E_con)
      
      # 5. Store beta_k_M
      K_M_current  <- length(beta_k_M)
      beta_k_M_flat <- c(K_M_current, unlist(beta_k_M))
      writeLines(paste(beta_k_M_flat, collapse = " "), con = beta_k_M_con)
      
      # 6. Store beta_k_D
      K_D_current  <- length(beta_k_D)
      beta_k_D_flat <- c(K_D_current, unlist(beta_k_D))
      writeLines(paste(beta_k_D_flat, collapse = " "), con = beta_k_D_con)
      
      # 7. Store alpha_E, alpha_M, alpha_D
      cadena$alpha_E[sample_count] <- alpha_E
      cadena$alpha_M[sample_count] <- alpha_M
      cadena$alpha_D[sample_count] <- alpha_D
      
      # 8. Store mu_E, mu_M, mu_D (as matrix with cols variables and rows samples)
      cadena$mu_E[sample_count, ] <- mu_E
      cadena$mu_M[sample_count, ] <- mu_M
      cadena$mu_D[sample_count, ] <- mu_D
      
      # 9. Store sigma2_E_base, sigma2_M_base, sigma2_D_base
      cadena$sigma2_E_base[sample_count] <- sigma2_E_base
      cadena$sigma2_M_base[sample_count] <- sigma2_M_base
      cadena$sigma2_D_base[sample_count] <- sigma2_D_base
      
      # 10. Store beta_int
      cadena$beta_int[sample_count] <- beta_int
      
      # 11. Store sigma2_beta
      cadena$sigma2_beta[sample_count] <- sigma2_beta
      
      # 12. Store kappa2_jq (municipality variances)
      for (i in seq_len(n_mun_total)) {
        q <- mun_map[[i]]$q
        j <- mun_map[[i]]$j
        cadena$kappa2_jq[[q]][[j]][sample_count] <- mun_kappa[i]
      }
      
      # 13. Store kappa2_q (department variances)
      for (q in seq_len(m)) {
        cadena$kappa2_q[[q]][sample_count] <- kappa2_q[q]
      }
    }
    
    # Progress indicator
    if (b %% 100 == 0) cat("Iteración", b, "de", total_iter, "completada\n")
  }
  
  #---------------------
  # CLOSE FILE CONNECTIONS
  #---------------------
  
  close(xi_E_con)
  close(xi_M_con)
  close(xi_D_con)
  close(beta_k_E_con)
  close(beta_k_M_con)
  close(beta_k_D_con)
  
  cat("\nFiles closed successfully!\n")
  cat("  xi_E samples stored:     ", sample_count, "rows x", length(xi_E), "columns\n")
  cat("  xi_M samples stored:     ", sample_count, "rows x", length(xi_M), "columns\n")
  cat("  xi_D samples stored:     ", sample_count, "rows x", length(xi_D), "columns\n")
  cat("  beta_k_E samples stored: ", sample_count, "rows (K_E varies per row)\n")
  cat("  beta_k_M samples stored: ", sample_count, "rows (K_M varies per row)\n")
  cat("  beta_k_D samples stored: ", sample_count, "rows (K_D varies per row)\n\n")
  
  #==================    ADDITIONAL INFO    ===================#
  
  cadena$info <- list(
    total_iterations = total_iter,
    burn_in          = burn_in,
    thin             = thin,
    samples_stored   = sample_count,
    samples_requested = B
  )
  
  return(cadena)
}



# =====================================================================================

#@                               3.3 INFERENCE

# =====================================================================================

# This block of functions are designed for model inference:
# - Filter the chain by k_mode
# - 

#===============================================

#           FILTER THE CHAIN BY k_mode

#===============================================

#-------------------------------------------------------------------------------
# This function filters the resulting chain in the GIBBS sampler to keep only 
#iterations that have as number of used clusters the mode of k (k_mode)
#-------------------------------------------------------------------------------

filter_cadena_by_k_3DP <- function(cadena, k_E_mode, k_M_mode, k_D_mode) {
  
  # -------------------------------------------------------
  # Step 1: Identify iterations where number of clusters == k_mode
  # Use both xi and beta_k and cross-check for consistency
  # -------------------------------------------------------
  
  n <- nrow(cadena$xi_E)
  
  # From xi matrix: unique clusters per row (iteration)
  k_E_from_xi      <- apply(cadena$xi_E, 1, function(row) length(unique(row))) #individual level
  k_M_from_xi      <- apply(cadena$xi_M, 1, function(row) length(unique(row))) #municipal level
  k_D_from_xi      <- apply(cadena$xi_D, 1, function(row) length(unique(row))) #departament level
  
  # From beta_k list: length of each sublist
  k_E_from_betak <- sapply(cadena$beta_k_E, length)
  k_M_from_betak <- sapply(cadena$beta_k_M, length)
  k_D_from_betak <- sapply(cadena$beta_k_D, length)
  
  # Cross-check: both should agree on k per iteration. 
  if (!all(k_E_from_xi == k_E_from_betak)) {
    # Gives a warning if elements don't agree on some k iteration
    warning("Mismatch between K from xi and K from beta_k. ",
            "Number of mismatched iterations: ",
            sum(k_E_from_xi != k_E_from_betak),
            ". Filtering will use k_from_xi as reference.")
  }
  
  # Index of the iterations to keep per level
  keep_idx_E <- which(k_E_from_xi == k_E_mode)
  keep_idx_M <- which(k_M_from_xi == k_M_mode)
  keep_idx_D <- which(k_D_from_xi == k_D_mode)
  
  # Index of iterations to keep (intersection of the tree vectors)
  keep_idx   <- Reduce(intersect, list(keep_idx_E, keep_idx_M, keep_idx_D))
  
  # Number of iterations to keep
  n_kept   <- length(keep_idx)
  
  cat("Filtering cadena by k_mode = (", k_E_mode, k_M_mode , k_D_mode, ")\n")
  cat("  Total iterations:  ", n, "\n")
  cat("  Iterations kept:   ", n_kept, "\n")
  cat("  Percentage kept:   ", round(100 * n_kept / n, 1), "%\n")
  
  # -------------------------------------------------------
  # Step 2: Build filtered cadena
  # -------------------------------------------------------
  
  # Creates a new chan to store
  cadena_filtered <- list()
  
  # 1. xi: keep rows in keep_idx
  cadena_filtered$xi_E <- cadena$xi_E[keep_idx, , drop = FALSE]
  cadena_filtered$xi_M <- cadena$xi_M[keep_idx, , drop = FALSE]
  cadena_filtered$xi_D <- cadena$xi_D[keep_idx, , drop = FALSE]
  
  # 2. beta_k: keep sublists in keep_idx
  cadena_filtered$beta_k_E <- cadena$beta_k_E[keep_idx]
  cadena_filtered$beta_k_M <- cadena$beta_k_M[keep_idx]
  cadena_filtered$beta_k_D <- cadena$beta_k_D[keep_idx]
  
  # 3. Scalar parameters (vectors of length B)
  cadena_filtered$alpha_E       <- cadena$alpha_E[keep_idx]
  cadena_filtered$alpha_M       <- cadena$alpha_M[keep_idx]
  cadena_filtered$alpha_D       <- cadena$alpha_D[keep_idx]
  
  cadena_filtered$sigma2_E_base <- cadena$sigma2_E_base[keep_idx]
  cadena_filtered$sigma2_M_base <- cadena$sigma2_M_base[keep_idx]
  cadena_filtered$sigma2_D_base <- cadena$sigma2_D_base[keep_idx]
  
  cadena_filtered$beta_int    <- cadena$beta_int[keep_idx]
  cadena_filtered$sigma2_beta <- cadena$sigma2_beta[keep_idx]
  
  # 4. mu_vec: matrix (B x p), keep rows
  cadena_filtered$mu_E <- cadena$mu_E[keep_idx, , drop = FALSE]
  cadena_filtered$mu_M <- cadena$mu_M[keep_idx, , drop = FALSE]
  cadena_filtered$mu_D <- cadena$mu_D[keep_idx, , drop = FALSE]
  
  # 5. kappa2_jq: nested list [dept][mun], each element is vector of length B
  cadena_filtered$kappa2_jq <- lapply(cadena$kappa2_jq, function(dept) {
    lapply(dept, function(mun_chain) mun_chain[keep_idx])
  })
  
  # 6. kappa2_q: list of length m, each element is vector of length B
  cadena_filtered$kappa2_q <- lapply(cadena$kappa2_q, function(q_chain) {
    q_chain[keep_idx]
  })
  
  # 7. Store filtering metadata
  cadena_filtered$filter_info <- list(
    k_E_mode          = k_E_mode,
    k_M_mode          = k_M_mode,
    k_D_mode          = k_D_mode,
    keep_idx_E        = keep_idx_E,
    keep_idx_M        = keep_idx_M,
    keep_idx_D        = keep_idx_D,
    kept_iterations   = keep_idx,
    b_kept            = n_kept,
    b_total           = n,
    pct_kept        = round(100 * n_kept / n, 1)
  )
  
  return(cadena_filtered)
}



#===============================================

#       SOLVE THE LABEL_SWITCHING PROBLEM

#===============================================

#-------------------------------------------------------------------------------
# This function solves the label_switching problem as described in Gelman (2013).
#It permutes the betas for each cluster over each iteration of the filtered chain 
#(by k_mode), and identifies the permutation that minimises the MSE. In that way,
#identifies the beta vector of each cluster
#-------------------------------------------------------------------------------


# Replaced by its C++ equivalent function

solve_label_switching_3DP <- function(cadena_filtered, y, X_E_full, Z_full, W_full,
                                      k_E_mode, k_M_mode, k_D_mode) {
  
  n_iter <- nrow(cadena_filtered$xi_E)
  
  # ----------------------------------------------------------------
  # Helper: generate permutations for a given k_mode
  # Exact same fallback rule as solve_label_switching
  # ----------------------------------------------------------------
  get_permutations <- function(k) {
    if (k <= 7) {
      permu <- gtools::permutations(n = k, r = k)
    } else {
      n_permu <- factorial(7)   # 5040
      permu   <- matrix(nrow = 0, ncol = k)
      while (nrow(permu) < n_permu) {
        needed <- n_permu - nrow(permu)
        batch  <- t(replicate(needed, sample.int(k)))
        permu  <- unique(rbind(permu, batch))
      }
      permu <- permu[seq_len(n_permu), ]
    }
    return(permu)
  }
  
  # Generate permutations for each level independently
  permu_E   <- get_permutations(k_E_mode)
  permu_M   <- get_permutations(k_M_mode)
  permu_D   <- get_permutations(k_D_mode)
  n_permu_E <- nrow(permu_E)
  n_permu_M <- nrow(permu_M)
  n_permu_D <- nrow(permu_D)
  
  cat("Solving label switching problem (3 DP - Sequential Greedy)\n")
  cat("  Filtered iterations:", n_iter,   "\n")
  cat("  k_E_mode:", k_E_mode, "| Permutations:", n_permu_E, "\n")
  cat("  k_M_mode:", k_M_mode, "| Permutations:", n_permu_M, "\n")
  cat("  k_D_mode:", k_D_mode, "| Permutations:", n_permu_D, "\n\n")
  
  # Storage for corrected beta_k (same structure as cadena_filtered$beta_k_E/M/D)
  beta_k_E_correct_order <- vector("list", n_iter)
  beta_k_M_correct_order <- vector("list", n_iter)
  beta_k_D_correct_order <- vector("list", n_iter)
  
  # Storage for corrected xi
  xi_E_correct_order <- matrix(NA_integer_, nrow = n_iter, ncol = ncol(cadena_filtered$xi_E))
  xi_M_correct_order <- matrix(NA_integer_, nrow = n_iter, ncol = ncol(cadena_filtered$xi_M))
  xi_D_correct_order <- matrix(NA_integer_, nrow = n_iter, ncol = ncol(cadena_filtered$xi_D))
  
  for (b in seq_len(n_iter)) {
    
    if (b %% 100 == 0) cat("  Processing iteration", b, "/", n_iter, "\n")
    
    # --------------------------------------------------
    # Extract iteration-specific quantities
    # --------------------------------------------------
    xi_E_b     <- cadena_filtered$xi_E[b, ]       # row b of length n
    xi_M_b     <- cadena_filtered$xi_M[b, ]       # row b of length n
    xi_D_b     <- cadena_filtered$xi_D[b, ]       # row b of length n
    beta_int_b <- cadena_filtered$beta_int[b]     # scalar for iteration b
    beta_k_E_b <- cadena_filtered$beta_k_E[[b]]   # list of k_E_mode vectors for iter b
    beta_k_M_b <- cadena_filtered$beta_k_M[[b]]   # list of k_M_mode vectors for iter b
    beta_k_D_b <- cadena_filtered$beta_k_D[[b]]   # list of k_D_mode vectors for iter b
    
    # Build beta matrices (k_mode x p): row k = beta_k_b[[k]]
    beta_E_matrix <- do.call(rbind, beta_k_E_b)   # k_E_mode x p_E
    beta_M_matrix <- do.call(rbind, beta_k_M_b)   # k_M_mode x p_M
    beta_D_matrix <- do.call(rbind, beta_k_D_b)   # k_D_mode x p_D
    
    # Pre-compute original (uncorrected) M and D contributions
    # These are fixed during Step 1, D is fixed during Step 2
    contrib_M_orig <- rowSums(Z_full * beta_M_matrix[xi_M_b, , drop = FALSE])  # length n
    contrib_D_orig <- rowSums(W_full * beta_D_matrix[xi_D_b, , drop = FALSE])  # length n
    
    # ==============================================================
    # STEP 1: Correct E-level
    # Fix M (original order) and D (original order)
    # Criterion: full fitted = beta + X*perm_E + Z*beta_M_orig + W*beta_D_orig
    # Equivalent to backfitting E against residual r_E = y - beta - M_orig - D_orig
    # ==============================================================
    
    mse_E <- numeric(n_permu_E)    # MSE for every permutation of E
    
    for (s in seq_len(n_permu_E)) {        # Retrieve permutation s
      sigma_E           <- permu_E[s, ]
      permuted_E_matrix <- beta_E_matrix[sigma_E, , drop = FALSE]       # apply the permutation to the matrix. dim = k_E_mode x p_E
      contrib_E_perm    <- rowSums(X_E_full * permuted_E_matrix[xi_E_b, , drop = FALSE])  # computes x_^t beta_k_E. length = n
      y_hat             <- beta_int_b + contrib_E_perm + contrib_M_orig + contrib_D_orig  # computes y_hat with M and D level fixed
      mse_E[s]          <- mean((y - y_hat)^2)        # Computes MSE
    }
    
    best_sigma_E                <- permu_E[which.min(mse_E), ]      # Keeps the permutation that minimises mse_E
    beta_k_E_correct_order[[b]] <- beta_k_E_b[best_sigma_E]         # Keeps its beta_k s
    
    # best_sigma_E maps NEW labels -> OLD labels, so we invert it using match()
    xi_E_correct_order[b, ] <- match(xi_E_b, best_sigma_E)
    
    # Compute corrected E contribution (used as fixed term in Step 2)
    beta_E_matrix_corrected <- beta_E_matrix[best_sigma_E, , drop = FALSE]   # Corrects beta matrix with best permutation
    # compute x_^t beta_k_E with the best beta_k_E permutation
    contrib_E_corrected     <- rowSums(X_E_full * beta_E_matrix_corrected[xi_E_b, , drop = FALSE]) 
    
    # ==============================================================
    # STEP 2: Correct M-level
    # Fix E (corrected) and D (original order)
    # Criterion: full fitted = beta + X*beta_E_corrected + Z*perm_M + W*beta_D_orig
    # Equivalent to backfitting M against residual r_M = y - beta - E_corrected - D_orig
    # ==============================================================
    
    # analog to the procedure for mse_E
    mse_M <- numeric(n_permu_M)
    
    for (s in seq_len(n_permu_M)) {
      sigma_M           <- permu_M[s, ]
      permuted_M_matrix <- beta_M_matrix[sigma_M, , drop = FALSE]              # k_M_mode x p_M
      contrib_M_perm    <- rowSums(Z_full * permuted_M_matrix[xi_M_b, , drop = FALSE])
      y_hat             <- beta_int_b + contrib_E_corrected + contrib_M_perm + contrib_D_orig
      mse_M[s]          <- mean((y - y_hat)^2)
    }
    
    best_sigma_M                <- permu_M[which.min(mse_M), ]
    beta_k_M_correct_order[[b]] <- beta_k_M_b[best_sigma_M]
    
    # relabel xi_M consistently with the selected permutation
    xi_M_correct_order[b, ] <- match(xi_M_b, best_sigma_M)
    
    # Compute corrected M contribution (used as fixed term in Step 3)
    beta_M_matrix_corrected <- beta_M_matrix[best_sigma_M, , drop = FALSE]
    contrib_M_corrected     <- rowSums(Z_full * beta_M_matrix_corrected[xi_M_b, , drop = FALSE])
    
    # ==============================================================
    # STEP 3: Correct D-level
    # Fix E (corrected) and M (corrected)
    # Criterion: full fitted = beta + X*beta_E_corrected + Z*beta_M_corrected + W*perm_D
    # Equivalent to backfitting D against residual r_D = y - beta - E_corrected - M_corrected
    # ==============================================================
    
    # analog to the procedure for mse_E
    mse_D <- numeric(n_permu_D)
    
    for (s in seq_len(n_permu_D)) {
      sigma_D           <- permu_D[s, ]
      permuted_D_matrix <- beta_D_matrix[sigma_D, , drop = FALSE]              # k_D_mode x p_D
      contrib_D_perm    <- rowSums(W_full * permuted_D_matrix[xi_D_b, , drop = FALSE])
      y_hat             <- beta_int_b + contrib_E_corrected + contrib_M_corrected + contrib_D_perm
      mse_D[s]          <- mean((y - y_hat)^2)
    }
    
    best_sigma_D                <- permu_D[which.min(mse_D), ]
    beta_k_D_correct_order[[b]] <- beta_k_D_b[best_sigma_D]
    
    # relabel xi_D consistently with the selected permutation
    xi_D_correct_order[b, ] <- match(xi_D_b, best_sigma_D)
    
  }
  
  # Add corrected beta_k to cadena_filtered (same pattern as solve_label_switching)
  cadena_filtered$beta_k_E_correct_order <- beta_k_E_correct_order
  cadena_filtered$beta_k_M_correct_order <- beta_k_M_correct_order
  cadena_filtered$beta_k_D_correct_order <- beta_k_D_correct_order
  
  
  # Add corrected xi to cadena_filtered
  cadena_filtered$xi_E_correct_order <- xi_E_correct_order
  cadena_filtered$xi_M_correct_order <- xi_M_correct_order
  cadena_filtered$xi_D_correct_order <- xi_D_correct_order
  
  cat("\ncadena_filtered$beta_k_E_correct_order successfully added.\n")
  cat("cadena_filtered$beta_k_M_correct_order successfully added.\n")
  cat("cadena_filtered$beta_k_D_correct_order successfully added.\n")
  
  cat("cadena_filtered$xi_E_correct_order successfully added.\n")
  cat("cadena_filtered$xi_M_correct_order successfully added.\n")
  cat("cadena_filtered$xi_D_correct_order successfully added.\n")
  
  return(cadena_filtered)
}


#===============================================

#       COMPUTE POSTERIOR SUMMARIES

#===============================================


compute_posterior_summaries_3DP <- function(cadena_filtered, k_E_mode, k_M_mode, k_D_mode) {
  
  n_iter <- length(cadena_filtered$beta_k_E_correct_order)
  p_E    <- length(cadena_filtered$beta_k_E_correct_order[[1]][[1]])
  p_M    <- length(cadena_filtered$beta_k_M_correct_order[[1]][[1]])
  p_D    <- length(cadena_filtered$beta_k_D_correct_order[[1]][[1]])
  
  cat("Computing posterior summaries\n")
  cat("  Iterations:", n_iter,   "\n")
  cat("  k_E_mode:  ", k_E_mode, "\n")
  cat("  k_M_mode:  ", k_M_mode, "\n")
  cat("  k_D_mode:  ", k_D_mode, "\n")
  cat("  p_E:       ", p_E,      "\n")
  cat("  p_M:       ", p_M,      "\n")
  cat("  p_D:       ", p_D,      "\n\n")
  
  # -------------------------------------------------
  # beta_k_E summaries (E-level, individual)
  # -------------------------------------------------
  
  # For each E-level cluster k, build a (n_iter x p_E) matrix where
  # row b = cadena_filtered$beta_k_E_correct_order[[b]][[k]]
  # Then compute quantiles and mean column-wise (across iterations for each beta)
  
  beta_k_E_summary <- vector("list", k_E_mode) # creates a list, size: number of E-level clusters
  
  # Iterate for E-level clusters
  for (k in seq_len(k_E_mode)) {
    
    # Stack iteration b, cluster k into matrix: (n_iter x p_E)
    beta_k_E_matrix <- do.call(rbind, lapply(cadena_filtered$beta_k_E_correct_order,
                                             function(iter_b) iter_b[[k]]))
    # dim(beta_k_E_matrix) = n_iter x p_E
    # Each column j contains the n_iter draws for coefficient j of E-level cluster k
    
    # Compute column-wise summaries (each col = one beta coefficient across iterations)
    # Computes: Mean, q025, q50, q975 (mean is the estimator and quantiles are used for
    # confidence intervals)
    beta_k_E_summary[[k]] <- list(
      mean = colMeans(beta_k_E_matrix),
      q025 = apply(beta_k_E_matrix, 2, quantile, probs = 0.025),  # Note: 2 for columns
      q50  = apply(beta_k_E_matrix, 2, quantile, probs = 0.5),
      q975 = apply(beta_k_E_matrix, 2, quantile, probs = 0.975)
    )
  }
  
  # -------------------------------------------------
  # beta_k_M summaries (M-level, municipal)
  # -------------------------------------------------
  
  # For each M-level cluster k, build a (n_iter x p_M) matrix where
  # row b = cadena_filtered$beta_k_M_correct_order[[b]][[k]]
  # Then compute quantiles and mean column-wise (across iterations for each beta)
  
  beta_k_M_summary <- vector("list", k_M_mode) # creates a list, size: number of M-level clusters
  
  # Iterate for M-level clusters
  for (k in seq_len(k_M_mode)) {
    
    # Stack iteration b, cluster k into matrix: (n_iter x p_M)
    beta_k_M_matrix <- do.call(rbind, lapply(cadena_filtered$beta_k_M_correct_order,
                                             function(iter_b) iter_b[[k]]))
    # dim(beta_k_M_matrix) = n_iter x p_M
    # Each column j contains the n_iter draws for coefficient j of M-level cluster k
    
    # Compute column-wise summaries (each col = one beta coefficient across iterations)
    # Computes: Mean, q025, q50, q975 (mean is the estimator and quantiles are used for
    # confidence intervals)
    beta_k_M_summary[[k]] <- list(
      mean = colMeans(beta_k_M_matrix),
      q025 = apply(beta_k_M_matrix, 2, quantile, probs = 0.025),  # Note: 2 for columns
      q50  = apply(beta_k_M_matrix, 2, quantile, probs = 0.5),
      q975 = apply(beta_k_M_matrix, 2, quantile, probs = 0.975)
    )
  }
  
  # -------------------------------------------------
  # beta_k_D summaries (D-level, departamental)
  # -------------------------------------------------
  
  # For each D-level cluster k, build a (n_iter x p_D) matrix where
  # row b = cadena_filtered$beta_k_D_correct_order[[b]][[k]]
  # Then compute quantiles and mean column-wise (across iterations for each beta)
  
  beta_k_D_summary <- vector("list", k_D_mode) # creates a list, size: number of D-level clusters
  
  # Iterate for D-level clusters
  for (k in seq_len(k_D_mode)) {
    
    # Stack iteration b, cluster k into matrix: (n_iter x p_D)
    beta_k_D_matrix <- do.call(rbind, lapply(cadena_filtered$beta_k_D_correct_order,
                                             function(iter_b) iter_b[[k]]))
    # dim(beta_k_D_matrix) = n_iter x p_D
    # Each column j contains the n_iter draws for coefficient j of D-level cluster k
    
    # Compute column-wise summaries (each col = one beta coefficient across iterations)
    # Computes: Mean, q025, q50, q975 (mean is the estimator and quantiles are used for
    # confidence intervals)
    beta_k_D_summary[[k]] <- list(
      mean = colMeans(beta_k_D_matrix),
      q025 = apply(beta_k_D_matrix, 2, quantile, probs = 0.025),  # Note: 2 for columns
      q50  = apply(beta_k_D_matrix, 2, quantile, probs = 0.5),
      q975 = apply(beta_k_D_matrix, 2, quantile, probs = 0.975)
    )
  }
  
  # ------------------------------------------------
  # beta_int summaries
  # ------------------------------------------------
  # Computes the mean and quantiles for the beta intercept
  # cadena_filtered$beta_int is a vector of length n_iter
  
  beta_int_summary <- list(
    mean = mean(cadena_filtered$beta_int),
    q025 = quantile(cadena_filtered$beta_int, probs = 0.025),
    q50  = quantile(cadena_filtered$beta_int, probs = 0.500),
    q975 = quantile(cadena_filtered$beta_int, probs = 0.975)
  )
  
  # -------------------------------------------------------
  # Store in cadena_filtered
  # -------------------------------------------------------
  
  # Store the results in the same chain (list)
  cadena_filtered$beta_k_E_summary  <- beta_k_E_summary
  cadena_filtered$beta_k_M_summary  <- beta_k_M_summary
  cadena_filtered$beta_k_D_summary  <- beta_k_D_summary
  cadena_filtered$beta_int_summary  <- beta_int_summary
  
  cat("--- Beta_int,Beta_E, Beta_M, Beta_D coefficients mean and quantiles computed ---\n")
  
  return(cadena_filtered)
}


# =====================================================================================

#@                          3.4  MODEL VALIDATION

# =====================================================================================


#===============================================

# EXTERNAL VALIDATION: MSE, MAE, R2, coverage

#===============================================

#-------------------------------------------------------------------------
# This function computes prediction testing metrics: MSE, MAE, and R2.
# It assigns to each observation a cluster randomly sampled considering the 
# proportion of observations per cluster in each iteration of the GIBBS 
# sampler. Adapted for the three-DP model (E, M, D levels).
#-------------------------------------------------------------------------

compute_test_metrics_3DP <- function(cadena_filtered, test_data) {
  
  k_E_mode <- cadena_filtered$filter_info$k_E_mode
  k_M_mode <- cadena_filtered$filter_info$k_M_mode
  k_D_mode <- cadena_filtered$filter_info$k_D_mode
  n_iter   <- nrow(cadena_filtered$xi_E_correct_order)
  n_test   <- length(test_data$y)
  
  # Fixed posterior mean estimates for beta_intercept (same across all iterations)
  beta_int_hat <- cadena_filtered$beta_int_summary$mean
  
  # Build (k_E_mode x p_E) matrix of posterior mean betas: row k = beta_k_E_summary[[k]]$mean
  # nrow = k_E, ncol = p_E (individual-level covariates)
  beta_k_E_hat <- do.call(rbind, lapply(seq_len(k_E_mode), function(k) {
    cadena_filtered$beta_k_E_summary[[k]]$mean
  }))
  # dim(beta_k_E_hat) = k_E_mode x p_E
  
  # Build (k_M_mode x p_M) matrix of posterior mean betas: row k = beta_k_M_summary[[k]]$mean
  # nrow = k_M, ncol = p_M (municipal-level covariates)
  beta_k_M_hat <- do.call(rbind, lapply(seq_len(k_M_mode), function(k) {
    cadena_filtered$beta_k_M_summary[[k]]$mean
  }))
  # dim(beta_k_M_hat) = k_M_mode x p_M
  
  # Build (k_D_mode x p_D) matrix of posterior mean betas: row k = beta_k_D_summary[[k]]$mean
  # nrow = k_D, ncol = p_D (departamental-level covariates)
  beta_k_D_hat <- do.call(rbind, lapply(seq_len(k_D_mode), function(k) {
    cadena_filtered$beta_k_D_summary[[k]]$mean
  }))
  # dim(beta_k_D_hat) = k_D_mode x p_D
  
  cat("Computing test predictions\n")
  cat("  Test observations:", n_test,   "\n")
  cat("  Iterations:       ", n_iter,   "\n")
  cat("  k_E_mode:         ", k_E_mode, "\n")
  cat("  k_M_mode:         ", k_M_mode, "\n")
  cat("  k_D_mode:         ", k_D_mode, "\n\n")
  
  # -------------------------------------------------------
  # Step 1: For each iteration b, compute y_hat for all
  # test observations using sampled cluster per level
  # -------------------------------------------------------
  
  # Matrix to store predictions: (n_iter x n_test)
  y_hat_matrix <- matrix(NA, nrow = n_iter, ncol = n_test)
  
  for (b in seq_len(n_iter)) { # Iterates over Gibbs iterations in filtered chain
    
    if (b %% 100 == 0) cat("  Processing iteration", b, "/", n_iter, "\n")
    
    # ---- E-level cluster sampling ----
    # Extract cluster assignments for iteration b from the corrected xi_E matrix
    xi_E_b <- cadena_filtered$xi_E_correct_order[b, ]
    # Compute the proportions of each E-level cluster k in row b from matrix xi_E
    pi_E_b <- as.numeric(table(factor(xi_E_b, levels = 1:k_E_mode))) / length(xi_E_b)
    # pi_E_b is a vector of length k_E_mode: pi_E_b[k] = proportion in cluster k
    
    # For each test observation i, sample an E-level cluster using pi_E_b as probabilities
    # size = n_test samples 1 cluster per observation in 1:k_E_mode with probability pi_E_b
    sampled_k_E <- sample.int(k_E_mode, size = n_test, replace = TRUE, prob = pi_E_b)
    # sampled_k_E[i] = E-level cluster assigned to observation i in iteration b
    
    # ---- M-level cluster sampling ----
    # Extract cluster assignments for iteration b from the corrected xi_M matrix
    xi_M_b <- cadena_filtered$xi_M_correct_order[b, ]
    # Compute the proportions of each M-level cluster k in row b from matrix xi_M
    pi_M_b <- as.numeric(table(factor(xi_M_b, levels = 1:k_M_mode))) / length(xi_M_b)
    # pi_M_b is a vector of length k_M_mode: pi_M_b[k] = proportion in cluster k
    
    # For each test observation i, sample an M-level cluster using pi_M_b as probabilities
    sampled_k_M <- sample.int(k_M_mode, size = n_test, replace = TRUE, prob = pi_M_b)
    # sampled_k_M[i] = M-level cluster assigned to observation i in iteration b
    
    # ---- D-level cluster sampling ----
    # Extract cluster assignments for iteration b from the corrected xi_D matrix
    xi_D_b <- cadena_filtered$xi_D_correct_order[b, ]
    # Compute the proportions of each D-level cluster k in row b from matrix xi_D
    pi_D_b <- as.numeric(table(factor(xi_D_b, levels = 1:k_D_mode))) / length(xi_D_b)
    # pi_D_b is a vector of length k_D_mode: pi_D_b[k] = proportion in cluster k
    
    # For each test observation i, sample a D-level cluster using pi_D_b as probabilities
    sampled_k_D <- sample.int(k_D_mode, size = n_test, replace = TRUE, prob = pi_D_b)
    # sampled_k_D[i] = D-level cluster assigned to observation i in iteration b
    
    # ---- Assign corresponding beta_k to each test observation ----
    
    # E-level: (n_test x p_E) - each observation gets the beta_k_E of its sampled E-cluster
    beta_E_per_obs <- beta_k_E_hat[sampled_k_E, , drop = FALSE]
    
    # M-level: (n_test x p_M) - each observation gets the beta_k_M of its sampled M-cluster
    beta_M_per_obs <- beta_k_M_hat[sampled_k_M, , drop = FALSE]
    
    # D-level: (n_test x p_D) - each observation gets the beta_k_D of its sampled D-cluster
    beta_D_per_obs <- beta_k_D_hat[sampled_k_D, , drop = FALSE]
    
    # Linear predictor: vartheta_i = beta_int + x_i^T beta_k_E + z_i^T beta_k_M + w_i^T beta_k_D
    # matrix nrow = n_iter, ncol = n_test
    y_hat_matrix[b, ] <- beta_int_hat +
      rowSums(test_data$X_E_full * beta_E_per_obs) +  # x_i^T beta_k_E
      rowSums(test_data$Z_full   * beta_M_per_obs) +  # z_i^T beta_k_M
      rowSums(test_data$W_full   * beta_D_per_obs)    # w_i^T beta_k_D
  }
  
  # -------------------------------------------------------
  # Step 2: Summarize across iterations for each observation
  # -------------------------------------------------------
  # Compute q0.025, mean and q0.975
  
  # Prediction matrix: 4 rows x n_test cols
  # Row 1 = mean, Row 2 = q0.025, Row 3 = median (q50), Row 4 = q0.975
  
  # The 2 applies by columns
  pred_matrix <- apply(y_hat_matrix, 2, function(obs_draws) {
    c(mean = mean(obs_draws),
      q025 = quantile(obs_draws, 0.025),
      q50  = quantile(obs_draws, 0.50),
      q975 = quantile(obs_draws, 0.975))
  })
  # Assign rownames
  rownames(pred_matrix) <- c("mean", "q025", "q50", "q975")
  # dim(pred_matrix) = 4 x n_test
  
  # -------------------------------------------------------
  # Step 3: Coverage — does true y fall in [q025, q975]?
  # -------------------------------------------------------
  
  y_test <- test_data$y    # observed y value (length n_test)
  # Assess whether y real value falls in the Confidence Interval
  # (TRUE if so, FALSE otherwise)
  in_CI    <- (y_test >= pred_matrix["q025", ]) & (y_test <= pred_matrix["q975", ])
  # Coverage is the percentage of observations whose real y falls in CI
  coverage <- mean(in_CI) * 100
  
  cat("\n  The % y_test in 95% CI:", round(coverage, 2), "%\n")
  
  # -------------------------------------------------------
  # Step 4: Compute MSE, MAE, R2 using posterior mean as point estimate
  # -------------------------------------------------------
  
  # Isolate the y_hat predictions (means)
  y_hat_mean <- pred_matrix["mean", ]
  # y observed - y predicted
  residuals  <- y_test - y_hat_mean
  
  # Compute metrics
  mse <- mean(residuals^2)
  mae <- mean(abs(residuals))
  r2  <- 1 - sum(residuals^2) / sum((y_test - mean(y_test))^2)
  
  cat("  MSE: ", round(mse, 4), "\n")
  cat("  RMSE: ", round(sqrt(mse), 4), "\n")
  cat("  MAE: ", round(mae, 4), "\n")
  cat("  R2:  ", round(r2,  4), "\n")
  
  
  # -----------------
  # Return results
  # -----------------
  
  return(list(
    pred_matrix = pred_matrix,    # 4 x n_test matrix (mean, q025, q50, q975)
    in_CI       = in_CI,          # logical vector length n_test
    coverage    = coverage,       # % of y_test inside CI
    mse         = mse,
    rmse        = sqrt(mse),
    mae         = mae,
    r2          = r2
  ))
}



#===============================================

# WAIC — MODEL 3 

#===============================================


# -------------------------------------------------------------------------------------
# Likelihood per observation i at posterior sample b:
#
#   log p(y_i | theta^(b)) = log N(y_i | vartheta_i^(b), kappa2_{j(i),q(i)}^(b))
#
#   where vartheta_i^(b) = beta_int^(b) + x_i' * beta_k_E(i)^(b)
#                                       + z_i' * beta_k_M(i)^(b)
#                                       + w_i' * beta_k_D(i)^(b)
#         k_E(i)^(b) = xi_E^(b)[i]   (E-level cluster of obs i at iter b)
#         k_M(i)^(b) = xi_M^(b)[i]   (M-level cluster of obs i at iter b)
#         k_D(i)^(b) = xi_D^(b)[i]   (D-level cluster of obs i at iter b)
#         kappa2      from cadena_filtered$kappa2_jq[[q]][[j]][b]
#
# Inputs to the function:
#   cadena_filtered : output of solve_label_switching_3DP()
#   y               : numeric vector (n_obs)     — observed scores
#   X_E_full        : matrix (n_obs x p_E)       — individual-level covariates
#   Z_full          : matrix (n_obs x p_M)       — municipal-level covariates
#   W_full          : matrix (n_obs x p_D)       — departamental-level covariates
#   obs_to_mun      : integer vector (n_obs)     — flat municipality index per obs
#   mun_map         : list (n_mun)               — maps flat index -> (q, j, start, end, n)
# ---------------------------------------------------------------------------------------

compute_WAIC_3DP <- function(cadena_filtered, y, X_E_full, Z_full, W_full,
                             obs_to_mun, mun_map) {
  
  cat("=== Computing WAIC — Model 3 (3-DP BNP) ===\n\n")
  
  # ------------------------------------------------------------------
  # 0. Dimensions
  # ------------------------------------------------------------------
  
  n_obs      <- length(y)                                        # number of observations
  b_filtered <- length(cadena_filtered$beta_k_E_correct_order)  # number of iterations after filtering by K_mode
  n_mun      <- length(mun_map)                                  # number of municipalities
  
  cat(sprintf("  n observations      : %d\n",   n_obs))
  cat(sprintf("  Posterior samples B : %d\n",   b_filtered))
  cat(sprintf("  Municipalities      : %d\n\n", n_mun))
  
  # ------------------------------------------------------------------
  # 1. Pre-extract kappa2_jq into matrix [n_obs x b_filtered]
  #
  #    kappa2_obs_mat[i, b] = kappa2_{j(i), q(i)}^(b)
  #
  #    All students in the same municipality share the same kappa2_jq.
  #    obs_to_mun[i] -> flat mun index -> mun_map[[.]]$q and $j
  # ------------------------------------------------------------------
  
  # Extract a kappa2_jq matrix with dim [n_obs x b_filtered]
  
  # For each observation, find its (q, j)
  obs_q <- integer(n_obs)   # initialise a vector with department index for each observation
  obs_j <- integer(n_obs)   # initialise a vector with municipality-within-dept index for each observation
  
  # Fill the vectors with department and municipality-within-department ID for observations
  for (i in seq_len(n_obs)) {
    mun_idx  <- obs_to_mun[i]       # flat municipality index
    obs_q[i] <- mun_map[[mun_idx]]$q
    obs_j[i] <- mun_map[[mun_idx]]$j
  }
  
  # Identify unique (q, j) pairs to avoid redundant list lookups
  unique_muns   <- unique(data.frame(q = obs_q, j = obs_j))  # one row per unique mun
  n_unique_muns <- nrow(unique_muns)                          # number of municipalities
  
  # Build a temporary mun-level kappa2 matrix [n_unique_muns x b_filtered]
  # then expand to observations
  kappa2_mun_mat <- matrix(NA_real_, nrow = n_unique_muns, ncol = b_filtered)  # dim [n_mun x b_filtered]
  
  for (r in seq_len(n_unique_muns)) {                                # iterate over municipalities
    q_r <- unique_muns$q[r]                                          # department q
    j_r <- unique_muns$j[r]                                          # mun index within department q
    kappa2_mun_mat[r, ] <- cadena_filtered$kappa2_jq[[q_r]][[j_r]]  # length b_filtered
  }
  
  # Map each observation to its row in kappa2_mun_mat
  obs_mun_row <- match(paste(obs_q, obs_j, sep = "_"),
                       paste(unique_muns$q, unique_muns$j, sep = "_"))
  
  # Expand: kappa2_obs_mat[i, b] = kappa2 for obs i at iteration b
  kappa2_obs_mat <- kappa2_mun_mat[obs_mun_row, , drop = FALSE]  # [n_obs x b_filtered]
  
  # ------------------------------------------------------------------
  # 2. Build log-likelihood matrix [n_obs x b_filtered]. Per iteration b
  #
  #    log_lik_mat[i, b] = log N(y_i | vartheta_i^(b), kappa2_{j(i),q(i)}^(b))
  #
  #    where vartheta_i^(b) = beta_int^(b) + x_i' * beta_k_E(i)^(b)
  #                                        + z_i' * beta_k_M(i)^(b)
  #                                        + w_i' * beta_k_D(i)^(b)
  #    Variance: kappa2_{j(i),q(i)}^(b)
  # ------------------------------------------------------------------
  
  # Compute log-likelihood matrix [n_obs x b_filtered]
  log_lik_mat <- matrix(NA_real_, nrow = n_obs, ncol = b_filtered)  # dim [n_obs x b_filtered]
  
  for (b in seq_len(b_filtered)) {  # iterate over filtered B iterations
    
    # Extract iteration-b quantities
    beta_int_b <- cadena_filtered$beta_int[b]                    # retrieve beta scalar
    xi_E_b     <- cadena_filtered$xi_E_correct_order[b, ]        # E-level cluster assignments [n_obs]
    xi_M_b     <- cadena_filtered$xi_M_correct_order[b, ]        # M-level cluster assignments [n_obs]
    xi_D_b     <- cadena_filtered$xi_D_correct_order[b, ]        # D-level cluster assignments [n_obs]
    beta_k_E_b <- cadena_filtered$beta_k_E_correct_order[[b]]    # E-level list of k_E_mode vectors
    beta_k_M_b <- cadena_filtered$beta_k_M_correct_order[[b]]    # M-level list of k_M_mode vectors
    beta_k_D_b <- cadena_filtered$beta_k_D_correct_order[[b]]    # D-level list of k_D_mode vectors
    
    # --- Fitted mean for all observations: vectorized ---
    
    # E-level: beta_E_per_obs[n_obs x p_E]
    # Row i = beta_k_E_b[[ xi_E_b[i] ]]  (beta_k_E vector for obs i's cluster at iter b)
    beta_E_per_obs <- do.call(rbind, beta_k_E_b[xi_E_b])         # dim [n_obs x p_E]
    
    # M-level: beta_M_per_obs[n_obs x p_M]
    # Row i = beta_k_M_b[[ xi_M_b[i] ]]  (beta_k_M vector for obs i's cluster at iter b)
    beta_M_per_obs <- do.call(rbind, beta_k_M_b[xi_M_b])         # dim [n_obs x p_M]
    
    # D-level: beta_D_per_obs[n_obs x p_D]
    # Row i = beta_k_D_b[[ xi_D_b[i] ]]  (beta_k_D vector for obs i's cluster at iter b)
    beta_D_per_obs <- do.call(rbind, beta_k_D_b[xi_D_b])         # dim [n_obs x p_D]
    
    # Compute the full fitted mean for every observation in iteration b:
    # vartheta_i^(b) = beta_int^(b) + x_i' * beta_k_E(i)^(b)
    #                               + z_i' * beta_k_M(i)^(b)
    #                               + w_i' * beta_k_D(i)^(b)
    mu_b <- beta_int_b +
      rowSums(X_E_full * beta_E_per_obs) +  # x_i' * beta_k_E(i)^(b), length [n_obs]
      rowSums(Z_full   * beta_M_per_obs) +  # z_i' * beta_k_M(i)^(b), length [n_obs]
      rowSums(W_full   * beta_D_per_obs)    # w_i' * beta_k_D(i)^(b), length [n_obs]
    
    # --- Log-likelihood for all observations at iteration b ---
    # log N(y_i | vartheta_i^(b), kappa2_{j(i),q(i)}^(b))
    # Extract variance per observation
    kappa2_b <- kappa2_obs_mat[, b]  # length [n_obs] — one kappa2 per obs
    
    # Compute the log N(mean^(b), variance^(b)) with the normal density formula:
    # = -0.5*log(2*pi) - 0.5*log(kappa2) - 0.5*(y - mu)^2 / kappa2
    # Fill the already created matrix [n_obs x b_filtered]
    log_lik_mat[, b] <- -0.5 * log(2 * pi)         -
      0.5 * log(kappa2_b)         -
      0.5 * (y - mu_b)^2 / kappa2_b  # [n_obs]
  }
  
  # ------------------------------------------------------------------
  # 3. WAIC components from log_lik_mat
  #
  #    lppd_i     = log( (1/B) * sum_b  exp(log_lik_mat[i,b]) )
  #               = log-sum-exp(row i) - log(B)          [numerical stability]
  #
  #    mean_ll_i  = (1/B) * sum_b  log_lik_mat[i,b]
  #               = rowMeans(log_lik_mat)
  #
  #    p_WAIC_i   = 2 * (lppd_i - mean_ll_i)
  #
  #    WAIC       = -2 * sum(lppd_i) + 2 * sum(p_WAIC_i)
  # ------------------------------------------------------------------
  
  # --------   lppd per observation       ---------
  # lppd_i = log( (1/B) * sum_b exp(log_lik_ib) )
  # Use log-sum-exp trick for numerical stability:
  #   log(mean(exp(x))) = max(x) + log( sum(exp(x - max(x))) ) - log(B)
  
  #---------
  # lppd_i = ll* + log(\sum_b exp{ll_ib - ll*})
  #---------
  lppd_vec <- apply(log_lik_mat, 1, function(ll) {
    ll_max <- max(ll)                                         # subtracts the max value -> [0,1]
    ll_max + log(sum(exp(ll - ll_max))) - log(b_filtered)    # length [n_obs]
  })
  
  # --- Posterior mean of log-likelihood per observation ---
  mean_loglik_vec <- rowMeans(log_lik_mat)                    # length [n_obs]
  
  # --- p_WAIC per observation ---
  # = 2 * (log(E[p(y|theta)]) - E[log p(y|theta)])
  # >= 0 always by Jensen's inequality
  p_WAIC_vec <- 2.0 * (lppd_vec - mean_loglik_vec)            # length [n_obs]
  
  # --- Aggregate ---
  lppd   <- sum(lppd_vec)                  # scalar: sum over all n_obs students
  p_WAIC <- sum(p_WAIC_vec)               # scalar: effective number of parameters
  WAIC   <- -2.0 * lppd + 2.0 * p_WAIC
  
  # ------------------------------------------------------------------
  # 4. Results
  # ------------------------------------------------------------------
  
  cat("\n=== WAIC Results — Model 3 (3-DP BNP) ===\n")
  cat(sprintf("  lppd                : %.4f\n", lppd))
  cat(sprintf("  p_WAIC              : %.4f\n", p_WAIC))
  cat(sprintf("  WAIC                : %.4f\n", WAIC))
  
  return(list(
    WAIC            = WAIC,
    lppd            = lppd,
    p_WAIC          = p_WAIC,
    lppd_vec        = lppd_vec,         # per-observation lppd     [n_obs]
    p_WAIC_vec      = p_WAIC_vec,       # per-observation p_WAIC   [n_obs]
    mean_loglik_vec = mean_loglik_vec,  # per-observation mean ll  [n_obs]
    log_lik_mat     = log_lik_mat,      # full matrix if needed    [n_obs x b_filtered]
    n_obs           = n_obs,
    b_filtered      = b_filtered
  ))
}




# =====================================================================================

#@                   3.5 SIMULATION: MODEL 3 DATABASE CREATION

# =====================================================================================


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




