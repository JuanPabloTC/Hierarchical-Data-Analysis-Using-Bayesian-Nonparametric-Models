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


# Load C++ function
Rcpp::sourceCpp("path/Samplers_M1.cpp")





##############################################################################################

#@--------                   2. NON-PARAMETRIC MODEL: ONE DP                 ---------------

##############################################################################################


# =====================================================================================

#@                                2.1  GIBBS

# =====================================================================================

#---------------------------------------------------------------------------------
# The following functions are the Full Conditional Distributions used in the GIBBS
#---------------------------------------------------------------------------------

# ------------------------------
# 1) Update cluster assignments (xi)
# ------------------------------
# Replaced by its equivalent C++ function
sample_xi_optimized <- function(y, X_full, xi, beta_k, sigma2_base, alpha, mu_vec,
                                mun_kappa, obs_to_mun, xTx_vec, sqrt_kappa, beta_int) {
  # Optimized version of sample_xi with pre-computed values
  #
  # INPUTS:
  # y           : numeric vector (n) - response variable
  # X_full      : numeric matrix n x p - covariates
  # xi          : integer vector (n) - current cluster assignments
  # beta_k      : list of vectors - cluster-specific coefficients
  # sigma2_base : scalar - variance of base distribution (sigma^2)
  # alpha       : scalar - DP concentration parameter
  # mu_vec      : numeric vector (p) - mean of base distribution
  # mun_kappa   : numeric vector - kappa^2_{j,q} per municipality
  # obs_to_mun  : integer vector (n) - pre-computed municipality indices
  # xTx_vec     : numeric vector (n) - pre-computed x_i^T x_i values
  # sqrt_kappa  : numeric vector - pre-computed sqrt(kappa^2) values
  # beta_int    : scalar - global intercept beta
  #
  # OUTPUTS:
  # xi          : list with (updated assignments) and beta_k (updated coefficients)
  
  n <- length(y)
  p <- ncol(X_full)
  
  # Pre-compute constants used in new cluster sampling
  sigma2_inv <- 1 / sigma2_base
  sigma2_diag <- sigma2_inv * diag(p)
  
  # Progress tracking for large datasets
  if (n > 1000) {
    cat("Starting sample_xi for", n, "observations...\n")
  }
  
  for (i in seq_len(n)) {
    
    # Progress indicator
    if (n > 10000 && i %% 10000 == 0) {
      cat("  Processed", i, "/", n, "observations (", 
          round(100*i/n, 1), "%)\n", sep = "")
    }
    
    # --- Get pre-computed values for observation i ---
    mun_idx <- obs_to_mun[i]          # lookup instead of search
    kappa2_i <- mun_kappa[mun_idx]
    sqrt_kappa_i <- sqrt_kappa[mun_idx]
    x_i <- X_full[i, ]
    xTx_i <- xTx_vec[i]
    
    # --- Get cluster assignments excluding observation i ---
    xi_minus_i <- xi[-i]
    
    # --- Get unique clusters and counts ---
    # Use table() instead of relabeling all observations
    cluster_table <- table(xi_minus_i)
    unique_clusters <- as.integer(names(cluster_table))
    cluster_counts <- as.integer(cluster_table)
    K_new <- length(unique_clusters)
    
    # --- Allocate log probability vector ---
    # new clusters
    log_probs <- numeric(K_new + 1) #last slot = new cluster
    
    # --- Probabilities for EXISTING clusters ---
    if (K_new > 0) {
      for (k_idx in seq_len(K_new)) {
        cluster_id <- unique_clusters[k_idx]
        n_k <- cluster_counts[k_idx]
        
        # Get beta for this cluster
        beta_k_vec <- beta_k[[cluster_id]]
        
        # Mean: beta_int + x_i^T beta_k
        mean_k <- beta_int + sum(x_i * beta_k_vec)
        
        # Log probability: log(n_{-i,k}) + log N(y_i | mean_k, kappa^2)
        log_probs[k_idx] <- log(n_k) + 
          dnorm(y[i], mean = mean_k, sd = sqrt_kappa_i, log = TRUE)
      }
    }
    
    # --- Probability for NEW cluster ---
    # Variance: kappa^2_{j,q} + sigma^2 * x_i^T x_i
    var_new <- kappa2_i + (sigma2_base * xTx_i)
    # Mean: beta_int + x_i^T mu
    mean_new <- beta_int + sum(x_i * mu_vec)
    # Log probability: log(alpha) + log N(y_i | mean_new, var_new)
    log_probs[K_new + 1] <- log(alpha) + 
      dnorm(y[i], mean = mean_new, sd = sqrt(var_new), log = TRUE)
    
    # --- Sample new cluster assignment ---
    # Normalize probabilities (log-sum-exp trick for numerical stability)
    max_log <- max(log_probs)
    probs <- exp(log_probs - max_log)   # subtracts max value for stability
    probs <- probs / sum(probs)
    
    # sample cluster index
    sampled_idx <- sample.int(K_new + 1, size = 1, prob = probs)
    
    # --- Update cluster assignment ---
    if (sampled_idx <= K_new) {
      # Assign to existing cluster (use original cluster ID)
      xi[i] <- unique_clusters[sampled_idx]
      
    } else {
      # Create new cluster
      # Find next available cluster ID
      if (length(beta_k) == 0) {
        new_cluster_id <- 1L
      } else {
        used_ids <- unique(xi)
        max_id <- max(used_ids)
        # Find first gap or use max+1
        all_ids <- seq_len(max_id + 1L)
        available <- setdiff(all_ids, used_ids)  # Looks for gapps, if exists uses min(gap)
        #If gap exists, uses the smallest missing label otherwise use max_id + 1
        new_cluster_id <- if (length(available) > 0) min(available) else (max_id + 1L)
      }
      
      # Sample new beta_{K+1} from posterior
      t_i <- y[i] - beta_int
      kappa2_inv <- 1 / kappa2_i
      
      # Precision matrix: A = (1/kappa_j,q^2) x_i x_i^T + (1/sigma^2) I_p
      A <- kappa2_inv * tcrossprod(x_i) + sigma2_diag
      
      # Optimisation: Use Cholesky throughout instead of solve()
      # Cholesky decomposition: A = R^T R
      R <- chol(A)
      
      # Mean vector: M = A^{-1} b
      # where b = (1/kappa^2) t_i x_i + (1/sigma^2) mu
      # Note: notice that the sum of kappa2_inv does not appear because we are sampling
      # based on observation i, the only observation in that cluster)
      b <- kappa2_inv * t_i * x_i + sigma2_inv * mu_vec
      
      # Solve A M = b using backsolve: M = R^{-1} (R^{-T} b)
      M <- backsolve(R, backsolve(R, b, transpose = TRUE))
      
      # Sample: beta_new ~ N(M, A^{-1}) = M + R^{-1} z
      z <- rnorm(p)
      #beta_new = M + R^{-1} z 
      beta_new <- as.numeric(M + backsolve(R, z))
      
      # Store new beta at new cluster ID
      beta_k[[new_cluster_id]] <- beta_new
      
      # Assign observation to new cluster
      xi[i] <- new_cluster_id
    }
  }
  
  # --- Final cleanup: relabel clusters consecutively ---
  #Only done ONCE at the end, not n times
  used_clusters <- sort(unique(xi))
  K_final <- length(used_clusters)
  
  # Extract only used beta_k
  beta_k_clean <- vector("list", K_final)
  for (idx in seq_along(used_clusters)) {
    beta_k_clean[[idx]] <- beta_k[[used_clusters[idx]]]
  }
  
  # Relabel xi to be 1, 2, 3, ..., K
  xi_final <- as.integer(factor(xi, levels = used_clusters, 
                                labels = seq_along(used_clusters)))
  
  if (n > 1000) {
    cat("Completed! Final number of clusters:", K_final, "\n")
  }
  
  return(list(xi = xi_final, beta_k = beta_k_clean))
}


# ------------------------------
# 2) Update atoms (beta_k)
# ------------------------------

sample_beta_k <- function(y, X_full, xi, beta_k, beta_int, sigma2_base, mu_vec,
                          kappa_inv, obs_to_mun) {
  # Update cluster-specific coefficient vectors beta_k
  #
  # INPUTS:
  # y           : numeric vector (n) - response variable
  # X_full      : numeric matrix n x p - covariates
  # xi          : integer vector (n) - current cluster assignments
  # beta_k      : list of vectors - current cluster-specific coefficients
  # beta_int    : scalar - global intercept beta
  # sigma2_base : scalar - variance of base distribution (sigma^2)
  # mu_vec      : numeric vector (p) - mean of base distribution
  # kappa_inv   : numeric vector - pre-computed 1/kappa^2 per municipality
  # obs_to_mun  : integer vector (n) - pre-computed municipality indices
  #
  # OUTPUTS:
  # beta_k      : list of updated coefficient vectors
  
  K <- length(beta_k)
  p <- ncol(X_full)
  
  # Pre-compute constant
  sigma2_inv <- 1 / sigma2_base
  sigma2_I <- sigma2_inv * diag(p)
  sigma2_mu <- sigma2_inv * mu_vec
  
  for (k in 1:K) {
    # Find observations in cluster k
    idx_k <- which(xi == k)
    n_k <- length(idx_k)
    
    # Skip empty clusters (shouldn't happen, but defensive programming)
    if (n_k == 0) next
    
    # Get data for cluster k
    X_k <- X_full[idx_k, , drop = FALSE]
    y_k <- y[idx_k]
    
    # Residuals: t_i = y_i - beta
    t_k <- y_k - beta_int
    
    # Get weights: 1/kappa^2 for each observation in cluster k
    kappa_inv_k <- kappa_inv[obs_to_mun[idx_k]]
    
    # Optimisation: Vectorized weighted sums
    # Weighted X^TX: sum_i (1/kappa^2_i) * x_i x_i^T
    # = X^T * diag(weights) * X
    # = (X .* weights)^T * X (more efficient, avoids creating diagonal matrix)
    X_weighted <- X_k * kappa_inv_k  # Element-wise multiply each row by weight
    XtX_weighted <- t(X_weighted) %*% X_k
    
    # Weighted X^T y: sum_i (1/kappa^2_i) * x_i * t_i
    Xty_weighted <- t(X_k) %*% (kappa_inv_k * t_k)
    
    # Precision matrix: A = sum_i (1/kappa^2_i) x_i x_i' + (1/sigma^2) I_p
    A <- XtX_weighted + sigma2_I
    
    # Optimisation: Use Cholesky decomposition instead of solve()
    R <- chol(A)
    
    # Mean vector: M = V * b
    # where b = sum_i (1/kappa^2_i) x_i (y_i - beta) + (1/sigma^2) mu
    b <- Xty_weighted + sigma2_mu
    M <- backsolve(R, backsolve(R, b, transpose = TRUE))
    
    # Sample: beta_k ~ N(M, V) where V = A^{-1}
    # Using: beta_k = M + R^{-1} * z where z ~ N(0, I)
    z <- rnorm(p)
    beta_k[[k]] <- as.numeric(M + backsolve(R, z))
  }
  
  return(beta_k)
}

# ------------------------------
# 3) Update concentration parameter (alpha)
# ------------------------------


sample_alpha <- function(alpha, xi, a_alpha, b_alpha) {
  
  # Update concentration parameter
  #
  # INPUTS:
  # alpha           : concentration parameter
  # xi              : vector of cluster asignments (1xn)
  # a_alpha         : parameter of alpha ~ Gamma(a_alpha, b_alpha)
  # a_alpha         : parameter of alpha ~ Gamma(a_alpha, b_alpha)
  #
  # OUTPUTS:
  # alpha           : updated concentration parameter
  
  K <- length(unique(xi))
  n <- length(xi)
  
  eta <- rbeta(1, shape1 = alpha + 1, shape2 = n)
  pi_eta <- (a_alpha + K - 1) / (a_alpha + K - 1 + n*(b_alpha - log(eta)))
  
  if (runif(1) < pi_eta) {
    return(rgamma(1, shape = a_alpha + K, rate = b_alpha - log(eta)))
  } else {
    return(rgamma(1, shape = a_alpha + K - 1, rate = b_alpha - log(eta)))
  }
}



# ------------------------------
# 4) Update mean of base distribution (mu_vec)
# ------------------------------

sample_mu <- function(beta_k, sigma2_base, eta_mu, nu2_mu_inv) {
  
  # Update mean (\mu) of base distribution mu_vec
  #
  # INPUTS:
  # beta_k      : coefficient vectors
  # sigma2_base : variance of beta_k ~ N_p(\mu, \sigma^2 I)
  # eta_mu      : mean of mu_vec ~ N_p(eta_mu, nu2_mu)
  # nu2_mu_inv  : inverse of variance of nu2_mu mu_vec ~ N_p(eta_mu, nu2_mu)
  #
  # OUTPUTS:
  # mu_vec      : updated means vector 
  
  K <- length(beta_k)
  
  # Stack all beta_k into matrix (K x p) and sum columns. Stack by rows
  beta_matrix <- do.call(rbind, beta_k)
  sum_beta <- colSums(beta_matrix)
  
  sigma2_inv <- 1 / sigma2_base
  v <- 1 / (K * sigma2_inv + nu2_mu_inv)
  # Posterior mean: M = v * (1/σ² * ∑β_k + 1/ν²_μ * η_μ)
  M <- v * (sigma2_inv * sum_beta + nu2_mu_inv * eta_mu)
  
  # Sample from N_p(M, v*I_p)
  # Since V = v*I_p is diagonal, we can sample each component independently:
  # μ = M + sqrt(v) * z where z ~ N(0, I_p)
  mu_vec <- M + sqrt(v) * rnorm(length(M))
  
  return(mu_vec)
}


# ------------------------------
# 5) Update variance of base distribution (sigma2_base)
# ------------------------------

sample_sigma2_base <- function(beta_k, mu_vec, p, a_sigma2, a_b_sigma2) {
  # Update variance (\sigma^2) of the base distribution
  #
  # INPUTS:
  # beta_k      : list of K vectors - cluster-specific coefficients
  # mu_vec      : numeric vector (p) - mean of base distribution
  # p           : integer - dimension of coefficient vectors
  # a_sigma2    : scalar - prior shape hyperparameter
  # a_b_sigma2  : scalar - pre-computed a_σ² * b_σ²
  #
  # OUTPUTS:
  # sigma2      : scalar - updated variance of base distribution
  
  K <- length(beta_k)
  
  # Posterior shape: (Kp + a_σ²)/2
  shape_post <- (K * p + a_sigma2) / 2
  
  # Optimisation: Vectorized computation of sum of squared deviations
  # Method 1: Stack into matrix and compute efficiently
  beta_matrix <- do.call(rbind, beta_k)  # K x p matrix
  
  # Center: subtract mu from each row
  centered <- sweep(beta_matrix, 2, mu_vec, "-")
  
  # Sum of (β_k - μ)ᵀ(β_k - μ) = sum of all squared elements
  sum_sq_dev <- sum(centered^2)
  
  # Posterior rate: (a_σ² b_σ² + Σ(β_k - μ)ᵀ(β_k - μ))/2
  rate_post <- (a_b_sigma2 + sum_sq_dev) / 2
  
  # Sample from Inverse Gamma
  # In R: rgamma gives Gamma, so IG(α,β) = 1/Gamma(α,β)
  sigma2 <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
  
  return(sigma2)
}

# ------------------------------
# 6) Update intercept \beta (beta_int)
# ------------------------------

sample_beta_int <- function(y, X_full, xi, beta_k, mu_beta, sigma2_beta, 
                            kappa_inv, obs_to_mun) {
  # Update global intercept β
  #
  # INPUTS:
  # y           : numeric vector (n) - response variable
  # X_full      : numeric matrix n x p - covariates
  # xi          : integer vector (n) - cluster assignments
  # beta_k      : list of vectors - cluster-specific coefficients
  # mu_beta     : scalar - prior mean for \beta (beta_int)
  # sigma2_beta : scalar - variance for \beta (beta_int)
  # kappa_inv   : numeric vector - pre-computed 1/\kappa2_jq per municipality
  # obs_to_mun  : integer vector (n) - municipality index for each observation
  #
  # OUTPUTS:
  # beta_int    : scalar - updated global intercept
  
  n <- length(y)
  
  # Optimisation: Vectorized computation of residuals
  # r_i = y_i - x_i^T \beta_{k(i)}
  
  # For each observation, get x_i^T β_{k(i)}
  fitted_values <- numeric(n)
  for (i in 1:n) {
    k_i <- xi[i]                            # Cluster of observation i
    x_i <- X_full[i, ]                      # Covariates for observation i
    beta_k_i <- beta_k[[k_i]]               # Coefficient vector for cluster k(i)
    fitted_values[i] <- sum(x_i * beta_k_i) # x_i^T β_{k(i)}
  }
  
  # Residuals: r_i^(β) = y_i - x_i^T β_{k(i)}
  residuals <- y - fitted_values
  
  # Get 1/kappa2_{j(i),q(i)} for each observation
  kappa_inv_i <- kappa_inv[obs_to_mun]
  
  # Posterior precision (inverse variance)
  precision_post <- sum(kappa_inv_i) + 1/sigma2_beta
  
  # Posterior variance: V = 1 / precision_post
  V <- 1 / precision_post
  
  # M
  # Weighted sum of residuals \sum(r_i/kappa2_jq):
  weighted_sum_r <- sum(residuals * kappa_inv_i)
  
  # Posterior mean: M = V · (Σ_i r_i/κ²_i + μ_β/σ²_β)
  M <- V * (weighted_sum_r + mu_beta / sigma2_beta)
  
  # Sample from N(M, V)
  beta_int <- rnorm(1, mean = M, sd = sqrt(V))
  
  return(beta_int)
}


# ------------------------------
# 7) Update sigma2_beta 
# ------------------------------


sample_sigma2_beta <- function(mu_beta, beta_int, nu_beta, gamma_beta){
  # Update beta_int variance sigma2_beta
  #
  # INPUTS:
  # mu_beta       : scalar - prior mean for \beta (beta_int)
  # beta_int      : scalar - global intercept 
  # nu_beta       : shape of GI in sigma2_beta
  # gamma_beta    : rate of GI in sigma2_beta
  #
  # OUTPUTS:
  # sigma2_beta   : Variance of global intercept
  
  a_sig2_beta <- 0.5 * (nu_beta + 1)
  b_sig2_beta <- 0.5 * (nu_beta * gamma_beta + (beta_int - mu_beta)^2)
  sigma2_beta <- 1 / rgamma(1, shape = a_sig2_beta, rate = b_sig2_beta)
  
  return(sigma2_beta)
}


# ------------------------------
# 8) Update mun_kappa (kappa2_{j,q})
# ------------------------------

sample_kappa2_jq <- function(y, X_full, xi, beta_k, beta_int, mun_map, 
                             kappa2_q, nu_kappa) {
  # Update municipality-level variances kappa_{j,q}
  #
  # INPUTS:
  # y           : numeric vector (n) - response variable
  # X_full      : numeric matrix n x p - covariates
  # xi          : integer vector (n) - cluster assignments
  # beta_k      : list - cluster-specific coefficients
  # beta_int    : scalar - global intercept β
  # mun_map     : list - municipality mapping (contains q, j, start, end, n)
  # kappa2_q    : numeric vector (m) - department-level hyperparameters
  # nu_kappa    : scalar - hyperparameter
  #
  # OUTPUTS:
  # kappa2_jq   : numeric vector - updated municipality variances κ²_{j,q}
  
  n_mun_total <- length(mun_map)
  kappa2_jq <- numeric(n_mun_total)
  
  for (mun_idx in seq_len(n_mun_total)) {
    # Get municipality info
    q <- mun_map[[mun_idx]]$q          # Department index
    j <- mun_map[[mun_idx]]$j          # Municipality index within department
    n_jq <- mun_map[[mun_idx]]$n       # Number of observations
    
    # Get observation indices I_{j,q} for this municipality
    idx_jq <- mun_map[[mun_idx]]$start:mun_map[[mun_idx]]$end
    
    # Extract data for municipality j,q
    y_jq <- y[idx_jq]                         # Responses
    X_jq <- X_full[idx_jq, , drop = FALSE]    # Covariates
    xi_jq <- xi[idx_jq]                       # Cluster assignments
    
    # Optimisation: Vectorized computation of fitted values
    # For each observation i in I_{j,q}: β + x_i^T β_{k(i)}
    
    # Stack beta_{k(i)} vectors for all observations in this municipality
    beta_matrix_jq <- do.call(rbind, beta_k[xi_jq])  # n_{j,q} x p matrix
    
    # Compute fitted values: β + X_{j,q} · beta_{k(i)}
    # rowSums(X_jq * beta_matrix_jq) computes x_i^T β_{k(i)} for each i
    fitted_jq <- beta_int + rowSums(X_jq * beta_matrix_jq)
    
    # Residuals: r_i = y_i - β - x_i^T β_{k(i)}
    residuals_jq <- y_jq - fitted_jq
    
    # Sum of squared residuals: \sum_{i \in I_{j,q}} (y_i - β - x_i^T β_{k(i)})^2
    sum_sq_resid <- sum(residuals_jq^2)
    
    # Posterior parameters for Inverse Gamma
    # shape: a = (ν_κ + n_{j,q}) / 2
    shape_post <- (nu_kappa + n_jq) / 2
    
    # Rate: β = (ν_κ kappa2_q + \sum residuals^2) / 2
    rate_post <- (nu_kappa * kappa2_q[q] + sum_sq_resid) / 2
    
    # Sample from Inverse Gamma: IG(α, β) = 1/Gamma(α, β)
    kappa2_jq[mun_idx] <- 1 / rgamma(1, shape = shape_post, rate = rate_post)
  }
  
  return(kappa2_jq)
}


# ------------------------------
# 9) Update kappa2_q (per department) 
# ------------------------------



sample_kappa2_q <- function(mun_kappa, dept_to_mun, nu_kappa, alpha_kappa, 
                            beta_kappa, m) {
  # Update department-level variance hyperparameters κ²_q
  #
  # INPUTS:
  # mun_kappa    : numeric vector - municipality variances κ²_{j,q}
  # dept_to_mun  : list (m) - dept_to_mun[[q]] = indices of municipalities in dept q
  #                           PRE-COMPUTED outside the function
  # nu_kappa     : hyperparameter ν_κ
  # alpha_kappa  : prior shape α_κ
  # beta_kappa   : prior rate β_κ
  # m            : number of departments
  #
  # OUTPUTS:
  # kappa2_q     : numeric vector (m) - updated department hyperparameters
  
  kappa2_q <- numeric(m)
  
  for (q in 1:m) {
    # Get pre-computed municipality indices for department q
    mun_indices_q <- dept_to_mun[[q]]
    
    # Number of municipalities in department q
    n_q <- length(mun_indices_q)
    
    # Sum of 1/kappa2_{j,q} for all municipalities j in department q
    sum_inv_kappa <- sum(1 / mun_kappa[mun_indices_q])
    
    # Posterior parameters for Gamma distribution
    # Shape: a = (ν_kappa n_q + alpha_kappa) / 2
    shape_post <- (nu_kappa * n_q + alpha_kappa) / 2
    
    # Rate: b = beta_kappa / 2 + (ν_kappa/2) \sum 1/kappa^2_{j,q}
    rate_post <- (beta_kappa / 2) + (nu_kappa / 2) * sum_inv_kappa
    
    # Sample from Gamma(α, β)
    kappa2_q[q] <- rgamma(1, shape = shape_post, rate = rate_post)
  }
  
  return(kappa2_q)
}



# ----------------------------------------------------------
# 10) Helper functions to read back xi and beta_k from files 
# ----------------------------------------------------------

# Read all xi samples from file
# Returns: matrix (n_samples x n_observations)
read_xi_samples <- function(file_path) {
  xi_matrix <- read.table(file_path, header = FALSE)
  return(as.matrix(xi_matrix))
}

# Read all beta_k samples from file
# Returns: list of length n_samples, each element is a list of K vectors of length p
read_beta_k_samples <- function(file_path, p) {
  beta_data <- read.table(file_path, header = FALSE)
  n_samples <- nrow(beta_data)
  
  beta_k_all <- vector("list", n_samples)
  
  for (i in 1:n_samples) {
    row <- as.numeric(beta_data[i, ])
    K <- as.integer(row[1])          # First element is K
    betas <- row[2:(K * p + 1)]      # Remaining elements are the beta vectors
    
    # Reconstruct list of K vectors of length p
    beta_k_all[[i]] <- vector("list", K)
    for (k in 1:K) {
      start <- (k - 1) * p + 1
      end   <- k * p
      beta_k_all[[i]][[k]] <- betas[start:end]
    }
  }
  
  return(beta_k_all)
}


#===============================================

#         GIBBS SAMPLER FUNCTION (MODEL 2)

#===============================================

#-------------------------------------------------------------------------------
# This function runs the GIBBS sampler using as input the functions created 
#above as the Full conditional distributions of each parameter in the model 2
#-------------------------------------------------------------------------------


MCMC2_BNP <- function(B, nu_beta, gamma_beta, nu_kappa, alpha_kappa, beta_kappa) {
  
  # --- sampling configuration ---- 
  
  # Configuración de burn-in y thinning
  burn_in <- ceiling(B * 1)         # 30% iterations as burn-in
  thin <- 7                           # store every 3 iterations
  
  # Number of total iterations to execute
  total_iter <- burn_in + (B * thin)
  # Number of samples to store afther thining and burn-in
  n_samples <- B
  
  
 
  #=====================     FILE CONNECTIONS FOR xi AND beta_k   =====================
  
  # Define file paths
  xi_file      <- file.path(path_text_files, "xi_samples.txt")
  beta_k_file  <- file.path(path_text_files, "beta_k_samples.txt")
  
  # Open connections (open = "w" truncates file if it exists)
  xi_con      <- file(xi_file,     open = "w")
  beta_k_con  <- file(beta_k_file, open = "w")
  
  cat("Files opened:\n")
  cat("  xi:    ", xi_file, "\n")
  cat("  beta_k:", beta_k_file, "\n")
  
  #=====================     STORAGE STRUCTURES   =====================
  
  # Pre-allocated list storage structure. Only for final samples. Parameter samples will be stored
  # as list except from xi and beta_k, which will be stored as text files.
  
  cadena <- list()
  
  # 1. For xi         ---> stored in xi_samples.txt
  
  # 2. For beta_k     ---> stored in beta_k_samples.txt
  
  
  # 3. For alpha 
  
  cadena$alpha <- rep(NA, n_samples)
  
  # 4. For mu_vec (as matrix, each column is a variable)
  
  cadena$mu_vec <- matrix(NA, nrow = n_samples, ncol = p)
  colnames(cadena$mu_vec) <- x_names
  
  # 5. For sigma2_base
  
  cadena$sigma2_base <- rep(NA, n_samples)
  
  # 6. For betas_int
  cadena$beta_int <- rep(NA, n_samples)
  
  # 7. For sigma2_beta
  
  cadena$sigma2_beta <- rep(NA, n_samples)
  
  # 8 - 9 For kappa_jq (mun_kappa), and for kappa_q (departamental var)
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
  
  #sampler count for samples stored
  sample_count <- 0L
  
  # Start sampling iterations
  for (b in 1:total_iter) {
    
    #=====================     PARAMETER UPDATING   =====================
    
    # ------------------------------
    #    NON PARAMETRIC PARAMETERS 
    # ------------------------------
    
    # 1. Update xi (Using C++ function)
    
    result_xi <- sample_xi_optimized_cpp(y, X_full, xi, beta_k, sigma2_base, alpha, 
                                         mu_vec, mun_kappa, obs_to_mun, xTx_vec, 
                                         sqrt_kappa, beta_int)
    xi <- result_xi$xi                  # xi
    beta_k <- result_xi$beta_k          # beta_k for new clusters
    
    
    # 2. Update beta_k
    
    kappa_inv <- 1 / mun_kappa  # Avoid repeated division in the loop
    
    beta_k <- sample_beta_k(y, X_full, xi, beta_k, beta_int, sigma2_base, mu_vec, 
                           kappa_inv, obs_to_mun)
    
    # 3. Update alpha
    
    alpha <- sample_alpha(alpha, xi, a_alpha, b_alpha) 
    
    # 4. Update mu_vec
    
    mu_vec <- sample_mu(beta_k, sigma2_base, eta_mu, nu2_mu_inv)
    
    # 5. Update sigma2_base
    
    sigma2_base = sample_sigma2_base(beta_k, mu_vec, p, a_sigma2, a_b_sigma2)
    
    # ------------------------------
    #    PARAMETRIC PARAMETERS 
    # ------------------------------
    
    # 6. Update beta_int
    
    beta_int <- sample_beta_int(y, X_full, xi, beta_k, mu_beta, sigma2_beta, 
                                kappa_inv, obs_to_mun)
    
    # 7. Update sigma2_beta
    
    sigma2_beta <- sample_sigma2_beta(mu_beta, beta_int, nu_beta, gamma_beta)
    
    # 8. Update kappa2_{j,q} = mun_kappa
    
    mun_kappa <- sample_kappa2_jq(y, X_full, xi, beta_k, beta_int, mun_map, 
                                  kappa2_q, nu_kappa)
    
    mun_kappa <- pmax(mun_kappa, 1e-4)  #this prevents the vector to go to cero
    
    # These are used in sample_xi (1) and sample_beta_k (2) next iteration
    sqrt_kappa <- sqrt(mun_kappa)
    
    # 9. Update kappa2_q
    
    kappa2_q <- sample_kappa2_q(mun_kappa, dept_to_mun, nu_kappa, alpha_kappa, beta_kappa, m)
    
    
    #=====================     SAMPLE STORAGE   =====================
    
    # ------------------------------
    # Save results after burn-in & thinning
    # ------------------------------
    if (b > burn_in && (b - burn_in) %% thin == 0) {
      sample_count <- sample_count + 1L
      
      # 1. Store xi
      
      # One row: n integers separated by spaces
      writeLines(paste(xi, collapse = " "), con = xi_con)
      
      # 2. Store beta_k
      
      # One row: K followed by all beta vectors concatenated:
      # |K|  beta_k[[1]]   beta_k[[2]]   beta_k[[3]]
      # Example row: "3 0.5 1.2 ... 0.8 -0.3 ... 1.1 0.4 ..."
      K_current <- length(beta_k)
      beta_k_flat <- c(K_current, unlist(beta_k))  # K, then all betas in order
      writeLines(paste(beta_k_flat, collapse = " "), con = beta_k_con)
      
      # 3. Store alpha 
      
      cadena$alpha[sample_count] <- alpha
      
      # 4. Store mu_vec (as matrix with cols variables and rows samples)
      
      cadena$mu_vec[sample_count, ] <- mu_vec
      
      # 5. Store sigma2_base
      
      cadena$sigma2_base[sample_count] <- sigma2_base
      
      # 6. Store betas_int
      cadena$beta_int[sample_count] <- beta_int
      
      # 7. Store sigma2_beta
      
      cadena$sigma2_beta[sample_count] <- sigma2_beta
      
      # 8. Store kappa_jq (mun_kappa)
      
      for (i in seq_len(n_mun_total)) {
        q <- mun_map[[i]]$q
        j <- mun_map[[i]]$j
        cadena$kappa2_jq[[q]][[j]][sample_count] <- mun_kappa[i]
      }
      
      # 9. Store kappa_q (departamental var)
      
      for (q in seq_len(m)) {
        cadena$kappa2_q[[q]][sample_count] <- kappa2_q[q]
      }
      
      
    }
    # Progress indicator
    if (b %% 10 == 0) cat("Iteración", b, "de", total_iter, "completada\n")
  }
  
  #---------------------
  #CLOSE FILE CONNECTIONS   
  #---------------------
  
  close(xi_con)
  close(beta_k_con)
  
  cat("\nFiles closed successfully!\n")
  cat("  xi samples stored:    ", sample_count, "rows x", length(xi), "columns\n")
  cat("  beta_k samples stored:", sample_count, "rows (K varies per row)\n\n")
  
  
  #==================    ADITIONAL INFO IN THE SAMPLER ===================#
  
  cadena$info <- list(
    total_iterations = total_iter,
    burn_in = burn_in,
    thin = thin,
    samples_stored = sample_count,
    samples_requested = B
  )
  
  return(cadena)
}


# =====================================================================================

#@                               2.2 INFERENCE

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

filter_cadena_by_k <- function(cadena, k_mode) {
  
  # -------------------------------------------------------
  # Step 1: Identify iterations where number of clusters == k_mode
  # Use both xi and beta_k and cross-check for consistency
  # -------------------------------------------------------
  
  # From xi matrix: unique clusters per row (iteration)
  k_from_xi <- apply(cadena$xi, 1, function(row) length(unique(row)))
  
  # From beta_k list: length of each sublist
  k_from_betak <- sapply(cadena$beta_k, length)
  
  # Cross-check: both should agree on k per iteration. 
  if (!all(k_from_xi == k_from_betak)) {
    # Gives a warning if elements don't agree on some k iteration
    warning("Mismatch between K from xi and K from beta_k. ",
            "Number of mismatched iterations: ",
            sum(k_from_xi != k_from_betak),
            ". Filtering will use k_from_xi as reference.")
  }
  
  # Index of the iterations to keep
  keep_idx <- which(k_from_xi == k_mode)
  # Number of iterations to keep
  n_kept   <- length(keep_idx)
  
  # if (n_kept == 0) {
  #   stop("No iterations found with k_mode = ", k_mode, 
  #        ". Available K values: ", 
  #        paste(sort(unique(k_from_xi)), collapse = ", "))
  # }
  
  cat("Filtering cadena by k_mode =", k_mode, "\n")
  cat("  Total iterations:  ", length(k_from_xi), "\n")
  cat("  Iterations kept:   ", n_kept, "\n")
  cat("  Percentage kept:   ", round(100 * n_kept / length(k_from_xi), 1), "%\n")
  
  # -------------------------------------------------------
  # Step 2: Build filtered cadena
  # -------------------------------------------------------
  
  # Creates a new chan to store
  cadena_filtered <- list()
  
  # 1. xi: keep rows in keep_idx
  cadena_filtered$xi <- cadena$xi[keep_idx, , drop = FALSE]
  
  # 2. beta_k: keep sublists in keep_idx
  cadena_filtered$beta_k <- cadena$beta_k[keep_idx]
  
  # 3. Scalar parameters (vectors of length B)
  cadena_filtered$alpha       <- cadena$alpha[keep_idx]
  cadena_filtered$sigma2_base <- cadena$sigma2_base[keep_idx]
  cadena_filtered$beta_int    <- cadena$beta_int[keep_idx]
  cadena_filtered$sigma2_beta <- cadena$sigma2_beta[keep_idx]
  
  # 4. mu_vec: matrix (B x p), keep rows
  cadena_filtered$mu_vec <- cadena$mu_vec[keep_idx, , drop = FALSE]
  
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
    k_mode          = k_mode,
    kept_iterations = keep_idx,
    b_kept          = n_kept,
    b_total         = length(k_from_xi),
    pct_kept        = round(100 * n_kept / length(k_from_xi), 1)
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

solve_label_switching <- function(cadena_filtered, y, X_full, k_mode) {
  
  set.seed(777)
  n_iter  <- nrow(cadena_filtered$xi)
  
  #--------------------------------
  # Define the posible permutations
  #--------------------------------
  
  # If k_mode > 7 is necesary to subsample the posible permutations to 7!. This is
  #because the computational cost of k>7 is not feasible
  if (k_mode <= 7) {
    # Generate all k_mode! permutations of {1, ..., k_mode}
    permu   <- gtools::permutations(n = k_mode, r = k_mode)
    n_permu <- nrow(permu)
  } else {  #If k_mode > 7
    n_permu <- factorial(7)   # <- this is 7!, i.e 5040
    permu <- matrix(nrow = 0, ncol = k_mode)
    
    # This while loop avoids to have duplicated permutations
    while (nrow(permu) < n_permu) {
      needed <- n_permu - nrow(permu)
      batch  <- t(replicate(needed, sample.int(k_mode)))
      permu  <- unique(rbind(permu, batch))   # unique() drops duplicate rows
    }
    # For defensive programming, in case permu > 5040
    permu <- permu[seq_len(n_permu), ]   # trim to exactly 5,040
  }
  
  cat("Solving label switching problem\n")
  cat("  Filtered iterations:", n_iter, "\n")
  cat("  k_mode:            ", k_mode, "\n")
  cat("  Permutations:      ", n_permu, "\n\n")
  
  # Storage for corrected beta_k (same structure as cadena_filtered$beta_k)
  beta_k_correct_order <- vector("list", n_iter)
  
  xi_correct_order_list <- vector("list", n_iter)
  
  for (b in seq_len(n_iter)) {
    
    if (b %% 100 == 0) cat("  Processing iteration", b, "/", n_iter, "\n")
    
    # --------------------------------------------------
    # Extract iteration-specific quantities
    # --------------------------------------------------
    xi_b       <- cadena_filtered$xi[b, ]     # cluster assignments, length n
    beta_int_b <- cadena_filtered$beta_int[b] # scalar intercept for iteration b
    beta_k_b   <- cadena_filtered$beta_k[[b]] # list of k_mode vectors, each length p
    
    # Build beta matrix (k_mode x p): row k = beta_k_b[[k]] ; col = p.
    beta_matrix <- do.call(rbind, beta_k_b)   # dim = num_k x p
    
    # --------------------------------------------------
    # Try all permutations, compute MSE for each
    # --------------------------------------------------
    
    # create a vector of length = num permutations
    mse_permu <- numeric(n_permu)
    
    # iterate over permutations
    for (s in seq_len(n_permu)) {
      
      # permutation s
      sigma <- permu[s, ]
      
      # Under permutation sigma:
      # observation i with cluster xi_b[i] uses beta_k_b[[sigma[xi_b[i]]]]
      # Vectorized: permuted_beta_matrix[xi_b[i], ] gives the right beta vector
      #reorders the matrix according to the current permutation s (sigma)
      permuted_beta_matrix <- beta_matrix[sigma, ]   # k_mode x p: row k = beta_k_b[[sigma[k]]]
      
      # Assign permuted beta to each observation: (n x p)
      #creates a matrix rows = n, cols = p. Each row corresponds to each observation i, and its values
      #are the corresponding beta_k (from its cluster xi) in the permuted_beta_matrix
      beta_per_obs <- permuted_beta_matrix[xi_b, , drop = FALSE]
      
      # Fitted values: y_hat_i = beta_int_b + x_i^T beta_k_b[[sigma[xi_b[i]]]]
      #length(y_hat) = n
      y_hat <- beta_int_b + rowSums(X_full * beta_per_obs)
      
      # COMPUTE MSE
      #stores the MSE of each permutation s. length(mse-permu) = n of permutations (n_permu)
      mse_permu[s] <- mean((y - y_hat)^2)
    }
    
    # ---------------------------------------
    # Select best permutation (minimum MSE)
    # ---------------------------------------
    # Selects the permutation that have minimum MSE
    best_sigma <- permu[which.min(mse_permu), ]
    
    # Reorder beta_k according to best_sigma:
    # beta_k_correct_order[[b]][[k]] = beta_k_b[[best_sigma[k]]]
    beta_k_correct_order[[b]] <- beta_k_b[best_sigma]
    
    xi_correct_order_list[[b]] <- best_sigma[xi_b]
  }
  
  # Add corrected beta_k to cadena_filtered (same structure as beta_k)
  cadena_filtered$beta_k_correct_order <- beta_k_correct_order
  # Add corrected xi
  cadena_filtered$xi_correct_order <- do.call(rbind, xi_correct_order_list)
  
  cat("\n cadena_filtered$beta_k_correct_order successfully added.\n")
  cat(" cadena_filtered$xi_correct_order successfully added.\n")
  
  return(cadena_filtered)
}



#===============================================

#       COMPUTE POSTERIOR SUMMARIES

#===============================================


compute_posterior_summaries <- function(cadena_filtered, k_mode) {
  
  n_iter <- length(cadena_filtered$beta_k_correct_order)
  p      <- length(cadena_filtered$beta_k_correct_order[[1]][[1]])
  
  cat("Computing posterior summaries\n")
  cat("  Iterations:", n_iter, "\n")
  cat("  k_mode:    ", k_mode, "\n")
  cat("  p:         ", p, "\n\n")
  
  # -------------------------------------------------
  # beta_k summaries
  # -------------------------------------------------
  
  # For each cluster k, build a (n_iter x p) matrix where
  # row b = cadena_filtered$beta_k_correct_order[[b]][[k]]
  # Then compute quantiles and mean column-wise (across iterations for each beta)
  
  beta_k_summary <- vector("list", k_mode) # creates a list, size: number of clusters 
  
  # Iterate for clusters
  
  for (k in seq_len(k_mode)) {
    
    # Stack iteration b, cluster k into matrix: (n_iter x p)
    beta_k_matrix <- do.call(rbind, lapply(cadena_filtered$beta_k_correct_order, 
                                           function(iter_b) iter_b[[k]]))
    # dim(beta_k_matrix) = n_iter x p
    # Each column j contains the n_iter draws for coefficient j of cluster k
    
    # Compute column-wise summaries (each col = one beta coefficient across iterations)
    #computes: Mean, q025, q50, q975 (mean is the estimator and quantiles are used for
    #confidence intervals)
    beta_k_summary[[k]] <- list(
      mean = colMeans(beta_k_matrix),
      q025 = apply(beta_k_matrix, 2, quantile, probs = 0.025),  #Note: 2 for columns
      q50  = apply(beta_k_matrix, 2, quantile, probs = 0.5),
      q975 = apply(beta_k_matrix, 2, quantile, probs = 0.975)
    )
    
  }
  
  # ------------------------------------------------
  # beta_int summaries
  # -----------------------------------------------
  # Computes te mean and quantiles for the beta intercept
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
  
  # store the results in the same chain (list)
  cadena_filtered$beta_k_summary    <- beta_k_summary
  cadena_filtered$beta_int_summary  <- beta_int_summary
  
  cat("--- Beta coefficients mean and quantiles computed --- \n")
  
  return(cadena_filtered)
}



# =====================================================================================

#@                   2.3 SIMULATION: MODEL 2 DATABASE CREATION

# =====================================================================================


#=============================================

#             DATABASE CREATION

#=============================================

#-------------------------------------------------------------------------------
# This function generates the simulated database keeping the hierarchies of 
#individual, municipality and department (3 hierarchies).Generates clusters of
#individuals by their differentiated beta coefficients
#-------------------------------------------------------------------------------

generate_hierarchical_data <- function(
    m,                           # number of departments
    n_mun_per_dept,             # vector of length m: municipalities per department
    n_obs_range,                # c(min, max): range for observations per municipality
    p_d,                        # number of department-level covariates
    p_m,                        # number of municipality-level covariates
    p_i,                        # number of individual-level covariates
    k_true,                     # true number of clusters/atoms
    mix_probs,                  # vector of length k_true: mixing probabilities (must sum to 1)
    mu_centers,                 # vector of length k_true: centers for beta_k atoms
    sigma_k = 1.0,              # sd for sampling beta_k coefficients around centers
    kappa_min = 0.5,            # minimum municipality variance
    kappa_max = 1.5            # maximum municipality variance
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
  
  # Check that mix_probs has length k_true
  if (length(mix_probs) != k_true) {
    stop("Error: length(mix_probs) must equal k_true. ",
         "You provided k_true = ", k_true, " but mix_probs has length ", 
         length(mix_probs))
  }
  
  # Check that mix_probs sums to 1 (with small tolerance for floating point)
  if (abs(sum(mix_probs) - 1.0) > 1e-10) {
    stop("Error: mix_probs must sum to 1. ",
         "Current sum = ", sum(mix_probs))
  }
  
  # Check that mu_centers has length k_true
  if (length(mu_centers) != k_true) {
    stop("Error: length(mu_centers) must equal k_true. ",
         "You provided k_true = ", k_true, " but mu_centers has length ", 
         length(mu_centers))
  }
  
  # Check that n_obs_range has length 2
  if (length(n_obs_range) != 2) {
    stop("Error: n_obs_range must be a vector of length 2: c(min, max)")
  }
  
  # Check that n_obs_range[1] >= 1
  if (n_obs_range[1] < 1) {
    stop("Error: minimum observations per municipality must be >= 1")
  }
  
  set.seed(777)
  
  
  # ===================================================================
  #                    Derived quantities
  # ===================================================================
  
  p <- p_d + p_m + p_i                    # total covariates
  n_mun_total <- sum(n_mun_per_dept)      # total municipalities
  
  # Sample observations per municipality from discrete uniform
  n_obs_per_mun_vec <- sample(n_obs_range[1]:n_obs_range[2], 
                              size = n_mun_total, 
                              replace = TRUE)
  
  n <- sum(n_obs_per_mun_vec)             # total observations
  
  cat("Generating dataset with:\n")
  cat("  Departments (m):", m, "\n")
  cat("  Total municipalities:", n_mun_total, "\n")
  cat("  Total observations (n):", n, "\n")
  cat("  Covariates: p_d =", p_d, ", p_m =", p_m, ", p_i =", p_i, 
      ", total p =", p, "\n")
  cat("  True clusters (k_true):", k_true, "\n\n")
  
  # ===================================================================
  #                   Build hierarchical structure
  # ===================================================================
  
  dept_id <- integer(n)
  mun_id  <- integer(n)
  obs_idx <- 1L:n
  
  mun_map <- vector("list", n_mun_total)
  dept_to_mun <- vector("list", m)
  obs_to_mun <- integer(n)
  
  # Build mapping (order by department -> municipality -> observations)
  obs_counter <- 1L
  mun_counter <- 1L
  
  for (q in seq_len(m)) {
    dept_to_mun[[q]] <- integer(0)
    n_mun_q <- n_mun_per_dept[q]
    
    for (j in seq_len(n_mun_q)) {
      n_jq <- n_obs_per_mun_vec[mun_counter]  # observations for this municipality
      start_i <- obs_counter
      end_i   <- obs_counter + n_jq - 1L
      
      mun_map[[mun_counter]] <- list(
        q = q, 
        j = j, 
        start = start_i, 
        end = end_i, 
        n = n_jq
      )
      dept_to_mun[[q]] <- c(dept_to_mun[[q]], mun_counter)
      
      # Fill IDs
      dept_id[start_i:end_i] <- q
      mun_id[start_i:end_i]  <- mun_counter
      obs_to_mun[start_i:end_i] <- mun_counter
      
      obs_counter <- obs_counter + n_jq
      mun_counter <- mun_counter + 1L
    }
  }
  
  # ===================================================================
  #                      Generate covariates
  # ===================================================================
  
  
  # 1) Individual-level covariates (p_i columns, varies per observation)
  if (p_i > 0) {
    ind_cov <- matrix(rnorm(n * p_i, mean = 0, sd = 1), 
                      nrow = n, ncol = p_i)
  } else {
    ind_cov <- matrix(nrow = n, ncol = 0)
  }
  colnames(ind_cov) <- paste0("Individual", seq_len(p_i))
  
  # 2) Municipality-level covariates (p_m columns, constant per municipality)
  if (p_m > 0) {
    mun_cov_vals <- matrix(rnorm(n_mun_total * p_m, mean = 0, sd = 1), 
                           nrow = n_mun_total, ncol = p_m)
    mun_cov <- mun_cov_vals[mun_id, , drop = FALSE]
  } else {
    mun_cov <- matrix(nrow = n, ncol = 0)
  }
  colnames(mun_cov) <- paste0("Municipal", seq_len(p_m))
  
  # 3) Department-level covariates (p_d columns, constant per department)
  if (p_d > 0) {
    dept_cov_vals <- matrix(rnorm(m * p_d, mean = 0, sd = 1), 
                            nrow = m, ncol = p_d)
    dept_cov <- dept_cov_vals[dept_id, , drop = FALSE]
  } else {
    dept_cov <- matrix(nrow = n, ncol = 0)
  }
  colnames(dept_cov) <- paste0("Departmental_", seq_len(p_d))
  
  
  # Assemble X_full (n x p)
  X_full <- cbind(ind_cov, mun_cov, dept_cov)
  x_names <- colnames(X_full)
  
  # Precompute xTx_vec (used by sampler)
  xTx_vec <- rowSums(X_full^2)
  
  # ===================================================================
  #                   Generate true parameters
  # ===================================================================
  
  # True intercept from N(0, 1)
  beta_int_true <- rnorm(1, mean = 0, sd = 1)
  
  # True atoms (beta_k_true) with offset centers
  # Each cluster k has coefficients centered around mu_centers[k]
  beta_k_true <- vector("list", k_true)
  for (k in seq_len(k_true)) {
    beta_k_true[[k]] <- rnorm(p, mean = mu_centers[k], sd = sigma_k)
  }
  
  # True municipality variances (uniform on [kappa_min, kappa_max])
  mun_kappa <- runif(n_mun_total, min = kappa_min, max = kappa_max)
  
  sqrt_kappa = sqrt(mun_kappa)  
  
  # ===================================================================
  #                   Assign clusters and generate Y
  # ===================================================================
  
  # Assign cluster for each observation based on mix_probs
  xi_true <- sample.int(k_true, size = n, replace = TRUE, prob = mix_probs)
  
  # Build matrix with beta per observation
  beta_mat <- do.call(rbind, beta_k_true)    # k_true x p
  beta_per_obs <- beta_mat[xi_true, ]        # n x p
  
  # Linear predictor: beta_int + X^T beta_{k(i)}
  linpred <- beta_int_true + rowSums(X_full * beta_per_obs)
  
  # Add observation-level noise with variance = mun_kappa[obs_to_mun]
  y <- linpred + rnorm(n, mean = 0, sd = sqrt(mun_kappa[obs_to_mun]))
  
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
    mun  = mun_id,
    y = y
  )
  datos <- cbind(datos, X_full)
  
  # ===================================================================
  #                   Print summary
  # ===================================================================
  
  cat("Summary:\n")
  cat("  Municipality sizes (min, median, max):", 
      min(n_obs_per_mun_vec), median(n_obs_per_mun_vec), max(n_obs_per_mun_vec), "\n")
  cat("  True intercept (beta_int_true):", round(beta_int_true, 3), "\n")
  cat("  mun_kappa range:", round(range(mun_kappa), 3), "\n")
  cat("  Cluster sizes:\n")
  print(table(xi_true))
  cat("\n")
  
  # ===================================================================
  #                     Return objects
  # ===================================================================
  
  #===================================
  # Objects to initialize
  #===================================
  # --- xi
  xi_false <- sample.int(k_true, size = n, replace = TRUE, prob = mix_probs)
  
  # In case some cluster is empty, it deletes it from the betas and from xi
  used <- sort(unique(xi_false))
  xi_false <- as.integer(factor(xi_false, levels = used, labels = seq_along(used)))
  
  # --- beta_k
  beta_k_false <- vector("list", k_true)
  for (k in seq_len(k_true)) {
    beta_k_false[[k]] <- rnorm(p, mean = mu_centers[k], sd = sigma_k)
  }
  # in case of empty clusters, deletes its beta
  beta_k_false <- beta_k_false[used]
  
  #kappa2_q
  kappa2_q <- sapply(seq_len(m), function(q) {
    mun_indices_q <- dept_to_mun[[q]]
    mean(mun_kappa[mun_indices_q])   # mean of true kappa2_{j,q} within dept q
  })
  
  beta_int_false <- rnorm(1, mean = 0, sd = 1)
  
  return(list(
    # Response and covariates
    y = y,
    X_full = X_full,
    x_names = x_names,
    
    # Dimensions
    n = n,
    m = m,
    p = p,
    p_d = p_d,
    p_m = p_m,
    p_i = p_i,
    k_true = k_true,
    n_mun_total = n_mun_total,
    
    # Hierarchical structure
    dept_id = dept_id,
    mun_id = mun_id,
    obs_to_mun = obs_to_mun,
    mun_map = mun_map,
    dept_to_mun = dept_to_mun,
    
    # Precomputed values (used in functions)
    xTx_vec = xTx_vec,
    sqrt_kappa = sqrt_kappa,
    
    # True parameters
    beta_int_true = beta_int_true,
    beta_k_true = beta_k_true,
    mun_kappa = mun_kappa,
    xi_true = xi_true,
    mix_probs = mix_probs,
    
    # Nested Y structure
    Y = Y,
    
    # Variable names
    x_names = x_names,
    
    # Objects to initialize
    xi_false = xi_false,
    beta_k_false = beta_k_false,
    beta_int_false = beta_int_false,
    kappa2_q = kappa2_q,
    
    # Data frame
    datos = datos
  ))
}





# =====================================================================================

#@                          2.4  MODEL VALIDATION

# =====================================================================================

#===============================================

# EXTERNAL VALIDATION: MSE, MAE, R2, coverage

#===============================================

#-------------------------------------------------------------------------
# This function computes prediction testing metrics: MSE, MAE, and R2.
#It assigns to each observation a cluster randomly sampled considering the 
#proportion of observations per cluster in each iteration of the GIBBS 
#sampler
#-------------------------------------------------------------------------

compute_test_metrics <- function(cadena_filtered, test_data) {
  
  k_mode <- cadena_filtered$filter_info$k_mode
  n_iter  <- nrow(cadena_filtered$xi_correct_order)
  n_test  <- nrow(test_data$X_full)
  
  # Fixed posterior mean estimates for beta_intercept (same across all iterations)
  beta_int_hat <- cadena_filtered$beta_int_summary$mean
  
  # Build (k_mode x p) matrix of posterior mean betas: row k = beta_k_summary[[k]]$mean
  #nrow = k, ncol = p (covariates)
  beta_k_hat <- do.call(rbind, lapply(seq_len(k_mode), function(k) {
    cadena_filtered$beta_k_summary[[k]]$mean
  }))
  # dim(beta_k_hat) = k_mode x p
  
  cat("Computing test predictions\n")
  cat("  Test observations:", n_test, "\n")
  cat("  Iterations:       ", n_iter, "\n")
  cat("  k_mode:           ", k_mode, "\n\n")
  
  # -------------------------------------------------------
  # Step 1: For each iteration b, compute y_hat for all
  # test observations using sampled cluster
  # -------------------------------------------------------
  
  # Matrix to store predictions: (n_iter x n_test)
  y_hat_matrix <- matrix(NA, nrow = n_iter, ncol = n_test)
  
  for (b in seq_len(n_iter)) { #Iterates over Gibbs iterations in filtered chain
    
    if (b %% 100 == 0) cat("  Processing iteration", b, "/", n_iter, "\n")
    
    # Cluster proportions from training xi at iteration b: Extract the for w
    #from the matrix xi that contains the cluster assignments for all iterations 
    #in the gibbs.
    xi_b    <- cadena_filtered$xi_correct_order[b, ]
    # Compute the proportions of each cluster k in the row b from matrix xi
    pi_b    <- as.numeric(table(factor(xi_b, levels = 1:k_mode))) / length(xi_b)
    # pi_b is a vector of length k_mode: pi_b[k] = proportion in cluster k
    
    # For each test observation i, sample a cluster k using pi_b as probabilities(sum1)
    #size = n_test samples 1 cluster per observation in 1:k_mode with probability pi_b
    sampled_k <- sample.int(k_mode, size = n_test, replace = TRUE, prob = pi_b)
    # sampled_k[i] = cluster assigned to observation i in iteration b
    
    # Assign corresponding beta_k to each test observation: (n_test x p)
    #nrow = n_test, ncol = p. Each observation have the beta_k corresponding to its 
    #sampled cluster
    beta_per_obs <- beta_k_hat[sampled_k, , drop = FALSE]
    
    # Linear predictor: y_hat_i = beta_int + x_i^T beta_k[sampled_k[i]]
    #matrix nrow = n_iter, ncol = n_test
    y_hat_matrix[b, ] <- beta_int_hat + rowSums(test_data$X_full * beta_per_obs)
  }
  
  # -------------------------------------------------------
  # Step 2: Summarize across iterations for each observation
  # -------------------------------------------------------
  # Compute q0.025, mean and q0.0975
  
  # Prediction matrix: 4 rows x n_test cols
  # Row 1 = mean, Row 2= q0.025, Row 3 = median (q50), Row 4 = q0.975
  
  #The 2 applies by columns
  pred_matrix <- apply(y_hat_matrix, 2, function(obs_draws) {
    c(mean  = mean(obs_draws),
      q025  = quantile(obs_draws, 0.025),
      q50   = quantile(obs_draws, 0.50),
      q975  = quantile(obs_draws, 0.975))
  })
  #assign rownames
  rownames(pred_matrix) <- c("mean", "q025", "q50" , "q975")
  # dim(pred_matrix) = 4 x n_test
  
  # -------------------------------------------------------
  # Step 3: Coverage — does true y fall in [q025, q975]?
  # -------------------------------------------------------
  
  y_test   <- test_data$y    #observed y value
  # assess wether y real value falls in the Confidence Interval 
  #(TRUE if so, FALSE otherwise)
  in_CI    <- (y_test >= pred_matrix["q025", ]) & (y_test <= pred_matrix["q975", ])
  # Coverage is the percentage of observations which real y falls in CI
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
  
  cat("  MSE:", round(mse, 4), "\n")
  cat("  MAE:", round(mae, 4), "\n")
  cat("  R2: ", round(r2,  4), "\n")
  
  # -----------------
  # Return results
  # -----------------
  
  return(list(
    pred_matrix = pred_matrix,      # 4 x n_test matrix (mean, q025,q50, q975)
    in_CI       = in_CI,            # logical vector length n_test
    coverage    = coverage,         # % of y_test inside CI
    mse         = mse,
    rmse        = sqrt(mse),      
    mae         = mae,
    r2          = r2
  ))
}


#===============================================

# WAIC — MODEL 2 (BNP)

#===============================================

# -------------------------------------------------------------------------------------
# Likelihood per observation i at posterior sample b:
#
#   log p(y_i | theta^(b)) = log N(y_i | mu_i^(b), kappa2_{j(i),q(i)}^(b))
#
#   where mu_i^(b) = beta_int^(b) + x_i' * beta_k(i)^(b)
#         k(i)^(b) = xi^(b)[i]               (cluster of obs i at iter b)
#         kappa2    from cadena_filtered$kappa2_jq[[q]][[j]][b]
#
# Inputs to the function:
#   cadena_filtered : output of solve_label_switching()
#   y               : numeric vector (n_obs) — observed scores
#   X_full          : matrix (n_obs x p)     — covariate matrix
#   obs_to_mun      : integer vector (n_obs) — flat municipality index per obs
#   mun_map         : list (n_mun)           — maps flat index -> (q, j, start, end, n)
# ---------------------------------------------------------------------------------------


compute_WAIC_model2 <- function(cadena_filtered, y, X_full, obs_to_mun, mun_map) {
  
  cat("=== Computing WAIC — Model 2 (BNP) ===\n\n")
  
  # ------------------------------------------------------------------
  # 0. Dimensions
  # ------------------------------------------------------------------
  
  n_obs      <- length(y)                                    # number of observations
  b_filtered <- length(cadena_filtered$beta_k_correct_order) # number of iterations after filtering by K_mode
  n_mun      <- length(mun_map)                              # number of municipalities
  
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
  unique_muns    <- unique(data.frame(q = obs_q, j = obs_j))  # one row per unique mun
  n_unique_muns  <- nrow(unique_muns)                         # Number of municipalities
  
  # Build a temporary mun-level kappa2 matrix [n_unique_muns x b_filtered]
  # then expand to observations
  kappa2_mun_mat <- matrix(NA_real_, nrow = n_unique_muns, ncol = b_filtered)  # dim [n_mun x b_filtered]
  
  for (r in seq_len(n_unique_muns)) {                               # Iterate over municipalities
    q_r <- unique_muns$q[r]                                         # department q
    j_r <- unique_muns$j[r]                                         # mun index within department q
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
  #    log_lik_mat[i, b] = log N(y_i | beta_int^(b) + x_i' beta_{k(i)}^(b),
  #                                    kappa2_{j(i),q(i)}^(b))
  #
  #    Mean= beta_int^(b) + x_i' beta_{k(i)}^(b)
  #    Variance: kappa2_{j(i),q(i)}^(b)
  # ------------------------------------------------------------------
  
  #Compute log-likelihood matrix [n_obs x b_filtered]
  log_lik_mat <- matrix(NA_real_, nrow = n_obs, ncol = b_filtered)   # dim [n_obs x b_filtered]
  
  for (b in seq_len(b_filtered)) {                           # Iterate over filtered B iterations
    
    #     Extract iteration-b quantities 
    beta_int_b <- cadena_filtered$beta_int[b]                # retrieve beta scalar
    xi_b       <- cadena_filtered$xi_correct_order[b, ]      # retrieve beta integer vector [n_obs]
    beta_k_b   <- cadena_filtered$beta_k_correct_order[[b]]  # retrieve beta list of k_mode vectors
    
    # --- Fitted mean for all observations: vectorized ---
    # beta_per_obs: [n_obs x p] matrix
    # Row i = beta_k_b[[ xi_b[i] ]]  (the beta vector for obs i's cluster at iter b)
    beta_per_obs <- do.call(rbind, beta_k_b[xi_b])            # dim [n_obs x p]
    
    # Compute the mean of every observation in this iteration b
    # mu_i^(b) = beta_int^(b) + x_i' * beta_{k(i)}^(b)
    mu_b <- beta_int_b + rowSums(X_full * beta_per_obs)       # length [n_obs]
    
    # --- Log-likelihood for all observations at iteration b ---
    # log N(y_i | mu_i^(b), kappa2_{j(i),q(i)}^(b))
    # Extract variance per observation
    kappa2_b <- kappa2_obs_mat[, b]                           # length [n_obs] — one kappa2 per obs
    
    # Compute the log N(mean^(b), variance^(b)) with the normal density formula
    # = -0.5*log(2*pi) - 0.5*log(kappa2) - 0.5*(y - mu)^2 / kappa2
    # Fill the already created  matrix [n_obs x b_filtered]
    log_lik_mat[, b] <- -0.5 * log(2 * pi)          -
      0.5 * log(kappa2_b)          -
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
    ll_max <- max(ll)                                       # Subtracts the max value -> [0,1]
    ll_max + log(sum(exp(ll - ll_max))) - log(b_filtered)   # length [n_obs]
  })                                                    
  
  # --- Posterior mean of log-likelihood per observation ---
  mean_loglik_vec <- rowMeans(log_lik_mat)                  # length [n_obs]
  
  # --- p_WAIC per observation ---
  # = 2 * (log(E[p(y|theta)]) - E[log p(y|theta)])
  # >= 0 always by Jensen's inequality
  p_WAIC_vec <- 2.0 * (lppd_vec - mean_loglik_vec)          # length [n_obs]
  
  # --- Aggregate ---
  lppd   <- sum(lppd_vec)                 # scalar: sum over all n_obs students
  p_WAIC <- sum(p_WAIC_vec)               # scalar: effective number of parameters
  WAIC   <- -2.0 * lppd + 2.0 * p_WAIC
  
  # ------------------------------------------------------------------
  # 4. Results
  # ------------------------------------------------------------------
  
  cat("\n=== WAIC Results — Model 2 ===\n")
  cat(sprintf("  lppd                : %.4f\n", lppd))
  cat(sprintf("  p_WAIC              : %.4f\n", p_WAIC))
  cat(sprintf("  WAIC                : %.4f\n", WAIC))
  
  return(list(
    WAIC            = WAIC,
    lppd            = lppd,
    p_WAIC          = p_WAIC,
    lppd_vec        = lppd_vec,        # per-observation lppd     [n_obs]
    p_WAIC_vec      = p_WAIC_vec,      # per-observation p_WAIC   [n_obs]
    mean_loglik_vec = mean_loglik_vec, # per-observation mean ll  [n_obs]
    log_lik_mat     = log_lik_mat,     # full matrix if needed    [n_obs x b_filtered]
    n_obs           = n_obs,
    b_filtered      = b_filtered
  ))
}
















