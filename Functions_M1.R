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



##############################################################################################

#@--------                         1. PARAMETRIC MODEL                         ---------------

##############################################################################################


# =====================================================================================

#@                                1.1  GIBBS

# =====================================================================================


#--------------------------

# GIBBS updating functions

#_------------------------

#---------------------------------------------------------------------------------
# The following functions are the Full Conditional Distributions used in the GIBBS
#---------------------------------------------------------------------------------


# ------------------------------
# 1. Update scalar beta (intercept)
# ------------------------------
update_beta <- function(y_vec, x_ijq_full, beta_E, z_jq_full, beta_M, 
                        w_q_full, beta_D, mun_kappa, mun_sizes_local,
                        nu_beta, sigma2_beta, mu_beta) {
  #' Update the scalar intercept parameter beta
  # 
  # INPUTS:
  # y_vec             : Vector of all observations
  # x_ijq_full        : Full design matrix for individual-level covariates
  # beta_E            : Current vector of individual-level coefficients
  # z_jq_full         : Full design matrix for municipality-level covariates
  # beta_M            : Current vector of municipality-level coefficients
  # w_q_full          : Full design matrix for department-level covariates
  # beta_D            : Current vector of department-level coefficients
  # mun_kappa         : Vector of municipality-level variance parameters
  # mun_sizes_local   : Vector with number of observations per municipality
  # nu_beta           : Prior hyperparameter
  # sigma2_beta       : Current variance for beta
  # mu_beta           : Prior mean for beta
  #
  # OUTPUTS:
  # Sampled value of beta from its full conditional distribution
  
  # Compute residual r_beta quickly using full design matrices
  r_beta <- y_vec - as.numeric(x_ijq_full %*% beta_E + z_jq_full %*% beta_M + w_q_full %*% beta_D)
  
  # Build observation-level kappa vector using rep() on flat mun_kappa
  kappa2_jq_vector <- rep(mun_kappa, times = mun_sizes_local)
  
  # Compute sum_v_beta (sum or var) and sum_r_sobre_kappa (vectorized)
  # sum_v_beta = sum_j (n_j / mun_kappa_j)
  sum_v_beta <- sum(mun_sizes_local / mun_kappa)
  sum_r_sobre_kappa <- sum(r_beta / kappa2_jq_vector)
  
  # Posterior variance and mean
  v_beta <- 1 / (sum_v_beta + (1 / sigma2_beta))
  if (!is.finite(v_beta) || v_beta <= 0) v_beta <- 1e-9   # numerical stability
  
  numerador_media <- sum_r_sobre_kappa + (nu_beta / sigma2_beta)
  m_beta <- numerador_media * v_beta
  
  # Sample from posterior
  beta_new <- rnorm(1, mean = m_beta, sd = sqrt(v_beta))
  
  return(beta_new)
}


# ------------------------------
# 2. Update beta_E vector (individual-level coefficients)
# ------------------------------

update_beta_E <- function(n_mun_total, mun_sizes_local, mun_x, mun_z, mun_w, 
                          mun_y, mun_kappa, beta, beta_M, beta_D, mun_xtx,
                          numero_Betas_E, sigma2_E, mu_E, ones_vector_E) {
  # Update the vector of individual-level coefficients beta_E
  #
  # INPUTS
  # n_mun_total       : Total number of municipalities
  # mun_sizes_local   : Vector with number of observations per municipality
  # mun_x             : List of individual-level covariate matrices per municipality
  # mun_z             : List of municipality-level covariate matrices per municipality
  # mun_w             : List of department-level covariate matrices per municipality
  # mun_y             : List of outcome vectors per municipality
  # mun_kappa         : Vector of municipality-level variance parameters
  # beta              : Current scalar intercept
  # beta_M            : Current vector of municipality-level coefficients
  # beta_D            : Current vector of department-level coefficients
  # mun_xtx           : List of precomputed X'X matrices per municipality
  # numero_Betas_E    : Number of individual-level coefficients
  # sigma2_E          : Current variance for beta_E
  # mu_E              : Prior mean vector for beta_E
  # ones_vector_E     : Vector of ones (length = numero_Betas_E)
  #
  # OUTPUT:
  # Sampled vector beta_E from its full conditional distribution
  
  # Initialize sufficient statistics compute sum_var and sum_med
  sum_var_E <- matrix(0, nrow = numero_Betas_E, ncol = numero_Betas_E)
  sum_med_E <- matrix(0, nrow = numero_Betas_E, ncol = 1)
  
  # Loop over municipalities to accumulate sufficient statistics
  for (j in seq_len(n_mun_total)) {
    n_j <- mun_sizes_local[j]
    Xj <- mun_x[[j]]
    Zj <- mun_z[[j]]
    Wj <- mun_w[[j]]
    Yj <- mun_y[[j]]
    kval <- mun_kappa[j]
    if (!is.finite(kval) || kval <= 0) kval <- 1e-9  # numerical stability
    
    # Compute residual r_E: y - beta*1 - Z*beta_M - W*beta_D
    r_E <- Yj - beta * rep(1, n_j) - Zj %*% beta_M - Wj %*% beta_D
    
    # Accumulate weighted cross-products
    # Crossproduct more efficient than t(Xi) %*% r_E
    sum_med_E <- sum_med_E + (1 / kval) * crossprod(Xj, r_E)
    sum_var_E <- sum_var_E + (1 / kval) * mun_xtx[[j]]
  }
  
  # Compute posterior precision matrix
  # Precision matrix: sum_var_E + (1/sigma2_E) * I
  precision_E <- sum_var_E + (1 / sigma2_E) * diag(numero_Betas_E)
  
  # Cholesky decomposition for efficient sampling
  # Cholesky of precision (upper triangular)
  R_E <- chol(precision_E)
  
  # Compute posterior mean
  # RHS for posterior mean
  rhs_E <- sum_med_E + (1 / sigma2_E) * (mu_E * ones_vector_E)
  # Solve for posterior mean m_beta_E: m = solve(precision, rhs) using backsolve (faster)
  tmp <- backsolve(R_E, rhs_E, transpose = TRUE)
  m_beta_E <- backsolve(R_E, tmp)
  
  # Sample from multivariate normal using Cholesky
  # Sample: m + solve(R, z) where z ~ N(0, I)
  z_E <- rnorm(numero_Betas_E)
  sample_noise_E <- backsolve(R_E, z_E)
  beta_E_new <- as.numeric(m_beta_E + sample_noise_E)
  
  return(beta_E_new)
}

# ------------------------------
# 3. Update beta_M vector (Same as beta_e)
# ------------------------------

update_beta_M <- function(n_mun_total, mun_sizes_local, mun_x, mun_z, mun_w, 
                          mun_y, mun_kappa, beta, beta_E, beta_D, mun_ztz,
                          numero_Betas_M, sigma2_M, mu_M, ones_vector_M) {
  # Update the vector of municipality-level coefficients beta_M
  #
  # INPUTS
  # n_mun_total         : Total number of municipalities
  # mun_sizes_local     : Vector with number of observations per municipality
  # mun_x               : List of individual-level covariate matrices per municipality
  # mun_z               : List of municipality-level covariate matrices per municipality
  # mun_w               : List of department-level covariate matrices per municipality
  # mun_y               : List of outcome vectors per municipality
  # mun_kappa           : Vector of municipality-level variance parameters
  # beta                : Current scalar intercept
  # beta_E              : Current vector of individual-level coefficients
  # beta_D              : Current vector of department-level coefficients
  # mun_ztz             : List of precomputed Z'Z matrices per municipality
  # numero_Betas_M      : Number of municipality-level coefficients
  # sigma2_M            : Current variance for beta_M
  # mu_M                : Prior mean vector for beta_M
  # ones_vector_M       : Vector of ones (length = numero_Betas_M)
  # 
  # OUTPUT:
  #Sampled vector beta_M from its full conditional distribution
  
  # Initialize sufficient statistics
  sum_var_M <- matrix(0, nrow = numero_Betas_M, ncol = numero_Betas_M)
  sum_med_M <- matrix(0, nrow = numero_Betas_M, ncol = 1)
  
  # Loop over municipalities to accumulate sufficient statistics
  for (j in seq_len(n_mun_total)) {
    n_j <- mun_sizes_local[j]
    Xj <- mun_x[[j]]
    Zj <- mun_z[[j]]
    Wj <- mun_w[[j]]
    Yj <- mun_y[[j]]
    kval <- mun_kappa[j]
    if (!is.finite(kval) || kval <= 0) kval <- 1e-9 #numerical stability
    
    # Compute residual: y - beta*1 - X*beta_E - W*beta_D
    r_M <- Yj - beta * rep(1, n_j) - Xj %*% beta_E - Wj %*% beta_D
    
    # Accumulate weighted cross-products
    sum_med_M <- sum_med_M + (1 / kval) * crossprod(Zj, r_M)
    sum_var_M <- sum_var_M + (1 / kval) * mun_ztz[[j]]
  }
  
  # Compute posterior precision matrix
  precision_M <- sum_var_M + (1 / sigma2_M) * diag(numero_Betas_M)
  
  # Cholesky decomposition for efficient sampling
  R_M <- chol(precision_M)
  
  # Compute posterior mean
  rhs_M <- sum_med_M + (1 / sigma2_M) * (mu_M * ones_vector_M)
  tmpM <- backsolve(R_M, rhs_M, transpose = TRUE)
  m_beta_M <- backsolve(R_M, tmpM)
  
  # Sample from multivariate normal using Cholesky
  z_M <- rnorm(numero_Betas_M)
  beta_M_new <- as.numeric(m_beta_M + backsolve(R_M, z_M))
  
  return(beta_M_new)
}


# ------------------------------
# 4. Update beta_D vector (same as beta_E)
# ------------------------------
update_beta_D <- function(n_mun_total, mun_sizes_local, mun_x, mun_z, mun_w, 
                          mun_y, mun_kappa, beta, beta_E, beta_M, mun_wtw,
                          numero_Betas_D, sigma2_D, mu_D, ones_vector_D) {
  # Update the vector of department-level coefficients beta_D
  #
  # INPUTS:
  # n_mun_total         : Total number of municipalities
  # mun_sizes_local     : Vector with number of observations per municipality
  # mun_x               : List of individual-level covariate matrices per municipality
  # mun_z               : List of municipality-level covariate matrices per municipality
  # mun_w               : List of department-level covariate matrices per municipality
  # mun_y               : List of outcome vectors per municipality
  # mun_kappa           : Vector of municipality-level variance parameters
  # beta                : Current scalar intercept
  # beta_E              : Current vector of individual-level coefficients
  # beta_M              : Current vector of municipality-level coefficients
  # mun_wtw             : List of precomputed W'W matrices per municipality
  # numero_Betas_D      : Number of department-level coefficients
  # sigma2_D            : Current variance for beta_D
  # mu_D                : Prior mean vector for beta_D
  # ones_vector_D       : Vector of ones (length = numero_Betas_D)
  # 
  # OUTPUT:
  #Sampled vector beta_D from its full conditional distribution
  
  # Initialize sufficient statistics
  sum_var_D <- matrix(0, nrow = numero_Betas_D, ncol = numero_Betas_D)
  sum_med_D <- matrix(0, nrow = numero_Betas_D, ncol = 1)
  
  # Loop over municipalities to accumulate sufficient statistics
  for (j in seq_len(n_mun_total)) {
    n_j <- mun_sizes_local[j]
    Xj <- mun_x[[j]]
    Zj <- mun_z[[j]]
    Wj <- mun_w[[j]]
    Yj <- mun_y[[j]]
    kval <- mun_kappa[j]
    if (!is.finite(kval) || kval <= 0) kval <- 1e-9
    
    # Compute residual: y - beta*1 - X*beta_E - Z*beta_M
    r_D <- Yj - beta * rep(1, n_j) - Xj %*% beta_E - Zj %*% beta_M
    
    # Accumulate weighted cross-products
    sum_med_D <- sum_med_D + (1 / kval) * crossprod(Wj, r_D)
    sum_var_D <- sum_var_D + (1 / kval) * mun_wtw[[j]]
  }
  
  # Compute posterior precision matrix
  precision_D <- sum_var_D + (1 / sigma2_D) * diag(numero_Betas_D)
  
  # Cholesky decomposition for efficient sampling
  R_D <- chol(precision_D)
  
  # Compute posterior mean
  rhs_D <- sum_med_D + (1 / sigma2_D) * (mu_D * ones_vector_D)
  tmpD <- backsolve(R_D, rhs_D, transpose = TRUE)
  m_beta_D <- backsolve(R_D, tmpD)
  
  # Sample from multivariate normal using Cholesky
  z_D <- rnorm(numero_Betas_D)
  beta_D_new <- as.numeric(m_beta_D + backsolve(R_D, z_D))
  
  return(beta_D_new)
}


# ------------------------------
# 5. Update sigma2_beta (variance for scalar beta)
# ------------------------------

update_sigma2_beta <- function(beta, mu_beta, nu_beta, gamma_beta) {
  # Update the variance parameter for the scalar intercept beta
  #
  # INPUT:
  # beta                : Current value of scalar intercept
  # mu_beta             : Prior mean for beta
  # nu_beta             : Prior hyperparameter (degrees of freedom)
  # gamma_beta          : Prior hyperparameter (scale)
  # 
  # OUTPUT:
  # Sampled value of sigma2_beta from its full conditional distribution
  
  a_sig2_beta <- 0.5 * (nu_beta + 1)
  b_sig2_beta <- 0.5 * (nu_beta * gamma_beta + (beta - mu_beta)^2)
  sigma2_beta_new <- 1 / rgamma(1, shape = a_sig2_beta, rate = b_sig2_beta)
  
  return(sigma2_beta_new)
}


# ------------------------------
# 6. Update sigma2_E (variance for beta_E vector)
# ------------------------------
update_sigma2_E <- function(beta_E, mu_E, ones_vector_E, nu_E, gamma_E, numero_Betas_E) {
  # Update the variance parameter for the individual-level coefficients beta_E
  #
  # INPUTS:
  # beta_E            : Current vector of individual-level coefficients
  # mu_E              : Prior mean for beta_E (scalar, applied to all elements)
  # ones_vector_E     : Vector of ones (length = numero_Betas_E)
  # nu_E              : Prior hyperparameter (degrees of freedom)
  # gamma_E           : Prior hyperparameter (scale)
  # numero_Betas_E    : Number of individual-level coefficients
  # 
  # OUTPUT
  #Sampled value of sigma2_E from its full conditional distribution
  
  p_E <- numero_Betas_E
  a_sig2_E <- 0.5 * (nu_E + p_E)
  b_sig2_E <- 0.5 * (nu_E * gamma_E + t(beta_E - mu_E * ones_vector_E) %*% (beta_E - mu_E * ones_vector_E))
  sigma2_E_new <- 1 / rgamma(1, shape = a_sig2_E, rate = b_sig2_E)
  
  return(sigma2_E_new)
}


# ------------------------------
# 7. Update sigma2_M (variance for beta_M vector)
# ------------------------------
update_sigma2_M <- function(beta_M, mu_M, ones_vector_M, nu_M, gamma_M, numero_Betas_M) {
  # Update the variance parameter for the municipality-level coefficients beta_M
  #
  # INPUT:
  # beta_M          : Current vector of municipality-level coefficients
  # mu_M            : Prior mean for beta_M (scalar, applied to all elements)
  # ones_vector_M   : Vector of ones (length = numero_Betas_M)
  # nu_M            : Prior hyperparameter (degrees of freedom)
  # gamma_M         : Prior hyperparameter (scale)
  # numero_Betas_M  : Number of municipality-level coefficients
  # 
  # OUTPUT:
  # Sampled value of sigma2_M from its full conditional distribution
  
  p_M <- numero_Betas_M
  a_sig2_M <- 0.5 * (nu_M + p_M)
  b_sig2_M <- 0.5 * (nu_M * gamma_M + t(beta_M - mu_M * ones_vector_M) %*% (beta_M - mu_M * ones_vector_M))
  sigma2_M_new <- 1 / rgamma(1, shape = a_sig2_M, rate = b_sig2_M)
  
  return(sigma2_M_new)
}


# ------------------------------
# 8. Update sigma2_D (variance for beta_D vector)
# ------------------------------
update_sigma2_D <- function(beta_D, mu_D, ones_vector_D, nu_D, gamma_D, numero_Betas_D) {
  # Update the variance parameter for the department-level coefficients beta_D
  #
  # INPUT:
  # beta_D         : Current vector of department-level coefficients
  # mu_D           : Prior mean for beta_D (scalar, applied to all elements)
  # ones_vector_D  : Vector of ones (length = numero_Betas_D)
  # nu_D           : Prior hyperparameter (degrees of freedom)
  # gamma_D        : Prior hyperparameter (scale)
  # numero_Betas_D : Number of department-level coefficients
  # 
  # OUTPUT:
  # Sampled value of sigma2_D from its full conditional distribution
  
  p_D <- numero_Betas_D
  a_sig2_D <- 0.5 * (nu_D + p_D)
  b_sig2_D <- 0.5 * (nu_D * gamma_D + t(beta_D - mu_D * ones_vector_D) %*% (beta_D - mu_D * ones_vector_D))
  sigma2_D_new <- 1 / rgamma(1, shape = a_sig2_D, rate = b_sig2_D)
  
  return(sigma2_D_new)
}


# ------------------------------
# 9. Update mun_kappa (municipality-level variances)
# ------------------------------
update_mun_kappa <- function(n_mun_total, mun_y, mun_x, mun_z, mun_w, 
                             mun_sizes_local, beta, beta_E, beta_M, beta_D,
                             nu_kappa, kappa2_q, mun_map) {
  # Update municipality-level variance parameters (kappa2_jq)
  #
  # INPUT:
  # n_mun_total       : Total number of municipalities
  # mun_y             : List of outcome vectors per municipality
  # mun_x             : List of individual-level covariate matrices per municipality
  # mun_z             : List of municipality-level covariate matrices per municipality
  # mun_w             : List of department-level covariate matrices per municipality
  # mun_sizes_local   : Vector with number of observations per municipality
  # beta              : Current scalar intercept
  # beta_E            : Current vector of individual-level coefficients
  # beta_M            : Current vector of municipality-level coefficients
  # beta_D            : Current vector of department-level coefficients
  # nu_kappa          : Prior hyperparameter
  # kappa2_q          : Vector of department-level variance parameters
  # mun_map           : List mapping municipalities to departments
  # 
  # OUTPUT:
  # Updated vector of municipality-level variances
  
  mun_kappa_new <- numeric(n_mun_total)
  
  for (j in seq_len(n_mun_total)) {
    Yj <- mun_y[[j]]
    Xj <- mun_x[[j]]
    Zj <- mun_z[[j]]
    Wj <- mun_w[[j]]
    n_j <- mun_sizes_local[j]
    
    # Predicted mean for municipal observations
    zeta <- beta * rep(1, n_j) + Xj %*% beta_E + Zj %*% beta_M + Wj %*% beta_D
    # residual
    y_minus_zeta <- Yj - zeta
    
    # Posterior parameters
    a_kappa2_jq <- 0.5 * (nu_kappa + n_j)
    b_kappa2_jq <- 0.5 * ((nu_kappa * kappa2_q[mun_map[[j]]$q]) + crossprod(y_minus_zeta))
    
    # Sample and store
    mun_kappa_new[j] <- max(1e-9, 1 / rgamma(1, shape = a_kappa2_jq, rate = b_kappa2_jq))
  }
  
  return(mun_kappa_new)
}


# ------------------------------
# 10. Update kappa2_q (department-level variances)
# ------------------------------
update_kappa2_q <- function(Y, mun_kappa, nu_kappa, alpha_kappa, beta_kappa, mun_map) {
  # Update department-level variance parameters (kappa2_q)
  #
  # INPUT:
  # Y                   : Nested list structure (departments -> municipalities)
  # mun_kappa           : Vector of municipality-level variance parameters
  # nu_kappa            : Prior hyperparameter
  # alpha_kappa         : Prior hyperparameter (shape)
  # beta_kappa          : Prior hyperparameter (rate)
  # mun_map             : List mapping municipalities to departments
  # 
  # OUTPUT:
  # Updated vector of department-level variances
  
  m <- length(Y)  # number of departments
  kappa2_q_new <- numeric(m)
  # For each department q, compute sum_inv = sum(1 / mun_kappa[municipalities in q])
  # We can find municipalities in department q via mun_map
  for (q in seq_along(Y)) {
    # collect indices of municipalities that belong to dept q
    inds_q <- which(sapply(mun_map, function(mm) mm$q) == q)
    sum_inv <- sum(1 / mun_kappa[inds_q])
    n_q <- length(inds_q)
    
    # Posterior parameters
    a_kappa2_q <- 0.5 * (nu_kappa * n_q + alpha_kappa)
    b_kappa2_q <- 0.5 * (beta_kappa + (nu_kappa * sum_inv))
    
    # Sample and store
    kappa2_q_new[q] <- max(1e-9, rgamma(1, shape = a_kappa2_q, rate = b_kappa2_q))
  }
  
  return(kappa2_q_new)
}


#===============================================

#         GIBBS SAMPLER FUNCTION (MODEL 1)

#===============================================

# -----------------------------
# MCMC1 function
# -----------------------------
MCMC1_optimized <- function(B, Y, mu_beta, mu_E, mu_M, mu_D, nu_beta,
                            nu_E, nu_M, nu_D, gamma_beta, gamma_E,
                            gamma_M, gamma_D, nu_kappa, alpha_kappa, beta_kappa) {
  
  # --- sampling config 
  
  # Configuración de burn-in y thinning
  burn_in <- ceiling(B * 0.5)         # 50% iterations as burn-in
  thin <- 4                           # store every 3 iterations
  
  # Number of total iterations to execute
  total_iter <- burn_in + (B * thin)
  # Number of samples to store afther thining and burn-in
  n_samples <- B
  
  #=====================     INITIALIZATION   =====================
  
  # --- Initialize priors and starting values 
  # Sufficient statistics
  
  # mean
  theta <- yb_m
  #var
  sigma2_q <- s2_m
  # departamental var
  kappa2_q <- sapply(Y, function(depto) {
    mean(sapply(depto, function(mun) var(mun, na.rm = TRUE) + 0.01), na.rm = TRUE)
  })
  mu <- mean(yb_m)
  tau2 <- var(yb_m)
  sigma2 <- mean(s2_m)
  
  # sigma initialization
  sigma2_beta <- 2
  sigma2_E <- 2
  sigma2_M <- 2
  sigma2_D <- 2
  
  # beta initialization in OLS
  beta <- mu_beta
  beta_E <- matrix(mu_E, nrow = numero_Betas_E, ncol = 1)
  beta_M <- matrix(mu_M, nrow = numero_Betas_M, ncol = 1)
  beta_D <- matrix(mu_D, nrow = numero_Betas_D, ncol = 1)
  
  # initialize kappa2_jq nested and flat representation: Sufficient statistics
  # Var per municipality
  kappa2_jq <- lapply(Y, function(depto) lapply(depto, function(x) var(x) + 0.01))
  # flattened mun_kappa vector in same order as 'mun_map'
  mun_kappa <- numeric(n_mun_total)
  idx <- 1L
  for (i in seq_len(n_mun_total)) {
    q <- mun_map[[i]]$q
    j <- mun_map[[i]]$j
    mun_kappa[i] <- kappa2_jq[[q]][[j]]
  }
  
  # Precompute mun_sizes vector
  mun_sizes_local <- mun_sizes
  
  
  #=====================     STORAGE STRUCTURES   =====================
  
  # Pre-allocated storage structure. Only for final samples
  cadena <- list()
  
  #For kappa
  cadena$kappa2_jq <- vector("list", length = m)
  cadena$kappa2_q  <- vector("list", length = m)
  
  #For betas
  cadena$beta = rep(NA, n_samples)
  cadena$beta_E = vector("list", length = numero_Betas_E)
  names(cadena$beta_E) <- x_ijq_names
  cadena$beta_M = vector("list", length = numero_Betas_M)
  names(cadena$beta_M) <- z_jq_names
  cadena$beta_D = vector("list", length = numero_Betas_D)
  names(cadena$beta_D) <- w_q_names
  
  # For beta's var (sigmas)
  cadena$sigma2_beta = rep(NA, n_samples)
  cadena$sigma2_E = rep(NA, n_samples)
  cadena$sigma2_M = rep(NA, n_samples)
  cadena$sigma2_D = rep(NA, n_samples)
  
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
    # 1) Update scalar beta
    # ------------------------------
    
    beta <- update_beta(y_vec = y_vec, x_ijq_full = x_ijq_full, beta_E = beta_E,
                        z_jq_full = z_jq_full, beta_M = beta_M, w_q_full = w_q_full, beta_D = beta_D,
                        mun_kappa = mun_kappa, mun_sizes_local = mun_sizes_local, nu_beta = nu_beta, 
                        sigma2_beta = sigma2_beta, mu_beta = mu_beta
    )
    
    # ------------------------------
    # 2) Update beta_E (vector) (precision-cholesky approach)
    # ------------------------------
    
    beta_E <- update_beta_E( n_mun_total = n_mun_total, mun_sizes_local = mun_sizes_local,
                             mun_x = mun_x, mun_z = mun_z, mun_w = mun_w, mun_y = mun_y, mun_kappa = mun_kappa,
                             beta = beta, beta_M = beta_M, beta_D = beta_D, mun_xtx = mun_xtx,
                             numero_Betas_E = numero_Betas_E, sigma2_E = sigma2_E, mu_E = mu_E,
                             ones_vector_E = ones_vector_E
    )
    
    # ------------------------------
    # 3) Update beta_M (vector) (same as beta_E)
    # ------------------------------
    
    beta_M <- update_beta_M(n_mun_total, mun_sizes_local, mun_x, mun_z, mun_w, 
                            mun_y, mun_kappa, beta, beta_E, beta_D, mun_ztz,
                            numero_Betas_M, sigma2_M, mu_M, ones_vector_M)
    
    # ------------------------------
    # 4) Update beta_D (vector) (same as beta_E)
    # ------------------------------
    
    beta_D <- update_beta_D(n_mun_total, mun_sizes_local, mun_x, mun_z, mun_w, 
                            mun_y, mun_kappa, beta, beta_E, beta_M, mun_wtw,
                            numero_Betas_D, sigma2_D, mu_D, ones_vector_D)
    
    # ------------------------------
    # 5) Update sigma2_beta 
    # ------------------------------
    
    sigma2_beta <- update_sigma2_beta(beta, mu_beta, nu_beta, gamma_beta)
    
    # ------------------------------
    # 6) Update sigma2_E 
    # ------------------------------
    
    sigma2_E <- update_sigma2_E(beta_E, mu_E, ones_vector_E, nu_E, gamma_E, numero_Betas_E)
    
    # ------------------------------
    # 7) Update sigma2_M
    # ------------------------------
    
    sigma2_M <- update_sigma2_M(beta_M, mu_M, ones_vector_M, nu_M, gamma_M, numero_Betas_M)
    
    # ------------------------------
    # 8) Update sigma2_D 
    # ------------------------------
    
    sigma2_D <- update_sigma2_D(beta_D, mu_D, ones_vector_D, nu_D, gamma_D, numero_Betas_D)
    
    # ------------------------------
    # 9) Update mun_kappa (flat loop over municipalities)
    # ------------------------------
    
    mun_kappa <- update_mun_kappa(n_mun_total, mun_y, mun_x, mun_z, mun_w, 
                                  mun_sizes_local, beta, beta_E, beta_M, beta_D,
                                  nu_kappa, kappa2_q, mun_map)
    
    # ------------------------------
    # 10) Update kappa2_q (per department) using mun_kappa
    # ------------------------------
    
    kappa2_q <- update_kappa2_q(Y, mun_kappa, nu_kappa, alpha_kappa, beta_kappa, mun_map)
    
    #=====================     SAMPLE STORAGE   =====================
    
    # ------------------------------
    # Save results after burn-in & thinning
    # ------------------------------
    if (b > burn_in && (b - burn_in) %% thin == 0) {
      sample_count <- sample_count + 1L
      
      # store beta scalar
      
      cadena$beta[sample_count] <- beta
      
      # Store betas vectors
      
      for (B_i in seq_len(numero_Betas_E)) {
        cadena$beta_E[[ x_ijq_names[B_i] ]][sample_count] <- beta_E[B_i]
      }
      for (B_i in seq_len(numero_Betas_M)) {
        cadena$beta_M[[ z_jq_names[B_i] ]][sample_count] <- beta_M[B_i]
      }
      for (B_i in seq_len(numero_Betas_D)) {
        cadena$beta_D[[ w_q_names[B_i] ]][sample_count] <- beta_D[B_i]
      }
      
      # Store sigmas
      
      cadena$sigma2_beta[sample_count] <- sigma2_beta
      cadena$sigma2_E[sample_count]    <- sigma2_E
      cadena$sigma2_M[sample_count]    <- sigma2_M
      cadena$sigma2_D[sample_count]    <- sigma2_D
      
      # Copy department-level kappas
      
      for (q in seq_len(m)) {
        cadena$kappa2_q[[q]][sample_count] <- kappa2_q[q]
      }
      
      # Copy municipality kappas from flat mun_kappa into nested storage
      # The mapping from flat i to nested q/j is in mun_map
      for (i in seq_len(n_mun_total)) {
        q <- mun_map[[i]]$q
        j <- mun_map[[i]]$j
        cadena$kappa2_jq[[q]][[j]][sample_count] <- mun_kappa[i]
      }
    }
    # Progress indicator
    if (b %% 1000 == 0) cat("Iteración", b, "de", total_iter, "completada\n")
  }
  
  #Aditional information of the sampler
  
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

#@                                1.2 INFERENCE

# =====================================================================================



# ===================================================================
#         FUNCTION 1: compute_posterior_summaries_M1
# ===================================================================
# Computes posterior summaries (mean, q025, q50, q975) for
# beta (intercept), beta_E, beta_M, beta_D from the MCMC chain.
#
# INPUT:
#   Resultados  : List returned by MCMC1_optimized containing the stored chains
#
# OUTPUT:
#   Resultados  : Same list with added elements:
#                   $beta_int_summary  — list(mean, q025, q50, q975) [scalar each]
#                   $beta_E_summary    — list(mean, q025, q50, q975) [vector length p_E]
#                   $beta_M_summary    — list(mean, q025, q50, q975) [vector length p_M]
#                   $beta_D_summary    — list(mean, q025, q50, q975) [vector length p_D]
# ===================================================================

compute_posterior_summaries_M1 <- function(Resultados) {
  
  # ------------------------------------------------------------------
  # 1. beta intercept  (scalars)
  # ------------------------------------------------------------------
  
  beta_int_summary <- list(
    mean = mean(Resultados$beta),
    q025 = quantile(Resultados$beta, probs = 0.025),
    q50  = quantile(Resultados$beta, probs = 0.500),
    q975 = quantile(Resultados$beta, probs = 0.975)
  )
  
  # ------------------------------------------------------------------
  # 2. beta_E  (list)
  #    each column is one coefficient's chain across iterations
  # ------------------------------------------------------------------
  
  # do.call(cbind, ...) stacks them into a [B x p_E] matrix so that
  beta_E_mat <- do.call(cbind, Resultados$beta_E)    # dim [B x p_E]
  
  # Compute the summary of beta_E
  beta_E_summary <- list(
    mean = colMeans(beta_E_mat),
    q025 = apply(beta_E_mat, 2, quantile, probs = 0.025),
    q50  = apply(beta_E_mat, 2, quantile, probs = 0.500),
    q975 = apply(beta_E_mat, 2, quantile, probs = 0.975)
  )
  
  # ------------------------------------------------------------------
  # 3. beta_M  (list)
  # ------------------------------------------------------------------
  
  beta_M_mat <- do.call(cbind, Resultados$beta_M)   # [B x p_M]
  
  beta_M_summary <- list(
    mean = colMeans(beta_M_mat),
    q025 = apply(beta_M_mat, 2, quantile, probs = 0.025),
    q50  = apply(beta_M_mat, 2, quantile, probs = 0.500),
    q975 = apply(beta_M_mat, 2, quantile, probs = 0.975)
  )
  
  
  # ------------------------------------------------------------------
  # 4. beta_D  (list)
  # ------------------------------------------------------------------
  
  beta_D_mat <- do.call(cbind, Resultados$beta_D)   # [B x p_D]
  
  beta_D_summary <- list(
    mean = colMeans(beta_D_mat),
    q025 = apply(beta_D_mat, 2, quantile, probs = 0.025),
    q50  = apply(beta_D_mat, 2, quantile, probs = 0.500),
    q975 = apply(beta_D_mat, 2, quantile, probs = 0.975)
  )
  
  
  # ------------------------------------------------------------------
  # 5. Store summaries inside Resultados and return
  # ------------------------------------------------------------------
  
  Resultados$beta_int_summary <- beta_int_summary
  Resultados$beta_E_summary   <- beta_E_summary
  Resultados$beta_M_summary   <- beta_M_summary
  Resultados$beta_D_summary   <- beta_D_summary
  
  return(Resultados)
}



# =====================================================================================

#@                                1.3 MODEL VALIDATION

# =====================================================================================



#===============================================

# EXTERNAL VALIDATION: MSE, MAE, R2, coverage

#===============================================


# ===================================================================
#         FUNCTION 2: compute_test_metrics_M1
# ===================================================================
# Computes MSE, MAE, R2 and 95% credibility interval coverage on a
# held-out test dataset.
#
# MSE / MAE / R2  — computed using posterior mean point predictions
#   y_hat_i = beta_int_mean + x_i' beta_E_mean + z_i' beta_M_mean + w_i' beta_D_mean
#
# Coverage — for each test observation i, the credibility interval
#   [q_0.025, q_0.975] is built from the B linear predictors
#   zeta_i^(b) = beta^(b) + x_i' beta_E^(b) + z_i' beta_M^(b) + w_i' beta_D^(b)
#   Coverage = % of y_test observations that fall inside their interval.
#
# ASSUMPTIONS:
#   - compute_posterior_summaries_M1(Resultados) has already been called,
#     so Resultados$beta_int_summary, $beta_E_summary, etc. exist.
#   - datos_testeo has the same column structure as the training data.
#   - x_ijq_names, z_jq_names, w_q_names are the character vectors that
#     select the covariate columns (same ones used during training).
#
# INPUTS:
#   Resultados    : MCMC chain list (with posterior summaries already added)
#   datos_testeo  : Raw test data frame (same structure as training data)
#   x_ijq_names   : Character vector — individual-level covariate names
#   z_jq_names    : Character vector — municipal-level covariate names
#   w_q_names     : Character vector — department-level covariate names
#
# OUTPUT: list with
#   $y_hat    — point predictions (posterior mean), length n_test
#   $y_test   — observed scores, length n_test
#   $q025     — lower credibility bound per observation
#   $q975     — upper credibility bound per observation
#   $in_CI    — logical vector: TRUE if y_test[i] in [q025[i], q975[i]]
#   $coverage — % of observations covered (scalar)
#   $mse, $rmse, $mae, $r2
# ===================================================================

compute_test_metrics_M1 <- function(Resultados, datos_testeo,
                                    x_ijq_names, z_jq_names, w_q_names) {
  
  # ------------------------------------------------------------------
  # 1. Prepare test data (mirrors the preprocessing in the main script)
  # ------------------------------------------------------------------
  
  n_test <- nrow(datos_testeo)
  y_test <- datos_testeo$punt_global
  
  # Build full covariate matrices  [n_test x p_*]
  X_test <- as.matrix(subset(datos_testeo, select = x_ijq_names))  # individual variables
  Z_test <- as.matrix(subset(datos_testeo, select = z_jq_names))   # municipal variables
  W_test <- as.matrix(subset(datos_testeo, select = w_q_names))    # departmental variables
  
  # ------------------------------------------------------------------
  # 2. Extract full posterior chains as matrices  [B x p_*]
  # ------------------------------------------------------------------
  
  # Extract number of samples in the MCMC
  B <- Resultados$info$samples_stored
  
  beta_samples <- Resultados$beta                         # length B       intercept 
  beta_E_mat   <- do.call(cbind, Resultados$beta_E)       # dim [B x p_E]  individual variables
  beta_M_mat   <- do.call(cbind, Resultados$beta_M)       # dim [B x p_M]  municipal variables
  beta_D_mat   <- do.call(cbind, Resultados$beta_D)       # dim [B x p_D]  departamental variables
  
  # ------------------------------------------------------------------
  # 3. Build zeta matrix for coverage  [n_test x B] (y_hat for iteration)
  #
  #   y_i^(b) = zeta_i^(b) = beta^(b) * 1  +  x_i' beta_E^(b)
  #                              +  z_i' beta_M^(b)
  #                              +  w_i' beta_D^(b)
  # ------------------------------------------------------------------
  
  # Compute a matrix with rows = number of observations and cols = number of iterations
  # it has the computed y_hat for each iteration
  
  zeta_matrix <- outer(rep(1.0, n_test), beta_samples) +    # dim [n_test x B]
    X_test %*% t(beta_E_mat)                             +  # dim [n_test x B]
    Z_test %*% t(beta_M_mat)                             +  # dim [n_test x B]
    W_test %*% t(beta_D_mat)                                # dim [n_test x B]
  
  # dim(zeta_matrix) = n_test x B
  
  # ------------------------------------------------------------------
  # 4. Credibility intervals and coverage
  #    For each observation i, compute q025 and q975 across B columns.
  #    Coverage = % of y_test[i] that fall in [q025[i], q975[i]].
  # ------------------------------------------------------------------
  
  # apply over rows (each row = one observation, B draws across columns)
  # 1 to compute per row
  q025 <- apply(zeta_matrix, 1, quantile, probs = 0.025)  # length n_test
  q975 <- apply(zeta_matrix, 1, quantile, probs = 0.975)  # length n_test
  
  # Compute coverage:
  # Compute as bolean if the y_test observation falls in the CI
  in_CI    <- (y_test >= q025) & (y_test <= q975)
  # Compute the percentage of observations that falls in the CI
  coverage <- mean(in_CI)
  
  # ------------------------------------------------------------------
  # 5. Point predictions using posterior means  ->  MSE, MAE, R2
  #    y_hat_i = beta_int_mean + x_i' beta_E_mean + z_i' beta_M_mean
  #                            + w_i' beta_D_mean
  # ------------------------------------------------------------------
  
  b0 <- Resultados$beta_int_summary$mean
  bE <- as.numeric(Resultados$beta_E_summary$mean)  # length p_E
  bM <- as.numeric(Resultados$beta_M_summary$mean)  # length p_M
  bD <- as.numeric(Resultados$beta_D_summary$mean)  # length p_D
  
  y_hat     <- as.numeric(b0 + X_test %*% bE + Z_test %*% bM + W_test %*% bD)
  residuals <- y_test - y_hat
  
  mse <- mean(residuals^2)
  mae <- mean(abs(residuals))
  r2  <- 1 - sum(residuals^2) / sum((y_test - mean(y_test))^2)
  
  # ------------------------------------------------------------------
  # 6. Print and return results
  # ------------------------------------------------------------------
  
  cat("=== Test Metrics — Model 1 ===\n")
  cat(sprintf("  n_test     : %d\n",     n_test))
  cat(sprintf("  B samples  : %d\n",     B))
  cat(sprintf("  MSE        : %.4f\n",   mse))
  cat(sprintf("  RMSE       : %.4f\n",   sqrt(mse)))
  cat(sprintf("  MAE        : %.4f\n",   mae))
  cat(sprintf("  R2         : %.4f\n",   r2))
  cat(sprintf("  Coverage   : %.2f\n", coverage))
  
  return(list(
    y_hat    = y_hat,        # point predictions (posterior mean), length n_test
    y_test   = y_test,       # observed scores,                    length n_test
    q025     = q025,         # lower 2.5% credibility bound,       length n_test
    q975     = q975,         # upper 97.5% credibility bound,      length n_test
    in_CI    = in_CI,        # logical: is y_test[i] in interval?  length n_test
    coverage = coverage,     # % observations covered              scalar
    mse      = mse,
    rmse     = sqrt(mse),
    mae      = mae,
    r2       = r2
  ))
}


#===============================================

# WAIC

#===============================================



# ===================================================================
#                     WAIC - WATANABE-AKAIKE INFORMATION CRITERION
# ===================================================================
# Formula (from image):
#   WAIC = -2*lppd + 2*p_WAIC
#
#   lppd   = sum_i log( (1/B) * sum_b p(y_i | theta^(b)) )
#
#   p_WAIC = 2 * sum_i [ log((1/B)*sum_b p(y_i|theta^b)) - (1/B)*sum_b log p(y_i|theta^b) ]
# ===================================================================

compute_WAIC <- function(Resultados, mun_map, mun_y, mun_x, mun_z, mun_w) {
  
  cat("=== Computing WAIC ===\n")
  
  B           <- Resultados$info$samples_stored
  n_mun_total <- length(mun_map)
  
  # ------------------------------------------------------------------
  # 1. Pre-extract posterior samples as matrices 
  # ------------------------------------------------------------------
  
  # beta scalar: length-B vector
  beta_samples <- Resultados$beta                       # [B]
  
  # beta_E: B x p_E matrix (each column = one coefficient's chain)
  # row = iterations, columns = covariates p_e
  beta_E_mat   <- do.call(cbind, Resultados$beta_E)    # [B x p_E]
  beta_M_mat   <- do.call(cbind, Resultados$beta_M)    # [B x p_M]
  beta_D_mat   <- do.call(cbind, Resultados$beta_D)    # [B x p_D]
  
  # kappa2 per municipality: n_mun x B matrix
  # Flatten nested list kappa2_jq[[q]][[j]] [1:B] into rows
  # rows = j municipalities , columns = B iterations
  kappa2_mat <- matrix(NA_real_, nrow = n_mun_total, ncol = B)
  for (mun in seq_len(n_mun_total)) {   #iterate over mun
    # Maps municipality department index and mun index inside the department
    q <- mun_map[[mun]]$q
    j <- mun_map[[mun]]$j
    # Retrieves samples from the municipality and stores in the kappa2_mat matrix
    kappa2_mat[mun, ] <- Resultados$kappa2_jq[[q]][[j]]  # B samples for mun (q,j)
  }
  
  # ------------------------------------------------------------------
  # 2. Loop over municipalities, accumulate lppd and p_WAIC
  # ------------------------------------------------------------------
  
  #Initialises count
  lppd   <- 0.0
  p_WAIC <- 0.0
  
  # Progress tracking alert
  progress_step <- max(1L, ceiling(n_mun_total / 10))
  
  for (mun_idx in seq_len(n_mun_total)) {
    
    # ---- Extract municipality data ----
    Yj    <- mun_y[[mun_idx]]          # observations: length n_j 
    Xj    <- mun_x[[mun_idx]]          # individual covariates: dim[n_j x p_E]
    Zj    <- mun_z[[mun_idx]]          # municipal covariates:  dim[n_j x p_M]
    Wj    <- mun_w[[mun_idx]]          # dept covariates:       dim[n_j x p_D]
    n_j   <- length(Yj)                # num of observations
    
    kappa2_j <- kappa2_mat[mun_idx, ]  # municipality variance: B samples
    
    # ---- Compute predicted mean zeta: [n_j x B] ----
    # zeta_{ij}^(b) = beta^(b) + X_i * beta_E^(b) + Z_i * beta_M^(b) + W_i * beta_D^(b)
    # Each term:
    #   outer(1_nj, beta_samples)    ->  n_j x B  (replicate scalar beta across rows)
    #   Xj %*% t(beta_E_mat)         ->  n_j x B  (X: n_j x p_E, beta_E: B x p_E)
    
    # Computes the matrix of means zeta for observation inside the municipality with its samples.
    # zeta_mat has dim (n_j x B)
    zeta_mat <- outer(rep(1.0, n_j), beta_samples) +   # n_j x B
      Xj %*% t(beta_E_mat)                 + # n_j x B
      Zj %*% t(beta_M_mat)                 + # n_j x B
      Wj %*% t(beta_D_mat)                   # n_j x B
    
    # ---- Compute pointwise log-likelihood: [n_j x B] ----
    # log p(y_i | theta^(b)) = log N(y_i | zeta^(b), kappa2^(b))
    # Vectorized normal log-density:
    #   -0.5*log(2*pi) - 0.5*log(kappa2) - 0.5*(y - zeta)^2 / kappa2
    
    # Expand kappa2_j [B] -> [n_j x B] (same variance for all students in municipality)
    kappa2_mat_j <- matrix(kappa2_j, nrow = n_j, ncol = B, byrow = TRUE)  # n_j x B with identical rows
    
    # Residuals [n_j x B]: y_i - zeta_{ij}^(b)
    resid_mat <- matrix(Yj, nrow = n_j, ncol = B) - zeta_mat              # n_j x B
    
    # Log-likelihood matrix [n_j x B]  Computed from the normal directly for all values
    #instead of using rnorm for each individual value
    log_lik_mat <- -0.5 * log(2 * pi) -
      0.5 * log(kappa2_mat_j)  -
      0.5 * resid_mat^2 / kappa2_mat_j                       # n_j x B
    
    # ---- lppd per observation ----
    # lppd_i = log( (1/B) * sum_b exp(log_lik_ib) )
    # Use log-sum-exp trick for numerical stability:
    #   log(mean(exp(x))) = max(x) + log( sum(exp(x - max(x))) ) - log(B)
    
    #---------
    # lppd_i = ll* + log(\sum_b exp{ll_ib - ll*})
    #---------
    
    lppd_vec <- apply(log_lik_mat, 1, function(ll) {
      ll_max <- max(ll)                               # Subtracts the max value -> [0,1]
      ll_max + log(sum(exp(ll - ll_max))) - log(B)    # length n_j
    })                                                                      
    
    #---- Compute llpd with the ean of lppd_i (lppd per observation i)
    # Mean log-likelihood per observation (1/B) * sum_b log p(y_i | theta^b) ----
    mean_loglik_vec <- rowMeans(log_lik_mat)          # length n_j
    
    # ---- p_WAIC per observation ----
    # p_WAIC_i = 2 * [ log(mean_b p(y_i|theta^b)) - mean_b(log p(y_i|theta^b)) ]
    #          = 2 * [ lppd_i - mean_loglik_i ]
    p_WAIC_vec <- 2.0 * (lppd_vec - mean_loglik_vec)                       # length n_j
    
    # ---- Accumulate ----
    lppd   <- lppd   + sum(lppd_vec)
    p_WAIC <- p_WAIC + sum(p_WAIC_vec)
    
    # Progress message
    if (mun_idx %% progress_step == 0 || mun_idx == n_mun_total) {
      cat(sprintf("  Municipality %d of %d (%.0f%%)\n",
                  mun_idx, n_mun_total, 100 * mun_idx / n_mun_total))
    }
  }
  
  # ------------------------------------------------------------------
  # 3. Final WAIC
  # ------------------------------------------------------------------
  
  WAIC <- -2.0 * lppd + 2.0 * p_WAIC
  
  results <- list(
    WAIC   = WAIC,
    lppd   = lppd,
    p_WAIC = p_WAIC,
    B      = B,
    n_obs  = sum(sapply(mun_y, length))
  )
  
  cat("\n=== WAIC Results ===\n")
  cat(sprintf("  n observations : %d\n",    results$n_obs))
  cat(sprintf("  B samples used : %d\n",    B))
  cat(sprintf("  lppd           : %.4f\n",  lppd))
  cat(sprintf("  p_WAIC         : %.4f\n",  p_WAIC))
  cat(sprintf("  WAIC           : %.4f\n",  WAIC))
  
  return(results)
}


