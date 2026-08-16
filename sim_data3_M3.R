#@ simulations data: 3. Model: M3

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
ruta <- setwd('C:/Users/jtorres/Desktop/From Documents/JP/Tesis/15-04-26/Simulaciones/sim_data3_M3')

# Path for text files
path_text_files <- paste0(ruta,'/text_files')
# Create dir (folder) if does not exist
if (!dir.exists(path_text_files)) dir.create(path_text_files, recursive = TRUE, showWarnings = FALSE)

# Path for results
path_to_results <- paste0(ruta,'/Resultados')
# Create dir (folder) if does not exist
if (!dir.exists(path_to_results)) dir.create(path_to_results, recursive = TRUE, showWarnings = FALSE)


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

B = 20000

# For alpha_ell
a_alpha_E         <-   1
b_alpha_E         <-   10
a_alpha_M         <-   1
b_alpha_M         <-   5
a_alpha_D         <-   1
b_alpha_D         <-   5

# For beta_k_E base distribution
nu_E              <-   1

# For beta_k_M base distribution
nu_M              <-   1

# For beta_k_D base distribution
nu_D              <-   1

# For beta_int 
nu_beta           <-   1

alpha_E           <- 1
alpha_M           <- 1
alpha_D           <- 1


for (s in 1:18) {
  
  simulated_df_name <- sprintf("simulated_data_%d", s)
  obj <- get(simulated_df_name, envir = .GlobalEnv)
  
  #------------------------------------------------
  #      Retrieve information from the database
  #------------------------------------------------
  
  #Dataframe extraction
  df          <- obj$datos
  m           <- obj$m
  beta_int    <- obj$beta_int_false
  X_E_names   <- obj$X_E_names
  Z_M_names   <- obj$Z_M_names
  W_D_names   <- obj$W_D_names
  Y           <- obj$Y
  y           <- obj$y
  
  p_E         <- obj$p_E
  p_M         <- obj$p_M
  p_D         <- obj$p_D
  
  X_E_full     <- obj$X_E_full
  xi_E         <- obj$xi_E_false
  beta_k_E     <- obj$beta_k_E_false
  
  Z_full       <- obj$Z_full
  xi_M         <- obj$xi_M_false
  beta_k_M     <- obj$beta_k_M_false
  
  W_full       <- obj$W_full
  xi_D         <- obj$xi_D_false
  beta_k_D     <- obj$beta_k_D_false
  
  mun_kappa    <- obj$mun_kappa
  obs_to_mun   <- obj$obs_to_mun
  sqrt_kappa   <- obj$sqrt_kappa
  mun_map      <- obj$mun_map
  dept_to_mun  <- obj$dept_to_mun
  n_mun_total  <- obj$n_mun_total
  kappa2_q     <- obj$kappa2_q
  
  xTx_vec      <- obj$xTx_vec
  zTz_vec      <- obj$zTz_vec
  wTw_vec      <- obj$wTw_vec
  
  mu_E         <- unlist(obj$beta_k_E_false[1])
  mu_M         <- unlist(obj$beta_k_M_false[1])
  mu_D         <- unlist(obj$beta_k_D_false[1])
  
  # VARIANCE COMPONENT HIPERPARAMETERS
  # using CV and 
  alpha_kappa      <-   max(0.1, 2 / (1^2))
  beta_kappa       <-   max(0.1, alpha_kappa / var(y))
  
  nu_kappa         <-   4 + (2 / 1^2)
  
  
  # Unitary priors hiperparameter specification
  lm_previa_int        <-  lm(y ~ 1)
  lm_previa_unitaria_E <-  lm(y ~ X_E_full -1) 
  lm_previa_unitaria_M <-  lm(y ~ Z_full   -1)    
  lm_previa_unitaria_D <-  lm(y ~ W_full   -1)    
  
  
  gamma2_beta       <-   summary(lm_previa_int)$sigma
  
  gamma2_E          <-   summary(lm_previa_unitaria_E)$sigma
  nu_gamma2_E       <-   nu_E*gamma2_E
  
  gamma2_M          <-   summary(lm_previa_unitaria_M)$sigma 
  nu_gamma2_M       <-   nu_M*gamma2_M
  
  gamma2_D          <-   summary(lm_previa_unitaria_D)$sigma
  nu_gamma2_D       <-   nu_D*gamma2_D
  
  mu_beta           <-  coef(lm_previa_int)
  eta_mu_E          <-  coef(lm_previa_unitaria_E)
  eta_mu_M          <-  coef(lm_previa_unitaria_M)
  eta_mu_D          <-  coef(lm_previa_unitaria_D)
  
  nu2_mu_E          <-   summary(lm_previa_unitaria_E)$sigma
  nu2_mu_E_inv      <-   1 / nu2_mu_E
  
  nu2_mu_M          <-   summary(lm_previa_unitaria_M)$sigma 
  nu2_mu_M_inv      <-   1 / nu2_mu_M
  
  nu2_mu_D          <-   summary(lm_previa_unitaria_D)$sigma
  nu2_mu_D_inv      <-   1 / nu2_mu_D
  
  sigma2_beta       <-   summary(lm_previa_int)$sigma
  sigma2_E_base     <-   summary(lm_previa_unitaria_E)$sigma
  sigma2_M_base     <-   summary(lm_previa_unitaria_M)$sigma 
  sigma2_D_base     <-   summary(lm_previa_unitaria_D)$sigma
  
  #------------------------------------------------
  #      Execute the GIBBS sampler
  #------------------------------------------------
  tictoc::tic()
  cadena <- MCMC3_BNP(B,
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
                      nu_kappa, alpha_kappa, beta_kappa)
  
  tiempo = tictoc::toc()
  cadena$info$Tiempo = tiempo$callback_msg
  
  # ======================================
  # LOAD TEXT FILES
  # ======================================
  
  #-------------------------------
  # load xi samples
  #-------------------------------
  
  # File names
  xi_E_file_path <- paste0(path_text_files, "/xi_E_samples.txt")
  xi_M_file_path <- paste0(path_text_files, "/xi_M_samples.txt")
  xi_D_file_path <- paste0(path_text_files, "/xi_D_samples.txt")
  # Load files
  xi_E_samples   <- read_xi_samples(xi_E_file_path)
  xi_M_samples   <- read_xi_samples(xi_M_file_path)
  xi_D_samples   <- read_xi_samples(xi_D_file_path)
  
  cat("xi_E, xi_M, and xi_D loaded: ", nrow(xi_E_samples), "samples", 
      ncol(xi_E_samples), "observations\n\n")
  
  #-------------------------------
  # Load beta_k samples line by line
  #-------------------------------
  
  # File names
  beta_k_E_file_path <- paste0(path_text_files, "/beta_k_E_samples.txt")
  beta_k_M_file_path <- paste0(path_text_files, "/beta_k_M_samples.txt")
  beta_k_D_file_path <- paste0(path_text_files, "/beta_k_D_samples.txt")
  # Load files
  beta_k_E_samples   <- read_beta_k_samples(beta_k_E_file_path, p_E)
  beta_k_M_samples   <- read_beta_k_samples(beta_k_M_file_path, p_M)
  beta_k_D_samples   <- read_beta_k_samples(beta_k_D_file_path, p_D)
  
  
  cat("beta_k (E,M,D) loaded:", length(beta_k_E_samples), "samples\n")
  
  
  # ==================================
  # Add to the chain results
  # ==================================
  
  # Add xi
  cadena$xi_E     <- xi_E_samples
  cadena$xi_M     <- xi_M_samples
  cadena$xi_D     <- xi_D_samples
  
  # Add beta_k 
  cadena$beta_k_E <- beta_k_E_samples
  cadena$beta_k_M <- beta_k_M_samples
  cadena$beta_k_D <- beta_k_D_samples
  
  
  # ==============================================================================
  #@                            LOG-LIKELIHOOD COMPUTATION
  # ==============================================================================
  
  cat("Computing log-likelihood for", length(beta_k_E_samples), "samples...\n")

  # ===================================================================
  #    VALIDATE DATA STRUCTURES
  # ===================================================================
  
  cat("xi_E dimensions   :", nrow(cadena$xi_E), "samples x", ncol(cadena$xi_E), "observations\n")
  cat("beta_k_E samples  :", length(cadena$beta_k_E), "\n")
  cat("beta_int samples  :", length(cadena$beta_int),  "\n\n")
  
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
    # Extract parameters for this sample
    # ------------------------------------------
    
    xi_E_sample     <- cadena$xi_E[sample_count, ]          # E-level cluster assignments (length n)
    xi_M_sample     <- cadena$xi_M[sample_count, ]          # M-level cluster assignments (length n)
    xi_D_sample     <- cadena$xi_D[sample_count, ]          # D-level cluster assignments (length n)
    
    beta_k_E_sample <- cadena$beta_k_E[[sample_count]]      # List of K_E beta_E vectors
    beta_k_M_sample <- cadena$beta_k_M[[sample_count]]      # List of K_M beta_M vectors
    beta_k_D_sample <- cadena$beta_k_D[[sample_count]]      # List of K_D beta_D vectors
    
    beta_int_sample <- cadena$beta_int[sample_count]        # Global intercept beta
    
    # -----------------------------------------------------------------------
    # Compute fitted values for all observations as vector
    # vartheta_{i,j,q} = beta + x_i^T beta_k_E + z_{j,q}^T beta_k_M + w_q^T beta_k_D
    # -----------------------------------------------------------------------
    
    # --- E-level contribution: x_i^T beta_k_E for each observation ---
    # Stack beta_k_E vectors into matrix (K_E x p_E), then index by xi_E_sample
    beta_E_mat     <- do.call(rbind, beta_k_E_sample)                # K_E x p_E
    beta_E_per_obs <- beta_E_mat[xi_E_sample, , drop = FALSE]        # n x p_E
    contrib_E      <- rowSums(X_E_full * beta_E_per_obs)             # length n
    
    # --- M-level contribution: z_{j,q}^T beta_k_M for each observation ---
    # Stack beta_k_M vectors into matrix (K_M x p_M), then index by xi_M_sample
    beta_M_mat     <- do.call(rbind, beta_k_M_sample)                # K_M x p_M
    beta_M_per_obs <- beta_M_mat[xi_M_sample, , drop = FALSE]        # n x p_M
    contrib_M      <- rowSums(Z_full * beta_M_per_obs)               # length n
    
    # --- D-level contribution: w_q^T beta_k_D for each observation ---
    # Stack beta_k_D vectors into matrix (K_D x p_D), then index by xi_D_sample
    beta_D_mat     <- do.call(rbind, beta_k_D_sample)                # K_D x p_D
    beta_D_per_obs <- beta_D_mat[xi_D_sample, , drop = FALSE]        # n x p_D
    contrib_D      <- rowSums(W_full * beta_D_per_obs)               # length n
    
    # --- Full fitted values: vartheta_{i,j,q} ---
    fitted_values <- beta_int_sample + contrib_E + contrib_M + contrib_D  # length n
    
    # ---------------------------------------
    # Compute log-likelihood by municipality
    # ---------------------------------------
    
    total_ll <- 0  # Initialize total log-likelihood for this sample
    
    for (q in 1:m) {
      
      # Number of municipalities in department q
      n_mun_q <- length(Y[[q]])
      
      for (j_local in 1:n_mun_q) {   # Iterate over municipalities in department q
        
        # Get global municipality index
        mun_idx <- dept_to_mun[[q]][j_local]
        
        # Get observation indices for this municipality
        idx_jq <- mun_map[[mun_idx]]$start:mun_map[[mun_idx]]$end
        
        # Extract response and fitted values for municipality (j, q)
        y_jq      <- y[idx_jq]
        fitted_jq <- fitted_values[idx_jq]
        
        # Get kappa2_{j,q} for this sample from the nested list
        kappa2_jq_sample <- cadena$kappa2_jq[[q]][[j_local]][sample_count]
        
        # Compute log-likelihood contribution for this municipality:
        # ll_{j,q} = sum_{i in I_{j,q}} log N(y_i | vartheta_{i,j,q}, kappa2_{j,q})
        ll_jq <- sum(dnorm(x    = y_jq,
                           mean = fitted_jq,
                           sd   = sqrt(kappa2_jq_sample),
                           log  = TRUE))
        
        # Accumulate total log-likelihood
        total_ll <- total_ll + ll_jq
      }
    }
    
    # Store log-likelihood for this sample
    log_lik_vec[sample_count] <- total_ll
  }
  
  # ===================
  # Store to the chain
  # ===================
  
  cadena$log_likelihood <- log_lik_vec
  
  cat("\nLog-likelihood computation completed!\n")

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
  
  # ===============================================
  #@             LOGLIKELIHOOD PLOTS
  # ===============================================
  
  B <- nrow(cadena$xi_E)
  
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
  
  # The current MCMC chain
  cadena <- get(paste0("cadena_", s))
  data   <- get(paste0("simulated_data_", s))
  
  # True number of clusters per level
  k_true_E <- data$k_true_E
  k_true_M <- data$k_true_M
  k_true_D <- data$k_true_D
  
  # Compute number of clusters K at each iteration for each level
  # Note: apply with 1 applies function row-wise (one row = one iteration)
  k_values_E <- apply(cadena$xi_E, 1, function(x) length(unique(x)))
  k_values_M <- apply(cadena$xi_M, 1, function(x) length(unique(x)))
  k_values_D <- apply(cadena$xi_D, 1, function(x) length(unique(x)))
  
  # Compute posterior distribution of K for each level
  k_table_E <- table(k_values_E) / length(k_values_E)
  k_table_M <- table(k_values_M) / length(k_values_M)
  k_table_D <- table(k_values_D) / length(k_values_D)
  
  # ----------------------------------------------------------------
  # Plot 1: Posterior distribution of K — no true K highlighted
  # One figure with 3 subplots (E, M, D)
  # ----------------------------------------------------------------
  x11()
  par(mfrow = c(1, 3), mar = c(3, 3, 2, 1.4), mgp = c(1.75, 0.75, 0))
  
  # --- E-level ---
  plot(as.numeric(names(k_table_E)), as.numeric(k_table_E),
       type = "h", lwd = 3, col = "dodgerblue4",
       main = paste0("K_E posterior: Sim", s),
       xlab = "Número de clusters", ylab = "Densidad",
       yaxt = "n")
  axis(2, at = pretty(k_table_E), labels = pretty(k_table_E))
  
  # --- M-level ---
  plot(as.numeric(names(k_table_M)), as.numeric(k_table_M),
       type = "h", lwd = 3, col = "dodgerblue4",
       main = paste0("K_M posterior: Sim", s),
       xlab = "Número de clusters", ylab = "Densidad",
       yaxt = "n")
  axis(2, at = pretty(k_table_M), labels = pretty(k_table_M))
  
  # --- D-level ---
  plot(as.numeric(names(k_table_D)), as.numeric(k_table_D),
       type = "h", lwd = 3, col = "dodgerblue4",
       main = paste0("K_D posterior: Sim", s),
       xlab = "Número de clusters", ylab = "Densidad",
       yaxt = "n")
  axis(2, at = pretty(k_table_D), labels = pretty(k_table_D))
  
  # ----------------------------------------------------------------
  # Plot 2: Posterior distribution of K — true K highlighted in red
  # One figure with 3 subplots (E, M, D)
  # ----------------------------------------------------------------
  x11()
  par(mfrow = c(1, 3), mar = c(3, 3, 2, 1.4), mgp = c(1.75, 0.75, 0))
  
  # --- E-level ---
  x_vals_E <- as.numeric(names(k_table_E))
  y_vals_E <- as.numeric(k_table_E)
  cols_E   <- ifelse(x_vals_E == k_true_E, "firebrick", "dodgerblue4")
  lwds_E   <- ifelse(x_vals_E == k_true_E, 8, 3)
  
  plot(x_vals_E, y_vals_E, type = "h", lwd = lwds_E, col = cols_E,
       main = paste0("K_E posterior: Sim", s),
       xlab = "Número de clústers", ylab = "Densidad",
       yaxt = "n")
  axis(2, at = pretty(y_vals_E), labels = pretty(y_vals_E))
  
  # --- M-level ---
  x_vals_M <- as.numeric(names(k_table_M))
  y_vals_M <- as.numeric(k_table_M)
  cols_M   <- ifelse(x_vals_M == k_true_M, "firebrick", "dodgerblue4")
  lwds_M   <- ifelse(x_vals_M == k_true_M, 8, 3)
  
  plot(x_vals_M, y_vals_M, type = "h", lwd = lwds_M, col = cols_M,
       main = paste0("K_M posterior: Sim", s),
       xlab = "Número de clústers", ylab = "Densidad",
       yaxt = "n")
  axis(2, at = pretty(y_vals_M), labels = pretty(y_vals_M))
  
  # --- D-level ---
  x_vals_D <- as.numeric(names(k_table_D))
  y_vals_D <- as.numeric(k_table_D)
  cols_D   <- ifelse(x_vals_D == k_true_D, "firebrick", "dodgerblue4")
  lwds_D   <- ifelse(x_vals_D == k_true_D, 8, 3)
  
  plot(x_vals_D, y_vals_D, type = "h", lwd = lwds_D, col = cols_D,
       main = paste0("K_D posterior: Sim", s),
       xlab = "Número de clústers", ylab = "Densidad",
       yaxt = "n")
  axis(2, at = pretty(y_vals_D), labels = pretty(y_vals_D))
}


# ==============================================================================================
#@                        ADJUSTED RAND INDEX (ARI) — Per level E, M, D
# ==============================================================================================

# Storage: one ARI mean per simulation per level
ARI_mean_E <- numeric(length(cadenas))
ARI_mean_M <- numeric(length(cadenas))
ARI_mean_D <- numeric(length(cadenas))

for (idx in seq_along(cadenas)) {
  
  s <- cadenas[idx]
  
  current_chain   <- get(paste0("cadena_",         s))
  current_sampled <- get(paste0("simulated_data_", s))
  
  n_iter <- nrow(current_chain$xi_E)   # number of stored iterations
  
  ARI_E <- numeric(n_iter)
  ARI_M <- numeric(n_iter)
  ARI_D <- numeric(n_iter)
  
  # Compute the ARI per iteration for each level
  for (i in seq_len(n_iter)) {
    
    # E-level: compare true E-level assignments vs sampled E-level assignments at iteration i
    ARI_E[i] <- adjustedRandIndex(
      current_sampled$xi_E_true,
      as.vector(current_chain$xi_E[i, ])
    )
    
    # M-level: compare true M-level assignments vs sampled M-level assignments at iteration i
    ARI_M[i] <- adjustedRandIndex(
      current_sampled$xi_M_true,
      as.vector(current_chain$xi_M[i, ])
    )
    
    # D-level: compare true D-level assignments vs sampled D-level assignments at iteration i
    ARI_D[i] <- adjustedRandIndex(
      current_sampled$xi_D_true,
      as.vector(current_chain$xi_D[i, ])
    )
  }
  
  # Mean ARI across iterations for simulation s
  ARI_mean_E[idx] <- mean(ARI_E)
  ARI_mean_M[idx] <- mean(ARI_M)
  ARI_mean_D[idx] <- mean(ARI_D)
  
  cat("Simulation", s,
      "| ARI_E:", round(ARI_mean_E[idx], 4),
      "| ARI_M:", round(ARI_mean_M[idx], 4),
      "| ARI_D:", round(ARI_mean_D[idx], 4), "\n")
}

# Print results
cat("\n=== ARI Results per level ===\n")
names(ARI_mean_E) <- paste0("sim_", cadenas)
names(ARI_mean_M) <- paste0("sim_", cadenas)
names(ARI_mean_D) <- paste0("sim_", cadenas)
cat("E-level:\n"); print(round(ARI_mean_E, 4))
cat("M-level:\n"); print(round(ARI_mean_M, 4))
cat("D-level:\n"); print(round(ARI_mean_D, 4))



# ==============================================================================================
#@                        ADJUSTED RAND INDEX (ARI) — Joint E_M_D clusters
# ==============================================================================================

# Storage: one joint ARI mean per simulation
ARI_mean_joint <- numeric(length(cadenas))

for (idx in seq_along(cadenas)) {
  
  s <- cadenas[idx]
  
  current_chain   <- get(paste0("cadena_",         s))
  current_sampled <- get(paste0("simulated_data_", s))
  
  # ----------------------------------------------------------------
  # Build TRUE joint cluster assignment: xi_E_M_D_true
  # Each unique combination of (xi_E_true, xi_M_true, xi_D_true)
  # gets a consecutive integer label
  # ----------------------------------------------------------------
  joint_key_true <- paste(current_sampled$xi_E_true,
                          current_sampled$xi_M_true,
                          current_sampled$xi_D_true,
                          sep = "_")
  
  # Relabel unique combinations consecutively: 1, 2, 3, ...
  xi_E_M_D_true <- as.integer(factor(joint_key_true,
                                     levels = unique(joint_key_true)))
  
  n_iter    <- nrow(current_chain$xi_E)   # number of stored iterations
  ARI_joint <- numeric(n_iter)
  
  # Compute the ARI per iteration comparing true joint vs sampled joint
  for (i in seq_len(n_iter)) {
    
    # Build SAMPLED joint cluster assignment at iteration i:
    # combine the three sampled levels into one joint key
    joint_key_i <- paste(as.vector(current_chain$xi_E[i, ]),
                         as.vector(current_chain$xi_M[i, ]),
                         as.vector(current_chain$xi_D[i, ]),
                         sep = "_")
    
    # Relabel unique combinations consecutively: 1, 2, 3, ...
    xi_E_M_D_i <- as.integer(factor(joint_key_i,
                                    levels = unique(joint_key_i)))
    
    # ARI between true joint assignment and sampled joint assignment at iteration i
    ARI_joint[i] <- adjustedRandIndex(xi_E_M_D_true, xi_E_M_D_i)
  }
  
  # Mean ARI across iterations for simulation s
  ARI_mean_joint[idx] <- mean(ARI_joint)
  
  cat("Simulation", s, "| ARI joint E_M_D:", round(ARI_mean_joint[idx], 4), "\n")
}

# Print results
cat("\n=== ARI Results joint E_M_D ===\n")
names(ARI_mean_joint) <- paste0("sim_", cadenas)
print(round(ARI_mean_joint, 4))





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
  
  # Compute the number of clusters K for all levels at each iteration. 
  k_E_values <- apply(cadena$xi_E, 1, function(x) length(unique(x))) #note: The 1 applies for rows
  k_M_values <- apply(cadena$xi_M, 1, function(x) length(unique(x))) 
  k_D_values <- apply(cadena$xi_D, 1, function(x) length(unique(x))) #note: The 1 applies for rows
  
  # Compute the posterior distribution of K (proportion of each k)
  k_E_table  <- table(k_E_values) / length(k_E_values)
  k_M_table  <- table(k_M_values) / length(k_M_values)
  k_D_table  <- table(k_D_values) / length(k_D_values)
  
  # Find the mode of k (most repeated k) 
  k_E_mode   <- as.numeric(names(k_E_table)[which.max(k_E_table)])
  k_M_mode   <- as.numeric(names(k_M_table)[which.max(k_M_table)])
  k_D_mode   <- as.numeric(names(k_D_table)[which.max(k_D_table)])
  
  
  # Apply filter_cadena_by_k 
  # 1) Execute the function, then 2) store the filtered chain as cadena_1_filtered
  # Function filter_cadena_by_k is stored in funciones.R
  
  assign(
    filtered_chain_name,
    filter_cadena_by_k_3DP(cadena, k_E_mode, k_M_mode, k_D_mode)
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


for (s in cadenas) {
  
  cat("--- PERMUTATING SIMULATION", s, " --- \n")
  
  # Get the cadena_s_filtered
  filtered_chain_name    <- paste0("cadena_", s, "_filtered")      
  current_filtered_chain <- get(filtered_chain_name)
  
  # Get its corresponding simulated data
  simulated_data_name    <- paste0("simulated_data_",s)
  current_simulated_data <- get(simulated_data_name)
  
  
  result <- solve_label_switching_3DP_cpp(
    cadena_filtered = current_filtered_chain,                      # Filtered MCMC chain
    y               = current_simulated_data$y,                    # The real y data
    X_E_full        = current_simulated_data$X_E_full,             # The real X data
    Z_full          = current_simulated_data$Z_full,               # The real Z data
    W_full          = current_simulated_data$W_full,               # The real W data
    k_E_mode        = current_filtered_chain$filter_info$k_E_mode, # The mode of k_E clusters
    k_M_mode        = current_filtered_chain$filter_info$k_M_mode, # The mode of k_M clusters
    k_D_mode        = current_filtered_chain$filter_info$k_D_mode  # The mode of _D clusters
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
  
  result <- compute_posterior_summaries_3DP(
    cadena_filtered      = current_filtered_chain, 
    k_E_mode             = current_filtered_chain$filter_info$k_E_mode, 
    k_M_mode             = current_filtered_chain$filter_info$k_M_mode, 
    k_D_mode             = current_filtered_chain$filter_info$k_D_mode
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
  
  result <- compute_test_metrics_3DP(
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
#                         EXECUTE WAIC — MODEL 3
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
  y          <- current_simulated_data$y
  X_E_full   <- current_simulated_data$X_E_full
  Z_full     <- current_simulated_data$Z_full
  W_full     <- current_simulated_data$W_full
  obs_to_mun <- current_simulated_data$obs_to_mun
  mun_map    <- current_simulated_data$mun_map
  
  # Compute WAIC
  result_waic <- compute_WAIC_3DP(cadena_filtered = current_filtered_chain,
                                  y               = y,
                                  X_E_full        = X_E_full,
                                  Z_full          = Z_full,
                                  W_full          = W_full,
                                  obs_to_mun      = obs_to_mun,
                                  mun_map         = mun_map)
  
  # Store WAIC inside the filtered chain 
  current_filtered_chain$WAIC_summary <- result_waic
  
  # Assigns the waic to the cadena_s_filtered
  assign(filtered_chain_name, current_filtered_chain)
  
  # Pint execution
  cat("--- WAIC computation for simulation", s, "Completed --- \n")
  
}

































