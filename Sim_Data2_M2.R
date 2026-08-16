#@ simulations

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
ruta <- setwd('D:/Actualizado/Maestria estadistica/Tesis/Modelos BNP/23-04-26/Simulaciones/sim_data2_M2')

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
# Generate 21 simulated scenarios using the function generate_hierarchical_data.
# 
#-----------------------------------------------------------------------------

# =============================
# 1. VARIANCE SCENARIOS
# =============================

# In this variance scenarios we keep the number of departments moderated 5, the number 
#of municipalities per department and the number of observations per department homogeneous.
#Number of covariates is also moderated 12. Mixing probabilities homogeneous in the 5 clusters

# The only change along scenarios is the variance. We process three scenarios: low, moderated and
#high variance:

# The variance is relative to the cluster centers

# Low: 4
# Moderated: 10 
# High: 20


# Constant parameters
m_1 = 15                                         # Number of departments 
n_mun_per_dept_1 = rep(8, 15)                   # varying municipalities per dept
n_obs_range_1 = c(15, 20)                       # observations per municipality
p_d_1 = 4                                       # department-level covariates
p_m_1 = 4                                       # municipality-level covariates
p_i_1 = 4                                       # individual-level covariates
k_true_1 = 3                                    # # clusters
mix_probs_1 = c(1/3, 1/3, 1/3)                  # mixing probabilities
mu_centers_1 = c(10, 25, 35)              # well-separated cluster centers
sigma_k_1 = 1                                   # moderate within-cluster variation

# ----------------------------------
# 1. Low municipal variance scenario
# ----------------------------------
kappa_min_1_1 = 3                  # min value of municipal variance (runif)                    
kappa_max_1_1 = 5                  # max value of municipal variance (runif)


simulated_data_1 <- generate_hierarchical_data(
  m                = m_1,
  n_mun_per_dept   = n_mun_per_dept_1,        
  n_obs_range      = n_obs_range_1,               
  p_d              = p_d_1,                                
  p_m              = p_m_1,                                
  p_i              = p_i_1,                                
  k_true           = k_true_1,                          
  mix_probs        = mix_probs_1,          
  mu_centers       = mu_centers_1,              
  sigma_k          = sigma_k_1,                         
  kappa_min        = kappa_min_1_1,
  kappa_max        = kappa_max_1_1
)

# ----------------------------------
# 2. Moderated municipal variance scenario
# ----------------------------------
kappa_min_1_2 = 8               # min value of municipal variance (runif)                     
kappa_max_1_2 = 10              # max value of municipal variance (runif)

simulated_data_2 <- generate_hierarchical_data(
  m                = m_1,
  n_mun_per_dept   = n_mun_per_dept_1,        
  n_obs_range      = n_obs_range_1,               
  p_d              = p_d_1,                                
  p_m              = p_m_1,                                
  p_i              = p_i_1,                                
  k_true           = k_true_1,                          
  mix_probs        = mix_probs_1,          
  mu_centers       = mu_centers_1,              
  sigma_k          = sigma_k_1,                         
  kappa_min        = kappa_min_1_2,
  kappa_max        = kappa_max_1_2
)

# ----------------------------------
# 3. High municipal variance scenario
# ----------------------------------
kappa_min_1_3 = 14                         # min value of municipal variance (runif)              
kappa_max_1_3 = 20                        # max value of municipal variance (runif)

simulated_data_3 <- generate_hierarchical_data(
  m                = m_1,
  n_mun_per_dept   = n_mun_per_dept_1,        
  n_obs_range      = n_obs_range_1,               
  p_d              = p_d_1,                                
  p_m              = p_m_1,                                
  p_i              = p_i_1,                                
  k_true           = k_true_1,                          
  mix_probs        = mix_probs_1,          
  mu_centers       = mu_centers_1,              
  sigma_k          = sigma_k_1,                         
  kappa_min        = kappa_min_1_3,
  kappa_max        = kappa_max_1_3
)


# =============================
# 2. NUMBER OF CLUSTER K SCENARIOS
# =============================

# In this three scenarios we maintain number of departments moderated (5), homogeneous
#number of municipalities per department and observations per municipality. Also, 
#moderated number of covariates (12), and homogeneous mixing probabilities, 
#and low variance (8-10).

# The only change is the number of clusters:
# - Low clusters K = 4
# - Moderated    K = 10
# - High         K = 50

# Constant parameters
m_2 = 15                                         # Number of departments 
n_mun_per_dept_2 = rep(8, 15)             # varying municipalities per dept
n_obs_range_2 = c(15, 22)                       # observations per municipality
p_d_2 = 4                                       # department-level covariates
p_m_2 = 4                                       # municipality-level covariates
p_i_2 = 4                                       # individual-level covariates
sigma_k_2 = 1                                   # moderate within-cluster variation
kappa_min_2 = 2                                 # min value of municipal variance (runif)              
kappa_max_2 = 4                                # max value of municipal variance (runif)

# ----------------------------------
# 4. Low number of clusters k scenario
# ----------------------------------
k_true_2_1 = 2                                   # number of clusters
mix_probs_2_1 = c(1/2, 1/2)                  # mixing probabilities
mu_centers_2_1 = c(8, 19)                # well-separated cluster centers

simulated_data_4 <- generate_hierarchical_data(
  m                = m_2,
  n_mun_per_dept   = n_mun_per_dept_2,        
  n_obs_range      = n_obs_range_2,               
  p_d              = p_d_2,                                
  p_m              = p_m_2,                                
  p_i              = p_i_2,                                
  k_true           = k_true_2_1,                          
  mix_probs        = mix_probs_2_1,          
  mu_centers       = mu_centers_2_1,              
  sigma_k          = sigma_k_2,                         
  kappa_min        = kappa_min_2,
  kappa_max        = kappa_max_2
)


# ----------------------------------
# 5. Moderated number of clusters k scenario
# ----------------------------------
k_true_2_2 = 5                                                        # number of clusters
mix_probs_2_2 = rep(1/k_true_2_2, k_true_2_2)                                   # mixing probabilities
mu_centers_2_2 = seq(from = -12, by = 5, length.out = k_true_2_2)      # well-separated cluster centers

simulated_data_5 <- generate_hierarchical_data(
  m                = m_2,
  n_mun_per_dept   = n_mun_per_dept_2,        
  n_obs_range      = n_obs_range_2,               
  p_d              = p_d_2,                                
  p_m              = p_m_2,                                
  p_i              = p_i_2,                                
  k_true           = k_true_2_2,                          
  mix_probs        = mix_probs_2_2,          
  mu_centers       = mu_centers_2_2,              
  sigma_k          = sigma_k_2,                         
  kappa_min        = kappa_min_2,
  kappa_max        = kappa_max_2
)

# ----------------------------------
# 6. High number of clusters k scenario
# ----------------------------------
k_true_2_3 = 10                                                          # number of clusters
mix_probs_2_3 = rep(1/k_true_2_3, k_true_2_3)                            # mixing probabilities
mu_centers_2_3 = seq(from = -14, by = 3, length.out = k_true_2_3)        # well-separated cluster centers


simulated_data_6 <- generate_hierarchical_data(
  m                = m_2,
  n_mun_per_dept   = n_mun_per_dept_2,        
  n_obs_range      = n_obs_range_2,               
  p_d              = p_d_2,                                
  p_m              = p_m_2,                                
  p_i              = p_i_2,                                
  k_true           = k_true_2_3,                          
  mix_probs        = mix_probs_2_3,          
  mu_centers       = mu_centers_2_3,              
  sigma_k          = sigma_k_2,                         
  kappa_min        = kappa_min_2,
  kappa_max        = kappa_max_2
)


# =============================
# 3. NUMBER OF COVARIATES p
# =============================

# In this three scenarios we maintain number of departments moderated (5), 
#moderated number of clusters (5), and homogeneous number of municipalities 
#per department and observations per municipality. 
#Also homogeneous mixing probabilities, and low variance (8-10).


# The only change is the number of covariates:
# - Low          p = 6
# - Moderated    p = 45
# - High         p = 90


# Constant parameters
m_3 = 15                                         # Number of departments 
n_mun_per_dept_3 = rep(8, 15)                     # varying municipalities per dept
n_obs_range_3 = c(18, 22)                       # observations per municipality
k_true_3 = 3                                    # # clusters
mix_probs_3 = c(1/3, 1/3, 1/3)             # mixing probabilities
mu_centers_3 = c(-2, 4, 12)             # well-separated cluster centers
sigma_k_3 = 1                                   # moderate within-cluster variation
kappa_min_3 = 2                                 # min value of municipal variance (runif)              
kappa_max_3 = 4                                # max value of municipal variance (runif)

# ----------------------------------
# 7. Low number of covariates p scenario
# ----------------------------------
p_d_3_1 = 2                                       # department-level covariates
p_m_3_1 = 2                                       # municipality-level covariates
p_i_3_1 = 2                                       # individual-level covariates

simulated_data_7 <- generate_hierarchical_data(
  m                = m_3,
  n_mun_per_dept   = n_mun_per_dept_3,        
  n_obs_range      = n_obs_range_3,               
  p_d              = p_d_3_1,                                
  p_m              = p_m_3_1,                                
  p_i              = p_i_3_1,                                
  k_true           = k_true_3,                          
  mix_probs        = mix_probs_3,          
  mu_centers       = mu_centers_3,              
  sigma_k          = sigma_k_3,                         
  kappa_min        = kappa_min_3,
  kappa_max        = kappa_max_3
)

# ----------------------------------
# 8. Moderated number of covariates p scenario
# ----------------------------------
p_d_3_2 = 8                                       # department-level covariates
p_m_3_2 = 8                                      # municipality-level covariates
p_i_3_2 = 8                                       # individual-level covariates

simulated_data_8 <- generate_hierarchical_data(
  m                = m_3,
  n_mun_per_dept   = n_mun_per_dept_3,        
  n_obs_range      = n_obs_range_3,               
  p_d              = p_d_3_2,                                
  p_m              = p_m_3_2,                                
  p_i              = p_i_3_2,                                
  k_true           = k_true_3,                          
  mix_probs        = mix_probs_3,          
  mu_centers       = mu_centers_3,              
  sigma_k          = sigma_k_3,                         
  kappa_min        = kappa_min_3,
  kappa_max        = kappa_max_3
)

# ----------------------------------
# 9. High number of covariates p scenario
# ----------------------------------
p_d_3_3 = 14                                       # department-level covariates
p_m_3_3 = 14                                       # municipality-level covariates
p_i_3_3 = 14                                       # individual-level covariates

simulated_data_9 <- generate_hierarchical_data(
  m                = m_3,
  n_mun_per_dept   = n_mun_per_dept_3,        
  n_obs_range      = n_obs_range_3,               
  p_d              = p_d_3_3,                                
  p_m              = p_m_3_3,                                
  p_i              = p_i_3_3,                                
  k_true           = k_true_3,                          
  mix_probs        = mix_probs_3,          
  mu_centers       = mu_centers_3,              
  sigma_k          = sigma_k_3,                         
  kappa_min        = kappa_min_3,
  kappa_max        = kappa_max_3
)


# =============================
# 4. NUMBER OF DEPARTMENTS
# =============================

# In this scenario we keep a moderated number of covariates and homogeneous number of
#municipalities per department and observations per municipality. Also, keep a moderated 
#number of clusters (5) with homogeneous mixing probabilities. 

# The only change is the number of departments:
# - Low: 5
# - Moderated: 20
# - High: 40

# Constant parameters
n_obs_range_4 = c(18, 20)                       # observations per municipality
p_d_4 = 4                                       # department-level covariates
p_m_4 = 4                                       # municipality-level covariates
p_i_4 = 4                                       # individual-level covariates
k_true_4 = 3                                    # # clusters
mix_probs_4 = c(1/3, 1/3, 1/3)                  # mixing probabilities
mu_centers_4 = c(-5, 7, 12)             # well-separated cluster centers
sigma_k_4 = 1                                   # moderate within-cluster variation
kappa_min_4 = 2                                 # min value of municipal variance (runif)              
kappa_max_4 = 4                                # max value of municipal variance (runif)


# ----------------------------------
# 10. Low number of departments m scenario
# ----------------------------------
m_4_1 = 5                                         # Number of departments 
n_mun_per_dept_4_1 = c(8, 8, 8, 8, 8)             # varying municipalities per dept

simulated_data_10 <- generate_hierarchical_data(
  m                = m_4_1,
  n_mun_per_dept   = n_mun_per_dept_4_1,        
  n_obs_range      = n_obs_range_4,               
  p_d              = p_d_4,                                
  p_m              = p_m_4,                                
  p_i              = p_i_4,                                
  k_true           = k_true_4,                          
  mix_probs        = mix_probs_4,          
  mu_centers       = mu_centers_4,              
  sigma_k          = sigma_k_4,                         
  kappa_min        = kappa_min_4,
  kappa_max        = kappa_max_4
)


# ----------------------------------
# 11. Moderated number of departments m scenario
# ----------------------------------
m_4_2 = 25                                   # Number of departments 
n_mun_per_dept_4_2 = rep(7, m_4_2)             # varying municipalities per dept

simulated_data_11 <- generate_hierarchical_data(
  m                = m_4_2,
  n_mun_per_dept   = n_mun_per_dept_4_2,        
  n_obs_range      = n_obs_range_4,               
  p_d              = p_d_4,                                
  p_m              = p_m_4,                                
  p_i              = p_i_4,                                
  k_true           = k_true_4,                          
  mix_probs        = mix_probs_4,          
  mu_centers       = mu_centers_4,              
  sigma_k          = sigma_k_4,                         
  kappa_min        = kappa_min_4,
  kappa_max        = kappa_max_4
)

# ----------------------------------
# 12. High number of departments m scenario
# ----------------------------------
m_4_3 = 40                                     # Number of departments 
n_mun_per_dept_4_3 = rep(6, m_4_3)             # varying municipalities per dept

simulated_data_12 <- generate_hierarchical_data(
  m                = m_4_3,
  n_mun_per_dept   = n_mun_per_dept_4_3,        
  n_obs_range      = n_obs_range_4,               
  p_d              = p_d_4,                                
  p_m              = p_m_4,                                
  p_i              = p_i_4,                                
  k_true           = k_true_4,                          
  mix_probs        = mix_probs_4,          
  mu_centers       = mu_centers_4,              
  sigma_k          = sigma_k_4,                         
  kappa_min        = kappa_min_4,
  kappa_max        = kappa_max_4
)

# =============================
# 5. NUMBER OF MUNICIPALITIES
# =============================

# In this scenario we keep a moderated number departments (5) and of covariates
#we keep an homogeneous number observations per municipality. Also, keep a moderated 
#number of clusters (5) with homogeneous mixing probabilities. 

# The only change is the number of Municipalities per department:
# - Low: 5
# - Moderated: 25
# - High: 50

# Constant parameters
m_5 = 15                                         # Number of departments 
n_obs_range_5 = c(18, 22)                       # observations per municipality
p_d_5 = 4                                       # department-level covariates
p_m_5 = 4                                       # municipality-level covariates
p_i_5 = 4                                       # individual-level covariates
k_true_5 = 5                                    # # clusters
mix_probs_5 = c(0.2, 0.2, 0.2, 0.2, 0.2)        # mixing probabilities
mu_centers_5 = c(-25, -8, 0, 5, 15)              # well-separated cluster centers
sigma_k_5 = 1                                   # moderate within-cluster variation
kappa_min_5 = 3                                 # min value of municipal variance (runif)              
kappa_max_5 = 5                                # max value of municipal variance (runif)

# -------------------------------------------------
# 13. Low number of municipalities j per department 
# -------------------------------------------------
n_mun_per_dept_5_1 = rep(3, m_5)             # varying municipalities per dept

simulated_data_13 <- generate_hierarchical_data(
  m                = m_5,
  n_mun_per_dept   = n_mun_per_dept_5_1,        
  n_obs_range      = n_obs_range_5,               
  p_d              = p_d_5,                                
  p_m              = p_m_5,                                
  p_i              = p_i_5,                                
  k_true           = k_true_5,                          
  mix_probs        = mix_probs_5,          
  mu_centers       = mu_centers_5,              
  sigma_k          = sigma_k_5,                         
  kappa_min        = kappa_min_5,
  kappa_max        = kappa_max_5
)


# -------------------------------------------------
# 14. Moderated number of municipalities j per department 
# -------------------------------------------------
n_mun_per_dept_5_2 = rep(10, m_5)             # varying municipalities per dept

simulated_data_14 <- generate_hierarchical_data(
  m                = m_5,
  n_mun_per_dept   = n_mun_per_dept_5_2,        
  n_obs_range      = n_obs_range_5,               
  p_d              = p_d_5,                                
  p_m              = p_m_5,                                
  p_i              = p_i_5,                                
  k_true           = k_true_5,                          
  mix_probs        = mix_probs_5,          
  mu_centers       = mu_centers_5,              
  sigma_k          = sigma_k_5,                         
  kappa_min        = kappa_min_5,
  kappa_max        = kappa_max_5
)


# -------------------------------------------------
# 15. High number of municipalities j per department 
# -------------------------------------------------
n_mun_per_dept_5_3 = rep(24, m_5)             # varying municipalities per dept

simulated_data_15 <- generate_hierarchical_data(
  m                = m_5,
  n_mun_per_dept   = n_mun_per_dept_5_3,        
  n_obs_range      = n_obs_range_5,               
  p_d              = p_d_5,                                
  p_m              = p_m_5,                                
  p_i              = p_i_5,                                
  k_true           = k_true_5,                          
  mix_probs        = mix_probs_5,          
  mu_centers       = mu_centers_5,              
  sigma_k          = sigma_k_5,                         
  kappa_min        = kappa_min_5,
  kappa_max        = kappa_max_5
)

# =============================
# 6. MIXING PROBABILITIES
# =============================

# In this scenario we keep a moderated number departments (5) and of covariates
#we keep an homogeneous number of municipalities per department, and number observations per municipality
#. Also, keep a moderated number of clusters (6). 

# The only change is the mixing probabilities:
# - 1 cluster with high probability (0.8, 0.04, 0.04, 0.04, 0.04, 0.04)
# - 1 cluster with low probability (0.02, 0.196, 0.196, 0.196, 0.196, 0.196)
# - 2 clusters with high probability, 2 with low (0.34, 0.34, 0.14, 0.14, 0.02, 0.02)


# Constant parameters
m_6 = 15                                         # Number of departments 
n_mun_per_dept_6 = rep(10, 15)                  # varying municipalities per dept
n_obs_range_6 = c(15, 20)                       # observations per municipality
p_d_6 = 4                                       # department-level covariates
p_m_6 = 4                                       # municipality-level covariates
p_i_6 = 4                                       # individual-level covariates
sigma_k_6 = 1                                   # moderate within-cluster variation
kappa_min_6 = 2                                 # min value of municipal variance (runif)              
kappa_max_6 = 4                                # max value of municipal variance (runif)



# -------------------------------------------------
# 16. One cluster with high probability
# -------------------------------------------------
k_true_6_1 = 3                                                        # number clusters
mix_probs_6_1 = c(0.45, 0.45, 0.1)                                    # mixing probabilities
mu_centers_6_1 = seq(from = -3, by = 9, length.out = k_true_6_1)       # well-separated cluster centers


simulated_data_16 <- generate_hierarchical_data(
  m                = m_6,
  n_mun_per_dept   = n_mun_per_dept_6,        
  n_obs_range      = n_obs_range_6,               
  p_d              = p_d_6,                                
  p_m              = p_m_6,                                
  p_i              = p_i_6,                                
  k_true           = k_true_6_1,                          
  mix_probs        = mix_probs_6_1,          
  mu_centers       = mu_centers_6_1,              
  sigma_k          = sigma_k_6,                         
  kappa_min        = kappa_min_6,
  kappa_max        = kappa_max_6
)

# -------------------------------------------------
# 17. One cluster with low probability
# -------------------------------------------------
k_true_6_2 = 3                                                          # number clusters
mix_probs_6_2 = c(0.8, 0.1, 0.1)                                        # mixing probabilities
mu_centers_6_2 = seq(from = -3, by = 9, length.out = k_true_6_2)       # well-separated cluster centers


simulated_data_17 <- generate_hierarchical_data(
  m                = m_6,
  n_mun_per_dept   = n_mun_per_dept_6,        
  n_obs_range      = n_obs_range_6,               
  p_d              = p_d_6,                                
  p_m              = p_m_6,                                
  p_i              = p_i_6,                                
  k_true           = k_true_6_2,                          
  mix_probs        = mix_probs_6_2,          
  mu_centers       = mu_centers_6_2,              
  sigma_k          = sigma_k_6,                         
  kappa_min        = kappa_min_6,
  kappa_max        = kappa_max_6
)

# -------------------------------------------------
# 18.Two clusters with high and two with low probability
# -------------------------------------------------
k_true_6_3 = 3                                                          # number clusters
mix_probs_6_3 = c(0.4, 0.25, 0.35)                                      # mixing probabilities
mu_centers_6_3 = seq(from = 2, by = 11, length.out = k_true_6_3)        # well-separated cluster centers


simulated_data_18 <- generate_hierarchical_data(
  m                = m_6,
  n_mun_per_dept   = n_mun_per_dept_6,        
  n_obs_range      = n_obs_range_6,               
  p_d              = p_d_6,                                
  p_m              = p_m_6,                                
  p_i              = p_i_6,                                
  k_true           = k_true_6_3,                          
  mix_probs        = mix_probs_6_3,          
  mu_centers       = mu_centers_6_3,              
  sigma_k          = sigma_k_6,                         
  kappa_min        = kappa_min_6,
  kappa_max        = kappa_max_6
)

# =============================
# 7. Number of observations n
# =============================

# In this scenario we keep a moderated number departments (5) and of covariates
#we keep an homogeneous number of municipalities per department, and number observations per municipality
#. Also, keep a moderated number of clusters (6). Mixing probabilities are homogeneous across clusters

# The only change is the number of observations per municipality n_obs_range:
# - Low: 5-8
# - Moderated: 30-40
# - High: 140-160


# Constant parameters
# Constant parameters
m_7 = 15                                        # Number of departments 
n_mun_per_dept_7 = rep(8, 15)                   # varying municipalities per dept
p_d_7 = 4                                       # department-level covariates
p_m_7 = 4                                       # municipality-level covariates
p_i_7 = 4                                       # individual-level covariates
k_true_7 = 3                                    # # clusters
mix_probs_7 = c(1/3, 1/3, 1/3)                  # mixing probabilities
mu_centers_7 = c(2, 15, 34)                     # well-separated cluster centers
sigma_k_7 = 1                                   # moderate within-cluster variation
kappa_min_7 = 2                                 # min value of municipal variance (runif)              
kappa_max_7 = 4                                 # max value of municipal variance (runif)

# -------------------------------------------------
# 19.Low number of observations n
# -------------------------------------------------
n_obs_range_7_1 = c(5, 8)                       # observations per municipality


simulated_data_19 <- generate_hierarchical_data(
  m                = m_7,
  n_mun_per_dept   = n_mun_per_dept_7,        
  n_obs_range      = n_obs_range_7_1,               
  p_d              = p_d_7,                                
  p_m              = p_m_7,                                
  p_i              = p_i_7,                                
  k_true           = k_true_7,                          
  mix_probs        = mix_probs_7,          
  mu_centers       = mu_centers_7,              
  sigma_k          = sigma_k_7,                         
  kappa_min        = kappa_min_7,
  kappa_max        = kappa_max_7
)

# -------------------------------------------------
# 20.Moderated number of observations n
# -------------------------------------------------
n_obs_range_7_2 = c(30, 40)                       # observations per municipality


simulated_data_20 <- generate_hierarchical_data(
  m                = m_7,
  n_mun_per_dept   = n_mun_per_dept_7,        
  n_obs_range      = n_obs_range_7_2,               
  p_d              = p_d_7,                                
  p_m              = p_m_7,                                
  p_i              = p_i_7,                                
  k_true           = k_true_7,                          
  mix_probs        = mix_probs_7,          
  mu_centers       = mu_centers_7,              
  sigma_k          = sigma_k_7,                         
  kappa_min        = kappa_min_7,
  kappa_max        = kappa_max_7
)


# -------------------------------------------------
# 21. High number of observations n
# -------------------------------------------------
n_obs_range_7_3 = c(80, 100)                       # observations per municipality


simulated_data_21 <- generate_hierarchical_data(
  m                = m_7,
  n_mun_per_dept   = n_mun_per_dept_7,        
  n_obs_range      = n_obs_range_7_3,               
  p_d              = p_d_7,                                
  p_m              = p_m_7,                                
  p_i              = p_i_7,                                
  k_true           = k_true_7,                          
  mix_probs        = mix_probs_7,          
  mu_centers       = mu_centers_7,              
  sigma_k          = sigma_k_7,                         
  kappa_min        = kappa_min_7,
  kappa_max        = kappa_max_7
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
# Training: n_obs_range = c(15, 20) 

# ----------------------------------
# Test 1. Low municipal variance
# ----------------------------------
test_data_1 <- generate_hierarchical_data(
  m                = m_1,
  n_mun_per_dept   = n_mun_per_dept_1,
  n_obs_range      = ceiling(n_obs_range_1 * 0.3),              # 30% aprox
  p_d              = p_d_1,
  p_m              = p_m_1,
  p_i              = p_i_1,
  k_true           = k_true_1,
  mix_probs        = mix_probs_1,
  mu_centers       = mu_centers_1,
  sigma_k          = sigma_k_1,
  kappa_min        = kappa_min_1_1,
  kappa_max        = kappa_max_1_1
)

# ----------------------------------
# Test 2. Moderated municipal variance
# ----------------------------------
test_data_2 <- generate_hierarchical_data(
  m                = m_1,
  n_mun_per_dept   = n_mun_per_dept_1,
  n_obs_range      = ceiling(n_obs_range_1 * 0.3),              # 30% aprox
  p_d              = p_d_1,
  p_m              = p_m_1,
  p_i              = p_i_1,
  k_true           = k_true_1,
  mix_probs        = mix_probs_1,
  mu_centers       = mu_centers_1,
  sigma_k          = sigma_k_1,
  kappa_min        = kappa_min_1_2,
  kappa_max        = kappa_max_1_2
)

# ----------------------------------
# Test 3. High municipal variance
# ----------------------------------
test_data_3 <- generate_hierarchical_data(
  m                = m_1,
  n_mun_per_dept   = n_mun_per_dept_1,
  n_obs_range      = ceiling(n_obs_range_1 * 0.3),              # 30% aprox
  p_d              = p_d_1,
  p_m              = p_m_1,
  p_i              = p_i_1,
  k_true           = k_true_1,
  mix_probs        = mix_probs_1,
  mu_centers       = mu_centers_1,
  sigma_k          = sigma_k_1,
  kappa_min        = kappa_min_1_3,
  kappa_max        = kappa_max_1_3
)


# =============================
# 2. NUMBER OF CLUSTERS K SCENARIOS
# =============================
# Training: n_obs_range = c(20, 25) → mean 22.5 → ~900 obs (40 mun)
# Test:     n_obs_range = c(6,  8)  → mean  7.0 → ~280 obs (~31%)

# ----------------------------------
# Test 4. Low number of clusters k
# ----------------------------------
test_data_4 <- generate_hierarchical_data(
  m                = m_2,
  n_mun_per_dept   = n_mun_per_dept_2,
  n_obs_range      = ceiling(n_obs_range_2 * 0.3),              # 30% aprox
  p_d              = p_d_2,
  p_m              = p_m_2,
  p_i              = p_i_2,
  k_true           = k_true_2_1,
  mix_probs        = mix_probs_2_1,
  mu_centers       = mu_centers_2_1,
  sigma_k          = sigma_k_2,
  kappa_min        = kappa_min_2,
  kappa_max        = kappa_max_2
)

# ----------------------------------
# Test 5. Moderated number of clusters k
# ----------------------------------
test_data_5 <- generate_hierarchical_data(
  m                = m_2,
  n_mun_per_dept   = n_mun_per_dept_2,
  n_obs_range      = ceiling(n_obs_range_2 * 0.3),              # 30% aprox
  p_d              = p_d_2,
  p_m              = p_m_2,
  p_i              = p_i_2,
  k_true           = k_true_2_2,
  mix_probs        = mix_probs_2_2,
  mu_centers       = mu_centers_2_2,
  sigma_k          = sigma_k_2,
  kappa_min        = kappa_min_2,
  kappa_max        = kappa_max_2
)

# ----------------------------------
# Test 6. High number of clusters k
# ----------------------------------
test_data_6 <- generate_hierarchical_data(
  m                = m_2,
  n_mun_per_dept   = n_mun_per_dept_2,
  n_obs_range      = ceiling(n_obs_range_2 * 0.3),              # 30% aprox
  p_d              = p_d_2,
  p_m              = p_m_2,
  p_i              = p_i_2,
  k_true           = k_true_2_3,
  mix_probs        = mix_probs_2_3,
  mu_centers       = mu_centers_2_3,
  sigma_k          = sigma_k_2,
  kappa_min        = kappa_min_2,
  kappa_max        = kappa_max_2
)


# =============================
# 3. NUMBER OF COVARIATES p SCENARIOS
# =============================
# Training: n_obs_range = c(17, 22) 

# ----------------------------------
# Test 7. Low number of covariates p
# ----------------------------------
test_data_7 <- generate_hierarchical_data(
  m                = m_3,
  n_mun_per_dept   = n_mun_per_dept_3,
  n_obs_range      = ceiling(n_obs_range_3 * 0.3),              # 30% aprox
  p_d              = p_d_3_1,
  p_m              = p_m_3_1,
  p_i              = p_i_3_1,
  k_true           = k_true_3,
  mix_probs        = mix_probs_3,
  mu_centers       = mu_centers_3,
  sigma_k          = sigma_k_3,
  kappa_min        = kappa_min_3,
  kappa_max        = kappa_max_3
)

# ----------------------------------
# Test 8. Moderated number of covariates p
# ----------------------------------
test_data_8 <- generate_hierarchical_data(
  m                = m_3,
  n_mun_per_dept   = n_mun_per_dept_3,
  n_obs_range      = ceiling(n_obs_range_3 * 0.3),              # 30% aprox
  p_d              = p_d_3_2,
  p_m              = p_m_3_2,
  p_i              = p_i_3_2,
  k_true           = k_true_3,
  mix_probs        = mix_probs_3,
  mu_centers       = mu_centers_3,
  sigma_k          = sigma_k_3,
  kappa_min        = kappa_min_3,
  kappa_max        = kappa_max_3
)

# ----------------------------------
# Test 9. High number of covariates p
# ----------------------------------
test_data_9 <- generate_hierarchical_data(
  m                = m_3,
  n_mun_per_dept   = n_mun_per_dept_3,
  n_obs_range      = ceiling(n_obs_range_3 * 0.3),              # 30% aprox
  p_d              = p_d_3_3,
  p_m              = p_m_3_3,
  p_i              = p_i_3_3,
  k_true           = k_true_3,
  mix_probs        = mix_probs_3,
  mu_centers       = mu_centers_3,
  sigma_k          = sigma_k_3,
  kappa_min        = kappa_min_3,
  kappa_max        = kappa_max_3
)


# =============================
# 4. NUMBER OF DEPARTMENTS SCENARIOS
# =============================
# Training: n_obs_range = c(15, 20) → mean 17.5


# ----------------------------------
# Test 10. Low number of departments m
# ----------------------------------
test_data_10 <- generate_hierarchical_data(
  m                = m_4_1,
  n_mun_per_dept   = n_mun_per_dept_4_1,
  n_obs_range      = ceiling(n_obs_range_4 * 0.3),              # 30% aprox
  p_d              = p_d_4,
  p_m              = p_m_4,
  p_i              = p_i_4,
  k_true           = k_true_4,
  mix_probs        = mix_probs_4,
  mu_centers       = mu_centers_4,
  sigma_k          = sigma_k_4,
  kappa_min        = kappa_min_4,
  kappa_max        = kappa_max_4
)

# ----------------------------------
# Test 11. Moderated number of departments m
# ----------------------------------
test_data_11 <- generate_hierarchical_data(
  m                = m_4_2,
  n_mun_per_dept   = n_mun_per_dept_4_2,
  n_obs_range      = ceiling(n_obs_range_4 * 0.3),              # 30% aprox
  p_d              = p_d_4,
  p_m              = p_m_4,
  p_i              = p_i_4,
  k_true           = k_true_4,
  mix_probs        = mix_probs_4,
  mu_centers       = mu_centers_4,
  sigma_k          = sigma_k_4,
  kappa_min        = kappa_min_4,
  kappa_max        = kappa_max_4
)

# ----------------------------------
# Test 12. High number of departments m
# ----------------------------------
test_data_12 <- generate_hierarchical_data(
  m                = m_4_3,
  n_mun_per_dept   = n_mun_per_dept_4_3,
  n_obs_range      = ceiling(n_obs_range_4 * 0.3),              # 30% aprox
  p_d              = p_d_4,
  p_m              = p_m_4,
  p_i              = p_i_4,
  k_true           = k_true_4,
  mix_probs        = mix_probs_4,
  mu_centers       = mu_centers_4,
  sigma_k          = sigma_k_4,
  kappa_min        = kappa_min_4,
  kappa_max        = kappa_max_4
)


# =============================
# 5. NUMBER OF MUNICIPALITIES SCENARIOS
# =============================
# Training: n_obs_range = c(15, 20) → mean 17.5


# ----------------------------------
# Test 13. Low municipalities per department
# ----------------------------------
test_data_13 <- generate_hierarchical_data(
  m                = m_5,
  n_mun_per_dept   = n_mun_per_dept_5_1,
  n_obs_range      = ceiling(n_obs_range_5 * 0.3),              # 30% aprox
  p_d              = p_d_5,
  p_m              = p_m_5,
  p_i              = p_i_5,
  k_true           = k_true_5,
  mix_probs        = mix_probs_5,
  mu_centers       = mu_centers_5,
  sigma_k          = sigma_k_5,
  kappa_min        = kappa_min_5,
  kappa_max        = kappa_max_5
)

# ----------------------------------
# Test 14. Moderated municipalities per department
# ----------------------------------
test_data_14 <- generate_hierarchical_data(
  m                = m_5,
  n_mun_per_dept   = n_mun_per_dept_5_2,
  n_obs_range      = ceiling(n_obs_range_5 * 0.3),              # 30% aprox
  p_d              = p_d_5,
  p_m              = p_m_5,
  p_i              = p_i_5,
  k_true           = k_true_5,
  mix_probs        = mix_probs_5,
  mu_centers       = mu_centers_5,
  sigma_k          = sigma_k_5,
  kappa_min        = kappa_min_5,
  kappa_max        = kappa_max_5
)

# ----------------------------------
# Test 15. High municipalities per department
# ----------------------------------
test_data_15 <- generate_hierarchical_data(
  m                = m_5,
  n_mun_per_dept   = n_mun_per_dept_5_3,
  n_obs_range      = ceiling(n_obs_range_5 * 0.3),              # 30% aprox
  p_d              = p_d_5,
  p_m              = p_m_5,
  p_i              = p_i_5,
  k_true           = k_true_5,
  mix_probs        = mix_probs_5,
  mu_centers       = mu_centers_5,
  sigma_k          = sigma_k_5,
  kappa_min        = kappa_min_5,
  kappa_max        = kappa_max_5
)


# =============================
# 6. MIXING PROBABILITIES SCENARIOS
# =============================
# Training: n_obs_range = c(15, 20)

# ----------------------------------
# Test 16. One cluster with high probability
# ----------------------------------
test_data_16 <- generate_hierarchical_data(
  m                = m_6,
  n_mun_per_dept   = n_mun_per_dept_6,
  n_obs_range      = ceiling(n_obs_range_6 * 0.3),              # 30% aprox
  p_d              = p_d_6,
  p_m              = p_m_6,
  p_i              = p_i_6,
  k_true           = k_true_6_1,
  mix_probs        = mix_probs_6_1,
  mu_centers       = mu_centers_6_1,
  sigma_k          = sigma_k_6,
  kappa_min        = kappa_min_6,
  kappa_max        = kappa_max_6
)

# ----------------------------------
# Test 17. One cluster with low probability
# ----------------------------------
test_data_17 <- generate_hierarchical_data(
  m                = m_6,
  n_mun_per_dept   = n_mun_per_dept_6,
  n_obs_range      = ceiling(n_obs_range_6 * 0.3),              # 30% aprox
  p_d              = p_d_6,
  p_m              = p_m_6,
  p_i              = p_i_6,
  k_true           = k_true_6_2,
  mix_probs        = mix_probs_6_2,
  mu_centers       = mu_centers_6_2,
  sigma_k          = sigma_k_6,
  kappa_min        = kappa_min_6,
  kappa_max        = kappa_max_6
)

# ----------------------------------
# Test 18. Two clusters high, two clusters low probability
# ----------------------------------
test_data_18 <- generate_hierarchical_data(
  m                = m_6,
  n_mun_per_dept   = n_mun_per_dept_6,
  n_obs_range      = ceiling(n_obs_range_6 * 0.3),              # 30% aprox
  p_d              = p_d_6,
  p_m              = p_m_6,
  p_i              = p_i_6,
  k_true           = k_true_6_3,
  mix_probs        = mix_probs_6_3,
  mu_centers       = mu_centers_6_3,
  sigma_k          = sigma_k_6,
  kappa_min        = kappa_min_6,
  kappa_max        = kappa_max_6
)


# =============================
# 7. NUMBER OF OBSERVATIONS n SCENARIOS
# =============================
# Here n_obs_range itself is the varying factor, so we scale each
# training range independently by 0.3.
#
# Test 19: c(5,8)     × 0.3 
# Test 20: c(30,40)   × 0.3 
# Test 21: c(140,160) × 0.3 

# ----------------------------------
# Test 19. Low number of observations n
# ----------------------------------
test_data_19 <- generate_hierarchical_data(
  m                = m_7,
  n_mun_per_dept   = n_mun_per_dept_7,
  n_obs_range      = ceiling(n_obs_range_7_1 * 0.3),              # 30% aprox
  p_d              = p_d_7,
  p_m              = p_m_7,
  p_i              = p_i_7,
  k_true           = k_true_7,
  mix_probs        = mix_probs_7,
  mu_centers       = mu_centers_7,
  sigma_k          = sigma_k_7,
  kappa_min        = kappa_min_7,
  kappa_max        = kappa_max_7
)

# ----------------------------------
# Test 20. Moderated number of observations n
# ----------------------------------
test_data_20 <- generate_hierarchical_data(
  m                = m_7,
  n_mun_per_dept   = n_mun_per_dept_7,
  n_obs_range      = ceiling(n_obs_range_7_2 * 0.3),              # 30% aprox
  p_d              = p_d_7,
  p_m              = p_m_7,
  p_i              = p_i_7,
  k_true           = k_true_7,
  mix_probs        = mix_probs_7,
  mu_centers       = mu_centers_7,
  sigma_k          = sigma_k_7,
  kappa_min        = kappa_min_7,
  kappa_max        = kappa_max_7
)

# ----------------------------------
# Test 21. High number of observations n
# ----------------------------------
test_data_21 <- generate_hierarchical_data(
  m                = m_7,
  n_mun_per_dept   = n_mun_per_dept_7,
  n_obs_range      = ceiling(n_obs_range_7_3 * 0.3),              # 30% aprox
  p_d              = p_d_7,
  p_m              = p_m_7,
  p_i              = p_i_7,
  k_true           = k_true_7,
  mix_probs        = mix_probs_7,
  mu_centers       = mu_centers_7,
  sigma_k          = sigma_k_7,
  kappa_min        = kappa_min_7,
  kappa_max        = kappa_max_7
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





#nu_beta         = 0.001
#gamma_beta      = 0.001
#nu_kappa        = 0.001
#alpha_kappa     = 0.001
#beta_kappa      = 0.001
#sigma2_base     = 5
#alpha           = 1
#a_alpha         = 0.001
#b_alpha         = 0.001
#nu2_mu_inv      = 1 / 10
#a_sigma2        = 0.0001
#b_sigma2        = 0.0001
#a_b_sigma2      = a_sigma2*b_sigma2
#sigma2_beta     = 5


# ===================================================================
#@                       GIBBS EXECUTION
# ===================================================================


#------------------------------------------------
#      ITERATE OVER SIMULATED DATAFRAMES
#_-----------------------------------------------

for (s in 1:21) {

  simulated_df_name <- sprintf("simulated_data_%d", s)
  obj <- get(simulated_df_name, envir = .GlobalEnv)
  
  #------------------------------------------------
  #      Retrieve information from the database
  #------------------------------------------------
  
  #Dataframe extraction
  df               <- obj$datos
  x_names          <- obj$x_names
  Y                <- obj$Y
  y                <- obj$y
  m                <- obj$m
  X_full           <- obj$X_full
  mun_kappa        <- obj$mun_kappa
  sqrt_kappa       <- obj$sqrt_kappa
  xTx_vec          <- obj$xTx_vec
  xi               <- obj$xi_false       #xi to initialize (Not the real xi cluster assignments)
  beta_k           <- obj$beta_k_false   #beta_k to initialize (Not the real beta_k atoms)
  beta_int         <- obj$beta_int_false
  
  #mu_vec           <- unlist(obj$beta_k_false[1])
  #eta_mu           <- unlist(obj$beta_k_false[2])
  p                <- obj$p
  obs_to_mun       <- obj$obs_to_mun
  mun_map          <- obj$mun_map
  kappa2_q         <- obj$kappa2_q
  dept_to_mun      <- obj$dept_to_mun
  n_mun_total      <- obj$n_mun_total
  
  mu_beta          <- beta_int
  
  # Unitary priors hiperparameter specification
  lm_previa_int       <-  lm(y ~ 1)
  lm_previa_unitaria  <-  lm(y ~ X_full -1)
  
  sigma2_beta         <-  summary(lm_previa_int)$sigma
  gamma_beta          <-  summary(lm_previa_int)$sigma
  b_sigma2            <-  summary(lm_previa_int)$sigma
  
  a_b_sigma2      <- a_sigma2*b_sigma2
  
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

chains_executed <- 21


for (s in 1:chains_executed) {
  
  results_file_name <- file.path(path_to_results, paste0("Chain_Simulation_", s, ".RData"))
  
  # Load into a temporary environment to avoid overwriting
  temp_env <- new.env()
  load(file = results_file_name, envir = temp_env)
  
  # Assign each chain to a uniquely named object
  assign(paste0("cadena_", s), temp_env$cadena)
  
  cat("Loaded:", results_file_name, "\n")
}


##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

# ==============================================================================================
#@                                     PLOTS
# ==============================================================================================

cadenas <- (1:21)

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
  
  
  
  # Option 2: with the real number of clusters
  # Convert once for convenience
  x_vals <- as.numeric(names(k_table))
  y_vals <- as.numeric(k_table)
  ymax  <- max(y_vals)
  
  x11()
  par(mfrow = c(1,1), mar = c(3,3,1.4,1.4), mgp = c(1.75,0.75,0))
  
  cols <- ifelse(x_vals == k_true, "firebrick", "dodgerblue4")
  lwds <- ifelse(x_vals == k_true, 8, 3)
  
  plot(x_vals, y_vals, type = "h", lwd = lwds, col = cols,
       main = paste0("Distribución posterior número clústers: Sim", s),
       xlab = "Número de clústers", ylab = "Densidad",
       yaxt = "n")
  
  axis(2, at = pretty(y_vals), labels = pretty(y_vals))
  
}




# ================================================================
#@                 CO-CLUSTERING MATRIX 
# ================================================================

# Compute co-clustering matrix A
compute_co_clustering_matrix <- function(chain) {
  n <- ncol(chain$xi)
  B <- nrow(chain$xi)
  A <- matrix(0, nrow = n, ncol = n)
  
  for (b in 1:B) {
    xi_b <- chain$xi[b, ]
    for (k in sort(unique(xi_b))) {
      cluster_members <- which(xi_b == k)
      A[cluster_members, cluster_members] <- A[cluster_members, cluster_members] + 1
    }
  }
  
  A <- A / B
  diag(A) <- 1
  
  return(A)
}

# Plot the co_clustering matrix for each chain indepently

# Compute and reorder A based on true partition
A <- compute_co_clustering_matrix(cadena_3)

# Convert A to dissimilarity
D <- 1 - A

# Apply classical multidimensional scaling
mds_coords <- cmdscale(as.dist(D), k = 4)

# Run Mclust on the embedded coordinates
clustering_result <- Mclust(mds_coords, G = 8, verbose = FALSE)

# Optimal partition
components <- clustering_result$classification

# Sort observations by true partition
A <- A[order(components), order(components)]

x11()
# Visualize co-clustering matrix
par(mar = c(2.75, 2.75, 0.5, 0.5), mgp = c(1.7, 0.7, 0))
colorscale <- c("white", rev(heat.colors(100)))
corrplot::corrplot(A, is.corr = FALSE, method = "color", col = colorscale,
                   tl.pos = "n", addgrid.col = NA)



# ==============================================================================================
#@                             ADJUSTED RAND INDEX (ARI)
# ==============================================================================================

# Compute the ARI for all simulations (all iterations within a simulation)
ARI_mean <- numeric(length(cadenas))

for (s in cadenas) {
  
  current_chain <- get(paste0("cadena_", s))
  
  current_sampled <- get(paste0("simulated_data_", s))
  
  n_iter <- nrow(current_chain$xi)
  
  ARI <- numeric(n_iter)
  
  # Compute the ARI per iteration
  for (i in seq(n_iter)) {
    ARI[i] <- adjustedRandIndex(
      current_sampled$xi_true, 
      as.vector(current_chain$xi[i, ])
    )
  }
  # Mean of ARI iterations
  ARI_mean[s] <- mean(ARI)
}
ARI_mean
round(ARI_mean,2)



# ==============================================================================================
#@                                    INFERENCE
# ==============================================================================================

cadenas <- c(1:21)

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


for (s in cadenas) {
  
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
























