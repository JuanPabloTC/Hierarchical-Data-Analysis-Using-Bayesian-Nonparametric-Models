#@ Modelo 3

# ===================================================================
#                             PREPARATION
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
suppressMessages(suppressWarnings(library(mclust)))

# Set working directory 
ruta <- setwd('D:/Actualizado/Maestria estadistica/Tesis/Modelos BNP/23-05-26')

# Path for results
path_to_results <- paste0(ruta,'/Resultados/Modelo_3/MCMC3')
# Path for text files
path_text_files <- paste0(ruta,'/Resultados/Modelo_3/text_files')


# Create dir (folder) if does not exist for text files
if (!dir.exists(path_text_files)) dir.create(path_text_files, recursive = TRUE, showWarnings = FALSE)

# Create dir (folder) if does not exist for MCMC2
if (!dir.exists(path_to_results)) dir.create(path_to_results, recursive = TRUE, showWarnings = FALSE)

# Load training data 
datos <- read.csv(paste0(ruta,'/Datos_Procesados/Examen_Saber_11_2022_2_Datos_train.csv'))

# Load testing data 
datos_testeo <- read.csv(paste0(ruta,'/Datos_Procesados/Examen_Saber_11_2022_2_Datos_test.csv'))


# ==============================================================================================
#                                  DATA MANIPULATION
# ==============================================================================================


cat("\n=== DATA ===\n")
cat("Total students:", nrow(datos), "\n")

dim(datos)


# Frequencies distribution 
table(datos$cole_depto_ubicacion)
# ===================================================================
#                     COVARIATES CREATION / MATRIXES
# ===================================================================


depto_ids <- sort(unique(datos$cole_cod_depto_ubicacion))
mcpio_ids <- sort(unique(datos$cole_cod_mcpio_ubicacion))

# m : number of departments
m <- length(depto_ids)   # q : departments
# n_j : number of municipalities
n_j <- length(mcpio_ids) # j : municipalities
# n: number of students
n_jq <- nrow(datos)      # i : observations

# Data manipulation
# y      : students score (c)
# Y      : students score as list (list)
# g_dep  : sequential departmental id (c)
# g_dep  : sequential municipal id (c)
# n_j    : number of students per department (c)
# yb     : average score per department (c)
# s2     : departamental variance (c)

y <- datos$punt_global
# y as vector
y_vec <- as.matrix(y)


# ID each observation per its corresponding municipality
g_mcpio_global <- rep(NA_integer_, nrow(datos))
mun_ID         <- interaction(datos$cole_cod_depto_ubicacion,
                        datos$cole_cod_mcpio_ubicacion,
                        drop = TRUE)

g_mcpio_global <- as.integer(factor(mun_ID))


#==============================================================
#@ =====         Y AN COVARIATES ORDERING                ======
#==============================================================

# ==== Prepare covariates matrix ====
# y       : y_{j,q} Scores
# x_jq    : x_{j,q} Individual level covariates
# Z_M    : z_{j,q} Municipal level covariates
# W_D     : w_{q}   Departamental level covariates

X_E_names <- c("fami_educacionmadre_modif", "computador", "internet", "etnia",
                 "libros_11_25", "libros_26_100", "libros_mas100", 
                 "estrato_1", "estrato_2", "estrato_3", "estrato_4", 
                 "estrato_5", "estrato_6", "Genero_mujer", "Calendario_A", 
                 "Calendario_B", "cole_privado", "trabaja_menos_de_10_horas",  
                 "trabaja__11_a_20_horas", "trabaja__21_a_30_horas", 
                 "trabaja_mas_de_30_horas", "cole_urbano")

Z_M_names <-  c("docenttotal_alumtotal", "RISK_VICTIM_2022", "Homi_x_100k_habitantes",
                 "porc_alumn_col_publico", "terrorismot", "hurto", "secuestros", "discapital",
                 "regalias_educa_per_capita")

W_D_names <-  c("PIB_percapita_DPTO", "proporcio_pob_rural", "perc_municipios_con_riesgo", 
                "Homicidios_ponderado_x_100k")

X_E_full <- as.matrix(subset(datos, select = X_E_names))
X_M_full  <- as.matrix(subset(datos, select = Z_M_names))
X_D_full   <- as.matrix(subset(datos, select = W_D_names))



# ==== ESTANDARIZACIÓN DE COVARIABLES CONTINUAS (M y D) ====
# Razón: las columnas de X_M_full tienen escalas que difieren hasta 10 órdenes
# de magnitud (regalias_educa_per_capita ~1e6 vs RISK_VICTIM_2022 ~0.1), lo
# que hace que ZtZ_weighted en sample_beta_k_M sea numéricamente singular y
# chol() falle. X_E_full no se toca: son dummies 0/1.
# El intercepto global beta_int absorbe los corrimientos por centrado.

# Guardar centros y escalas (necesarios para des-estandarizar coeficientes
# después y para transformar datos_testeo en predicción)
X_M_center <- colMeans(X_M_full)
X_M_scale  <- apply(X_M_full, 2, sd)
X_D_center <- colMeans(X_D_full)
X_D_scale  <- apply(X_D_full, 2, sd)

# Guarda contra columnas constantes (no debería pasar con tus datos, pero
# por seguridad)
if (any(X_M_scale == 0) || any(X_D_scale == 0)) {
  stop("Columna con sd = 0 en X_M_full o X_D_full. Revisar datos.")
}

# Aplicar estandarización (Z-score por columna)
X_M_full <- scale(X_M_full, center = X_M_center, scale = X_M_scale)
X_M_full <- unname(as.matrix(X_M_full))   # quitar atributos de scale()
colnames(X_M_full) <- Z_M_names           # restaurar nombres

X_D_full <- scale(X_D_full, center = X_D_center, scale = X_D_scale)
X_D_full <- unname(as.matrix(X_D_full))
colnames(X_D_full) <- W_D_names

cat("Z_full estandarizado. Nuevos rangos por columna:\n")
print(apply(X_M_full, 2, function(x) round(c(min = min(x), max = max(x)), 3)))

#==============================



# ==== CREATE NESTED LISTS =====

Y       <-  vector(mode = "list", length = m)
X_E     <-  vector(mode = "list", length = m)
Z_M     <-  vector(mode = "list", length = m)
W_D     <-  vector(mode = "list", length = m)

g_dep   <-  rep(NA, nrow(datos))
g_mcpio <-  rep(NA, nrow(datos))

for (q in 1:m) {  #loop over deparments
  d_q            <- depto_ids[q] # get department code
  # create boolean. TRUE for students from department q
  idx_dep        <- datos$cole_cod_depto_ubicacion == d_q 
  g_dep[idx_dep] <- q
  
  # Identify and sort unique municipalities inside each department
  mcpio_ids <- sort(unique(datos$cole_cod_mcpio_ubicacion[idx_dep]))
  n_mcpios <- length(mcpio_ids)
  
  # Initialise lists over departments
  
  # Empty lists creation
  Y[[q]]     <-  vector("list", length = n_mcpios)
  X_E[[q]]   <-  vector("list", length = n_mcpios)
  Z_M[[q]]   <-  vector("list", length = n_mcpios)
  W_D[[q]]   <-  vector("list", length = n_mcpios)
  
  for (j in 1:n_mcpios) { # loop for municipalities
    mcpio_j            <- mcpio_ids[j]
    # Create boolean for students in both department and municipality
    idx_mcpio          <- idx_dep & datos$cole_cod_mcpio_ubicacion == mcpio_j
    g_mcpio[idx_mcpio] <- j
    
    # Store observed variables and covariates in lists dept - mun
    Y[[q]][[j]]      <-  datos$punt_global[idx_mcpio]
    X_E[[q]][[j]]    <-  X_E_full[idx_mcpio, , drop = FALSE]
    Z_M[[q]][[j]]    <-  X_M_full[idx_mcpio, , drop = FALSE]
    W_D[[q]][[j]]    <-  X_D_full[idx_mcpio, , drop = FALSE]
  }
}


#==== Table with sufficient statistics per department and municipality ====

# Score

# Per departament
estadisticos_departamentos <- datos %>% 
  group_by(cole_cod_depto_ubicacion) %>% 
  summarise(
    codigo = first(cole_cod_depto_ubicacion),
    nombre = first(cole_depto_ubicacion),
    nq_x_depto = n(),                         # Number of observations per department
    yb = mean(punt_global),                   # Mean
    s2 = var(punt_global)                     # Variance
  ) %>% 
  ungroup() %>% 
  arrange(codigo)

# Per municipality
estadisticos_municipios <- datos %>% 
  group_by(cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion) %>% 
  summarise(
    codigo = first(cole_cod_mcpio_ubicacion),
    nombre = first(cole_mcpio_ubicacion),
    nj_x_mun = n(),                           # Number of observations per municipality
    yb = mean(punt_global, na.rm = TRUE),     # Mean  
    s2 = var(punt_global, na.rm = TRUE),      # Variance
    .groups = "drop") %>%
  select(-codigo)

# Sample sizes
nq_x_depto <- estadisticos_departamentos$nq_x_depto
nj_x_mun   <- estadisticos_municipios$nj_x_mun

# Sufficient statistics
yb_m <- estadisticos_departamentos$yb
yb_j <- estadisticos_municipios$yb
s2_m <- estadisticos_departamentos$s2
s2_j <- estadisticos_municipios$s2




#==============================================================
#@ =====        PARAMETERS E HIPERPARAMETERS             ======
#==============================================================

# Number of parameters
numero_Betas_E  <-  ncol(X_E_full)
numero_Betas_M  <-  ncol(X_M_full)
numero_Betas_D  <-  ncol(X_D_full)
numero_Betas    <-  1 + numero_Betas_E + numero_Betas_M + numero_Betas_D


#=================     DEFINE UNITARY PRIOR     =================

# Join covariate-level matrixes
X_full             <-  cbind(X_E_full, X_M_full, X_D_full)
# MCO regression
lm_previa_unitaria <-  lm(y_vec ~ X_full)
# Beta coefficients
betas_lm           <-  coef(lm_previa_unitaria)

# Intercept
mu_beta  <-  betas_lm[1]
# Beta estimators: Indidivual level
mu_E     <-  betas_lm[ 2 : (1 + numero_Betas_E) ]

# Beta estimators: Municipal level
start_M  <-  2 + numero_Betas_E
end_M    <-  1 + numero_Betas_E + numero_Betas_M
mu_M     <-  betas_lm[ start_M : end_M ]

# Beta estimators: Departamental level
start_D  <-  end_M + 1
end_D    <-  length(betas_lm)
mu_D     <-  betas_lm[ start_D : end_D ]


# If variance was needed
#sigma2_ols_E <- summary(lm_previa_unitaria)$sigma^2
#Sigma0_E <- sigma2_ols_E * solve(t(X_full)%*%X_full) * nrow(y_vec)



#=================     SOME QUANTITIES     =================

# 1s observation vector
ones_vector_beta <- matrix(1, nrow = n_jq, ncol = 1)
# 1S betas vectors
ones_vector_E <- matrix(1, nrow = numero_Betas_E, ncol = 1)
ones_vector_M <- matrix(1, nrow = numero_Betas_M, ncol = 1)
ones_vector_D <- matrix(1, nrow = numero_Betas_D, ncol = 1)

# Product x^tx, z^tz, w^tw
xtx_list = vector("list",m)
ztz_list = vector("list",m)
wtw_list = vector("list",m)
for (q in 1:m) {
  for (j in 1:length(Y[[q]])) {
    n_j_q <- length(Y[[q]][[j]])
    # Sum of each of the addends of the mean
    xtx_list[[q]][[j]] = t(X_E[[q]][[j]]) %*% X_E[[q]][[j]]
    ztz_list[[q]][[j]] = t(Z_M[[q]][[j]]) %*% Z_M[[q]][[j]]
    wtw_list[[q]][[j]] = t(W_D[[q]][[j]]) %*% W_D[[q]][[j]]
  }
}  

# -----------------------------
# Precompute flattened municipality map & flat lists 
# -----------------------------
# Build a flattened map that enumerates municipalities
mun_map <- list()
mun_y    <- list()
mun_x    <- list()
mun_z    <- list()
mun_w    <- list()
mun_xtx  <- list()
mun_ztz  <- list()
mun_wtw  <- list()
mun_sizes <- integer(0)
start_idx <- 1L

for (q in seq_along(Y)) {
  for (j in seq_along(Y[[q]])) {
    n_jq_local <- length(Y[[q]][[j]])
    end_idx <- start_idx + n_jq_local - 1L
    mun_map[[length(mun_map) + 1L]] <- list(q = q, j = j, 
                                            start = start_idx, end = end_idx, n = n_jq_local)
    mun_y[[length(mun_y) + 1L]]   <- as.numeric(Y[[q]][[j]])
    mun_x[[length(mun_x) + 1L]]   <- X_E[[q]][[j]]
    mun_z[[length(mun_z) + 1L]]   <- Z_M[[q]][[j]]
    mun_w[[length(mun_w) + 1L]]   <- W_D[[q]][[j]]
    mun_xtx[[length(mun_xtx) + 1L]] <- xtx_list[[q]][[j]]
    mun_ztz[[length(mun_ztz) + 1L]] <- ztz_list[[q]][[j]]
    mun_wtw[[length(mun_wtw) + 1L]] <- wtw_list[[q]][[j]]
    mun_sizes <- c(mun_sizes, n_jq_local)
    start_idx <- end_idx + 1L
  }
}
n_mun_total <- length(mun_map)       # should equal sum(sapply(Y,length)) 

# initialize kappa2_jq nested and flat representation: Sufficient statistics
# Var per municipality
kappa2_jq <- lapply(Y, function(depto) {
  lapply(depto, function(x) {
    v <- var(x)
    if (is.na(v) || v < 0.1) 0.1 else v  # 0.1 for municipalities with 1 observation. Avoids NA
  })
})
# flattened mun_kappa vector in same order as 'mun_map'
mun_kappa <- numeric(n_mun_total)
idx <- 1L
for (i in seq_len(n_mun_total)) {
  q <- mun_map[[i]]$q
  j <- mun_map[[i]]$j
  mun_kappa[i] <- kappa2_jq[[q]][[j]]
}


# Pre-compute observation-to-municipality mapping
obs_to_mun <- integer(n_jq)
for (idx in seq_along(mun_map)) {
  obs_to_mun[mun_map[[idx]]$start:mun_map[[idx]]$end] <- idx
}

# Pre-compute x_i^T x_i for all observations
# Avoid computing sum(x_i^2) B times
cat("Pre-computing XtX, ZtZ, and WtW values...\n")
xTx_vec <- rowSums(X_E_full^2)
zTz_vec <- rowSums(X_M_full^2)
wTw_vec <- rowSums(X_D_full^2)

# Pre-compute sqrt(kappa^2) for all municipalities
# Avoid repeated sqrt() calls
sqrt_kappa <- sqrt(mun_kappa)



#==============================================================================
#@                       GIBBS SAMPLER DP FUNCTIONS              
#==============================================================================

#=====================     INITIALIZATION   =====================

Z_full <- X_M_full
W_full <- X_D_full

# -----------------------------
# Initialization for beta_k and xi
# -----------------------------
# OLS for each department
p_E <- ncol(X_E_full)
p_M <- ncol(Z_full)
p_D <- ncol(W_full)


# Initialize xi: cluster = department index q
xi_E <- integer(nrow(datos))      # Individual level
xi_M <- integer(nrow(datos))      # Municipal level
xi_D <- integer(nrow(datos))      # Departamental level



# We will initialize with m clusters, m: the number of departments
# and in j clusters: j: number of municipalities 

# Initialize xi with department and municipality indicators
xi_E <- g_dep                 #<---- change g_mcpio_global: Observation's departamental ID
xi_M <- g_dep                                # g_dep: Observation's departamental ID
xi_D <- g_dep


# Initialize beta_k_E, beta_k_M, beta_k_D Using OLS for department and municipality

initialise_beta_k <- function(y, level, g_dep){
  set.seed(777)
  # Retrieve information about the level = {E,M,D}
  level_size <- length(unique(g_dep))                         # level size
  beta_k_ell <- vector("list", level_size)                    # Initialises beta_k_ell
  X_ell  <- get(paste0("X_", level, "_full"))                 # Covariate matrix
  p_ell  <- ncol(X_ell)                                       # Number of parameters
  mu_ell <- get(paste0("mu_", level))                         # OLS vector for NA cases
  
  # Iterate over the unique elements in the level ell
  for (ell in 1:level_size) {
    # Get indices for students in level ell
    idx_q <- which(g_dep == ell)            
    # Extract y and X for department q
    y_q   <- y[idx_q]
    X_q   <- X_ell[idx_q, , drop = FALSE]
    
    # Try to fit OLS without intercept
    tryCatch({
      lm_q       <- lm(y_q ~ X_q)
      beta_coefs <- as.numeric(coef(lm_q)[-1])  
      # Commonly we are going to have NA because multicolinearity
      if (any(is.na(beta_coefs))) {
        na_indices <- which(is.na(beta_coefs))
        for (j in na_indices) {
          # Replaces the value for a N(beta_OLS, 1)
          beta_coefs[j] <- rnorm(1, mean = mu_ell[j], sd = 1)   #Samples from N(beta_OLS, 1)
        }
      }
      
      beta_k_ell[[ell]] <- beta_coefs        
      
    }, error = function(e) {                                     # Runs if tryCatch returns an error
      # If OLS completely fails, sample all coefficients from the prior
      beta_k_ell[[ell]] <<- rnorm(p_ell, mean = mu_ell, sd = 1)  
    })
  }
  # Return the beta_k vector
  return(beta_k_ell)                                                        
}

# Initialise list of vectors beta_k_M assuming the clusters are the municipalities
beta_k_E <- initialise_beta_k(y     = y,
                              level = "E",
                              g_dep = g_dep)   #<----- Considerar cambiarlo inicializar en departamentos
                                                        #       g_mcpio_global

# Initialise list of vectors beta_k_E assuming the clusters are the departments
beta_k_M <- initialise_beta_k(y     = y,
                              level = "M",
                              g_dep = g_dep)

# Initialise list of vectors beta_k_D assuming the clusters are the departments
beta_k_D <- initialise_beta_k(y     = y,
                              level = "D",
                              g_dep = g_dep)


# Initialise parameters with OLS estimators
lm_previa_int = lm(y_vec ~ 1)
summary(lm_previa_int)
sigma2_ols_beta <- summary(lm_previa_int)$sigma^2

# To initialise sigma_base's, compute OLS by level
lm_previa_unitaria_E <-  lm(y_vec ~ X_E_full)           # Individual level covariates
sigma2_E_ols <- summary(lm_previa_unitaria_E)$sigma^2

lm_previa_unitaria_M <-  lm(y_vec ~ X_M_full)           # Municipal level covariates
sigma2_M_ols <- summary(lm_previa_unitaria_M)$sigma^2  

lm_previa_unitaria_D <-  lm(y_vec ~ X_D_full)           # Departamental level covariates
sigma2_D_ols <- summary(lm_previa_unitaria_D)$sigma^2

#=====================           HIPERPARAMETERS           =====================

# For beta_int 
mu_beta = mu_beta

# For sigma2_beta
nu_beta          <-   1
gamma2_beta       <-   sigma2_ols_beta


#For alpha_E
a_alpha_E          <-   1
b_alpha_E          <-   1

#For alpha_M
a_alpha_M          <-   1
b_alpha_M          <-   1

#For alpha_D
a_alpha_D          <-   1
b_alpha_D          <-   1

# For beta_k_E base distribution
  # For mu_E
eta_mu_E           <-   betas_lm[ 2 : (1 + numero_Betas_E) ]
nu2_mu_E           <-   sigma2_E_ols
nu2_mu_E_inv       <-   1 / nu2_mu_E
  # For sigma2_E
nu_E               <-   1
gamma2_E            <-   sigma2_E_ols
nu_gamma2_E         <-   nu_E*gamma2_E


# For beta_k_M base distribution
  # For mu_M
eta_mu_M           <-  betas_lm[ start_M : end_M ]
nu2_mu_M           <-   sigma2_M_ols
nu2_mu_M_inv       <-   1 / nu2_mu_M
  # For sigma2_M
nu_M               <-   1
gamma2_M            <-   sigma2_M_ols
nu_gamma2_M         <-   nu_M*gamma2_M

# For beta_k_D base distribution
  # For mu_D
eta_mu_D           <-  betas_lm[ start_D : end_D ]
nu2_mu_D           <-   sigma2_D_ols
nu2_mu_D_inv       <-   1 / nu2_mu_D
  # For sigma2_D
nu_D               <-   1
gamma2_D            <-   sigma2_D_ols
nu_gamma2_D         <-   nu_D*gamma2_D


# VARIANCE COMPONENT HIPERPARAMETERS
alpha_kappa      <-   2 / (1^2)               # sd = 50 <- expected (m)
beta_kappa       <-   alpha_kappa / 50        

nu_kappa         <-   4 + (2 / 1^2)         # CV = 1


#=====================     PARAMETERS INITIALISATION       =====================

beta_int         <-  mu_beta
mu_E             <-  betas_lm[ 2 : (1 + numero_Betas_E) ]
mu_M             <-  betas_lm[ start_M : end_M ]
mu_D             <-  betas_lm[ start_D : end_D ]

alpha_E          <- 1
alpha_M          <- 1
alpha_D          <- 1
sigma2_beta      <-   vcov(lm_previa_unitaria)[1, 1]  # Variance of intercept from OLS
sigma2_E_base    <- sigma2_E_ols
sigma2_M_base    <- sigma2_M_ols
sigma2_D_base    <- sigma2_D_ols
kappa_inv        <-   1 / mun_kappa  # Avoid repeated division in the loop


# Departamental variance
kappa2_q <- sapply(Y, function(depto) {
  mean(sapply(depto, function(mun) var(mun, na.rm = TRUE) + 0.01), na.rm = TRUE)
})



#  QUANTITIES
mun_sizes_local <- mun_sizes
n_mun_total     <- length(mun_map)

# Pre-compute which municipalities belong to each department
dept_to_mun <- vector("list", m)
for (mun_idx in seq_along(mun_map)) {
  q <- mun_map[[mun_idx]]$q  # Department of this municipality
  dept_to_mun[[q]] <- c(dept_to_mun[[q]], mun_idx)
}

# ===================================================================
#                       EXECUTE THE GIBBS SAMPLER
# ===================================================================

B  = 25000

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



# ===================================================================
#               XI AND BETA_K EXTRACTION FROM TEXT FILES
# ===================================================================


# Retrieve B from cadena
B = cadena$info$samples_stored

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


#-------------------------------
# Verify
#-------------------------------
cat("Summary:\n")
cat("  Number of samples:", nrow(xi_E_samples), "\n")
cat("  Number of observations:", ncol(xi_E_samples), "\n")
cat("  Number of covariates:", p_E, "\n")
cat("  Sample 1 has", length(beta_k_E_samples[[1]]), "clusters\n")
cat("  Sample", B, "has", length(beta_k_E_samples[[B]]), "clusters\n\n")

# Example access
sample_idx <- 1
cat("Example - Sample", sample_idx, ":\n")
cat("  K_E =", length(beta_k_E_samples[[sample_idx]]), "clusters\n")
cat("  First cluster beta (first 5 values):", beta_k_E_samples[[sample_idx]][[1]][1:5], "\n")



# ==============================================================================
#@                     LOG-LIKELIHOOD COMPUTATION 
# ==============================================================================

cat("Computing log-likelihood for", B, "samples...\n")

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
cat("Range: [", round(min(log_lik_vec), 2), ",", round(max(log_lik_vec), 2), "]\n")
cat("Mean :", round(mean(log_lik_vec), 2), "\n")



# ===================================================================
#                   SAVE OR LOAD THE CHAIN
# ===================================================================

#---------------------------------------------------------------------------------------
# Save or load the filtered chain after permutating to solve the label switching problem
#---------------------------------------------------------------------------------------

# Save chain
#save(cadena, file = paste0(path_to_results, "/cadena_Ifces_MCMC3.RData"))


# Load the chain 
load(file = paste0(path_to_results, "/cadena_Ifces_MCMC3.RData"))



# ===================================================================
#                        PLOTS
# ===================================================================


# ===============================================
#@             LOGLIKELIHOOD PLOT
# ===============================================


cat("Summary statistics:\n")
cat("  Min log-lik:  ", sprintf("%.2f", min(cadena$log_likelihood)), "\n")
cat("  Max log-lik:  ", sprintf("%.2f", max(cadena$log_likelihood)), "\n")
cat("  Mean log-lik: ", sprintf("%.2f", mean(cadena$log_likelihood)), "\n")
cat("  Median log-lik:", sprintf("%.2f", median(cadena$log_likelihood)), "\n")
cat("  SD log-lik:   ", sprintf("%.2f", sd(cadena$log_likelihood)), "\n\n")

# Traceplot 

df_loglik <- data.frame(
  iteration = 1:B,
  log_likelihood = cadena$log_likelihood
)

ll_trace <- ggplot(df_loglik, aes(x = iteration, y = log_likelihood)) +
  geom_line(color = "steelblue", alpha = 0.4) +
  geom_point(size = 1.1, alpha = 0.55) +
  labs(title = paste0("Gráfico de traza log-verosimilitud Modelo M3"),
       x = "Iteración MCMC",
       y = "Log-verosimilitud") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
x11()
print(ll_trace)

# ==============================================================================
# THINNED LOG-LIKELIHOOD TRACEPLOTS (subsample every 5, 10, and 25 iterations)
# ==============================================================================

#-----------------------------------------
# The chain suffers from autocorrelation
#Thinned loglikelihood traceplots would 
#help to observe whether higher thinning 
#would aliviate the problem
#-----------------------------------------

# Thinning levels to use
thin_levels <- c(5, 10, 25)

# Storage list for the three plots
ll_trace_thinned <- vector("list", length(thin_levels))
names(ll_trace_thinned) <- paste0("thin_", thin_levels)

for (k in thin_levels) {
  
  # Systematic subsampling: keep iterations 1, 1+k, 1+2k, ...
  idx_keep <- seq(from = 1, to = B, by = k)
  
  # Build data frame for ggplot
  # x-axis = subsample position (1, 2, 3, ..., length(idx_keep))
  # so the range shrinks as k grows
  df_thin <- data.frame(
    iteration      = seq_along(idx_keep),              # subsample position
    log_likelihood = cadena$log_likelihood[idx_keep]
  )
  
  # Build traceplot (same style as ll_trace)
  p <- ggplot(df_thin, aes(x = iteration, y = log_likelihood)) +
    geom_line(color = "steelblue", alpha = 0.4) +
    geom_point(size = 1.1, alpha = 0.55) +
    labs(title = paste0("Gráfico de traza log-verosimilitud Modelo M3 (submuestreo cada ", k, " iteraciones)"),
         x = "Iteración MCMC (submuestra)",
         y = "Log-verosimilitud") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Store and display
  ll_trace_thinned[[paste0("thin_", k)]] <- p
  x11()
  print(p)
  
  # Brief console summary
  cat(sprintf("Thinning cada %2d : %d muestras conservadas de %d originales\n",
              k, length(idx_keep), B))
}


# ==============================================================================
#               EFFECTIVE SAMPLE SIZE (ESS) — Modelo M3
# ==============================================================================

# Helper: ESS for a single numeric chain
ess_of <- function(x) {
  suppressWarnings(as.numeric(effectiveSize(x)))
}

# Total stored samples
B_total <- length(cadena$log_likelihood)

# ------------------------------------------------------------------
# 1. Collect global scalar/coordinate chains into a named list
# ------------------------------------------------------------------

chains_named <- list()

# --- DP concentration parameters ---
chains_named[["alpha_E"]] <- cadena$alpha_E
chains_named[["alpha_M"]] <- cadena$alpha_M
chains_named[["alpha_D"]] <- cadena$alpha_D

# --- Base-measure means (one entry per covariate coordinate) ---
for (j in seq_len(ncol(cadena$mu_E))) {
  nm <- paste0("mu_E[", colnames(cadena$mu_E)[j], "]")
  chains_named[[nm]] <- cadena$mu_E[, j]
}
for (j in seq_len(ncol(cadena$mu_M))) {
  nm <- paste0("mu_M[", colnames(cadena$mu_M)[j], "]")
  chains_named[[nm]] <- cadena$mu_M[, j]
}
for (j in seq_len(ncol(cadena$mu_D))) {
  nm <- paste0("mu_D[", colnames(cadena$mu_D)[j], "]")
  chains_named[[nm]] <- cadena$mu_D[, j]
}

# --- Base-measure variances ---
chains_named[["sigma2_E_base"]] <- cadena$sigma2_E_base
chains_named[["sigma2_M_base"]] <- cadena$sigma2_M_base
chains_named[["sigma2_D_base"]] <- cadena$sigma2_D_base

# --- Intercept ---
chains_named[["beta_int"]]    <- cadena$beta_int
chains_named[["sigma2_beta"]] <- cadena$sigma2_beta

# --- Log-likelihood (global mixing diagnostic) ---
chains_named[["log_likelihood"]] <- cadena$log_likelihood

# ------------------------------------------------------------------
# 2. Compute ESS for each global parameter
# ------------------------------------------------------------------

ess_vals <- sapply(chains_named, ess_of)

ess_df <- data.frame(
  parametro = names(ess_vals),
  ESS       = round(ess_vals, 1),
  pct_of_B  = round(100 * ess_vals / B_total, 2),
  row.names = NULL,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# 3. Hierarchical variances — kappa2_q (33 depts) and kappa2_jq
#    Compute ESS per parameter, then report only min/median/max
# ------------------------------------------------------------------

# kappa2_q : list of m vectors (one per department)
ess_kappa2_q <- sapply(cadena$kappa2_q, ess_of)

# kappa2_jq : nested list [dept][mun] — flatten across all municipalities
ess_kappa2_jq <- unlist(lapply(cadena$kappa2_jq, function(dept) {
  sapply(dept, ess_of)
}))

# Summary rows
summary_kappa_q <- data.frame(
  parametro = c(sprintf("kappa2_q  [min, n=%d]",    length(ess_kappa2_q)),
                sprintf("kappa2_q  [median, n=%d]", length(ess_kappa2_q)),
                sprintf("kappa2_q  [max, n=%d]",    length(ess_kappa2_q))),
  ESS       = round(c(min(ess_kappa2_q),
                      median(ess_kappa2_q),
                      max(ess_kappa2_q)), 1),
  pct_of_B  = round(100 * c(min(ess_kappa2_q),
                            median(ess_kappa2_q),
                            max(ess_kappa2_q)) / B_total, 2),
  stringsAsFactors = FALSE
)

summary_kappa_jq <- data.frame(
  parametro = c(sprintf("kappa2_jq [min, n=%d]",    length(ess_kappa2_jq)),
                sprintf("kappa2_jq [median, n=%d]", length(ess_kappa2_jq)),
                sprintf("kappa2_jq [max, n=%d]",    length(ess_kappa2_jq))),
  ESS       = round(c(min(ess_kappa2_jq),
                      median(ess_kappa2_jq),
                      max(ess_kappa2_jq)), 1),
  pct_of_B  = round(100 * c(min(ess_kappa2_jq),
                            median(ess_kappa2_jq),
                            max(ess_kappa2_jq)) / B_total, 2),
  stringsAsFactors = FALSE
)

ess_df <- rbind(ess_df, summary_kappa_q, summary_kappa_jq)

# ------------------------------------------------------------------
# 4. Print results
# ------------------------------------------------------------------

cat("\n=================================================================\n")
cat(" Effective Sample Size (ESS) - Modelo M3\n")
cat(" Total stored samples B =", B_total, "\n")
cat("=================================================================\n")
print(ess_df, row.names = FALSE)
cat("=================================================================\n")
cat(sprintf(" Resumen general (parametros globales):\n"))
cat(sprintf("   ESS minimo : %.1f  (%.2f%% de B)\n",
            min(ess_vals), 100 * min(ess_vals) / B_total))
cat(sprintf("   ESS mediano: %.1f  (%.2f%% de B)\n",
            median(ess_vals), 100 * median(ess_vals) / B_total))
cat(sprintf("   ESS maximo : %.1f  (%.2f%% de B)\n",
            max(ess_vals), 100 * max(ess_vals) / B_total))
cat("=================================================================\n\n")


# ==============================================================================================
#@                             ADJUSTED RAND INDEX (ARI)
# ==============================================================================================


n_iter <- nrow(cadena$xi_E)

ARI_E <- numeric(n_iter)
ARI_M <- numeric(n_iter)
ARI_D <- numeric(n_iter)

# Compute the ARI per iteration
for (b in seq(n_iter)) {
  ARI_E[b] <- adjustedRandIndex(
    g_dep, 
    as.vector(cadena$xi_E[b, ])
  )
  
  ARI_M[b] <- adjustedRandIndex(
    g_dep, 
    as.vector(cadena$xi_M[b, ])
  )
  
  ARI_D[b] <- adjustedRandIndex(
    g_dep, 
    as.vector(cadena$xi_D[b, ])
  )
}

# Mean of ARI_E iterations
ARI_E_mean <- mean(ARI_E)
ARI_M_mean <- mean(ARI_M)
ARI_D_mean <- mean(ARI_D)


# There is no evidence that the model retrieves a geographic aggrupation




# ==============================================================================================
#@                                    INFERENCE
# ==============================================================================================


# ==============================================================================
#@                 INFERENCE ON THE NUMBER OF CLUSTERS
# ==============================================================================

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
# Plot: Posterior distribution of K (E, M, D) — Modelo M3
# One figure with 3 subplots
# ----------------------------------------------------------------
x11()
par(mfrow = c(1, 3), mar = c(3, 3, 2, 1.4), mgp = c(1.75, 0.75, 0))

# --- E-level ---
plot(as.numeric(names(k_table_E)), as.numeric(k_table_E),
     type = "h", lwd = 3, col = "dodgerblue4",
     main = "Distribución posterior de K_E - Modelo M3",
     xlab = "Número de clusters", ylab = "Densidad",
     yaxt = "n")
axis(2, at = pretty(k_table_E), labels = pretty(k_table_E))

# --- M-level ---
plot(as.numeric(names(k_table_M)), as.numeric(k_table_M),
     type = "h", lwd = 3, col = "dodgerblue4",
     main = "Distribución posterior de K_M - Modelo M3",
     xlab = "Número de clusters", ylab = "Densidad",
     yaxt = "n")
axis(2, at = pretty(k_table_M), labels = pretty(k_table_M))

# --- D-level ---
plot(as.numeric(names(k_table_D)), as.numeric(k_table_D),
     type = "h", lwd = 3, col = "dodgerblue4",
     main = "Distribución posterior de K_D - Modelo M3",
     xlab = "Número de clusters", ylab = "Densidad",
     yaxt = "n")
axis(2, at = pretty(k_table_D), labels = pretty(k_table_D))

# ----------------------------------------------------------------
# Posterior summaries of K (mode, mean, 95% credible interval)
# ----------------------------------------------------------------
cat("\n=== Posterior summary of K - Modelo M3 ===\n")

summarise_K <- function(k_values, level) {
  k_tab  <- table(k_values)
  k_mode <- as.integer(names(k_tab)[which.max(k_tab)])
  cat(sprintf("  K_%s : moda = %d | media = %.2f | IC95%% = [%d, %d] | P(K=moda) = %.3f\n",
              level,
              k_mode,
              mean(k_values),
              as.integer(quantile(k_values, 0.025)),
              as.integer(quantile(k_values, 0.975)),
              max(k_tab) / length(k_values)))
}

summarise_K(k_values_E, "E")
summarise_K(k_values_M, "M")
summarise_K(k_values_D, "D")



# ===============================================
#@     FILTER THE CHAIN BY k_mode ITERATIONS
# ===============================================

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

cadena_filtered <- filter_cadena_by_k_3DP(cadena, k_E_mode, k_M_mode, k_D_mode)




# ===============================================
#@          SOLVE LABEL SWITCHING
# ===============================================


#cadena_filtered <- solve_label_switching_3DP(cadena_filtered, y, X_E_full, Z_full, W_full,
#                                             k_E_mode, k_M_mode, k_D_mode)


cadena_filtered <- solve_label_switching_3DP_cpp(cadena_filtered, y, X_E_full, Z_full, W_full,
                                             k_E_mode, k_M_mode, k_D_mode)




# ===================================================================
#                   SAVE OR LOAD THE FILTERED CHAIN
# ===================================================================

#---------------------------------------------------------------------------------------
# Save or load the filtered chain after permutating to solve the label switching problem
#---------------------------------------------------------------------------------------

# Save chain
save(cadena_filtered, file = paste0(path_to_results, "/cadena_filtered_Ifces_MCMC3.RData"))


# Load the chain 
#load(file = paste0(path_to_results, "/cadena_filtered_Ifces_MCMC3.RData"))



# ===============================================
# COMPUTE EACH CLUSTER BETAs MEAN AND QUANTILES
# ===============================================

cadena_filtered <- compute_posterior_summaries_3DP(cadena_filtered, k_E_mode, 
                                                   k_M_mode, k_D_mode)



# =====================================================================================
# =====================================================================================
# =====================================================================================
# =====================================================================================




# ==============================================================================================
#@                                EXTERNAL VALIDATION (PREDICTION)
# ==============================================================================================

#---------------------------------------------------------------------
# Using the testing samples, predictive capacity metrics will be 
#computed (MSE, MAE, R2)
#---------------------------------------------------------------------

prepare_test_data <- function(datos_testeo, X_E_names, 
                              Z_M_names, W_D_names) {
  
  # ---- Response variable ----
  y <- datos_testeo$punt_global
  
  # ---- Covariate matrix (no intercept column, consistent with training) ----
  # Individual level
  X_full <- as.matrix(subset(datos_testeo, select = X_E_names))
  # Municipal level
  Z_full <- as.matrix(subset(datos_testeo, select = Z_M_names))
  # Departamental level
  W_full <- as.matrix(subset(datos_testeo, select = W_D_names))
  
  # ---- Return list ----
  return(list(
    X_full = X_full,   # matrix: n_test x p_E
    Z_full = Z_full,   # matrix: n_test x p_M
    W_full = W_full,   # matrix: n_test x p_D
    y      = y         # vector: length n_test
  ))
}

test_data <- prepare_test_data(datos_testeo, X_E_names, 
                               Z_M_names, W_D_names)

# ===================================================================
#   COMPUTE MSE, RMSE, MAE, R2
# ===================================================================

cadena_filtered$external_validation <- compute_test_metrics_3DP(cadena_filtered, 
                                                                test_data)

# ===================================================================
#   EXECUTE WAIC 
# ===================================================================

cadena_filtered$WAIC <- compute_WAIC_3DP(cadena_filtered, y, X_E_full, Z_full, 
                                         W_full, obs_to_mun, mun_map)



# =====================================================================================
# =====================================================================================
# =====================================================================================
# =====================================================================================




# =====================================================================================
# =====================================================================================
# =====================================================================================
# =====================================================================================















