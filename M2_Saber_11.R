#@ Modelo 2

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
ruta <- setwd('path')

# Path for results
path_to_results <- paste0(ruta,'/Resultados/Modelo_2/MCMC2')
# Path for text files
path_text_files <- paste0(ruta,'/Resultados/Modelo_2/text_files')

# Create dir (folder) if does not exist for text files
if (!dir.exists(path_text_files)) dir.create(path_text_files, recursive = TRUE, showWarnings = FALSE)

# Create dir (folder) if does not exist for MCMC2
if (!dir.exists(path_to_results)) dir.create(path_to_results, recursive = TRUE, showWarnings = FALSE)



# Load train data 
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
m <- length(depto_ids)              # q : departments
# n_j : number of municipalities
n_j <- length(mcpio_ids)            # j : municipalities
# n: number of students
n_jq <- nrow(datos)                 # i : observations

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

# Departments and municipalities ID
g_dep <- rep(NA, m)
g_mcpio <- rep(NA, n_j)


#==============================================================
#@ =====         Y AN COVARIATES ORDERING                ======
#==============================================================

# ==== Prepare covariates matrix ====
# y   :   y_{j,q} Scores
# x   :   X       Covariates


x_names <- c("fami_educacionmadre_modif", "computador", "internet", "etnia",
             "libros_11_25", "libros_26_100", "libros_mas100", 
             "estrato_1", "estrato_2", "estrato_3", "estrato_4", 
             "estrato_5", "estrato_6", "Genero_mujer", "Calendario_A", 
             "Calendario_B", "cole_privado", "trabaja_menos_de_10_horas",  
             "trabaja__11_a_20_horas", "trabaja__21_a_30_horas", 
             "trabaja_mas_de_30_horas", "cole_urbano", "docenttotal_alumtotal", 
             "RISK_VICTIM_2022", "Homi_x_100k_habitantes",
             "porc_alumn_col_publico", "terrorismot", "hurto", "secuestros", "discapital",
             "regalias_educa_per_capita", "PIB_percapita_DPTO", "proporcio_pob_rural", 
             "perc_municipios_con_riesgo", "Homicidios_ponderado_x_100k")

X_full <- as.matrix(subset(datos, select = x_names))
# Note: The intercept is modeled apart from the X matrix, that is why is not included
#in X_full, it is later include as a 1s column


# ==== Crear listas anidadas =====
Y      <- vector(mode = "list", length = m)
X      <- vector(mode = "list", length = m)

g_dep   <- rep(NA, nrow(datos))
g_mcpio <- rep(NA, nrow(datos))

for (q in 1:m) {  #loop over deparments
  d_q <- depto_ids[q] # get department code
  # create boolean. TRUE for students from department q
  idx_dep <- datos$cole_cod_depto_ubicacion == d_q 
  g_dep[idx_dep] <- q
  
  # Municipios únicos dentro del departamento q
  mcpio_ids <- sort(unique(datos$cole_cod_mcpio_ubicacion[idx_dep]))
  n_mcpios <- length(mcpio_ids)
  
  # Inicializar listas por departamento
  
  # Crear listas vacias
  Y[[q]]     <- vector("list", length = n_mcpios)
  X[[q]]     <- vector("list", length = n_mcpios)
  
  for (j in 1:n_mcpios) { # loop for municipalities
    mcpio_j <- mcpio_ids[j]
    # Create boolean for students in both department and municipality
    idx_mcpio <- idx_dep & datos$cole_cod_mcpio_ubicacion == mcpio_j
    g_mcpio[idx_mcpio] <- j
    
    # Guardar variables observadas y covariables en listas depto - municipio
    Y[[q]][[j]]     <- datos$punt_global[idx_mcpio]
    X[[q]][[j]]     <- X_full[idx_mcpio, , drop = FALSE]
  }
}

#==== Tabla con estadisticas suficientes por departamento y municipio ====

#puntaje

# Por departamento
estadisticos_departamentos <- datos %>% 
  group_by(cole_cod_depto_ubicacion) %>% 
  summarise(
    codigo = first(cole_cod_depto_ubicacion),
    nombre = first(cole_depto_ubicacion),
    nq_x_depto = n(),  # se calcula el numero de observaciones por departamento
    yb = mean(punt_global), #media
    s2 = var(punt_global) #varianza
  ) %>% 
  ungroup() %>% 
  arrange(codigo)

# Por municipio
estadisticos_municipios <- datos %>% 
  group_by(cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion) %>% 
  summarise(
    codigo = first(cole_cod_mcpio_ubicacion),
    nombre = first(cole_mcpio_ubicacion),
    nj_x_mun = n(), # se calcula el numero de datos por cada municipio
    yb = mean(punt_global, na.rm = TRUE), #media  
    s2 = var(punt_global, na.rm = TRUE), #varianza
    .groups = "drop") %>%
  select(-codigo)

### Tamaños de muestra
nq_x_depto <- estadisticos_departamentos$nq_x_depto
nj_x_mun <- estadisticos_municipios$nj_x_mun

### Estadísticos suficientes
yb_m <- estadisticos_departamentos$yb
yb_j <- estadisticos_municipios$yb
s2_m <- estadisticos_departamentos$s2
s2_j <- estadisticos_municipios$s2

#==============================================================
#@ =====        PARAMETERS E HIPERPARAMETERS             ======
#==============================================================

# Numero parametros (Sin parametros de cluster)
numero_Betas   = ncol(X_full)  
numero_parametros = numero_Betas + 4 + m + n_j 

num_tot_parametros = numero_Betas + numero_parametros

###@ hiperparámetros No informativos
nu_beta = 1
#nu_E = 1
#nu_M = 1
#nu_D = 1
gamma_beta = 1
#gamma_E = 1
#gamma_M = 1
#gamma_D = 1
nu_kappa = 1
alpha_kappa = 5
beta_kappa = 5

#=================     DEFINE UNITARY PRIOR     =================

# Define a MCO regresion
lm_previa_unitaria = lm(y_vec ~ X_full)
betas_lm = coef(lm_previa_unitaria)

# Intercept
mu_beta      <- betas_lm[1]
# MCO Coeficient estimators
mu_cov_betas <- betas_lm[-1]

# If variance was needed

#sigma2_ols <- summary(lm_previa_unitaria)$sigma^2

#temp_var <- summary(lm_previa_unitaria)$coefficients["(Intercept)", "Std. Error"]
# n \sigma^2_{ols} (X^T X)^{-1}
#Sigma0_E <- sigma2_ols * solve(t(X_full)%*%X_full) * nrow(y_vec)


#==============================================================
#@ =====             SOME QUANTITIES                     ======
#==============================================================

# Vector de 1s de individuos
ones_vector_beta <- matrix(1, nrow = n_jq, ncol = 1)
# Vector de 1s de betas
ones_vector_covariables <- matrix(1, nrow = numero_Betas, ncol = 1)


# Producto x^tx, z^tz, w^tw
xtx_list = vector("list",m)

for (q in 1:m) {
  for (j in 1:length(Y[[q]])) {
    n_j_q <- length(Y[[q]][[j]])
    # Sumatoria de cada uno de los sumandos de la media
    xtx_list[[q]][[j]] = t(X[[q]][[j]]) %*% X[[q]][[j]]
  }
}  

# -----------------------------
# Precompute flattened municipality map & flat lists 
# -----------------------------
# Build a flattened map that enumerates municipalities
mun_map <- list()
mun_y    <- list()
mun_x    <- list()
mun_xtx  <- list()
mun_sizes <- integer(0)
start_idx <- 1L 

for (q in seq_along(Y)) {
  for (j in seq_along(Y[[q]])) {
    n_jq_local <- length(Y[[q]][[j]])
    end_idx <- start_idx + n_jq_local - 1L
    mun_map[[length(mun_map) + 1L]] <- list(q = q, j = j, start = start_idx, end = end_idx, n = n_jq_local)
    mun_y[[length(mun_y) + 1L]]   <- as.numeric(Y[[q]][[j]])
    mun_x[[length(mun_x) + 1L]]   <- X[[q]][[j]]
    mun_xtx[[length(mun_xtx) + 1L]] <- xtx_list[[q]][[j]]
    mun_sizes <- c(mun_sizes, n_jq_local)
    start_idx <- end_idx + 1L
  }
}


n_mun_total <- length(mun_map)  # should equal sum(sapply(Y,length))
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


# Pre-compute observation-to-municipality mapping
obs_to_mun <- integer(n_jq)
for (idx in seq_along(mun_map)) {
  obs_to_mun[mun_map[[idx]]$start:mun_map[[idx]]$end] <- idx
}

# Pre-compute x_i^T x_i for all observations
# Avoid computing sum(x_i^2) B times
cat("Pre-computing x^T x values...\n")
xTx_vec <- rowSums(X_full^2)

# Pre-compute sqrt(kappa^2) for all municipalities
# Avoid repeated sqrt() calls
sqrt_kappa <- sqrt(mun_kappa)




#==============================================================================
#@                       GIBBS SAMPLER DP FUNCTIONS              
#==============================================================================

#=====================     INITIALIZATION   =====================

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

# -----------------------------
# Initialization for beta_k and xi
# -----------------------------
# OLS for each department
p <- ncol(X_full)

# Initialize beta_k: one OLS per department
beta_k <- vector("list", m)

# Initialize xi: cluster = department index q
xi <- integer(nrow(datos))



# We will initialize with m clusters, m: the number of departments
# Initialize xi with department indicators
xi <- g_dep

# Initialize beta_k: one OLS per department 
beta_k <- vector("list", m)

#Initialize beta_k: one OLS per department 
beta_k <- vector("list", m)

# Compute the initial beta_k with defensive programming for NA
# that are present for departamental variables

for (q in 1:m) {
  # Get indices for students in department q
  idx_q <- which(g_dep == q)
  
  # Extract y and X for department q
  y_q <- y[idx_q]
  X_q <- X_full[idx_q, , drop = FALSE]
  
  # Try to fit OLS without intercept
  tryCatch({
    lm_q <- lm(y_q ~ X_q)
    beta_coefs <- coef(lm_q)[-1]  # Remove intercept
    
    # Check for NA values and replace them using the OLS for all data
    if (any(is.na(beta_coefs))) {
      
      # Replace each NA with a sample from N(mu_cov_betas[j], 1)
      na_indices <- which(is.na(beta_coefs))
      for (j in na_indices) {
        # Replaces the value for a N(beta_OLS, 1)
        beta_coefs[j] <- rnorm(1, mean = mu_cov_betas[j], sd = 1)
      }
    }
    
    beta_k[[q]] <- beta_coefs
    
  }, error = function(e) {
    # If OLS completely fails, sample all coefficients from the prior
    beta_k[[q]] <<- rnorm(p, mean = mu_cov_betas, sd = 1)
  })
}


#=====================           HIPERPARAMETERS           =====================

lm_previa_intento = lm(y_vec ~ 1)
summary(lm_previa_intento)
sigma2_ols_beta <- summary(lm_previa_intento)$sigma^2

sigma2_ols <- summary(lm_previa_unitaria)$sigma^2

#   HIPERPARAMETERS

# For beta_int
mu_beta = mu_beta

#temp_var <- summary(lm_previa_unitaria)$coefficients["(Intercept)", "Std. Error"]
nu_beta     <-   1
#gamma_beta <-   1
gamma_beta  <-   sigma2_ols_beta
alpha       <-   1


# For sigma2_beta
nu_beta          <-   1

# For beta_k  base distribution
eta_mu           <-   mu_cov_betas
nu2_mu           <-   sigma2_ols
nu2_mu_inv       <-   1 / nu2_mu

a_sigma2         <-   1
b_sigma2         <-   sigma2_ols
#b_sigma2        <-   10
a_b_sigma2  <-   a_sigma2*b_sigma2

#For alpha
a_alpha          <-   1
b_alpha          <-   1

# VARIANCE COMPONENT HIPERPARAMETERS
alpha_kappa      <-   2 / (1^2)               # sd = 50 <- expected (m)
beta_kappa       <-   alpha_kappa / 50        

nu_kappa         <-   4 + (2 / 1^2)         # CV = 1


#        PARAMETERS INITIALISATION
beta_int         <-   mu_beta
#sigma2_base     <-   5
sigma2_base      <-   sigma2_ols
sigma2_beta      <-   vcov(lm_previa_unitaria)[1, 1]  # Variance of intercept from OLS
mu_vec           <-   mu_cov_betas 
alpha            <-   1
kappa_inv        <-   1 / mun_kappa  # Avoid repeated division in the loop

# IMPROPER PRIOR
mu_vec = mu_cov_betas
eta_mu = mu_cov_betas
nu2_mu = sigma2_ols

nu2_mu_inv <- 1 / nu2_mu  # Constant, compute once

# Departamental variance
kappa2_q <- sapply(Y, function(depto) {
  mean(sapply(depto, function(mun) var(mun, na.rm = TRUE) + 0.01), na.rm = TRUE)
})

#  QUANTITIES
mun_sizes_local <- mun_sizes
n_mun_total <- length(mun_map)

# Pre-compute which municipalities belong to each department
dept_to_mun <- vector("list", m)
for (mun_idx in seq_along(mun_map)) {
  q <- mun_map[[mun_idx]]$q  # Department of this municipality
  dept_to_mun[[q]] <- c(dept_to_mun[[q]], mun_idx)
}


# ===================================================================
#                       EXECUTE THE GIBBS SAMPLER
# ===================================================================

tictoc::tic()
cadena <- MCMC2_BNP(
  B = 25000,
  nu_beta = nu_beta,
  gamma_beta = gamma_beta,
  nu_kappa = nu_kappa,
  alpha_kappa = alpha_kappa,
  beta_kappa = beta_kappa
)

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

# Example access
sample_idx <- 1
cat("Example - Sample", sample_idx, ":\n")
cat("  K =", length(beta_k_samples[[sample_idx]]), "clusters\n")
cat("  First cluster beta (first 5 values):", beta_k_samples[[sample_idx]][[1]][1:5], "\n")



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



# ===================================================================
#                       SAVE OR LOAD THE CHAIN 
# ===================================================================
# Save chain
save(cadena, file = paste0(path_to_results, "/cadena_Ifces_MCMC2.RData"))


# Load the chain 
load(file = paste0(path_to_results, "/cadena_Ifces_MCMC2.RData"))



# ===================================================================
#                        PLOTS
# ===================================================================


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
  labs(title = paste0("Gráfico de traza log-verosimilitud Modelo 2"),
       x = "Iteración MCMC",
       y = "Log-verosimilitud") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
x11()
print(ll_trace)


# ==============================================================================================
#@                             ADJUSTED RAND INDEX (ARI)
# ==============================================================================================


n_iter <- nrow(cadena$xi)

ARI <- numeric(n_iter)

# Compute the ARI per iteration
for (b in seq(n_iter)) {
  ARI[b] <- adjustedRandIndex(
    g_dep, 
    as.vector(cadena$xi[b, ])
  )
}

# Mean of ARI iterations
ARI_mean <- mean(ARI)

# There is no evidence that the model retrieves a geographic aggrupation

# ==============================================================================================
#@                                    INFERENCE
# ==============================================================================================

# ==============================================================================
#@                 INFERENCE ON THE NUMBER OF CLUSTERS
# ==============================================================================


# --- 1) Calcular K en cada iteración ---
# MARGIN = 1 -> aplica por filas (cada fila = una iteración del MCMC)
k_values <- apply(cadena$xi, 1, function(x) length(unique(x)))

# --- 2) Distribución posterior de K ---
k_table <- table(k_values) / length(k_values)

# --- 3) Gráfico ---
x11()
par(mfrow = c(1, 1), mar = c(3, 3, 1.4, 1.4), mgp = c(1.75, 0.75, 0))
plot(as.numeric(names(k_table)), k_table,
     type = "h", lwd = 3, col = "dodgerblue4",
     main = "Distribución posterior del número de clusters del Modelo 2",
     xlab = "Número de clusters", ylab = "Densidad",
     yaxt = "n")
axis(2, at = pretty(k_table), labels = pretty(k_table))

# --- 4) (Opcional) Resumen numérico ---
cat("Moda posterior de K:", as.numeric(names(which.max(k_table))), "\n")
cat("Media posterior de K:", round(mean(k_values), 2), "\n")
cat("Mediana posterior de K:", median(k_values), "\n")
cat("IC 95% para K: [", quantile(k_values, 0.025), ",",
    quantile(k_values, 0.975), "]\n")




# ===============================================
#@     FILTER THE CHAIN BY k_mode ITERATIONS
# ===============================================

# Compute the number of clusters K at each iteration. 
k_values <- apply(cadena$xi, 1, function(x) length(unique(x))) #note: The 1 applies for rows

# Compute the posterior distribution of K (proportion of each k)
k_table <- table(k_values) / length(k_values)

# Find the mode of k (most repeated k) 
k_mode <- as.numeric(names(k_table)[which.max(k_table)])


# Apply filter_cadena_by_k 

cadena_filtered <- filter_cadena_by_k(cadena, k_mode)


# ===============================================
#@          SOLVE LABEL SWITCHING
# ===============================================


cadena_filtered <- solve_label_switching_cpp(
                      cadena_filtered = cadena_filtered,                   # Filtered MCMC chain
                      y               = y,                                 # The real y data
                      X_full          = X_full,                            # The real X data
                      k_mode          = cadena_filtered$filter_info$k_mode # The mode of k clusters
                    )

# ===================================================================
#                   SAVE OR LOAD THE FILTERED CHAIN
# ===================================================================

#---------------------------------------------------------------------------------------
# Save or load the filtered chain after permutating to solve the label switching problem
#---------------------------------------------------------------------------------------

# Save chain
#save(cadena_filtered, file = paste0(path_to_results, "/cadena_filtered_Ifces_MCMC2.RData"))


# Load the chain 
load(file = paste0(path_to_results, "/cadena_filtered_Ifces_MCMC2.RData"))



# ===============================================
# COMPUTE EACH CLUSTER BETAs MEAN AND QUANTILES
# ===============================================

cadena_filtered <- compute_posterior_summaries(
            cadena_filtered      <- cadena_filtered,                    # Filtered MCMC chain
            k_mode               <- cadena_filtered$filter_info$k_mode  # K_mode
          )




# ==============================================================================================
#@                                EXTERNAL VALIDATION (PREDICTION)
# ==============================================================================================

#---------------------------------------------------------------------
# Using the testing samples, predictive capacity metrics will be 
#computed (MSE, MAE, R2)
#---------------------------------------------------------------------

prepare_test_data <- function(datos_testeo, x_names) {
  
  # ---- Response variable ----
  y <- datos_testeo$punt_global
  
  # ---- Covariate matrix (no intercept column, consistent with training) ----
  X_full <- as.matrix(subset(datos_testeo, select = x_names))
  
  # ---- Return list ----
  return(list(
    X_full = X_full,   # matrix: n_test x p
    y      = y         # vector: length n_test
  ))
}

test_data <- prepare_test_data(datos_testeo, x_names)

# ===================================================================
#   COMPUTE MSE, RMSE, MAE, R2
# ===================================================================

cadena_filtered$external_validation <- compute_test_metrics(
  cadena_filtered      <- cadena_filtered,
  test_data            <- test_data
)

# ===================================================================
#                         EXECUTE WAIC — MODEL 2
# ===================================================================

cadena_filtered$WAIC <- compute_WAIC_model2(cadena_filtered = cadena_filtered,
                                   y               = y,
                                   X_full          = X_full,
                                   obs_to_mun      = obs_to_mun,
                                   mun_map         = mun_map)




# ============================================================================
# Heatmap posterior coeficients per cluster (significance)
# ============================================================================

# ---------------------------------------------------------------------------
# 1. Configuration
# ---------------------------------------------------------------------------

normalize <- TRUE       # TRUE: normalises each row by its absolute max

n_clusters <- 7
n_vars     <- length(x_names)

# Labels (diferent as x_names)
var_labels <- c(
  "Educación de la madre",
  "Computador (hogar)",
  "Internet (hogar)",
  "Etnia (minoría)",
  "Libros en el hogar (11-25)",
  "Libros en el hogar (26-100)",
  "Libros en el hogar (>100)",
  "Estrato 1",
  "Estrato 2",
  "Estrato 3",
  "Estrato 4",
  "Estrato 5",
  "Estrato 6",
  "Género (mujer)",
  "Calendario A",
  "Calendario B",
  "Colegio privado",
  "Trabaja <10 h/sem",
  "Trabaja 11-20 h/sem",
  "Trabaja 21-30 h/sem",
  "Trabaja >30 h/sem",
  "Colegio urbano",
  "Razón alumnos/docente",
  "Riesgo de victimización",
  "Homicidios por 100k hab.",
  "% alumnos en col. público",
  "Terrorismo",
  "Hurto",
  "Secuestros",
  "Distancia a capital",
  "Regalías educ. per cápita",
  "PIB per cápita (dpto.)",
  "Proporción pob. rural",
  "% municipios en riesgo",
  "Homicidios pond. x 100k"
)

# thematic clustering
bloque <- c(
  rep("Individuales y hogar",   14),
  rep("Situación laboral",       4),
  rep("Colegio",                 5),
  rep("Contexto municipal",      7),
  rep("Contexto departamental",  5)
)

# ---------------------------------------------------------------------------
# 2. DF construction with means and CI
# ---------------------------------------------------------------------------

df_list <- list()
for (k in 1:n_clusters) {
  df_k <- data.frame(
    variable = var_labels,
    bloque   = bloque,
    cluster  = paste0("C", k),
    mean     = cadena_filtered$beta_k_summary[[k]]$mean,
    q025     = cadena_filtered$beta_k_summary[[k]]$q025,
    q975     = cadena_filtered$beta_k_summary[[k]]$q975
  )
  df_k$excluye_cero <- with(df_k, (q025 > 0) | (q975 < 0))
  df_list[[k]] <- df_k
}
df_long <- do.call(rbind, df_list)

# ---------------------------------------------------------------------------
# 3. Normalisation per covariables (preserves sign + -)
# ---------------------------------------------------------------------------

if (normalize) {
  df_long <- df_long %>%
    group_by(variable) %>%
    mutate(mean_norm = mean / max(abs(mean), na.rm = TRUE)) %>%
    ungroup()
} else {
  df_long$mean_norm <- df_long$mean
}

# ---------------------------------------------------------------------------
# 4. Order of the factors 
# ---------------------------------------------------------------------------

df_long$variable <- factor(df_long$variable, levels = rev(var_labels))
df_long$cluster  <- factor(df_long$cluster,  levels = paste0("C", 1:n_clusters))

# DIvisor lines to separate blocks 
sep_positions <- cumsum(rev(rle(rev(bloque))$lengths))
sep_positions <- sep_positions[-length(sep_positions)] + 0.5

# ---------------------------------------------------------------------------
# 5. Heatmap ploting
# ---------------------------------------------------------------------------

p <- ggplot(df_long, aes(x = cluster, y = variable, fill = mean_norm)) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  geom_point(
    data = subset(df_long, excluye_cero),
    aes(x = cluster, y = variable),
    color = "black", size = 1.4, inherit.aes = FALSE
  ) +
  geom_hline(
    yintercept = sep_positions,
    color = "grey30", linewidth = 0.5
  ) +
  scale_fill_gradient2(
    low      = "#B2182B",     
    mid      = "white",
    high     = "#2166AC",     
    midpoint = 0,
    limits   = if (normalize) c(-1, 1) else NULL,
    name     = if (normalize) "Coef.\nnormalizado" else "Media\nposterior"
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title = "Mapa de calor de coeficientes posteriores por clúster - Modelo M2",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 13, hjust = 0.5,
                                     margin = margin(b = 10)),
    axis.text.x       = element_text(face = "bold", size = 11),
    axis.text.y       = element_text(size = 9),
    panel.grid        = element_blank(),
    legend.position   = "right",
    legend.key.height = unit(1.2, "cm")
  )

X11()
print(p)



# Build wide-format table of normalized values
df_wide_norm <- df_long %>%
  select(variable, cluster, mean_norm) %>%
  tidyr::pivot_wider(names_from = cluster, values_from = mean_norm) %>%
  arrange(match(variable, var_labels))   # reorder rows to match var_labels

print(df_wide_norm, n = Inf)


table(cadena_filtered$xi)



#===================
#  Cluster analysis
#===================

# Posterior mode of the cluster assignments
# (Most frequent assignment per observation across iterations)

# Option 1: cluster sizes at a single posterior draw (e.g., the last sample)
last_xi <- cadena_filtered$xi[, ncol(cadena_filtered$xi)]
cat("Cluster sizes at last MCMC iteration:\n")
print(table(last_xi))

# Option 2: average cluster sizes across iterations
mean_sizes <- rowMeans(apply(cadena_filtered$xi, 2, function(xi_b) {
  tabulate(xi_b, nbins = 7)
}))
cat("\nMean cluster sizes across iterations:\n")
print(round(mean_sizes, 1))
cat("\nRelative proportions (%):\n")
print(round(100 * mean_sizes / sum(mean_sizes), 2))



  
# Compute ARI between latent partition and observable categorical variables
library(mclust)  # for adjustedRandIndex

# Use a representative posterior assignment (e.g., last iteration or modal assignment)
xi_post <- cadena_filtered$xi[, ncol(cadena_filtered$xi)]

# Adapt these variable names to your dataset
cat("ARI between latent partition and observable categorical variables:\n")
cat("  Estrato:           ", adjustedRandIndex(xi_post, datos$estrato),           "\n")
cat("  Calendario:        ", adjustedRandIndex(xi_post, datos$calendario),        "\n")
cat("  Colegio privado:   ", adjustedRandIndex(xi_post, datos$cole_privado),      "\n")
cat("  Colegio urbano:    ", adjustedRandIndex(xi_post, datos$cole_urbano),       "\n")
cat("  Género:            ", adjustedRandIndex(xi_post, datos$Genero_mujer),      "\n")
cat("  Etnia:             ", adjustedRandIndex(xi_post, datos$etnia),             "\n")
cat("  Trabaja:           ", adjustedRandIndex(xi_post, datos$trabaja_categoria), "\n")










