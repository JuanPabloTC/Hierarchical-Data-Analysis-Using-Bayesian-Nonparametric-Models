#@ Modelo 1

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

# Set working directory 
ruta <- setwd('path')

# Path for results
path_to_results <- paste0(ruta,'/Resultados/Modelo_1')

# Create dir (folder) if does not exist
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

# Departments and municipalities ID
g_dep <- rep(NA, m)
g_mcpio <- rep(NA, n_j)


#==============================================================
#@ =====         Y AN COVARIATES ORDERING                ======
#==============================================================

# ==== Prepare covariates matrix ====
# y       : y_{j,q} Scores
# x_jq    : x_{j,q} Individual level covariates
# z_jq    : z_{j,q} Municipal level covariates
# w_q     : w_{q}   Departamental level covariates

x_ijq_names <- c("fami_educacionmadre_modif", "computador", "internet", "etnia",
                 "libros_11_25", "libros_26_100", "libros_mas100", 
                 "estrato_1", "estrato_2", "estrato_3", "estrato_4", 
                 "estrato_5", "estrato_6", "Genero_mujer", "Calendario_A", 
                 "Calendario_B", "cole_privado", "trabaja_menos_de_10_horas",  
                 "trabaja__11_a_20_horas", "trabaja__21_a_30_horas", 
                 "trabaja_mas_de_30_horas", "cole_urbano")

z_jq_names <-  c("docenttotal_alumtotal", "RISK_VICTIM_2022", "Homi_x_100k_habitantes",
                 "porc_alumn_col_publico", "terrorismot", "hurto", "secuestros", "discapital",
                 "regalias_educa_per_capita")

w_q_names <-  c("PIB_percapita_DPTO", "proporcio_pob_rural", "perc_municipios_con_riesgo", 
                "Homicidios_ponderado_x_100k")

x_ijq_full <- as.matrix(subset(datos, select = x_ijq_names))
z_jq_full  <- as.matrix(subset(datos, select = z_jq_names))
w_q_full   <- as.matrix(subset(datos, select = w_q_names))

# ==== CREATE NESTED LISTS =====

Y      <- vector(mode = "list", length = m)
x_ijq  <- vector(mode = "list", length = m)
z_jq   <- vector(mode = "list", length = m)
w_q    <- vector(mode = "list", length = m)

g_dep   <- rep(NA, nrow(datos))
g_mcpio <- rep(NA, nrow(datos))

for (q in 1:m) {  #loop over deparments
  d_q <- depto_ids[q] # get department code
  # create boolean. TRUE for students from department q
  idx_dep <- datos$cole_cod_depto_ubicacion == d_q 
  g_dep[idx_dep] <- q
  
  # Identify and sort unique municipalities inside each department
  mcpio_ids <- sort(unique(datos$cole_cod_mcpio_ubicacion[idx_dep]))
  n_mcpios <- length(mcpio_ids)
  
  # Initialise lists over departments
  
  # Empty lists creation
  Y[[q]]     <- vector("list", length = n_mcpios)
  x_ijq[[q]] <- vector("list", length = n_mcpios)
  z_jq[[q]]  <- vector("list", length = n_mcpios)
  w_q[[q]]   <- vector("list", length = n_mcpios)
  
  for (j in 1:n_mcpios) { # loop for municipalities
    mcpio_j <- mcpio_ids[j]
    # Create boolean for students in both department and municipality
    idx_mcpio <- idx_dep & datos$cole_cod_mcpio_ubicacion == mcpio_j
    g_mcpio[idx_mcpio] <- j
    
    # Store observed variables and covariates in lists dept - mun
    Y[[q]][[j]]      <- datos$punt_global[idx_mcpio]
    x_ijq[[q]][[j]]  <- x_ijq_full[idx_mcpio, , drop = FALSE]
    z_jq[[q]][[j]]   <- z_jq_full[idx_mcpio, , drop = FALSE]
    w_q[[q]][[j]]    <- w_q_full[idx_mcpio, , drop = FALSE]
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
nj_x_mun <- estadisticos_municipios$nj_x_mun

# Sufficient statistics
yb_m <- estadisticos_departamentos$yb
yb_j <- estadisticos_municipios$yb
s2_m <- estadisticos_departamentos$s2
s2_j <- estadisticos_municipios$s2

#==============================================================
#@ =====      PARAMETERS AND HIPERPARAMETERS             ======
#==============================================================

# Number of parameters
numero_Betas_E = ncol(x_ijq_full)
numero_Betas_M = ncol(z_jq_full)
numero_Betas_D = ncol(w_q_full)
numero_Betas = 1 + numero_Betas_E + numero_Betas_M + numero_Betas_D
numero_parametros = 4 + m + n_j 

num_tot_parametros = numero_Betas + numero_parametros

### Non informative hiperparameters
nu_beta = 10
nu_E = 1.5
nu_M = 1.5
nu_D = 1.5
gamma_beta = 1.5
gamma_E = 1.5
gamma_M = 1.5
gamma_D = 1.5

# VARIANCE COMPONENT HIPERPARAMETERS
nu_kappa = 4 + (2 / 0.5^2)
alpha_kappa = 2 / (0.5^2)
beta_kappa = alpha_kappa / 250

#=================     DEFINE UNITARY PRIOR     =================

# Join covariate-level matrixes
X_full = cbind(x_ijq_full, z_jq_full, w_q_full)
lm_previa_unitaria = lm(y_vec ~ X_full)
betas_lm = coef(lm_previa_unitaria)

# Intercept
mu_beta = betas_lm[1]
# Observation stimators
mu_E <- betas_lm[ 2 : (1 + numero_Betas_E) ]

# Municipality stimators
start_M <- 2 + numero_Betas_E
end_M   <- 1 + numero_Betas_E + numero_Betas_M
mu_M     <- betas_lm[ start_M : end_M ]

# Department stimators
start_D <- end_M + 1
end_D   <- length(betas_lm)
mu_D     <- betas_lm[ start_D : end_D ]


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
    xtx_list[[q]][[j]] = t(x_ijq[[q]][[j]]) %*% x_ijq[[q]][[j]]
    ztz_list[[q]][[j]] = t(z_jq[[q]][[j]]) %*% z_jq[[q]][[j]]
    wtw_list[[q]][[j]] = t(w_q[[q]][[j]]) %*% w_q[[q]][[j]]
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
    mun_map[[length(mun_map) + 1L]] <- list(q = q, j = j, start = start_idx, end = end_idx, n = n_jq_local)
    mun_y[[length(mun_y) + 1L]]   <- as.numeric(Y[[q]][[j]])
    mun_x[[length(mun_x) + 1L]]   <- x_ijq[[q]][[j]]
    mun_z[[length(mun_z) + 1L]]   <- z_jq[[q]][[j]]
    mun_w[[length(mun_w) + 1L]]   <- w_q[[q]][[j]]
    mun_xtx[[length(mun_xtx) + 1L]] <- xtx_list[[q]][[j]]
    mun_ztz[[length(mun_ztz) + 1L]] <- ztz_list[[q]][[j]]
    mun_wtw[[length(mun_wtw) + 1L]] <- wtw_list[[q]][[j]]
    mun_sizes <- c(mun_sizes, n_jq_local)
    start_idx <- end_idx + 1L
  }
}
n_mun_total <- length(mun_map)  # should equal sum(sapply(Y,length)) 


# ===================================================================
#                       EXECUTE GIBBS SAMPLER
# ===================================================================

tictoc::tic()
Resultados <- MCMC1_optimized(25000, Y, mu_beta, mu_E, mu_M, mu_D, nu_beta,
                    nu_E, nu_M, nu_D, gamma_beta, gamma_E,
                    gamma_M, gamma_D, nu_kappa, alpha_kappa, beta_kappa)
tiempo = tictoc::toc()
# Retrieve data over iterations
burn_in    <- Resultados$info$burn_in
thin       <- Resultados$info$thin
B          <- Resultados$info$samples_stored
total_iter <- Resultados$info$total_iterations

Resultados$info$Tiempo = tiempo$callback_msg

# ===================================================================
#                       LOG-LIKELIHOOD TRACE PLOT
# ===================================================================

# Initialise a vector with log-likelihood computation
log_lik_vec <- numeric(B)


for (sample_count in seq_len(B)) {
  
  # Progress alerts
  if (sample_count %% 10 == 0 || sample_count == 1 || sample_count == B) {
    cat("Iteración", sample_count, "de", B, "\n")
    flush.console()
  }
  
  # Initialise the log-likelihood count 
  total_ll <- 0 # re starts loglikelihood count for each iteration
  for (q in 1:m) {
    for (j in 1:length(Y[[q]])) {
      n_j_q <- length(Y[[q]][[j]])
      # retrieve y_{j,q} and kappa^2_{j,q}
      y_vec      <- Y[[q]][[j]]
      k2_b       <- Resultados$kappa2_jq[[q]][[j]][sample_count]
      
      # Calculate the mean zeta_{j,q}
      
      # beta_E, beta_M, beta_D as vectors 
      beta_E_vec <- as.matrix(sapply(Resultados$beta_E, `[`, sample_count))
      beta_M_vec <- as.matrix(sapply(Resultados$beta_M, `[`, sample_count))
      beta_D_vec <- as.matrix(sapply(Resultados$beta_D, `[`, sample_count))
      
      Zeta_jq_b = Resultados$beta[sample_count]*rep(1, n_j_q) + x_ijq[[q]][[j]]%*%beta_E_vec + 
        z_jq[[q]][[j]]%*%beta_M_vec + w_q[[q]][[j]]%*%beta_D_vec
      
      
      # Sum loglikelihood of municipality, and then sum them for iteration b
      total_ll  <- total_ll + 
        sum(dnorm(x   = y_vec,
                  mean= Zeta_jq_b,
                  sd  = sqrt(k2_b),
                  log = TRUE))
    }
  }
  # Store loglikelihood as the vector already created
  log_lik_vec[sample_count] <- total_ll
} 

# Store loglikelihood in results
Resultados$log_likelihood <- log_lik_vec


# ===================================================================
#                 LOG-LIKELIHOOD CONVERGENCE PLOT
# ===================================================================

df_loglik <- data.frame(
  iteration = 1:B,
  log_likelihood = Resultados$log_likelihood
)

ll_trace <- ggplot(df_loglik, aes(x = iteration, y = log_likelihood)) +
  geom_line(color = "steelblue", alpha = 0.4) +
  geom_point(size = 1.1, alpha = 0.55) +
  labs(title = paste0("Gráfico de traza log-verosimilitud Modelo M1"),
       x = "Iteración MCMC",
       y = "Log-verosimilitud") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
x11()
print(ll_trace)



#===============================================================
#@========          INFORMATION CRITERIA               =========
#===============================================================



# ===================================================================
#                        EXECUTE WAIC
# ===================================================================

Resultados$WAIc_list <- compute_WAIC(Resultados, mun_map, mun_y, mun_x, mun_z, mun_w)



# ===================================================================
#                       SAVE OR LOAD THE CHAIN 
# ===================================================================
# Save chain
#save(Resultados, file = paste0(path_to_results, "/Resultados_Ifces_MCMC1.RData"))


# Load the chain 
load(file = paste0(path_to_results, "/Resultados_Ifces_MCMC1.RData"))


#===========================================================================
#                       EFFECTIVE SAMPLE SIZES    
#===========================================================================

ess_beta          = effectiveSize(Resultados$beta)
ess_sigma2_beta   = effectiveSize(Resultados$sigma2_beta)



kappa2_q_matrix <- do.call(cbind, Resultados$kappa2_q)
ess_kappa2_q <- effectiveSize(as.mcmc(kappa2_q_matrix))


beta_E_matrix <- do.call(cbind, Resultados$beta_E)
ess_beta_E <- effectiveSize(as.mcmc(beta_E_matrix))

beta_M_matrix <- do.call(cbind, Resultados$beta_M)
ess_beta_M <- effectiveSize(as.mcmc(beta_M_matrix))


beta_D_matrix <- do.call(cbind, Resultados$beta_D)
ess_beta_D <- effectiveSize(as.mcmc(beta_D_matrix))



ess_kappa2_jq <- lapply(Resultados$kappa2_jq, function(dep) { 
  #dep es una funcion que se refiere a cada elemento a nivel departamental
  lapply(dep, function(mun) {
    # Se refiere a los elementos del nivel municipal
    effectiveSize(as.mcmc(mun))
  })
})

ess_summary <- list(
  beta        = summary(ess_beta),
  sigma2_beta = summary(ess_sigma2_beta),
  #sigma2_E    = summary(ess_sigma2_E),
  #sigma2_M    = summary(ess_sigma2_M),
  #sigma2_D    = summary(ess_sigma2_D),
  beta_E      = summary(ess_beta_E),
  beta_M      = summary(ess_beta_M),
  beta_D      = summary(ess_beta_D),
  kappa2_q    = summary(ess_kappa2_q),
  kappa2_jq   = summary(unlist(ess_kappa2_jq))
)



#=========================================
#    INFERENCE
#=========================================

#--------------------------------------------------------------------------------
# Creates lists: beta_int_summary, beta_E_summary, beta_M_summary, beta_D_summary
# 1. mean 
# 2. q0025
# 3. q50
# 4. q975
#--------------------------------------------------------------------------------

Resultados <- compute_posterior_summaries_M1(Resultados)


#metricas_M1 <- compute_test_metrics_M1(
#  Resultados    = Resultados,
#  datos_testeo  = datos,
#  x_ijq_names   = x_ijq_names,
#  z_jq_names    = z_jq_names,
#  w_q_names     = w_q_names
#)




# ===================================================================
#                COEFICIENTS SIGNIFICANCE PLOTS
# ===================================================================


# ------------------------------------------------------------------
# 1. Build df for quantiles from Results.
# ------------------------------------------------------------------

df_quantiles_betas <- rbind(
  data.frame(
      param = names(Resultados$beta_E_summary$mean),
    q025  = as.numeric(Resultados$beta_E_summary$q025),
    q50   = as.numeric(Resultados$beta_E_summary$q50),
    q975  = as.numeric(Resultados$beta_E_summary$q975),
    nivel = "Individual (β_E)",
    stringsAsFactors = FALSE
  ),
  data.frame(
    param = names(Resultados$beta_M_summary$mean),
    q025  = as.numeric(Resultados$beta_M_summary$q025),
    q50   = as.numeric(Resultados$beta_M_summary$q50),
    q975  = as.numeric(Resultados$beta_M_summary$q975),
    nivel = "Municipal (β_M)",
    stringsAsFactors = FALSE
  ),
  data.frame(
    param = names(Resultados$beta_D_summary$mean),
    q025  = as.numeric(Resultados$beta_D_summary$q025),
    q50   = as.numeric(Resultados$beta_D_summary$q50),
    q975  = as.numeric(Resultados$beta_D_summary$q975),
    nivel = "Departamental (β_D)",
    stringsAsFactors = FALSE
  )
)

# Flag of significance: True if the IC 95% does not cross 0
df_quantiles_betas$signif <- with(df_quantiles_betas,
                                  ifelse(q025 * q975 > 0, "Sí", "No"))

# Ordered factors: Hierarchical order and within levels
# Order paramiters with magnitude (q50)
df_quantiles_betas$nivel <- factor(
  df_quantiles_betas$nivel,
  levels = c("Individual (β_E)", "Municipal (β_M)", "Departamental (β_D)")
)

df_quantiles_betas <- df_quantiles_betas[
  order(df_quantiles_betas$nivel, df_quantiles_betas$q50), ]

df_quantiles_betas$param <- factor(df_quantiles_betas$param,
                                   levels = df_quantiles_betas$param)


# ------------------------------------------------------------------
# Mapping Name -> Label
# ------------------------------------------------------------------
var_labels_map <- c(
  # β_E  — Individual
  "fami_educacionmadre_modif"   = "Educación madre",
  "computador"                  = "Computador (hogar)",
  "internet"                    = "Internet (hogar)",
  "etnia"                       = "Etnia (minoría)",
  "libros_11_25"                = "Libros en hogar (11-25)",
  "libros_26_100"               = "Libros en hogar (26-100)",
  "libros_mas100"               = "Libros en hogar (>100)",
  "estrato_1"                   = "Estrato 1",
  "estrato_2"                   = "Estrato 2",
  "estrato_3"                   = "Estrato 3",
  "estrato_4"                   = "Estrato 4",
  "estrato_5"                   = "Estrato 5",
  "estrato_6"                   = "Estrato 6",
  "Genero_mujer"                = "Género (mujer)",
  "Calendario_A"                = "Calendario A",
  "Calendario_B"                = "Calendario B",
  "cole_privado"                = "Colegio privado",
  "trabaja_menos_de_10_horas"   = "Trabaja <10 h/sem",
  "trabaja__11_a_20_horas"      = "Trabaja 11-20 h/sem",
  "trabaja__21_a_30_horas"      = "Trabaja 21-30 h/sem",
  "trabaja_mas_de_30_horas"     = "Trabaja >30 h/sem",
  # β_M — Municipal
  "cole_urbano"                 = "Colegio urbano",
  "docenttotal_alumtotal"       = "Razón alumnos/docente",
  "RISK_VICTIM_2022"            = "Riesgo victimización",
  "Homi_x_100k_habitantes"      = "Homicidios por 100k hab.",
  "porc_alumn_col_publico"      = "% alumnos col. público",
  "terrorismot"                 = "Terrorismo",
  "hurto"                       = "Hurto",
  "secuestros"                  = "Secuestros",
  # β_D — Departamental
  "discapital"                  = "Distancia a capital",
  "regalias_educa_per_capita"   = "Regalías educ. per cápita",
  "PIB_percapita_DPTO"          = "PIB per cápita (dpto.)",
  "proporcio_pob_rural"         = "Proporción pob. rural",
  "perc_municipios_con_riesgo"  = "% municipios en riesgo",
  "Homicidios_ponderado_x_100k" = "Homicidios pond. x 100k"
)


# ------------------------------------------------------------------
# 2. Plot
# ------------------------------------------------------------------

paleta_niveles <- c(
  "Individual"    = "#1F77B4",
  "Municipal"     = "#2CA02C",
  "Departamental" = "#D62728"
)

df_quantiles_betas$nivel <- factor(
  gsub(" \\(.*\\)$", "", as.character(df_quantiles_betas$nivel)),
  levels = c("Individual", "Municipal", "Departamental")
)

x11()

p_betas <- ggplot(df_quantiles_betas,
                  aes(x = q50, y = param, xmin = q025, xmax = q975,
                      color = nivel, shape = signif, linetype = signif)) +
  # reference line in 0
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  # CI and mean
  geom_pointrange(size = 0.6, linewidth = 0.75) +
  # Color por nivel jerárquico
  scale_color_manual(values = paleta_niveles, name = "Nivel") +
  # Significance: solid circle vs empty circle (no fill) 
  scale_shape_manual(values = c("Sí" = 19, "No" = 1),
                     name = "IC 95%",
                     labels = c("Sí" = "no incluye 0", "No" = "incluye 0")) +
  scale_linetype_manual(values = c("Sí" = "solid", "No" = "dotted"),
                        name = "IC 95%",
                        labels = c("Sí" = "no incluye 0", "No" = "incluye 0")) +
  labs(
    #title    = "Efectos y significancia por nivel jerárquico",
    #subtitle = "Media posterior con intervalo de credibilidad al 95%",
    x        = "Efecto estimado",
    y        = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title         = element_text(size = 16, face = "bold"),
    plot.subtitle      = element_text(size = 12, color = "grey30",
                                      margin = margin(b = 10)),
    axis.title.x       = element_text(size = 13, margin = margin(t = 8)),
    axis.text.y        = element_text(size = 11),
    axis.text.x        = element_text(size = 11),
    legend.position    = "top",
    legend.box         = "horizontal",
    legend.title       = element_text(size = 11, face = "bold"),
    legend.text        = element_text(size = 10),
    legend.margin      = margin(0, 0, 0, 0),
    legend.box.spacing = unit(4, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "grey90"),
    plot.margin        = margin(15, 20, 15, 15)
  ) +
  guides(
    color    = guide_legend(order = 1,
                            override.aes = list(shape = 15, linetype = 0, size = 0.8)),
    shape    = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  ) +
  scale_y_discrete(labels = var_labels_map) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)))

print(p_betas)


# ============================================================================
# STANDARDIZED COEFFICIENTS - MODEL M1
# ============================================================================
# Computation of standardized coefficients:
#   beta_std_j = beta_j * sd(x_j)
#
# Interpretation: expected effect on the Saber 11 score associated with a
# one-standard-deviation change in covariate x_j.
#
# For binary dummies, sd(x_j) = sqrt(p*(1-p)) where p = proportion of 1s.
# The interpretation becomes: effect associated with a typical change
# (one sd) in the prevalence of the category in the population.
# ============================================================================

# --- 1. Compute empirical standard deviations of the covariates ---

sd_x_E <- apply(x_ijq_full, 2, sd)   # sd of individual-level covariates
sd_x_M <- apply(z_jq_full,  2, sd)   # sd of municipal-level covariates
sd_x_D <- apply(w_q_full,   2, sd)   # sd of departmental-level covariates

# --- 2. Compute standardized coefficients ---
# For each level, compute the standardized mean and its percentiles

# === Individual level (beta_E) ===

beta_E_mat <- do.call(cbind, Resultados$beta_E)   # [B x p_E]

# Multiply each column of the chain by the corresponding sd
# sweep(matrix, MARGIN=2, vector, FUN="*") multiplies column-wise
beta_E_std_mat <- sweep(beta_E_mat, MARGIN = 2, STATS = sd_x_E, FUN = "*")

beta_E_std_summary <- list(
  mean = colMeans(beta_E_std_mat),
  q025 = apply(beta_E_std_mat, 2, quantile, probs = 0.025),
  q50  = apply(beta_E_std_mat, 2, quantile, probs = 0.500),
  q975 = apply(beta_E_std_mat, 2, quantile, probs = 0.975)
)

# === Municipal level (beta_M) ===

beta_M_mat <- do.call(cbind, Resultados$beta_M)   # [B x p_M]
beta_M_std_mat <- sweep(beta_M_mat, MARGIN = 2, STATS = sd_x_M, FUN = "*")

beta_M_std_summary <- list(
  mean = colMeans(beta_M_std_mat),
  q025 = apply(beta_M_std_mat, 2, quantile, probs = 0.025),
  q50  = apply(beta_M_std_mat, 2, quantile, probs = 0.500),
  q975 = apply(beta_M_std_mat, 2, quantile, probs = 0.975)
)

# === Departmental level (beta_D) ===

beta_D_mat <- do.call(cbind, Resultados$beta_D)   # [B x p_D]
beta_D_std_mat <- sweep(beta_D_mat, MARGIN = 2, STATS = sd_x_D, FUN = "*")

beta_D_std_summary <- list(
  mean = colMeans(beta_D_std_mat),
  q025 = apply(beta_D_std_mat, 2, quantile, probs = 0.025),
  q50  = apply(beta_D_std_mat, 2, quantile, probs = 0.500),
  q975 = apply(beta_D_std_mat, 2, quantile, probs = 0.975)
)

# --- 3. Store in Resultados ---

Resultados$beta_E_std_summary <- beta_E_std_summary
Resultados$beta_M_std_summary <- beta_M_std_summary
Resultados$beta_D_std_summary <- beta_D_std_summary

# --- 4. Display results ---

cat("\n=== Standardized coefficients - Individual level ===\n")
df_E <- data.frame(
  variable   = names(Resultados$beta_E_summary$mean),
  beta_raw   = Resultados$beta_E_summary$mean,
  sd_x       = sd_x_E,
  beta_std   = Resultados$beta_E_std_summary$mean
)
df_E[, -1] <- round(df_E[, -1], 3)
print(df_E)

cat("\n=== Standardized coefficients - Municipal level ===\n")
df_M <- data.frame(
  variable   = names(Resultados$beta_M_summary$mean),
  beta_raw   = Resultados$beta_M_summary$mean,
  sd_x       = sd_x_M,
  beta_std   = Resultados$beta_M_std_summary$mean
)
df_M[, -1] <- round(df_M[, -1], 4)
print(df_M)

cat("\n=== Standardized coefficients - Departmental level ===\n")
df_D <- data.frame(
  variable   = names(Resultados$beta_D_summary$mean),
  beta_raw   = Resultados$beta_D_summary$mean,
  sd_x       = sd_x_D,
  beta_std   = Resultados$beta_D_std_summary$mean
)
df_D[, -1] <- round(df_D[, -1], 4)
print(df_D)

# --- 5. Sort and display covariates by relative importance ---
# Useful to validate the analysis in the body of the document

all_std <- data.frame(
  variable = c(
    names(Resultados$beta_E_std_summary$mean),
    names(Resultados$beta_M_std_summary$mean),
    names(Resultados$beta_D_std_summary$mean)
  ),
  level = c(
    rep("Individual",    length(Resultados$beta_E_std_summary$mean)),
    rep("Municipal",     length(Resultados$beta_M_std_summary$mean)),
    rep("Departmental",  length(Resultados$beta_D_std_summary$mean))
  ),
  beta_std = c(
    Resultados$beta_E_std_summary$mean,
    Resultados$beta_M_std_summary$mean,
    Resultados$beta_D_std_summary$mean
  )
)
all_std$abs_beta_std <- abs(all_std$beta_std)
all_std <- all_std[order(-all_std$abs_beta_std), ]

cat("\n=== Ranking of relative importance by absolute standardized magnitude ===\n")
print(round_df <- data.frame(
  variable = all_std$variable,
  level    = all_std$level,
  beta_std = round(all_std$beta_std, 3)
))


