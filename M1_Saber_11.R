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
ruta <- setwd('D:/Actualizado/Maestria estadistica/Tesis/Modelos BNP/15-05-26')

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
#












































































































#=====================     DIC   =====================

# 1.  Recuperar las muestras de logverosimilitud calculada en GIBBS
LL_samples = Resultados$log_likelihood

# 2. Calcular los valores esperados E[zeta_{jq}] y E[kappa2_{jq}]

# a) Hallar la media posterior de Zeta_jq[[q]][[j]] para cada estudiante i

# Se inicializa la lista
zeta_hat = vector("list", length = m)

for(q in seq_len(m)) {
  # numero de municipios en departamento q
  n_muns = length(Y[[q]])
  # Se inicializa las sublistas por depto
  zeta_hat[[q]] = vector("list", length = n_muns)
  cat("Procesando departamento", q, "of", m, "...\n")
  for(j in 1:n_muns) {
    n_j_q <- length(Y[[q]][[j]])
    
    # Pre‑alocamos la matriz n_j_q x B donde pondremos los zeta_{j,q} para las B iteraciones
    mat <- matrix(NA_real_, nrow = n_j_q, ncol = B)
    
    # Calcular zeta_{j,q} para las B iteraciones
    
    # Iterar sobre todas las iteraciones del muestreador
    for (sample_count in seq_len(B)) {
      # Calculat zeta_{j,q}
      
      #Formar beta_E, beta_M, beta_D como vectores
      beta0       = Resultados$beta[sample_count]
      beta_E_vec  = as.matrix(sapply(Resultados$beta_E, `[`, sample_count))
      beta_M_vec  = as.matrix(sapply(Resultados$beta_M, `[`, sample_count))
      beta_D_vec  = as.matrix(sapply(Resultados$beta_D, `[`, sample_count))
      
      # calcular zeta_{j,q} en la iteracion sample_count
      Zeta_jq_b = beta0*rep(1, n_j_q) + x_ijq[[q]][[j]]%*%beta_E_vec + 
        z_jq[[q]][[j]]%*%beta_M_vec + w_q[[q]][[j]]%*%beta_D_vec
      
      # Guardamos como columna sample_count
      mat[, sample_count] <- as.vector(Zeta_jq_b)
      
    }
    zeta_hat[[q]][[j]] <- rowMeans(mat)
  }
}

# b) Hallar media posterior de kappa2_jq[[q]][[j]] para cada municipio

# Se inicializa la lista
kappa2_hat = vector("list", length = m)
for(q in 1:m) {
  # numero de municipios en el dpto
  n_muns = length(Y[[q]])
  # se unicializa las sublistas
  kappa2_hat[[q]] = numeric(n_muns)
  for(j in 1:n_muns) {
    # se calcula la media de las B iteraciones por departamento
    kappa2_hat[[q]][j] = mean(Resultados$kappa2_jq[[q]][[j]])
  }
}

# 3) Se evalua la log-verosimilitud con las medias posteriores halladas en 2.

# Con los valores de 2. se calcula la log-verosimilitud se calculara como 
#\sum_{q}\sum_{j}\sum_{i} log N(y_{i,j,q} | \hat{\zeta_{i,j,q}}, \hat{\kappa^2_{j,q}})

# Se inicializa el valor de la log-verosimilitud
lpy_hat <- 0 #log-verosimilitud 
for(q in 1:m) {
  for(j in 1:length(Y[[q]])) {
    # observaciones puntaje en el municipio j del depto q
    y_vec_jq = Y[[q]][[j]]
    # vector de medias de estudiantes i en Zeta_jq de n_jq x 1
    mu_vec = zeta_hat[[q]][[j]]
    # Varianza esperada del municipio j (escalar)
    sigma_jq = sqrt(kappa2_hat[[q]][j])
    # Se calcula la log-verosimilitud (se esta sumando en lpy_hat)
    lpy_hat  = lpy_hat + sum(dnorm(y_vec_jq, mean = mu_vec, sd   = sigma_jq,
                                   log  = TRUE)) #log = TRUE, para que sea escala log
  }
}

# 4) Numero efectivo de parmaetros

# log-verosimilitud posterior - media de log-verosimilitud de las muestras de Gibbs
p_DIC = 2 * (lpy_hat - mean(LL_samples))

# 5.) Calculo de DIC.

DIC = (-2 * lpy_hat) + (2*p_DIC)



# ===================================================================
#                       GUARDAR RESULTADOS
# ===================================================================

# save(Resultados, file = paste0(ruta, "/Resultados/Modelo_2/Resultados_MCMC2.RData"))

#=====================     Cargar Resultados   =====================

if (run_gibbs == 0) {
  load(paste0(ruta, "/Resultados/Modelo 2/Nuevos_Resultados_sin_error_MCM1.RData"))
}






#===============================================================
#@========                Valores PPP                  =========
#===============================================================

# Extraer parámetros de los resultados
B <- length(Resultados$sigma2_beta)  # Número de iteraciones MCMC
n_dept <- length(Y)  # Número de departamentos (33)

# Crear vector con todos los datos observados originales
Y_obs <- c()
for (q in 1:n_dept) {
  for (j in 1:length(Y[[q]])) {
    Y_obs <- c(Y_obs, Y[[q]][[j]])
  }
}

# Convertir a vector numérico
Y_obs_vector <- unlist(Y_obs)
Y_obs_vector <- as.numeric(Y_obs_vector)

# Calcular estadísticos observados
media_obs <- mean(Y_obs_vector, na.rm = TRUE)
sd_obs <- sd(Y_obs_vector, na.rm = TRUE)

cat("Estadísticos observados:\n")
cat("Media observada:", round(media_obs, 4), "\n")
cat("Desviación estándar observada:", round(sd_obs, 4), "\n\n")

# Matriz para almacenar estadísticos de prueba
TS <- matrix(NA, B, 2)
colnames(TS) <- c("media", "desv_std")

cat("Generando datos replicados y calculando estadísticos...\n")
cat("NOTA: Calculando zeta sobre la marcha para ahorrar memoria\n\n")
set.seed(123)

# Generar datos replicados para cada iteración MCMC - SIN ALMACENAR ZETA_FULL
for (b in 1:B) {
  if (b %% 500 == 0) cat("Iteración", b, "de", B, "\n")
  
  # Extraer parámetros de la iteración b UNA SOLA VEZ
  beta0 <- Resultados$beta[b]
  beta_E_vec <- as.matrix(sapply(Resultados$beta_E, `[`, b))
  beta_M_vec <- as.matrix(sapply(Resultados$beta_M, `[`, b))
  beta_D_vec <- as.matrix(sapply(Resultados$beta_D, `[`, b))
  
  # Vector para almacenar datos replicados de la iteración b
  Y_rep <- c()
  
  # Generar datos replicados para cada departamento
  for (q in 1:n_dept) {
    n_munic_q <- length(Y[[q]])
    
    # Para cada municipio en el departamento q
    for (j in 1:n_munic_q) {
      n_students_jq <- length(Y[[q]][[j]])
      
      # CALCULAR ZETA SOLO PARA ESTA ITERACIÓN Y ESTE MUNICIPIO
      Zeta_jq_b <- beta0 * rep(1, n_students_jq) + 
        x_ijq[[q]][[j]] %*% beta_E_vec + 
        z_jq[[q]][[j]] %*% beta_M_vec + 
        w_q[[q]][[j]] %*% beta_D_vec
      
      # Varianza específica del municipio para esta iteración
      kappa2_jq <- Resultados$kappa2_jq[[q]][[j]][b]
      
      # Generar puntajes replicados para todos los estudiantes de este municipio
      Y_rep_jq <- rnorm(n_students_jq, 
                        mean = as.vector(Zeta_jq_b), 
                        sd = sqrt(kappa2_jq))
      
      # Agregar al vector de datos replicados
      Y_rep <- c(Y_rep, Y_rep_jq)
    }
  }
  
  # Calcular estadísticos de los datos replicados para la iteración b
  TS[b, 1] <- mean(Y_rep, na.rm = TRUE)    # Media
  TS[b, 2] <- sd(Y_rep, na.rm = TRUE)      # Desviación estándar
}

cat("\nCálculo completado. Calculando valores ppp...\n")

# Calcular valores p-posteriores predictivos (ppp)
ppp_media <- mean(TS[, 1] >= media_obs, na.rm = TRUE)
ppp_sd <- mean(TS[, 2] >= sd_obs, na.rm = TRUE)

# Interpretación
cat("INTERPRETACIÓN:\n")
cat("- ppp cerca de 0.5 indica buen ajuste del modelo\n")
cat("- ppp < 0.05 o ppp > 0.95 sugiere falta de ajuste\n\n")

# Evaluación del ajuste
ajuste_media <- ifelse(ppp_media > 0.05 & ppp_media < 0.95, "Bueno", "Cuestionable")
ajuste_sd <- ifelse(ppp_sd > 0.05 & ppp_sd < 0.95, "Bueno", "Cuestionable")

cat("EVALUACIÓN DEL AJUSTE:\n")
cat("Media: ", ajuste_media, " (ppp = ", round(ppp_media, 4), ")\n", sep = "")
cat("Desv. estándar: ", ajuste_sd, " (ppp = ", round(ppp_sd, 4), ")\n\n", sep = "")

# Gráficos de diagnóstico
par(mfrow = c(1, 2))

# Histograma para la media
hist(TS[, 1], main = "Distribución predictiva\nposterior de la media - Modelo 2", 
     xlab = "Media", col = "grey80", breaks = 30, 
     freq = FALSE, probability = TRUE)
abline(v = media_obs, col = "red", lwd = 3)
legend("topright", 
       legend = c(paste("Valor observado =", round(media_obs, 3)), paste("ppp =", round(ppp_media, 3))), 
       col = c("red", "black"), lty = c(1, NA), lwd = c(3, NA),
       cex = 0.8)

# Histograma para la desviación estándar
hist(TS[, 2], main = "Distribución predictiva\nposterior de la desv. estándar - Modelo 2", 
     xlab = "Desviación estándar", col = "grey80", breaks = 30,
     freq = FALSE, probability = TRUE)
abline(v = sd_obs, col = "red", lwd = 3)
legend("topright", 
       legend = c(paste("Valor observado =", round(sd_obs, 3)), paste("ppp =", round(ppp_sd, 3))), 
       col = c("red", "black"), lty = c(1, NA), lwd = c(3, NA),
       cex = 0.8)

par(mfrow = c(1, 1))

# Tabla resumen
resultados_ppp <- data.frame(
  Estadistico = c("Media", "Desviación estándar"),
  Valor_observado = round(c(media_obs, sd_obs), 4),
  ppp = round(c(ppp_media, ppp_sd), 4),
  Ajuste = c(ajuste_media, ajuste_sd),
  stringsAsFactors = FALSE
)

cat("TABLA RESUMEN:\n")
print(resultados_ppp, row.names = FALSE)


#===============================================================
#@========         QUANTILES 0.025, 0.5, 0.975         =========
#===============================================================
numero_Betas_E = ncol(x_ijq_full)
numero_Betas_M = ncol(z_jq_full)
numero_Betas_D = ncol(w_q_full)


# Almacenar cuantiles  (0.025, 0.5, 0.975)

# Copiamos la estructura de resultados
quantiles_parametros <- Resultados

# Eliminamos las que no vamos a usar
quantiles_parametros$info          <- NULL
quantiles_parametros$log_likelihood <- NULL
quantiles_parametros$proposal_sd_alpha <- NULL


# Para los kappa2_jq y kappa2_q 
for (q in 1:m) {
  cuantiles_kappa_q = quantile(
    Resultados$kappa2_q[[q]],
    probs = c(0.025, 0.5, 0.975), names = FALSE
  )
  quantiles_parametros$kappa2_q[[q]] = cuantiles_kappa_q
  for (j in seq_along(Resultados$kappa2_jq[[q]])) {
    cuantiles_kappa_jq = quantile(
      Resultados$kappa2_jq[[q]][[j]],
      probs = c(0.025, 0.5, 0.975), names = FALSE
    )
    quantiles_parametros$kappa2_jq[[q]][[j]] = cuantiles_kappa_jq
  }
}

# Para los betas
cuantiles_beta = quantile(
  Resultados$beta,
  probs = c(0.025, 0.5, 0.975), names = FALSE
)
quantiles_parametros$beta = cuantiles_beta

for (p in 1:numero_Betas_E) {
  cuantiles_beta_E <- quantile(
    Resultados$beta_E[[p]],
    probs = c(0.025, 0.5, 0.975),
    names = FALSE
  )
  quantiles_parametros$beta_E[[p]] <- cuantiles_beta_E
}

for (p in 1:numero_Betas_M) {
  cuantiles_beta_M <- quantile(
    Resultados$beta_M[[p]],
    probs = c(0.025, 0.5, 0.975),
    names = FALSE
  )
  quantiles_parametros$beta_M[[p]] <- cuantiles_beta_M
}

for (p in 1:numero_Betas_D) {
  cuantiles_beta_D <- quantile(
    Resultados$beta_D[[p]],
    probs = c(0.025, 0.5, 0.975),
    names = FALSE
  )
  quantiles_parametros$beta_D[[p]] <- cuantiles_beta_D
}

# Para los sigma2
cuantiles_sigma2_beta = quantile(
  Resultados$sigma2_beta,
  probs = c(0.025, 0.5, 0.975), names = FALSE
)
quantiles_parametros$sigma2_beta = cuantiles_sigma2_beta

cuantiles_sigma2_E = quantile(
  Resultados$sigma2_E,
  probs = c(0.025, 0.5, 0.975), names = FALSE
)
quantiles_parametros$sigma2_E = cuantiles_sigma2_E

cuantiles_sigma2_M = quantile(
  Resultados$sigma2_M,
  probs = c(0.025, 0.5, 0.975), names = FALSE
)
quantiles_parametros$sigma2_M = cuantiles_sigma2_M

cuantiles_sigma2_D = quantile(
  Resultados$sigma2_D,
  probs = c(0.025, 0.5, 0.975), names = FALSE
)
quantiles_parametros$sigma2_D = cuantiles_sigma2_D



#=====================     Cuantiles betas como dataframe   =====================


# Extraer cada bloque de cuantiles
qe <- quantiles_parametros$beta_E   # lista de length p_e, cada elemento num[3]
qm <- quantiles_parametros$beta_M   # lista de length p_m, cada elemento num[3]
qd <- quantiles_parametros$beta_D   # lista de length p_d, cada elemento num[3]

# Convertir cada uno en data.frame
df_e <- as.data.frame(do.call(rbind, qe))
df_m <- as.data.frame(do.call(rbind, qm))
df_d <- as.data.frame(do.call(rbind, qd))

# Poner nombres de cuantiles a las colmnas
colnames(df_e) <- c("q025", "q50", "q975")
colnames(df_m) <- c("q025", "q50", "q975")
colnames(df_d) <- c("q025", "q50", "q975")

# Apilar todo en uno solo
df_quantiles_betas <- rbind(df_e, df_m, df_d)


#=====================     Grafico significancia betas   =====================


# Preparar el dataframe para graficar con ggplot2
df_quantiles_betas$param <- rownames(df_quantiles_betas)

# Flag de significancia: TRUE si el intervalo no cruza 0. 
#Se evalua que la multiplicación de los cuantiles tenga mismo signo, si no es porque
#cruza por el cero
df_quantiles_betas$signif <- with(df_quantiles_betas, ifelse(q025 * q975 > 0, "sí", "no"))

# Graficar
x11()
ggplot(df_quantiles_betas, aes(x = param, y = q50, 
                               ymin = q025, ymax = q975,
                               color = signif)) +
  geom_pointrange(size = 0.45) +
  scale_color_manual(
    values = c("sí" = "blue",    # no toca cero
               "no" = "black")   # toca cero
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey50") +
  labs(x = "Parámetro", 
       y = "Efecto (q50 ± intervalo 95%)",
       color = "Significativo al 95%") +
  theme_minimal(base_size = 16) +    # Tamaño base más grande
  theme(
    axis.title   = element_text(size = 16),  # Títulos de ejes
    axis.text    = element_text(size = 14),  # Texto de ticks
    legend.title = element_text(size = 16),  # Título de la leyenda
    legend.text  = element_text(size = 12),  # Texto de la leyenda
    plot.title   = element_text(size = 14, face = "bold"),  # Si añades un título
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )


#===============================================================
#@========       Medias (Estimadores Bayesianos)       =========
#===============================================================

# Copiamos la estructura de resultados
Estimadores_Bayesianos <- Resultados

# Eliminamos las que no vamos a usar
Estimadores_Bayesianos$info             <- NULL
Estimadores_Bayesianos$log_likelihood   <- NULL
Estimadores_Bayesianos$proposal_sd_alpha <- NULL

# Para los kappa2_q y kappa2_jq 
for (q in 1:m) {
  # media de kappa2_q en el departamento q
  media_kappa_q <- mean(Resultados$kappa2_q[[q]], na.rm = TRUE)
  Estimadores_Bayesianos$kappa2_q[[q]] <- media_kappa_q
  
  # media de kappa2_jq en cada municipio j de q
  for (j in seq_along(Resultados$kappa2_jq[[q]])) {
    media_kappa_jq <- mean(Resultados$kappa2_jq[[q]][[j]], na.rm = TRUE)
    Estimadores_Bayesianos$kappa2_jq[[q]][[j]] <- media_kappa_jq
  }
}

# Para el beta global
media_beta <- mean(Resultados$beta, na.rm = TRUE)
Estimadores_Bayesianos$beta <- media_beta

# Para los beta_E
for (p in 1:numero_Betas_E) {
  media_beta_E <- mean(Resultados$beta_E[[p]], na.rm = TRUE)
  Estimadores_Bayesianos$beta_E[[p]] <- media_beta_E
}

# Para los beta_M
for (p in 1:numero_Betas_M) {
  media_beta_M <- mean(Resultados$beta_M[[p]], na.rm = TRUE)
  Estimadores_Bayesianos$beta_M[[p]] <- media_beta_M
}

# Para los beta_D
for (p in 1:numero_Betas_D) {
  media_beta_D <- mean(Resultados$beta_D[[p]], na.rm = TRUE)
  Estimadores_Bayesianos$beta_D[[p]] <- media_beta_D
}

# Para las varianzas sigma2
media_sigma2_beta <- mean(Resultados$sigma2_beta, na.rm = TRUE)
Estimadores_Bayesianos$sigma2_beta <- media_sigma2_beta

media_sigma2_E <- mean(Resultados$sigma2_E, na.rm = TRUE)
Estimadores_Bayesianos$sigma2_E <- media_sigma2_E

media_sigma2_M <- mean(Resultados$sigma2_M, na.rm = TRUE)
Estimadores_Bayesianos$sigma2_M <- media_sigma2_M

media_sigma2_D <- mean(Resultados$sigma2_D, na.rm = TRUE)
Estimadores_Bayesianos$sigma2_D <- media_sigma2_D


#===========================================================================
#@========         BONDAD DE AJUSTE / HABILIDAD PREDICTIVA         =========
#===========================================================================


#=====================     Cargar datos testeo   =====================

datos_testeo <- read.csv(paste0(ruta,'/Datos_Procesados/Examen_Saber_11_2022_2_Datos_test.csv'))
datos_testeo <- datos_testeo[!(datos_testeo$estu_cod_reside_depto == 99999), ]


#datos_testeo <- datos_testeo[, !names(datos_testeo) %in% c("homicidios", "codprovincia"
#                                      , "ano")]  # drop the column
anyna(datos_testeo)
#datos_testeo <- na.omit(datos_testeo)                             # drop rows with any NA
#datos_testeo <- datos_testeo[!(datos_testeo$cole_cod_depto_ubicacion == 68 & datos_testeo$cole_cod_mcpio_ubicacion == 68264), ]

#=====================     Cargar datos testeo   =====================


depto_ids_testeo <- sort(unique(datos_testeo$cole_cod_depto_ubicacion))
mcpio_ids_testeo <- sort(unique(datos_testeo$cole_cod_mcpio_ubicacion))

# m : numero de departamentos
m <- length(depto_ids_testeo) # m : número de grupos (departamentos)
# n_j : numero de municipios
n_j <- length(mcpio_ids_testeo) # j : número de sub-grupos (municipios)
# n: numero estudiantes
n_jq <- nrow(datos_testeo) # n: número de individuos (estudiantes)


y_test <- datos_testeo$punt_global
# y como vector
y_test_vec <- as.matrix(y_test)

g_dep <- rep(NA, m)
g_mcpio <- rep(NA, n_j)


#=====================     Ordenar datos testeo   =====================

# ==== Preparar matrices de covariables ====
# y: y_{j,q} Puntajes
# x_jq: x_{j,q} Variables a nivel individual
# z_jq: z_{j,q} Variables a nivel municipal
# w_q: w_{q} Variables a nivel departamental


x_ijq_test_full <- as.matrix(subset(datos_testeo, select = x_ijq_names))
z_jq_test_full  <- as.matrix(subset(datos_testeo, select = z_jq_names))
w_q_test_full   <- as.matrix(subset(datos_testeo, select = w_q_names))

# ==== Crear listas anidadas =====
Y_test      <- vector(mode = "list", length = m)
x_ijq_test  <- vector(mode = "list", length = m)
z_jq_test   <- vector(mode = "list", length = m)
w_q_test    <- vector(mode = "list", length = m)

g_dep   <- rep(NA, nrow(datos_testeo))
g_mcpio <- rep(NA, nrow(datos_testeo))

for (q in 1:m) {
  d_q <- depto_ids_testeo[q]
  idx_dep <- datos_testeo$cole_cod_depto_ubicacion == d_q
  g_dep[idx_dep] <- q
  
  # Municipios únicos dentro del departamento q
  mcpio_ids_testeo <- sort(unique(datos_testeo$cole_cod_mcpio_ubicacion[idx_dep]))
  n_mcpios <- length(mcpio_ids_testeo)
  
  # Inicializar listas por departamento
  Y_test[[q]]          <- vector("list", length = n_mcpios)
  x_ijq_test[[q]] <- vector("list", length = n_mcpios)
  z_jq_test[[q]]  <- vector("list", length = n_mcpios)
  w_q_test[[q]]   <- vector("list", length = n_mcpios)
  
  for (j in 1:n_mcpios) {
    mcpio_j <- mcpio_ids_testeo[j]
    idx_mcpio <- idx_dep & datos_testeo$cole_cod_mcpio_ubicacion == mcpio_j
    g_mcpio[idx_mcpio] <- j
    
    # Guardar variables observadas y covariables
    Y_test[[q]][[j]]          <- datos_testeo$punt_global[idx_mcpio]
    x_ijq_test[[q]][[j]] <- x_ijq_test_full[idx_mcpio, , drop = FALSE]
    z_jq_test[[q]][[j]]  <- z_jq_test_full[idx_mcpio, , drop = FALSE]
    w_q_test[[q]][[j]]   <- w_q_test_full[idx_mcpio, , drop = FALSE]
  }
}


#=====================     PREDICCION, MSE, MAE Y R2   =====================

# 1) Extraer estimadores puntuales
b0   <- unlist(Estimadores_Bayesianos$beta)
bE   <- as.matrix(unlist(Estimadores_Bayesianos$beta_E))   
bM   <- as.matrix(unlist(Estimadores_Bayesianos$beta_M))
bD   <- as.matrix(unlist(Estimadores_Bayesianos$beta_D))


# 2) Loop sobre departamentos q y municipios j
m <- length(Y_test)

zeta_list <- vector("list", m)

for (q in 1:m) {
  cat(sprintf("Iniciando iteración q = %d de %d\n", q, m))
  flush.console()
  
  Jk <- length(Y_test[[q]])
  for (j in seq_len(Jk)) {
    # datos observados
    y_vec_jq = Y_test[[q]][[j]]               # vector de longitud n_{j,q}
    n_j_q    = length(y_vec_jq)
    X        = x_ijq_test[[q]][[j]]           # matriz n_{j,q} × p_e
    Z        = z_jq_test[[q]][[j]]            # matriz n_{j,q} × p_m
    W        = w_q_test[[q]][[j]]             # matriz n_{j,q} × p_d
    
    # predicción puntual para cada i
    
    zeta_i <- b0*rep(1, n_j_q) + X %*% bE + Z %*% bM + W %*% bD
    
    zeta_list[[q]][[j]] <- zeta_i
  }
}




# 2) Inicializar vectores para observados y predichos
y_obs <- c()
y_hat <- c()

n_samples = 100

set.seed(777)
for (q in 1:m) {
  cat(sprintf("Generando muestras q = %d de %d\n", q, m))
  flush.console()
  
  Jk <- length(Y_test[[q]])
  
  for (j in seq_len(Jk)) {
    
    if (j > length(Estimadores_Bayesianos$kappa2_jq[[q]]) || 
        is.null(Estimadores_Bayesianos$kappa2_jq[[q]][[j]])) {
      cat(sprintf("Índice j = %d no encontrado en q = %d, pasando al siguiente q\n", j, q))
      break  # Sale del loop j y pasa al siguiente q
    }
    
    # Get stored data
    y_vec_jq = Y_test[[q]][[j]]               # vector de longitud n_{j,q}
    n_j_q    = length(y_vec_jq)
    zeta_i   = zeta_list[[q]][[j]]            # Get the stored zeta_i
    
    # Generate multivariate normal sample
    #kappa2_i_temp = Estimadores_Bayesianos$kappa2_jq[[q]][[j]] * diag(n_j_q)
    #y_hat_jq <- c(mvtnorm::rmvnorm(1, zeta_i, kappa2_i_temp))
    
    kappa2_escalar = Estimadores_Bayesianos$kappa2_jq[[q]][[j]]  
    
    #Promediar muestras
    
    y_hat_samples <- matrix(0, nrow = n_j_q , ncol = n_samples)
    
    for (s in 1:n_samples) {
      y_hat_samples[, s] <- rnorm(n_j_q, mean = zeta_i, sd = sqrt(kappa2_escalar))
    }
    
    # Promediar para cada individuo
    y_hat_jq <- rowMeans(y_hat_samples)
    
    # acumular en los vectores globales
    y_obs <- c(y_obs, y_vec_jq)
    y_hat <- c(y_hat, y_hat_jq)
  }
}


# 4) Calcular metricas
n = length(y_obs)
mse = mean((y_obs - y_hat)^2)
mae = mean(abs(y_obs - y_hat))
r2  =  1 - sum((y_obs - y_hat)^2) / sum((y_obs - mean(y_obs))^2)

# Almacenamos las metricas
metricas = data.frame(
  MSE = mse,
  MAE = mae,
  R2  = r2
)






### Estimaciones con zeta



# 2) Inicializar vectores para observados y predichos
y_obs <- c()
y_hat <- c()


set.seed(777)
for (q in 1:m) {
  cat(sprintf("Generando muestras q = %d de %d\n", q, m))
  flush.console()
  
  Jk <- length(Y_test[[q]])
  
  for (j in seq_len(Jk)) {
    
    if (j > length(Estimadores_Bayesianos$kappa2_jq[[q]]) || 
        is.null(Estimadores_Bayesianos$kappa2_jq[[q]][[j]])) {
      cat(sprintf("Índice j = %d no encontrado en q = %d, pasando al siguiente q\n", j, q))
      break  # Sale del loop j y pasa al siguiente q
    }
    
    # Get stored data
    y_vec_jq = Y_test[[q]][[j]]               # vector de longitud n_{j,q}
    n_j_q    = length(y_vec_jq)
    zeta_i   = zeta_list[[q]][[j]]            # Get the stored zeta_i
    
    # Generate multivariate normal sample
    #kappa2_i_temp = Estimadores_Bayesianos$kappa2_jq[[q]][[j]] * diag(n_j_q)
    #y_hat_jq <- c(mvtnorm::rmvnorm(1, zeta_i, kappa2_i_temp))
    y_hat_jq <- zeta_i
    
    
    # acumular en los vectores globales
    y_obs <- c(y_obs, y_vec_jq)
    y_hat <- c(y_hat, y_hat_jq)
  }
}


# 4) Calcular metricas
n = length(y_obs)
mse = mean((y_obs - y_hat)^2)
mae = mean(abs(y_obs - y_hat))
r2  =  1 - sum((y_obs - y_hat)^2) / sum((y_obs - mean(y_obs))^2)

# Almacenamos las metricas
metricas = data.frame(
  MSE = mse,
  MAE = mae,
  R2  = r2
)

summary(y_obs)
summary(y_hat)

#===========================================================================
#@========         Tamaños Efectivos de Muestra         =========
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

ess_summary$beta_D

#




calcular_mcse <- function(cadena) {
  if (length(cadena) == 0 || all(is.na(cadena))) {
    return(NA)
  }
  tryCatch({
    mcse_mean(as_draws_matrix(matrix(cadena, ncol = 1)))
  }, error = function(e) {
    warning(paste("Error calculando MCSE:", e$message))
    return(NA)
  })
}


mcse_beta <-  calcular_mcse(Resultados$beta)

beta_E_matrix <- do.call(cbind, Resultados$beta_E)
mcse_beta_E <- apply(beta_E_matrix, 2, calcular_mcse)


beta_M_matrix <- do.call(cbind, Resultados$beta_M)
mcse_beta_M <- apply(beta_M_matrix, 2, calcular_mcse)

beta_D_matrix <- do.call(cbind, Resultados$beta_D)
mcse_beta_D <- apply(beta_D_matrix, 2, calcular_mcse)


kappa2_q_matrix <- do.call(cbind, Resultados$kappa2_q)
mcse_kappa2_q <- apply(kappa2_q_matrix, 2, calcular_mcse)


kappa2_jq_list <- unlist(Resultados$kappa2_jq, recursive = FALSE)
mcse_kappa2_jq <- sapply(kappa2_jq_list, calcular_mcse)




resumen_5_numeros <- function(x) {
  x <- na.omit(as.numeric(x))
  c(
    min     = min(x),
    q1      = quantile(x, 0.25),
    median  = median(x),
    mean    = mean(x),
    q3      = quantile(x, 0.75),
    max     = max(x)
  )
}



mcse_summary <- list(
  beta_E = resumen_5_numeros(mcse_beta_E),
  beta_M = resumen_5_numeros(mcse_beta_M),
  beta_D = resumen_5_numeros(mcse_beta_D),
  kappa2_q = resumen_5_numeros(mcse_kappa2_q),
  kappa2_jq = resumen_5_numeros(mcse_kappa2_jq)
)






mcse_sigma2  <- calcular_mcse(Resultados$sigma2)
mcse_alpha   <- calcular_mcse(Resultados$alpha_kappa)
mcse_beta    <- calcular_mcse(Resultados$beta_kappa)

## MCSE para matrices: aplicar por columna
# mcse_theta <- apply(as_draws_matrix(Resultados$theta), 2, calcular_mcse)
# mcse_sigma2_q <- apply(as_draws_matrix(Resultados$sigma2_q), 2, calcular_mcse)
mcse_kappa2_q <- apply(as_draws_matrix(Resultados$kappa2_q), 2, calcular_mcse)

## MCSE para estructura jerárquica: zeta (departamentos -> municipios)
mcse_zeta <- unlist(lapply(1:length(Resultados$zeta), function(dept) {
  lapply(1:length(Resultados$zeta[[dept]]), function(mun) {
    calcular_mcse(Resultados$zeta[[dept]][[mun]])
  })
}))

## MCSE para estructura jerárquica: kappa2_jq
mcse_kappa2 <- unlist(lapply(1:length(Resultados$kappa2_jq), function(dept) {
  lapply(1:length(Resultados$kappa2_jq[[dept]]), function(mun) {
    calcular_mcse(Resultados$kappa2_jq[[dept]][[mun]])
  })
}))



summary(mcse_kappa2)










