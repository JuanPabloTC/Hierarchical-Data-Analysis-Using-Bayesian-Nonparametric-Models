#@ Datos
## Datos ICFES

# ===================================================================
#                             PREPARACION
# ===================================================================

library(readr)
library(tidyverse)
library(readxl)
library(stringi)
library(VIM)  # for KNN imputation
library(mice) # alternative imputation methods if needed

path <- c("C:/Users/jtorres/Desktop/From Documents/JP/Tesis/23-05-26")
carpeta_datos_originales <- paste0(path,"/Datos_originales/")
setwd(path)
path_to_save <- paste0(path,'/Datos_Procesados')
# Crea la carpeta si no existe
if (!dir.exists(path_to_save)) dir.create(path_to_save, recursive = TRUE, showWarnings = FALSE)

# ===================================================================
#                       CREACION BASE DE DATOS
# ===================================================================

# ========================================
# 1. Datos Prueba Saber 11
# ========================================

#Examen_Saber_11_20222_Orig <- read_delim(paste0(carpeta_datos_originales,"Examen_Saber_11_20222.txt"), 
#                                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

Examen_Saber_11_20222 <- read_delim(paste0(carpeta_datos_originales,"Examen_Saber_11_20222.txt"), 
                                    delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  filter(
    estu_nacionalidad == "COLOMBIA",
    estu_pais_reside == "COLOMBIA",
    cole_cod_depto_ubicacion != 88,
    estu_cod_reside_depto != 99999
  ) %>%
  select(
    periodo, estu_consecutivo, estu_cod_reside_depto, estu_depto_reside, 
    estu_cod_reside_mcpio, estu_mcpio_reside,
    fami_educacionmadre, fami_tienecomputador, fami_tieneinternet, 
    fami_numlibros, fami_estratovivienda, estu_tieneetnia, 
    estu_genero, cole_naturaleza, cole_calendario, estu_horassemanatrabaja, 
    cole_cod_depto_ubicacion, cole_cod_mcpio_ubicacion, 
    cole_mcpio_ubicacion, punt_global, cole_area_ubicacion
  ) %>%
  mutate(
    # Just factor the categorical variables - NO DUMMIES YET
    fami_numlibros = factor(fami_numlibros,
                            levels = c("0 A 10 LIBROS",
                                       "11 A 25 LIBROS",
                                       "26 A 100 LIBROS",
                                       "MÁS DE 100 LIBROS")
    ),
    fami_estratovivienda = factor(
      fami_estratovivienda,
      levels = c(
        "Sin Estrato",
        "Estrato 1",   
        "Estrato 2",
        "Estrato 3",
        "Estrato 4",
        "Estrato 5",
        "Estrato 6"
      )
    ),
    cole_calendario = factor(
      cole_calendario,
      levels = c(
        "OTRO",   
        "A",
        "B"
      )
    ),
    estu_horassemanatrabaja = factor(
      estu_horassemanatrabaja,
      levels = c(
        "0",
        "Menos de 10 horas",
        "Entre 11 y 20 horas", 
        "Entre 21 y 30 horas",
        "Más de 30 horas"
      )
    )
  )

#Correccion de municipio
Examen_Saber_11_20222$cole_cod_mcpio_ubicacion[Examen_Saber_11_20222$cole_mcpio_ubicacion == "BELÉN DE BAJIRÁ"] <- "27493"
Examen_Saber_11_20222$cole_cod_mcpio_ubicacion[Examen_Saber_11_20222$cole_cod_mcpio_ubicacion == 27086] <- 27493

# ========================================
# 2. Datos CEDE
# ========================================

PANEL_DE_EDUCACION_2022_ <- read_excel(paste0(carpeta_datos_originales,"PANEL_DE_EDUCACION(2022).xlsx"))
CEDE_Educacion <- PANEL_DE_EDUCACION_2022_ %>%
                  filter(ano == 2021) %>%
                  select(codmpio, ano, docen_total, alumn_total, alumn_oficial) %>%
                  mutate(docenttotal_alumtotal  = alumn_total/docen_total) %>%      # Alumnos por Docente
                  mutate(porc_alumn_col_publico = (alumn_oficial/alumn_total)*100)  # Porcentaje alumnos colegio publico

CEDE_Educacion$codmpio[CEDE_Educacion$codmpio == 27086] <- 27493

PANEL_DE_VIOLENCIA_2022_ <- read_excel(paste0(carpeta_datos_originales,"PANEL_CONFLICTO_Y_VIOLENCIA(2022).xlsx"))
CEDE_Violencia <- PANEL_DE_VIOLENCIA_2022_ %>%
                  filter(ano == 2019) %>%
                  select(codmpio, homicidios, terrorismot, ataq_instpol, hurto, secuestros)

CEDE_Violencia$codmpio[CEDE_Violencia$codmpio == 27086] <- 27493

PANEL_DE_GENERAL_2022_ <- read_excel(paste0(carpeta_datos_originales,"PANEL_CARACTERISTICAS_GENERALES(2022).xlsx"))
CEDE_General <- PANEL_DE_GENERAL_2022_ %>%
              filter(ano == 2021) %>%
              select(coddepto, codprovincia, codmpio, ano, pobl_rur, pobl_tot, discapital) #%>%
              #mutate(proporcio_pob_rural = (100*pobl_rur)/pobl_tot)

CEDE_General <- CEDE_General %>%
  group_by(coddepto) %>%
  mutate(
    total_pobl_rur_depto = sum(pobl_rur, na.rm = TRUE),
    total_pobl_depto     = sum(pobl_tot, na.rm = TRUE),
    proporcio_pob_rural =  100 * total_pobl_rur_depto / total_pobl_depto
  ) %>%
  ungroup()

CEDE_General$codmpio[CEDE_General$codmpio == 27086] <- 27493

CEDE_Data <- CEDE_General %>%
              left_join(CEDE_Violencia, by = c("codmpio")) %>%
              left_join(CEDE_Educacion, by = c("codmpio", "ano")) %>%
              mutate("Homi_x_100k_habitantes" = (homicidios/pobl_tot)*100000)

# ========================================
# 3. Datos DIVIPOLA
# ========================================

DANE_DIVIPOLA <- read_excel(paste0(carpeta_datos_originales,"DIVIPOLA_Municipios.xlsx")) %>%
                 drop_na(`...2`) %>%
                  {
                    header <- as.character(.[1, ])
                    .[-1, ] %>% 
                      setNames(header)
                  }

DANE_DIVIPOLA <- DANE_DIVIPOLA[, !is.na(names(DANE_DIVIPOLA)) & names(DANE_DIVIPOLA) != ""]
names(DANE_DIVIPOLA)[which(names(DANE_DIVIPOLA) == "Nombre")[1]] <- "nombre_departamento"
names(DANE_DIVIPOLA)[4] <- "nombre_municipio"
names(DANE_DIVIPOLA)[which(names(DANE_DIVIPOLA) == "Código")[1]] <- "codigo_departamento"
names(DANE_DIVIPOLA)[3] <- "codigo_municipio"

DANE_DIVIPOLA <- DANE_DIVIPOLA %>%
  mutate(
    departamento_sin_tilde = stri_trans_general(nombre_departamento, "Latin-ASCII"),
    municipio_sin_tilde = stri_trans_general(nombre_municipio, "Latin-ASCII"),
    codigo_departamento = sub("^0+", "", codigo_departamento)  # Remove leading zeros
  ) %>%
  mutate(
    nombre_departamento = trimws(nombre_departamento),
    nombre_municipio = trimws(nombre_municipio)
  )


#== Create new column named cole_depto_ubicacion

# Create lookup table
dept_lookup <- DANE_DIVIPOLA %>%
  distinct(codigo_departamento, .keep_all = TRUE) %>%
  select(codigo_departamento, nombre_departamento)

# Join with type conversion 
Examen_Saber_11_20222 <- Examen_Saber_11_20222 %>%
  mutate(cole_cod_depto_ubicacion = as.character(cole_cod_depto_ubicacion)) %>%  # Convert numeric to character
  left_join(dept_lookup, 
            by = c("cole_cod_depto_ubicacion" = "codigo_departamento")) %>%
  rename(cole_depto_ubicacion = nombre_departamento)


# ========================================
# 3. Datos LAFT
# ========================================


Data_LAFT_Orig <- read_excel(paste0(carpeta_datos_originales,"LAFT.xlsx"))
                
Data_LAFT <- Data_LAFT_Orig %>%
             select(CODIGO, RISK_VICTIM_2022) %>%
             separate(CODIGO, into = c("departamento", "municipio", "resto"), sep = "-", fill = "right") %>%
             mutate(
              across(c(departamento, municipio), trimws),
              departamento_sin_tilde = stri_trans_general(departamento, "Latin-ASCII"),
              municipio_sin_tilde    = stri_trans_general(municipio, "Latin-ASCII")
             )

Data_LAFT <- Data_LAFT %>%
             left_join(DANE_DIVIPOLA, by = c("departamento_sin_tilde", "municipio_sin_tilde")) %>%
             select(departamento, municipio, codigo_departamento, codigo_municipio, RISK_VICTIM_2022) %>%
             mutate(codigo_departamento = sub("^0", "", as.character(codigo_departamento)),
                    codigo_municipio = sub("^0", "", as.character(codigo_municipio)))

Data_LAFT$codigo_departamento[Data_LAFT$municipio == "PIENDAMO"] <- "19"
Data_LAFT$codigo_departamento[Data_LAFT$municipio == "SOTARA PAISPAMBA"] <- "19"
Data_LAFT$codigo_departamento[Data_LAFT$municipio == "MIRITI"] <- "91"
Data_LAFT$codigo_municipio[Data_LAFT$municipio == "SOTARA PAISPAMBA"] <- "19760"
Data_LAFT$codigo_municipio[Data_LAFT$municipio == "PIENDAMO"] <- "19548"
Data_LAFT$codigo_municipio[Data_LAFT$municipio == "MIRITI"] <- "91460"

# ========================================
# 5. Datos DANE (PIB)
# ========================================

DANE_PIB_percapita_Orig <- read_excel(paste0(carpeta_datos_originales,"DANE - PIB.xlsx"), sheet = 4)
DANE_PIB_percapita <- read_excel(paste0(carpeta_datos_originales,"DANE - PIB.xlsx"), sheet = 4) %>%
                            drop_na(`...2`) %>%
                            {
                              header <- as.character(.[1, ])
                              .[-1, ] %>% 
                                setNames(header)
                            } %>% 
                            select(`Código Departamento (DIVIPOLA)`, `DEPARTAMENTOS`, `2022`) %>% 
                            rename('PIB_percapita_DPTO' = '2022') %>%   
                            mutate(`PIB_percapita_DPTO` = `PIB_percapita_DPTO`/1000000,   
                                   `Código Departamento (DIVIPOLA)` = sub("^0", "", as.character(`Código Departamento (DIVIPOLA)`)))                          

# ========================================
# 5. Datos MOE
# ========================================

MOE_mun_riesgo_orig <- read_excel(paste0(carpeta_datos_originales,"MOE.xlsx")) 
MOE_mun_riesgo <- MOE_mun_riesgo_orig %>% 
                  select(Depto,
                         perc_municipios_con_riesgo = `% municipios con riesgo`) %>%
                  mutate(departamento_sin_tilde = trimws(stri_trans_general(toupper(Depto), "Latin-ASCII")))

DANE_departamentos <- DANE_DIVIPOLA %>%
                      distinct(departamento_sin_tilde, .keep_all = TRUE) %>%
                      select(codigo_departamento, departamento_sin_tilde)

MOE_mun_riesgo <- MOE_mun_riesgo %>%
                  left_join(DANE_departamentos, by = "departamento_sin_tilde") %>%
                  select(departamento_sin_tilde, codigo_departamento, `perc_municipios_con_riesgo`) %>%
                  mutate(codigo_departamento = sub("^0", "", as.character(codigo_departamento)))

MOE_mun_riesgo$codigo_departamento[MOE_mun_riesgo$departamento_sin_tilde == "ARCHIPIELAGO DE SAN ANDRES"] <- "88"
MOE_mun_riesgo$codigo_departamento[MOE_mun_riesgo$departamento_sin_tilde == "BOGOTA D.C."] <- "11"


# ========================================
# 7. Datos SICODIS/Depto Nal de Planeacion
# ========================================

regalias_mun_Orig <- read_excel(paste0(carpeta_datos_originales,"ResumenDsitribucionSGP_2022.xlsx")) 

regalias_mun <- regalias_mun_Orig %>%
                select("Codigo DANE Entidad", "Nombre Entidad", "Tipo Entidad", "Educación", "Salud") %>%
                rename('codmpio' = 'Codigo DANE Entidad',
                       'Educacion' = 'Educación') %>% 
                mutate(codmpio = as.character(codmpio), 
                       codmpio = str_replace(codmpio, "^0+", "")) %>% #Remover 0 cuando inician con 0
                left_join(CEDE_General %>%                            #Juntar con poblacion total
                          mutate(codmpio = as.character(codmpio)) %>%
                          select(codmpio, pobl_tot), 
                          by = "codmpio") %>%
                mutate("regalias_educa_per_capita" = Educacion / pobl_tot) %>% # Regalias por habitante
                select(- pobl_tot)
# ========================================
# 7. Unificar datos
# ========================================

Examen_Saber_11_20222$cole_cod_mcpio_ubicacion <- as.character(Examen_Saber_11_20222$cole_cod_mcpio_ubicacion)
Examen_Saber_11_20222$cole_cod_depto_ubicacion <- as.character(Examen_Saber_11_20222$cole_cod_depto_ubicacion)
CEDE_Data$codmpio <- as.character(CEDE_Data$codmpio)

Datos_agrupados <- Examen_Saber_11_20222 %>%
                   left_join(CEDE_Data, by = c("cole_cod_mcpio_ubicacion" = "codmpio")) %>%
                   left_join(Data_LAFT, by = c("cole_cod_mcpio_ubicacion" = "codigo_municipio")) %>%
                   left_join(DANE_PIB_percapita, by = c("cole_cod_depto_ubicacion" = "Código Departamento (DIVIPOLA)")) %>%
                   left_join(MOE_mun_riesgo, by = c("cole_cod_depto_ubicacion" = "codigo_departamento")) %>%
                   left_join(regalias_mun, by = c("cole_cod_mcpio_ubicacion" = "codmpio"))

# Creacion variable homicidios ponderados por depto
Homicidios_ponderado <- Datos_agrupados %>%
  select(cole_cod_depto_ubicacion, Homi_x_100k_habitantes, pobl_tot) %>%
  group_by(cole_cod_depto_ubicacion) %>%
  summarise(
    Homicidios_ponderado_x_100k = sum(Homi_x_100k_habitantes * pobl_tot, na.rm = TRUE) / 
      sum(pobl_tot, na.rm = TRUE)
  )
Datos_agrupados <- Datos_agrupados %>%
  left_join(Homicidios_ponderado, by = "cole_cod_depto_ubicacion")

Datos_agrupados <- Datos_agrupados %>%
                   select(
                     # ID and geographic variables
                     "periodo", "estu_consecutivo", 
                     "estu_cod_reside_depto", "estu_depto_reside", 
                     "estu_cod_reside_mcpio", "estu_mcpio_reside",
                     "cole_cod_depto_ubicacion", "cole_cod_mcpio_ubicacion", "cole_mcpio_ubicacion",
                     "cole_depto_ubicacion",
                     
                     # Outcome variable
                     "punt_global",
                     
                     # Individual-level categorical variables (before dummies)
                     "fami_educacionmadre", "fami_tienecomputador", "fami_tieneinternet", 
                     "fami_numlibros", "fami_estratovivienda", "estu_tieneetnia",
                     "estu_genero", "cole_naturaleza", "cole_calendario", 
                     "estu_horassemanatrabaja", "cole_area_ubicacion",
                     
                     # CEDE variables (municipality/department level)
                     "coddepto", "codprovincia", "ano", 
                     "pobl_rur", "pobl_tot", "proporcio_pob_rural",
                     "homicidios", "Homi_x_100k_habitantes",
                     "docen_total", "alumn_total", "docenttotal_alumtotal",
                     "porc_alumn_col_publico", "terrorismot", "ataq_instpol",
                     "hurto", "secuestros", "discapital", "regalias_educa_per_capita",
                     
                     # LAFT variable
                     "RISK_VICTIM_2022",
                     
                     # PIB variable
                     "PIB_percapita_DPTO",
                     
                     # MOE variable
                     "perc_municipios_con_riesgo",
                     
                     # Homicidios ponderado
                     "Homicidios_ponderado_x_100k"
                   )

# --- Nos aseguramos que no haya valores faltantes para proporcio_pob_rural

Datos_agrupados <- Datos_agrupados %>%
  group_by(cole_cod_depto_ubicacion) %>%
  mutate(
    proporcio_pob_rural = if_else(
      is.na(proporcio_pob_rural),
      first(proporcio_pob_rural[!is.na(proporcio_pob_rural)]),
      proporcio_pob_rural
    )
  ) %>%
  ungroup()

# --- Columnas a usar 

#x_ijk_names <- c("fami_educacionmadre_modif", "computador", "internet", "etnia",
#                 "libros_11_25", "libros_26_100", "libros_mas100", 
#                 "estrato_1", "estrato_2", "estrato_3", "estrato_4", 
#                 "estrato_5", "estrato_6", "Genero_mujer", "Calendario_A", 
#                 "Calendario_B", "cole_privado", "trabaja_menos_de_10_horas",  
#                 "trabaja__11_a_20_horas", "trabaja__21_a_30_horas", 
#                 "trabaja_mas_de_30_horas", "cole_urbano")

#z_jk_names <-  c("docenttotal_alumtotal", "RISK_VICTIM_2022", "Homi_x_100k_habitantes",
#                 "porc_alumn_col_publico", "terrorismot", "hurto", "secuestros", "discapital",
#                 "regalias_educa_per_capita")

#w_k_names <-  c("PIB_percapita_DPTO", "proporcio_pob_rural", "perc_municipios_con_riesgo", 
#                "Homicidios_ponderado_x_100k")

# Eliminar municipio Mapiripana, tiene muchos datos faltantes y es problemático.

Datos_agrupados <- Datos_agrupados %>%
  filter(estu_cod_reside_mcpio != "94663")


#load("Datos_Procesados/Datos_pre_imputacion.RData")

# ===================================================================
#                     DEJAR DATOS REGIÓN CENTRAL
# ===================================================================

# Se va a dejar los datos de la región central y Bogotá acorde con el DANE, la cual está 
#compuesta por los siguientes departamentos y su respectivo codigo Divipola:
# - Caldas       (17)
# - Risaralda    (66)
# - Quindío      (63)
# - Tolima       (73)
# - Huila        (41)
# - Caquetá      (18)
# - Antioquia    (05)
# - Bogotá       (11)

# # Filter the column "cole_cod_depto_ubicacion"
# 
# codes_region_central <- c("17","66","63","73","41","18","5","11")
# 
# Datos_filtrado <- Datos_agrupados %>%
#   filter(cole_cod_depto_ubicacion %in% codes_region_central)
# 
# 
# cat("The central region sample have ", nrow(Datos_filtrado), "\n")
# 
# cat("which means that ", nrow(Datos_agrupados) - nrow(Datos_filtrado), 
#     "or ", 100 *(nrow(Datos_agrupados) - nrow(Datos_filtrado))/ nrow(Datos_agrupados), "%", 
#     " were eliminated", "\n")
# 
# 
# count_11_Bog <- Datos_agrupados %>%
#   filter(cole_cod_depto_ubicacion == "11") %>%
#   nrow()
# 
# cat("There are", count_11_Bog, "samples only in Bogotá")
# 
# # overwrite the database
# Datos_agrupados <- Datos_filtrado
# 
# #eliminate the temporal database
# rm(Datos_filtrado)

# ===================================================================
#                     IMPUTACION DATOS FALTANTES
# ===================================================================

# ========================================
# 1. Eliminar obs con muchos faltantes
# ========================================

# Se eliminan individuos con más de cierto porcentaje de datos faltantes dependiendo del tamaño de su municipio
# y se eliminan municipios con menos de 10 individuos

#------------------  DIAGNOSTICOS -----------------------

# --- Conteo variables

# Individuales
individual_vars_to_impute <- c(
  "fami_educacionmadre", "fami_tienecomputador", "fami_tieneinternet", 
  "estu_tieneetnia", "fami_numlibros", "fami_estratovivienda",
  "estu_genero", "cole_calendario", "cole_naturaleza", 
  "estu_horassemanatrabaja", "cole_area_ubicacion"
)
#Municipales
municipality_vars_to_impute <- c(
  "docenttotal_alumtotal", "RISK_VICTIM_2022", "Homi_x_100k_habitantes",
  "porc_alumn_col_publico", "terrorismot", "hurto", "secuestros", "discapital",
  "regalias_educa_per_capita"
)
#Departamentales
department_vars_to_impute <- c(
  "PIB_percapita_DPTO", "proporcio_pob_rural", "perc_municipios_con_riesgo", 
  "Homicidios_ponderado_x_100k"
)

# Numero total de variables
all_vars_to_impute <- c(individual_vars_to_impute, 
                        municipality_vars_to_impute, 
                        department_vars_to_impute)

total_vars <- length(all_vars_to_impute)


missing_summary <- Datos_agrupados %>%
       summarise(across(everything(), ~sum(is.na(.)))) %>%
       pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
       mutate(
             total_obs = nrow(Datos_agrupados),
             pct_missing = (n_missing / total_obs) * 100
         ) %>%
       arrange(desc(pct_missing))

# No hay variables con más de 10% de datos faltantes.
cat("Variable with most missing values:", missing_summary[[1, 1]], 
    "with", missing_summary[[1, "pct_missing"]], "% missing values")

# --- numero individuos, municipios y departamento

cat("\n=== INITIAL DATA ===\n")
cat("Total students:", nrow(Datos_agrupados), "\n")
cat("Total municipalities:", n_distinct(Datos_agrupados$cole_cod_mcpio_ubicacion), "\n")
cat("Total departments:", n_distinct(Datos_agrupados$cole_cod_depto_ubicacion), "\n")

# --- Faltantes a nivel departamental 

# A nivel departamental no hay que imputar 
anyNA(Datos_agrupados$PIB_percapita_DPTO)
anyNA(Datos_agrupados$proporcio_pob_rural)
anyNA(Datos_agrupados$perc_municipios_con_riesgo)
anyNA(Datos_agrupados$Homicidios_ponderado_x_100k)



#------------------  CALCULO PORCENTAJE DATOS FALTANTES -----------------------

# To observe the sizes of the municipalities
resumen_concentracion <- Datos_agrupados %>%
  count(cole_cod_mcpio_ubicacion) %>%
  arrange(desc(n)) %>%
  mutate(
    pct = 100 * n / sum(n),
    cum_pct = cumsum(pct)
  )

#view(resumen_concentracion)


# Computar tamaños de muestras en municipio y asignarles un tamaño (pequeño, mediano, grande)

# municipal sizes are:
# small if < 64
# medium if > 65 but < 299
# large if > 300
mun_size <- Datos_agrupados %>%
  count(cole_cod_mcpio_ubicacion, name = "n_mun") %>%
  mutate(
    size_class = case_when(
      n_mun < 65  ~ "small",
      n_mun < 300 ~ "medium",
      TRUE        ~ "large"
    )
  )

Datos_agrupados <- Datos_agrupados %>%
  left_join(mun_size, by = "cole_cod_mcpio_ubicacion")

cat("\n=== MUNICIPALITY SIZE DISTRIBUTION ===\n")
print(table(mun_size$size_class))



# --- Calcular datos faltantes por estudiantes

# The threshold to eliminate observations for missing variables % according to their mun_size:
# small  < 35%
# medium < 25%
# large  < 10%

Datos_agrupados <- Datos_agrupados %>%
  mutate(
    n_missing_vars = rowSums(is.na(select(., all_of(all_vars_to_impute)))),
    pct_missing = (n_missing_vars / total_vars) * 100,
    # Define threshold based on municipality size
    threshold = case_when(  #Se crea umbral 
      size_class == "large"  ~ 5,
      size_class == "medium" ~ 10,
      size_class == "small"  ~ 20
    )
  )

# Note que en general todas las observaciones tienen datos
cat("\n=== MISSING DATA PER STUDENT (before filtering) ===\n")
print(summary(Datos_agrupados$pct_missing))

cat("\nStudents by missingness category:\n")
print(table(cut(Datos_agrupados$pct_missing, 
                breaks = c(-Inf, 0, 10, 20, 30, 40, 50, Inf),
                labels = c("0%", "1-10%", "11-20%", "21-30%", "31-40%", "41-50%", ">50%"))))


#------------------  FILTRAR INDIVIDUOS -----------------------

# Se filtran individuos de acuerdo al umbral
#      size_class == "large"  ~ 25%
#      size_class == "medium" ~ 35%
#      size_class == "small"  ~ 45%

Datos_agrupados_filtered <- Datos_agrupados %>%
  filter(pct_missing <= threshold)

# --- Diagnostico de cuantos se filtraron
cat("\n=== STUDENTS REMOVED ===\n")
cat("Students before filtering:", nrow(Datos_agrupados), "\n")
cat("Students after filtering:", nrow(Datos_agrupados_filtered), "\n")
cat("Students removed:", nrow(Datos_agrupados) - nrow(Datos_agrupados_filtered), "\n")
cat("% removed:", round((1 - nrow(Datos_agrupados_filtered)/nrow(Datos_agrupados)) * 100, 2), "%\n")

# --- Recalcular tamaños de municipios despues del filtrado
mun_size_after <- Datos_agrupados_filtered %>%
  count(cole_cod_mcpio_ubicacion, name = "n_mun_after")

print(table(mun_size$size_class))

#------------------  FILTRAR MUNICIPIOS -----------------------

# Se filtran municipios con menos de 50 observaciones

municipalities_to_keep <- mun_size_after %>%
  filter(n_mun_after >= 50) %>%
  pull(cole_cod_mcpio_ubicacion)

Datos_agrupados_filtered <- Datos_agrupados_filtered %>%
  filter(cole_cod_mcpio_ubicacion %in% municipalities_to_keep)

cat("\n=== MUNICIPALITIES REMOVED ===\n")
cat("Municipalities before:", n_distinct(Datos_agrupados$cole_cod_mcpio_ubicacion), "\n")
cat("Municipalities after:", n_distinct(Datos_agrupados_filtered$cole_cod_mcpio_ubicacion), "\n")
cat("Municipalities removed:", 
    n_distinct(Datos_agrupados$cole_cod_mcpio_ubicacion) - 
      n_distinct(Datos_agrupados_filtered$cole_cod_mcpio_ubicacion), "\n")

#------------------  DIAGNOSTICO DATOS DESPUES DEL FILTRO -----------------------

cat("\n=== REMAINING MISSING DATA (in variables to impute) ===\n")
missing_summary_after <- Datos_agrupados_filtered %>%
  select(all_of(all_vars_to_impute)) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(
    total_obs = nrow(Datos_agrupados_filtered),
    pct_missing = (n_missing / total_obs) * 100
  ) %>%
  filter(n_missing > 0) %>%
  arrange(desc(pct_missing))

print(missing_summary_after)

cat("Final number of observations:", nrow(Datos_agrupados_filtered))

# Update the main dataset
Datos_agrupados <- Datos_agrupados_filtered
# Eliminamos el dataframe por memoria
rm(Datos_agrupados_filtered)

# ========================================
# 2. Crear funcion para Dummys
# ========================================

# Se va a crear una funcion para crear las variables dummy a partir de los datos categoricos.
# Esta función se va a usar después de la imputacion

## Guia de dummies:
# - fami_educacionmadre_modif:  {0: no tiene educacion superior} {1: tiene educacion superior}
# - computador:                 {0: no tiene computador}         {1: tiene educacion computador}
# - internet:                   {0: no tiene internet}           {1: tiene educacion internet}
# - Genero_mujer:               {0: Estudiante hombre}           {1: Estudiante mujer}
# - cole_privado:               {0: Colegio publico}             {1: Colegio privado}
# - cole_urbano:                {0: Colegio Rural}               {1: Colegio urbano}


# - fami_numlibros:             CATEGORIA BASE "0 A 10 LIBROS"
#                               Columnas: "libros_11_25","libros_26_100", "libros_mas100"
# - fami_estratovivienda:       CATEGORIA BASE "Sin Estrato"
#                               Columnas: "estrato_1","estrato_2", "estrato_3", "estrato_4","estrato_5", "estrato_6"
# - estu_horassemanatrabaja:    CATEGORIA BASE: "0"
#                               Columnas: "trabaja_menos_de_10_horas", "trabaja__11_a_20_horas", "trabaja__21_a_30_horas", "trabaja_mas_de_30_horas"
# - cole_calendario:            CATEGORIA BASE: "OTRO"
#                               Columnas: "Calendario_A", Calendario_B"

create_dummy_variables <- function(data) {
  data %>%
    mutate(
      # Educacion Madre (from fami_educacionmadre)
      fami_educacionmadre_modif = ifelse(
        fami_educacionmadre %in% c("Educación profesional completa", "Postgrado"),
        1, 0
      ),
      
      # Libros (from fami_numlibros) - CATEGORIA BASE "0 A 10 LIBROS"
      libros_11_25  = as.integer(fami_numlibros == "11 A 25 LIBROS"),
      libros_26_100 = as.integer(fami_numlibros == "26 A 100 LIBROS"),
      libros_mas100 = as.integer(fami_numlibros == "MÁS DE 100 LIBROS"),
      
      # Estrato (from fami_estratovivienda) - CATEGORIA BASE "Sin Estrato"
      estrato_1 = as.integer(fami_estratovivienda == "Estrato 1"),
      estrato_2 = as.integer(fami_estratovivienda == "Estrato 2"),
      estrato_3 = as.integer(fami_estratovivienda == "Estrato 3"),
      estrato_4 = as.integer(fami_estratovivienda == "Estrato 4"),
      estrato_5 = as.integer(fami_estratovivienda == "Estrato 5"),
      estrato_6 = as.integer(fami_estratovivienda == "Estrato 6"),
      
      # Calendario (from cole_calendario) - CATEGORIA BASE "OTRO"
      Calendario_A = as.integer(cole_calendario == "A"),
      Calendario_B = as.integer(cole_calendario == "B"),
      
      # Horas trabajadas (from estu_horassemanatrabaja) - CATEGORIA BASE "0"
      trabaja_menos_de_10_horas = as.integer(estu_horassemanatrabaja == "Menos de 10 horas"),
      trabaja__11_a_20_horas    = as.integer(estu_horassemanatrabaja == "Entre 11 y 20 horas"),
      trabaja__21_a_30_horas    = as.integer(estu_horassemanatrabaja == "Entre 21 y 30 horas"),
      trabaja_mas_de_30_horas   = as.integer(estu_horassemanatrabaja == "Más de 30 horas"),
      
      # Binary variables (from various sources)
      computador   = as.integer(fami_tienecomputador == "Si"),
      internet     = as.integer(fami_tieneinternet   == "Si"),
      etnia        = as.integer(estu_tieneetnia      == "Si"),
      Genero_mujer = as.integer(estu_genero          == "F"),
      cole_privado = as.integer(cole_naturaleza == "NO OFICIAL"),
      cole_urbano  = as.integer(cole_area_ubicacion == "URBANO")
    )
}



# ========================================
# 3. Imputar por KNN
# ========================================

# --- Imputar variables individuales.

# Check missing data before imputation
cat("\n=== MISSING DATA BEFORE IMPUTATION ===\n")
missing_before <- Datos_agrupados %>%
  select(all_of(all_vars_to_impute)) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0) %>%
  arrange(desc(n_missing))
print(missing_before)


# Todas las variables individuales son categoricas


cat("\n=== IMPUTING INDIVIDUAL-LEVEL VARIABLES ===\n")

# Get list of departments
departments <- unique(Datos_agrupados$cole_cod_depto_ubicacion)
cat("Number of departments:", length(departments), "\n")

# Create a list to store imputed data by department
imputed_data_list <- list()

# Loop through each department and apply KNN imputation

set.seed(777)

for (dept in departments) {
  
  cat("\nProcessing department:", dept, "\n")
  
  # Subset data for this department
  dept_data <- Datos_agrupados %>%
    filter(cole_cod_depto_ubicacion == dept)
  
  cat("  Students in department:", nrow(dept_data), "\n")
  
  # Check if there's any missing data in individual variables for this department
  missing_in_dept <- dept_data %>%
    select(all_of(individual_vars_to_impute)) %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    sum()
  
  if (missing_in_dept > 0) {
    cat("  Missing values to impute:", missing_in_dept, "\n")
    
    # Select variables for imputation (individual level + identifiers)
    vars_for_knn <- c(individual_vars_to_impute, "punt_global")
    
    # Apply KNN imputation
    
    # k = min(5, floor(nrow(dept_data) * 0.1)) ensures we don't use too many neighbors in small departments
    # El numero de vecinos está entre 3 y 5:
      # nrow(dept_data) * 0.05): Se calcula el 5% de las filas
      # floor redondea al numero entero por debajo
    k_neighbors <- min(5, max(3, floor(nrow(dept_data) * 0.05))) 
    
    cat("  Using k =", k_neighbors, "neighbors\n")
    
    # Se imputan las variables por cada depto
    dept_data_imputed <- kNN(
      data = dept_data,
      variable = individual_vars_to_impute, #variables where missing values should be imputed
      k = k_neighbors, #number of Nearest Neighbours used
      dist_var = vars_for_knn, #names or variables to be used for distance calculation
      imp_var = FALSE  # Don't create indicator variables for imputation
    )
    
    imputed_data_list[[dept]] <- dept_data_imputed
    
  } else {
    cat("  No missing data - skipping imputation\n")
    imputed_data_list[[dept]] <- dept_data
  }
}

# Combine all departments back together
cat("\n=== COMBINING IMPUTED DATA ===\n")
Datos_agrupados_imputed <- bind_rows(imputed_data_list)
cat("Total students after individual-level imputation:", nrow(Datos_agrupados_imputed), "\n")


#save(Datos_agrupados_imputed, file = paste0(path_to_save, "/Datos_agrupados_imputados.RData"))

# --- Verificar que haya funcionado la imputacion

load("Datos_Procesados/Datos_agrupados_imputados.RData")
# Se verifica que haya funcionado 
Datos_agrupados_imputed %>%
  summarise(across(all_of(individual_vars_to_impute), ~ sum(is.na(.))))

# Todas las variables han sido imputadas, no se necesita el bloque de municipio.
Datos_agrupados_imputed %>%
  summarise(across(all_of(all_vars_to_impute), ~ sum(is.na(.))))

anyNA(Datos_agrupados_imputed)

cat("\n=== FINAL MISSING DATA CHECK ===\n")
missing_after <- Datos_agrupados_imputed %>%
  select(all_of(all_vars_to_impute)) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0)

if (nrow(missing_after) > 0) {
  cat("WARNING: Some missing data remains\n")
  print(missing_after)
} else {
  cat("Success: All missing data has been imputed\n")
}

# Update main dataset
Datos_agrupados <- Datos_agrupados_imputed
rm(Datos_agrupados_imputed)


# ========================================
# 4. Crear Dummys y finalizar DataFrame
# ========================================


cat("\n\n===============================================\n")
cat("CREATING DUMMY VARIABLES\n")
cat("===============================================\n")

# Apply the dummy variable creation function we defined in Block 2
Datos_agrupados <- create_dummy_variables(Datos_agrupados)

cat("Dummy variables created successfully!\n")

cat("\n=== SELECTING FINAL VARIABLES ===\n")

Datos_agrupados <- Datos_agrupados %>%
  select(
    # ID and geographic variables
    "periodo", "estu_consecutivo", 
    "estu_cod_reside_depto", "estu_depto_reside", 
    "estu_cod_reside_mcpio", "estu_mcpio_reside",
    "cole_cod_depto_ubicacion", "cole_cod_mcpio_ubicacion", "cole_mcpio_ubicacion",
    "cole_depto_ubicacion",
    
    # Outcome variable
    "punt_global",
    
    # Individual-level DUMMY variables (for modeling)
    "fami_educacionmadre_modif", "computador", "internet", "etnia",
    "libros_11_25", "libros_26_100", "libros_mas100",
    "estrato_1", "estrato_2", "estrato_3", "estrato_4", "estrato_5", "estrato_6",
    "Genero_mujer", "Calendario_A", "Calendario_B", 
    "cole_privado", "trabaja_menos_de_10_horas", "trabaja__11_a_20_horas",
    "trabaja__21_a_30_horas", "trabaja_mas_de_30_horas", "cole_urbano",
    
    # Municipality-level variables (for modeling)
    "docenttotal_alumtotal", "RISK_VICTIM_2022", "Homi_x_100k_habitantes",
    "porc_alumn_col_publico", "terrorismot", "hurto", "secuestros", "discapital",
    "regalias_educa_per_capita",
    
    # Department-level variables (for modeling)
    "PIB_percapita_DPTO", "proporcio_pob_rural", "perc_municipios_con_riesgo", 
    "Homicidios_ponderado_x_100k"
  )

cat("Final dataset dimensions:", nrow(Datos_agrupados), "rows x", ncol(Datos_agrupados), "columns\n")

#x_ijk_names <- c("fami_educacionmadre_modif", "computador", "internet", "etnia",
#                 "libros_11_25", "libros_26_100", "libros_mas100", 
#                 "estrato_1", "estrato_2", "estrato_3", "estrato_4", 
#                 "estrato_5", "estrato_6", "Genero_mujer", "Calendario_A", 
#                 "Calendario_B", "cole_privado", "trabaja_menos_de_10_horas",  
#                 "trabaja__11_a_20_horas", "trabaja__21_a_30_horas", 
#                 "trabaja_mas_de_30_horas", "cole_urbano")

#z_jk_names <-  c("docenttotal_alumtotal", "RISK_VICTIM_2022", "Homi_x_100k_habitantes",
#                 "porc_alumn_col_publico", "terrorismot", "hurto", "secuestros", "discapital",
#                 "regalias_educa_per_capita")

#w_k_names <-  c("PIB_percapita_DPTO", "proporcio_pob_rural", "perc_municipios_con_riesgo", 
#                "Homicidios_ponderado_x_100k")


# ===================================================================
#                     REMUESTREAR
# ===================================================================

# Se remuestrea la base de datos para dejar alrededor de 100.000 observaciones. Se considera 
# la diferencia de tamaño de los municipios para que todos tengan representación.

# Function: stratified_subsample_simple
# - Datos: full dataframe
# - mun_id: string name of municipality id column
# - outcome_var: variable to compute variance: "punt_global"
# - n_target: desired overall sample size (approx.)
# - m_min: minimum sample per municipality (integer)
# - m_max: maximum sample per municipality: Inf for no cap
# - seed: random seed for reproducibility
# - return_test_set: TRUE to return the test_set
#
# Returns a list with:
#  $train_df      -> sampled dataframe (with sample_weight and sample_prob)
#  $test_df       -> test set (if return_test_set = TRUE)
#  $alloc_table   -> dataframe with N_h, S_h, n_alloc, and final sample fraction
#  $meta          -> metadata (requested n_target, actual n_sampled, sum_n_alloc)

# IMPLEMENTATION:
#   1) Neyman-allocate a combined pool of size 2 * n_target, with a
#      per-municipality floor of 2 * m_min. This guarantees each
#      municipality enters the combined pool with >= 2 * m_min observations
#      (or all of N_h if N_h < 2 * m_min).
#   2) Stratified 50/50 split within each municipality:
#      - n_test = floor(n_h / 2) -> test
#      - n_h - n_test = ceiling(n_h / 2) -> train
#      So train always has ceiling((2 * m_min) / 2) = m_min observations
#      per municipality (or all of N_h if smaller).

subsample_stratified_neyman <- function(Datos,
                                        mun_id = "cole_cod_mcpio_ubicacion",
                                        outcome_var = "punt_global",
                                        n_target = 10000,
                                        m_min = 10,
                                        m_max = Inf,
                                        seed = 777,
                                        return_test_set = TRUE) {
  set.seed(seed)
  
  # Internal pool target and per-municipality floor
  pool_target <- 2L * n_target
  m_min_pool  <- 2L * m_min
  
  # Temporary Row ID for tracking
  Datos_with_id <- Datos %>%
    mutate(.row_id_temp = row_number())
  
  # 1) Compute N_h and S_h (standard deviation) for each municipality
  N_by_mun <- Datos %>%
    group_by(.data[[mun_id]]) %>%
    summarise(
      N_h = n(),
      S_h = sd(.data[[outcome_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(N_h))
  
  # Defensive: municipalities with N_h = 1 have S_h = NA; treat as 0 so
  # Neyman gives them their floor allocation rather than NA.
  N_by_mun <- N_by_mun %>%
    mutate(S_h = ifelse(is.na(S_h), 0, S_h))
  
  N_total <- sum(N_by_mun$N_h)
  H <- nrow(N_by_mun)
  cat("Total observations N = ", N_total, "; Municipalities H = ", H, "\n")
  cat("Training target = ", n_target, "; pool target = ", pool_target, "\n")
  cat("Min per municipality (train) = ", m_min,
      "; min per municipality (pool) = ", m_min_pool, "\n")
  
  # 2) Neyman allocation toward the POOL: n_h = pool_target * (N_h * S_h) / sum(N_i * S_i)
  sum_N_S <- sum(N_by_mun$N_h * N_by_mun$S_h)
  
  N_by_mun <- N_by_mun %>%
    mutate(
      n_neyman = round(pool_target * (N_h * S_h) / sum_N_S),
      # Enforce pool-level min/max, and cannot exceed N_h
      n_alloc  = pmin(N_h, pmax(n_neyman, m_min_pool)),
      n_alloc  = pmin(n_alloc, m_max)
    )
  
  # 3) Compute sums
  planned_total <- sum(N_by_mun$n_alloc)
  cat("Planned pool total after Neyman + min/max: ", planned_total,
      " (pool target: ", pool_target, ")\n")
  
  rel_dev <- abs(planned_total - pool_target) / pool_target
  
  # 4) Draw stratified sample: SRSWOR within each municipality
  alloc_vec <- setNames(N_by_mun$n_alloc, N_by_mun[[mun_id]])
  mun_ids <- unique(Datos_with_id[[mun_id]])
  
  sampled_list <- vector("list", length(mun_ids))
  names(sampled_list) <- as.character(mun_ids)
  
  for (i in seq_along(mun_ids)) {
    mun_code <- mun_ids[i]
    mun_data <- Datos_with_id[Datos_with_id[[mun_id]] == mun_code, ]
    n_req    <- as.integer(alloc_vec[as.character(mun_code)])
    
    if (is.na(n_req) || n_req >= nrow(mun_data)) {
      sampled_list[[i]] <- mun_data
    } else if (n_req <= 0) {
      sampled_list[[i]] <- mun_data[0, ]
    } else {
      sampled_list[[i]] <- mun_data[
        sample.int(n = nrow(mun_data), size = n_req, replace = FALSE),
      ]
    }
  }
  
  sampled_df_with_id <- dplyr::bind_rows(sampled_list)
  actual_pool_n <- nrow(sampled_df_with_id)
  cat("Actual sampled rows (pool): ", actual_pool_n, "\n")
  
  # 5) Compute inclusion probabilities and sample weights (POOL-level)
  alloc_table <- N_by_mun %>%
    rename(!!mun_id := .data[[mun_id]]) %>%
    select(all_of(mun_id), N_h, S_h, n_neyman, n_alloc)
  
  sampled_df_with_id <- sampled_df_with_id %>%
    left_join(alloc_table, by = mun_id) %>%
    mutate(
      sample_prob   = ifelse(N_h > 0, n_alloc / N_h, NA_real_),
      sample_weight = ifelse(sample_prob > 0, 1 / sample_prob, NA_real_)
    )
  
  # 6) Tidy allocation table
  alloc_table <- alloc_table %>%
    mutate(
      frac          = n_alloc / N_h,
      prop_of_total = N_h / N_total
    )
  
  # 7) Stratified 50/50 split within each municipality
  #    n_test = floor(n_h / 2); n_train = ceiling(n_h / 2)
  #    Guarantees train has >= m_min per municipality (or all of N_h).
  if (return_test_set) {
    sampled_df_with_id <- sampled_df_with_id %>%
      group_by(.data[[mun_id]]) %>%
      mutate(.in_test = {
        n_h    <- dplyr::n()
        n_test <- floor(n_h / 2)
        sample(c(rep(TRUE, n_test), rep(FALSE, n_h - n_test)))
      }) %>%
      ungroup()
    
    test_df <- sampled_df_with_id %>%
      filter(.in_test) %>%
      select(-.row_id_temp, -.in_test)
    
    sampled_df <- sampled_df_with_id %>%
      filter(!.in_test) %>%
      select(-.row_id_temp, -.in_test)
    
    cat("Test set size:  ", nrow(test_df), "\n")
    cat("Train set size: ", nrow(sampled_df), "\n")
    
    # Diagnostic: verify the m_min guarantee on train (base R, robust)
    train_sizes <- as.integer(table(sampled_df[[mun_id]]))
    cat("Train per-municipality sizes: min =", min(train_sizes),
        ", max =", max(train_sizes),
        ", # municipalities =", length(train_sizes), "\n")
    
    if (min(train_sizes) < m_min) {
      n_below <- sum(train_sizes < m_min)
      cat("NOTE:", n_below, "municipalities have < m_min =", m_min,
          "in train (because their N_h < 2 * m_min =", m_min_pool, ").\n")
    }
    
  } else {
    test_df    <- NULL
    sampled_df <- sampled_df_with_id %>% select(-.row_id_temp)
    cat("No test split requested. Returning full pool as train.\n")
    cat("Train set size: ", nrow(sampled_df), "\n")
  }
  
  # Full data without temporary ID
  full_df <- Datos_with_id %>% select(-.row_id_temp)
  
  # 8) Return
  meta <- list(
    outcome_variable     = outcome_var,
    requested_n_train    = n_target,
    pool_target          = pool_target,
    m_min_train          = m_min,
    m_min_pool           = m_min_pool,
    planned_pool_total   = planned_total,
    actual_pool_sampled  = actual_pool_n,
    actual_train         = nrow(sampled_df),
    actual_test          = ifelse(return_test_set, nrow(test_df), NA),
    relative_deviation   = rel_dev
  )
  
  cat("\n=== SUMMARY ===\n")
  cat("Allocation method: Neyman (based on variance of", outcome_var, ")\n")
  cat("Full dataset:     ", nrow(full_df), "observations\n")
  cat("Training set:     ", nrow(sampled_df), "observations (target:", n_target, ")\n")
  if (return_test_set) {
    cat("Test set:         ", nrow(test_df), "observations\n")
    cat("Train + Test =", nrow(sampled_df) + nrow(test_df),
        "(should equal", actual_pool_n, ")\n")
  }
  
  return(list(
    full_df     = full_df,
    train_df    = sampled_df,
    test_df     = test_df,
    alloc_table = alloc_table,
    meta        = meta
  ))
}
# ---------------------------
# Define final number of samples (n) 
# ---------------------------

# sample size formula for estimating a population mean with the finite-population 
#correction (FPC).

# confidence level 95%
z <- 1.96
# Standard error of y (global score)
s <- sd(Datos_agrupados$punt_global)
# Total number of samples N 
N <- nrow(Datos_agrupados)
# Margin of error of 5 points
e <- 0.6

# Apply the formula

num <- (z^2) * N * (s^2)
den <- (z^2) * (s^2) + (N - 1) * (e^2)
n <- ceiling(num / den)

cat("Optimal number of samples:", n)

# ---------------------------
# Sampling Function Usage 
# ---------------------------
# params

#Note: put a value in n_target that is x2 of the desired training dataset size.
#test is ging to have half observations of the training df
n_target <- 50000  # training set size (test will be ~ same size)
m_min    <- 12     # minimum observations per municipality in train
m_max    <- Inf    # cap to large municipalities
seed     <- 777

res <- subsample_stratified_neyman(Datos = Datos_agrupados,
                                   mun_id = "cole_cod_mcpio_ubicacion",
                                   n_target = n_target,
                                   m_min = m_min,
                                   m_max = m_max,
                                   seed = seed,
                                   return_test_set = TRUE)

# Extract the 3 datasets
Datos_completos <- res$full_df       # Original complete data (around 200k obs)
Datos_train <- res$train_df          # Training sample (10k obs) 
Datos_test <- res$test_df            # Test set (10k obs)

# Metadata
Meta_info <- res$meta
Allocation_table <- res$alloc_table

# Verify
cat("\nVerification:\n")
cat("Full:", nrow(Datos_completos), "\n")
cat("Train:", nrow(Datos_train), "\n")
cat("Test:", nrow(Datos_test), "\n")
cat("Train + Test:", nrow(Datos_train) + nrow(Datos_test), "\n")
#cat("Match:", nrow(Datos_train) + nrow(Datos_test) == nrow(Datos_completos), "\n")


write.csv(Datos_completos, 
          file = paste0(path_to_save, "/Examen_Saber_11_2022_2_Datos_completos.csv"), 
          row.names = FALSE)
save(Datos_completos, file = paste0(path_to_save, "/Examen_Saber_11_2022_2_Datos_completos.RData"))

write.csv(Datos_train, 
          file = paste0(path_to_save, "/Examen_Saber_11_2022_2_Datos_train.csv"), 
          row.names = FALSE)
save(Datos_train, file = paste0(path_to_save, "/Examen_Saber_11_2022_2_Datos_train.RData"))

write.csv(Datos_test, 
          file = paste0(path_to_save, "/Examen_Saber_11_2022_2_Datos_test.csv"), 
          row.names = FALSE)
save(Datos_test, file = paste0(path_to_save, "/Examen_Saber_11_2022_2_Datos_test.RData"))


