

# dados de clima por fase com o script

setwd("C:\\Users\\elisv\\Downloads\\novo")

# --- climate data extract --- #
rm(list=ls()) 

if(!require(pacman)) install.packages("pacman")

pacman::p_load(tidyverse, nasapower, lubridate, furrr, future)


# ---- data load and processing ----
setwd("C:\\Users\\elisv\\Downloads\\novo")

#load
data <- read.csv("C:\\Users\\elisv\\Downloads\\novo\\dados_arroz.csv")

#processing
source("data_processing_elis2.R") 


data1<- data_processing_elis2(data)

dados <- data1

unique(dados$ST)
unique(dados$GEN)
unique(dados$LOCATION)
unique(dados$TYPE)
unique(dados$TRIAL)

table(dados$ST)
table(dados$GEN)
table(dados$LOCATION)
table(dados$TYPE)
table(dados$TRIAL)

#preparation for climate data
data_clim<- data1 %>%
  distinct(LATITUDE, LONGITUDE, DATE, PI_data, DTF_data, PM_data)


# ------------ functions ------------


variables <- c(                                                          #são 16.
  "T2M", "T2M_MIN", "T2M_MAX",   # Temperatura
  "PRECTOTCORR",                  # Precipitação
  "RH2M",                         # Umidade relativa
  "GWETROOT",                      # Umidade no solo (raízes)
  "ALLSKY_SFC_PAR_TOT",           # Radiação fotossinteticamente ativa
  "EVPTRNS",                       # Evapotranspiração
  "WS2M",                           # Velocidade do vento
  "ALLSKY_SFC_SW_DWN",  # Radiação solar total
  "ALLSKY_SRF_ALB",     # Albedo do solo
  "GWETPROF",            # Umidade no solo mais profundo
  "TSOIL1",              # Temperatura do solo superficial
  "CDD0", "HDD0",        # Dias de calor e frio
  "FROST_DAYS"           # Geadas
)


# saber quais variáveis pode puxar
nomes_variaveis_totais <- names(nasapower::query_parameters(community = "AG", temporal_api = "daily"))

#calculate stats function
calc_stats <- function(clim_data, stage_suffix) { 
  stats_df <- clim_data %>%
    summarise(across(all_of(variables), 
                     list(mean = ~mean(., na.rm = TRUE),
                          median = ~median(., na.rm = TRUE),
                          iqr = ~IQR(., na.rm = TRUE),
                          cumulative = ~sum(mean(., na.rm = TRUE)),
                          p90_prop = ~mean(. > quantile(., 0.90, na.rm = TRUE), na.rm = TRUE),
                          p10_prop = ~mean(. < quantile(., 0.10, na.rm = TRUE), na.rm = TRUE)),
                     .names = "{stage_suffix}_{.col}_{.fn}"))
  
  return(stats_df)
}

#get climate data function
get_climate_data <- function(lat, lon, start_date, end_date) { 
  
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  clim_data <- get_power(
    community = "AG",
    pars = variables,
    temporal_api = "daily",
    lonlat = c(lon, lat),
    dates = c(format(start_date, "%Y%m%d"), format(end_date, "%Y%m%d"))
  )
  
  return(as_tibble(clim_data))
}


# Função principal para processar cada linha do data_clim
process_climate_stats <- function(row) {
  lat <- row$LATITUDE
  lon <- row$LONGITUDE
  date <- row$DATE
  pi <- row$PI_data
  pf <- row$DTF_data
  pm <- row$PM_data
  
  tryCatch({
    #get data for each stage
    
    #vegetative (DATE to PI)
    veg_data <- get_climate_data(lat, lon, date, pi)
    
    #reproductive (PI to PF)
    repro_data <- get_climate_data(lat, lon, pi, pf)
    
    #grain_filling (PF to PM)
    gf_data <- get_climate_data(lat, lon, pf, pm)
    
    #stats calc for each stage
    veg_stats <- calc_stats(veg_data, "veg")
    repro_stats <- calc_stats(repro_data, "repro")
    gf_stats <- calc_stats(gf_data, "gf")
    
    #bind results in a single row
    result_row <- bind_cols(
      tibble(LATITUDE = lat, LONGITUDE = lon, DATE = date, PI = pi, PF = pf, PM = pm),
      veg_stats,
      repro_stats,
      gf_stats
    )
    
    return(result_row)
    
  }, error = function(e) {
    message(paste("Error in lat:", lat, "lon:", lon, "-", e$message))
    return(NULL)
  })
}



# ------ processamento em paralelo -----
parallel::detectCores()
plan(multisession, workers = 2)


#aplicar a função a cada linha de data_clim
climate_results <- future_map_dfr(
  1:nrow(data_clim),
  ~process_climate_stats(data_clim[.x, ]),
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)


# ----------- savar os dados -----------

# Salvar dataset completo
write_csv(climate_results, "climate_data_dezembro.csv") 


message(paste("Total de observações processadas:", nrow(climate_results)))
message(paste("Observações com erro:", nrow(data_clim) - nrow(climate_results)))

###############################################################################################################################

setwd("C:\\Users\\elisv\\Downloads\\novo")


#resultado da obtenção dos dados climáticos
dados_climaticos_por_fase <- read.csv("climate_data_dezembro.csv")

dados <- dados_climaticos_por_fase

unique(dados$ST)
unique(dados$GEN)
unique(dados$LOCATION)
unique(dados$TYPE)
unique(dados$TRIAL)

table(dados$ST)
table(dados$GEN)
table(dados$LOCATION)
table(dados$TYPE)
table(dados$TRIAL)


# juntarei posteriormente.
